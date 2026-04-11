//! Kekulization round-trip checks against committed fixtures and full `PubChem`
//! aromatic corpora.

use std::{
    env, fs,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use flate2::read::GzDecoder;
use geometric_traits::traits::SparseValuedMatrixRef;
use rayon::prelude::*;
use serde::Deserialize;
use smiles_parser::bond::Bond;
use smiles_parser::prelude::{AromaticityPolicy, AromaticityStatus, Smiles};

const PUBCHEM_AROMATIC_CASES_PATH: &str =
    "tests/fixtures/aromaticity/corpus/pubchem_aromaticity_cases.json";
const DEFAULT_CORPUS_PATH: &str =
    "target/pubchem_aromaticity/default/pubchem_aromaticity_default.jsonl.gz";
const MDL_CORPUS_PATH: &str = "target/pubchem_aromaticity/mdl/pubchem_aromaticity_mdl.jsonl.gz";
const SIMPLE_CORPUS_PATH: &str =
    "target/pubchem_aromaticity/simple/pubchem_aromaticity_simple.jsonl.gz";
const DEFAULT_VALIDATION_BATCH_SIZE: usize = 4_096;
const DEFAULT_PROGRESS_EVERY: usize = 100_000;

#[derive(Debug, Deserialize)]
struct AromaticCorpus {
    cases: Vec<AromaticCorpusCase>,
}

#[derive(Debug, Deserialize)]
struct AromaticCorpusCase {
    id: String,
    smiles: String,
}

#[derive(Debug, Deserialize)]
struct PubChemAromaticRecord {
    cid: u64,
    smiles: String,
    atoms: usize,
    bonds: usize,
    aromatic_atom_ids: Vec<usize>,
    aromatic_bond_edges: Vec<[usize; 2]>,
}

#[test]
fn kekulization_roundtrips_pubchem_aromatic_cases_under_default_perception() {
    let corpus_path = Path::new(PUBCHEM_AROMATIC_CASES_PATH);
    let corpus = serde_json::from_str::<AromaticCorpus>(
        &fs::read_to_string(corpus_path).expect("pubchem aromatic corpus should be readable"),
    )
    .expect("pubchem aromatic corpus should parse");

    let mut failures = Vec::new();

    for case in corpus.cases {
        let smiles: Smiles = case.smiles.parse().expect("fixture smiles should parse");
        let perceived = match smiles.perceive_aromaticity() {
            Ok(perception) => perception,
            Err(error) => {
                failures.push(format!("{}: initial perception failed: {error}", case.id));
                continue;
            }
        };
        let expected_atom_ids = perceived.assignment().atom_ids().to_vec();
        let expected_bond_edges = perceived.assignment().bond_edges().to_vec();

        let kekulized = match perceived.aromaticized().kekulize() {
            Ok(kekulized) => kekulized,
            Err(error) => {
                failures.push(format!("{}: kekulization failed: {error}", case.id));
                continue;
            }
        };

        let reperceived = match kekulized.perceive_aromaticity() {
            Ok(perception) => perception,
            Err(error) => {
                failures.push(format!("{}: reperception failed: {error}", case.id));
                continue;
            }
        };

        if reperceived.assignment().atom_ids() != expected_atom_ids
            || reperceived.assignment().bond_edges() != expected_bond_edges
        {
            failures.push(format!(
                "{}: roundtrip mismatch atoms {:?} != {:?}, bonds {:?} != {:?}",
                case.id,
                reperceived.assignment().atom_ids(),
                expected_atom_ids,
                reperceived.assignment().bond_edges(),
                expected_bond_edges,
            ));
        }
    }

    assert!(failures.is_empty(), "kekulization roundtrip mismatches:\n{}", failures.join("\n"));
}

#[test]
#[ignore = "This test streams the whole-PubChem aromatic corpus and validates perceive->kekulize->perceive roundtrips for RDKit Default."]
fn validate_rdkit_default_pubchem_kekulization_roundtrip_corpus()
-> Result<(), Box<dyn std::error::Error>> {
    validate_pubchem_kekulization_roundtrip_corpus(
        AromaticityPolicy::RdkitDefault,
        "PUBCHEM_KEKULIZATION_DEFAULT",
        DEFAULT_CORPUS_PATH,
    )
}

#[test]
#[ignore = "This test streams the whole-PubChem aromatic corpus and validates perceive->kekulize->perceive roundtrips for RDKit MDL."]
fn validate_rdkit_mdl_pubchem_kekulization_roundtrip_corpus()
-> Result<(), Box<dyn std::error::Error>> {
    validate_pubchem_kekulization_roundtrip_corpus(
        AromaticityPolicy::RdkitMdl,
        "PUBCHEM_KEKULIZATION_MDL",
        MDL_CORPUS_PATH,
    )
}

#[test]
#[ignore = "This test streams the whole-PubChem aromatic corpus and validates perceive->kekulize->perceive roundtrips for RDKit Simple."]
fn validate_rdkit_simple_pubchem_kekulization_roundtrip_corpus()
-> Result<(), Box<dyn std::error::Error>> {
    validate_pubchem_kekulization_roundtrip_corpus(
        AromaticityPolicy::RdkitSimple,
        "PUBCHEM_KEKULIZATION_SIMPLE",
        SIMPLE_CORPUS_PATH,
    )
}

fn validate_pubchem_kekulization_roundtrip_corpus(
    policy: AromaticityPolicy,
    env_prefix: &str,
    default_corpus_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    require_release_build()?;

    let corpus_path = env::var(format!("{env_prefix}_CORPUS"))
        .map_or_else(|_| PathBuf::from(default_corpus_path), PathBuf::from);
    let limit = env_usize(&format!("{env_prefix}_VALIDATE_LIMIT"));
    let progress_every =
        env_usize(&format!("{env_prefix}_PROGRESS_EVERY")).unwrap_or(DEFAULT_PROGRESS_EVERY);
    let batch_size =
        env_usize(&format!("{env_prefix}_BATCH_SIZE")).unwrap_or(DEFAULT_VALIDATION_BATCH_SIZE);
    let mismatch_limit = env_usize(&format!("{env_prefix}_MAX_MISMATCHES"));
    let rejects_path = env::var(format!("{env_prefix}_REJECTS_PATH")).ok().map(PathBuf::from);
    let mut rejects_writer =
        rejects_path.as_ref().map(|path| File::create(path).map(BufWriter::new)).transpose()?;

    let reader = BufReader::new(GzDecoder::new(File::open(&corpus_path)?));
    let started = Instant::now();
    let mut validated = 0_usize;
    let mut mismatch_count = 0_usize;
    let mut stopped_for_mismatch_limit = false;
    let mut batch = Vec::with_capacity(batch_size);
    let mut next_progress_mark = progress_every;

    'records: for line in reader.lines() {
        batch.push(line?);
        if batch.len() == batch_size {
            mismatch_count += handle_batch_result(
                &batch,
                policy,
                &mut rejects_writer,
                mismatch_limit,
                mismatch_count,
                &mut stopped_for_mismatch_limit,
            )?;
            validated += batch.len();
            batch.clear();
            report_progress(validated, progress_every, started, &mut next_progress_mark);
            if stopped_for_mismatch_limit {
                break 'records;
            }
            if limit.is_some_and(|limit| validated >= limit) {
                break 'records;
            }
        }
    }

    if !batch.is_empty()
        && !stopped_for_mismatch_limit
        && limit.is_none_or(|limit| validated < limit)
    {
        mismatch_count += handle_batch_result(
            &batch,
            policy,
            &mut rejects_writer,
            mismatch_limit,
            mismatch_count,
            &mut stopped_for_mismatch_limit,
        )?;
        validated += batch.len();
        report_progress(validated, progress_every, started, &mut next_progress_mark);
    }

    if let Some(writer) = rejects_writer.as_mut() {
        writer.flush()?;
    }

    if mismatch_count > 0 {
        let rejects_suffix = rejects_path
            .as_ref()
            .map(|path| format!("; see {}", path.display()))
            .unwrap_or_default();
        return Err(io::Error::other(format!(
            "found {mismatch_count} PubChem kekulization roundtrip mismatches after validating {validated} records{rejects_suffix}"
        ))
        .into());
    }

    if let Some(limit) = limit {
        assert!(
            validated >= limit || stopped_for_mismatch_limit,
            "validation limit requested {limit}, but corpus only yielded {validated} records"
        );
    }

    let elapsed = started.elapsed();
    let rate = records_per_second(validated, elapsed);
    eprintln!(
        "completed validated={validated} elapsed={:.1}s rate={rate}/s",
        elapsed.as_secs_f64()
    );

    Ok(())
}

fn handle_batch_result(
    batch: &[String],
    policy: AromaticityPolicy,
    rejects_writer: &mut Option<BufWriter<File>>,
    mismatch_limit: Option<usize>,
    mismatches_seen_so_far: usize,
    stopped_for_mismatch_limit: &mut bool,
) -> Result<usize, Box<dyn std::error::Error>> {
    let mismatches = batch
        .par_iter()
        .map(|line| validate_roundtrip_record_line_with_policy(line, policy).err())
        .collect::<Vec<_>>();

    if rejects_writer.is_none() {
        if let Some(error) = mismatches.into_iter().flatten().next() {
            return Err(io::Error::other(error).into());
        }
        return Ok(0);
    }

    if let Some(writer) = rejects_writer.as_mut() {
        let mut written = 0_usize;
        for mismatch in mismatches.into_iter().flatten() {
            writer.write_all(mismatch.as_bytes())?;
            writer.write_all(b"\n")?;
            written += 1;
            if mismatch_limit.is_some_and(|limit| mismatches_seen_so_far + written >= limit) {
                *stopped_for_mismatch_limit = true;
                break;
            }
        }
        Ok(written)
    } else {
        Ok(0)
    }
}

fn validate_roundtrip_record_line_with_policy(
    line: &str,
    policy: AromaticityPolicy,
) -> Result<(), String> {
    let record = serde_json::from_str::<PubChemAromaticRecord>(line)
        .map_err(|error| format!("failed to parse corpus record JSON: {error}\nline={line}"))?;

    let smiles = record.smiles.parse::<Smiles>().map_err(|error| {
        format!("cid={} failed to parse SMILES:\n{}", record.cid, error.render(&record.smiles))
    })?;

    if smiles.nodes().len() != record.atoms {
        return Err(format!(
            "cid={} atom-count mismatch: parser={} corpus={} smiles={}",
            record.cid,
            smiles.nodes().len(),
            record.atoms,
            record.smiles
        ));
    }

    if smiles.number_of_bonds() != record.bonds {
        return Err(format!(
            "cid={} bond-count mismatch: parser={} corpus={} smiles={}",
            record.cid,
            smiles.number_of_bonds(),
            record.bonds,
            record.smiles
        ));
    }

    let perceived = smiles.perceive_aromaticity_for(policy).map_err(|error| {
        format!(
            "cid={} initial aromaticity perception failed: {error} smiles={}",
            record.cid, record.smiles
        )
    })?;
    if perceived.assignment().status() == AromaticityStatus::Unsupported {
        return Err(format!(
            "cid={} initial aromaticity status is unsupported: parser={:?} smiles={}",
            record.cid,
            perceived.assignment().status(),
            record.smiles
        ));
    }
    if perceived.assignment().atom_ids() != record.aromatic_atom_ids
        || perceived.assignment().bond_edges() != record.aromatic_bond_edges
    {
        return Err(format!(
            "cid={} initial perception mismatch: atoms {:?} != {:?}, bonds {:?} != {:?}, smiles={}",
            record.cid,
            perceived.assignment().atom_ids(),
            record.aromatic_atom_ids,
            perceived.assignment().bond_edges(),
            record.aromatic_bond_edges,
            record.smiles
        ));
    }

    let kekulized = perceived.aromaticized().kekulize().map_err(|error| {
        format!("cid={} kekulization failed: {error} smiles={}", record.cid, record.smiles)
    })?;
    let residual_aromatic_bonds = kekulized
        .bond_matrix()
        .sparse_entries()
        .filter(|((row, column), entry)| *row < *column && entry.bond() == Bond::Aromatic)
        .count();
    if residual_aromatic_bonds > 0 {
        return Err(format!(
            "cid={} kekulized graph still contains {residual_aromatic_bonds} aromatic bonds: smiles={}",
            record.cid, record.smiles
        ));
    }

    let reperceived = kekulized.perceive_aromaticity_for(policy).map_err(|error| {
        format!(
            "cid={} reperception failed after kekulization: {error} smiles={}",
            record.cid, record.smiles
        )
    })?;
    if reperceived.assignment().status() == AromaticityStatus::Unsupported {
        return Err(format!(
            "cid={} reperception status is unsupported after kekulization: parser={:?} smiles={}",
            record.cid,
            reperceived.assignment().status(),
            record.smiles
        ));
    }

    if reperceived.assignment().atom_ids() != record.aromatic_atom_ids
        || reperceived.assignment().bond_edges() != record.aromatic_bond_edges
    {
        return Err(format!(
            "cid={} roundtrip mismatch: atoms {:?} != {:?}, bonds {:?} != {:?}, smiles={}",
            record.cid,
            reperceived.assignment().atom_ids(),
            record.aromatic_atom_ids,
            reperceived.assignment().bond_edges(),
            record.aromatic_bond_edges,
            record.smiles
        ));
    }

    Ok(())
}

fn require_release_build() -> Result<(), Box<dyn std::error::Error>> {
    if cfg!(debug_assertions) {
        return Err(io::Error::other(
            "full-PubChem kekulization roundtrip validation must run with `cargo test --release`",
        )
        .into());
    }
    Ok(())
}

fn env_usize(name: &str) -> Option<usize> {
    env::var(name).ok().and_then(|value| value.parse::<usize>().ok())
}

fn report_progress(
    validated: usize,
    progress_every: usize,
    started: Instant,
    next_progress_mark: &mut usize,
) {
    if progress_every == 0 {
        return;
    }
    while validated >= *next_progress_mark {
        let elapsed = started.elapsed();
        let rate = records_per_second(validated, elapsed);
        eprintln!("validated={validated} elapsed={:.1}s rate={rate}/s", elapsed.as_secs_f64());
        *next_progress_mark += progress_every;
    }
}

fn records_per_second(validated: usize, elapsed: std::time::Duration) -> usize {
    let elapsed_millis = elapsed.as_millis();
    if elapsed_millis == 0 {
        0
    } else {
        let validated = u128::try_from(validated).unwrap_or(u128::MAX);
        let rounded_rate =
            validated.saturating_mul(1000).saturating_add(elapsed_millis / 2) / elapsed_millis;
        usize::try_from(rounded_rate).unwrap_or(usize::MAX)
    }
}
