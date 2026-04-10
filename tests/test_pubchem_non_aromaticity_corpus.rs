//! Validate `smiles-parser` against an extracted full-PubChem `RDKit` default
//! non-aromatic corpus.

use std::{
    env, fs,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use flate2::read::GzDecoder;
use rayon::prelude::*;
use serde::Deserialize;
use smiles_parser::prelude::{AromaticityPolicy, AromaticityStatus, Smiles};

const DEFAULT_CORPUS_PATH: &str =
    "target/pubchem_aromaticity/default/pubchem_non_aromaticity_default.jsonl.gz";
const DEFAULT_MANIFEST_PATH: &str =
    "target/pubchem_aromaticity/default/pubchem_non_aromaticity_default.manifest.json";
const MDL_CORPUS_PATH: &str = "target/pubchem_aromaticity/mdl/pubchem_non_aromaticity_mdl.jsonl.gz";
const MDL_MANIFEST_PATH: &str =
    "target/pubchem_aromaticity/mdl/pubchem_non_aromaticity_mdl.manifest.json";
const SIMPLE_CORPUS_PATH: &str =
    "target/pubchem_aromaticity/simple/pubchem_non_aromaticity_simple.jsonl.gz";
const SIMPLE_MANIFEST_PATH: &str =
    "target/pubchem_aromaticity/simple/pubchem_non_aromaticity_simple.manifest.json";
const DEFAULT_VALIDATION_BATCH_SIZE: usize = 4_096;
const DEFAULT_PROGRESS_EVERY: usize = 100_000;

#[derive(Debug, Deserialize)]
struct PubChemNonAromaticManifest {
    non_aromatic_records: usize,
}

#[derive(Debug, Deserialize)]
struct PubChemNonAromaticRecord {
    cid: u64,
    smiles: String,
    atoms: usize,
    bonds: usize,
    aromatic_atom_ids: Vec<usize>,
    aromatic_bond_edges: Vec<[usize; 2]>,
    aromatic_atom_count: usize,
    aromatic_bond_count: usize,
}

#[test]
#[ignore = "This test streams a whole-PubChem non-aromatic corpus and checks that aromatic assignments stay empty."]
fn validate_rdkit_default_pubchem_non_aromaticity_corpus() -> Result<(), Box<dyn std::error::Error>>
{
    validate_pubchem_non_aromaticity_corpus(
        AromaticityPolicy::RdkitDefault,
        "PUBCHEM_NON_AROMATICITY",
        DEFAULT_CORPUS_PATH,
        DEFAULT_MANIFEST_PATH,
    )
}

#[test]
#[ignore = "This test streams a whole-PubChem non-aromatic corpus and checks that aromatic assignments stay empty for RDKit MDL."]
fn validate_rdkit_mdl_pubchem_non_aromaticity_corpus() -> Result<(), Box<dyn std::error::Error>> {
    validate_pubchem_non_aromaticity_corpus(
        AromaticityPolicy::RdkitMdl,
        "PUBCHEM_MDL_NON_AROMATICITY",
        MDL_CORPUS_PATH,
        MDL_MANIFEST_PATH,
    )
}

#[test]
#[ignore = "This test streams a whole-PubChem non-aromatic corpus and checks that aromatic assignments stay empty for RDKit Simple."]
fn validate_rdkit_simple_pubchem_non_aromaticity_corpus() -> Result<(), Box<dyn std::error::Error>>
{
    validate_pubchem_non_aromaticity_corpus(
        AromaticityPolicy::RdkitSimple,
        "PUBCHEM_SIMPLE_NON_AROMATICITY",
        SIMPLE_CORPUS_PATH,
        SIMPLE_MANIFEST_PATH,
    )
}

fn validate_pubchem_non_aromaticity_corpus(
    policy: AromaticityPolicy,
    env_prefix: &str,
    default_corpus_path: &str,
    default_manifest_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    require_release_build()?;

    let corpus_path = non_aromatic_corpus_path(env_prefix, default_corpus_path)?;
    let manifest_path = non_aromatic_manifest_path(env_prefix, default_manifest_path);
    let expected_total = manifest_path
        .exists()
        .then(|| load_manifest(&manifest_path))
        .transpose()?
        .map(|m| m.non_aromatic_records);

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
            "found {mismatch_count} PubChem non-aromaticity mismatches after validating {validated} records{rejects_suffix}"
        ))
        .into());
    }

    if let Some(limit) = limit {
        assert!(
            validated >= limit || stopped_for_mismatch_limit,
            "validation limit requested {limit}, but corpus only yielded {validated} records"
        );
    } else if let Some(expected_total) = expected_total {
        assert_eq!(
            validated, expected_total,
            "validated record count does not match the corpus manifest"
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
        .map(|line| validate_record_line_with_policy(line, policy).err())
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

fn validate_record_line_with_policy(line: &str, policy: AromaticityPolicy) -> Result<(), String> {
    let record = serde_json::from_str::<PubChemNonAromaticRecord>(line)
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

    if !record.aromatic_atom_ids.is_empty() || record.aromatic_atom_count != 0 {
        return Err(format!(
            "cid={} corpus non-aromatic record carries aromatic atoms: ids={:?} count={} smiles={}",
            record.cid, record.aromatic_atom_ids, record.aromatic_atom_count, record.smiles
        ));
    }

    if !record.aromatic_bond_edges.is_empty() || record.aromatic_bond_count != 0 {
        return Err(format!(
            "cid={} corpus non-aromatic record carries aromatic bonds: edges={:?} count={} smiles={}",
            record.cid, record.aromatic_bond_edges, record.aromatic_bond_count, record.smiles
        ));
    }

    let aromaticity = smiles.aromaticity_assignment_for(policy);
    if aromaticity.status() == AromaticityStatus::Unsupported {
        return Err(format!(
            "cid={} aromaticity status is unsupported: parser={:?} smiles={}",
            record.cid,
            aromaticity.status(),
            record.smiles
        ));
    }

    if !aromaticity.atom_ids().is_empty() {
        return Err(format!(
            "cid={} parser produced unexpected aromatic atoms: parser={:?} smiles={}",
            record.cid,
            aromaticity.atom_ids(),
            record.smiles
        ));
    }

    if !aromaticity.bond_edges().is_empty() {
        return Err(format!(
            "cid={} parser produced unexpected aromatic bonds: parser={:?} smiles={}",
            record.cid,
            aromaticity.bond_edges(),
            record.smiles
        ));
    }

    Ok(())
}

fn report_progress(
    validated: usize,
    progress_every: usize,
    started: Instant,
    next_progress_mark: &mut usize,
) {
    if progress_every == 0 || validated == 0 {
        return;
    }
    if validated < *next_progress_mark {
        return;
    }

    let elapsed = started.elapsed();
    let rate = records_per_second(validated, elapsed);
    eprintln!("validated={validated} elapsed={:.1}s rate={rate}/s", elapsed.as_secs_f64());
    while *next_progress_mark <= validated {
        *next_progress_mark = next_progress_mark.saturating_add(progress_every);
    }
}

fn require_release_build() -> Result<(), Box<dyn std::error::Error>> {
    if cfg!(debug_assertions) {
        return Err(io::Error::other(
            "PubChem non-aromaticity validation must be run with --release; debug builds are intentionally rejected",
        )
        .into());
    }
    Ok(())
}

fn non_aromatic_corpus_path(
    env_prefix: &str,
    default_path: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let path = env::var(format!("{env_prefix}_CORPUS_PATH"))
        .map_or_else(|_| PathBuf::from(default_path), PathBuf::from);

    if path.exists() {
        Ok(path)
    } else {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "PubChem non-aromatic corpus not found at {}. Run the extractor first or set PUBCHEM_NON_AROMATICITY_CORPUS_PATH.",
                path.display()
            ),
        )
        .into())
    }
}

fn non_aromatic_manifest_path(env_prefix: &str, default_path: &str) -> PathBuf {
    env::var(format!("{env_prefix}_MANIFEST_PATH"))
        .map_or_else(|_| PathBuf::from(default_path), PathBuf::from)
}

fn load_manifest(path: &Path) -> Result<PubChemNonAromaticManifest, Box<dyn std::error::Error>> {
    let bytes = fs::read(path)?;
    Ok(serde_json::from_slice(&bytes)?)
}

fn env_usize(name: &str) -> Option<usize> {
    env::var(name).ok().and_then(|value| value.parse().ok())
}

fn records_per_second(validated: usize, elapsed: std::time::Duration) -> usize {
    let elapsed_millis = elapsed.as_millis();
    if elapsed_millis == 0 {
        return validated;
    }

    let elapsed_millis = usize::try_from(elapsed_millis).unwrap_or(usize::MAX);
    validated.saturating_mul(1_000) / elapsed_millis.max(1)
}

#[test]
fn exact_pubchem_non_aromatic_record_validation_accepts_matching_empty_assignment() {
    let line = r#"{"cid":6324,"smiles":"CC","atoms":2,"bonds":1,"aromatic_atom_ids":[],"aromatic_bond_edges":[],"aromatic_atom_count":0,"aromatic_bond_count":0}"#;
    validate_record_line_with_policy(line, AromaticityPolicy::RdkitDefault)
        .expect("matching non-aromatic assignment should validate");
}

#[test]
fn exact_pubchem_non_aromatic_record_validation_rejects_unexpected_aromatic_atoms() {
    let line = r#"{"cid":241,"smiles":"c1ccccc1","atoms":6,"bonds":6,"aromatic_atom_ids":[],"aromatic_bond_edges":[],"aromatic_atom_count":0,"aromatic_bond_count":0}"#;
    let error = validate_record_line_with_policy(line, AromaticityPolicy::RdkitDefault)
        .expect_err("benzene should not validate as non-aromatic");
    assert!(error.contains("unexpected aromatic atoms"), "{error}");
}

#[test]
fn exact_pubchem_mdl_non_aromatic_record_validation_accepts_imidazole() {
    let line = r#"{"cid":795,"smiles":"c1ncc[nH]1","atoms":5,"bonds":5,"aromatic_atom_ids":[],"aromatic_bond_edges":[],"aromatic_atom_count":0,"aromatic_bond_count":0}"#;
    validate_record_line_with_policy(line, AromaticityPolicy::RdkitMdl)
        .expect("imidazole should validate as non-aromatic under MDL");
}

#[test]
fn exact_pubchem_simple_non_aromatic_record_validation_accepts_tropylium() {
    let line = r#"{"cid":7066,"smiles":"[CH+]1C=CC=CC=C1","atoms":7,"bonds":7,"aromatic_atom_ids":[],"aromatic_bond_edges":[],"aromatic_atom_count":0,"aromatic_bond_count":0}"#;
    validate_record_line_with_policy(line, AromaticityPolicy::RdkitSimple)
        .expect("tropylium should validate as non-aromatic under Simple");
}

#[test]
fn exact_pubchem_simple_non_aromatic_record_validation_accepts_cid_672() {
    let line = r#"{"cid":672,"smiles":"CC1=C(N(C2=NC(=O)NC(=O)C2=N1)CC(C(C(CO)O)O)O)C","atoms":23,"bonds":24,"aromatic_atom_ids":[],"aromatic_bond_edges":[],"aromatic_atom_count":0,"aromatic_bond_count":0}"#;
    validate_record_line_with_policy(line, AromaticityPolicy::RdkitSimple)
        .expect("CID 672 should validate as non-aromatic under Simple");
}
