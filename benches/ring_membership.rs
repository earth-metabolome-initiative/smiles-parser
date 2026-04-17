//! Benchmarks for [`Smiles::ring_membership`].
//!
//! The benchmark keeps a few single-molecule control cases and
//! opportunistically adds broader corpus-backed datasets when local fixtures
//! are available.

use core::str::FromStr;
use std::{
    fs::File,
    hint::black_box,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use flate2::read::GzDecoder;
use serde::Deserialize;
use smiles_parser::Smiles;

const SAMPLE_COUNT: usize = 9;
const MIN_BATCH_TIME: Duration = Duration::from_millis(100);
const DEFAULT_PUBCHEM_CORPUS_PATH: &str = "tests/CID-SMILES";
const DEFAULT_PUBCHEM_CORPUS_GZ_PATH: &str = "tests/CID-SMILES.gz";
const AROMATICITY_FIXTURE_PATH: &str =
    "tests/fixtures/aromaticity/corpus/pubchem_aromaticity_cases.json.gz";
const PUBCHEM_MIXED_SIZE: usize = 4_096;
const PUBCHEM_ACYCLIC_SIZE: usize = 2_048;
const PUBCHEM_RINGED_SIZE: usize = 2_048;
const PUBCHEM_LARGE_RINGED_SIZE: usize = 512;
const PUBCHEM_LARGE_RINGED_MIN_ATOMS: usize = 24;

const LARGE_POLYCYCLE_FRONTIER_CASE: &str = "CC1=C2CC3=C(C4=C5C=C3[C@@H]6C2=CC7=C1CC8=C(C9=C1C=C8[C@@H]7CCC[C@@H]2C3=C7CC8=C2C=C2C%10CCCC%11C%12=CC%13=C%14CC(=C7C)C(=C3)[C@@H]%13CCC[C@H]3C7=C%13CC%15=C3C=C3C(CCCC%16C%17=CC(=C(C4)C(=C%17CC4=C%16C=C%16[C@H](CCC6)C(=C7)C(=C%13C)CC%16=C4C)C)C5CCCC1C1=CC%10=C(CC2=C8C)C(=C1C9)C)C1=CC%11=C(CC%12=C%14C)C(=C1CC3=C%15C)C)C)C";

#[derive(Debug)]
struct BenchmarkCase {
    name: String,
    smiles_set: Vec<Smiles>,
    molecule_count: usize,
    atom_count: usize,
    bond_count: usize,
    ring_atom_count: usize,
    ring_bond_count: usize,
}

#[derive(Debug, Default)]
struct PubchemSamples {
    mixed: Vec<Smiles>,
    acyclic: Vec<Smiles>,
    ringed: Vec<Smiles>,
    large_ringed: Vec<Smiles>,
}

impl PubchemSamples {
    fn is_complete(&self) -> bool {
        self.mixed.len() >= PUBCHEM_MIXED_SIZE
            && self.acyclic.len() >= PUBCHEM_ACYCLIC_SIZE
            && self.ringed.len() >= PUBCHEM_RINGED_SIZE
            && self.large_ringed.len() >= PUBCHEM_LARGE_RINGED_SIZE
    }
}

#[derive(Debug, Deserialize)]
struct AromaticityCorpus {
    cases: Vec<AromaticityCorpusCase>,
}

#[derive(Debug, Deserialize)]
struct AromaticityCorpusCase {
    #[serde(default)]
    id: String,
    smiles: String,
    #[serde(default)]
    milestones: Vec<String>,
}

fn main() {
    let filter = std::env::args().skip(1).find(|argument| !argument.starts_with('-'));
    let cases = benchmark_cases();
    let mut matched_case_count = 0_usize;

    println!(
        "{:<32} {:>8} {:>10} {:>10} {:>10} {:>10} {:>16} {:>16}",
        "case", "mols", "atoms", "bonds", "ring_a", "ring_b", "median ns/mol", "best ns/mol",
    );

    for case in &cases {
        if let Some(filter) = &filter
            && !case.name.contains(filter)
        {
            continue;
        }
        matched_case_count += 1;
        print_case_result(case);
    }

    if matched_case_count == 0 {
        eprintln!("no benchmark case matched filter {filter:?}");
        std::process::exit(1);
    }
}

fn benchmark_cases() -> Vec<BenchmarkCase> {
    let mut cases = vec![
        single_case("acyclic_chain", "CCCCCCCCCCCC"),
        single_case("benzene", "C1=CC=CC=C1"),
        single_case("naphthalene", "Cc1cccc2ccccc12"),
        single_case("cubane", "C12C3C4C1C5C2C3C45"),
        single_case("polycycle_frontier", LARGE_POLYCYCLE_FRONTIER_CASE),
    ];

    if let Some(case) = load_ring_membership_fixture_case() {
        cases.push(case);
    }

    cases.extend(load_optional_pubchem_corpus_cases());
    cases
}

fn single_case(name: &str, source: &str) -> BenchmarkCase {
    benchmark_case(name.to_owned(), vec![parse_smiles_or_panic(name, source)])
}

fn benchmark_case(name: String, smiles_set: Vec<Smiles>) -> BenchmarkCase {
    assert!(!smiles_set.is_empty(), "benchmark case {name} must not be empty");

    let mut atom_count = 0_usize;
    let mut bond_count = 0_usize;
    let mut ring_atom_count = 0_usize;
    let mut ring_bond_count = 0_usize;

    for smiles in &smiles_set {
        atom_count += smiles.nodes().len();
        bond_count += smiles.number_of_bonds();
        let ring_membership = smiles.ring_membership();
        ring_atom_count += ring_membership.atom_ids().len();
        ring_bond_count += ring_membership.bond_edges().len();
    }

    BenchmarkCase {
        name,
        molecule_count: smiles_set.len(),
        atom_count,
        bond_count,
        ring_atom_count,
        ring_bond_count,
        smiles_set,
    }
}

fn load_ring_membership_fixture_case() -> Option<BenchmarkCase> {
    let path = Path::new(AROMATICITY_FIXTURE_PATH);
    if !path.is_file() {
        return None;
    }

    let reader = GzDecoder::new(File::open(path).ok()?);
    let corpus = serde_json::from_reader::<_, AromaticityCorpus>(reader).ok()?;
    let smiles_set = corpus
        .cases
        .into_iter()
        .filter(|case| case.milestones.iter().any(|milestone| milestone == "ring-membership"))
        .map(|case| parse_smiles_or_panic(case.id.as_str(), case.smiles.as_str()))
        .collect::<Vec<_>>();

    if smiles_set.is_empty() {
        return None;
    }

    Some(benchmark_case(format!("ring_fixture_{}", smiles_set.len()), smiles_set))
}

fn load_optional_pubchem_corpus_cases() -> Vec<BenchmarkCase> {
    let Some(path) = pubchem_corpus_path() else {
        return Vec::new();
    };

    let samples = collect_pubchem_samples(path.as_path());
    let mut cases = Vec::new();

    if !samples.mixed.is_empty() {
        cases.push(benchmark_case(format!("pubchem_mixed_{}", samples.mixed.len()), samples.mixed));
    }
    if !samples.acyclic.is_empty() {
        cases.push(benchmark_case(
            format!("pubchem_acyclic_{}", samples.acyclic.len()),
            samples.acyclic,
        ));
    }
    if !samples.ringed.is_empty() {
        cases.push(benchmark_case(
            format!("pubchem_ringed_{}", samples.ringed.len()),
            samples.ringed,
        ));
    }
    if !samples.large_ringed.is_empty() {
        cases.push(benchmark_case(
            format!("pubchem_large_ringed_{}", samples.large_ringed.len()),
            samples.large_ringed,
        ));
    }

    cases
}

fn pubchem_corpus_path() -> Option<PathBuf> {
    if let Some(path) = std::env::var_os("RING_MEMBERSHIP_BENCH_CORPUS").map(PathBuf::from) {
        assert!(
            path.is_file(),
            "RING_MEMBERSHIP_BENCH_CORPUS does not point to a file: {}",
            path.display()
        );
        return Some(path);
    }

    [DEFAULT_PUBCHEM_CORPUS_PATH, DEFAULT_PUBCHEM_CORPUS_GZ_PATH]
        .into_iter()
        .map(PathBuf::from)
        .find(|path| path.is_file())
}

fn collect_pubchem_samples(path: &Path) -> PubchemSamples {
    let mut samples = PubchemSamples::default();
    let reader = open_text_reader(path);

    for line in reader.lines() {
        let Ok(line) = line else {
            continue;
        };
        let Some((_, smiles_text)) = line.split_once('\t') else {
            continue;
        };
        let Ok(smiles) = Smiles::from_str(smiles_text) else {
            continue;
        };

        let atom_count = smiles.nodes().len();
        let ring_membership = smiles.ring_membership();
        let has_ring = !ring_membership.bond_edges().is_empty();

        if samples.mixed.len() < PUBCHEM_MIXED_SIZE {
            samples.mixed.push(smiles.clone());
        }
        if has_ring {
            if samples.ringed.len() < PUBCHEM_RINGED_SIZE {
                samples.ringed.push(smiles.clone());
            }
            if atom_count >= PUBCHEM_LARGE_RINGED_MIN_ATOMS
                && samples.large_ringed.len() < PUBCHEM_LARGE_RINGED_SIZE
            {
                samples.large_ringed.push(smiles.clone());
            }
        } else if samples.acyclic.len() < PUBCHEM_ACYCLIC_SIZE {
            samples.acyclic.push(smiles.clone());
        }

        if samples.is_complete() {
            break;
        }
    }

    samples
}

fn open_text_reader(path: &Path) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap_or_else(|error| {
        panic!("failed to open benchmark corpus {}: {error}", path.display())
    });

    if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

fn parse_smiles_or_panic(name: &str, source: &str) -> Smiles {
    Smiles::from_str(source).unwrap_or_else(|error| {
        panic!("failed to parse benchmark input {name}:\n{}", error.render(source))
    })
}

fn print_case_result(case: &BenchmarkCase) {
    let batch_iterations = calibrate_iterations(case);
    let sample_ns = sample_batch_ns(case, batch_iterations);
    let per_molecule_divisor = u128::from(batch_iterations)
        * u128::try_from(case.molecule_count)
            .unwrap_or_else(|_| unreachable!("usize always fits into u128"));
    let median_ns = sample_ns[SAMPLE_COUNT / 2] / per_molecule_divisor;
    let best_ns = sample_ns[0] / per_molecule_divisor;

    println!(
        "{:<32} {:>8} {:>10} {:>10} {:>10} {:>10} {:>16} {:>16}",
        case.name,
        case.molecule_count,
        case.atom_count,
        case.bond_count,
        case.ring_atom_count,
        case.ring_bond_count,
        median_ns,
        best_ns,
    );
}

fn calibrate_iterations(case: &BenchmarkCase) -> u64 {
    let mut iterations = 1_u64;
    loop {
        if measure_batch(case, iterations) >= MIN_BATCH_TIME {
            return iterations;
        }
        if iterations >= u64::MAX / 2 {
            return iterations;
        }
        iterations *= 2;
    }
}

fn sample_batch_ns(case: &BenchmarkCase, iterations: u64) -> [u128; SAMPLE_COUNT] {
    let mut sample_ns = [0_u128; SAMPLE_COUNT];
    for slot in &mut sample_ns {
        *slot = measure_batch(case, iterations).as_nanos();
    }
    sample_ns.sort_unstable();
    sample_ns
}

fn measure_batch(case: &BenchmarkCase, iterations: u64) -> Duration {
    let started_at = Instant::now();
    let mut checksum = 0_usize;

    for _ in 0..iterations {
        for smiles in &case.smiles_set {
            let ring_membership = smiles.ring_membership();
            checksum = checksum.wrapping_add(ring_membership.atom_ids().len());
            checksum = checksum.wrapping_add(ring_membership.bond_edges().len());
        }
    }

    black_box(checksum);
    started_at.elapsed()
}
