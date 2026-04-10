//! Criterion benchmarks for SMILES parsing, ring analysis, rendering, and
//! aromaticity assignment.
//!
//! This bench uses two corpora:
//! - `tests/fixtures/benchmark_smiles_corpus.txt` for the legacy parse/display
//!   and implicit-hydrogen groups
//! - `tests/fixtures/aromaticity/corpus/pubchem_aromaticity_cases.json` for the
//!   PubChem parse/ring/SSSR/aromaticity groups
//!
//! The benchmarked operations are kept separate:
//! - parse-only: `Smiles::from_str(...)`
//! - display-only: `Smiles::to_string()` on already parsed molecules
//! - implicit-only: `Smiles::implicit_hydrogen_counts()` on already parsed
//!   molecules
//!
//! The aromaticity-oriented groups use the frozen PubChem fixture corpus in
//! `tests/fixtures/aromaticity/corpus/pubchem_aromaticity_cases.json`.

use std::{hint::black_box, str::FromStr};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use serde::Deserialize;
use smiles_parser::smiles::Smiles;

const CORPUS: &str = include_str!("../tests/fixtures/benchmark_smiles_corpus.txt");
const PUBCHEM_CORPUS: &str =
    include_str!("../tests/fixtures/aromaticity/corpus/pubchem_aromaticity_cases.json");

#[derive(Debug, Deserialize)]
struct PubChemCorpus {
    cases: Vec<PubChemCase>,
}

#[derive(Debug, Deserialize)]
struct PubChemCase {
    smiles: String,
}

fn corpus_smiles() -> Vec<&'static str> {
    CORPUS.lines().map(str::trim).filter(|line| !line.is_empty()).collect()
}

fn pubchem_smiles() -> Vec<String> {
    serde_json::from_str::<PubChemCorpus>(PUBCHEM_CORPUS)
        .expect("PubChem benchmark corpus should deserialize")
        .cases
        .into_iter()
        .map(|case| case.smiles)
        .collect()
}

fn bench_parse_only(c: &mut Criterion) {
    let corpus = corpus_smiles();
    let mut group = c.benchmark_group("smiles_parse");
    group.throughput(Throughput::Elements(corpus.len() as u64));

    group.bench_function("from_str", |b| {
        b.iter(|| {
            let mut atom_total = 0usize;
            for smiles in corpus.iter().copied() {
                let parsed = Smiles::from_str(black_box(smiles)).expect("fixture should parse");
                atom_total += parsed.nodes().len();
            }
            black_box(atom_total)
        });
    });

    group.finish();
}

fn bench_pubchem_parse_only(c: &mut Criterion) {
    let corpus = pubchem_smiles();
    let mut group = c.benchmark_group("pubchem_parse");
    group.throughput(Throughput::Elements(corpus.len() as u64));

    group.bench_function("from_str", |b| {
        b.iter(|| {
            let mut atom_total = 0usize;
            for smiles in &corpus {
                let parsed = Smiles::from_str(black_box(smiles)).expect("fixture should parse");
                atom_total += parsed.nodes().len();
            }
            black_box(atom_total)
        });
    });

    group.finish();
}

fn bench_implicit_only(c: &mut Criterion) {
    let corpus = corpus_smiles();
    let parsed = corpus
        .iter()
        .map(|smiles| Smiles::from_str(smiles).expect("fixture should parse"))
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("implicit_hydrogens");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("counts", |b| {
        b.iter(|| {
            let mut hydrogen_total = 0usize;
            for molecule in &parsed {
                let counts = black_box(molecule).implicit_hydrogen_counts();
                let per_molecule_total = counts.into_iter().map(usize::from).sum::<usize>();
                hydrogen_total += black_box(per_molecule_total);
            }
            black_box(hydrogen_total)
        });
    });

    group.finish();
}

fn bench_display_only(c: &mut Criterion) {
    let corpus = corpus_smiles();
    let parsed = corpus
        .iter()
        .map(|smiles| Smiles::from_str(smiles).expect("fixture should parse"))
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("smiles_display");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("to_string", |b| {
        b.iter(|| {
            let mut rendered_total = 0usize;
            for molecule in &parsed {
                let rendered = black_box(molecule).to_string();
                rendered_total += black_box(rendered.len());
            }
            black_box(rendered_total)
        });
    });

    group.finish();
}

fn bench_ring_membership(c: &mut Criterion) {
    let parsed = pubchem_smiles()
        .into_iter()
        .map(|smiles| Smiles::from_str(&smiles).expect("PubChem fixture should parse"))
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("pubchem_ring_membership");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("ring_membership", |b| {
        b.iter(|| {
            let mut ring_member_total = 0usize;
            for molecule in &parsed {
                let membership = black_box(molecule).ring_membership();
                ring_member_total += black_box(membership.atom_ids().len());
                ring_member_total += black_box(membership.bond_edges().len());
            }
            black_box(ring_member_total)
        });
    });

    group.finish();
}

fn bench_symm_sssr(c: &mut Criterion) {
    let parsed = pubchem_smiles()
        .into_iter()
        .map(|smiles| Smiles::from_str(&smiles).expect("PubChem fixture should parse"))
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("pubchem_symm_sssr");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("symm_sssr", |b| {
        b.iter(|| {
            let mut cycle_total = 0usize;
            for molecule in &parsed {
                let result = black_box(molecule).symm_sssr_result();
                cycle_total += black_box(result.cycles().len());
            }
            black_box(cycle_total)
        });
    });

    group.finish();
}

fn bench_aromaticity_assignment(c: &mut Criterion) {
    let parsed = pubchem_smiles()
        .into_iter()
        .map(|smiles| Smiles::from_str(&smiles).expect("PubChem fixture should parse"))
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("pubchem_aromaticity");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("assignment", |b| {
        b.iter(|| {
            let mut aromatic_total = 0usize;
            for molecule in &parsed {
                let assignment = black_box(molecule).aromaticity_assignment();
                aromatic_total += black_box(assignment.atom_ids().len());
                aromatic_total += black_box(assignment.bond_edges().len());
            }
            black_box(aromatic_total)
        });
    });

    group.finish();
}

criterion_group!(
    smiles_benches,
    bench_parse_only,
    bench_pubchem_parse_only,
    bench_display_only,
    bench_implicit_only,
    bench_ring_membership,
    bench_symm_sssr,
    bench_aromaticity_assignment
);
criterion_main!(smiles_benches);
