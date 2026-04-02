//! Criterion benchmarks for SMILES parsing and implicit-hydrogen counting.
//!
//! This bench uses the corpus in `tests/fixtures/benchmark_smiles_corpus.txt`.
//! The two benchmarked operations are kept separate:
//! - parse-only: `Smiles::from_str(...)`
//! - implicit-only: `Smiles::implicit_hydrogen_counts()` on already parsed
//!   molecules

use std::{hint::black_box, str::FromStr};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use smiles_parser::smiles::Smiles;

const CORPUS: &str = include_str!("../tests/fixtures/benchmark_smiles_corpus.txt");

fn corpus_smiles() -> Vec<&'static str> {
    CORPUS.lines().map(str::trim).filter(|line| !line.is_empty()).collect()
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

criterion_group!(smiles_benches, bench_parse_only, bench_implicit_only);
criterion_main!(smiles_benches);
