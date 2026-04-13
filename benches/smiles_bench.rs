//! Criterion benchmarks for SMILES parsing, canonicalization, ring analysis,
//! rendering, and aromaticity assignment.
//!
//! This bench uses two corpora:
//! - `tests/fixtures/benchmark_smiles_corpus.txt` for the legacy parse/display
//!   and implicit-hydrogen groups
//! - `tests/fixtures/aromaticity/corpus/pubchem_benchmark_cases.json.gz` for
//!   the PubChem parse/ring/SSSR/aromaticity groups
//!
//! The benchmarked operations are kept separate:
//! - parse-only: `Smiles::from_str(...)`
//! - canonical-labeling only: `Smiles::canonical_labeling()` on already parsed
//!   molecules
//! - canonicalize-only: `Smiles::canonicalize()` on already parsed molecules
//! - display-only: `Smiles::to_string()` on already parsed molecules
//! - implicit-only: `Smiles::implicit_hydrogen_counts()` on already parsed
//!   molecules
//!
//! The aromaticity-oriented groups use the frozen PubChem benchmark corpus in
//! `tests/fixtures/aromaticity/corpus/pubchem_benchmark_cases.json.gz`.

use std::{hint::black_box, str::FromStr};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use flate2::read::GzDecoder;
use serde::Deserialize;
use smiles_parser::smiles::Smiles;

const CORPUS: &str = include_str!("../tests/fixtures/benchmark_smiles_corpus.txt");
const PUBCHEM_CORPUS: &[u8] =
    include_bytes!("../tests/fixtures/aromaticity/corpus/pubchem_benchmark_cases.json.gz");
const HEAVY_PUBCHEM_CASES: &[(&str, &str)] = &[
    ("cid_1", "CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C"),
    ("naphthalene", "c1ccccc1c2ccccc2"),
    (
        "cid_158446604",
        "B(C1=CC=CC=C1)(O)O.C.CCOC(=O)[C@@]12C[C@H]1/C=C\\CCCCC[C@@H](C(=O)N3C[C@@H](C[C@H]3C(=O)C2)OC4=NC5=C(C=CC=C5N4C(C)C)C6=NC(=CS6)C7CCCC7)CC8=CC=CC=C8.CCOC(=O)[C@@]12C[C@H]1/C=C\\CCCCC[C@@H](C(=O)N3C[C@@H](C[C@H]3C(=O)C2)OC4=NC5=C(C=CC=C5N4C(C)C)C6=NC(=CS6)C7CCCC7)N.CC#CC#CC#CC#CC#CC#CC#CC#CC#CC.C[C@H]1CCCCC/C=C\\[C@@H]2C[C@]2(CC(=O)[C@@H]3C[C@H](CN3C1=O)OC4=NC5=C(C=CC=C5N4C(C)C)C6=NC(=CS6)C7CCCC7)C(=O)O.CC(C)N1C2=CC=CC(=C2N=C1O[C@@H]3C[C@H]4C(=O)C[C@@]5(C[C@H]5/C=C\\CCCCC[C@@H](C(=O)N4C3)CC6=CC=CC=C6)C(=O)O)C7=NC(=CS7)C8CCCC8.CF.CF.CF.[OH-].[Na+]",
    ),
    (
        "cid_158446616",
        "CC(=O)C1=CC=CC(=C1)C2=CN=C3N2C=CN=C3NCC4=CC=C(C=C4)S(=O)(=O)N.COC1=C(C=C(C=C1)C2=CN=C3N2C=CN=C3NCC4=CC=C(C=C4)S(=O)(=O)C)OC.COC1=CC=CC(=C1)C2=CN=C3N2C=CN=C3NCC4=CC=C(C=C4)S(=O)(=O)N.CS(=O)(=O)C1=CC=C(C=C1)CNC2=NC=CN3C2=NC=C3C4=CC=C(C=C4)O.CS(=O)(=O)NC1=CC=CC(=C1)C2=CN=C3N2C=CN=C3NCCC4=CC=NC=C4.C1=CC(=CC=C1CNC2=NC=CN3C2=NC=C3C4=CC(=C(C=C4)O)F)S(=O)(=O)N",
    ),
    (
        "cid_158446617",
        "CCC1=CC(=C2C=C3C4=C5C2=C1C=CC5=C(C=C4CCC3(C)C)N(C6=CC=CC=C6)C7=C(C=CC=C7C)C)N(C8=CC=CC=C8)C9=C(C=CC=C9C)C.CC1=CC=C(C=C1)N(C2=CC(=C3C=CC4=C(C=C5CCC(C6=C5C4=C3C2=C6)(C)C)N(C7=CC=C(C=C7)C)C8=C(C=CC=C8C)C)C(C)C)C9=C(C=CC=C9C)C.CC1=C(C(=CC=C1)C)N(C2=CC=C(C=C2)C(C)C)C3=C4C=C5C6=C7C4=C(C=CC7=C(C=C6CCC5(C)C)N(C8=CC=C(C=C8)C(C)C)C9=C(C=CC=C9C)C)C(=C3)C",
    ),
    (
        "cid_158446619",
        "[CH3-].[CH3-].[CH3-].CC(C)(C)CC(=C)CN(C)CCOCCC1=CC(=C2N1N=CN=C2OC3=C(C=C(C=C3)NC(=O)NC(=O)CC4=CC=C(C=C4)F)F)CN5CCC(CC5)N.C=CC1=NC=CN1CCCCCC2=CC(=C3N2N=CN=C3OC4=C(C=C(C=C4)NC(=O)NC(=O)CC5=CC=C(C=C5)F)F)CN6CCC(CC6)N.C1CN(CCC1N)CC2=C3C(=NC=NN3C(=C2)CCCCCC4=NSC(=N4)Cl)OC5=C(C=C(C=C5)NC(=O)NC(=O)CC6=CC=C(C=C6)F)F.[V].[V].[V]",
    ),
    (
        "cid_158446620",
        "C.C1C=NC2=C(C1=O)N=CN2[C@H]3[C@H]([C@H]4[C@H](O3)COP(=S)(O[C@@H]5[C@H]([C@@H](COP(=S)(O4)O)O[C@H]5N6C7=C(C(=O)NC(=N7)N)N=N6)F)O)F.C1C=NC2=C(C1=O)N=NN2[C@H]3[C@H]([C@H]4[C@H](O3)COP(=S)(O[C@@H]5[C@@H]([C@@H](COP(=O)(O4)O)O[C@H]5N6C=NC7=C6N=C(NC7=O)N)F)O)F.C1C=NC2=C(C1=O)N=NN2[C@H]3[C@H]([C@H]4[C@H](O3)COP(=O)(O[C@@H]5[C@@H]([C@@H](COP(=O)(O4)O)O[C@H]5N6C=NC7=C6N=C(NC7=O)N)F)S)F",
    ),
    (
        "cid_158446625",
        "CC1(C2=CC=CC=C2C3=C1C=C(C=C3)N(C4=CC=C(C=C4)C5=CC6=C(C=C5)OC7=CC=CC=C76)C8=CC=CC9=C8C1=CC=CC=C1C=C9)C.CC1(C2=CC=CC=C2C3=C1C=C(C=C3)N(C4=CC=C(C=C4)C5=CC=CC6=C5OC7=CC=CC=C67)C8=CC=CC9=C8C1=CC=CC=C1C=C9)C.C1=CC=C(C=C1)C2(C3=CC=CC=C3C4=C2C=C(C=C4)N(C5=CC=C(C=C5)C6=CC7=C(C=C6)N(C8=CC=CC=C87)C9=CC=CC=C9)C1=CC=CC2=C1C1=CC=CC=C1C=C2)C1=CC=CC=C1",
    ),
];

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

fn parse_corpus<'a, I>(smiles: I, context: &str) -> Vec<Smiles>
where
    I: IntoIterator<Item = &'a str>,
{
    smiles
        .into_iter()
        .map(|smiles| {
            Smiles::from_str(smiles).unwrap_or_else(|error| {
                panic!("{context} fixture should parse: {}", error.render(smiles))
            })
        })
        .collect()
}

fn pubchem_smiles() -> Vec<String> {
    serde_json::from_reader::<_, PubChemCorpus>(GzDecoder::new(PUBCHEM_CORPUS))
        .expect("PubChem benchmark corpus should deserialize")
        .cases
        .into_iter()
        .map(|case| case.smiles)
        .collect()
}

fn parse_heavy_pubchem_cases() -> Vec<(&'static str, Smiles)> {
    HEAVY_PUBCHEM_CASES
        .iter()
        .map(|(cid, smiles)| {
            (
                *cid,
                Smiles::from_str(smiles).unwrap_or_else(|error| {
                    panic!("heavy PubChem fixture {cid} should parse: {}", error.render(smiles))
                }),
            )
        })
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
    let parsed = parse_corpus(corpus_smiles(), "benchmark corpus");

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
    let parsed = parse_corpus(corpus_smiles(), "benchmark corpus");

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

fn bench_canonicalization(c: &mut Criterion) {
    let parsed = parse_corpus(corpus_smiles(), "benchmark corpus");

    let mut group = c.benchmark_group("smiles_canonicalization");
    group.throughput(Throughput::Elements(parsed.len() as u64));

    group.bench_function("canonical_labeling", |b| {
        b.iter(|| {
            let mut position_total = 0usize;
            for molecule in &parsed {
                let labeling = black_box(molecule).canonical_labeling();
                position_total += black_box(labeling.order().len());
                position_total += black_box(labeling.new_index_of_old_node().len());
            }
            black_box(position_total)
        });
    });

    group.bench_function("canonicalize", |b| {
        b.iter(|| {
            let mut node_total = 0usize;
            for molecule in &parsed {
                let canonicalized = black_box(molecule).canonicalize();
                node_total += black_box(canonicalized.nodes().len());
            }
            black_box(node_total)
        });
    });

    group.finish();
}

fn bench_heavy_pubchem_canonicalization(c: &mut Criterion) {
    let parsed = parse_heavy_pubchem_cases();

    let mut group = c.benchmark_group("pubchem_heavy_canonicalization");
    group.sample_size(30);

    for (cid, molecule) in &parsed {
        group.bench_with_input(BenchmarkId::new("canonical_labeling", cid), molecule, |b, m| {
            b.iter(|| {
                let labeling = black_box(m).canonical_labeling();
                black_box(labeling.order().len() + labeling.new_index_of_old_node().len())
            });
        });

        group.bench_with_input(BenchmarkId::new("canonicalize", cid), molecule, |b, m| {
            b.iter(|| {
                let canonicalized = black_box(m).canonicalize();
                black_box(canonicalized.nodes().len())
            });
        });
    }

    group.finish();
}

#[cfg(feature = "benchmarking")]
fn bench_heavy_pubchem_canonicalization_stages(c: &mut Criterion) {
    let parsed = parse_heavy_pubchem_cases();

    let mut group = c.benchmark_group("pubchem_heavy_canonicalization_stages");
    group.sample_size(20);

    for (cid, molecule) in &parsed {
        group.bench_with_input(
            BenchmarkId::new("canonicalization_normal_form", cid),
            molecule,
            |b, m| {
                b.iter(|| black_box(m).benchmark_canonicalization_normal_form());
            },
        );
        group.bench_with_input(BenchmarkId::new("collapse_explicit_h", cid), molecule, |b, m| {
            b.iter(|| black_box(m).benchmark_collapse_removable_explicit_hydrogens());
        });
        group.bench_with_input(BenchmarkId::new("stereo_normal_form", cid), molecule, |b, m| {
            b.iter(|| black_box(m).benchmark_stereo_normal_form());
        });
        group.bench_with_input(
            BenchmarkId::new("exact_canonical_labeling", cid),
            molecule,
            |b, m| {
                b.iter(|| {
                    let labeling = black_box(m).benchmark_exact_canonical_labeling();
                    black_box(labeling.order().len() + labeling.new_index_of_old_node().len())
                });
            },
        );
        group.bench_with_input(BenchmarkId::new("exact_canonicalize", cid), molecule, |b, m| {
            b.iter(|| black_box(m).benchmark_exact_canonicalize());
        });
        group.bench_with_input(BenchmarkId::new("canonicalization_step", cid), molecule, |b, m| {
            b.iter(|| black_box(m).benchmark_canonicalization_step());
        });
        group.bench_with_input(BenchmarkId::new("canonicalize", cid), molecule, |b, m| {
            b.iter(|| black_box(m).canonicalize());
        });
    }

    group.finish();
}

#[cfg(not(feature = "benchmarking"))]
fn bench_heavy_pubchem_canonicalization_stages(_c: &mut Criterion) {}

fn bench_ring_membership(c: &mut Criterion) {
    let pubchem = pubchem_smiles();
    let parsed = parse_corpus(pubchem.iter().map(String::as_str), "PubChem benchmark corpus");

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
    let pubchem = pubchem_smiles();
    let parsed = parse_corpus(pubchem.iter().map(String::as_str), "PubChem benchmark corpus");

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
    let pubchem = pubchem_smiles();
    let parsed = parse_corpus(pubchem.iter().map(String::as_str), "PubChem benchmark corpus");

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
    bench_canonicalization,
    bench_heavy_pubchem_canonicalization,
    bench_heavy_pubchem_canonicalization_stages,
    bench_implicit_only,
    bench_ring_membership,
    bench_symm_sssr,
    bench_aromaticity_assignment
);
criterion_main!(smiles_benches);
