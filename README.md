# smiles-rs
[![crates.io](https://img.shields.io/crates/v/smiles-rs.svg)](https://crates.io/crates/smiles-rs)
[![docs.rs](https://img.shields.io/docsrs/smiles-rs)](https://docs.rs/smiles-rs)
[![downloads](https://img.shields.io/crates/d/smiles-rs.svg)](https://crates.io/crates/smiles-rs)
[![Rust CI](https://github.com/earth-metabolome-initiative/smiles-parser/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/earth-metabolome-initiative/smiles-parser/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/earth-metabolome-initiative/smiles-parser/graph/badge.svg)](https://codecov.io/gh/earth-metabolome-initiative/smiles-parser)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/earth-metabolome-initiative/smiles-parser/blob/main/LICENSE)
[![MSRV](https://img.shields.io/badge/rustc-1.92%2B-orange.svg)](https://blog.rust-lang.org/)

Parses SMILES strings into molecular graphs, following the [OpenSMILES specification](http://opensmiles.org/opensmiles.html). `no_std` with `alloc`, no unsafe code.

## What it does

- **Canonicalization**: `canonicalize` and `canonical_labeling` give a canonical SMILES and a canonical atom ordering.
- **Aromaticity perception**: `perceive_aromaticity` under the RDKit default, MDL and simple models, plus `kekulize` back to alternating bonds.
- **Maximum common edge subgraph**: `mces` compares two molecules through [`geometric-traits`](https://crates.io/crates/geometric-traits).
- **Ring analysis**: symmetrized SSSR, ring membership, and fragment and connected component decomposition.
- **Atom environments**: `atom_environment` yields radius-bounded neighbourhoods, the basis for MAP4-style fingerprints.
- **Stereochemistry**: tetrahedral and double bond configuration, preserved across canonicalization.
- **Molecular formulas**: conversion into [`molecular-formulas`](https://crates.io/crates/molecular-formulas) types.
- **Wildcard SMILES**: `WildcardSmiles` accepts `*` atoms, and conversion back into `Smiles` is fallible.
- **Public corpora**: PubChem, ZINC20, COCONUT, LOTUS and MassSpecGym stream from a local cache behind the `datasets` feature.

## Example

```rust
use core::str::FromStr;

use molecular_formulas::prelude::ChemicalFormula;
use smiles_rs::prelude::Smiles;

let ethanol = Smiles::from_str("CCO")?;

assert_eq!(ethanol.nodes().len(), 3);
assert_eq!(ethanol.number_of_bonds(), 2);
assert_eq!(ethanol.render(), "CCO");

let formula: ChemicalFormula<u32, i32> = ChemicalFormula::from(&ethanol);
assert_eq!(formula.to_string(), "C₂H₆O");
# Ok::<(), smiles_rs::SmilesErrorWithSpan>(())
```

## Features

| Feature | Effect |
| ------- | ------ |
| `datasets` | Fetches and streams public SMILES corpora, requires `std`. See [`datasets`](https://docs.rs/smiles-rs/latest/smiles_rs/datasets/). |
| `fuzzing` | Exposes the parser internals the fuzz targets drive. |
