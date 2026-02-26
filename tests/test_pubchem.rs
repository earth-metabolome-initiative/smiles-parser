//! Test suite for validating SMILES parsing against PubChem Data.
//!
//! # Running Tests
//!
//! To run this test (validates SMILES in the PubChem dataset), ensure that:
//!
//! ```
//! cargo test --release --test test_pubchem_inchi_validation -- --ignored --nocapture
//! ```

use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Write},
    path::Path,
    result,
};

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use smiles_parser::{smiles::Smiles, token::{Token, TokenWithSpan}};

/// Structure representing a PubChem compound as a SMILES string
#[derive(Debug, Deserialize)]
struct SmilesPubChemCompound {
    /// The id for the SMILES
    id: u64,
    /// Smiles String
    smiles: String,
}

fn validate_pubchem_smiles(file_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut csv_reader =
        ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(reader);
    let pb = ProgressBar::new(123_458_626);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {post}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    let start = std::time::Instant::now();
    for result in csv_reader.deserialize::<SmilesPubChemCompound>() {
        let result = result?;
        pb.inc(1);

        let smiles_str = &result.smiles;
        match smiles_str.parse::<Smiles>() {
            Ok(_) => todo!(),
            Err(_) => todo!(),
        }
    }
    Ok(())
}
