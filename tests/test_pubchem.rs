//! Test suite for validating SMILES parsing against PubChem Data.
//!
//! # Running Tests
//!
//! To run this test (validates SMILES in the PubChem dataset), ensure that:
//!
//! ```
//! cargo test --release --test validate_pubchem_smiles -- --ignored --nocapture
//! ```

use std::{fs::File, io::BufReader};

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use smiles_parser::smiles::Smiles;

/// Structure representing a PubChem compound as a SMILES string
#[derive(Debug, Deserialize)]
struct SmilesPubChemCompound {
    /// The id for the SMILES
    id: u64,
    /// Smiles String
    smiles: String,
}

#[test]
#[ignore = "This test downloads a ~6.79 GB file and is time-consuming."]
fn validate_pubchem_smiles() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("tests/CID-SMILES.gz")?;
    // let file = File::open(file_path)?;
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
    let mut count: u64 = 0;
    for result in csv_reader.deserialize::<SmilesPubChemCompound>() {
        count += 1;
        let result = result?;
       if count % 10_000 == 0 {
        pb.set_position(count);
       }

        let smiles_str = &result.smiles;
        if let Err(err) = smiles_str.parse::<Smiles>() {
            panic!("{}\n{}", result.id, err.render(smiles_str));
        }
    }
    Ok(())
}
