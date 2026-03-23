//! Test suite for validating SMILES parsing against PubChem Data.
//!
//! # Running Tests
//!
//! To run this test (validates SMILES in the PubChem dataset), ensure that:
//!
//! ```
//! cargo test --release --test validate_pubchem_smiles
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
    /// SMILES String
    smiles: String,
}

#[test]
#[ignore = "This test downloads a ~6.79 GB file and is time-consuming."]
fn validate_pubchem_smiles() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("tests/CID-SMILES.gz")?;
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
        if count.is_multiple_of(10_000) {
            pb.set_position(count);
        }

        let smiles_str = &result.smiles;
        let smiles = smiles_str.parse::<Smiles>();
        match smiles {
            Ok(smiles) => {
                let rendered = smiles.render().unwrap_or_else(|e| {
                    panic!(
                        "Failed to render SMILES.\nRecord: {}\nPubChem ID: {}\nOriginal:\n{}\n{:?}",
                        count, result.id, smiles_str, e
                    )
                });

                let reparsed = rendered.parse::<Smiles>().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse rendered SMILES.\nPubChem ID: {}\nOriginal:\n{}\nRendered:\n{}\n{}",
                        result.id,
                        smiles_str,
                        rendered,
                        e.render(&rendered)
                    )
                });
                assert_eq!(smiles.nodes().len(), reparsed.nodes().len());
                assert_eq!(smiles.edges().len(), reparsed.edges().len());
            }
            Err(err) => panic!("{}\n{}", result.id, err.render(smiles_str)),
        }
    }
    Ok(())
}
