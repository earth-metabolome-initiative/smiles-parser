//! Test suite for validating SMILES parsing against PubChem Data.
//!
//! # Running Tests
//!
//! To run this test (validates SMILES in the PubChem dataset), ensure that:
//!
//! ```
//! cargo test --release --test validate_pubchem_smiles
//! ```

use std::{
    fmt::Write as _,
    fs::File,
    io::{BufReader, Write},
};

use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use smiles_parser::{SmilesError, smiles::Smiles};

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
    let file = File::open("tests/CID-SMILES")?;
    let reader = BufReader::new(file);

    let mut csv_reader =
        ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(reader);
    let pb = ProgressBar::new(123_458_626);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}]
{post}/{len} ({eta})",
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
                        "Failed to render SMILES.\nRecord: {}\nPubChem ID:
{}\nOriginal:\n{}\n{:?}",
                        count, result.id, smiles_str, e
                    )
                });

                let reparsed = rendered.parse::<Smiles>().unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse rendered SMILES.\nPubChem ID:
{}\nOriginal:\n{}\nRendered:\n{}\n{}",
                        result.id,
                        smiles_str,
                        rendered,
                        e.render(&rendered)
                    )
                });
                assert_eq!(smiles.nodes().len(), reparsed.nodes().len());
                assert_eq!(smiles.number_of_bonds(), reparsed.number_of_bonds());
            }
            Err(err) => {
                panic!("id: {}\n SMILES:\n{}", result.id, err.render(smiles_str));
            }
        }
    }
    Ok(())
}

#[test]
#[ignore = "This test uses a ~6.79 GB file and is time-consuming."]
fn collect_pubchem_fails() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("tests/CID-SMILES")?;
    let reader = BufReader::new(file);

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
    let mut rejects = String::new();
    for result in csv_reader.deserialize::<SmilesPubChemCompound>() {
        count += 1;
        let result = result?;
        if count.is_multiple_of(10_000) {
            pb.set_position(count);
        }
        let smiles_str = &result.smiles;
        let smiles = smiles_str.parse::<Smiles>();
        match smiles {
            Ok(_) => {}
            Err(err) => {
                if err.smiles_error() == SmilesError::InvalidHydrogenWithExplicitHydrogensFound {
                    let id = result.id;
                    let replacement = replace_hydrogen_hcount(smiles_str);
                    writeln!(&mut rejects, "{id}, {smiles_str}, {replacement}")
                        .expect("writing a string has failed and should not have");
                }
            }
        }
    }
    let mut file = File::create("invalid_hydrogens.csv")?;
    file.write_all(rejects.as_bytes())?;
    Ok(())
}

fn replace_hydrogen_hcount(smiles: &str) -> String {
    let bytes = smiles.as_bytes();
    let mut out = String::new();
    let mut i = 0;

    while i < bytes.len() {
        if bytes[i] == b'[' {
            // Find the closing ']'
            if let Some(rel_end) = smiles[i + 1..].find(']') {
                let end = i + 1 + rel_end;
                let inner = &smiles[i + 1..end];

                if let Some(expanded) = expand_bracket_h_atom(inner) {
                    out.push_str(&expanded);
                } else {
                    out.push_str(&smiles[i..=end]);
                }

                i = end + 1;
            }
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }

    out
}

fn expand_bracket_h_atom(inner: &str) -> Option<String> {
    let b = inner.as_bytes();
    let mut i = 0;

    // isotope?
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }

    // symbol must be exactly H
    if i >= b.len() || b[i] != b'H' {
        return None;
    }

    // Avoid matching elements like He, Hf, Hg, Ho, Hs...
    if i + 1 < b.len() && b[i + 1].is_ascii_lowercase() {
        return None;
    }

    let first_h_end = i + 1;

    // Need an hcount marker H immediately after the symbol H
    if first_h_end >= b.len() || b[first_h_end] != b'H' {
        return None;
    }

    // Parse optional digits after that second H
    let mut j = first_h_end + 1;
    while j < b.len() && b[j].is_ascii_digit() {
        j += 1;
    }

    let extra_h_count =
        if j == first_h_end + 1 { 1 } else { inner[first_h_end + 1..j].parse::<usize>().ok()? };

    // Keep everything except the hcount marker and its digits
    // Example:
    // inner = "3HH2-2:2"
    // prefix = "3H"
    // suffix = "-2:2"
    let prefix = &inner[..first_h_end];
    let suffix = &inner[j..];

    let mut result = String::new();
    result.push('[');
    result.push_str(prefix);
    result.push_str(suffix);
    result.push(']');

    for _ in 0..extra_h_count {
        result.push_str("[H]");
    }

    Some(result)
}
