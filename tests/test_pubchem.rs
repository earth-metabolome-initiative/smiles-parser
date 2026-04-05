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
    env,
    fmt::Write as _,
    fs::{self, File},
    io::{self, BufReader, Read, Write},
    path::{Path, PathBuf},
};

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use reqwest::blocking::Client;
use serde::Deserialize;
use smiles_parser::{SmilesError, smiles::Smiles};

const PUBCHEM_CID_SMILES_PATH: &str = "tests/CID-SMILES";
const PUBCHEM_CID_SMILES_GZ_PATH: &str = "tests/CID-SMILES.gz";
const PUBCHEM_CID_SMILES_URL: &str =
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz";
const PUBCHEM_RECORD_COUNT: u64 = 123_458_626;
const PUBCHEM_BATCH_SIZE: usize = 16_384;

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
    require_release_build()?;
    let reader = open_pubchem_reader()?;
    let mut csv_reader =
        ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(reader);
    let pb = records_progress_bar();

    while let Some(batch) = next_pubchem_batch(&mut csv_reader)? {
        let batch_len = batch.len() as u64;
        if let Some(error) = batch.par_iter().find_map_any(validate_pubchem_record) {
            return Err(io::Error::other(error).into());
        }
        pb.inc(batch_len);
    }
    pb.finish_and_clear();
    Ok(())
}

#[test]
#[ignore = "This test uses a ~6.79 GB file and is time-consuming."]
fn collect_pubchem_fails() -> Result<(), Box<dyn std::error::Error>> {
    require_release_build()?;
    let reader = open_pubchem_reader()?;
    let mut csv_reader =
        ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(reader);
    let pb = records_progress_bar();
    let mut rejects = String::new();
    while let Some(batch) = next_pubchem_batch(&mut csv_reader)? {
        let batch_len = batch.len() as u64;
        for reject in batch.par_iter().filter_map(collect_pubchem_reject).collect::<Vec<_>>() {
            writeln!(&mut rejects, "{reject}")
                .expect("writing a string has failed and should not have");
        }
        pb.inc(batch_len);
    }
    pb.finish_and_clear();
    let mut file = File::create("invalid_hydrogens.csv")?;
    file.write_all(rejects.as_bytes())?;
    Ok(())
}

fn records_progress_bar() -> ProgressBar {
    let pb = ProgressBar::new(PUBCHEM_RECORD_COUNT);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    pb
}

fn require_release_build() -> Result<(), Box<dyn std::error::Error>> {
    if cfg!(debug_assertions) {
        return Err(io::Error::other(
            "PubChem tests must be run with --release; debug builds are intentionally rejected",
        )
        .into());
    }
    Ok(())
}

fn next_pubchem_batch<R: Read>(
    csv_reader: &mut csv::Reader<R>,
) -> Result<Option<Vec<SmilesPubChemCompound>>, csv::Error> {
    let mut batch = Vec::with_capacity(PUBCHEM_BATCH_SIZE);
    for result in csv_reader.deserialize::<SmilesPubChemCompound>() {
        batch.push(result?);
        if batch.len() == PUBCHEM_BATCH_SIZE {
            return Ok(Some(batch));
        }
    }
    if batch.is_empty() { Ok(None) } else { Ok(Some(batch)) }
}

fn validate_pubchem_record(result: &SmilesPubChemCompound) -> Option<String> {
    let smiles_str = &result.smiles;
    let smiles = match smiles_str.parse::<Smiles>() {
        Ok(smiles) => smiles,
        Err(err) => {
            return Some(format!("id: {}\n SMILES:\n{}", result.id, err.render(smiles_str)));
        }
    };

    let rendered = match smiles.render() {
        Ok(rendered) => rendered,
        Err(error) => {
            return Some(format!(
                "Failed to render SMILES.\nPubChem ID:\n{}\nOriginal:\n{}\n{:?}",
                result.id, smiles_str, error
            ));
        }
    };

    let reparsed = match rendered.parse::<Smiles>() {
        Ok(reparsed) => reparsed,
        Err(error) => {
            return Some(format!(
                "Failed to parse rendered SMILES.\nPubChem ID:\n{}\nOriginal:\n{}\nRendered:\n{}\n{}",
                result.id,
                smiles_str,
                rendered,
                error.render(&rendered)
            ));
        }
    };

    if smiles.nodes().len() != reparsed.nodes().len()
        || smiles.number_of_bonds() != reparsed.number_of_bonds()
    {
        return Some(format!(
            "Round-trip topology mismatch.\nPubChem ID:\n{}\nOriginal:\n{}\nRendered:\n{}",
            result.id, smiles_str, rendered
        ));
    }

    None
}

fn collect_pubchem_reject(result: &SmilesPubChemCompound) -> Option<String> {
    let smiles_str = &result.smiles;
    let err = smiles_str.parse::<Smiles>().err()?;
    if err.smiles_error() != SmilesError::InvalidHydrogenWithExplicitHydrogensFound {
        return None;
    }
    let replacement = replace_hydrogen_hcount(smiles_str);
    Some(format!("{}, {}, {}", result.id, smiles_str, replacement))
}

fn open_pubchem_reader() -> Result<BufReader<Box<dyn Read>>, Box<dyn std::error::Error>> {
    let path = pubchem_fixture_path()?;
    let reader: Box<dyn Read> = Box::new(File::open(path)?);
    Ok(BufReader::new(reader))
}

fn pubchem_fixture_path() -> Result<PathBuf, Box<dyn std::error::Error>> {
    if let Ok(path) = env::var("PUBCHEM_CID_SMILES_PATH") {
        return Ok(PathBuf::from(path));
    }

    let plain = Path::new(PUBCHEM_CID_SMILES_PATH);
    if plain.exists() {
        return Ok(plain.to_path_buf());
    }

    let gz = Path::new(PUBCHEM_CID_SMILES_GZ_PATH);
    if gz.exists() {
        expand_pubchem_fixture(gz, plain)?;
        return Ok(plain.to_path_buf());
    }

    download_pubchem_fixture(gz)?;
    expand_pubchem_fixture(gz, plain)?;
    Ok(plain.to_path_buf())
}

fn download_pubchem_fixture(destination: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let client = Client::builder().build()?;
    let mut response = client.get(PUBCHEM_CID_SMILES_URL).send()?.error_for_status()?;
    let total_bytes = response.content_length().unwrap_or(0);

    let pb =
        if total_bytes > 0 { ProgressBar::new(total_bytes) } else { ProgressBar::new_spinner() };
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} downloading CID-SMILES.gz [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let mut file = File::create(destination)?;
    let mut buffer = vec![0_u8; 64 * 1024];
    loop {
        let read = response.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        file.write_all(&buffer[..read])?;
        pb.inc(read as u64);
    }
    file.flush()?;
    pb.finish_with_message("downloaded CID-SMILES.gz");
    Ok(())
}

fn expand_pubchem_fixture(
    compressed: &Path,
    destination: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    if destination.exists() {
        return Ok(());
    }

    let total_bytes = fs::metadata(compressed)?.len();
    let pb = ProgressBar::new(total_bytes);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} expanding CID-SMILES [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let partial = destination.with_extension("part");
    let mut input = GzDecoder::new(File::open(compressed)?);
    let mut output = File::create(&partial)?;
    let mut buffer = vec![0_u8; 64 * 1024];
    loop {
        let read = input.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        output.write_all(&buffer[..read])?;
        pb.inc(read as u64);
    }
    output.flush()?;
    fs::rename(partial, destination)?;
    pb.finish_with_message("expanded CID-SMILES");
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
