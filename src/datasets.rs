//! Programmatic dataset downloads with local caching.
//!
//! This module is available behind the `datasets` cargo feature.
//!
//! The download API is intentionally path-oriented:
//!
//! - dataset definitions provide stable metadata
//! - downloads are cached on disk
//! - callers can choose whether gzip archives stay compressed, are
//!   decompressed, or both
//!
//! ```no_run
//! use smiles_parser::datasets::{DatasetFetchOptions, DatasetSource, GzipMode, PUBCHEM_SMILES};
//!
//! let artifact = PUBCHEM_SMILES.fetch_with_options(&DatasetFetchOptions {
//!     gzip_mode: GzipMode::Decompress,
//!     ..DatasetFetchOptions::default()
//! })?;
//!
//! println!("{}", artifact.path().display());
//! # Ok::<(), smiles_parser::DatasetError>(())
//! ```
//!
//! ```no_run
//! use smiles_parser::datasets::{PUBCHEM_SMILES, SmilesDatasetSource};
//!
//! let mut smiles = PUBCHEM_SMILES.iter_smiles()?;
//! if let Some(first) = smiles.next() {
//!     println!("{}", first?);
//! }
//! # Ok::<(), smiles_parser::DatasetError>(())
//! ```

use alloc::{borrow::ToOwned, boxed::Box, string::String};
use std::{
    env,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, IsTerminal, Read, Write},
    path::{Path, PathBuf},
    time::Duration,
};

use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::blocking::Client;
use thiserror::Error;

const DOWNLOAD_USER_AGENT: &str = concat!("smiles-parser/", env!("CARGO_PKG_VERSION"));
const PROGRESS_TICK_INTERVAL: Duration = Duration::from_millis(100);
const BYTE_PROGRESS_TEMPLATE: &str = "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})";
const SPINNER_PROGRESS_TEMPLATE: &str =
    "{spinner:.green} {msg} [{elapsed_precise}] {bytes} ({bytes_per_sec})";

/// Compression used by the upstream dataset artifact.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum DatasetCompression {
    /// The upstream file is not compressed.
    None,
    /// The upstream file is gzip-compressed.
    Gzip,
}

/// Cache behavior for dataset fetches.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum CacheMode {
    /// Reuse an existing cached artifact when present.
    #[default]
    UseCache,
    /// Always redownload the upstream artifact.
    Redownload,
}

/// How gzip-compressed datasets should be materialized locally.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum GzipMode {
    /// Keep the cached gzip file as-is.
    #[default]
    KeepCompressed,
    /// Materialize and return a decompressed copy.
    Decompress,
    /// Keep the gzip cache and also materialize a decompressed copy.
    KeepBoth,
}

/// Options controlling dataset fetch and cache behavior.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetFetchOptions {
    /// Override the cache directory. When `None`, a per-user cache directory is
    /// selected automatically.
    pub cache_dir: Option<PathBuf>,
    /// Whether to reuse existing cached files or redownload them.
    pub cache_mode: CacheMode,
    /// How gzip-compressed datasets should be stored locally.
    pub gzip_mode: GzipMode,
}

impl Default for DatasetFetchOptions {
    fn default() -> Self {
        Self {
            cache_dir: None,
            cache_mode: CacheMode::UseCache,
            gzip_mode: GzipMode::KeepCompressed,
        }
    }
}

/// A materialized dataset artifact on disk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetArtifact {
    dataset_id: &'static str,
    path: PathBuf,
    compressed_path: Option<PathBuf>,
    decompressed_path: Option<PathBuf>,
    was_downloaded: bool,
    was_decompressed: bool,
}

impl DatasetArtifact {
    /// Returns the dataset identifier that produced this artifact.
    #[must_use]
    pub fn dataset_id(&self) -> &'static str {
        self.dataset_id
    }

    /// Returns the primary path callers should consume.
    ///
    /// This is the decompressed path when one was requested, otherwise the
    /// cached upstream file path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Returns the cached compressed path, when present.
    #[must_use]
    pub fn compressed_path(&self) -> Option<&Path> {
        self.compressed_path.as_deref()
    }

    /// Returns the cached decompressed path, when present.
    #[must_use]
    pub fn decompressed_path(&self) -> Option<&Path> {
        self.decompressed_path.as_deref()
    }

    /// Returns whether the upstream artifact was downloaded during this call.
    #[must_use]
    pub fn was_downloaded(&self) -> bool {
        self.was_downloaded
    }

    /// Returns whether a decompressed copy was created during this call.
    #[must_use]
    pub fn was_decompressed(&self) -> bool {
        self.was_decompressed
    }
}

/// Errors raised while fetching and materializing datasets.
#[derive(Debug, Error)]
pub enum DatasetError {
    /// The HTTP download failed.
    #[error("failed to download dataset from {url}: {source}")]
    Download {
        /// The upstream URL.
        url: &'static str,
        /// The underlying transport or HTTP error.
        #[source]
        source: reqwest::Error,
    },
    /// The on-disk materialization step failed.
    #[error("failed to access dataset path {path}: {source}")]
    Io {
        /// The path being read or written.
        path: PathBuf,
        /// The underlying filesystem error.
        #[source]
        source: io::Error,
    },
    /// The dataset contents did not match the expected record layout.
    #[error("failed to parse dataset {dataset_id} at line {line_number}: {message}")]
    Format {
        /// The dataset identifier.
        dataset_id: &'static str,
        /// The 1-based line number within the materialized dataset file.
        line_number: usize,
        /// A human-readable explanation of the malformed record.
        message: String,
    },
}

/// Metadata for a downloadable dataset source.
pub trait DatasetSource {
    /// Stable dataset identifier used for cache subdirectories and diagnostics.
    fn id(&self) -> &'static str;

    /// Upstream URL for the dataset artifact.
    fn url(&self) -> &'static str;

    /// File name of the cached upstream artifact.
    fn file_name(&self) -> &'static str;

    /// File name of the decompressed artifact.
    ///
    /// For uncompressed sources this should equal [`Self::file_name`].
    fn extracted_file_name(&self) -> &'static str {
        self.file_name()
    }

    /// Compression used by the upstream dataset artifact.
    fn compression(&self) -> DatasetCompression {
        DatasetCompression::None
    }

    /// Fetches the dataset using default options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the download fails or the cached artifact
    /// cannot be materialized on disk.
    fn fetch(&self) -> Result<DatasetArtifact, DatasetError> {
        self.fetch_with_options(&DatasetFetchOptions::default())
    }

    /// Fetches the dataset using explicit options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the download fails or the requested cached
    /// artifact layout cannot be materialized on disk.
    fn fetch_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetArtifact, DatasetError> {
        fetch_dataset(self, options)
    }
}

/// A dataset source that can stream SMILES strings directly.
pub trait SmilesDatasetSource: DatasetSource {
    /// Opens a streaming iterator over the dataset SMILES using default fetch
    /// options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the dataset cannot be fetched or if the
    /// materialized dataset cannot be opened for streaming.
    fn iter_smiles(&self) -> Result<DatasetSmilesIter, DatasetError> {
        self.iter_smiles_with_options(&DatasetFetchOptions::default())
    }

    /// Opens a streaming iterator over the dataset SMILES using explicit fetch
    /// options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the dataset cannot be fetched or if the
    /// materialized dataset cannot be opened for streaming.
    fn iter_smiles_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesIter, DatasetError>;
}

/// A streaming iterator over SMILES strings extracted from a dataset file.
pub struct DatasetSmilesIter {
    dataset_id: &'static str,
    path: PathBuf,
    reader: Box<dyn BufRead>,
    parser: DatasetSmilesParser,
    line_number: usize,
    line_buffer: String,
}

/// The official PubChem `CID-SMILES.gz` bulk download.
///
/// Source:
/// `https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz`
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq)]
pub struct PubChemSmiles;

impl DatasetSource for PubChemSmiles {
    fn id(&self) -> &'static str {
        "pubchem-smiles"
    }

    fn url(&self) -> &'static str {
        "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
    }

    fn file_name(&self) -> &'static str {
        "CID-SMILES.gz"
    }

    fn extracted_file_name(&self) -> &'static str {
        "CID-SMILES"
    }

    fn compression(&self) -> DatasetCompression {
        DatasetCompression::Gzip
    }
}

impl SmilesDatasetSource for PubChemSmiles {
    fn iter_smiles_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesIter, DatasetError> {
        let artifact = self.fetch_with_options(options)?;
        DatasetSmilesIter::for_pubchem(&artifact)
    }
}

/// The official MassSpecGym benchmark TSV containing a `smiles` column.
///
/// Source:
/// `https://huggingface.co/datasets/roman-bushuiev/MassSpecGym/resolve/main/data/MassSpecGym.tsv`
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq)]
pub struct MassSpecGymSmiles;

impl DatasetSource for MassSpecGymSmiles {
    fn id(&self) -> &'static str {
        "massspecgym-smiles"
    }

    fn url(&self) -> &'static str {
        "https://huggingface.co/datasets/roman-bushuiev/MassSpecGym/resolve/main/data/MassSpecGym.tsv"
    }

    fn file_name(&self) -> &'static str {
        "MassSpecGym.tsv"
    }
}

impl SmilesDatasetSource for MassSpecGymSmiles {
    fn iter_smiles_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesIter, DatasetError> {
        let artifact = self.fetch_with_options(options)?;
        DatasetSmilesIter::for_mass_spec_gym(&artifact)
    }
}

/// Convenient constant handle for the PubChem SMILES dataset.
pub const PUBCHEM_SMILES: PubChemSmiles = PubChemSmiles;

/// Convenient constant handle for the MassSpecGym SMILES dataset.
pub const MASS_SPEC_GYM_SMILES: MassSpecGymSmiles = MassSpecGymSmiles;

/// Returns the default cache directory used by dataset fetches.
///
/// The selection order is:
///
/// 1. `XDG_CACHE_HOME/smiles-parser/datasets`
/// 2. `LOCALAPPDATA/smiles-parser/datasets`
/// 3. `HOME/.cache/smiles-parser/datasets`
/// 4. `${TMPDIR}/smiles-parser/datasets`
#[must_use]
pub fn default_dataset_cache_dir() -> PathBuf {
    if let Some(path) = env::var_os("XDG_CACHE_HOME") {
        return PathBuf::from(path).join("smiles-parser").join("datasets");
    }
    if let Some(path) = env::var_os("LOCALAPPDATA") {
        return PathBuf::from(path).join("smiles-parser").join("datasets");
    }
    if let Some(path) = env::var_os("HOME") {
        return PathBuf::from(path).join(".cache").join("smiles-parser").join("datasets");
    }
    env::temp_dir().join("smiles-parser").join("datasets")
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum DatasetSmilesParser {
    PubChem,
    MassSpecGym { smiles_column: usize },
}

impl DatasetSmilesIter {
    fn for_pubchem(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        Self::from_artifact(artifact, DatasetSmilesParser::PubChem)
    }

    fn for_mass_spec_gym(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        let dataset_id = artifact.dataset_id;
        let path = artifact.path.clone();
        let mut reader = open_text_reader(&path)?;
        let mut header = String::new();
        let bytes_read = reader
            .read_line(&mut header)
            .map_err(|source| DatasetError::Io { path: path.clone(), source })?;
        if bytes_read == 0 {
            return Err(DatasetError::Format {
                dataset_id,
                line_number: 1,
                message: "expected a TSV header row with a smiles column".into(),
            });
        }

        let smiles_column = header
            .trim_end_matches(['\r', '\n'])
            .split('\t')
            .position(|field| field.eq_ignore_ascii_case("smiles"))
            .ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number: 1,
                    message: "expected a TSV header containing a smiles column".into(),
                }
            })?;

        Ok(Self {
            dataset_id,
            path,
            reader,
            parser: DatasetSmilesParser::MassSpecGym { smiles_column },
            line_number: 1,
            line_buffer: String::new(),
        })
    }

    fn from_artifact(
        artifact: &DatasetArtifact,
        parser: DatasetSmilesParser,
    ) -> Result<Self, DatasetError> {
        let path = artifact.path.clone();
        Ok(Self {
            dataset_id: artifact.dataset_id,
            path: path.clone(),
            reader: open_text_reader(&path)?,
            parser,
            line_number: 0,
            line_buffer: String::new(),
        })
    }
}

impl Iterator for DatasetSmilesIter {
    type Item = Result<String, DatasetError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line_number += 1;
                }
                Err(source) => {
                    return Some(Err(DatasetError::Io { path: self.path.clone(), source }));
                }
            }

            let line = self.line_buffer.trim_end_matches(['\r', '\n']);
            if line.is_empty() {
                continue;
            }

            return Some(parse_smiles_line(self.dataset_id, self.line_number, self.parser, line));
        }
    }
}

fn parse_smiles_line(
    dataset_id: &'static str,
    line_number: usize,
    parser: DatasetSmilesParser,
    line: &str,
) -> Result<String, DatasetError> {
    match parser {
        DatasetSmilesParser::PubChem => {
            let (_, smiles) = line.split_once('\t').ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a CID<TAB>SMILES record".into(),
                }
            })?;
            Ok(smiles.to_owned())
        }
        DatasetSmilesParser::MassSpecGym { smiles_column } => {
            let smiles = tsv_field(line, smiles_column).ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a TSV row with a smiles column value".into(),
                }
            })?;
            Ok(smiles.to_owned())
        }
    }
}

fn tsv_field(line: &str, column_index: usize) -> Option<&str> {
    line.split('\t').nth(column_index)
}

fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>, DatasetError> {
    let file =
        File::open(path).map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })?;

    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn fetch_dataset<D: DatasetSource + ?Sized>(
    dataset: &D,
    options: &DatasetFetchOptions,
) -> Result<DatasetArtifact, DatasetError> {
    let cache_root = options.cache_dir.clone().unwrap_or_else(default_dataset_cache_dir);
    let dataset_dir = cache_root.join(dataset.id());
    create_dir_all(&dataset_dir)?;

    let compressed_path = dataset_dir.join(dataset.file_name());
    let decompressed_path = dataset_dir.join(dataset.extracted_file_name());

    match dataset.compression() {
        DatasetCompression::None => {
            let was_downloaded = ensure_downloaded(dataset, &compressed_path, options.cache_mode)?;
            Ok(DatasetArtifact {
                dataset_id: dataset.id(),
                path: compressed_path.clone(),
                compressed_path: Some(compressed_path),
                decompressed_path: None,
                was_downloaded,
                was_decompressed: false,
            })
        }
        DatasetCompression::Gzip => {
            match options.gzip_mode {
                GzipMode::KeepCompressed => {
                    let was_downloaded =
                        ensure_downloaded(dataset, &compressed_path, options.cache_mode)?;
                    Ok(DatasetArtifact {
                        dataset_id: dataset.id(),
                        path: compressed_path.clone(),
                        compressed_path: Some(compressed_path),
                        decompressed_path: None,
                        was_downloaded,
                        was_decompressed: false,
                    })
                }
                GzipMode::Decompress | GzipMode::KeepBoth => {
                    let (was_downloaded, was_decompressed) = ensure_decompressed(
                        dataset,
                        &compressed_path,
                        &decompressed_path,
                        options.cache_mode,
                    )?;
                    Ok(DatasetArtifact {
                        dataset_id: dataset.id(),
                        path: decompressed_path.clone(),
                        compressed_path: compressed_path.is_file().then_some(compressed_path),
                        decompressed_path: Some(decompressed_path),
                        was_downloaded,
                        was_decompressed,
                    })
                }
            }
        }
    }
}

fn ensure_downloaded<D: DatasetSource + ?Sized>(
    dataset: &D,
    target_path: &Path,
    cache_mode: CacheMode,
) -> Result<bool, DatasetError> {
    if cache_mode == CacheMode::UseCache && target_path.is_file() {
        return Ok(false);
    }

    download_to_path(dataset.url(), target_path)?;
    Ok(true)
}

fn ensure_decompressed<D: DatasetSource + ?Sized>(
    dataset: &D,
    compressed_path: &Path,
    decompressed_path: &Path,
    cache_mode: CacheMode,
) -> Result<(bool, bool), DatasetError> {
    if cache_mode == CacheMode::UseCache && decompressed_path.is_file() {
        return Ok((false, false));
    }

    let was_downloaded = ensure_downloaded(dataset, compressed_path, cache_mode)?;
    let was_decompressed = gunzip_file(compressed_path, decompressed_path)?;
    Ok((was_downloaded, was_decompressed))
}

fn download_to_path(url: &'static str, target_path: &Path) -> Result<(), DatasetError> {
    let client = Client::builder()
        .user_agent(DOWNLOAD_USER_AGENT)
        .build()
        .map_err(|source| DatasetError::Download { url, source })?;
    let response = client
        .get(url)
        .send()
        .and_then(reqwest::blocking::Response::error_for_status)
        .map_err(|source| DatasetError::Download { url, source })?;
    let file_label = progress_label("downloading", target_path);
    let progress_bar = new_byte_progress_bar(response.content_length(), &file_label);
    let mut response = ProgressReader::new(response, progress_bar.clone());

    let temporary_path = temporary_download_path(target_path);
    write_parent_dir(target_path)?;
    let file = File::create(&temporary_path)
        .map_err(|source| DatasetError::Io { path: temporary_path.clone(), source })?;
    let mut writer = BufWriter::new(file);
    if let Err(source) = io::copy(&mut response, &mut writer) {
        progress_bar.abandon();
        return Err(DatasetError::Io { path: target_path.to_path_buf(), source });
    }
    writer
        .flush()
        .map_err(|source| DatasetError::Io { path: target_path.to_path_buf(), source })?;
    progress_bar.finish_and_clear();

    if target_path.exists() {
        fs::remove_file(target_path)
            .map_err(|source| DatasetError::Io { path: target_path.to_path_buf(), source })?;
    }
    fs::rename(&temporary_path, target_path)
        .map_err(|source| DatasetError::Io { path: target_path.to_path_buf(), source })?;
    Ok(())
}

fn gunzip_file(compressed_path: &Path, decompressed_path: &Path) -> Result<bool, DatasetError> {
    write_parent_dir(decompressed_path)?;
    let source_file = File::open(compressed_path)
        .map_err(|source| DatasetError::Io { path: compressed_path.to_path_buf(), source })?;
    let progress_bar = new_byte_progress_bar(
        source_file.metadata().ok().map(|metadata| metadata.len()),
        &progress_label("decompressing", decompressed_path),
    );
    let source_file = ProgressReader::new(source_file, progress_bar.clone());
    let mut decoder = GzDecoder::new(source_file);
    let temporary_path = temporary_download_path(decompressed_path);
    let target_file = File::create(&temporary_path)
        .map_err(|source| DatasetError::Io { path: temporary_path.clone(), source })?;
    let mut writer = BufWriter::new(target_file);
    if let Err(source) = io::copy(&mut decoder, &mut writer) {
        progress_bar.abandon();
        return Err(DatasetError::Io { path: decompressed_path.to_path_buf(), source });
    }
    writer
        .flush()
        .map_err(|source| DatasetError::Io { path: decompressed_path.to_path_buf(), source })?;
    progress_bar.finish_and_clear();

    if decompressed_path.exists() {
        fs::remove_file(decompressed_path)
            .map_err(|source| DatasetError::Io { path: decompressed_path.to_path_buf(), source })?;
    }
    fs::rename(&temporary_path, decompressed_path)
        .map_err(|source| DatasetError::Io { path: decompressed_path.to_path_buf(), source })?;
    Ok(true)
}

fn temporary_download_path(target_path: &Path) -> PathBuf {
    let file_name = target_path
        .file_name()
        .map_or_else(|| "download".into(), |name| name.to_string_lossy().into_owned());
    target_path.with_file_name(format!("{file_name}.part"))
}

fn write_parent_dir(path: &Path) -> Result<(), DatasetError> {
    if let Some(parent) = path.parent() {
        create_dir_all(parent)?;
    }
    Ok(())
}

fn create_dir_all(path: &Path) -> Result<(), DatasetError> {
    fs::create_dir_all(path).map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })
}

fn progress_label(action: &str, path: &Path) -> String {
    let file_name = path
        .file_name()
        .map_or_else(|| "dataset".into(), |name| name.to_string_lossy().into_owned());
    format!("{action} {file_name}")
}

fn new_byte_progress_bar(total_bytes: Option<u64>, message: &str) -> ProgressBar {
    if !io::stderr().is_terminal() {
        return ProgressBar::hidden();
    }

    if let Some(total_bytes) = total_bytes.filter(|&bytes| bytes > 0) {
        let progress_bar = ProgressBar::new(total_bytes);
        progress_bar.set_style(
            ProgressStyle::with_template(BYTE_PROGRESS_TEMPLATE)
                .unwrap_or_else(|_| unreachable!("progress template is static and valid")),
        );
        progress_bar.set_message(message.to_owned());
        progress_bar
    } else {
        let progress_bar = ProgressBar::new_spinner();
        progress_bar.set_style(
            ProgressStyle::with_template(SPINNER_PROGRESS_TEMPLATE)
                .unwrap_or_else(|_| unreachable!("progress template is static and valid")),
        );
        progress_bar.set_message(message.to_owned());
        progress_bar.enable_steady_tick(PROGRESS_TICK_INTERVAL);
        progress_bar
    }
}

struct ProgressReader<R> {
    inner: R,
    progress_bar: ProgressBar,
}

impl<R> ProgressReader<R> {
    fn new(inner: R, progress_bar: ProgressBar) -> Self {
        Self { inner, progress_bar }
    }
}

impl<R: Read> Read for ProgressReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, io::Error> {
        let bytes_read = self.inner.read(buf)?;
        if bytes_read > 0 {
            self.progress_bar.inc(u64::try_from(bytes_read).unwrap_or(u64::MAX));
        }
        Ok(bytes_read)
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs::{self, File},
        io::Write,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
        vec::Vec,
    };

    use flate2::{Compression, write::GzEncoder};

    use super::{
        CacheMode, DatasetArtifact, DatasetCompression, DatasetFetchOptions, DatasetSmilesIter,
        DatasetSource, GzipMode, MASS_SPEC_GYM_SMILES, PUBCHEM_SMILES, PubChemSmiles,
        default_dataset_cache_dir, gunzip_file,
    };

    fn temporary_directory(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_else(|_| unreachable!("system time is after unix epoch"))
            .as_nanos();
        std::env::temp_dir().join(format!("smiles-parser-{name}-{unique}"))
    }

    #[test]
    fn pubchem_smiles_metadata_matches_current_upstream_layout() {
        let dataset = PubChemSmiles;

        assert_eq!(dataset.id(), "pubchem-smiles");
        assert_eq!(dataset.file_name(), "CID-SMILES.gz");
        assert_eq!(dataset.extracted_file_name(), "CID-SMILES");
        assert_eq!(dataset.compression(), DatasetCompression::Gzip);
        assert!(dataset.url().contains("pubchem/Compound/Extras/CID-SMILES.gz"));
    }

    #[test]
    fn massspecgym_smiles_metadata_matches_current_upstream_layout() {
        assert_eq!(MASS_SPEC_GYM_SMILES.id(), "massspecgym-smiles");
        assert_eq!(MASS_SPEC_GYM_SMILES.file_name(), "MassSpecGym.tsv");
        assert_eq!(MASS_SPEC_GYM_SMILES.extracted_file_name(), "MassSpecGym.tsv");
        assert_eq!(MASS_SPEC_GYM_SMILES.compression(), DatasetCompression::None);
        assert!(MASS_SPEC_GYM_SMILES.url().contains("/MassSpecGym.tsv"));
    }

    #[test]
    fn default_fetch_options_keep_compressed_cache_behavior() {
        let options = DatasetFetchOptions::default();

        assert_eq!(options.cache_mode, CacheMode::UseCache);
        assert_eq!(options.gzip_mode, GzipMode::KeepCompressed);
        assert!(options.cache_dir.is_none());
    }

    #[test]
    fn default_dataset_cache_dir_has_stable_suffix() {
        let cache_dir = default_dataset_cache_dir();

        assert!(cache_dir.ends_with(PathBuf::from("smiles-parser").join("datasets")));
    }

    #[test]
    fn gunzip_file_materializes_plaintext_copy() {
        let directory = temporary_directory("datasets-gunzip");
        fs::create_dir_all(&directory).unwrap();
        let compressed_path = directory.join("sample.txt.gz");
        let decompressed_path = directory.join("sample.txt");

        {
            let file = File::create(&compressed_path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(b"cid\tsmiles\n1\tCCO\n").unwrap();
            encoder.finish().unwrap();
        }

        gunzip_file(&compressed_path, &decompressed_path).unwrap();
        assert_eq!(fs::read(&decompressed_path).unwrap(), b"cid\tsmiles\n1\tCCO\n");

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn pubchem_and_massspecgym_constants_are_usable_dataset_handles() {
        assert_eq!(PUBCHEM_SMILES.id(), "pubchem-smiles");
        assert_eq!(MASS_SPEC_GYM_SMILES.id(), "massspecgym-smiles");
    }

    #[test]
    fn pubchem_smiles_iterator_streams_smiles_from_gzip_records() {
        let directory = temporary_directory("datasets-pubchem-iter");
        fs::create_dir_all(&directory).unwrap();
        let compressed_path = directory.join("CID-SMILES.gz");

        {
            let file = File::create(&compressed_path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(b"1\tCCO\n2\tc1ccccc1\n").unwrap();
            encoder.finish().unwrap();
        }

        let artifact = DatasetArtifact {
            dataset_id: "pubchem-smiles",
            path: compressed_path.clone(),
            compressed_path: Some(compressed_path),
            decompressed_path: None,
            was_downloaded: false,
            was_decompressed: false,
        };
        let smiles = DatasetSmilesIter::for_pubchem(&artifact)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(smiles, ["CCO", "c1ccccc1"]);

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn massspecgym_smiles_iterator_uses_smiles_tsv_column() {
        let directory = temporary_directory("datasets-massspecgym-iter");
        fs::create_dir_all(&directory).unwrap();
        let dataset_path = directory.join("MassSpecGym.tsv");
        fs::write(&dataset_path, "spec_id\tname\tsmiles\n1\tethanol\tCCO\n2\tbenzene\tc1ccccc1\n")
            .unwrap();

        let artifact = DatasetArtifact {
            dataset_id: "massspecgym-smiles",
            path: dataset_path,
            compressed_path: None,
            decompressed_path: None,
            was_downloaded: false,
            was_decompressed: false,
        };
        let smiles = DatasetSmilesIter::for_mass_spec_gym(&artifact)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(smiles, ["CCO", "c1ccccc1"]);

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn massspecgym_smiles_iterator_requires_smiles_header_column() {
        let directory = temporary_directory("datasets-massspecgym-header");
        fs::create_dir_all(&directory).unwrap();
        let dataset_path = directory.join("MassSpecGym.tsv");
        fs::write(&dataset_path, "spec_id\tname\n1\tethanol\n").unwrap();

        let artifact = DatasetArtifact {
            dataset_id: "massspecgym-smiles",
            path: dataset_path,
            compressed_path: None,
            decompressed_path: None,
            was_downloaded: false,
            was_decompressed: false,
        };
        match DatasetSmilesIter::for_mass_spec_gym(&artifact) {
            Ok(_) => panic!("expected a missing smiles header to fail"),
            Err(super::DatasetError::Format {
                dataset_id: "massspecgym-smiles",
                line_number: 1,
                ..
            }) => {}
            Err(error) => panic!("unexpected error: {error}"),
        }

        fs::remove_dir_all(directory).unwrap();
    }
}
