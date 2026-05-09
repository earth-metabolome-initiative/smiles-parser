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

use alloc::{borrow::ToOwned, boxed::Box, string::String, vec::Vec};
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
use tar::Archive;
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
    /// The upstream file is a gzip-compressed tar archive.
    TarGzip,
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
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetSource, MASS_SPEC_GYM_SMILES};
    ///
    /// let artifact = MASS_SPEC_GYM_SMILES.fetch()?;
    /// assert_eq!(artifact.dataset_id(), "massspecgym-smiles");
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn dataset_id(&self) -> &'static str {
        self.dataset_id
    }

    /// Returns the primary path callers should consume.
    ///
    /// This is the decompressed path when one was requested, otherwise the
    /// cached upstream file path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetSource, MASS_SPEC_GYM_SMILES};
    ///
    /// let artifact = MASS_SPEC_GYM_SMILES.fetch()?;
    /// println!("{}", artifact.path().display());
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Returns the cached compressed path, when present.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetFetchOptions, DatasetSource, GzipMode, PUBCHEM_SMILES};
    ///
    /// let artifact = PUBCHEM_SMILES.fetch_with_options(&DatasetFetchOptions {
    ///     gzip_mode: GzipMode::KeepBoth,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// assert!(artifact.compressed_path().is_some());
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn compressed_path(&self) -> Option<&Path> {
        self.compressed_path.as_deref()
    }

    /// Returns the cached decompressed path, when present.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetFetchOptions, DatasetSource, GzipMode, PUBCHEM_SMILES};
    ///
    /// let artifact = PUBCHEM_SMILES.fetch_with_options(&DatasetFetchOptions {
    ///     gzip_mode: GzipMode::Decompress,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// assert!(artifact.decompressed_path().is_some());
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn decompressed_path(&self) -> Option<&Path> {
        self.decompressed_path.as_deref()
    }

    /// Returns whether the upstream artifact was downloaded during this call.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{
    ///     CacheMode, DatasetFetchOptions, DatasetSource, MASS_SPEC_GYM_SMILES,
    /// };
    ///
    /// let artifact = MASS_SPEC_GYM_SMILES.fetch_with_options(&DatasetFetchOptions {
    ///     cache_mode: CacheMode::UseCache,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// let _downloaded = artifact.was_downloaded();
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn was_downloaded(&self) -> bool {
        self.was_downloaded
    }

    /// Returns whether a decompressed copy was created during this call.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetFetchOptions, DatasetSource, GzipMode, PUBCHEM_SMILES};
    ///
    /// let artifact = PUBCHEM_SMILES.fetch_with_options(&DatasetFetchOptions {
    ///     gzip_mode: GzipMode::Decompress,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// let _decompressed = artifact.was_decompressed();
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    #[must_use]
    pub fn was_decompressed(&self) -> bool {
        self.was_decompressed
    }
}

/// A materialized multi-file dataset collection on disk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetCollectionArtifact {
    dataset_id: &'static str,
    paths: Vec<PathBuf>,
    compressed_paths: Vec<PathBuf>,
    was_downloaded: bool,
    was_extracted: bool,
}

impl DatasetCollectionArtifact {
    /// Returns the dataset identifier that produced this artifact collection.
    #[must_use]
    pub fn dataset_id(&self) -> &'static str {
        self.dataset_id
    }

    /// Returns the primary paths callers should consume.
    ///
    /// For archive-based datasets these are extracted directories when
    /// extraction was requested, otherwise the cached archive paths.
    #[must_use]
    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    /// Returns the cached compressed archive paths.
    #[must_use]
    pub fn compressed_paths(&self) -> &[PathBuf] {
        &self.compressed_paths
    }

    /// Returns whether any upstream artifact was downloaded during this call.
    #[must_use]
    pub fn was_downloaded(&self) -> bool {
        self.was_downloaded
    }

    /// Returns whether any upstream archive was extracted during this call.
    #[must_use]
    pub fn was_extracted(&self) -> bool {
        self.was_extracted
    }
}

/// Metadata for one file within a downloadable dataset collection.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct DatasetFile {
    url: &'static str,
    file_name: &'static str,
    extracted_file_name: &'static str,
    compression: DatasetCompression,
}

impl DatasetFile {
    /// Creates static metadata for one downloadable dataset file.
    #[must_use]
    pub const fn new(
        url: &'static str,
        file_name: &'static str,
        extracted_file_name: &'static str,
        compression: DatasetCompression,
    ) -> Self {
        Self { url, file_name, extracted_file_name, compression }
    }

    /// Returns the upstream URL for this dataset file.
    #[must_use]
    pub fn url(&self) -> &'static str {
        self.url
    }

    /// Returns the cached upstream file name.
    #[must_use]
    pub fn file_name(&self) -> &'static str {
        self.file_name
    }

    /// Returns the extracted file or directory name.
    #[must_use]
    pub fn extracted_file_name(&self) -> &'static str {
        self.extracted_file_name
    }

    /// Returns the compression used by this dataset file.
    #[must_use]
    pub fn compression(&self) -> DatasetCompression {
        self.compression
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
    /// The requested dataset subset is not valid.
    #[error("invalid dataset selection for {dataset_id}: {message}")]
    InvalidSelection {
        /// The dataset identifier.
        dataset_id: &'static str,
        /// A human-readable explanation of the invalid selection.
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
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetSource, MASS_SPEC_GYM_SMILES};
    ///
    /// let artifact = MASS_SPEC_GYM_SMILES.fetch()?;
    /// assert!(artifact.path().ends_with("MassSpecGym.tsv"));
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    fn fetch(&self) -> Result<DatasetArtifact, DatasetError> {
        self.fetch_with_options(&DatasetFetchOptions::default())
    }

    /// Fetches the dataset using explicit options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the download fails or the requested cached
    /// artifact layout cannot be materialized on disk.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{DatasetFetchOptions, DatasetSource, GzipMode, PUBCHEM_SMILES};
    ///
    /// let artifact = PUBCHEM_SMILES.fetch_with_options(&DatasetFetchOptions {
    ///     gzip_mode: GzipMode::Decompress,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// assert!(artifact.path().ends_with("CID-SMILES"));
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    fn fetch_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetArtifact, DatasetError> {
        fetch_dataset(self, options)
    }
}

/// Metadata for a dataset that is distributed across multiple files.
pub trait DatasetCollectionSource {
    /// Stable dataset identifier used for cache subdirectories and diagnostics.
    fn id(&self) -> &'static str;

    /// Upstream files that make up this dataset collection.
    fn files(&self) -> Vec<DatasetFile>;

    /// Fetches the dataset collection using default options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if any file cannot be downloaded or
    /// materialized on disk.
    fn fetch_collection(&self) -> Result<DatasetCollectionArtifact, DatasetError> {
        self.fetch_collection_with_options(&DatasetFetchOptions::default())
    }

    /// Fetches the dataset collection using explicit options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if any file cannot be downloaded or
    /// materialized on disk.
    fn fetch_collection_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetCollectionArtifact, DatasetError> {
        fetch_dataset_collection(self, options)
    }
}

/// A dataset source that can stream SMILES strings directly.
pub trait SmilesDatasetSource {
    /// Opens a streaming iterator over the dataset SMILES using default fetch
    /// options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the dataset cannot be fetched or if the
    /// materialized dataset cannot be opened for streaming.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{PUBCHEM_SMILES, SmilesDatasetSource};
    ///
    /// let mut smiles = PUBCHEM_SMILES.iter_smiles()?;
    /// if let Some(first) = smiles.next() {
    ///     let _ = first?;
    /// }
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
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
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use smiles_parser::datasets::{
    ///     DatasetFetchOptions, GzipMode, PUBCHEM_SMILES, SmilesDatasetSource,
    /// };
    ///
    /// let mut smiles = PUBCHEM_SMILES.iter_smiles_with_options(&DatasetFetchOptions {
    ///     gzip_mode: GzipMode::KeepCompressed,
    ///     ..DatasetFetchOptions::default()
    /// })?;
    /// if let Some(first) = smiles.next() {
    ///     let _ = first?;
    /// }
    /// # Ok::<(), smiles_parser::DatasetError>(())
    /// ```
    fn iter_smiles_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesIter, DatasetError>;
}

/// A dataset source that can stream SMILES records with dataset identifiers.
pub trait SmilesDatasetRecordSource {
    /// Opens a streaming iterator over dataset records using default fetch
    /// options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the dataset cannot be fetched or opened.
    fn iter_records(&self) -> Result<DatasetSmilesRecordIter, DatasetError> {
        self.iter_records_with_options(&DatasetFetchOptions::default())
    }

    /// Opens a streaming iterator over dataset records using explicit fetch
    /// options.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError`] if the dataset cannot be fetched or opened.
    fn iter_records_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesRecordIter, DatasetError>;
}

/// One SMILES record from a dataset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetSmilesRecord {
    id: String,
    smiles: String,
}

impl DatasetSmilesRecord {
    /// Creates a dataset SMILES record.
    #[must_use]
    pub fn new(id: String, smiles: String) -> Self {
        Self { id, smiles }
    }

    /// Returns the dataset-specific record identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the SMILES string.
    #[must_use]
    pub fn smiles(&self) -> &str {
        &self.smiles
    }

    /// Consumes the record and returns its SMILES string.
    #[must_use]
    pub fn into_smiles(self) -> String {
        self.smiles
    }
}

/// A streaming iterator over SMILES strings extracted from a dataset file.
pub struct DatasetSmilesIter {
    inner: DatasetSmilesRecordIter,
}

/// A streaming iterator over SMILES records extracted from a dataset file.
pub struct DatasetSmilesRecordIter {
    dataset_id: &'static str,
    paths: Vec<PathBuf>,
    next_path_index: usize,
    current: Option<DatasetReader>,
    parser: DatasetSmilesParser,
    line_number: usize,
    line_buffer: String,
}

struct DatasetReader {
    path: PathBuf,
    reader: Box<dyn BufRead + Send>,
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
        Ok(DatasetSmilesIter::from_records(self.iter_records_with_options(options)?))
    }
}

impl SmilesDatasetRecordSource for PubChemSmiles {
    fn iter_records_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesRecordIter, DatasetError> {
        let artifact = self.fetch_with_options(options)?;
        DatasetSmilesRecordIter::for_pubchem(&artifact)
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
        Ok(DatasetSmilesIter::from_records(self.iter_records_with_options(options)?))
    }
}

impl SmilesDatasetRecordSource for MassSpecGymSmiles {
    fn iter_records_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesRecordIter, DatasetError> {
        let artifact = self.fetch_with_options(options)?;
        DatasetSmilesRecordIter::for_mass_spec_gym(&artifact)
    }
}

/// Number of rows reported by the ZINC20-ML `smiles_count.txt` manifest.
pub const ZINC20_EXPECTED_RECORD_COUNT: usize = 1_006_651_037;

macro_rules! zinc20_chunk_file {
    ($chunk:literal) => {
        DatasetFile::new(
            concat!(
                "https://files.docking.org/zinc20-ML/smiles/ZINC20_smiles_chunk_",
                $chunk,
                ".tar.gz"
            ),
            concat!("ZINC20_smiles_chunk_", $chunk, ".tar.gz"),
            concat!("ZINC20_smiles_chunk_", $chunk),
            DatasetCompression::TarGzip,
        )
    };
}

const ZINC20_CHUNK_FILES: [DatasetFile; 20] = [
    zinc20_chunk_file!(1),
    zinc20_chunk_file!(2),
    zinc20_chunk_file!(3),
    zinc20_chunk_file!(4),
    zinc20_chunk_file!(5),
    zinc20_chunk_file!(6),
    zinc20_chunk_file!(7),
    zinc20_chunk_file!(8),
    zinc20_chunk_file!(9),
    zinc20_chunk_file!(10),
    zinc20_chunk_file!(11),
    zinc20_chunk_file!(12),
    zinc20_chunk_file!(13),
    zinc20_chunk_file!(14),
    zinc20_chunk_file!(15),
    zinc20_chunk_file!(16),
    zinc20_chunk_file!(17),
    zinc20_chunk_file!(18),
    zinc20_chunk_file!(19),
    zinc20_chunk_file!(20),
];

/// The ZINC20-ML SMILES dataset distributed as 20 `tar.gz` chunks.
///
/// Source:
/// `https://files.docking.org/zinc20-ML/smiles/`
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Zinc20Smiles {
    first_chunk: u8,
    last_chunk: u8,
}

impl Default for Zinc20Smiles {
    fn default() -> Self {
        Self::all()
    }
}

impl Zinc20Smiles {
    /// First available ZINC20-ML SMILES chunk.
    pub const FIRST_CHUNK: u8 = 1;
    /// Last available ZINC20-ML SMILES chunk.
    pub const LAST_CHUNK: u8 = 20;

    /// Returns a handle spanning all ZINC20-ML SMILES chunks.
    #[must_use]
    pub const fn all() -> Self {
        Self { first_chunk: Self::FIRST_CHUNK, last_chunk: Self::LAST_CHUNK }
    }

    /// Returns a handle for one ZINC20-ML SMILES chunk.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError::InvalidSelection`] if `chunk` is outside
    /// `1..=20`.
    pub fn chunk(chunk: u8) -> Result<Self, DatasetError> {
        Self::chunk_range(chunk, chunk)
    }

    /// Returns a handle for an inclusive ZINC20-ML chunk range.
    ///
    /// # Errors
    ///
    /// Returns [`DatasetError::InvalidSelection`] if the range is empty or
    /// outside `1..=20`.
    pub fn chunk_range(first_chunk: u8, last_chunk: u8) -> Result<Self, DatasetError> {
        if first_chunk < Self::FIRST_CHUNK
            || last_chunk > Self::LAST_CHUNK
            || first_chunk > last_chunk
        {
            return Err(DatasetError::InvalidSelection {
                dataset_id: "zinc20-smiles",
                message: format!(
                    "expected an inclusive chunk range within {}..={}, got {first_chunk}..={last_chunk}",
                    Self::FIRST_CHUNK,
                    Self::LAST_CHUNK
                ),
            });
        }

        Ok(Self { first_chunk, last_chunk })
    }

    /// Returns the first selected chunk.
    #[must_use]
    pub fn first_chunk(&self) -> u8 {
        self.first_chunk
    }

    /// Returns the last selected chunk.
    #[must_use]
    pub fn last_chunk(&self) -> u8 {
        self.last_chunk
    }
}

impl DatasetCollectionSource for Zinc20Smiles {
    fn id(&self) -> &'static str {
        "zinc20-smiles"
    }

    fn files(&self) -> Vec<DatasetFile> {
        let start = usize::from(self.first_chunk - 1);
        let end = usize::from(self.last_chunk);
        ZINC20_CHUNK_FILES[start..end].to_vec()
    }
}

impl SmilesDatasetSource for Zinc20Smiles {
    fn iter_smiles_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesIter, DatasetError> {
        Ok(DatasetSmilesIter::from_records(self.iter_records_with_options(options)?))
    }
}

impl SmilesDatasetRecordSource for Zinc20Smiles {
    fn iter_records_with_options(
        &self,
        options: &DatasetFetchOptions,
    ) -> Result<DatasetSmilesRecordIter, DatasetError> {
        let mut options = options.clone();
        if options.gzip_mode == GzipMode::KeepCompressed {
            options.gzip_mode = GzipMode::Decompress;
        }
        let artifact = self.fetch_collection_with_options(&options)?;
        DatasetSmilesRecordIter::for_zinc20(&artifact)
    }
}

/// Convenient constant handle for the PubChem SMILES dataset.
pub const PUBCHEM_SMILES: PubChemSmiles = PubChemSmiles;

/// Convenient constant handle for the MassSpecGym SMILES dataset.
pub const MASS_SPEC_GYM_SMILES: MassSpecGymSmiles = MassSpecGymSmiles;

/// Convenient constant handle for the full ZINC20-ML SMILES dataset.
pub const ZINC20_SMILES: Zinc20Smiles = Zinc20Smiles::all();

/// Returns the default cache directory used by dataset fetches.
///
/// The selection order is:
///
/// 1. `XDG_CACHE_HOME/smiles-parser/datasets`
/// 2. `LOCALAPPDATA/smiles-parser/datasets`
/// 3. `HOME/.cache/smiles-parser/datasets`
/// 4. `${TMPDIR}/smiles-parser/datasets`
///
/// # Examples
///
/// ```
/// use smiles_parser::datasets::default_dataset_cache_dir;
///
/// assert!(default_dataset_cache_dir().ends_with("smiles-parser/datasets"));
/// ```
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
    Zinc20,
}

impl DatasetSmilesIter {
    #[cfg(test)]
    fn for_pubchem(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        Ok(Self::from_records(DatasetSmilesRecordIter::for_pubchem(artifact)?))
    }

    #[cfg(test)]
    fn for_mass_spec_gym(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        Ok(Self::from_records(DatasetSmilesRecordIter::for_mass_spec_gym(artifact)?))
    }

    fn from_records(inner: DatasetSmilesRecordIter) -> Self {
        Self { inner }
    }
}

impl Iterator for DatasetSmilesIter {
    type Item = Result<String, DatasetError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|record| record.map(DatasetSmilesRecord::into_smiles))
    }
}

impl DatasetSmilesRecordIter {
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
            paths: Vec::new(),
            next_path_index: 0,
            current: Some(DatasetReader { path, reader }),
            parser: DatasetSmilesParser::MassSpecGym { smiles_column },
            line_number: 1,
            line_buffer: String::new(),
        })
    }

    fn for_zinc20(artifact: &DatasetCollectionArtifact) -> Result<Self, DatasetError> {
        let mut paths = Vec::new();
        for path in artifact.paths() {
            collect_zinc20_smiles_paths(artifact.dataset_id, path, &mut paths)?;
        }
        if paths.is_empty() {
            return Err(DatasetError::Format {
                dataset_id: artifact.dataset_id,
                line_number: 1,
                message: "expected at least one extracted ZINC20 smiles_all_*.txt file".into(),
            });
        }

        Ok(Self {
            dataset_id: artifact.dataset_id,
            paths,
            next_path_index: 0,
            current: None,
            parser: DatasetSmilesParser::Zinc20,
            line_number: 0,
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
            paths: Vec::new(),
            next_path_index: 0,
            current: Some(DatasetReader { path: path.clone(), reader: open_text_reader(&path)? }),
            parser,
            line_number: 0,
            line_buffer: String::new(),
        })
    }

    fn open_next_reader(&mut self) -> Option<Result<(), DatasetError>> {
        let path = self.paths.get(self.next_path_index)?.clone();
        self.next_path_index += 1;
        match open_text_reader(&path) {
            Ok(reader) => {
                self.current = Some(DatasetReader { path, reader });
                self.line_number = 0;
                Some(Ok(()))
            }
            Err(error) => Some(Err(error)),
        }
    }
}

impl Iterator for DatasetSmilesRecordIter {
    type Item = Result<DatasetSmilesRecord, DatasetError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.current.is_none() {
                match self.open_next_reader()? {
                    Ok(()) => {}
                    Err(error) => return Some(Err(error)),
                }
            }

            self.line_buffer.clear();
            let read_result = {
                let current =
                    self.current.as_mut().unwrap_or_else(|| unreachable!("current reader is open"));
                current.reader.read_line(&mut self.line_buffer)
            };
            match read_result {
                Ok(0) => {
                    self.current = None;
                    continue;
                }
                Ok(_) => {
                    self.line_number += 1;
                }
                Err(source) => {
                    let path = self
                        .current
                        .as_ref()
                        .map_or_else(PathBuf::new, |current| current.path.clone());
                    return Some(Err(DatasetError::Io { path, source }));
                }
            }

            let line = self.line_buffer.trim_end_matches(['\r', '\n']);
            if line.is_empty() {
                continue;
            }

            return Some(parse_smiles_record(self.dataset_id, self.line_number, self.parser, line));
        }
    }
}

fn parse_smiles_record(
    dataset_id: &'static str,
    line_number: usize,
    parser: DatasetSmilesParser,
    line: &str,
) -> Result<DatasetSmilesRecord, DatasetError> {
    match parser {
        DatasetSmilesParser::PubChem => {
            let (id, smiles) = line.split_once('\t').ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a CID<TAB>SMILES record".into(),
                }
            })?;
            Ok(DatasetSmilesRecord::new(id.to_owned(), smiles.to_owned()))
        }
        DatasetSmilesParser::MassSpecGym { smiles_column } => {
            let smiles = tsv_field(line, smiles_column).ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a TSV row with a smiles column value".into(),
                }
            })?;
            let id = tsv_field(line, 0).unwrap_or("");
            Ok(DatasetSmilesRecord::new(id.to_owned(), smiles.to_owned()))
        }
        DatasetSmilesParser::Zinc20 => {
            let mut fields = line.split_whitespace();
            let smiles = fields.next().ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a ZINC20 SMILES record".into(),
                }
            })?;
            let id = fields.next().ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a ZINC20 SMILES and identifier record".into(),
                }
            })?;
            if fields.next().is_some() {
                return Err(DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected exactly two whitespace-separated ZINC20 fields".into(),
                });
            }
            Ok(DatasetSmilesRecord::new(id.to_owned(), smiles.to_owned()))
        }
    }
}

fn tsv_field(line: &str, column_index: usize) -> Option<&str> {
    line.split('\t').nth(column_index)
}

fn collect_zinc20_smiles_paths(
    dataset_id: &'static str,
    path: &Path,
    output: &mut Vec<PathBuf>,
) -> Result<(), DatasetError> {
    if path.is_file() {
        if is_zinc20_smiles_file(path) {
            output.push(path.to_path_buf());
            return Ok(());
        }
        return Err(DatasetError::Format {
            dataset_id,
            line_number: 1,
            message: format!("expected an extracted ZINC20 directory, got {}", path.display()),
        });
    }

    let entries = fs::read_dir(path)
        .map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })?;
    let mut paths = Vec::new();
    for entry in entries {
        let entry =
            entry.map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })?;
        let entry_path = entry.path();
        if entry_path.is_dir() {
            collect_zinc20_smiles_paths(dataset_id, &entry_path, &mut paths)?;
        } else if is_zinc20_smiles_file(&entry_path) {
            paths.push(entry_path);
        }
    }
    paths.sort();
    output.extend(paths);
    Ok(())
}

fn is_zinc20_smiles_file(path: &Path) -> bool {
    path.file_name().and_then(|name| name.to_str()).is_some_and(|name| {
        name.starts_with("smiles_all_")
            && Path::new(name)
                .extension()
                .and_then(|extension| extension.to_str())
                .is_some_and(|extension| extension.eq_ignore_ascii_case("txt"))
    })
}

fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead + Send>, DatasetError> {
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
        DatasetCompression::TarGzip => {
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
                    let (was_downloaded, was_extracted) = ensure_extracted_tar_gzip(
                        dataset.url(),
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
                        was_decompressed: was_extracted,
                    })
                }
            }
        }
    }
}

fn fetch_dataset_collection<D: DatasetCollectionSource + ?Sized>(
    dataset: &D,
    options: &DatasetFetchOptions,
) -> Result<DatasetCollectionArtifact, DatasetError> {
    let cache_root = options.cache_dir.clone().unwrap_or_else(default_dataset_cache_dir);
    let dataset_dir = cache_root.join(dataset.id());
    create_dir_all(&dataset_dir)?;

    let mut paths = Vec::new();
    let mut compressed_paths = Vec::new();
    let mut was_downloaded = false;
    let mut was_extracted = false;

    for file in dataset.files() {
        let compressed_path = dataset_dir.join(file.file_name());
        let extracted_path = dataset_dir.join(file.extracted_file_name());
        match file.compression() {
            DatasetCompression::None => {
                was_downloaded |=
                    ensure_downloaded_url(file.url(), &compressed_path, options.cache_mode)?;
                paths.push(compressed_path.clone());
                compressed_paths.push(compressed_path);
            }
            DatasetCompression::Gzip => {
                match options.gzip_mode {
                    GzipMode::KeepCompressed => {
                        was_downloaded |= ensure_downloaded_url(
                            file.url(),
                            &compressed_path,
                            options.cache_mode,
                        )?;
                        paths.push(compressed_path.clone());
                        compressed_paths.push(compressed_path);
                    }
                    GzipMode::Decompress | GzipMode::KeepBoth => {
                        let (downloaded, decompressed) = ensure_decompressed_url(
                            file.url(),
                            &compressed_path,
                            &extracted_path,
                            options.cache_mode,
                        )?;
                        was_downloaded |= downloaded;
                        was_extracted |= decompressed;
                        paths.push(extracted_path);
                        if compressed_path.is_file() {
                            compressed_paths.push(compressed_path);
                        }
                    }
                }
            }
            DatasetCompression::TarGzip => {
                match options.gzip_mode {
                    GzipMode::KeepCompressed => {
                        was_downloaded |= ensure_downloaded_url(
                            file.url(),
                            &compressed_path,
                            options.cache_mode,
                        )?;
                        paths.push(compressed_path.clone());
                        compressed_paths.push(compressed_path);
                    }
                    GzipMode::Decompress | GzipMode::KeepBoth => {
                        let (downloaded, extracted) = ensure_extracted_tar_gzip(
                            file.url(),
                            &compressed_path,
                            &extracted_path,
                            options.cache_mode,
                        )?;
                        was_downloaded |= downloaded;
                        was_extracted |= extracted;
                        paths.push(extracted_path);
                        if compressed_path.is_file() {
                            compressed_paths.push(compressed_path);
                        }
                    }
                }
            }
        }
    }

    Ok(DatasetCollectionArtifact {
        dataset_id: dataset.id(),
        paths,
        compressed_paths,
        was_downloaded,
        was_extracted,
    })
}

fn ensure_downloaded<D: DatasetSource + ?Sized>(
    dataset: &D,
    target_path: &Path,
    cache_mode: CacheMode,
) -> Result<bool, DatasetError> {
    ensure_downloaded_url(dataset.url(), target_path, cache_mode)
}

fn ensure_downloaded_url(
    url: &'static str,
    target_path: &Path,
    cache_mode: CacheMode,
) -> Result<bool, DatasetError> {
    if cache_mode == CacheMode::UseCache && target_path.is_file() {
        return Ok(false);
    }

    download_to_path(url, target_path)?;
    Ok(true)
}

fn ensure_decompressed<D: DatasetSource + ?Sized>(
    dataset: &D,
    compressed_path: &Path,
    decompressed_path: &Path,
    cache_mode: CacheMode,
) -> Result<(bool, bool), DatasetError> {
    ensure_decompressed_url(dataset.url(), compressed_path, decompressed_path, cache_mode)
}

fn ensure_decompressed_url(
    url: &'static str,
    compressed_path: &Path,
    decompressed_path: &Path,
    cache_mode: CacheMode,
) -> Result<(bool, bool), DatasetError> {
    if cache_mode == CacheMode::UseCache && decompressed_path.is_file() {
        return Ok((false, false));
    }

    let was_downloaded = ensure_downloaded_url(url, compressed_path, cache_mode)?;
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

fn ensure_extracted_tar_gzip(
    url: &'static str,
    compressed_path: &Path,
    extracted_path: &Path,
    cache_mode: CacheMode,
) -> Result<(bool, bool), DatasetError> {
    if cache_mode == CacheMode::UseCache && extracted_path.is_dir() {
        return Ok((false, false));
    }

    let was_downloaded = ensure_downloaded_url(url, compressed_path, cache_mode)?;
    let was_extracted = untar_gzip_file(compressed_path, extracted_path)?;
    Ok((was_downloaded, was_extracted))
}

fn untar_gzip_file(compressed_path: &Path, extracted_path: &Path) -> Result<bool, DatasetError> {
    write_parent_dir(extracted_path)?;
    let source_file = File::open(compressed_path)
        .map_err(|source| DatasetError::Io { path: compressed_path.to_path_buf(), source })?;
    let progress_bar = new_byte_progress_bar(
        source_file.metadata().ok().map(|metadata| metadata.len()),
        &progress_label("extracting", extracted_path),
    );
    let source_file = ProgressReader::new(source_file, progress_bar.clone());
    let decoder = GzDecoder::new(source_file);
    let mut archive = Archive::new(decoder);

    let temporary_path = temporary_download_path(extracted_path);
    remove_path_if_exists(&temporary_path)?;
    create_dir_all(&temporary_path)?;
    if let Err(source) = archive.unpack(&temporary_path) {
        progress_bar.abandon();
        remove_path_if_exists(&temporary_path)?;
        return Err(DatasetError::Io { path: extracted_path.to_path_buf(), source });
    }
    progress_bar.finish_and_clear();

    let extracted_name = extracted_path
        .file_name()
        .unwrap_or_else(|| unreachable!("extracted path has a file name"));
    let unpacked_path = temporary_path.join(extracted_name);
    if !unpacked_path.exists() {
        remove_path_if_exists(&temporary_path)?;
        return Err(DatasetError::Io {
            path: extracted_path.to_path_buf(),
            source: io::Error::new(
                io::ErrorKind::NotFound,
                "tar archive did not contain the expected top-level directory",
            ),
        });
    }

    remove_path_if_exists(extracted_path)?;
    fs::rename(&unpacked_path, extracted_path)
        .map_err(|source| DatasetError::Io { path: extracted_path.to_path_buf(), source })?;
    remove_path_if_exists(&temporary_path)?;
    Ok(true)
}

fn remove_path_if_exists(path: &Path) -> Result<(), DatasetError> {
    if !path.exists() {
        return Ok(());
    }
    let result = if path.is_dir() { fs::remove_dir_all(path) } else { fs::remove_file(path) };
    result.map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })
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
        path::{Path, PathBuf},
        time::{SystemTime, UNIX_EPOCH},
        vec::Vec,
    };

    use flate2::{Compression, write::GzEncoder};

    use super::{
        CacheMode, DatasetArtifact, DatasetCollectionArtifact, DatasetCollectionSource,
        DatasetCompression, DatasetFetchOptions, DatasetSmilesIter, DatasetSmilesRecordIter,
        DatasetSource, GzipMode, MASS_SPEC_GYM_SMILES, PUBCHEM_SMILES, PubChemSmiles,
        ZINC20_EXPECTED_RECORD_COUNT, ZINC20_SMILES, Zinc20Smiles, default_dataset_cache_dir,
        gunzip_file, untar_gzip_file,
    };

    fn temporary_directory(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_else(|_| unreachable!("system time is after unix epoch"))
            .as_nanos();
        std::env::temp_dir().join(format!("smiles-parser-{name}-{unique}"))
    }

    fn write_zinc20_tar_gzip(path: &Path, chunk_dir: &str, contents: &[u8]) {
        let file = File::create(path).unwrap();
        let encoder = GzEncoder::new(file, Compression::default());
        let mut builder = tar::Builder::new(encoder);
        let mut header = tar::Header::new_gnu();
        header.set_size(
            u64::try_from(contents.len())
                .unwrap_or_else(|_| unreachable!("fixture length fits into u64")),
        );
        header.set_mode(0o644);
        header.set_cksum();
        builder
            .append_data(&mut header, format!("{chunk_dir}/smiles_all_01.txt"), contents)
            .unwrap();
        let encoder = builder.into_inner().unwrap();
        encoder.finish().unwrap();
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
    fn zinc20_smiles_metadata_matches_current_upstream_layout() {
        assert_eq!(ZINC20_SMILES.id(), "zinc20-smiles");
        assert_eq!(ZINC20_SMILES.first_chunk(), 1);
        assert_eq!(ZINC20_SMILES.last_chunk(), 20);
        assert_eq!(ZINC20_EXPECTED_RECORD_COUNT, 1_006_651_037);

        let files = ZINC20_SMILES.files();
        assert_eq!(files.len(), 20);
        assert_eq!(files[0].file_name(), "ZINC20_smiles_chunk_1.tar.gz");
        assert_eq!(files[0].extracted_file_name(), "ZINC20_smiles_chunk_1");
        assert_eq!(files[0].compression(), DatasetCompression::TarGzip);
        assert!(files[0].url().contains("zinc20-ML/smiles/ZINC20_smiles_chunk_1.tar.gz"));
        assert_eq!(files[19].file_name(), "ZINC20_smiles_chunk_20.tar.gz");
    }

    #[test]
    fn zinc20_chunk_range_validates_selection() {
        let chunk = Zinc20Smiles::chunk(7).unwrap();
        assert_eq!(chunk.first_chunk(), 7);
        assert_eq!(chunk.last_chunk(), 7);

        let range = Zinc20Smiles::chunk_range(3, 5).unwrap();
        assert_eq!(range.files().len(), 3);

        assert!(Zinc20Smiles::chunk(0).is_err());
        assert!(Zinc20Smiles::chunk_range(5, 3).is_err());
        assert!(Zinc20Smiles::chunk_range(1, 21).is_err());
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
        assert_eq!(ZINC20_SMILES.id(), "zinc20-smiles");
    }

    #[test]
    fn dataset_smiles_iterator_is_send() {
        fn assert_send<T: Send>() {}

        assert_send::<DatasetSmilesIter>();
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
    fn pubchem_record_iterator_streams_identifiers_and_smiles() {
        let directory = temporary_directory("datasets-pubchem-record-iter");
        fs::create_dir_all(&directory).unwrap();
        let compressed_path = directory.join("CID-SMILES.gz");

        {
            let file = File::create(&compressed_path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            encoder.write_all(b"123\tCCO\n456\tc1ccccc1\n").unwrap();
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
        let records = DatasetSmilesRecordIter::for_pubchem(&artifact)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(records[0].id(), "123");
        assert_eq!(records[0].smiles(), "CCO");
        assert_eq!(records[1].id(), "456");
        assert_eq!(records[1].smiles(), "c1ccccc1");

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
    fn zinc20_record_iterator_streams_records_from_extracted_chunk() {
        let directory = temporary_directory("datasets-zinc20-iter");
        let extracted_path = directory.join("ZINC20_smiles_chunk_1");
        fs::create_dir_all(&extracted_path).unwrap();
        fs::write(
            extracted_path.join("smiles_all_01.txt"),
            "CCO ZINC000000000001_1\nc1ccccc1 ZINC000000000002_1\n",
        )
        .unwrap();

        let artifact = DatasetCollectionArtifact {
            dataset_id: "zinc20-smiles",
            paths: vec![extracted_path],
            compressed_paths: Vec::new(),
            was_downloaded: false,
            was_extracted: false,
        };
        let records = DatasetSmilesRecordIter::for_zinc20(&artifact)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(records[0].smiles(), "CCO");
        assert_eq!(records[0].id(), "ZINC000000000001_1");
        assert_eq!(records[1].smiles(), "c1ccccc1");
        assert_eq!(records[1].id(), "ZINC000000000002_1");

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn zinc20_record_iterator_rejects_malformed_rows() {
        let directory = temporary_directory("datasets-zinc20-malformed");
        let extracted_path = directory.join("ZINC20_smiles_chunk_1");
        fs::create_dir_all(&extracted_path).unwrap();
        fs::write(extracted_path.join("smiles_all_01.txt"), "CCO\n").unwrap();

        let artifact = DatasetCollectionArtifact {
            dataset_id: "zinc20-smiles",
            paths: vec![extracted_path],
            compressed_paths: Vec::new(),
            was_downloaded: false,
            was_extracted: false,
        };
        match DatasetSmilesRecordIter::for_zinc20(&artifact).unwrap().next() {
            Some(Err(super::DatasetError::Format {
                dataset_id: "zinc20-smiles",
                line_number: 1,
                ..
            })) => {}
            other => panic!("unexpected result: {other:?}"),
        }

        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn untar_gzip_file_materializes_zinc20_chunk_directory() {
        let directory = temporary_directory("datasets-zinc20-untar");
        fs::create_dir_all(&directory).unwrap();
        let compressed_path = directory.join("ZINC20_smiles_chunk_1.tar.gz");
        let extracted_path = directory.join("ZINC20_smiles_chunk_1");
        write_zinc20_tar_gzip(
            &compressed_path,
            "ZINC20_smiles_chunk_1",
            b"CCO ZINC000000000001_1\n",
        );

        untar_gzip_file(&compressed_path, &extracted_path).unwrap();
        assert_eq!(
            fs::read_to_string(extracted_path.join("smiles_all_01.txt")).unwrap(),
            "CCO ZINC000000000001_1\n"
        );

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
