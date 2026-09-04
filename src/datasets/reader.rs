use alloc::{
    borrow::ToOwned,
    boxed::Box,
    string::{String, ToString},
    vec::Vec,
};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use flate2::read::GzDecoder;

use super::{
    types::{DatasetArtifact, DatasetCollectionArtifact, DatasetError, DatasetSmilesRecord},
    zinc20::collect_zinc20_smiles_paths,
};

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
    csv: Option<CsvReader>,
}

struct DatasetReader {
    path: PathBuf,
    reader: Box<dyn BufRead + Send>,
}

struct CsvReader {
    path: PathBuf,
    reader: csv::Reader<File>,
    record: csv::StringRecord,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum DatasetSmilesParser {
    PubChem,
    MassSpecGym { smiles_column: usize },
    Zinc20,
    Coconut { id_column: usize, smiles_column: usize },
    Lotus,
}

impl DatasetSmilesIter {
    pub(crate) fn from_records(inner: DatasetSmilesRecordIter) -> Self {
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
    pub(crate) fn for_pubchem(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        Self::from_artifact(artifact, DatasetSmilesParser::PubChem)
    }

    pub(crate) fn for_mass_spec_gym(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        let dataset_id = artifact.dataset_id();
        let path = artifact.path();
        let mut reader = open_text_reader(path)?;
        let mut header = String::new();
        let bytes_read = reader
            .read_line(&mut header)
            .map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })?;
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
            current: Some(DatasetReader { path: path.to_path_buf(), reader }),
            parser: DatasetSmilesParser::MassSpecGym { smiles_column },
            line_number: 1,
            line_buffer: String::new(),
            csv: None,
        })
    }

    pub(crate) fn for_zinc20(artifact: &DatasetCollectionArtifact) -> Result<Self, DatasetError> {
        let mut paths = Vec::new();
        for path in artifact.paths() {
            collect_zinc20_smiles_paths(artifact.dataset_id(), path, &mut paths)?;
        }
        if paths.is_empty() {
            return Err(DatasetError::Format {
                dataset_id: artifact.dataset_id(),
                line_number: 1,
                message: "expected at least one extracted ZINC20 smiles_all_*.txt file".into(),
            });
        }

        Ok(Self {
            dataset_id: artifact.dataset_id(),
            paths,
            next_path_index: 0,
            current: None,
            parser: DatasetSmilesParser::Zinc20,
            line_number: 0,
            line_buffer: String::new(),
            csv: None,
        })
    }

    pub(crate) fn for_coconut(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        let dataset_id = artifact.dataset_id();
        let path = artifact.path();
        let file = File::open(path)
            .map_err(|source| DatasetError::Io { path: path.to_path_buf(), source })?;
        let mut reader = csv::ReaderBuilder::new().has_headers(true).from_reader(file);
        let (id_column, smiles_column) = {
            let headers =
                reader.headers().map_err(|error| csv_error(dataset_id, path, 1, error))?;
            let id_column = headers
                .iter()
                .position(|field| field.eq_ignore_ascii_case("identifier"))
                .ok_or_else(|| {
                    DatasetError::Format {
                        dataset_id,
                        line_number: 1,
                        message: "expected a CSV header containing an identifier column".into(),
                    }
                })?;
            let smiles_column = headers
                .iter()
                .position(|field| field.eq_ignore_ascii_case("canonical_smiles"))
                .ok_or_else(|| {
                    DatasetError::Format {
                        dataset_id,
                        line_number: 1,
                        message: "expected a CSV header containing a canonical_smiles column"
                            .into(),
                    }
                })?;
            (id_column, smiles_column)
        };
        Ok(Self {
            dataset_id,
            paths: Vec::new(),
            next_path_index: 0,
            current: None,
            parser: DatasetSmilesParser::Coconut { id_column, smiles_column },
            line_number: 1,
            line_buffer: String::new(),
            csv: Some(CsvReader {
                path: path.to_path_buf(),
                reader,
                record: csv::StringRecord::new(),
            }),
        })
    }

    fn next_coconut_record(
        &mut self,
        id_column: usize,
        smiles_column: usize,
    ) -> Option<Result<DatasetSmilesRecord, DatasetError>> {
        let fallback_line = self.line_number.saturating_add(1);
        let state =
            self.csv.as_mut().unwrap_or_else(|| unreachable!("COCONUT CSV reader is initialized"));
        match state.reader.read_record(&mut state.record) {
            Ok(false) => None,
            Ok(true) => {
                let line_number = state
                    .record
                    .position()
                    .and_then(|position| usize::try_from(position.line()).ok())
                    .unwrap_or(fallback_line);
                self.line_number = line_number;
                let id = state.record.get(id_column).filter(|field| !field.is_empty()).ok_or_else(
                    || {
                        DatasetError::Format {
                            dataset_id: self.dataset_id,
                            line_number,
                            message: "expected a COCONUT CSV row with an identifier".into(),
                        }
                    },
                );
                let smiles =
                    state.record.get(smiles_column).filter(|field| !field.is_empty()).ok_or_else(
                        || {
                            DatasetError::Format {
                                dataset_id: self.dataset_id,
                                line_number,
                                message: "expected a COCONUT CSV row with a canonical_smiles value"
                                    .into(),
                            }
                        },
                    );
                match (id, smiles) {
                    (Ok(id), Ok(smiles)) => {
                        Some(Ok(DatasetSmilesRecord::new(id.to_owned(), smiles.to_owned())))
                    }
                    (Err(error), _) | (_, Err(error)) => Some(Err(error)),
                }
            }

            Err(error) => {
                let line_number = error
                    .position()
                    .and_then(|position| usize::try_from(position.line()).ok())
                    .unwrap_or(fallback_line);
                self.line_number = line_number;
                Some(Err(csv_error(self.dataset_id, &state.path, line_number, error)))
            }
        }
    }

    pub(crate) fn for_lotus(artifact: &DatasetArtifact) -> Result<Self, DatasetError> {
        Self::from_artifact(artifact, DatasetSmilesParser::Lotus)
    }

    fn from_artifact(
        artifact: &DatasetArtifact,
        parser: DatasetSmilesParser,
    ) -> Result<Self, DatasetError> {
        let path = artifact.path();
        Ok(Self {
            dataset_id: artifact.dataset_id(),
            paths: Vec::new(),
            next_path_index: 0,
            current: Some(DatasetReader {
                path: path.to_path_buf(),
                reader: open_text_reader(path)?,
            }),
            parser,
            line_number: 0,
            line_buffer: String::new(),
            csv: None,
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
        if let DatasetSmilesParser::Coconut { id_column, smiles_column } = self.parser {
            return self.next_coconut_record(id_column, smiles_column);
        }
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
        DatasetSmilesParser::Coconut { .. } => {
            unreachable!("COCONUT records are parsed by csv::Reader")
        }
        DatasetSmilesParser::Lotus => {
            // still use split_once to avoid iterator use (only one space)
            let (smiles, id) = line.split_once(' ').ok_or_else(|| {
                DatasetError::Format {
                    dataset_id,
                    line_number,
                    message: "expected a SMILES ID record".into(),
                }
            })?;
            Ok(DatasetSmilesRecord::new(id.to_owned(), smiles.to_owned()))
        }
    }
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

fn tsv_field(line: &str, column_index: usize) -> Option<&str> {
    line.split('\t').nth(column_index)
}

fn csv_error(
    dataset_id: &'static str,
    path: &Path,
    fallback_line: usize,
    error: csv::Error,
) -> DatasetError {
    let line_number = error
        .position()
        .and_then(|position| usize::try_from(position.line()).ok())
        .unwrap_or(fallback_line);
    let message = error.to_string();
    match error.into_kind() {
        csv::ErrorKind::Io(source) => DatasetError::Io { path: path.to_path_buf(), source },
        _ => DatasetError::Format { dataset_id, line_number, message },
    }
}
