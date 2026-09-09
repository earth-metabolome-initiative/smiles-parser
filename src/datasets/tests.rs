use std::{
    fs::{self, File},
    io::Write,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
    vec::Vec,
};

use flate2::{Compression, write::GzEncoder};

use super::{
    CacheMode, DatasetFetchOptions, GzipMode, ZINC20_EXPECTED_RECORD_COUNT, Zinc20Smiles,
    fetch::{default_dataset_cache_dir, gunzip_file, untar_gzip_file},
    massspecgym::MASS_SPEC_GYM_SMILES,
    pubchem::{PUBCHEM_SMILES, PubChemSmiles},
    reader::{DatasetSmilesIter, DatasetSmilesRecordIter},
    source::{DatasetCollectionSource, DatasetSource},
    types::{DatasetArtifact, DatasetCollectionArtifact, DatasetCompression, DatasetError},
    zinc20::ZINC20_SMILES,
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
    builder.append_data(&mut header, format!("{chunk_dir}/smiles_all_01.txt"), contents).unwrap();
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
    let smiles =
        DatasetSmilesIter::from_records(DatasetSmilesRecordIter::for_pubchem(&artifact).unwrap())
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
    let smiles = DatasetSmilesIter::from_records(
        DatasetSmilesRecordIter::for_mass_spec_gym(&artifact).unwrap(),
    )
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
        Some(Err(DatasetError::Format { dataset_id: "zinc20-smiles", line_number: 1, .. })) => {}
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
    write_zinc20_tar_gzip(&compressed_path, "ZINC20_smiles_chunk_1", b"CCO ZINC000000000001_1\n");

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
    match DatasetSmilesRecordIter::for_mass_spec_gym(&artifact) {
        Ok(_) => panic!("expected a missing smiles header to fail"),
        Err(DatasetError::Format { dataset_id: "massspecgym-smiles", line_number: 1, .. }) => {}
        Err(error) => panic!("unexpected error: {error}"),
    }

    fs::remove_dir_all(directory).unwrap();
}
