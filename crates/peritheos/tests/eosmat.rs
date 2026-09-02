use std::fs;
use std::path::{Path, PathBuf};

use peritheos::{
    load_eosmat, load_eosmat_str, save_eosmat, serialize_eosmat, validate_eosmat_document,
    EosmatError, EosmatErrorKind,
};

fn materials_directory() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../peritheos/data/materials")
}

fn simple_document() -> &'static str {
    r#"{
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": "test_material",
        "name": "Test material",
        "formula": "X",
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell"
        },
        "extension": {"preserved": true},
        "eos_records": [{
            "identifier": "test_bm3",
            "label": "Test BM3",
            "default": true,
            "reference": "Synthetic test record",
            "eos": {
                "type": "BM3",
                "model": "birch_murnaghan_3",
                "parameters": {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0}
            },
            "parameter_errors": {},
            "fixed_parameters": [],
            "scientific_validation": {"status": "primary_source_validated"},
            "record_extension": 42
        }]
    }"#
}

#[test]
fn canonical_documents_validate_serialize_and_round_trip_without_losing_extensions() {
    let document = serde_json::from_str(simple_document()).unwrap();
    validate_eosmat_document(&document).unwrap();

    let serialized = serialize_eosmat(&document).unwrap();
    assert!(serialized.ends_with('\n'));
    let reloaded = load_eosmat_str(&serialized).unwrap();
    assert_eq!(reloaded.document, document);
    assert_eq!(reloaded.document["extension"]["preserved"], true);
    assert_eq!(reloaded.eos_records[0].document["record_extension"], 42);

    let unique = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "peritheos-eosmat-roundtrip-{}-{unique}.eosmat",
        std::process::id()
    ));
    save_eosmat(&path, &document).unwrap();
    let saved = load_eosmat(&path).unwrap();
    fs::remove_file(path).unwrap();
    assert_eq!(saved.document, document);
}

#[test]
fn public_validation_rejects_ambiguous_defaults_and_covariance_shapes() {
    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    let duplicate = document["eos_records"][0].clone();
    document["eos_records"]
        .as_array_mut()
        .unwrap()
        .push(duplicate);
    let error = validate_eosmat_document(&document).unwrap_err();
    assert_eq!(error.kind(), EosmatErrorKind::InvalidDocument);
    assert_eq!(error.code(), "eosmat.invalid_document");
    assert!(error
        .to_string()
        .contains("duplicate EOS record identifier"));

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["eos_records"][0]["parameter_covariance"] = serde_json::json!({
        "parameter_order": ["K0", "K0_prime"],
        "matrix": [[1.0, 0.0]]
    });
    let error = validate_eosmat_document(&document).unwrap_err();
    assert!(error.to_string().contains("square matrix"));
}

#[test]
fn canonical_document_loads_as_an_executable_model_and_preserves_extensions() {
    let material = load_eosmat_str(simple_document()).unwrap();
    assert_eq!(material.identifier, "test_material");
    assert_eq!(material.document["extension"]["preserved"], true);

    let record = material.default_record().unwrap();
    assert_eq!(record.identifier, "test_bm3");
    assert_eq!(
        record.eos.isothermal_model_identifier(),
        "birch_murnaghan_3"
    );
    assert_eq!(record.eos.thermal_model_identifier(), None);
    assert_eq!(record.document["record_extension"], 42);
    assert!((record.reference_volume() - 10.0).abs() < f64::EPSILON);
    assert!(record.pressure(10.0, 900.0).unwrap().abs() < f64::EPSILON);
    assert!(record.pressure(9.0, 300.0).unwrap() > 0.0);

    let pressure = record.pressure(8.5, 300.0).unwrap();
    let recovered = record.volume(pressure, 300.0).unwrap();
    assert!((recovered - 8.5).abs() < 1.0e-10);
}

#[test]
fn legacy_dioptas_format_two_is_accepted_and_identifiers_are_generated() {
    let source = r#"{
        "format_version": 2,
        "name": "Legacy Test",
        "formula": "X",
        "eos_records": [{
            "label": "Legacy Vinet",
            "eos": {
                "type": "Vinet",
                "parameters": {"V0": 12.0, "K0": 150.0, "K0_prime": 4.2}
            },
            "thermal": {
                "type": "LinearThermalPressure",
                "parameters": {"Tr": 300.0, "alpha_KT": 0.005}
            }
        }]
    }"#;
    let material = load_eosmat_str(source).unwrap();
    assert_eq!(material.identifier, "legacy_test");
    assert_eq!(material.eos_records[0].identifier, "legacy_vinet");
    assert_eq!(
        material.eos_records[0].eos.thermal_model_identifier(),
        Some("linear_thermal_pressure")
    );
    assert!((material.eos_records[0].reference_volume() - 12.0).abs() < f64::EPSILON);
    assert!((material.eos_records[0].pressure(12.0, 500.0).unwrap() - 1.0).abs() < 1.0e-12);
}

#[test]
fn wrong_type_model_pair_is_rejected() {
    let source = simple_document().replace("\"type\": \"BM3\"", "\"type\": \"Vinet\"");
    let error = load_eosmat_str(&source).unwrap_err();
    assert!(matches!(error, EosmatError::InvalidRecord { .. }));
    assert!(error.to_string().contains("requires type \"BM3\""));
}

#[test]
fn all_bundled_material_records_load_and_round_trip_through_rust() {
    let materials_directory = materials_directory();
    if !materials_directory.is_dir() {
        // The crates.io archive deliberately contains the reusable loader but
        // not the Python distribution's complete material catalog.
        return;
    }
    let mut paths = fs::read_dir(materials_directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "eosmat")
        })
        .collect::<Vec<_>>();
    paths.sort();

    let mut records = 0;
    let mut thermal_records = 0;
    for path in &paths {
        let material = load_eosmat(path)
            .unwrap_or_else(|error| panic!("{} failed to load: {error}", path.display()));
        let serialized = material
            .to_json()
            .unwrap_or_else(|error| panic!("{} failed to serialize: {error}", path.display()));
        let round_tripped = load_eosmat_str(&serialized)
            .unwrap_or_else(|error| panic!("{} failed to round trip: {error}", path.display()));
        assert_eq!(round_tripped.document, material.document);
        records += material.eos_records.len();
        thermal_records += material
            .eos_records
            .iter()
            .filter(|record| record.eos.is_thermal())
            .count();
        for record in &material.eos_records {
            let reference_volume = record.reference_volume();
            let pressure = record
                .pressure(reference_volume, record.reference_temperature)
                .unwrap_or_else(|error| {
                    panic!(
                        "{} / {} failed at its reference state: {error}",
                        material.identifier, record.identifier
                    )
                });
            if record.eos.thermal_model_identifier() == Some("double_debye_helmholtz") {
                assert!(pressure.is_finite());
                assert!(pressure > 0.0);
                let pressure = record
                    .pressure(8.0 * 4.654_270_411_587_497, 3000.0)
                    .unwrap();
                assert!((pressure - 150.0).abs() < 1.0e-7);
            } else {
                assert!(pressure.abs() < 1.0e-8, "reference pressure was {pressure}");
            }
        }
    }

    assert_eq!(paths.len(), 115);
    assert_eq!(records, 147);
    assert_eq!(thermal_records, 28);
}
