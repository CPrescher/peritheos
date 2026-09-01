use std::fs;
use std::path::{Path, PathBuf};

use peritheos_core::{load_eosmat, load_eosmat_str, EosmatError};

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
fn all_bundled_material_records_load_through_the_rust_registry() {
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
