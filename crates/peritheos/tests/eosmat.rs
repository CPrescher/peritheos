use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use peritheos::{
    load_eosmat, load_eosmat_str, save_eosmat, serialize_eosmat, validate_eosmat_document,
    EosmatError, EosmatErrorKind, Material,
};

fn materials_directory() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../peritheos/data/materials")
}

fn load_bundled_material(filename: &str) -> Option<Material> {
    let path = materials_directory().join(filename);
    if !path.is_file() {
        eprintln!(
            "skipping bundled-material integration test because {} is not packaged with the Rust crate",
            path.display()
        );
        return None;
    }
    Some(load_eosmat(path).unwrap())
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

    saved.validate().unwrap();
    assert!(saved.to_json().unwrap().ends_with('\n'));
    let path = std::env::temp_dir().join(format!(
        "peritheos-eosmat-method-save-{}-{unique}.eosmat",
        std::process::id()
    ));
    saved.save(&path).unwrap();
    fs::remove_file(path).unwrap();
}

#[test]
fn static_isothermal_record_accepts_zero_kelvin_reference_temperature() {
    let source = simple_document().replace(
        "\"default\": true,",
        "\"default\": true, \"equation_kind\": \"isothermal\", \"temperature_ref\": 0.0,",
    );
    let material = load_eosmat_str(&source).unwrap();
    let record = material.record("test_bm3").unwrap();

    assert!(record.reference_temperature.abs() < f64::EPSILON);
    assert!(record.pressure(9.0, 0.0).unwrap().is_finite());

    let error =
        load_eosmat_str(&source.replace("\"temperature_ref\": 0.0", "\"temperature_ref\": -1.0"))
            .unwrap_err();
    assert!(error.to_string().contains("static 0 K isothermal"));
}

#[test]
fn loaded_thermal_record_exposes_dac_forward_state_in_cell_units() {
    let Some(material) = load_bundled_material("diamond.eosmat") else {
        return;
    };
    let record = material
        .record("diamond_benedict_2014_double_debye_4")
        .unwrap();
    let cold_pressure = 60.0;
    let temperature = 2000.0;
    let f_dac = 0.25;

    let heated_volume = record
        .volume_with_dac_confinement(cold_pressure, temperature, f_dac)
        .unwrap();
    let thermal_increment = record
        .thermal_pressure_increment(heated_volume, temperature)
        .unwrap();
    let confinement_pressure = record
        .dac_thermal_pressure(heated_volume, temperature, f_dac)
        .unwrap();

    assert_close(confinement_pressure, f_dac * thermal_increment, 1.0e-12);
    assert_close(
        record.pressure(heated_volume, temperature).unwrap(),
        cold_pressure + confinement_pressure,
        1.0e-10,
    );

    let isothermal = material.record("diamond_dewaele_2008_vinet_2").unwrap();
    assert!(isothermal
        .thermal_pressure_increment(40.0, temperature)
        .is_err());
    assert!(isothermal
        .dac_thermal_pressure(40.0, temperature, f_dac)
        .is_err());
    assert!(isothermal
        .volume_with_dac_confinement(cold_pressure, temperature, f_dac)
        .is_err());
}

#[test]
fn bundled_suzuki_epsilon_feooh_uses_reference_temperature_expansivity() {
    let Some(material) = load_bundled_material("e_feooh.eosmat") else {
        return;
    };
    let record = material
        .record("e_feooh_suzuki_2016_bm3_thermal_2")
        .unwrap();

    assert_close(record.pressure(66.278, 300.0).unwrap(), 0.0, 1.0e-12);
    assert_close(
        record.pressure(62.63, 700.0).unwrap(),
        10.833_562_803_310_029,
        1.0e-11,
    );
}

#[test]
fn qin_2023_calcium_ferrite_records_load_and_reproduce_high_pressure_states() {
    let cases = [
        (
            "na093al102si100o4_calcium_ferrite.eosmat",
            "na093al102si100o4_calcium_ferrite_qin_2023_bm3_1",
            207.3,
            40.609_612_068_923_79,
        ),
        (
            "na088al099fe013si094o4_calcium_ferrite.eosmat",
            "na088al099fe013si094o4_calcium_ferrite_qin_2023_bm3_1",
            206.6,
            43.209_378_809_850_9,
        ),
    ];

    for (filename, identifier, volume, expected_pressure) in cases {
        let Some(material) = load_bundled_material(filename) else {
            return;
        };
        let record = material.record(identifier).unwrap();

        assert_eq!(
            record.eos.isothermal_model_identifier(),
            "birch_murnaghan_3"
        );
        assert_close(
            record.pressure(volume, 293.0).unwrap(),
            expected_pressure,
            1.0e-12,
        );
        assert_close(
            record.volume(expected_pressure, 293.0).unwrap(),
            volume,
            1.0e-12,
        );
    }
}

#[test]
fn luo_mgo_eosmat_preserves_absolute_thermal_pressure() {
    let Some(material) = load_bundled_material("mgo.eosmat") else {
        return;
    };
    let record = material.record("mgo_b1_luo_2023_vinet_thermal_5").unwrap();
    let ambient_cell_volume = 74.569_767_758_6;
    let table_volume = ambient_cell_volume * (1.0 - 0.42);

    assert_eq!(record.eos.isothermal_model_identifier(), "vinet");
    assert_eq!(
        record.eos.thermal_model_identifier(),
        Some("second_order_taylor_thermal_pressure")
    );
    assert_close(
        record.pressure(table_volume, 8500.0).unwrap(),
        342.01,
        1.5 / 342.01,
    );

    let cold_reference_volume = 74.074_102_512_3;
    let baseline = record.pressure(cold_reference_volume, 300.0).unwrap();
    assert!(baseline.abs() > 1.0e-3);
    assert_close(
        record
            .thermal_pressure_increment(table_volume, 300.0)
            .unwrap(),
        0.0,
        1.0e-12,
    );
}

#[test]
fn double_debye_eosmat_reference_temperature_anchors_the_dewaele_isotherm() {
    let Some(material) = load_bundled_material("diamond.eosmat") else {
        return;
    };
    let reference = material.record("diamond_dewaele_2008_vinet_2").unwrap();
    let cases = [
        (
            "diamond_benedict_2014_double_debye_4",
            "diamond_benedict_2014_dewaele_anchored",
        ),
        (
            "diamond_correa_2008_double_debye_log_moment_5",
            "diamond_correa_2008_dewaele_anchored",
        ),
    ];
    let volume = 40.0;
    let temperature = 3000.0;

    for (absolute_identifier, anchored_identifier) in cases {
        let absolute = material.record(absolute_identifier).unwrap();
        let anchored = material.record(anchored_identifier).unwrap();
        let reference_pressure = reference.pressure(volume, 298.0).unwrap();
        let expected = reference_pressure + absolute.pressure(volume, temperature).unwrap()
            - absolute.pressure(volume, 298.0).unwrap();

        assert_close(
            anchored.pressure(volume, 298.0).unwrap(),
            reference_pressure,
            1.0e-12,
        );
        assert_close(
            anchored.pressure(volume, temperature).unwrap(),
            expected,
            1.0e-12,
        );
    }
}

fn assert_close(actual: f64, expected: f64, relative_tolerance: f64) {
    assert!(
        (actual - expected).abs() <= relative_tolerance * expected.abs().max(1.0),
        "actual {actual:.17e} differs from expected {expected:.17e}"
    );
}

#[test]
fn decoding_and_io_errors_retain_kinds_codes_and_sources() {
    let json_error = load_eosmat_str("{").unwrap_err();
    assert_eq!(json_error.kind(), EosmatErrorKind::Json);
    assert_eq!(json_error.code(), "eosmat.json");
    assert!(json_error.source().is_some());
    assert_eq!(json_error.record_identifier(), None);
    assert!(json_error.to_string().contains("invalid eosmat JSON"));

    let missing = std::env::temp_dir().join(format!(
        "peritheos-missing-{}-never-created.eosmat",
        std::process::id()
    ));
    let io_error = load_eosmat(missing).unwrap_err();
    assert_eq!(io_error.kind(), EosmatErrorKind::Io);
    assert_eq!(io_error.code(), "eosmat.io");
    assert!(io_error.source().is_some());
    assert_eq!(io_error.record_identifier(), None);
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

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["eos_records"][0]["record_kind"] = serde_json::json!("refit");
    let error = validate_eosmat_document(&document).unwrap_err();
    assert!(error.to_string().contains("refit records require"));
}

#[test]
fn public_validation_reports_each_document_section_with_context() {
    let base: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    let mut invalid_documents = Vec::new();

    for (pointer, replacement, expected) in [
        ("/format", serde_json::json!("unknown"), "supported formats"),
        ("/name", serde_json::Value::Null, "name must be a string"),
        ("/identifier", serde_json::json!(""), "non-empty identifier"),
        (
            "/eos_records",
            serde_json::json!({}),
            "eos_records must be a JSON array",
        ),
        (
            "/eos_records/0/label",
            serde_json::Value::Null,
            "label must be a string",
        ),
        (
            "/eos_records/0/identifier",
            serde_json::json!(""),
            "identifier must be a non-empty string",
        ),
        (
            "/eos_records/0/reference",
            serde_json::json!(7),
            "reference must be a string or object",
        ),
        (
            "/eos_records/0/parameter_errors",
            serde_json::json!([]),
            "parameter_errors must be a JSON object",
        ),
        (
            "/eos_records/0/fixed_parameters",
            serde_json::json!({}),
            "fixed_parameters must be a JSON array",
        ),
        (
            "/eos_records/0/scientific_validation/status",
            serde_json::json!("unknown"),
            "scientific_validation.status is invalid",
        ),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).unwrap() = replacement;
        invalid_documents.push((document, expected));
    }

    assert!(matches!(
        validate_eosmat_document(&serde_json::Value::Null),
        Err(EosmatError::InvalidDocument(_))
    ));
    for (document, expected) in invalid_documents {
        let error = validate_eosmat_document(&document).unwrap_err();
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?} in {error}"
        );
    }
}

#[test]
fn public_validation_reports_optional_material_and_record_sections() {
    let base: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    let mut invalid_documents = Vec::new();

    for (field, replacement, expected) in [
        (
            "formula_units_per_cell",
            serde_json::json!(0.0),
            "must be greater than zero",
        ),
        ("aliases", serde_json::json!({}), "must be a JSON array"),
        ("lattice", serde_json::json!({}), "lattice.a is required"),
    ] {
        let mut document = base.clone();
        document[field] = replacement;
        invalid_documents.push((document, expected));
    }

    for (field, replacement, expected) in [
        (
            "parameter_error_confidence",
            serde_json::json!(1.0),
            "must lie between zero and one",
        ),
        (
            "volume",
            serde_json::json!([]),
            "volume must be a JSON object",
        ),
        (
            "parameter_covariance",
            serde_json::json!({"parameter_order": ["K0"]}),
            "matrix is required",
        ),
        (
            "thermal",
            serde_json::json!([]),
            "thermal must be a JSON object",
        ),
        (
            "experimental_pressure_range_gpa",
            serde_json::json!([2.0, 1.0]),
            "must be ordered",
        ),
        (
            "validity",
            serde_json::json!([]),
            "validity must be a JSON object",
        ),
    ] {
        let mut document = base.clone();
        document["eos_records"][0][field] = replacement;
        invalid_documents.push((document, expected));
    }

    for (document, expected) in invalid_documents {
        let error = validate_eosmat_document(&document).unwrap_err();
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?} in {error}"
        );
    }

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["phase"] = serde_json::json!(7);
    assert!(validate_eosmat_document(&document)
        .unwrap_err()
        .to_string()
        .contains("phase must be a string"));

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["eos_records"][0]["validity"] = serde_json::json!({
        "pressure_gpa": [0.0],
        "temperature_k": [300.0, 300.0],
        "notes": []
    });
    assert!(validate_eosmat_document(&document)
        .unwrap_err()
        .to_string()
        .contains("must contain two values"));

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["eos_records"][0]
        .as_object_mut()
        .unwrap()
        .remove("parameter_errors");
    assert!(validate_eosmat_document(&document)
        .unwrap_err()
        .to_string()
        .contains("parameter_errors is required"));

    let mut document: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    document["eos_records"][0]["fixed_parameters"] = serde_json::json!([7]);
    assert!(validate_eosmat_document(&document)
        .unwrap_err()
        .to_string()
        .contains("fixed_parameters must contain strings"));
}

#[test]
fn public_validation_reports_pressure_calibration_errors() {
    let mut base: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    base["eos_records"][0]["pressure_calibration"] = serde_json::json!({
        "status": "resolved",
        "methods": [{
            "kind": "ruby_fluorescence",
            "reference": "Synthetic ruby calibration",
            "reference_calibration_record": "ruby_mao_1986",
            "source_location": "Methods"
        }],
        "recalculation": {"status": "missing_calibrant_observations"}
    });
    validate_eosmat_document(&base).unwrap();

    let mut cases = Vec::new();
    for (pointer, replacement, expected) in [
        (
            "/eos_records/0/pressure_calibration/status",
            serde_json::json!("invalid"),
            "status is invalid",
        ),
        (
            "/eos_records/0/pressure_calibration/methods",
            serde_json::json!({}),
            "methods must be a JSON array",
        ),
        (
            "/eos_records/0/pressure_calibration/methods/0/kind",
            serde_json::json!("invalid"),
            "kind is invalid",
        ),
        (
            "/eos_records/0/pressure_calibration/methods/0/source_location",
            serde_json::json!(""),
            "source_location must be a non-empty string",
        ),
        (
            "/eos_records/0/pressure_calibration/methods/0/reference",
            serde_json::json!([]),
            "reference must be a string or object",
        ),
        (
            "/eos_records/0/pressure_calibration/methods/0/reference_calibration_record",
            serde_json::json!(""),
            "reference_calibration_record must be a non-empty string",
        ),
        (
            "/eos_records/0/pressure_calibration/recalculation",
            serde_json::json!([]),
            "recalculation must be a JSON object",
        ),
        (
            "/eos_records/0/pressure_calibration/recalculation/status",
            serde_json::json!("invalid"),
            "recalculation.status is invalid",
        ),
    ] {
        let mut document = base.clone();
        *document.pointer_mut(pointer).unwrap() = replacement;
        cases.push((document, expected));
    }

    let mut document = base.clone();
    document["eos_records"][0]["pressure_calibration"]["methods"] = serde_json::json!([]);
    cases.push((document, "methods must not be empty"));

    let mut document = base.clone();
    let method = &mut document["eos_records"][0]["pressure_calibration"]["methods"][0];
    method["reference_eos_record"] = serde_json::json!("");
    cases.push((document, "reference_eos_record must be a non-empty string"));

    let mut document = base.clone();
    let method = &mut document["eos_records"][0]["pressure_calibration"]["methods"][0];
    method["reference_eos_record"] = serde_json::json!("test_bm3");
    cases.push((document, "requires an equation_of_state method"));

    let mut document = base.clone();
    let method = &mut document["eos_records"][0]["pressure_calibration"]["methods"][0];
    method["kind"] = serde_json::json!("ambient_pressure");
    cases.push((document, "requires a ruby_fluorescence method"));

    let mut document = base.clone();
    let method = &mut document["eos_records"][0]["pressure_calibration"]["methods"][0];
    method["kind"] = serde_json::json!("equation_of_state");
    method.as_object_mut().unwrap().remove("reference");
    method
        .as_object_mut()
        .unwrap()
        .remove("reference_calibration_record");
    cases.push((document, "reference is required for an EOS method"));

    for (document, expected) in cases {
        let error = validate_eosmat_document(&document).unwrap_err();
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?} in {error}"
        );
    }
}

#[test]
fn canonical_document_loads_as_an_executable_model_and_preserves_extensions() {
    let material = load_eosmat_str(simple_document()).unwrap();
    assert_eq!(material.identifier, "test_material");
    assert_eq!(material.document["extension"]["preserved"], true);

    let record = material.default_record().unwrap();
    assert_eq!(material.hugoniot_records().count(), 0);
    assert_eq!(material.equilibrium_records().count(), 1);
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
#[allow(clippy::too_many_lines)]
fn canonical_hugoniot_record_loads_as_an_executable_eos_path() {
    let source = r#"{
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": "shock_test",
        "name": "Shock test",
        "formula": "X",
        "phase": "alpha",
        "formula_units_per_cell": 1.0,
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell"
        },
        "eos_records": [{
            "identifier": "shock_alpha",
            "label": "Alpha principal Hugoniot",
            "equation_kind": "hugoniot",
            "default_for": "hugoniot",
            "loading_path": "principal",
            "branch_kind": "untransformed",
            "initial_state": {
                "phase": "alpha",
                "material_identifier": "shock_test",
                "temperature_k": 298.15,
                "pressure_gpa": 0.0,
                "density_g_cm3": 8.0
            },
            "volume_basis": {
                "kind": "formula_units",
                "formula_units": 1.0,
                "molar_mass_g_mol": 48.17712608
            },
            "branch_domain": {
                "particle_velocity_km_s": [0.0, 3.0],
                "kind": "phase_stability",
                "boundary_status": "reported_exactly"
            },
            "reference": "Synthetic test record",
            "eos": {
                "type": "LinearUsUpHugoniot",
                "model": "linear_us_up_hugoniot",
                "parameters": {
                    "V0": 10.0,
                    "rho0": 8.0,
                    "c0": 4.0,
                    "s": 1.5,
                    "P0": 0.0
                }
            },
            "parameter_errors": {},
            "fixed_parameters": [],
            "scientific_validation": {"status": "primary_source_validated"}
        }]
    }"#;
    let material = load_eosmat_str(source).unwrap();
    let record = material.default_record().unwrap();
    assert_eq!(material.hugoniot_records().count(), 1);
    assert_eq!(material.equilibrium_records().count(), 0);
    assert!(record.eos.is_hugoniot());
    assert!(!record.eos.is_thermal());
    assert_eq!(record.eos.model_identifier(), "linear_us_up_hugoniot");
    assert!((record.reference_temperature - 298.15).abs() < f64::EPSILON);
    let pressure = record.pressure(8.0, 298.15).unwrap();
    assert!(record.pressure(8.0, 300.0).is_err());
    assert!(pressure > 0.0);
    assert!((record.volume(pressure, 298.15).unwrap() - 8.0).abs() < 1.0e-10);
    assert!(record.bulk_modulus(8.0, 298.15).is_err());
    let typed = material.default_hugoniot_record().unwrap();
    assert_eq!(
        typed.metadata.loading_path,
        peritheos::HugoniotLoadingPath::Principal
    );
    assert!((typed.density(8.0).unwrap() - 10.0).abs() < 1.0e-12);
    assert!(typed.shock_velocity(8.0).unwrap() > 0.0);
    assert!(typed.particle_velocity(8.0).unwrap() > 0.0);
    assert!(typed.specific_internal_energy_change(8.0).unwrap() > 0.0);
    assert!(typed.tangent_modulus(8.0).unwrap() > 0.0);
    let state = typed.state_from_particle_velocity(1.0).unwrap();
    assert!((state.particle_velocity - 1.0).abs() < f64::EPSILON);
    assert!(typed.state_from_particle_velocity(3.1).is_err());

    let mut mixed: serde_json::Value = serde_json::from_str(simple_document()).unwrap();
    let hugoniot_document: serde_json::Value = serde_json::from_str(source).unwrap();
    mixed["identifier"] = "shock_test".into();
    mixed["phase"] = "alpha".into();
    mixed["formula_units_per_cell"] = 1.0.into();
    mixed["eos_records"]
        .as_array_mut()
        .unwrap()
        .push(hugoniot_document["eos_records"][0].clone());
    let mixed = load_eosmat_str(&mixed.to_string()).unwrap();
    assert!(!mixed.default_record().unwrap().eos.is_hugoniot());
    assert!(mixed.default_equilibrium_record().is_some());
    assert!(mixed.default_hugoniot_record().is_some());

    let mut derived = hugoniot_document.clone();
    derived["eos_records"][0]["record_kind"] = "derived".into();
    assert!(validate_eosmat_document(&derived).is_err());
    derived["eos_records"][0]["derivation"] = serde_json::json!({
        "source_kind": "sesame_table",
        "source_identifier": "SESAME 1234, release 2026-01",
        "method": "sampled the principal Hugoniot and fit linear Us-up coefficients"
    });
    validate_eosmat_document(&derived).unwrap();

    let mut precompressed_transformed = hugoniot_document;
    precompressed_transformed["identifier"] = "shock_beta".into();
    precompressed_transformed["phase"] = "beta".into();
    precompressed_transformed["eos_records"][0]["loading_path"] = "precompressed".into();
    precompressed_transformed["eos_records"][0]["branch_kind"] = "transformed".into();
    precompressed_transformed["eos_records"][0]["initial_state"]["material_identifier"] =
        "shock_alpha".into();
    precompressed_transformed["eos_records"][0]["initial_state"]["pressure_gpa"] = 5.0.into();
    precompressed_transformed["eos_records"][0]["eos"]["parameters"]["P0"] = 5.0.into();
    let precompressed_transformed =
        load_eosmat_str(&precompressed_transformed.to_string()).unwrap();
    let metadata = &precompressed_transformed
        .default_hugoniot_record()
        .unwrap()
        .metadata;
    assert_eq!(
        metadata.loading_path,
        peritheos::HugoniotLoadingPath::Precompressed
    );
    assert_eq!(
        metadata.branch_kind,
        peritheos::HugoniotBranchKind::Transformed
    );

    for (old, new, expected) in [
        (
            "\"density_g_cm3\": 8.0",
            "\"density_g_cm3\": 7.9",
            "must match eos.parameters.rho0",
        ),
        (
            "\"formula_units\": 1.0,",
            "\"formula_units\": 2.0,",
            "must match formula_units_per_cell",
        ),
        (
            "\"loading_path\": \"principal\"",
            "\"loading_path\": \"reshock\"",
            "supported loading_path",
        ),
    ] {
        let error = load_eosmat_str(&source.replace(old, new)).unwrap_err();
        assert!(error.to_string().contains(expected));
    }
}

#[test]
#[allow(clippy::too_many_lines)]
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
            let evaluation_volume = record.as_hugoniot().map_or(reference_volume, |hugoniot| {
                let [lower, upper] = hugoniot.metadata.branch_domain.particle_velocity_km_s;
                hugoniot
                    .state_from_particle_velocity(0.5 * (lower + upper))
                    .unwrap_or_else(|error| {
                        panic!(
                            "{} / {} failed inside its declared Hugoniot branch: {error}",
                            material.identifier, record.identifier
                        )
                    })
                    .volume
            });
            let pressure = record
                .pressure(evaluation_volume, record.reference_temperature)
                .unwrap_or_else(|error| {
                    panic!(
                        "{} / {} failed at its catalog validation state: {error}",
                        material.identifier, record.identifier
                    )
                });
            if record.as_hugoniot().is_some() {
                assert!(pressure.is_finite());
                assert!(pressure > 0.0);
                continue;
            }
            if record
                .document
                .pointer("/thermal/thermal_pressure_reference")
                .and_then(serde_json::Value::as_str)
                == Some("absolute_zero")
            {
                assert!(pressure.is_finite());
                assert!(pressure > 0.0);
                continue;
            }
            match (
                record.identifier.as_str(),
                record.eos.thermal_model_identifier(),
            ) {
                ("diamond_benedict_2014_double_debye_4", Some("double_debye_helmholtz")) => {
                    assert!(pressure.is_finite());
                    assert!(pressure > 0.0);
                    let pressure = record
                        .pressure(8.0 * 4.654_270_411_587_497, 3000.0)
                        .unwrap();
                    assert!((pressure - 150.0).abs() < 1.0e-7);
                }
                (
                    "diamond_correa_2008_double_debye_log_moment_5",
                    Some("double_debye_log_moment_helmholtz"),
                ) => {
                    assert!(pressure.is_finite());
                    assert!(pressure > 0.0);
                    let pressure = record.pressure(8.0 * 4.43, 5000.0).unwrap();
                    assert!((pressure - 202.628_115_197_741_86).abs() < 1.0e-7);
                }
                (
                    "mgo_b1_luo_2023_vinet_thermal_5",
                    Some("second_order_taylor_thermal_pressure"),
                ) => {
                    assert_close(pressure, 0.785_335_88, 1.0e-8);
                }
                (_, _) => {
                    assert!(pressure.abs() < 1.0e-8, "reference pressure was {pressure}");
                }
            }
        }
    }

    assert_eq!(paths.len(), 182);
    assert_eq!(records, 353);
    assert_eq!(thermal_records, 56);
}

#[test]
fn datchi_cbn_absolute_mgd_reproduces_table_vi_volume() {
    let Some(material) = load_bundled_material("boron_nitride.eosmat") else {
        return;
    };
    let record = material
        .eos_records
        .iter()
        .find(|record| record.identifier == "boron_nitride_datchi_2007_vinet_mgd_2")
        .unwrap();

    let atomic_volume = record.volume(0.0, 300.0).unwrap() / 8.0;
    assert!((atomic_volume - 5.9055).abs() < 1.0e-5);
    assert_eq!(
        record
            .document
            .pointer("/thermal/thermal_pressure_reference")
            .and_then(serde_json::Value::as_str),
        Some("absolute_zero")
    );
}
