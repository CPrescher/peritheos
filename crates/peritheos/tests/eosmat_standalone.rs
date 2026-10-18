use std::path::Path;

use peritheos::{export_eosmat_v3, load_eosmat, load_eosmat_str, material_from_value};
use serde_json::json;

#[test]
fn minimal_python_authoring_template_is_a_regular_lossless_v4_document() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/examples/simple-material.eosmat");
    // The standalone Rust package does not include repository documentation.
    if !path.is_file() {
        return;
    }
    let material = load_eosmat(path).unwrap();
    let record = &material.eos_records[0];
    assert_eq!(material.document["format_version"], 4);
    assert_eq!(record.document["model"], "BM3");
    assert_eq!(record.document["parameters"]["K0_prime"]["fixed"], true);
    // Independent BM3 checkpoint with the template's K0'=4.
    let ratio: f64 = 10.0 / 9.0;
    let expected = 240.0 * (ratio.powf(7.0 / 3.0) - ratio.powf(5.0 / 3.0));
    assert!((record.pressure(9.0, 300.0).unwrap() - expected).abs() < 1e-12);
    let reloaded = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(reloaded.document, material.document);
}

fn native() -> serde_json::Value {
    json!({
        "format_version": 4, "identifier": "native", "formula": "X",
        "eos_records": [{"identifier": "native_bm3", "model": "BM3", "temperature_ref": 300,
            "parameters": {"V0": 10, "K0": {"value":160,"error":0}, "K0_prime":{"value":4,"fixed":true}},
            "evidence": {"scientific_validation":{"status":"not_reproduced", "note":"Synthetic diagnostic"}},
            "extension": {"rows":[[3,null],[1,2]]}}]
    })
}

#[test]
fn native_annotations_and_evidence_survive_read_save_and_legacy_export() {
    let document = native();
    let material = material_from_value(document.clone()).unwrap();
    assert_eq!(material.document, document);
    let legacy = export_eosmat_v3(&document).unwrap();
    assert_eq!(legacy["eos_records"][0]["parameter_errors"]["K0"], 0);
    assert_eq!(
        legacy["eos_records"][0]["fixed_parameters"],
        json!(["K0_prime"])
    );
    assert_eq!(
        legacy["eos_records"][0]["scientific_validation"]["status"],
        "not_reproduced"
    );
    let old = material_from_value(legacy).unwrap();
    assert_eq!(
        old.eos_records[0].pressure(9., 300.).unwrap().to_bits(),
        material.eos_records[0]
            .pressure(9., 300.)
            .unwrap()
            .to_bits()
    );
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, document);
}

#[test]
fn native_thermal_law_and_cell_normalization_use_the_existing_full_evaluator() {
    let mut document = native();
    document["formula_units_per_cell"] = json!(4);
    document["eos_records"][0]["thermal"] = json!({"model":"MieGruneisenDebye",
        "parameters":{"Tr":300,"theta0":760,"gamma0":1.5,"q":2,"n":2},
        "debye_temperature_law":"variable_exponent"});
    let new = material_from_value(document.clone()).unwrap();
    let old = material_from_value(export_eosmat_v3(&document).unwrap()).unwrap();
    assert_eq!(
        new.eos_records[0].pressure(9., 1000.).unwrap().to_bits(),
        old.eos_records[0].pressure(9., 1000.).unwrap().to_bits()
    );
    assert_eq!(new.document, document);
}

#[test]
fn native_conflicts_and_unknown_models_are_rejected() {
    let mut document = native();
    document["eos_records"][0]["model"] = json!("Unknown");
    assert!(material_from_value(document).is_err());
    let mut document = native();
    document["eos_records"][0]["evidence"]["identifier"] = json!("other");
    assert!(material_from_value(document).is_err());
    let mut document = native();
    document["eos_records"][0]["evidence"]["default_for"] = json!("equilibrium");
    assert!(material_from_value(document).is_err());
    let mut document = native();
    document["eos_records"][0]["parameters"]["V0"] = json!({"value":10,"fixed":1});
    assert!(material_from_value(document).is_err());
}

#[test]
fn optional_complete_migration_audit() {
    let Ok(directory) = std::env::var("PERITHEOS_V4_AUDIT") else {
        return;
    };
    let mut count = 0;
    for file in std::fs::read_dir(directory).unwrap() {
        let path = file.unwrap().path();
        if path.extension().and_then(|s| s.to_str()) != Some("eosmat") {
            continue;
        }
        let native =
            load_eosmat(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let old = material_from_value(export_eosmat_v3(&native.document).unwrap()).unwrap();
        assert_eq!(native.eos_records.len(), old.eos_records.len());
        assert_eq!(
            load_eosmat_str(&native.to_json().unwrap())
                .unwrap()
                .document,
            native.document
        );
        count += 1;
    }
    assert_eq!(count, 232);
}
