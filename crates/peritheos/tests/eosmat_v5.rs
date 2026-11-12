use peritheos::{export_eosmat_v3, load_eosmat_str, material_from_value, validate_eosmat_document};
use serde_json::{json, Value};

fn fixture() -> Value {
    serde_json::from_str(include_str!("data/eosmat-v5-uncertainty.json")).unwrap()
}

#[test]
fn v5_is_lossless_and_exports_covariance_in_historical_model_units() {
    let document = fixture();
    let material = material_from_value(document.clone()).unwrap();
    let legacy = export_eosmat_v3(&document).unwrap();
    let scale = 0.060_221_407_6 / 4.0;
    let record = &legacy["eos_records"][0];
    assert!((record["parameter_errors"]["V0"].as_f64().unwrap() - 0.2 * scale).abs() < 1e-17);
    let matrix = &record["parameter_covariance"]["matrix"];
    assert!((matrix[0][0].as_f64().unwrap() - 0.04 * scale * scale).abs() < 1e-19);
    assert!((matrix[0][1].as_f64().unwrap() - 0.1 * scale).abs() < 1e-17);
    assert_eq!(matrix[1][1], 4);
    assert!(record.get("parameter_error_volume_basis").is_none());
    assert!(record["parameter_covariance"].get("volume_basis").is_none());
    let old = material_from_value(legacy.clone()).unwrap();
    for (v, t) in [(100.0, 300.0), (90.0, 1200.0)] {
        assert_eq!(
            material.eos_records[0].pressure(v, t).unwrap(),
            old.eos_records[0].pressure(v, t).unwrap()
        );
    }
    assert_eq!(export_eosmat_v3(&legacy).unwrap(), legacy);
    assert_eq!(
        load_eosmat_str(&material.to_json().unwrap())
            .unwrap()
            .document,
        document
    );
}

#[test]
fn v5_rejects_ambiguous_legacy_unit_flags_and_still_reads_v4() {
    let original = fixture();
    let mut wrong = original.clone();
    wrong["eos_records"][0]["parameter_error_volume_basis"] = json!("public");
    assert!(validate_eosmat_document(&wrong).is_err());
    let mut wrong = original.clone();
    wrong["eos_records"][0]["parameter_covariance"]["volume_basis"] = json!("public");
    assert!(validate_eosmat_document(&wrong).is_err());
    wrong["format_version"] = json!(4);
    validate_eosmat_document(&wrong).unwrap();
}
