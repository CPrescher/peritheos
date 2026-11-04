use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn original_sata_static_source_preserves_coefficients_and_existing_default() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/cementite_nonmagnetic.eosmat");
    if !path.is_file() {
        return; // Python catalog is not packaged in the standalone Rust crate.
    }
    let material = load_eosmat(path).unwrap();
    let record = material
        .record("cementite_nonmagnetic_sata_2010_bm3")
        .unwrap();
    assert_eq!(
        material.default_record().unwrap().identifier,
        "cementite_nonmagnetic_takahashi_2019_kbr_bm3_mgd"
    );
    assert_eq!(record.reference_temperature, 300.0);
    assert!((record.reference_volume() - 149.456).abs() < 1e-12);
    assert!((record.pressure(120.0, 300.0).unwrap() - 96.030_055_482_629_88).abs() < 1e-9);
    for pressure in [27.2, 50.0, 100.0, 186.6] {
        let volume = record.volume(pressure, 300.0).unwrap();
        assert!((record.pressure(volume, 300.0).unwrap() - pressure).abs() < 1e-8);
    }
    let raw = material.document["eos_records"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["identifier"] == record.identifier)
        .unwrap();
    assert_eq!(raw["model"], "BM3");
    assert_eq!(raw["parameters"]["K0"]["value"], 290.0);
    assert_eq!(raw["parameters"]["K0_prime"]["value"], 3.76);
    assert_eq!(
        raw["validity"]["pressure_gpa"],
        serde_json::json!([27.2, 186.6])
    );
    assert_eq!(
        raw["validity"]["temperature_k"],
        serde_json::json!([300, 300])
    );
    assert_eq!(raw["published_reference_fit"]["Pr_error_gpa"], 1.6);
    assert_eq!(
        raw["evidence"]["scientific_validation"]["data_agreement_status"],
        "poor"
    );
    assert_eq!(
        raw["evidence"]["scientific_validation"]["usage_recommendation"],
        "qualified_use"
    );
    assert_eq!(
        raw["evidence"]["reproduction"]["assessment"]["status"],
        "parity"
    );
    assert_eq!(
        raw["evidence"]["reproduction"]["assessment"]["scope"],
        "parameter_agreement"
    );
    assert_eq!(
        raw["evidence"]["reproduction"]["assessment"]["fit"]["data_agreement_status"],
        "poor"
    );
    assert!(raw.get("thermal").is_none());
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
