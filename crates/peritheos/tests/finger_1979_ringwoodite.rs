use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn installed_finger_static_curve_keeps_its_reference_and_metadata() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/ferrous_ringwoodite.eosmat");
    if !path.is_file() {
        // The standalone Rust crate does not package the Python catalog.
        return;
    }
    let material = load_eosmat(path).unwrap();
    let record = material.default_record().unwrap();
    assert_eq!(
        record.identifier,
        "ferrous_ringwoodite_finger_1979_bm2_ambient_cell_reconstruction"
    );
    for (volume, pressure) in [(558.7, 0.0), (553.5, 1.867_379_040_453_365)] {
        assert!((record.pressure(volume, 296.15).unwrap() - pressure).abs() < 1e-10);
    }
    for pressure in [0.0, 2.0, 4.0] {
        let volume = record.volume(pressure, 296.15).unwrap();
        assert!((record.pressure(volume, 296.15).unwrap() - pressure).abs() < 1e-9);
    }
    // Low-level Rust isothermal evaluators ignore temperature; the record's
    // source reference and validity metadata must survive loading unchanged.
    assert_eq!(record.reference_temperature, 296.15);
    let raw = &material.document["eos_records"][0];
    assert_eq!(
        material.document["eos_records"].as_array().unwrap().len(),
        1
    );
    assert_eq!(raw["record_kind"], "derived");
    assert_eq!(raw["temperature_ref"], 296.15);
    assert_eq!(
        raw["validity"]["pressure_gpa"],
        serde_json::json!([0.0, 4.0])
    );
    assert_eq!(
        raw["validity"]["temperature_k"],
        serde_json::json!([296.15, 296.15])
    );
    assert!(raw.get("thermal").is_none());
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
