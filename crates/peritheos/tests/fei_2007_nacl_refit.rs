use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn nacl_thermal_refit_matches_independent_quadrature_and_round_trips() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../peritheos/data/materials/nacl_b2.eosmat");
    if !path.is_file() {
        return; // The material bundle is not included in the standalone Rust crate.
    }
    let material = load_eosmat(path).unwrap();
    let refit_id = "nacl_b2_fei_2007_staged_thermal_refit";
    let record = material.record(refit_id).unwrap();
    assert_eq!(record.document["record_kind"], "refit");
    assert_eq!(record.document["derived_from_record"], "nacl_b2_fei_2007");
    assert_eq!(record.document["default"], false);
    assert!(record.document["parameter_covariance"].is_null());
    // Independent SI quadrature using the two-stage figure-data coefficients.
    for (volume, temperature, pressure) in [
        (27.0, 300.0, 33.004_773_420_467_62),
        (27.0, 1000.0, 35.312_252_856_742_79),
        (22.0, 1000.0, 80.323_604_484_309_83),
    ] {
        assert!((record.pressure(volume, temperature).unwrap() - pressure).abs() < 1e-7);
        assert!((record.volume(pressure, temperature).unwrap() - volume).abs() < 1e-8);
    }
    assert_eq!(
        material.record("nacl_b2_fei_2007").unwrap().document["thermal"]["parameters"]["q"],
        0.5
    );
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
