//! Native loading retains published `Fe3C` coefficients and scoped evidence.

use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn published_fe3c_fits_preserve_adoption_parity_and_uncertainty_qualifications() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/cementite.eosmat");
    if !path.is_file() {
        return; // The standalone Rust crate does not bundle the Python catalog.
    }
    let material = load_eosmat(path).unwrap();
    {
        let (suffix, k0, kp, error_k0, pressure) = ("full_range_bm3", 161.0, 5.9, 2.0, 40.0);
        let record = material
            .record(&format!("cementite_prescher_2012_{suffix}"))
            .unwrap();
        let doc = &record.document;
        assert_eq!(doc["parameters"]["K0"]["value"], k0);
        assert_eq!(doc["parameters"]["K0"]["error"], error_k0);
        assert_eq!(doc["parameters"]["K0_prime"]["value"], kp);
        assert_eq!(doc["parameters"]["V0"]["value"], 154.567_886_370_481_2);
        assert_eq!(
            doc["evidence"]["parameter_provenance"]["V0"]["value_kind"],
            "adopted"
        );
        assert_eq!(
            doc["evidence"]["reproduction"]["assessment"]["status"],
            "parity"
        );
        assert_eq!(
            doc["evidence"]["reproduction"]["assessment"]["scope"],
            "parameter_agreement"
        );
        let volume = record.volume(pressure, 293.0).unwrap();
        assert!((record.pressure(volume, 293.0).unwrap() - pressure).abs() < 1e-9);
    }
    for suffix in ["below_22_gpa_bm3", "22_47_gpa_bm3"] {
        assert!(material
            .record(&format!("cementite_prescher_2012_{suffix}"))
            .is_none());
    }
    let full = material
        .record("cementite_prescher_2012_full_range_bm3")
        .unwrap();
    for parameter in full.document["evidence"]["reproduction"]["uncertainty_comparisons"]
        .as_array()
        .unwrap()
    {
        assert_eq!(parameter["finding"], "matched_at_printed_precision");
    }
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
