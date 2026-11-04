use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn farla_b2_halides_load_preserve_metadata_and_execute_in_rust() {
    for (compound, source_ratio) in [("rbcl", 0.8149), ("rbbr", 0.7886), ("rbi", 0.7558)] {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{compound}.eosmat"));
        let material = load_eosmat(path).unwrap();
        let identifier = format!("{compound}_b2_farla_2025_bm3_mgd_cscl");
        let record = material.record(&identifier).unwrap();
        let reference_volume = record.document["parameters"]["V0"]["value"]
            .as_f64()
            .unwrap();
        let volume = record.volume(10.0, 1500.0).unwrap();
        // S5-S7: coefficient rounding can exceed the last printed ratio digit.
        assert!((volume / reference_volume - source_ratio).abs() < 0.00025);
        assert!((record.pressure(volume, 1500.0).unwrap() - 10.0).abs() < 1e-9);
        assert!(record.pressure(reference_volume, 300.0).unwrap().abs() < 1e-12);
        assert_eq!(
            record.document["thermal"]["debye_temperature_law"],
            "integrated_gruneisen"
        );
        let round_trip = load_eosmat_str(&material.to_json().unwrap()).unwrap();
        assert_eq!(round_trip.document, material.document);
    }
}
