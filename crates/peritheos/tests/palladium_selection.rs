use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str, serialize_eosmat};

#[test]
fn palladium_default_hold_preserves_separate_explicit_source_records() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/palladium.eosmat");
    if !path.is_file() {
        return;
    }
    let material = load_eosmat(path).unwrap();
    assert!(material.default_record().is_none());
    assert!(material.default_equilibrium_record().is_none());
    assert!(material.eos_records.iter().all(|r| !r.is_default));
    let frost = material.record("palladium_frost_2023_vinet_1").unwrap();
    // Independent printed pressure from supplementary Table S3.
    assert!((frost.pressure(56.0, 300.0).unwrap() - 10.04).abs() < 0.005);
    for identifier in ["palladium_frost_2023_vinet_1", "palladium_frost_2023_bm3_2"] {
        let record = material.record(identifier).unwrap();
        assert_eq!(
            record.document["evidence"]["scientific_validation"]["status"],
            "not_reproduced"
        );
        assert_eq!(record.document["catalog_access"], "explicit_selection");
    }

    let dft = material.record("palladium_baty_2024_bm3_dft_2").unwrap();
    assert!(!dft.is_default);
    assert_eq!(dft.reference_temperature.to_bits(), 0.0_f64.to_bits());
    assert!((dft.reference_volume() - 58.4668).abs() < 1e-12);
    assert!(dft.pressure(dft.reference_volume(), 0.0).unwrap().abs() < 1e-10);

    let experimental = material.record("palladium_baty_2024_bm3_1").unwrap();
    assert!(!experimental.is_default);
    assert_eq!(
        experimental.document["evidence"]["scientific_validation"]["status"],
        "not_reproduced"
    );
    assert_eq!(
        experimental.document["evidence"]["scientific_validation"]["usage_recommendation"],
        "not_recommended_for_quantitative_use"
    );

    let constrained = material.record("palladium_fedotenko_2020_bm2_2").unwrap();
    assert_eq!(
        constrained.document["evidence"]["scientific_validation"]["status"],
        "primary_source_validated"
    );
    assert_eq!(
        constrained.document["evidence"]["reproduction"]["assessment"]["status"],
        "parity"
    );
    assert_eq!(
        constrained.document["evidence"]["scientific_validation"]["data_agreement_status"],
        "poor"
    );
    assert_eq!(constrained.document["catalog_access"], "explicit_selection");
    assert!((constrained.pressure(52.411, 300.0).unwrap() - 29.761_441_011_4).abs() < 1e-9);

    let round_trip = load_eosmat_str(&serialize_eosmat(&material.document).unwrap()).unwrap();
    assert!(round_trip.default_record().is_none());
    assert!(round_trip.default_equilibrium_record().is_none());
    assert_eq!(
        round_trip
            .record("palladium_fedotenko_2020_bm2_2")
            .unwrap()
            .document,
        constrained.document
    );
}
