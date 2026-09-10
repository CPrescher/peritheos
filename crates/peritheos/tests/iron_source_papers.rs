use peritheos::isothermal::BM3;
use peritheos::thermal::Dewaele2006;

#[test]
fn yamazaki_zero_gamma_infinity_power_law_limit() {
    let reference = BM3::new(0.667, 202.0, 4.5).unwrap();
    let model = Dewaele2006::new(
        reference, 300.0, 1173.0, 3.2, 0.0, 0.8, 3.7e-5, 1.87, 1.95e-4, 1.339, 1.0,
    )
    .unwrap();
    let volume: f64 = 0.5;
    let expected_gamma = 3.2 * (volume / 0.667).powf(0.8);
    assert!((model.volume_gruneisen_parameter(volume).unwrap() - expected_gamma).abs() < 1e-12);
    assert!(Dewaele2006::new(
        reference, 300.0, 1173.0, 3.2, -0.1, 0.8, 3.7e-5, 1.87, 1.95e-4, 1.339, 1.0
    )
    .is_err());
}

#[test]
fn primary_iron_and_alloy_records_load_in_native_interchange() {
    for material in ["iron", "fe09ni01_hcp", "nacl_b2"] {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../peritheos/data/materials/{material}.eosmat"));
        if !path.exists() {
            continue;
        }
        let loaded = peritheos::load_eosmat(path).unwrap();
        let restored = peritheos::load_eosmat_str(&loaded.to_json().unwrap()).unwrap();
        assert_eq!(loaded.document, restored.document);
    }
}
