use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn diopside_published_records_load_evaluate_invert_and_preserve_metadata() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/diopside.eosmat");
    if !path.is_file() {
        // The Python resource tree is not included in the standalone Rust crate.
        return;
    }
    let material = load_eosmat(path).unwrap();
    material.validate().unwrap();
    assert_eq!(material.eos_records.len(), 2);
    assert_eq!(
        material.default_record().unwrap().identifier,
        "diopside_thompson_2008_bm3"
    );
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
    let static_record = material.record("diopside_thompson_2008_bm3").unwrap();
    assert!((static_record.pressure(407.16, 300.0).unwrap() - 10.16).abs() < 0.04);
    let thermal = material.record("diopside_zhao_1998_thermal_bm3").unwrap();
    assert!((thermal.pressure(421.43, 1272.0).unwrap() - 8.277_209_555_111_074).abs() < 1e-10);
    assert!((thermal.pressure(426.65, 1280.0).unwrap() - 6.88).abs() < 0.12);
    assert!(thermal.pressure(438.67, 300.0).unwrap().abs() < 1e-12);
    for temperature in [300.0, 800.0, 1280.0] {
        for pressure in [0.0, 4.0, 8.18] {
            let volume = thermal.volume(pressure, temperature).unwrap();
            assert!((thermal.pressure(volume, temperature).unwrap() - pressure).abs() < 1e-9);
        }
    }
    for pressure in [0.0001, 5.0, 10.16] {
        let volume = static_record.volume(pressure, 300.0).unwrap();
        assert!((static_record.pressure(volume, 300.0).unwrap() - pressure).abs() < 1e-9);
    }
    assert_eq!(
        thermal.document["scientific_validation"]["numerical_reproduction"]
            ["author_fit_reproduced"],
        false
    );
    assert_eq!(
        thermal.document["scientific_validation"]["parity_assessment"]
            ["parameter_parity_with_documented_table_error"],
        true
    );
}
