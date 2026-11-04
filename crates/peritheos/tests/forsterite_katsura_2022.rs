use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn forsterite_final_table_loads_evaluates_inverts_and_preserves_source_metadata() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/forsterite.eosmat");
    let material = load_eosmat(path).unwrap();
    assert_eq!(
        material.eos_records[0].identifier,
        "forsterite_finkelstein_2014_bm3_1"
    );
    let record = material.record("forsterite_katsura_2022_bm3_mgd").unwrap();
    // Independent translation of the deposited author's polynomial and constants
    // with the final Table 2 theta0=768 K, V/V0=.95, T=1700 K.
    let expected = 12.195_391_775_358_03;
    let pressure = record.pressure(289.58 * 0.95, 1700.0).unwrap();
    assert!((pressure - expected).abs() < 7e-5);
    for temperature in [310.0, 1300.0, 1900.0] {
        for pressure in [0.0, 5.0, 13.5] {
            let volume = record.volume(pressure, temperature).unwrap();
            assert!((record.pressure(volume, temperature).unwrap() - pressure).abs() < 1e-8);
        }
    }
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(material.document, restored.document);
    assert_eq!(
        record.document["thermal"]["parameters"]["theta0"]["value"],
        768.0
    );
    assert_eq!(record.document["thermal"]["parameters"]["n"]["value"], 7.0);
}

#[test]
fn forsterite_author_script_preserves_the_full_precision_coefficient_vector() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/forsterite.eosmat");
    let material = load_eosmat(path).unwrap();
    let record = material
        .record("forsterite_katsura_2022_author_script_bm3_mgd")
        .unwrap();
    assert_eq!(
        record.document["thermal"]["parameters"]["theta0"]["value"],
        760.0
    );
    assert_eq!(
        record.document["parameters"]["K0_prime"]["value"],
        4.211_168_969_885_74
    );
    // Literal geotherm-script polynomial, R and NA at V/V0=.95, T=1700 K.
    let pressure = record.pressure(289.58 * 0.95, 1700.0).unwrap();
    assert!((pressure - 12.180_116_341_180_696).abs() < 7e-5);
    let volume = record.volume(13.5, 1900.0).unwrap();
    assert!((record.pressure(volume, 1900.0).unwrap() - 13.5).abs() < 1e-8);
}
