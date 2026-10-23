use peritheos::{load_eosmat, serialize_eosmat};
use std::path::Path;

#[test]
fn source_table_iv_ag_preserves_fixed_volume_scale_and_execution() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dewaele_2008_ag.eosmat");
    let material = load_eosmat(path).unwrap();
    let record = material.record("silver_dewaele_2008_table4_vinet").unwrap();
    let volume: f64 = 4.0 * 12.3348;
    let x = (12.3348_f64 / 17.021).cbrt();
    let expected = 3.0 * 103.9 * (1.0 - x) / (x * x) * (1.5 * (5.78 - 1.0) * (1.0 - x)).exp();
    let pressure = record.pressure(volume, 300.0).unwrap();
    assert!((pressure - expected).abs() < 1e-10);
    assert!((record.volume(pressure, 300.0).unwrap() - volume).abs() < 1e-8);
    let encoded = serialize_eosmat(&material.document).unwrap();
    let reloaded = peritheos::load_eosmat_str(&encoded).unwrap();
    assert_eq!(reloaded.document, material.document);
    let raw = &material.document["eos_records"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["identifier"] == "silver_dewaele_2008_table4_vinet")
        .unwrap();
    assert_eq!(raw["parameters"]["V0"]["fixed"], true);
    assert_eq!(
        raw["pressure_calibration"]["methods"][0]["reference_calibration_record"],
        "ruby_dewaele_2004"
    );
}
