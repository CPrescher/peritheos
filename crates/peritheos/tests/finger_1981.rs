use peritheos::isothermal::SecondOrderMurnaghan;
use peritheos::{load_eosmat, IsothermalEos};
use std::path::Path;

#[test]
fn finger_native_material_and_quadratic_bulk_modulus() {
    let eos = SecondOrderMurnaghan::new(22.557, 2.3701, 6.97, -0.4, -0.10289).unwrap();
    assert!((eos.pressure(22.557).unwrap() + 0.10289).abs() < 1e-14);
    let v = 17.0;
    let step = 1e-5;
    let numerical =
        -v * (eos.pressure(v + step).unwrap() - eos.pressure(v - step).unwrap()) / (2.0 * step);
    assert!((eos.bulk_modulus(v).unwrap() - numerical).abs() < 1e-7);
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/argon_fcc.eosmat");
    let material = load_eosmat(path).unwrap();
    let record = material
        .record("argon_fcc_finger_1981_murnaghan2_debye")
        .unwrap();
    // Independent adaptive Debye quadrature, four-atom cell units.
    assert!((record.pressure(110.0, 293.0).unwrap() - 3.327_019_022_150_685_7).abs() < 1e-9);
    assert!((record.volume(3.327_019_022_150_685_7, 293.0).unwrap() - 110.0).abs() < 1e-7);
}

#[test]
fn zero_point_configuration_cannot_be_silently_ignored_by_other_models() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/argon_fcc.eosmat");
    let mut document: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let records = document["eos_records"].as_array_mut().unwrap();
    let xiao = records
        .iter_mut()
        .find(|r| r["identifier"] == "argon_fcc_xiao_2025_helmholtz")
        .unwrap();
    xiao["thermal"]["zero_point_pressure"] = "included".into();
    xiao["thermal"]["thermal_pressure_reference"] = "absolute_zero".into();
    assert!(peritheos::load_eosmat_str(&document.to_string()).is_err());
}

#[test]
fn second_order_murnaghan_inverse_does_not_cross_compression_pole() {
    let eos = SecondOrderMurnaghan::new(10.0, 2.0, 4.0, 1.0, 0.2).unwrap();
    for volume in [5.0, 7.0, 9.0, 10.0, 12.0] {
        assert!((eos.volume(eos.pressure(volume).unwrap()).unwrap() - volume).abs() < 1e-12);
    }
    assert!(eos.volume(-1.0).is_err());
}
