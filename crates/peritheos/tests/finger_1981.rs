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
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/argon_fcc.eosmat");
    if !path.is_file() {
        return;
    }
    let material = load_eosmat(path).unwrap();
    let record = material
        .record("argon_fcc_finger_1981_murnaghan2_debye")
        .unwrap();
    // Independent adaptive Debye quadrature, four-atom cell units.
    assert!((record.pressure(110.0, 293.0).unwrap() - 3.327_019_022_150_685_7).abs() < 1e-9);
    assert!((record.volume(3.327_019_022_150_685_7, 293.0).unwrap() - 110.0).abs() < 1e-7);
}
