use peritheos::isothermal::BM3;
use peritheos::thermal::DebyeQuadraticThermalPressure;
use peritheos::ThermalEos;

#[test]
fn fei_material_loads_and_round_trips_on_cell_volume_basis() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/iron.eosmat");
    if !path.exists() {
        return;
    }
    let material = peritheos::load_eosmat(path).unwrap();
    let record = material
        .eos_records
        .iter()
        .find(|r| r.identifier == "iron_fei_2016_bm3_debye_quadratic")
        .unwrap();
    let p = record.pressure(18.764_694, 1208.0).unwrap();
    assert!((p - 53.872_809_446_406_59).abs() < 1e-8);
    let text = material.to_json().unwrap();
    let restored = peritheos::load_eosmat_str(&text).unwrap();
    assert_eq!(material.document, restored.document);
}

#[test]
fn fei_2016_equation_and_source_state() {
    let factor = 6.022_140_76e23 * 1e-25 / 2.0;
    let reference = BM3::new(22.427_668_953_702_995 * factor, 172.7, 4.79).unwrap();
    let model = DebyeQuadraticThermalPressure::new(
        reference,
        300.0,
        422.0,
        1.74,
        0.78,
        1.0,
        5.788_65e-7,
        0.34,
    )
    .unwrap();
    let volume = 18.764_694 * factor;
    // Independent SI-unit Debye quadrature, and Table S2 row 3 (53.96 GPa).
    let pressure = model.pressure(volume, 1208.0).unwrap();
    assert!((pressure - 53.872_809_446_406_59).abs() < 1e-8);
    assert!((pressure - 53.96).abs() < 0.1);
    assert!(model.thermal_pressure(volume, 300.0).unwrap().abs() < 1e-12);
    assert!((model.volume(pressure, 1208.0).unwrap() - volume).abs() < 1e-9);
    for (v, t) in [(0.0, 300.0), (volume, 0.0), (f64::NAN, 300.0)] {
        assert!(model.pressure(v, t).is_err());
    }
    assert!(DebyeQuadraticThermalPressure::new(
        reference,
        300.0,
        422.0,
        1.74,
        0.78,
        1.0,
        f64::NAN,
        0.34,
    )
    .is_err());
}
