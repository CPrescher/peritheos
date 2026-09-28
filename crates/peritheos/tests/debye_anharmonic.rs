use peritheos::isothermal::NaturalStrain4;
use peritheos::thermal::DebyeAnharmonicHelmholtz;
use peritheos::{CaloricEos, ThermalEos};

fn model() -> DebyeAnharmonicHelmholtz<NaturalStrain4> {
    let k0 = 2.6565;
    let kp = 2.0 + 2.0 * 7.298 / k0;
    let kpp = (6.0 * 0.010 / k0 - kp * kp + 3.0 * kp - 3.0) / k0;
    let cold = NaturalStrain4::new(2.2555, k0, kp, kpp).unwrap();
    DebyeAnharmonicHelmholtz::new(cold, 70.0, 86.44, 2.68, 0.0024, 0.0128, 0.388, 7.85).unwrap()
}

#[test]
fn authors_workbook_and_thermodynamic_derivatives() {
    let eos = model();
    let volume = 2.3;
    let temperature = 70.0;
    let pressure = eos.pressure(volume, temperature).unwrap();
    let helmholtz = eos.helmholtz_free_energy(volume, temperature).unwrap();
    let energy = eos.internal_energy(volume, temperature).unwrap();
    let entropy = eos.entropy(volume, temperature).unwrap();
    let cv = eos.molar_heat_capacity_v(volume, temperature).unwrap();
    // Independently cached cells from the authors' supplementary Excel workbook.
    assert!((pressure - 0.078_351_737_584_362_9).abs() < 1e-7);
    assert!((helmholtz - -998.229_660_630_414_4).abs() < 1e-8);
    assert!((energy - 1_098.200_890_676_024_3).abs() < 1e-8);
    assert!((entropy - 29.949_007_875_806_267).abs() < 1e-10);
    assert!((cv - 22.854_917_560_882_143).abs() < 1e-10);
    assert!(
        (eos.bulk_modulus(volume, temperature, 1e-6).unwrap()
            - 0.001 / 0.000_444_851_552_096_565_1)
            .abs()
            < 2e-6
    );
    let dv = 1e-5;
    let dt = 1e-3;
    let dp = -(eos.helmholtz_free_energy(volume + dv, temperature).unwrap()
        - eos.helmholtz_free_energy(volume - dv, temperature).unwrap())
        / (2.0 * dv * 1e4);
    let ds = -(eos.helmholtz_free_energy(volume, temperature + dt).unwrap()
        - eos.helmholtz_free_energy(volume, temperature - dt).unwrap())
        / (2.0 * dt);
    let du = (eos.internal_energy(volume, temperature + dt).unwrap()
        - eos.internal_energy(volume, temperature - dt).unwrap())
        / (2.0 * dt);
    assert!((dp - pressure).abs() < 1e-8);
    assert!((ds - entropy).abs() < 1e-7);
    assert!((du - cv).abs() < 1e-7);
    assert!((eos.volume(pressure, temperature).unwrap() - volume).abs() < 1e-8);
    assert!(eos.thermal_pressure_increment(volume, 70.0).unwrap().abs() < 1e-14);
    assert!(eos.thermal_pressure(volume, 70.0).unwrap() > 0.0);
}

#[test]
fn invalid_states_and_parameters() {
    let eos = model();
    assert!(eos.pressure(-1.0, 70.0).is_err());
    assert!(eos.pressure(2.3, 0.0).is_err());
    assert!(eos.entropy(2.3, f64::NAN).is_err());
    assert!(DebyeAnharmonicHelmholtz::new(
        eos.rt_eos, 70.0, 86.44, 2.68, 0.0024, 0.0128, -1.0, 7.85
    )
    .is_err());
}

#[test]
fn native_eosmat_roundtrip_converts_cell_volume() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/argon_fcc.eosmat");
    let material = peritheos::load_eosmat(path).unwrap();
    let record = material.record("argon_fcc_xiao_2025_helmholtz").unwrap();
    // Four monatomic formula units per conventional fcc cell.
    let cell_volume = 23.0 * 4.0 / 0.602_214_076;
    let pressure = record.pressure(cell_volume, 70.0).unwrap();
    assert!((pressure - 0.078_351_737_584_362_9).abs() < 1e-10);
    assert!((record.volume(pressure, 70.0).unwrap() - cell_volume).abs() < 1e-8);
    let exported = peritheos::serialize_eosmat(&material.document).unwrap();
    let reloaded = peritheos::load_eosmat_str(&exported).unwrap();
    assert!(
        (reloaded
            .record("argon_fcc_xiao_2025_helmholtz")
            .unwrap()
            .pressure(cell_volume, 70.0)
            .unwrap()
            - pressure)
            .abs()
            < 1e-12
    );
}
