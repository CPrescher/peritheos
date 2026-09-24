use peritheos::isothermal::DensityPolynomial3;
use peritheos::{load_eosmat_str, IsothermalEos};

#[test]
fn density_polynomial_reproduces_source_bulk_modulus_and_roundtrip() {
    let eos =
        DensityPolynomial3::new(132.670_429_310_921_66, 2.0, 12.65, -11.43, 1.5, 0.68).unwrap();
    assert!((eos.pressure(eos.v0).unwrap() - 1.23).abs() < 1e-12);
    for rho in [2.017, 2.5, 3.0, 4.0, 4.8] {
        let v = eos.v0 * eos.rho0 / rho;
        let expected = 12.65 - 11.43 * rho + 1.5 * rho * rho + 0.68 * rho * rho * rho;
        assert!((eos.pressure(v).unwrap() - expected).abs() < 1e-12);
        let k = rho * (-11.43 + 3.0 * rho + 2.04 * rho * rho);
        assert!((eos.bulk_modulus(v).unwrap() - k).abs() < 1e-12);
        assert!((eos.volume(expected).unwrap() - v).abs() < 1e-10);
    }
    assert!(eos.pressure(eos.v0 * 2.0).is_err());
    assert!(eos.volume(0.0).is_err());
    assert!(DensityPolynomial3::new(100.0, 1.0, 12.65, -11.43, 1.5, 0.68).is_err());
}

#[test]
fn native_material_loader_accepts_published_density_refit() {
    let material = load_eosmat_str(include_str!(
        "../../../peritheos/data/materials/argon_fcc.eosmat"
    ))
    .unwrap();
    let record = material
        .eos_records
        .iter()
        .find(|r| r.document["identifier"] == "argon_fcc_grimsditch_1986_density_polynomial")
        .unwrap();
    assert_eq!(
        record.document["reproduction"]["fit_status"],
        "not_reproduced"
    );
}

#[test]
fn expansion_roots_remain_on_stable_density_branch() {
    let eos = DensityPolynomial3::new(100.0, 2.0, 12.65, -11.43, 1.5, 0.68).unwrap();
    for density in [1.7436, 1.76, 1.8, 1.9] {
        let volume = 200.0 / density;
        let pressure = eos.pressure(volume).unwrap();
        assert!((eos.volume(pressure).unwrap() - volume).abs() < 1e-8);
    }
    assert!(eos.volume(0.885).is_err());
}
