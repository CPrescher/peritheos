// Exact reference identities and zero-offset equivalence are deliberate checks.
#![allow(clippy::float_cmp)]

use peritheos::isothermal::{BM2Reference, ReferenceEnergyEos, BM2};
use peritheos::thermal::ReferenceStateEos;
use peritheos::IsothermalEos;

#[test]
fn finite_pressure_reference_and_independent_pressure_case() {
    let eos = BM2Reference::new(375.04, 254.0, 32.7).unwrap();
    assert_eq!(eos.pressure(375.04).unwrap(), 32.7);
    assert_eq!(eos.bulk_modulus(375.04).unwrap(), 254.0);
    // Vref/V = 8 => eta = 2, P-Pref = 3*Kref*(128-32)/2.
    assert!((eos.pressure(375.04 / 8.0).unwrap() - 36608.7).abs() < 1e-8);
    for p in [20.0, 32.7, 40.0, 70.4] {
        let v = eos.volume(p).unwrap();
        assert!((eos.pressure(v).unwrap() - p).abs() < 1e-7);
        let step = v * 1e-5;
        let numerical_bulk_modulus =
            -v * (eos.pressure(v + step).unwrap() - eos.pressure(v - step).unwrap()) / (2.0 * step);
        assert!((eos.bulk_modulus(v).unwrap() - numerical_bulk_modulus).abs() < 1e-6);
    }
    let h = 375.04 * 1e-5;
    let derivative = (eos.bulk_modulus(375.04 + h).unwrap()
        - eos.bulk_modulus(375.04 - h).unwrap())
        / (eos.pressure(375.04 + h).unwrap() - eos.pressure(375.04 - h).unwrap());
    assert!((derivative - 4.0).abs() < 1e-7);
    let shifted = eos.with_reference_state(360.0, 270.0).unwrap();
    assert_eq!(shifted.pressure(360.0).unwrap(), 32.7);
}

#[test]
fn reference_energy_integrates_absolute_pressure_and_zero_offset_is_bm2() {
    let eos = BM2Reference::new(2.0, 254.0, 32.7).unwrap();
    let v = 1.8;
    let h = 1e-6;
    let pressure_from_energy = -(eos.reference_energy(v + h).unwrap()
        - eos.reference_energy(v - h).unwrap())
        / (2.0 * h * 1e4);
    assert!((pressure_from_energy - eos.pressure(v).unwrap()).abs() < 1e-7);
    assert_eq!(eos.reference_energy(2.0).unwrap(), 0.0);
    let ambient = BM2Reference::new(2.0, 254.0, 0.0).unwrap();
    let legacy = BM2::new(2.0, 254.0).unwrap();
    for v in [1.5, 1.8, 2.0, 2.1] {
        assert_eq!(ambient.pressure(v).unwrap(), legacy.pressure(v).unwrap());
        assert_eq!(
            ambient.bulk_modulus(v).unwrap(),
            legacy.bulk_modulus(v).unwrap()
        );
        assert_eq!(
            ambient.reference_energy(v).unwrap(),
            legacy.reference_energy(v).unwrap()
        );
    }
    assert!(BM2Reference::new(0.0, 254.0, 32.7).is_err());
    assert!(BM2Reference::new(2.0, -1.0, 32.7).is_err());
    assert!(BM2Reference::new(2.0, 254.0, f64::NAN).is_err());
    assert!(eos.pressure(-1.0).is_err());
}
