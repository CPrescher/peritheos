use peritheos::isothermal::{BM2, BM3};
use peritheos::thermal::{ReferenceVolumeLaw, ThermalExpansionLaw, ThermalReferenceState};
use peritheos::{load_eosmat_str, IsothermalEos, ThermalEos};
use serde_json::json;

fn model() -> ThermalReferenceState<BM3> {
    ThermalReferenceState::new(
        BM3::new(40.917, 41.8, 6.6).unwrap(),
        300.0,
        6.69e-5,
        -0.0319,
        0.0,
        ThermalExpansionLaw::Constant,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap()
    .with_linear_pressure_derivative(0.00571)
    .unwrap()
}

#[test]
fn linear_kprime_matches_direct_equation_and_bulk_derivative() {
    let eos = model();
    for t in [300.0_f64, 473.0, 673.0, 873.0] {
        let delta = t - 300.0;
        let v0 = 40.917 * (6.69e-5 * delta).exp();
        let k0 = 41.8 - 0.0319 * delta;
        let kp = 6.6 + 0.00571 * delta;
        let v = 34.0;
        let eta = (v0 / v).powf(1.0 / 3.0);
        let expected = 1.5
            * k0
            * (eta.powi(7) - eta.powi(5))
            * (1.0 + 0.75 * (kp - 4.0) * (eta.powi(2) - 1.0));
        assert!((eos.pressure(v, t).unwrap() - expected).abs() < 1e-11);
        assert!(eos.pressure(v0, t).unwrap().abs() < 1e-12);
        assert!((eos.bulk_modulus(v0, t, 1e-6).unwrap() - k0).abs() < 1e-10);
        let h = 1e-4;
        let derivative =
            -v * (eos.pressure(v + h, t).unwrap() - eos.pressure(v - h, t).unwrap()) / (2.0 * h);
        assert!((eos.bulk_modulus(v, t, 1e-6).unwrap() - derivative).abs() < 1e-7);
        let p = eos.pressure(v, t).unwrap();
        assert!((eos.volume(p, t).unwrap() - v).abs() < 1e-9);
    }
    assert!((eos.pressure(34.0, 873.0).unwrap() - 14.051_454_515_344_124).abs() < 1e-10);
    assert!(eos.thermal_pressure(34.0, 300.0).unwrap().abs() < 1e-12);
}

#[test]
fn linear_kprime_composes_with_logarithmic_and_reciprocal_modulus_laws() {
    let eos = model().with_temperature_laws(None, 0.0003).unwrap();
    let t = 700.0_f64;
    let delta = t - 300.0;
    let expected = BM3::new(
        40.917 * (6.69e-5 * delta).exp(),
        41.8 - 0.0319 * delta,
        6.6 + 0.00571 * delta + 0.0003 * delta * (t / 300.0).ln(),
    )
    .unwrap();
    assert!((eos.pressure(35.0, t).unwrap() - expected.pressure(35.0).unwrap()).abs() < 1e-11);

    let cubic = ThermalReferenceState::new(
        BM3::new(40.917, 41.8, 6.6).unwrap(),
        300.0,
        6.69e-5,
        0.0,
        0.0,
        ThermalExpansionLaw::Constant,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap()
    .with_temperature_laws(Some([1e-6, 0.0, 0.0]), 0.0)
    .unwrap()
    .with_linear_pressure_derivative(0.00571)
    .unwrap();
    let expected = BM3::new(
        40.917 * (6.69e-5 * delta).exp(),
        1.0 / (1.0 / 41.8 + 1e-6 * delta),
        6.6 + 0.00571 * delta,
    )
    .unwrap();
    assert!((cubic.pressure(35.0, t).unwrap() - expected.pressure(35.0).unwrap()).abs() < 1e-11);
}

#[test]
fn linear_kprime_rejects_invalid_parameters_and_states() {
    let eos = model();
    for slope in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(eos.with_linear_pressure_derivative(slope).is_err());
    }
    for (v, t) in [(0.0, 300.0), (34.0, 0.0), (f64::NAN, 300.0), (34.0, 2000.0)] {
        assert!(eos.pressure(v, t).is_err());
    }
    let bm2 = ThermalReferenceState::new(
        BM2::new(40.917, 41.8).unwrap(),
        300.0,
        6.69e-5,
        -0.0319,
        0.0,
        ThermalExpansionLaw::Constant,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap();
    assert!(bm2.with_linear_pressure_derivative(0.001).is_err());
    assert!(bm2.with_linear_pressure_derivative(0.0).is_ok());
}

#[test]
fn linear_kprime_eosmat_load_roundtrip_and_wrong_model_rejection() {
    let doc = json!({
        "format": "peritheos.material", "format_version": 3,
        "identifier": "conditional_brucite", "name": "Conditional equation test",
        "formula": "Mg(OH)2",
        "units": {"pressure": "GPa", "temperature": "K", "volume": "angstrom^3/conventional_unit_cell"},
        "eos_records": [{
            "identifier": "conditional_brucite_linear_kprime", "label": "Conditional 300 K anchored audit surface",
            "record_kind": "diagnostic", "equation_kind": "thermal", "default": true,
            "parameter_errors": {}, "fixed_parameters": [],
            "reference": "Fukui (2003) slopes; intercepts assumed from 300 K slice, not author-workflow reproduction",
            "eos": {"type": "BM3", "model": "birch_murnaghan_3", "parameters": {"V0": 40.917, "K0": 41.8, "K0_prime": 6.6}},
            "thermal": {"type": "AlphaKT", "model": "thermal_reference_state", "parameters": {
                "Tr": 300.0, "alpha0": 6.69e-5, "dK_dT": -0.0319, "dKprime_dT": 0.00571
            }},
            "scientific_validation": {"status": "pending_primary_source_check"}
        }]
    });
    let material = load_eosmat_str(&doc.to_string()).unwrap();
    let record = material
        .record("conditional_brucite_linear_kprime")
        .unwrap();
    assert!((record.pressure(34.0, 873.0).unwrap() - 14.051_454_515_344_124).abs() < 1e-10);
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(material.document, restored.document);
    let mut wrong = doc.clone();
    wrong["eos_records"][0]["thermal"]["type"] = json!("LinearThermalPressure");
    wrong["eos_records"][0]["thermal"]["model"] = json!("linear_thermal_pressure");
    assert!(load_eosmat_str(&wrong.to_string()).is_err());
    // Explicit null must not silently default to zero.
    let mut missing = doc;
    missing["eos_records"][0]["thermal"]["parameters"]["dKprime_dT"] = json!(null);
    assert!(load_eosmat_str(&missing.to_string()).is_err());
}
