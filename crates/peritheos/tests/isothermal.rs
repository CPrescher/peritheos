use std::collections::HashMap;

use peritheos::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
use peritheos::{EosError, IsothermalEos};
use serde::Deserialize;

fn assert_close(actual: f64, expected: f64, relative_tolerance: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= relative_tolerance * scale,
        "actual {actual:.17e} differs from expected {expected:.17e}"
    );
}

#[derive(Debug, Deserialize)]
struct LiteratureCase {
    id: String,
    model: String,
    parameters: HashMap<String, f64>,
    volume: f64,
    quantity: String,
    expected: f64,
}

fn parameter(parameters: &HashMap<String, f64>, context: &str, name: &str) -> f64 {
    *parameters
        .get(name)
        .unwrap_or_else(|| panic!("{context} is missing parameter {name}"))
}

fn model_from_parameters(model: &str, parameters: &HashMap<String, f64>) -> Box<dyn IsothermalEos> {
    let p = |name| parameter(parameters, model, name);
    match model {
        "BM2" => Box::new(BM2::new(p("V0"), p("K0")).unwrap()),
        "BM3" => Box::new(BM3::new(p("V0"), p("K0"), p("K0_prime")).unwrap()),
        "BM4" => Box::new(BM4::new(p("V0"), p("K0"), p("K0_prime"), p("K0_double_prime")).unwrap()),
        "Murnaghan" => Box::new(Murnaghan::new(p("V0"), p("K0"), p("K0_prime")).unwrap()),
        "ModifiedTait" => Box::new(
            ModifiedTait::new(p("V0"), p("K0"), p("K0_prime"), p("K0_double_prime")).unwrap(),
        ),
        "NaturalStrain2" => Box::new(NaturalStrain2::new(p("V0"), p("K0")).unwrap()),
        "NaturalStrain3" => Box::new(NaturalStrain3::new(p("V0"), p("K0"), p("K0_prime")).unwrap()),
        "NaturalStrain4" => Box::new(
            NaturalStrain4::new(p("V0"), p("K0"), p("K0_prime"), p("K0_double_prime")).unwrap(),
        ),
        "Vinet" => Box::new(Vinet::new(p("V0"), p("K0"), p("K0_prime")).unwrap()),
        "Holzapfel" => {
            Box::new(Holzapfel::new(p("V0"), p("K0"), p("K0_prime"), p("n"), p("Z")).unwrap())
        }
        _ => panic!("unknown isothermal model {model}"),
    }
}

#[test]
fn shared_literature_pressure_cases_match() {
    let raw = include_str!("data/literature_reference_cases.json");
    let cases: Vec<LiteratureCase> = serde_json::from_str(raw).expect("valid reference fixture");

    for case in cases {
        let supported = matches!(
            case.model.as_str(),
            "BM2"
                | "BM3"
                | "BM4"
                | "Murnaghan"
                | "NaturalStrain3"
                | "ModifiedTait"
                | "Vinet"
                | "Holzapfel"
        );
        if !supported {
            continue;
        }
        let pressure = model_from_parameters(&case.model, &case.parameters)
            .pressure(case.volume)
            .unwrap_or_else(|error| panic!("{} failed: {error}", case.id));

        assert_eq!(case.quantity, "pressure", "unexpected case {}", case.id);
        assert_close(pressure, case.expected, 1.0e-13);
    }
}

#[derive(Debug, Deserialize)]
struct CompatibilityFixture {
    schema_version: u32,
    source_commit: String,
    relative_tolerance: f64,
    cases: Vec<CompatibilityCase>,
}

#[derive(Debug, Deserialize)]
struct CompatibilityCase {
    model: String,
    parameters: HashMap<String, f64>,
    volume_fractions: Vec<f64>,
    pressures: Vec<f64>,
    bulk_moduli: Vec<f64>,
}

#[test]
fn python_compatibility_grid_matches_pressure_and_bulk_modulus() {
    let raw = include_str!("data/isothermal_compatibility_cases.json");
    let fixture: CompatibilityFixture =
        serde_json::from_str(raw).expect("valid compatibility fixture");
    assert_eq!(fixture.schema_version, 1);
    assert_eq!(
        fixture.source_commit,
        "12d033378418cfd6c9ece6050c550fc748ffe02a"
    );

    for case in fixture.cases {
        let model = model_from_parameters(&case.model, &case.parameters);
        assert_eq!(case.volume_fractions.len(), case.pressures.len());
        assert_eq!(case.volume_fractions.len(), case.bulk_moduli.len());
        for ((fraction, expected_pressure), expected_modulus) in case
            .volume_fractions
            .iter()
            .zip(&case.pressures)
            .zip(&case.bulk_moduli)
        {
            let volume = fraction * model.reference_volume();
            assert_close(
                model.pressure(volume).unwrap(),
                *expected_pressure,
                fixture.relative_tolerance,
            );
            assert_close(
                model.bulk_modulus(volume).unwrap(),
                *expected_modulus,
                fixture.relative_tolerance,
            );
        }
    }
}

fn representative_models() -> Vec<Box<dyn IsothermalEos>> {
    vec![
        Box::new(BM2::new(10.0, 120.0).unwrap()),
        Box::new(BM3::new(10.0, 120.0, 4.3).unwrap()),
        Box::new(BM4::new(10.0, 120.0, 4.3, -0.02).unwrap()),
        Box::new(Murnaghan::new(10.0, 120.0, 4.3).unwrap()),
        Box::new(ModifiedTait::new(10.0, 120.0, 4.3, -0.02).unwrap()),
        Box::new(NaturalStrain2::new(10.0, 120.0).unwrap()),
        Box::new(NaturalStrain3::new(10.0, 120.0, 4.3).unwrap()),
        Box::new(NaturalStrain4::new(10.0, 120.0, 4.3, -0.02).unwrap()),
        Box::new(Vinet::new(10.0, 120.0, 4.3).unwrap()),
        Box::new(Holzapfel::new(0.3414, 441.5, 3.9, 1.0, 6.0).unwrap()),
    ]
}

#[test]
fn all_isothermal_models_match_reference_state_and_round_trip_compression() {
    for model in representative_models() {
        let v0 = model.reference_volume();
        assert_close(model.pressure(v0).unwrap(), 0.0, 1.0e-12);
        for fraction in [0.55, 0.7, 0.85, 1.0] {
            let expected_volume = fraction * v0;
            let pressure = model.pressure(expected_volume).unwrap();
            assert!(pressure.is_finite());
            assert!(model.bulk_modulus(expected_volume).unwrap().is_finite());
            assert_close(model.volume(pressure).unwrap(), expected_volume, 1.0e-10);
        }
    }
}

#[test]
fn analytical_bulk_modulus_matches_pressure_derivative() {
    for model in representative_models() {
        let volume = 0.8 * model.reference_volume();
        let step = 1.0e-6 * volume;
        let derivative = (model.pressure(volume + step).unwrap()
            - model.pressure(volume - step).unwrap())
            / (2.0 * step);
        let numerical_bulk_modulus = -volume * derivative;
        assert_close(
            model.bulk_modulus(volume).unwrap(),
            numerical_bulk_modulus,
            2.0e-9,
        );
    }
}

#[test]
fn validation_and_expansion_branch_failures_are_typed() {
    assert!(matches!(
        BM2::new(0.0, 100.0),
        Err(EosError::InvalidParameter { name: "V0", .. })
    ));
    let model = BM2::new(10.0, 100.0).unwrap();
    assert!(matches!(
        model.pressure(0.0),
        Err(EosError::InvalidState { name: "volume", .. })
    ));
    assert_eq!(model.volume(-100.0), Err(EosError::OutsideInvertibleRange));
}

#[test]
fn murnaghan_zero_derivative_limit_is_supported() {
    let model = Murnaghan::new(10.0, 100.0, 0.0).unwrap();
    let volume = 8.0;
    assert_close(
        model.pressure(volume).unwrap(),
        100.0 * (10.0_f64 / volume).ln(),
        1.0e-14,
    );
    assert_close(model.bulk_modulus(volume).unwrap(), 100.0, 1.0e-14);
}
