use peritheos::isothermal::{Holzapfel, ModifiedTait, Vinet, BM3};
use peritheos::thermal::{
    debye_function_3, AsymptoticPowerLawMieGruneisenDebye, DebyeTemperatureLaw,
    DoubleDebyeHelmholtz, LinearThermalPressure, LogVolumeThermalPressure, MieGruneisenDebye,
    MieGruneisenEinstein, MultiOscillatorGruneisen, ReferenceVolumeLaw, Sokolova2016,
    SokolovaParameters, ThermalExpansionLaw, ThermalModifiedTait, ThermalReferenceState,
    GAS_CONSTANT,
};
use peritheos::{CaloricEos, EosError, ThermalEos};
use serde_json::Value;

fn assert_close(actual: f64, expected: f64, relative_tolerance: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= relative_tolerance * scale,
        "actual {actual:.17e} differs from expected {expected:.17e}"
    );
}

#[test]
fn double_debye_helmholtz_matches_the_benedict_diamond_regression() {
    let model = DoubleDebyeHelmholtz::new(
        Vinet::new(0.343_466_776_105_840_03, 432.4, 3.793).unwrap(),
        0.335_493_461_739_6,
        1887.8,
        -5.247_303_452_269_356,
        0.913,
        1887.8,
        2.789_705_632_852_062_4,
        0.429,
        1887.8,
        1.404_816_050_829_074_1,
        0.499,
        1.0,
        3.79e-5,
        0.348_380_842_966,
        0.0,
        -874_736.021_029_928_6,
    )
    .unwrap();
    for (atomic_volume, temperature, expected_pressure) in [
        (5.7034, 300.0, 5.097_381_463_860_107),
        (5.4, 1000.0, 34.724_104_556_157_07),
        (5.0, 2000.0, 87.831_067_537_911_82),
        (4.654_270_411_587_497, 3000.0, 150.0),
    ] {
        let volume = atomic_volume * 0.060_221_407_6;
        assert_close(
            model.pressure(volume, temperature).unwrap(),
            expected_pressure,
            5.0e-10,
        );
    }

    let volume = 4.654_270_411_587_497 * 0.060_221_407_6;
    assert_close(model.volume(150.0, 3000.0).unwrap(), volume, 1.0e-10);
    assert_close(model.temperature(150.0, volume).unwrap(), 3000.0, 1.0e-10);
    assert!(model.zero_point_energy(volume).unwrap() > 0.0);
    assert!(model.molar_heat_capacity_v(volume, 3000.0).unwrap() > 0.0);
    assert!(model
        .helmholtz_free_energy(volume, 3000.0)
        .unwrap()
        .is_finite());
}

#[test]
fn shared_mie_gruneisen_literature_cases_match() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let debye = MieGruneisenDebye::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
    assert_close(
        debye.thermal_pressure(0.9, 1800.0).unwrap(),
        10.430_895_000_629_15,
        1.0e-10,
    );

    let einstein = MieGruneisenEinstein::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
    assert_close(
        einstein.thermal_pressure(0.9, 1800.0).unwrap(),
        9.961_068_596_665_895,
        1.0e-10,
    );
}

#[test]
fn debye_temperature_laws_are_explicit_and_distinct() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let integrated = MieGruneisenDebye::new(reference, 300.0, 800.0, 1.5, 1.2, 2.0).unwrap();
    let variable = MieGruneisenDebye::new_with_temperature_law(
        reference,
        300.0,
        800.0,
        1.5,
        1.2,
        2.0,
        DebyeTemperatureLaw::VariableExponent,
    )
    .unwrap();
    let volume = 0.8;
    let ratio = volume / reference.v0;
    let gamma = 1.5 * ratio.powf(1.2);
    assert_close(
        variable.characteristic_temperature(volume).unwrap(),
        800.0 * ratio.powf(-gamma),
        1.0e-14,
    );
    assert!(
        (variable.characteristic_temperature(volume).unwrap()
            - integrated.characteristic_temperature(volume).unwrap())
        .abs()
            > 1.0
    );
}

#[test]
fn simple_thermal_pressure_models_preserve_reference_state() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let linear = LinearThermalPressure::new(reference, 300.0, 0.002).unwrap();
    let logarithmic = LogVolumeThermalPressure::new(reference, 300.0, 0.002, -0.000_01).unwrap();

    for model_pressure in [
        linear.thermal_pressure(0.8, 300.0).unwrap(),
        logarithmic.thermal_pressure(0.8, 300.0).unwrap(),
    ] {
        assert_close(model_pressure, 0.0, 1.0e-15);
    }
    assert_close(linear.thermal_pressure(0.8, 800.0).unwrap(), 1.0, 1.0e-14);
    assert_close(
        logarithmic.thermal_pressure(0.8, 800.0).unwrap(),
        (0.002 - 0.000_01 * (1.0_f64 / 0.8).ln()) * 500.0,
        1.0e-14,
    );
}

#[test]
fn thermal_reference_state_supports_both_volume_laws_and_domains() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let integrated = ThermalReferenceState::new(
        reference,
        300.0,
        2.0e-5,
        -0.01,
        1.0e-8,
        ThermalExpansionLaw::LinearTemperature,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap();
    let linear = ThermalReferenceState::new(
        reference,
        300.0,
        2.0e-5,
        -0.01,
        0.0,
        ThermalExpansionLaw::Constant,
        ReferenceVolumeLaw::LinearTemperature,
    )
    .unwrap();

    assert_close(
        integrated.thermal_pressure(0.9, 300.0).unwrap(),
        0.0,
        1.0e-14,
    );
    assert_close(linear.thermal_pressure(0.9, 300.0).unwrap(), 0.0, 1.0e-14);
    assert!(integrated.pressure(0.9, 1200.0).unwrap().is_finite());
    assert!(linear.pressure(0.9, 1200.0).unwrap().is_finite());

    let invalid_modulus = ThermalReferenceState::new(
        reference,
        300.0,
        0.0,
        -1.0,
        0.0,
        ThermalExpansionLaw::Constant,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap();
    assert!(invalid_modulus.pressure(0.9, 500.0).is_err());
}

#[test]
fn thermal_configuration_and_state_errors_are_typed() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    assert!(matches!(
        ThermalReferenceState::new(
            reference,
            300.0,
            2.0e-5,
            -0.01,
            1.0e-8,
            ThermalExpansionLaw::Constant,
            ReferenceVolumeLaw::IntegratedExpansivity,
        ),
        Err(EosError::InvalidParameter { name: "alpha1", .. })
    ));
    assert!(matches!(
        ThermalReferenceState::new(
            reference,
            300.0,
            2.0e-5,
            -0.01,
            1.0e-8,
            ThermalExpansionLaw::LinearTemperature,
            ReferenceVolumeLaw::LinearTemperature,
        ),
        Err(EosError::InvalidParameter {
            name: "reference_volume_law",
            ..
        })
    ));
    assert!(matches!(
        AsymptoticPowerLawMieGruneisenDebye::new(reference, 300.0, 760.0, 1.5, 1.1, 2.5, 2.0,),
        Err(EosError::InvalidParameter { name: "a", .. })
    ));

    let double_debye = DoubleDebyeHelmholtz::new(
        Vinet::new(1.0, 160.0, 4.0).unwrap(),
        1.0,
        1000.0,
        0.0,
        1.0,
        1000.0,
        0.0,
        1.0,
        1000.0,
        0.0,
        1.0,
        1.0,
        -1.0,
        1.0,
        0.0,
        0.0,
    );
    assert!(matches!(
        double_debye,
        Err(EosError::InvalidParameter { name: "alpha0", .. })
    ));

    let model = DoubleDebyeHelmholtz::new(
        Vinet::new(1.0, 160.0, 4.0).unwrap(),
        1.0,
        1000.0,
        0.0,
        1.0,
        1000.0,
        0.0,
        1.0,
        1000.0,
        0.0,
        1.0,
        1.0,
        0.0,
        1.0,
        0.0,
        0.0,
    )
    .unwrap();
    assert!(matches!(
        model.dac_thermal_pressure(1.0, 300.0, 1.0),
        Err(EosError::InvalidState { name: "f_dac", .. })
    ));
    assert!(matches!(
        model.ion_helmholtz_free_energy(1.0, -1.0),
        Err(EosError::InvalidState {
            name: "temperature",
            ..
        })
    ));
    assert!(matches!(
        debye_function_3(0.0),
        Err(EosError::InvalidState { .. })
    ));
}

#[test]
fn asymptotic_power_law_model_round_trips_and_preserves_reference_state() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let model =
        AsymptoticPowerLawMieGruneisenDebye::new(reference, 300.0, 760.0, 1.5, 0.3, 2.5, 2.0)
            .unwrap();
    let volume = 0.83;
    let temperature = 1800.0;
    let pressure = model.pressure(volume, temperature).unwrap();

    assert_close(model.thermal_pressure(volume, 300.0).unwrap(), 0.0, 1.0e-14);
    assert_close(
        model.volume(pressure, temperature).unwrap(),
        volume,
        1.0e-10,
    );
    assert_close(
        model.temperature(pressure, volume).unwrap(),
        temperature,
        1.0e-9,
    );
}

#[test]
fn multi_oscillator_model_accepts_a_generic_reference_isotherm() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let parameters = SokolovaParameters::reduced(
        298.15, 684.0, 0.564, 1561.0, 2.436, -0.506, 1.085, 0.0, 0.0, 0.0, 0.0,
    );
    let model = MultiOscillatorGruneisen::new_with_atom_count(reference, parameters, 2.0).unwrap();

    assert_close(model.thermal_pressure(0.9, 298.15).unwrap(), 0.0, 1.0e-14);
    assert!(model.pressure(0.9, 1800.0).unwrap().is_finite());
}

#[test]
fn shared_holland_powell_literature_case_matches() {
    let reference = ModifiedTait::new(1.0, 160.0, 4.0, -0.01).unwrap();
    let model = ThermalModifiedTait::new(reference, 298.15, 700.0, 2.5e-5, 2.0).unwrap();
    assert_close(
        model.thermal_pressure(0.9, 1200.0).unwrap(),
        5.029_535_604_766_256,
        1.0e-13,
    );
}

#[test]
fn debye_function_matches_known_limits_and_reference_values() {
    assert_close(
        debye_function_3(1.0e-5).unwrap(),
        0.999_996_250_005,
        1.0e-14,
    );
    assert_close(
        debye_function_3(1.0).unwrap(),
        0.674_415_564_077_814_7,
        1.0e-13,
    );
    assert_close(
        debye_function_3(10.0).unwrap(),
        0.019_295_765_690_345_49,
        1.0e-13,
    );
    assert!(debye_function_3(200.0).unwrap().is_finite());
}

#[test]
fn thermal_models_round_trip_volume_and_temperature() {
    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let debye = MieGruneisenDebye::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
    let volume = 0.85;
    let temperature = 1800.0;
    let pressure = debye.pressure(volume, temperature).unwrap();
    assert_close(
        debye.volume(pressure, temperature).unwrap(),
        volume,
        1.0e-10,
    );
    assert_close(
        debye.temperature(pressure, volume).unwrap(),
        temperature,
        1.0e-9,
    );
}

#[test]
fn mie_gruneisen_caloric_identities_are_preserved() {
    let model = MieGruneisenEinstein::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();
    let volume = 0.9;
    let temperature = 1200.0;
    let energy = model.thermal_energy(volume, temperature).unwrap();
    let entropy = model.thermal_entropy(volume, temperature).unwrap();
    assert_close(
        model
            .thermal_helmholtz_free_energy(volume, temperature)
            .unwrap(),
        energy - temperature * entropy,
        1.0e-13,
    );
    let cv = model.molar_heat_capacity_v(volume, temperature).unwrap();
    let cp = model.molar_heat_capacity_p(volume, temperature).unwrap();
    assert!(cv > 0.0);
    assert!(cp > cv);
    assert_close(
        model.gruneisen_parameter(volume, temperature).unwrap(),
        model.volume_gruneisen_parameter(volume).unwrap(),
        1.0e-14,
    );
}

#[test]
fn high_temperature_debye_energy_reaches_dulong_petit_limit() {
    let model = MieGruneisenDebye::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        1.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();
    let temperature = 10_000.0;
    assert_close(
        model.thermal_energy(1.0, temperature).unwrap(),
        3.0 * model.n * GAS_CONSTANT * temperature,
        5.0e-5,
    );
}

fn reduced_sokolova() -> Sokolova2016 {
    let reference = Holzapfel::new(0.3414, 441.5, 3.9, 1.0, 6.0).unwrap();
    let parameters = SokolovaParameters::reduced(
        298.15, 684.0, 0.564, 1561.0, 2.436, -0.506, 1.085, 0.0, 0.0, 0.0, 0.0,
    );
    Sokolova2016::new(reference, parameters).unwrap()
}

#[test]
fn sokolova_reduced_pressure_and_integral_regressions_match() {
    let model = reduced_sokolova();
    assert_close(
        model
            .thermal_pressure(0.9 * model.rt_eos.v0, 3000.0)
            .unwrap(),
        14.860_881_778_954_129,
        1.0e-9,
    );
    assert_close(
        model.gruneisen_integral(0.9).unwrap().exp(),
        1.094_276_545_296_065_2,
        1.0e-9,
    );
    assert_close(
        model.gruneisen_integral(1.1).unwrap().exp(),
        0.910_574_594_040_097_6,
        1.0e-9,
    );
}

#[test]
fn sokolova_complete_pressure_terms_match_python_regressions() {
    let reference = Holzapfel::new(0.3414, 441.5, 3.9, 1.0, 6.0).unwrap();
    let mut parameters = SokolovaParameters::reduced(
        298.15, 684.0, 0.564, 1561.0, 2.436, -0.506, 1.085, 5.2, 1.3, 0.8, 2.7,
    );
    parameters.beta = 0.35;
    parameters.qbo = 480.0;
    parameters.d = 2.4;
    parameters.mb = 0.75;
    parameters.qb1o = 1120.0;
    parameters.d1 = 1.6;
    parameters.mb1 = 0.4;
    let model = Sokolova2016::new(reference, parameters).unwrap();

    for (volume_ratio, temperature, expected) in [
        (0.72, 900.0, 3.097_970_564_162_762_4),
        (0.88, 2400.0, 17.277_471_289_674_384),
        (1.05, 700.0, 2.622_165_296_289_313),
    ] {
        assert_close(
            model
                .thermal_pressure(volume_ratio * model.rt_eos.v0, temperature)
                .unwrap(),
            expected,
            1.0e-9,
        );
    }
}

#[test]
fn sokolova_round_trips_total_pressure_in_both_state_variables() {
    let model = reduced_sokolova();
    let volume = 0.8 * model.rt_eos.v0;
    let temperature = 2500.0;
    let pressure = model.pressure(volume, temperature).unwrap();
    assert_close(
        model.volume(pressure, temperature).unwrap(),
        volume,
        1.0e-10,
    );
    assert_close(
        model.temperature(pressure, volume).unwrap(),
        temperature,
        1.0e-9,
    );
}

fn fixture_number(state: &Value, name: &str) -> f64 {
    state["quantities"][name]
        .as_f64()
        .unwrap_or_else(|| panic!("missing fixture quantity {name}"))
}

fn assert_common_thermal<T: ThermalEos>(
    model: &T,
    state: &Value,
    tolerance: f64,
    derivative_tolerance: f64,
) {
    let volume = state["volume"].as_f64().unwrap();
    let temperature = state["temperature"].as_f64().unwrap();
    assert_close(
        model.thermal_pressure(volume, temperature).unwrap(),
        fixture_number(state, "thermal_pressure"),
        tolerance,
    );
    assert_close(
        model.pressure(volume, temperature).unwrap(),
        fixture_number(state, "pressure"),
        tolerance,
    );
    assert_close(
        model.bulk_modulus(volume, temperature, 1.0e-6).unwrap(),
        fixture_number(state, "bulk_modulus"),
        derivative_tolerance,
    );
    assert_close(
        model
            .isothermal_compressibility(volume, temperature)
            .unwrap(),
        fixture_number(state, "isothermal_compressibility"),
        derivative_tolerance,
    );
    assert_close(
        model
            .thermal_expansivity(volume, temperature, 1.0e-5)
            .unwrap(),
        fixture_number(state, "thermal_expansivity"),
        derivative_tolerance,
    );
}

fn assert_common_caloric<T: CaloricEos>(model: &T, state: &Value, tolerance: f64) {
    let volume = state["volume"].as_f64().unwrap();
    let temperature = state["temperature"].as_f64().unwrap();
    assert_close(
        model.molar_heat_capacity_v(volume, temperature).unwrap(),
        fixture_number(state, "molar_heat_capacity_v"),
        tolerance,
    );
    assert_close(
        model.molar_heat_capacity_p(volume, temperature).unwrap(),
        fixture_number(state, "molar_heat_capacity_p"),
        tolerance,
    );
    assert_close(
        model.gruneisen_parameter(volume, temperature).unwrap(),
        fixture_number(state, "gruneisen_parameter"),
        tolerance,
    );
    assert_close(
        model.adiabatic_bulk_modulus(volume, temperature).unwrap(),
        fixture_number(state, "adiabatic_bulk_modulus"),
        tolerance,
    );
}

fn assert_mie<const DEBYE: bool>(
    model: &peritheos::thermal::MieGruneisen<BM3, DEBYE>,
    state: &Value,
    tolerance: f64,
) {
    assert_common_thermal(model, state, tolerance, tolerance);
    assert_common_caloric(model, state, tolerance);
    let volume = state["volume"].as_f64().unwrap();
    let temperature = state["temperature"].as_f64().unwrap();
    for (name, value) in [
        (
            "characteristic_temperature",
            model.characteristic_temperature(volume).unwrap(),
        ),
        (
            "thermal_energy",
            model.thermal_energy(volume, temperature).unwrap(),
        ),
        (
            "thermal_entropy",
            model.thermal_entropy(volume, temperature).unwrap(),
        ),
        (
            "vibrational_pressure",
            model.vibrational_pressure(volume, temperature).unwrap(),
        ),
        (
            "thermal_helmholtz_free_energy",
            model
                .thermal_helmholtz_free_energy(volume, temperature)
                .unwrap(),
        ),
        (
            "thermal_enthalpy",
            model.thermal_enthalpy(volume, temperature).unwrap(),
        ),
        (
            "thermal_gibbs_free_energy",
            model
                .thermal_gibbs_free_energy(volume, temperature)
                .unwrap(),
        ),
    ] {
        assert_close(value, fixture_number(state, name), tolerance);
    }
}

#[test]
fn shared_python_thermal_compatibility_fixture_matches() {
    let raw = include_str!("data/thermal_compatibility_cases.json");
    let document: Value = serde_json::from_str(raw).expect("valid compatibility fixture");
    assert_eq!(document["schema_version"].as_u64(), Some(1));
    let tolerance = document["relative_tolerance"].as_f64().unwrap();
    let sokolova_derivative_tolerance = document["sokolova_derivative_relative_tolerance"]
        .as_f64()
        .unwrap();
    for case in document["debye_function_3"].as_array().unwrap() {
        assert_close(
            debye_function_3(case["argument"].as_f64().unwrap()).unwrap(),
            case["value"].as_f64().unwrap(),
            tolerance,
        );
    }

    let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
    let debye = MieGruneisenDebye::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
    let einstein = MieGruneisenEinstein::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
    let tait = ThermalModifiedTait::new(
        ModifiedTait::new(1.0, 160.0, 4.0, -0.01).unwrap(),
        298.15,
        700.0,
        2.5e-5,
        2.0,
    )
    .unwrap();
    let reduced = reduced_sokolova();
    let mut complete_parameters = reduced.parameters;
    complete_parameters.a_0 = 5.2;
    complete_parameters.m = 1.3;
    complete_parameters.g = 0.8;
    complete_parameters.e_0 = 2.7;
    complete_parameters.beta = 0.35;
    complete_parameters.qbo = 480.0;
    complete_parameters.d = 2.4;
    complete_parameters.mb = 0.75;
    complete_parameters.qb1o = 1120.0;
    complete_parameters.d1 = 1.6;
    complete_parameters.mb1 = 0.4;
    let complete = Sokolova2016::new(reduced.rt_eos, complete_parameters).unwrap();

    for case in document["cases"].as_array().unwrap() {
        for state in case["states"].as_array().unwrap() {
            match case["model"].as_str().unwrap() {
                "MieGruneisenDebye" => assert_mie(&debye, state, tolerance),
                "MieGruneisenEinstein" => assert_mie(&einstein, state, tolerance),
                "ThermalModifiedTait" => {
                    assert_common_thermal(&tait, state, tolerance, tolerance);
                    assert_common_caloric(&tait, state, tolerance);
                }
                "Sokolova2016Reduced" => {
                    assert_common_thermal(
                        &reduced,
                        state,
                        tolerance,
                        sokolova_derivative_tolerance,
                    );
                }
                "Sokolova2016Complete" => {
                    assert_common_thermal(
                        &complete,
                        state,
                        tolerance,
                        sokolova_derivative_tolerance,
                    );
                }
                model => panic!("unexpected thermal fixture model {model}"),
            }
        }
    }
}
