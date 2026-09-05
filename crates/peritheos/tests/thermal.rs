use peritheos::isothermal::{Holzapfel, ModifiedTait, Vinet, BM3};
use peritheos::thermal::{
    debye_function_3, AsymptoticPowerLawMieGruneisenDebye, DebyeTemperatureLaw,
    DorogokupetsOganov2007, DorogokupetsOganov2007Parameters, DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz, LinearThermalPressure, LogVolumeThermalPressure,
    MieGruneisenDebye, MieGruneisenEinstein, MultiOscillatorGruneisen, ReferenceVolumeLaw,
    SecondOrderTaylorThermalPressure, Sokolova2016, SokolovaParameters, ThermalExpansionLaw,
    ThermalModifiedTait, ThermalReferenceState, GAS_CONSTANT,
};
use peritheos::{CaloricEos, EosError, IsothermalEos, ThermalEos};
use serde_json::Value;

fn assert_close(actual: f64, expected: f64, relative_tolerance: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= relative_tolerance * scale,
        "actual {actual:.17e} differs from expected {expected:.17e}"
    );
}

#[test]
fn dorogokupets_oganov_pt_reproduces_table_vi_isochors() {
    let model = DorogokupetsOganov2007::new(
        Vinet::new(0.9091, 276.07, 5.30).unwrap(),
        DorogokupetsOganov2007Parameters {
            tr: 298.15,
            theta_b1: 95.2,
            d_b1: 8.199,
            m_b1: 0.329,
            theta_b2: 148.4,
            d_b2: 4.005,
            m_b2: 0.383,
            theta_e1: 214.6,
            m_e1: 1.211,
            theta_e2: 140.8,
            m_e2: 1.077,
            gamma0: 2.802,
            gamma_inf: 1.538,
            beta: 5.550,
            anharmonic_a: 160.9,
            anharmonic_m: 4.06,
            electronic_e: 260.0,
            electronic_g: 2.4,
            defect_h: 32572.0,
            defect_s: 0.631,
        },
        1.0,
    )
    .unwrap();
    let temperatures = [298.15, 1000.0, 2000.0, 3000.0];
    let rows = [
        (1.00, [0.0, 5.309, 12.864, 20.349]),
        (0.95, [16.207, 21.219, 28.485, 35.822]),
        (0.90, [38.307, 43.111, 50.199, 57.476]),
        (0.85, [68.389, 73.069, 80.079, 87.376]),
        (0.80, [109.388, 114.016, 121.043, 128.434]),
        (0.75, [165.473, 170.119, 177.253, 184.808]),
        (0.70, [242.676, 247.403, 254.730, 262.523]),
    ];
    for (ratio, expected_pressures) in rows {
        for (temperature, expected_pressure) in temperatures.into_iter().zip(expected_pressures) {
            let actual = model.pressure(ratio * 0.9091, temperature).unwrap();
            assert!(
                (actual - expected_pressure).abs() <= 0.02,
                "x={ratio}, T={temperature}: actual {actual}, expected {expected_pressure}"
            );
        }
    }
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
fn log_moment_double_debye_matches_the_correa_diamond_regression() {
    let model = DoubleDebyeLogMomentHelmholtz::new(
        Vinet::new(0.348_380_842_966, 368.2, 4.038).unwrap(),
        0.335_493_461_739_6,
        1887.8,
        -5.247_303_452_269_356,
        0.913,
        1887.8,
        2.789_705_632_852_062_4,
        0.429,
        1887.8,
        2.175_306_177_997_739,
        0.202,
        1.0,
        3.8e-5,
        -14_960_919.113_708_327,
    )
    .unwrap();
    for (atomic_volume, temperature, expected_pressure) in [
        (5.785, 300.0, 4.710_783_105_788_316),
        (5.571, 300.0, 20.018_830_109_151_665),
        (4.43, 5000.0, 202.628_115_197_741_86),
        (3.21, 9000.0, 729.946_307_877_533_7),
    ] {
        let volume = atomic_volume * 0.060_221_407_6;
        assert_close(
            model.pressure(volume, temperature).unwrap(),
            expected_pressure,
            5.0e-10,
        );
    }

    let (weight_a, weight_b) = model.double_debye_weights(4.0 * 0.060_221_407_6).unwrap();
    assert_close(weight_a + weight_b, 1.0, 1.0e-14);
    assert!(model
        .helmholtz_free_energy(4.0 * 0.060_221_407_6, 5000.0)
        .unwrap()
        .is_finite());
    assert!(
        model
            .molar_heat_capacity_v(4.0 * 0.060_221_407_6, 5000.0)
            .unwrap()
            > 0.0
    );
}

#[test]
fn double_debye_reference_temperatures_anchor_the_supplied_vinet_isotherm() {
    let reference = Vinet::new(5.6693 * 0.060_221_407_6, 444.5, 4.18).unwrap();
    let benedict_absolute = DoubleDebyeHelmholtz::new(
        reference,
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
        0.0,
    )
    .unwrap();
    let benedict = benedict_absolute.with_reference_temperature(298.0).unwrap();
    let correa_absolute = DoubleDebyeLogMomentHelmholtz::new(
        reference,
        0.335_493_461_739_6,
        1887.8,
        -5.247_303_452_269_356,
        0.913,
        1887.8,
        2.789_705_632_852_062_4,
        0.429,
        1887.8,
        2.175_306_177_997_739,
        0.202,
        1.0,
        3.8e-5,
        0.0,
    )
    .unwrap();
    let correa = correa_absolute.with_reference_temperature(298.0).unwrap();
    let volume = 5.0 * 0.060_221_407_6;
    let reference_pressure = reference.pressure(volume).unwrap();

    for (anchored_at_reference, anchored_hot, absolute_at_reference, absolute_hot) in [
        (
            benedict.pressure(volume, 298.0).unwrap(),
            benedict.pressure(volume, 5000.0).unwrap(),
            benedict_absolute.pressure(volume, 298.0).unwrap(),
            benedict_absolute.pressure(volume, 5000.0).unwrap(),
        ),
        (
            correa.pressure(volume, 298.0).unwrap(),
            correa.pressure(volume, 5000.0).unwrap(),
            correa_absolute.pressure(volume, 298.0).unwrap(),
            correa_absolute.pressure(volume, 5000.0).unwrap(),
        ),
    ] {
        assert_close(anchored_at_reference, reference_pressure, 1.0e-13);
        assert_close(
            anchored_hot,
            reference_pressure + absolute_hot - absolute_at_reference,
            1.0e-12,
        );
    }
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
fn dac_confined_volume_closes_forward_boundary_condition() {
    let model = MieGruneisenDebye::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();
    let cold_pressure = 40.0;
    let temperature = 2500.0;
    let f_dac = 0.25;

    let cold_volume = model.volume(cold_pressure, 300.0).unwrap();
    let isobaric_volume = model.volume(cold_pressure, temperature).unwrap();
    let confined_volume = model
        .volume_with_dac_confinement(cold_pressure, temperature, f_dac)
        .unwrap();
    let thermal_increment = model
        .thermal_pressure_increment(confined_volume, temperature)
        .unwrap();
    let confinement_pressure = model
        .dac_thermal_pressure(confined_volume, temperature, f_dac)
        .unwrap();
    let total_pressure = model.pressure(confined_volume, temperature).unwrap();

    assert!(isobaric_volume > confined_volume);
    assert!(confined_volume > cold_volume);
    assert_close(confinement_pressure, f_dac * thermal_increment, 1.0e-12);
    assert_close(
        total_pressure,
        cold_pressure + confinement_pressure,
        1.0e-11,
    );
    assert_close(
        model
            .temperature_from_volumes(cold_volume, confined_volume, f_dac)
            .unwrap(),
        temperature,
        1.0e-10,
    );
}

#[test]
fn dac_confined_volume_rejects_invalid_fraction() {
    let model =
        LinearThermalPressure::new(BM3::new(1.0, 160.0, 4.0).unwrap(), 300.0, 0.005).unwrap();

    assert!(matches!(
        model.volume_with_dac_confinement(40.0, 2000.0, 1.0),
        Err(EosError::InvalidState { name: "f_dac", .. })
    ));
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
fn second_order_taylor_thermal_pressure_preserves_its_absolute_baseline() {
    let reference = Vinet::new(74.074_102_512_3, 169.8, 4.501).unwrap();
    let model = SecondOrderTaylorThermalPressure::new(
        reference, 300.0, 0.02, 0.5096, -13.4246, 6.3295e-3, 36.2194, 5.4705e-8, 3.2238e-3,
    )
    .unwrap();
    let volume = reference.v0 * 0.8;
    let temperature = 3000.0;
    let delta_eta = 0.18;
    let delta_temperature = 2700.0;
    let expected = 0.5096 - 13.4246 * delta_eta
        + 6.3295e-3 * delta_temperature
        + 0.5 * 36.2194 * delta_eta * delta_eta
        + 0.5 * 5.4705e-8 * delta_temperature * delta_temperature
        + 0.5 * 3.2238e-3 * delta_eta * delta_temperature;

    assert_close(
        model.thermal_pressure(volume, temperature).unwrap(),
        expected,
        1.0e-14,
    );
    assert!(model.thermal_pressure(volume, 300.0).unwrap().abs() > 0.1);
    assert_close(
        model.thermal_pressure_increment(volume, 300.0).unwrap(),
        0.0,
        1.0e-14,
    );
    let pressure = model.pressure(volume, temperature).unwrap();
    assert_close(
        model.volume(pressure, temperature).unwrap(),
        volume,
        1.0e-10,
    );
    assert_close(
        model.temperature(pressure, volume).unwrap(),
        temperature,
        1.0e-10,
    );
}

#[test]
fn thermal_reference_state_supports_all_volume_laws_and_domains() {
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
    let berman = ThermalReferenceState::new(
        reference,
        298.0,
        1.94e-5,
        -0.008,
        5.73e-10,
        ThermalExpansionLaw::LinearTemperature,
        ReferenceVolumeLaw::Berman,
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
    assert!(berman.pressure(0.9, 1200.0).unwrap().is_finite());

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
fn thermal_reference_state_integrates_inverse_temperature_squared_expansivity() {
    let reference = BM3::new(1.0, 261.0, 4.0).unwrap();
    let model = ThermalReferenceState::new_with_inverse_square(
        reference,
        300.0,
        1.982e-5,
        -0.0280,
        0.818e-8,
        0.474,
        ThermalExpansionLaw::LinearTemperatureInverseSquare,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap();
    let temperature: f64 = 2000.0;
    let exponent = 1.982e-5 * (temperature - 300.0)
        + 0.5 * 0.818e-8 * (temperature.powi(2) - 300.0_f64.powi(2))
        + 0.474 * (temperature.recip() - 300.0_f64.recip());
    let shifted_volume = exponent.exp();

    assert_close(
        model.pressure(shifted_volume, temperature).unwrap(),
        0.0,
        1.0e-12,
    );
    assert_close(model.pressure(1.0, 300.0).unwrap(), 0.0, 1.0e-14);
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
