use peritheos::isothermal::{ModifiedTait, BM3};
use peritheos::thermal::{MieGruneisenDebye, ThermalModifiedTait};
use peritheos::{
    CaloricEos, CaloricEosBatch, IsothermalEos, IsothermalEosBatch, ThermalEos, ThermalEosBatch,
};

#[test]
fn isothermal_batches_match_ordered_scalar_evaluation_and_round_trip() {
    let model = BM3::new(10.0, 120.0, 4.3).unwrap();
    let volumes = [7.5, 8.5, 9.5, 10.0];
    let pressures = model.pressures(&volumes).unwrap();
    let expected = volumes
        .iter()
        .map(|&volume| model.pressure(volume).unwrap())
        .collect::<Vec<_>>();

    assert_eq!(pressures, expected);
    assert_eq!(model.bulk_moduli(&volumes).unwrap().len(), volumes.len());
    assert!(IsothermalEosBatch::pressures(&model, &[])
        .unwrap()
        .is_empty());
    for (recovered, expected) in model.volumes(&pressures).unwrap().iter().zip(volumes) {
        assert!((recovered - expected).abs() < 1.0e-9);
    }

    let error = IsothermalEosBatch::pressures(&model, &[10.0, 0.0, 9.0]).unwrap_err();
    assert!(error.to_string().contains("volume"));
}

#[test]
fn thermal_and_caloric_batches_match_scalar_methods() {
    let model = MieGruneisenDebye::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();
    let volumes = [0.75, 0.85, 0.95];
    let temperatures = [500.0, 1_500.0, 2_500.0];
    let pressures = model.pressures(&volumes, &temperatures).unwrap();

    for index in 0..volumes.len() {
        assert_eq!(
            pressures[index].to_bits(),
            model
                .pressure(volumes[index], temperatures[index])
                .unwrap()
                .to_bits()
        );
    }
    assert_eq!(
        model
            .molar_heat_capacities_v(&volumes, &temperatures)
            .unwrap(),
        volumes
            .iter()
            .zip(temperatures)
            .map(|(&volume, temperature)| {
                model.molar_heat_capacity_v(volume, temperature).unwrap()
            })
            .collect::<Vec<_>>()
    );
    let recovered = model.volumes(&pressures, &temperatures).unwrap();
    for (actual, expected) in recovered.iter().zip(volumes) {
        assert!((actual - expected).abs() < 1.0e-9);
    }
}

#[test]
fn remaining_thermal_and_caloric_batches_match_scalar_methods() {
    let model = MieGruneisenDebye::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();
    let volumes = [0.75, 0.85, 0.95];
    let temperatures = [500.0, 1_500.0, 2_500.0];
    let pressures = model.pressures(&volumes, &temperatures).unwrap();
    let thermal_pressures = model.thermal_pressures(&volumes, &temperatures).unwrap();
    let recovered_temperatures = model.temperatures(&pressures, &volumes).unwrap();
    let bulk_moduli = model.bulk_moduli(&volumes, &temperatures, 1.0e-6).unwrap();
    let compressibilities = model
        .isothermal_compressibilities(&volumes, &temperatures)
        .unwrap();
    let expansivities = model
        .thermal_expansivities(&volumes, &temperatures, 1.0e-5)
        .unwrap();
    let heat_capacities_p = model
        .molar_heat_capacities_p(&volumes, &temperatures)
        .unwrap();
    let gruneisen = model.gruneisen_parameters(&volumes, &temperatures).unwrap();
    let adiabatic_moduli = model
        .adiabatic_bulk_moduli(&volumes, &temperatures)
        .unwrap();

    for index in 0..volumes.len() {
        let volume = volumes[index];
        let temperature = temperatures[index];
        assert_eq!(
            thermal_pressures[index].to_bits(),
            model
                .thermal_pressure(volume, temperature)
                .unwrap()
                .to_bits()
        );
        assert!((recovered_temperatures[index] - temperature).abs() < 1.0e-8);
        assert_eq!(
            bulk_moduli[index].to_bits(),
            model
                .bulk_modulus(volume, temperature, 1.0e-6)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            compressibilities[index].to_bits(),
            model
                .isothermal_compressibility(volume, temperature)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            expansivities[index].to_bits(),
            model
                .thermal_expansivity(volume, temperature, 1.0e-5)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            heat_capacities_p[index].to_bits(),
            model
                .molar_heat_capacity_p(volume, temperature)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            gruneisen[index].to_bits(),
            model
                .gruneisen_parameter(volume, temperature)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            adiabatic_moduli[index].to_bits(),
            model
                .adiabatic_bulk_modulus(volume, temperature)
                .unwrap()
                .to_bits()
        );
    }
}

#[test]
fn temperature_from_volumes_batch_preserves_scalar_convention() {
    let model = ThermalModifiedTait::new(
        ModifiedTait::new(1.0, 160.0, 4.0, -0.01).unwrap(),
        298.15,
        700.0,
        2.5e-5,
        2.0,
    )
    .unwrap();
    let ambient = [0.98, 0.96];
    let heated = [0.99, 0.98];
    let batch = model
        .temperatures_from_volumes(&ambient, &heated, 0.2)
        .unwrap();

    for index in 0..ambient.len() {
        assert_eq!(
            batch[index].to_bits(),
            model
                .temperature_from_volumes(ambient[index], heated[index], 0.2)
                .unwrap()
                .to_bits()
        );
    }
}

#[test]
fn dac_confined_volume_batch_preserves_scalar_convention() {
    let model = ThermalModifiedTait::new(
        ModifiedTait::new(1.0, 160.0, 4.0, -0.01).unwrap(),
        298.15,
        700.0,
        2.5e-5,
        2.0,
    )
    .unwrap();
    let cold_pressures = [20.0, 40.0];
    let temperatures = [1000.0, 2000.0];
    let batch = model
        .volumes_with_dac_confinement(&cold_pressures, &temperatures, 0.2)
        .unwrap();
    let increments = model
        .thermal_pressure_increments(&batch, &temperatures)
        .unwrap();

    for index in 0..cold_pressures.len() {
        assert_eq!(
            batch[index].to_bits(),
            model
                .volume_with_dac_confinement(cold_pressures[index], temperatures[index], 0.2,)
                .unwrap()
                .to_bits()
        );
        assert_eq!(
            increments[index].to_bits(),
            model
                .thermal_pressure_increment(batch[index], temperatures[index])
                .unwrap()
                .to_bits()
        );
    }
}

#[test]
fn paired_batches_reject_mismatched_lengths_before_evaluation() {
    let model = MieGruneisenDebye::new(
        BM3::new(1.0, 160.0, 4.0).unwrap(),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    .unwrap();

    let error = model.pressures(&[0.9, 1.0], &[300.0]).unwrap_err();
    assert!(error.to_string().contains("matching lengths"));
}
