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
    for (recovered, expected) in model.volumes(&pressures).unwrap().iter().zip(volumes) {
        assert!((recovered - expected).abs() < 1.0e-9);
    }
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
