use peritheos::hugoniot::{Hugoniot, LinearUsUpHugoniot};
use peritheos::EosError;

fn assert_close(actual: f64, expected: f64, relative_tolerance: f64) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= relative_tolerance * scale,
        "actual {actual:.17e} differs from expected {expected:.17e}"
    );
}

#[test]
fn reference_state_and_jump_conditions_are_exact() {
    let model = LinearUsUpHugoniot::new(10.0, 8.0, 4.0, 1.5, 0.0).unwrap();
    assert_close(model.pressure(10.0).unwrap(), 0.0, f64::EPSILON);
    assert_close(model.particle_velocity(10.0).unwrap(), 0.0, f64::EPSILON);
    assert_close(model.shock_velocity(10.0).unwrap(), 4.0, f64::EPSILON);
    assert_close(model.density(10.0).unwrap(), 8.0, f64::EPSILON);
    assert_close(
        model.specific_internal_energy_change(10.0).unwrap(),
        0.0,
        f64::EPSILON,
    );

    let volume = 8.0;
    let pressure = model.pressure(volume).unwrap();
    let up = model.particle_velocity(volume).unwrap();
    let us = model.shock_velocity(volume).unwrap();
    assert_close(pressure, model.rho0 * us * up, 1.0e-14);
    assert_close(volume / model.v0, 1.0 - up / us, 1.0e-14);
    assert_close(us, model.c0 + model.s * up, 1.0e-14);
}

#[test]
fn pressure_volume_inverse_is_stable() {
    let model = LinearUsUpHugoniot::new(10.0, 8.0, 4.0, 1.5, 0.2).unwrap();
    for fraction in [1.0, 0.95, 0.8, 0.5] {
        let volume = fraction * model.v0;
        let pressure = model.pressure(volume).unwrap();
        assert_close(model.volume(pressure).unwrap(), volume, 1.0e-13);
    }
}

#[test]
fn energy_and_tangent_modulus_match_independent_identities() {
    let model = LinearUsUpHugoniot::new(10.0, 8.0, 4.0, 1.5, 0.2).unwrap();
    let volume = 8.0;
    let density = model.density(volume).unwrap();
    let pressure = model.pressure(volume).unwrap();
    let expected_energy = 0.5 * (pressure + model.p0) * (1.0 / model.rho0 - 1.0 / density);
    assert_close(
        model.specific_internal_energy_change(volume).unwrap(),
        expected_energy,
        1.0e-14,
    );

    let step = 1.0e-6 * volume;
    let derivative = (model.pressure(volume + step).unwrap()
        - model.pressure(volume - step).unwrap())
        / (2.0 * step);
    assert_close(
        model.tangent_modulus(volume).unwrap(),
        -volume * derivative,
        1.0e-9,
    );
}

#[test]
fn direct_particle_velocity_relations_are_consistent() {
    let model = LinearUsUpHugoniot::new(10.0, 8.0, 4.0, 1.5, 0.0).unwrap();
    let up = 1.0;
    let us = 5.5;
    assert_close(
        model.shock_velocity_from_particle_velocity(up).unwrap(),
        us,
        f64::EPSILON,
    );
    assert_close(
        model.pressure_from_particle_velocity(up).unwrap(),
        44.0,
        f64::EPSILON,
    );
    let volume = model.volume_from_particle_velocity(up).unwrap();
    assert_close(model.particle_velocity(volume).unwrap(), up, 1.0e-14);
}

#[test]
fn invalid_parameters_and_states_are_rejected() {
    assert!(matches!(
        LinearUsUpHugoniot::new(10.0, 8.0, 0.0, 1.5, 0.0),
        Err(EosError::InvalidParameter { name: "c0", .. })
    ));
    let model = LinearUsUpHugoniot::new(10.0, 8.0, 4.0, 1.5, 0.0).unwrap();
    assert!(matches!(
        model.pressure(10.1),
        Err(EosError::InvalidState { name: "volume", .. })
    ));
    assert_eq!(model.volume(-1.0), Err(EosError::OutsideInvertibleRange));
    assert_eq!(model.pressure(3.0), Err(EosError::OutsideInvertibleRange));
}
