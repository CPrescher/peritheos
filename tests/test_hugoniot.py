"""Tests for shock Hugoniot equations of state."""

import numpy as np
import pytest

from peritheos.errors import EosValidationError
from peritheos.hugoniot import LinearUsUpHugoniot


@pytest.fixture
def hugoniot():
    return LinearUsUpHugoniot(V0=10.0, rho0=8.0, c0=4.0, s=1.5)


def test_parameters_and_reconstruction(hugoniot):
    assert hugoniot.parameter_values() == {
        "V0": 10.0,
        "rho0": 8.0,
        "c0": 4.0,
        "s": 1.5,
        "P0": 0.0,
    }
    updated = hugoniot.with_parameters(s=1.6)
    assert updated.s == 1.6
    assert hugoniot.s == 1.5
    with pytest.raises(EosValidationError, match="Unknown parameters"):
        hugoniot.with_parameters(K0=100.0)


def test_scalar_array_and_round_trip(hugoniot):
    volumes = np.array([10.0, 9.0, 8.0, 7.0])
    pressures = hugoniot.pressure(volumes)
    assert pressures.shape == volumes.shape
    assert pressures[0] == 0.0
    assert np.all(np.diff(pressures) > 0.0)
    assert np.allclose(hugoniot.volume(pressures), volumes)


def test_rankine_hugoniot_identities(hugoniot):
    volume = 8.0
    pressure = hugoniot.pressure(volume)
    up = hugoniot.particle_velocity(volume)
    us = hugoniot.shock_velocity(volume)
    density = hugoniot.density(volume)

    assert np.isclose(pressure - hugoniot.P0, hugoniot.rho0 * us * up)
    assert np.isclose(volume / hugoniot.V0, 1.0 - up / us)
    assert np.isclose(us, hugoniot.c0 + hugoniot.s * up)
    assert np.isclose(
        hugoniot.specific_internal_energy_change(volume),
        0.5 * (pressure + hugoniot.P0) * (1.0 / hugoniot.rho0 - 1.0 / density),
    )


def test_state_collects_all_quantities(hugoniot):
    state = hugoniot.state(np.array([9.0, 8.0]))
    assert np.allclose(state.pressure, hugoniot.pressure(state.volume))
    assert np.allclose(state.shock_velocity, hugoniot.shock_velocity(state.volume))
    scalar = hugoniot.state(8.0)
    assert scalar.volume == 8.0


def test_state_from_particle_velocity_is_self_consistent(hugoniot):
    state = hugoniot.state_from_particle_velocity(1.0)
    assert state.particle_velocity == pytest.approx(1.0)
    assert state.shock_velocity == pytest.approx(5.5)
    assert state.pressure == pytest.approx(44.0)
    assert hugoniot.pressure(state.volume) == pytest.approx(state.pressure)


def test_direct_particle_velocity_and_tangent_relations(hugoniot):
    up = 1.0
    us = hugoniot.shock_velocity_from_particle_velocity(up)
    pressure = hugoniot.pressure_from_particle_velocity(up)
    volume = hugoniot.volume_from_particle_velocity(up)
    assert us == pytest.approx(hugoniot.c0 + hugoniot.s * up)
    assert pressure == pytest.approx(hugoniot.rho0 * us * up)
    assert hugoniot.pressure(volume) == pytest.approx(pressure)
    assert hugoniot.tangent_modulus(volume) > 0.0


def test_invalid_states(hugoniot):
    with pytest.raises(EosValidationError):
        hugoniot.pressure(11.0)
    with pytest.raises(EosValidationError):
        hugoniot.volume(-1.0)
    with pytest.raises(EosValidationError):
        hugoniot.pressure(3.0)
