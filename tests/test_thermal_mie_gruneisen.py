"""Tests for Mie-Gruneisen thermal equations of state."""

import numpy as np
import pytest
from scipy.constants import R

from peritheos.eos import EosBase
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import (
    AsymptoticPowerLawMieGruneisenDebyeExcess,
    MieGruneisenDebye,
    MieGruneisenEinstein,
)
from peritheos.eos.thermal.mie_gruneisen import _debye_function_3


@pytest.fixture(params=[MieGruneisenDebye, MieGruneisenEinstein])
def eos(request):
    rt_eos = BM3(V0=1.0, K0=160.0, K0_prime=4.0)
    return request.param(
        rt_eos=rt_eos,
        Tr=300.0,
        theta0=800.0,
        gamma0=1.5,
        q=1.0,
        n=2.0,
    )


def test_parameters(eos):
    assert eos.Tr == 300.0
    assert eos.theta0 == 800.0
    assert eos.gamma0 == 1.5
    assert eos.q == 1.0
    assert eos.n == 2.0


def test_thermal_pressure_is_zero_at_reference_temperature(eos):
    volumes = eos.rt_eos.V0 * np.array([0.8, 1.0, 1.2])

    assert np.all(eos.thermal_pressure(volumes, eos.Tr) == 0.0)


def test_thermal_pressure_sign(eos):
    assert eos.thermal_pressure(eos.rt_eos.V0, 2000.0) > 0.0
    assert eos.thermal_pressure(eos.rt_eos.V0, 200.0) < 0.0


def test_absolute_zero_pressure_reference_uses_full_vibrational_pressure():
    eos = MieGruneisenDebye(
        BM3(1.0, 160.0, 4.0),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
        thermal_pressure_reference="absolute_zero",
    )
    volume = 0.9

    assert eos.thermal_pressure(volume, 300.0) == pytest.approx(
        eos.vibrational_pressure(volume, 300.0)
    )
    assert eos.thermal_pressure_increment(volume, 300.0) == pytest.approx(0.0)
    assert eos.thermal_pressure_increment(volume, 1800.0) == pytest.approx(
        eos.vibrational_pressure(volume, 1800.0)
        - eos.vibrational_pressure(volume, 300.0)
    )
    assert eos.configuration_values()["thermal_pressure_reference"] == "absolute_zero"


@pytest.mark.parametrize("value", ["zero", "cold_curve", None])
def test_invalid_thermal_pressure_reference_is_rejected(value):
    with pytest.raises(ValueError, match="thermal_pressure_reference"):
        MieGruneisenDebye(
            BM3(1.0, 160.0, 4.0),
            300.0,
            800.0,
            1.5,
            1.0,
            2.0,
            thermal_pressure_reference=value,
        )


def test_pressure_broadcasting(eos):
    volumes = eos.rt_eos.V0 * np.array([[0.8], [1.0]])
    temperatures = np.array([[500.0, 1000.0, 2000.0]])

    pressure = eos.thermal_pressure(volumes, temperatures)

    assert pressure.shape == (2, 3)
    assert np.all(np.isfinite(pressure))


def test_pressure_volume_round_trip(eos):
    expected_volume = 0.85 * eos.rt_eos.V0
    temperature = 1800.0
    pressure = eos.pressure(expected_volume, temperature)

    assert np.isclose(eos.volume(pressure, temperature), expected_volume, rtol=1.0e-10)


def test_characteristic_temperature_obeys_gruneisen_definition(eos):
    volume = 0.9 * eos.rt_eos.V0
    step = 1.0e-6 * volume
    numerical_gamma = -(
        np.log(eos.characteristic_temperature(volume + step))
        - np.log(eos.characteristic_temperature(volume - step))
    ) / (np.log(volume + step) - np.log(volume - step))

    assert np.isclose(numerical_gamma, eos.gruneisen_parameter(volume))


def test_zero_q_limit():
    eos = MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 0.0, 2.0)
    volume = 0.8 * eos.rt_eos.V0

    assert np.isclose(eos.gruneisen_parameter(volume), eos.gamma0)
    assert np.isclose(
        eos.characteristic_temperature(volume),
        eos.theta0 * (volume / eos.rt_eos.V0) ** -eos.gamma0,
    )


def test_debye_high_temperature_limit():
    eos = MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 1.0, 1.5, 1.0, 2.0)
    temperature = 10000.0

    assert np.isclose(
        eos.thermal_energy(eos.rt_eos.V0, temperature),
        3.0 * eos.n * R * temperature,
        rtol=5.0e-5,
    )


def test_einstein_energy_at_reference_volume():
    eos = MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0)
    temperature = 1200.0
    expected = 3.0 * eos.n * R * eos.theta0 / np.expm1(eos.theta0 / temperature)

    assert np.isclose(eos.thermal_energy(eos.rt_eos.V0, temperature), expected)


@pytest.mark.parametrize("argument", [1.0e-5, 1.0, 200.0])
def test_debye_function_is_finite_in_all_numerical_branches(argument):
    assert np.isfinite(_debye_function_3(argument))


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf])
def test_invalid_temperature_is_rejected(eos, temperature):
    with pytest.raises(ValueError, match="Temperature"):
        eos.thermal_pressure(eos.rt_eos.V0, temperature)


def test_incompatible_shapes_are_rejected(eos):
    with pytest.raises(ValueError, match="broadcast-compatible"):
        eos.thermal_pressure(np.ones(2), np.ones(3) * 300.0)


@pytest.mark.parametrize(
    "keyword,value",
    [
        ("Tr", 0.0),
        ("theta0", -1.0),
        ("gamma0", np.nan),
        ("q", np.inf),
        ("n", 0.0),
    ],
)
def test_invalid_parameters_are_rejected(keyword, value):
    parameters = {
        "Tr": 300.0,
        "theta0": 800.0,
        "gamma0": 1.5,
        "q": 1.0,
        "n": 2.0,
    }
    parameters[keyword] = value

    with pytest.raises(ValueError):
        MieGruneisenDebye(BM3(1.0, 160.0, 4.0), **parameters)


def test_non_eos_reference_is_rejected():
    with pytest.raises(TypeError, match="rt_eos"):
        MieGruneisenDebye(object(), 300.0, 800.0, 1.5, 1.0, 2.0)


def test_user_defined_reference_eos_retains_python_fallback():
    class CustomReference(EosBase):
        def __init__(self):
            self.V0 = 1.0

        def pressure(self, V):
            return 160.0 * (1.0 / np.asarray(V) - 1.0)

        def bulk_modulus(self, V):
            return 160.0 / np.asarray(V)

    eos = MieGruneisenDebye(CustomReference(), 300.0, 800.0, 1.5, 1.0, 2.0)
    expected_volume = 0.9
    temperature = 1200.0

    assert not hasattr(eos, "_native")
    pressure = eos.pressure(expected_volume, temperature)
    assert np.isclose(eos.volume(pressure, temperature), expected_volume)


def test_debye_subclass_retains_debye_python_behavior():
    class CustomDebye(MieGruneisenDebye):
        pass

    parameters = (BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0)
    expected = MieGruneisenDebye(*parameters)
    custom = CustomDebye(*parameters)

    assert not hasattr(custom, "_native")
    assert custom.thermal_energy(1.0, 300.0) == pytest.approx(
        expected.thermal_energy(1.0, 300.0)
    )
    pressure = custom.pressure(0.9, 1200.0)
    assert custom.volume(pressure, 1200.0) == pytest.approx(0.9)


def test_asymptotic_debye_excess_pressure_and_thermodynamics():
    eos = AsymptoticPowerLawMieGruneisenDebyeExcess(
        BM3(1.0, 160.0, 4.0),
        Tr=300.0,
        theta0=761.0,
        gamma0=1.52,
        a=1.0,
        b=1.43,
        n=2.0,
        beta0=-8.061e-4,
        m=4.8,
    )
    volumes = np.array([0.8, 1.0, 1.2])
    temperatures = np.array([1000.0, 1500.0, 2000.0])

    assert np.all(eos.thermal_pressure(volumes, 300.0) == 0.0)
    assert eos.excess_pressure(0.8, 2000.0, referenced=True) == pytest.approx(
        eos.excess_pressure(0.8, 2000.0) - eos.excess_pressure(0.8, 300.0)
    )
    assert np.allclose(
        eos.thermal_helmholtz_free_energy(volumes, temperatures),
        eos.thermal_energy(volumes, temperatures)
        - temperatures * eos.thermal_entropy(volumes, temperatures),
    )
    assert np.all(np.isfinite(eos.molar_heat_capacity_v(volumes, temperatures)))


def test_asymptotic_debye_excess_round_trip_and_parameters():
    eos = AsymptoticPowerLawMieGruneisenDebyeExcess(
        BM3(1.0, 160.0, 4.0), 300.0, 230.0, 2.75, 0.39, 5.1, 1.0, 0.002145, 0.65
    )
    volume = 0.8
    temperature = 2000.0
    pressure = eos.pressure(volume, temperature)

    assert eos.volume(pressure, temperature) == pytest.approx(volume)
    assert eos.parameter_values()["beta0"] == pytest.approx(0.002145)
    changed = eos.with_parameters(beta0=0.0)
    assert changed.beta0 == 0.0
    assert changed.rt_eos.parameter_values() == eos.rt_eos.parameter_values()
