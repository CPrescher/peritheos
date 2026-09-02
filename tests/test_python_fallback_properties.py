"""Property-oriented checks for native and custom-model compatibility paths."""

import numpy as np
import pytest

from peritheos import EosValidationError
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import (
    LinearThermalPressure,
    LogVolumeThermalPressure,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    MultiOscillatorGruneisenThermalEOS,
    Tange2009Debye,
    ThermalReferenceStateEOS,
)


class PythonDebye(MieGruneisenDebye):
    """Force the inherited Python compatibility implementation."""


class PythonEinstein(MieGruneisenEinstein):
    """Force the inherited Python compatibility implementation."""


class PythonLinear(LinearThermalPressure):
    """Force the inherited Python compatibility implementation."""


class PythonLogVolume(LogVolumeThermalPressure):
    """Force the inherited Python compatibility implementation."""


class PythonReferenceState(ThermalReferenceStateEOS):
    """Force the inherited Python compatibility implementation."""


class PythonMultiOscillator(MultiOscillatorGruneisenThermalEOS):
    """Force the inherited Python compatibility implementation."""


class PythonTange(Tange2009Debye):
    """Force the inherited Python compatibility implementation."""


def _model_pair(model_class, python_class):
    reference = BM3(V0=1.0, K0=160.0, K0_prime=4.1)
    common = {"rt_eos": reference, "Tr": 300.0}
    if model_class is LinearThermalPressure:
        parameters = {"alpha_KT": 0.005}
    elif model_class is LogVolumeThermalPressure:
        parameters = {"alpha_KT_ref": 0.005, "dK_dT_V": -0.002}
    elif model_class is ThermalReferenceStateEOS:
        parameters = {"alpha0": 3.0e-5, "dK_dT": -0.02}
    elif model_class is MultiOscillatorGruneisenThermalEOS:
        parameters = {
            "QE1o": 1500.0,
            "mE1": 2.4,
            "QE2o": 700.0,
            "mE2": 0.6,
            "delta": -0.5,
            "t": 1.1,
            "a_0": 0.0,
            "m": 0.0,
            "g": 0.0,
            "e_0": 0.0,
            "n": 1.0,
        }
    elif model_class is Tange2009Debye:
        parameters = {
            "theta0": 760.0,
            "gamma0": 1.45,
            "a": 0.14,
            "b": 5.4,
            "n": 2.0,
        }
    else:
        parameters = {
            "theta0": 800.0,
            "gamma0": 1.5,
            "q": 1.0,
            "n": 2.0,
        }
    return model_class(**common, **parameters), python_class(**common, **parameters)


@pytest.mark.parametrize(
    ("model_class", "python_class"),
    [
        (LinearThermalPressure, PythonLinear),
        (LogVolumeThermalPressure, PythonLogVolume),
        (MieGruneisenDebye, PythonDebye),
        (MieGruneisenEinstein, PythonEinstein),
        (MultiOscillatorGruneisenThermalEOS, PythonMultiOscillator),
        (Tange2009Debye, PythonTange),
        (ThermalReferenceStateEOS, PythonReferenceState),
    ],
)
def test_python_compatibility_path_matches_native_on_deterministic_state_grid(
    model_class, python_class
):
    native, python = _model_pair(model_class, python_class)
    rng = np.random.default_rng(20260902)
    volumes = rng.uniform(0.78, 0.98, 8)
    temperatures = rng.uniform(350.0, 2500.0, 8)

    native_pressure = native.pressure(volumes, temperatures)
    python_pressure = python.pressure(volumes, temperatures)

    assert np.allclose(python_pressure, native_pressure, rtol=2.0e-10, atol=1.0e-11)
    assert np.allclose(
        python.volume(native_pressure, temperatures),
        volumes,
        rtol=2.0e-9,
        atol=1.0e-11,
    )
    assert np.allclose(
        python.temperature(native_pressure, volumes),
        temperatures,
        rtol=2.0e-9,
        atol=1.0e-7,
    )
    cold_pressures = np.linspace(10.0, 60.0, volumes.size)
    assert np.allclose(
        python.volume_with_dac_confinement(cold_pressures, temperatures, f_dac=0.25),
        native.volume_with_dac_confinement(cold_pressures, temperatures, f_dac=0.25),
        rtol=2.0e-9,
        atol=1.0e-11,
    )
    assert np.allclose(
        python.thermal_pressure_increment(volumes, temperatures),
        native.thermal_pressure_increment(volumes, temperatures),
        rtol=2.0e-10,
        atol=1.0e-11,
    )
    assert np.allclose(
        python.bulk_modulus(volumes, temperatures),
        native.bulk_modulus(volumes, temperatures),
        rtol=2.0e-5,
        atol=1.0e-7,
    )
    assert np.allclose(
        python.thermal_expansivity(volumes, temperatures),
        native.thermal_expansivity(volumes, temperatures),
        rtol=2.0e-5,
        atol=1.0e-10,
    )


@pytest.mark.parametrize("python_class", [PythonDebye, PythonEinstein])
def test_python_caloric_path_preserves_thermodynamic_ordering(python_class):
    model_class = (
        MieGruneisenDebye if python_class is PythonDebye else MieGruneisenEinstein
    )
    _, eos = _model_pair(model_class, python_class)
    volumes = np.array([0.82, 0.9, 0.98])
    temperatures = np.array([500.0, 1200.0, 2200.0])

    cv = eos.molar_heat_capacity_v(volumes, temperatures)
    cp = eos.molar_heat_capacity_p(volumes, temperatures)
    kt = eos.bulk_modulus(volumes, temperatures)
    ks = eos.adiabatic_bulk_modulus(volumes, temperatures)

    assert np.all(np.asarray(cv) > 0.0)
    assert np.all(np.asarray(cp) >= np.asarray(cv))
    assert np.all(np.asarray(ks) >= np.asarray(kt))


@pytest.mark.parametrize(
    ("model_class", "python_class"),
    [
        (MieGruneisenDebye, PythonDebye),
        (MieGruneisenEinstein, PythonEinstein),
    ],
)
def test_python_caloric_helpers_match_native_for_scalars_and_arrays(
    model_class, python_class
):
    native, python = _model_pair(model_class, python_class)
    states = [
        (0.87, 1400.0),
        (np.array([0.82, 0.9, 0.98]), np.array([500.0, 1200.0, 2200.0])),
    ]
    state_methods = (
        "thermal_energy",
        "thermal_entropy",
        "thermal_internal_energy",
        "vibrational_pressure",
        "thermal_helmholtz_free_energy",
        "thermal_enthalpy",
        "thermal_gibbs_free_energy",
        "molar_heat_capacity_v",
    )

    for volume, temperature in states:
        for method_name in state_methods:
            expected = getattr(native, method_name)(volume, temperature)
            actual = getattr(python, method_name)(volume, temperature)
            assert np.allclose(actual, expected, rtol=2.0e-8, atol=1.0e-8), method_name

        assert np.allclose(
            python.characteristic_temperature(volume),
            native.characteristic_temperature(volume),
            rtol=2.0e-12,
        )
        assert np.allclose(
            python.gruneisen_parameter(volume),
            native.gruneisen_parameter(volume),
            rtol=2.0e-12,
        )


@pytest.mark.parametrize("python_class", [PythonDebye, PythonEinstein])
def test_python_caloric_helpers_reject_invalid_state_shapes(python_class):
    model_class = (
        MieGruneisenDebye if python_class is PythonDebye else MieGruneisenEinstein
    )
    _, eos = _model_pair(model_class, python_class)

    with pytest.raises(EosValidationError, match="Temperature"):
        eos.thermal_energy(0.9, 0.0)
    with pytest.raises(EosValidationError, match="broadcast-compatible"):
        eos.thermal_energy(np.ones(2), np.ones(3) * 1000.0)


def test_python_debye_zero_q_and_tange_zero_b_limits_are_finite():
    reference = BM3(V0=1.0, K0=160.0, K0_prime=4.1)
    debye = PythonDebye(reference, 300.0, 800.0, 1.5, 0.0, 2.0)
    tange = PythonTange(reference, 300.0, 760.0, 1.45, 0.14, 0.0, 2.0)

    expected = 800.0 * (0.85 ** (-1.5))
    assert debye.characteristic_temperature(0.85) == pytest.approx(expected)
    assert np.isfinite(tange.characteristic_temperature(0.85))
