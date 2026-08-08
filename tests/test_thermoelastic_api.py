"""Tests for common temperature-aware thermoelastic properties."""

import numpy as np
import pytest

from peritheos.eos.rt import BM3, ModifiedTait
from peritheos.eos.thermal import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    ThermalModifiedTait,
)


@pytest.fixture(params=[MieGruneisenDebye, MieGruneisenEinstein])
def eos(request):
    return request.param(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0)


def test_common_thermoelastic_properties(eos):
    volume = 0.9
    temperature = 1500.0
    kt = eos.bulk_modulus(volume, temperature)
    alpha = eos.thermal_expansivity(volume, temperature)
    cv = eos.molar_heat_capacity_v(volume, temperature)
    cp = eos.molar_heat_capacity_p(volume, temperature)
    ks = eos.adiabatic_bulk_modulus(volume, temperature)

    assert kt > 0.0
    assert alpha > 0.0
    assert cv > 0.0
    assert cp > cv
    assert ks > kt
    assert np.isclose(eos.isothermal_compressibility(volume, temperature), 1.0 / kt)


def test_thermodynamic_potential_identities(eos):
    volume = 0.9
    temperature = 1200.0
    internal_energy = eos.thermal_internal_energy(volume, temperature)
    entropy = eos.thermal_entropy(volume, temperature)

    assert np.isclose(
        eos.thermal_helmholtz_free_energy(volume, temperature),
        internal_energy - temperature * entropy,
    )
    assert np.isclose(
        eos.thermal_enthalpy(volume, temperature),
        internal_energy
        + eos.vibrational_pressure(volume, temperature) * volume * 1.0e4,
    )


def test_gruneisen_parameter_matches_thermodynamic_definition(eos):
    volume = 0.9
    temperature = 1200.0
    expected = (
        eos.thermal_expansivity(volume, temperature)
        * eos.bulk_modulus(volume, temperature)
        * volume
        * 1.0e4
        / eos.molar_heat_capacity_v(volume, temperature)
    )

    assert np.isclose(
        eos.gruneisen_parameter(volume, temperature), expected, rtol=1.0e-5
    )


def test_heat_capacity_is_internal_energy_temperature_derivative(eos):
    volume = 0.9
    temperature = 1200.0
    step = 0.01
    numerical = (
        eos.thermal_internal_energy(volume, temperature + step)
        - eos.thermal_internal_energy(volume, temperature - step)
    ) / (2.0 * step)

    assert np.isclose(
        eos.molar_heat_capacity_v(volume, temperature), numerical, rtol=1.0e-7
    )


def test_free_energy_derivatives_recover_entropy_and_vibrational_pressure(eos):
    volume = 0.9
    temperature = 1200.0
    temperature_step = 0.01
    volume_step = 1.0e-6
    entropy_from_free_energy = -(
        eos.thermal_helmholtz_free_energy(volume, temperature + temperature_step)
        - eos.thermal_helmholtz_free_energy(volume, temperature - temperature_step)
    ) / (2.0 * temperature_step)
    pressure_from_free_energy = -(
        eos.thermal_helmholtz_free_energy(volume + volume_step, temperature)
        - eos.thermal_helmholtz_free_energy(volume - volume_step, temperature)
    ) / (2.0 * volume_step) / 1.0e4

    assert np.isclose(
        entropy_from_free_energy,
        eos.thermal_entropy(volume, temperature),
        rtol=1.0e-7,
    )
    assert np.isclose(
        pressure_from_free_energy,
        eos.vibrational_pressure(volume, temperature),
        rtol=1.0e-7,
    )


def test_thermal_modified_tait_reference_properties():
    eos = ThermalModifiedTait(
        ModifiedTait(1.0, 160.0, 4.0, -0.01),
        Tr=298.15,
        theta=700.0,
        alpha0=2.5e-5,
        n=2.0,
    )

    assert eos.thermal_pressure(eos.rt_eos.V0, eos.Tr) == 0.0
    assert np.isclose(eos.bulk_modulus(eos.rt_eos.V0, eos.Tr), eos.rt_eos.K0)
    assert np.isclose(
        eos.thermal_expansivity(eos.rt_eos.V0, eos.Tr), eos.alpha0, rtol=1.0e-5
    )
    assert np.isfinite(eos.gruneisen_parameter(eos.rt_eos.V0, eos.Tr))
    pressure = eos.pressure(0.9, 1200.0)
    assert np.isclose(eos.volume(pressure, 1200.0), 0.9)


def test_thermal_modified_tait_requires_modified_tait_reference():
    with pytest.raises(TypeError, match="ModifiedTait"):
        ThermalModifiedTait(BM3(1.0, 160.0, 4.0), 300.0, 700.0, 2.5e-5, 2.0)
