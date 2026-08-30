"""Deterministic stress grids across the supported EOS families."""

import numpy as np
import pytest

from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    Vinet,
)
from peritheos.eos.thermal import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Sokolova2016,
    ThermalModifiedTait,
)


@pytest.mark.parametrize(
    "eos",
    [
        BM2(10.0, 120.0),
        BM3(10.0, 120.0, 4.3),
        BM4(10.0, 120.0, 4.3, -0.02),
        Murnaghan(10.0, 120.0, 4.3),
        ModifiedTait(10.0, 120.0, 4.3, -0.02),
        NaturalStrain2(10.0, 120.0),
        NaturalStrain3(10.0, 120.0, 4.3),
        NaturalStrain4(10.0, 120.0, 4.3, -0.02),
        Vinet(10.0, 120.0, 4.3),
        Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
    ],
)
def test_room_temperature_compression_grid_is_finite_and_invertible(eos):
    volumes = eos.V0 * np.array([0.55, 0.7, 0.85, 1.0])
    pressures = np.asarray(eos.pressure(volumes))
    moduli = np.asarray(eos.bulk_modulus(volumes))

    assert np.all(np.isfinite(pressures))
    assert np.all(np.isfinite(moduli))
    assert np.allclose(eos.volume(pressures), volumes, rtol=1.0e-9)


@pytest.fixture(
    params=[
        MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
        MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
        ThermalModifiedTait(
            ModifiedTait(1.0, 160.0, 4.0, -0.01),
            298.15,
            700.0,
            2.5e-5,
            2.0,
        ),
        Sokolova2016(
            Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
            298.15,
            684.0,
            0.564,
            1561.0,
            2.436,
            -0.506,
            1.085,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
        Sokolova2016(
            Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
            298.15,
            684.0,
            0.564,
            1561.0,
            2.436,
            -0.506,
            1.085,
            5.2,
            1.3,
            0.8,
            2.7,
            beta=0.35,
            QBo=480.0,
            d=2.4,
            mb=0.75,
            QB1o=1120.0,
            d1=1.6,
            mb1=0.4,
        ),
    ]
)
def thermal_eos(request):
    return request.param


def test_thermal_state_grid_is_finite_and_representative_state_inverts(thermal_eos):
    volumes = thermal_eos.rt_eos.V0 * np.array([0.7, 0.85, 1.0])[:, None]
    temperatures = np.array([300.0, 1500.0, 4000.0])[None, :]

    pressures = np.asarray(thermal_eos.pressure(volumes, temperatures))
    moduli = np.asarray(thermal_eos.bulk_modulus(volumes, temperatures))

    assert pressures.shape == (3, 3)
    assert moduli.shape == (3, 3)
    assert np.all(np.isfinite(pressures))
    assert np.all(np.isfinite(moduli))

    expected_volume = 0.8 * thermal_eos.rt_eos.V0
    pressure = thermal_eos.pressure(expected_volume, 2500.0)
    assert np.isclose(
        thermal_eos.volume(pressure, 2500.0), expected_volume, rtol=1.0e-9
    )
    assert np.isclose(
        thermal_eos.temperature(pressure, expected_volume), 2500.0, rtol=1.0e-9
    )


def test_thermal_temperature_inversion_broadcasts(thermal_eos):
    volumes = thermal_eos.rt_eos.V0 * np.array([[0.75], [0.9]])
    expected_temperatures = np.array([[200.0, 1500.0, 3000.0]])
    pressures = thermal_eos.pressure(volumes, expected_temperatures)

    temperatures = thermal_eos.calculate_temperature(pressures, volumes)

    assert temperatures.shape == (2, 3)
    assert np.allclose(temperatures, expected_temperatures, rtol=1.0e-9)


def test_thermal_temperature_inversion_recovers_reference_temperature(thermal_eos):
    volume = 0.85 * thermal_eos.rt_eos.V0
    pressure = thermal_eos.pressure(volume, thermal_eos.Tr)

    assert thermal_eos.temperature(pressure, volume) == thermal_eos.Tr


def test_thermal_temperature_inversion_with_dac_pressure(thermal_eos):
    ambient_volume = 0.8 * thermal_eos.rt_eos.V0
    heated_volume = 0.80001 * thermal_eos.rt_eos.V0
    ambient_pressure = thermal_eos.rt_eos.pressure(ambient_volume)
    uncorrected_temperature = thermal_eos.temperature(ambient_pressure, heated_volume)
    corrected_temperature = thermal_eos.temperature_from_volumes(
        ambient_volume,
        heated_volume,
        f_dac=0.2,
    )
    thermal_pressure = thermal_eos.thermal_pressure(
        heated_volume, corrected_temperature
    )
    dac_pressure = thermal_eos.dac_thermal_pressure(
        heated_volume, corrected_temperature, f_dac=0.2
    )

    assert corrected_temperature > uncorrected_temperature
    assert np.isclose(
        thermal_pressure,
        (
            thermal_eos.rt_eos.pressure(ambient_volume)
            - thermal_eos.rt_eos.pressure(heated_volume)
        )
        / 0.8,
    )
    assert np.isclose(
        thermal_eos.pressure(heated_volume, corrected_temperature),
        ambient_pressure + dac_pressure,
    )

    temperature_at_larger_fraction = thermal_eos.temperature_from_volumes(
        ambient_volume,
        heated_volume,
        f_dac=0.3,
    )
    assert temperature_at_larger_fraction > corrected_temperature


@pytest.mark.parametrize("f_dac", [-0.1, 1.0, 1.1, np.nan, np.inf])
def test_thermal_temperature_inversion_rejects_invalid_dac_fraction(thermal_eos, f_dac):
    with pytest.raises(ValueError, match="f_dac"):
        thermal_eos.temperature_from_volumes(
            0.8 * thermal_eos.rt_eos.V0,
            thermal_eos.rt_eos.V0,
            f_dac=f_dac,
        )


def test_temperature_from_volumes_zero_dac_matches_isobaric_inversion(thermal_eos):
    ambient_volume = 0.8 * thermal_eos.rt_eos.V0
    heated_volume = 0.80001 * thermal_eos.rt_eos.V0
    ambient_pressure = thermal_eos.rt_eos.pressure(ambient_volume)

    assert np.isclose(
        thermal_eos.temperature_from_volumes(
            ambient_volume,
            heated_volume,
            f_dac=0.0,
        ),
        thermal_eos.temperature(ambient_pressure, heated_volume),
    )


def test_temperature_from_volumes_broadcasts(thermal_eos):
    ambient_volumes = thermal_eos.rt_eos.V0 * np.array([[0.80], [0.82]])
    heated_volumes = ambient_volumes + thermal_eos.rt_eos.V0 * np.array(
        [[1.0e-5, 2.0e-5, 3.0e-5]]
    )

    temperatures = thermal_eos.temperature_from_volumes(
        ambient_volumes, heated_volumes, f_dac=0.25
    )

    assert temperatures.shape == (2, 3)
    assert np.all(np.isfinite(temperatures))
    expected_thermal_pressures = (
        thermal_eos.rt_eos.pressure(ambient_volumes)
        - thermal_eos.rt_eos.pressure(heated_volumes)
    ) / 0.75
    assert np.allclose(
        thermal_eos.thermal_pressure(heated_volumes, temperatures),
        expected_thermal_pressures,
    )


def test_temperature_from_volumes_rejects_nonheated_state(thermal_eos):
    with pytest.raises(ValueError, match="below the reference temperature"):
        thermal_eos.temperature_from_volumes(
            0.8 * thermal_eos.rt_eos.V0,
            0.79 * thermal_eos.rt_eos.V0,
            f_dac=0.25,
        )


def test_temperature_from_volumes_rejects_noninvertible_thermal_pressure():
    eos = MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 0.0, 1.0, 2.0)

    with pytest.raises(ValueError, match="invertible temperature range"):
        eos.temperature_from_volumes(0.8, 0.80001, f_dac=0.25)


def test_temperature_from_volumes_rejects_root_below_reference_temperature():
    eos = MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, -1.5, 1.0, 2.0)

    with pytest.raises(ValueError, match="below the reference temperature"):
        eos.temperature_from_volumes(0.8, 0.80001, f_dac=0.25)


def test_temperature_from_volumes_rejects_incompatible_shapes(thermal_eos):
    with pytest.raises(ValueError, match="broadcast-compatible"):
        thermal_eos.temperature_from_volumes(np.ones(2), np.ones(3), f_dac=0.25)


@pytest.mark.parametrize("volume", [0.0, -1.0, np.nan, np.inf])
def test_thermal_temperature_inversion_rejects_invalid_volume(thermal_eos, volume):
    with pytest.raises(ValueError, match="Volume"):
        thermal_eos.temperature(1.0, volume)


@pytest.mark.parametrize("pressure", [np.nan, np.inf, -np.inf])
def test_thermal_temperature_inversion_rejects_invalid_pressure(thermal_eos, pressure):
    with pytest.raises(ValueError, match="Pressure"):
        thermal_eos.temperature(pressure, thermal_eos.rt_eos.V0)


def test_thermal_temperature_inversion_rejects_incompatible_shapes(thermal_eos):
    with pytest.raises(ValueError, match="broadcast-compatible"):
        thermal_eos.temperature(np.ones(2), np.ones(3))
