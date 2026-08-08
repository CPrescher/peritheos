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
