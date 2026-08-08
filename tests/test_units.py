"""Tests for explicit molar-volume and density conversions."""

import numpy as np
import pytest

from peritheos.units import (
    convert_density,
    convert_molar_volume,
    density_from_molar_volume,
    molar_volume_from_density,
)


def test_thermal_volume_convention_conversion():
    assert np.isclose(convert_molar_volume(10.0, "cm^3/mol", "J/bar/mol"), 1.0)
    assert np.isclose(convert_molar_volume(1.0, "J bar^-1 mol^-1", "cm3/mol"), 10.0)


def test_array_conversion():
    values = np.array([1.0, 2.0])
    assert np.allclose(convert_density(values, "g/cm3", "kg/m3"), values * 1000.0)


def test_density_volume_round_trip():
    density = density_from_molar_volume(40.304, 11.25)
    volume = molar_volume_from_density(40.304, density)

    assert np.isclose(volume, 11.25)


@pytest.mark.parametrize("unit", ["cubits/mol", "", "cm/mol"])
def test_unsupported_volume_unit(unit):
    with pytest.raises(ValueError, match="Unsupported unit"):
        convert_molar_volume(1.0, unit, "cm3/mol")


def test_nonpositive_density_and_volume_are_rejected():
    with pytest.raises(ValueError, match="Molar volume"):
        density_from_molar_volume(40.0, 0.0)
    with pytest.raises(ValueError, match="Density"):
        molar_volume_from_density(40.0, 0.0)
    with pytest.raises(ValueError, match="Molar mass"):
        density_from_molar_volume(0.0, 10.0)
