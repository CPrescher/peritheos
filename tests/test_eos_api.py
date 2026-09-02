import numpy as np
import pytest

from peritheos.constants import R
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
from peritheos.units import convert_pressure, convert_temperature
from peritheos.utils import (
    compressibility_factor,
    derivative,
)


@pytest.mark.parametrize(
    "eos",
    [
        BM2(10.0, 100.0),
        BM3(10.0, 100.0, 4.0),
        BM4(10.0, 100.0, 4.0, -0.01),
        Murnaghan(10.0, 100.0, 4.0),
        ModifiedTait(10.0, 100.0, 4.0, -0.01),
        NaturalStrain2(10.0, 100.0),
        NaturalStrain3(10.0, 100.0, 4.0),
        NaturalStrain4(10.0, 100.0, 4.0, -0.01),
        Vinet(10.0, 100.0, 4.0),
        Holzapfel(0.3414, 441.5, 3.9, 1, 6),
    ],
)
def test_room_temperature_pressure_volume_round_trip(eos):
    expected_volume = eos.V0 * 0.8
    pressure = eos.pressure(expected_volume)

    assert np.isclose(eos.volume(pressure), expected_volume, rtol=1e-10)


def test_subclass_volume_inversion_honors_overridden_pressure():
    class ShiftedBM3(BM3):
        def pressure(self, V):
            return np.asarray(super().pressure(V)) + 2.0

    eos = ShiftedBM3(10.0, 120.0, 4.0)
    expected_volume = eos.V0
    pressure = eos.pressure(expected_volume)

    recovered = eos.volume(pressure)

    assert recovered == pytest.approx(expected_volume)
    assert eos.pressure(recovered) == pytest.approx(pressure)


def test_volume_solver_rejects_pressure_outside_expansion_branch():
    eos = BM2(10.0, 100.0)

    with pytest.raises(ValueError, match="outside the invertible expansion range"):
        eos.calculate_volume(-100.0)


@pytest.mark.parametrize("volume", [0, -1, np.nan, np.inf])
def test_eos_rejects_invalid_volume(volume):
    with pytest.raises(ValueError):
        BM3(10.0, 100.0, 4.0).pressure(volume)


@pytest.mark.parametrize("V0,K0", [(0, 100), (10, 0), (-1, 100), (10, np.nan)])
def test_eos_rejects_invalid_reference_parameters(V0, K0):
    with pytest.raises(ValueError):
        BM2(V0, K0)


def test_pressure_conversion_supports_eos_units():
    assert convert_pressure(1, "GPa", "kbar") == 10
    assert convert_pressure(10, "kbar", "GPa") == 1


def test_pressure_conversion_reports_unsupported_units():
    with pytest.raises(ValueError, match="Unsupported pressure unit"):
        convert_pressure(1, "bogus", "GPa")
    with pytest.raises(ValueError, match="Unsupported pressure unit"):
        convert_pressure(1, "GPa", "bogus")


@pytest.mark.parametrize(
    "value,source,target,expected",
    [
        (273.15, "K", "C", 0.0),
        (0.0, "C", "K", 273.15),
        (32.0, "F", "C", 0.0),
        (100.0, "C", "F", 212.0),
        (32.0, "F", "K", 273.15),
        (273.15, "K", "F", 32.0),
    ],
)
def test_temperature_conversion(value, source, target, expected):
    assert np.isclose(convert_temperature(value, source, target), expected)


def test_temperature_conversion_rejects_unknown_units():
    with pytest.raises(ValueError, match="Unsupported temperature unit"):
        convert_temperature(1.0, "rankine", "K")
    with pytest.raises(ValueError, match="Unsupported temperature unit"):
        convert_temperature(1.0, "K", "rankine")


def test_general_thermodynamic_utilities():
    assert np.isclose(compressibility_factor(2.0 * R * 300.0, 1.0, 300.0, 2.0), 1.0)
    assert np.isclose(derivative(lambda value: value**2, 3.0), 6.0)


def test_legacy_unit_imports_warn_and_delegate():
    from peritheos import utils

    with pytest.deprecated_call(match="peritheos.units"):
        assert utils.convert_pressure(1.0, "GPa", "kbar") == 10.0
    with pytest.deprecated_call(match="peritheos.units"):
        assert utils.convert_temperature(273.15, "K", "C") == 0.0
