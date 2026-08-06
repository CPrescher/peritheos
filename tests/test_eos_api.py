import numpy as np
import pytest

from peritheos.eos.rt import BM2, BM3, BM4, Holzapfel, Vinet
from peritheos.utils import convert_pressure


@pytest.mark.parametrize(
    "eos",
    [
        BM2(10.0, 100.0),
        BM3(10.0, 100.0, 4.0),
        BM4(10.0, 100.0, 4.0, -0.01),
        Vinet(10.0, 100.0, 4.0),
        Holzapfel(0.3414, 441.5, 3.9, 1, 6),
    ],
)
def test_room_temperature_pressure_volume_round_trip(eos):
    expected_volume = eos.V0 * 0.8
    pressure = eos.pressure(expected_volume)

    assert np.isclose(eos.volume(pressure), expected_volume, rtol=1e-10)


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
