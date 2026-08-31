import numpy as np
import pytest

from peritheos.eos import EosBase
from peritheos.eos.rt import BM2, Vinet
from peritheos.eos.thermal import LinearThermalPressure, ThermalReferenceStateEOS


def test_linear_thermal_pressure_reference_identity_and_arrays():
    eos = LinearThermalPressure(Vinet(50.0, 20.0, 5.0), 300.0, 0.002)
    volumes = np.array([50.0, 45.0])

    assert eos.thermal_pressure(50.0, 300.0) == pytest.approx(0.0)
    assert np.allclose(eos.thermal_pressure(volumes, 800.0), [1.0, 1.0])
    assert np.allclose(eos.pressure(volumes, 800.0), eos.rt_eos.pressure(volumes) + 1.0)


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan])
def test_linear_thermal_pressure_rejects_invalid_temperature(temperature):
    eos = LinearThermalPressure(Vinet(50.0, 20.0, 5.0), 300.0, 0.002)
    with pytest.raises(ValueError):
        eos.thermal_pressure(40.0, temperature)


def test_linear_thermal_pressure_requires_broadcast_compatible_states():
    eos = LinearThermalPressure(Vinet(50.0, 20.0, 5.0), 300.0, 0.002)
    with pytest.raises(ValueError, match="broadcast-compatible"):
        eos.thermal_pressure(np.ones(2), np.ones(3) * 300.0)


def _bm2_pressure(volume, reference_volume, bulk_modulus):
    compression = reference_volume / np.asarray(volume, dtype=float)
    return (
        1.5 * bulk_modulus * (compression ** (7.0 / 3.0) - compression ** (5.0 / 3.0))
    )


def test_thermal_reference_state_matches_bezacier_equations_and_arrays():
    eos = ThermalReferenceStateEOS(
        BM2(V0=41.480265898, K0=20.15),
        Tr=300.0,
        alpha0=11.58e-5,
        dK_dT=0.0,
    )
    volumes = np.array([41.480265898, 32.406])
    temperatures = np.array([300.0, 300.6])
    shifted_volumes = eos.rt_eos.V0 * np.exp(eos.alpha0 * (temperatures - eos.Tr))
    shifted_moduli = eos.rt_eos.K0 + eos.dK_dT * (temperatures - eos.Tr)

    expected = _bm2_pressure(volumes, shifted_volumes, shifted_moduli)

    assert eos.pressure(eos.rt_eos.V0, eos.Tr) == pytest.approx(0.0, abs=1e-14)
    assert np.allclose(eos.pressure(volumes, temperatures), expected)
    assert np.allclose(
        eos.thermal_pressure(volumes, temperatures),
        expected - eos.rt_eos.pressure(volumes),
    )


def test_thermal_reference_state_volume_and_temperature_round_trips():
    eos = ThermalReferenceStateEOS(
        BM2(V0=235.2983858185, K0=14.05),
        Tr=300.0,
        alpha0=14.6e-5,
        dK_dT=0.0,
    )
    volumes = np.array([206.233, 214.814])
    temperatures = np.array([340.7, 298.7])
    pressures = eos.pressure(volumes, temperatures)

    assert np.allclose(eos.volume(pressures, temperatures), volumes)
    assert np.allclose(eos.temperature(pressures, volumes), temperatures)

    shifted_reference_volume = eos.rt_eos.V0 * np.exp(
        eos.alpha0 * (temperatures[0] - eos.Tr)
    )
    assert eos.bulk_modulus(shifted_reference_volume, temperatures[0]) == pytest.approx(
        eos.rt_eos.K0
    )


@pytest.mark.parametrize("relative_step", [0.0, -1.0, np.nan])
def test_thermal_reference_state_rejects_invalid_bulk_modulus_step(relative_step):
    eos = ThermalReferenceStateEOS(BM2(10.0, 20.0), 300.0, 1e-5, 0.0)
    with pytest.raises(ValueError):
        eos.bulk_modulus(9.0, 400.0, relative_step=relative_step)


def test_thermal_reference_state_requires_reconstructable_reference_parameters():
    class MissingBulkModulus(EosBase):
        def __init__(self, V0):
            self.V0 = V0

    with pytest.raises(ValueError, match="V0 and K0"):
        ThermalReferenceStateEOS(
            MissingBulkModulus(V0=10.0),
            Tr=300.0,
            alpha0=1e-5,
            dK_dT=0.0,
        )


def test_thermal_reference_state_rejects_nonpositive_shifted_modulus():
    eos = ThermalReferenceStateEOS(BM2(10.0, 20.0), 300.0, 1e-5, -0.1)
    with pytest.raises(ValueError, match="non-positive bulk modulus"):
        eos.pressure(9.0, 501.0)


def test_thermal_reference_state_rejects_nonfinite_shifted_volume():
    eos = ThermalReferenceStateEOS(BM2(10.0, 20.0), 300.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="non-positive reference volume"):
        eos.pressure(9.0, 1.0e4)
