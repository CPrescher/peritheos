import numpy as np
import pytest

from peritheos.eos import EosBase
from peritheos.eos.rt import BM2, BM3, Vinet
from peritheos.eos.thermal import (
    LinearThermalPressure,
    LogVolumeThermalPressure,
    SecondOrderTaylorThermalPressure,
    ThermalReferenceStateEOS,
)
from peritheos.fitting import fit_thermal_eos


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


def test_log_volume_thermal_pressure_matches_anderson_equation_and_arrays():
    reference = BM2(V0=67.79, K0=166.65)
    eos = LogVolumeThermalPressure(
        reference,
        Tr=300.0,
        alpha_KT_ref=0.00714,
        dK_dT_V=-0.0115,
    )
    volumes = reference.V0 * np.array([1.0, 0.8, 0.66])
    temperatures = np.array([3000.0, 1000.0, 300.0])
    expected = (0.00714 - 0.0115 * np.log(reference.V0 / volumes)) * (
        temperatures - 300.0
    )

    assert eos.thermal_pressure(reference.V0, 300.0) == pytest.approx(0.0)
    assert np.allclose(eos.thermal_pressure(volumes, temperatures), expected)
    assert np.allclose(
        eos.pressure(volumes, temperatures), reference.pressure(volumes) + expected
    )


def test_log_volume_thermal_pressure_round_trips_and_requires_reference_volume():
    eos = LogVolumeThermalPressure(
        BM2(V0=67.79, K0=166.65),
        Tr=300.0,
        alpha_KT_ref=0.00714,
        dK_dT_V=-0.0115,
    )
    volumes = np.array([60.0, 50.0])
    temperatures = np.array([1000.0, 2000.0])
    pressures = eos.pressure(volumes, temperatures)

    assert np.allclose(eos.volume(pressures, temperatures), volumes)
    assert np.allclose(eos.temperature(pressures, volumes), temperatures)

    class MissingReferenceVolume(EosBase):
        pass

    with pytest.raises(ValueError, match="expose V0"):
        LogVolumeThermalPressure(
            MissingReferenceVolume(),
            Tr=300.0,
            alpha_KT_ref=0.00714,
            dK_dT_V=-0.0115,
        )


def test_second_order_taylor_thermal_pressure_is_absolute_and_round_trips():
    reference = Vinet(V0=74.0741025123, K0=169.8, K0_prime=4.501)
    eos = SecondOrderTaylorThermalPressure(
        reference,
        Tr=300.0,
        eta0=0.02,
        c0=0.5096,
        c1=-13.4246,
        c2=6.3295e-3,
        c3=36.2194,
        c4=5.4705e-8,
        c5=3.2238e-3,
    )
    volume = reference.V0 * 0.8
    temperature = 3000.0
    delta_eta = 0.2 - 0.02
    delta_temperature = temperature - 300.0
    expected_thermal = (
        0.5096
        - 13.4246 * delta_eta
        + 6.3295e-3 * delta_temperature
        + 0.5 * 36.2194 * delta_eta**2
        + 0.5 * 5.4705e-8 * delta_temperature**2
        + 0.5 * 3.2238e-3 * delta_eta * delta_temperature
    )

    assert eos.thermal_pressure(volume, temperature) == pytest.approx(expected_thermal)
    assert eos.thermal_pressure(volume, 300.0) != pytest.approx(0.0)
    assert eos.thermal_pressure_increment(volume, 300.0) == pytest.approx(0.0)
    pressure = eos.pressure(volume, temperature)
    assert eos.volume(pressure, temperature) == pytest.approx(volume)
    assert eos.temperature(pressure, volume) == pytest.approx(temperature)


def test_second_order_taylor_thermal_pressure_broadcasts_and_validates():
    eos = SecondOrderTaylorThermalPressure(
        Vinet(50.0, 20.0, 5.0),
        300.0,
        0.02,
        0.5,
        -13.0,
        0.006,
        36.0,
        5.0e-8,
        0.003,
    )
    result = eos.thermal_pressure(np.array([50.0, 45.0]), 1000.0)
    assert result.shape == (2,)
    with pytest.raises(ValueError):
        eos.thermal_pressure(45.0, 0.0)


def test_second_order_taylor_thermal_pressure_can_be_refitted():
    reference = Vinet(1.0, 169.8, 4.501)
    expected = SecondOrderTaylorThermalPressure(
        reference,
        300.0,
        0.02,
        0.5096,
        -13.4246,
        0.0063295,
        36.2194,
        5.4705e-8,
        0.0032238,
    )
    volumes = np.tile(np.array([0.65, 0.75, 0.85, 0.95]), 3)
    temperatures = np.repeat(np.array([500.0, 3000.0, 7000.0]), 4)
    pressures = expected.pressure(volumes, temperatures)
    result = fit_thermal_eos(
        SecondOrderTaylorThermalPressure,
        reference,
        volumes,
        temperatures,
        pressures,
        initial={"c2": 0.005},
        fixed={
            "Tr": 300.0,
            "eta0": 0.02,
            "c0": 0.5096,
            "c1": -13.4246,
            "c3": 36.2194,
            "c4": 5.4705e-8,
            "c5": 0.0032238,
        },
    )

    assert result.success
    assert result.parameters["c2"] == pytest.approx(0.0063295)


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


def test_thermal_reference_state_integrates_linear_temperature_expansivity():
    eos = ThermalReferenceStateEOS(
        BM3(V0=227.14, K0=65.4, K0_prime=2.7),
        Tr=298.0,
        alpha0=6.5e-5,
        alpha1=3.1e-9,
        dK_dT=-0.013,
        thermal_expansion_law="linear_temperature",
    )
    temperatures = np.array([298.0, 673.0, 1273.0])
    exponent = 6.5e-5 * (temperatures - 298.0) + 0.5 * 3.1e-9 * (
        temperatures**2 - 298.0**2
    )
    shifted_volumes = 227.14 * np.exp(exponent)
    shifted_moduli = 65.4 - 0.013 * (temperatures - 298.0)
    expected = np.array(
        [
            BM3(volume, modulus, 2.7).pressure(220.0)
            for volume, modulus in zip(shifted_volumes, shifted_moduli)
        ]
    )

    assert eos.configuration_values() == {
        "thermal_expansion_law": "linear_temperature",
        "reference_volume_law": "integrated_expansivity",
    }
    assert eos.parameter_values(include_reference=False)["alpha1"] == pytest.approx(
        3.1e-9
    )
    assert np.allclose(eos.pressure(220.0, temperatures), expected)
    assert np.allclose(eos.pressure(shifted_volumes, temperatures), 0.0, atol=1e-13)

    reconstructed = eos.with_parameters(alpha1=3.0e-9)
    assert reconstructed.alpha1 == pytest.approx(3.0e-9)
    assert reconstructed.thermal_expansion_law == "linear_temperature"


def test_thermal_reference_state_integrates_inverse_temperature_squared_term():
    eos = ThermalReferenceStateEOS(
        BM3(V0=1.0, K0=261.0, K0_prime=4.0),
        Tr=300.0,
        alpha0=1.982e-5,
        alpha1=0.818e-8,
        dK_dT=-0.0280,
        thermal_expansion_law="linear_temperature_inverse_square",
        alpha_inverse_square=0.474,
    )
    temperatures = np.array([300.0, 1000.0, 2000.0])
    exponent = (
        eos.alpha0 * (temperatures - eos.Tr)
        + 0.5 * eos.alpha1 * (temperatures**2 - eos.Tr**2)
        + eos.alpha_inverse_square * (1.0 / temperatures - 1.0 / eos.Tr)
    )
    shifted_volumes = eos.rt_eos.V0 * np.exp(exponent)
    shifted_moduli = eos.rt_eos.K0 + eos.dK_dT * (temperatures - eos.Tr)
    expected = np.array(
        [
            BM3(volume, modulus, 4.0).pressure(0.92)
            for volume, modulus in zip(shifted_volumes, shifted_moduli)
        ]
    )

    assert eos.configuration_values() == {
        "thermal_expansion_law": "linear_temperature_inverse_square",
        "reference_volume_law": "integrated_expansivity",
    }
    assert eos.parameter_values(include_reference=False)[
        "alpha_inverse_square"
    ] == pytest.approx(0.474)
    assert np.allclose(eos.pressure(0.92, temperatures), expected)
    assert np.allclose(eos.pressure(shifted_volumes, temperatures), 0.0, atol=1e-12)

    reconstructed = eos.with_parameters(alpha_inverse_square=0.5)
    assert reconstructed.alpha_inverse_square == pytest.approx(0.5)
    assert reconstructed.thermal_expansion_law == "linear_temperature_inverse_square"


def test_thermal_reference_state_supports_direct_linear_reference_volume():
    eos = ThermalReferenceStateEOS(
        BM2(V0=227.5, K0=64.81),
        Tr=298.0,
        alpha0=6.5e-5,
        dK_dT=-0.018,
        reference_volume_law="linear_temperature",
    )
    temperatures = np.array([298.0, 673.0, 973.0])
    shifted_volumes = 227.5 * (1.0 + 6.5e-5 * (temperatures - 298.0))
    shifted_moduli = 64.81 - 0.018 * (temperatures - 298.0)
    expected = np.array(
        [
            BM2(volume, modulus).pressure(220.0)
            for volume, modulus in zip(shifted_volumes, shifted_moduli)
        ]
    )

    assert eos.configuration_values() == {
        "thermal_expansion_law": "constant",
        "reference_volume_law": "linear_temperature",
    }
    assert np.allclose(eos.pressure(220.0, temperatures), expected)
    assert np.allclose(eos.pressure(shifted_volumes, temperatures), 0.0, atol=1e-13)

    reconstructed = eos.with_parameters(alpha0=6.4e-5)
    assert reconstructed.alpha0 == pytest.approx(6.4e-5)
    assert reconstructed.reference_volume_law == "linear_temperature"


def test_thermal_reference_state_supports_berman_reference_volume():
    eos = ThermalReferenceStateEOS(
        BM3(V0=328.4, K0=221.0, K0_prime=3.3),
        Tr=298.0,
        alpha0=1.94e-5,
        alpha1=5.73e-10,
        dK_dT=-0.008,
        thermal_expansion_law="linear_temperature",
        reference_volume_law="berman",
    )
    temperature = 2023.0
    delta = temperature - 298.0
    shifted_volume = 328.4 * (1.0 + 1.94e-5 * delta + 0.5 * 5.73e-10 * delta**2)
    shifted_modulus = 221.0 - 0.008 * delta
    expected = BM3(shifted_volume, shifted_modulus, 3.3).pressure(289.1)

    assert eos.pressure(289.1, temperature) == pytest.approx(expected)
    assert eos.pressure(shifted_volume, temperature) == pytest.approx(0.0, abs=1e-13)


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


def test_thermal_reference_state_validates_expansion_law_configuration():
    with pytest.raises(ValueError, match="thermal_expansion_law"):
        ThermalReferenceStateEOS(
            BM2(10.0, 20.0),
            300.0,
            1e-5,
            0.0,
            thermal_expansion_law="quadratic",
        )
    with pytest.raises(ValueError, match="alpha1 must be zero"):
        ThermalReferenceStateEOS(BM2(10.0, 20.0), 300.0, 1e-5, 0.0, 1e-9)
    with pytest.raises(ValueError, match="reference_volume_law"):
        ThermalReferenceStateEOS(
            BM2(10.0, 20.0),
            300.0,
            1e-5,
            0.0,
            reference_volume_law="quadratic",
        )
    with pytest.raises(ValueError, match="requires constant thermal expansion"):
        ThermalReferenceStateEOS(
            BM2(10.0, 20.0),
            300.0,
            1e-5,
            0.0,
            alpha1=1e-9,
            thermal_expansion_law="linear_temperature",
            reference_volume_law="linear_temperature",
        )
    with pytest.raises(ValueError, match="berman reference volume requires"):
        ThermalReferenceStateEOS(
            BM2(10.0, 20.0),
            300.0,
            1e-5,
            0.0,
            reference_volume_law="berman",
        )
