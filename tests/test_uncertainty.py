"""Tests for EOS parameter and prediction uncertainty propagation."""

import numpy as np
import pytest

from peritheos.eos import EosBase
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
from peritheos.fitting import fit_rt_eos, fit_thermal_eos
from peritheos.uncertainty import EOSUncertainty, ParameterUncertainty


class BoundaryEOS(EosBase):
    """Tiny test EOS whose parameter and state must remain non-negative."""

    def __init__(self, coefficient):
        coefficient = float(coefficient)
        if coefficient < 0.0:
            raise ValueError("coefficient must be non-negative")
        self.coefficient = coefficient
        self.V0 = 1.0

    def pressure(self, V):
        values = np.asarray(V, dtype=float)
        if np.any(values <= 0.0):
            raise ValueError("V must be positive")
        result = self.coefficient * values
        return float(result) if result.ndim == 0 else result

    def bulk_modulus(self, V):
        values = np.asarray(V, dtype=float)
        result = np.broadcast_to(self.coefficient, values.shape)
        return float(result) if result.ndim == 0 else result


class LockedParameterEOS(BoundaryEOS):
    """Accept only its nominal parameter to exercise failed sampling."""

    def __init__(self, coefficient):
        if float(coefficient) != 1.0:
            raise ValueError("coefficient is locked")
        super().__init__(coefficient)


class NonFiniteEOS(BoundaryEOS):
    def pressure(self, V):
        return np.asarray(V, dtype=float) * np.nan


class ShapeChangingEOS(BoundaryEOS):
    def pressure(self, V):
        result = super().pressure(V)
        if self.coefficient == 1.0:
            return result
        return np.atleast_1d(result)


def test_parameter_values_and_safe_reconstruction():
    eos = BM3(10.0, 120.0, 4.3)

    assert eos.parameter_values() == {
        "V0": 10.0,
        "K0": 120.0,
        "K0_prime": 4.3,
    }
    modified = eos.with_parameters(K0=130.0)
    assert modified.parameter_values()["K0"] == 130.0
    assert eos.K0 == 120.0


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
        Holzapfel(0.3414, 441.5, 3.9, 1.0, 6),
    ],
)
def test_all_room_temperature_models_can_reconstruct(eos):
    rebuilt = eos.with_parameters()

    assert type(rebuilt) is type(eos)
    assert rebuilt.parameter_values() == eos.parameter_values()


@pytest.mark.parametrize(
    "eos",
    [
        MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
        MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
        ThermalModifiedTait(
            ModifiedTait(1.0, 160.0, 4.0, -0.01),
            300.0,
            800.0,
            3.0e-5,
            2.0,
        ),
        Sokolova2016(
            Holzapfel(0.3414, 441.5, 3.9, 1.0, 6),
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
    ],
)
def test_all_thermal_models_can_reconstruct(eos):
    rebuilt = eos.with_parameters()

    assert type(rebuilt) is type(eos)
    assert rebuilt.parameter_values() == eos.parameter_values()


def test_nested_thermal_parameter_reconstruction():
    thermal = MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0)

    parameters = thermal.parameter_values()
    assert parameters["gamma0"] == 1.5
    assert parameters["rt_eos.K0"] == 160.0
    modified = thermal.with_parameters(gamma0=1.7, **{"rt_eos.K0": 170.0})
    assert modified.gamma0 == 1.7
    assert modified.rt_eos.K0 == 170.0
    assert thermal.rt_eos.K0 == 160.0


def test_covariance_derives_errors_and_correlation():
    uncertainty = ParameterUncertainty(
        parameter_names=("V0", "K0"),
        covariance=np.array([[0.04, 0.1], [0.1, 4.0]]),
    )

    assert uncertainty.standard_errors == {"V0": 0.2, "K0": 2.0}
    assert np.allclose(uncertainty.correlation, [[1.0, 0.25], [0.25, 1.0]])


def test_parameter_errors_and_correlation_build_covariance():
    uncertainty = ParameterUncertainty(
        parameter_errors={"V0": 0.2, "K0": 2.0},
        correlation=np.array([[1.0, -0.5], [-0.5, 1.0]]),
    )

    assert np.allclose(uncertainty.covariance, [[0.04, -0.2], [-0.2, 4.0]])
    assert not uncertainty.assumptions


def test_individual_errors_explicitly_assume_independence():
    uncertainty = ParameterUncertainty(parameter_errors={"K0": 2.0})

    assert np.array_equal(uncertainty.covariance, [[4.0]])
    assert "mutually independent" in uncertainty.assumptions[0]


def test_invalid_uncertainty_inputs():
    eos = BM3(10.0, 120.0, 4.3)
    with pytest.raises(ValueError, match="correlation requires"):
        EOSUncertainty(eos, correlation=np.eye(3))
    with pytest.raises(ValueError, match="positive semidefinite"):
        EOSUncertainty(
            eos,
            covariance=[[1.0, 2.0], [2.0, 1.0]],
            parameter_names=("V0", "K0"),
        )
    with pytest.raises(ValueError, match="Unknown parameters"):
        EOSUncertainty(eos, parameter_errors={"missing": 1.0})


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (
            {"covariance": [[1.0]], "parameter_errors": {"K0": 1.0}},
            "cannot be combined",
        ),
        ({"covariance": [[1.0]]}, "parameter_names are required"),
        ({"covariance": [], "parameter_names": ()}, "At least one"),
        (
            {"covariance": [[np.nan]], "parameter_names": ("K0",)},
            "finite values",
        ),
        (
            {
                "covariance": [[1.0, 0.2], [0.1, 1.0]],
                "parameter_names": ("V0", "K0"),
            },
            "symmetric",
        ),
        ({"parameter_errors": {}}, "At least one"),
        (
            {"parameter_errors": {"K0": 1.0}, "parameter_names": ("V0",)},
            "must match",
        ),
        ({"parameter_errors": {"K0": 0.0}}, "greater than zero"),
        (
            {
                "parameter_errors": {"V0": 0.1, "K0": 1.0},
                "correlation": [[1.0, 0.0], [0.0, 0.9]],
            },
            "ones on its diagonal",
        ),
        (
            {
                "parameter_errors": {"V0": 0.1, "K0": 1.0},
                "correlation": [[1.0, 1.1], [1.1, 1.0]],
            },
            "positive semidefinite|lie in",
        ),
        (
            {"covariance": np.eye(2), "parameter_names": ("K0", "K0")},
            "unique",
        ),
        ({}, "Provide covariance"),
    ],
)
def test_parameter_uncertainty_validation(kwargs, match):
    with pytest.raises(ValueError, match=match):
        ParameterUncertainty(**kwargs)


def test_eos_uncertainty_infers_complete_covariance_order():
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), covariance=np.diag([0.01, 4.0]))

    assert uncertainty.parameter_names == ("V0", "K0")
    assert uncertainty.standard_errors == {"V0": 0.1, "K0": 2.0}
    assert np.array_equal(uncertainty.correlation, np.eye(2))


def test_eos_uncertainty_rejects_invalid_eos_and_ambiguous_covariance():
    with pytest.raises(TypeError, match="equation of state"):
        EOSUncertainty(object(), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="parameter_names are required"):
        EOSUncertainty(BM2(10.0, 120.0), covariance=[[1.0]])


def test_linear_pressure_uncertainty_matches_analytic_bm2_result():
    eos = BM2(10.0, 120.0)
    volumes = np.array([8.5, 9.0, 9.5])
    uncertainty = EOSUncertainty(eos, parameter_errors={"K0": 2.0})

    prediction = uncertainty.pressure(volumes, full_covariance=True)
    derivative = np.asarray(eos.pressure(volumes)) / eos.K0
    expected_errors = np.abs(derivative) * 2.0

    assert np.allclose(prediction.value, eos.pressure(volumes))
    assert np.allclose(prediction.standard_error, expected_errors, rtol=1.0e-7)
    assert np.allclose(prediction.covariance, 4.0 * np.outer(derivative, derivative))
    assert prediction.method == "linear"


def test_parameter_correlation_changes_prediction_uncertainty():
    eos = BM3(10.0, 120.0, 4.3)
    errors = {"V0": 0.02, "K0": 2.0}
    positive = EOSUncertainty(
        eos,
        parameter_errors=errors,
        correlation=[[1.0, 0.9], [0.9, 1.0]],
    )
    negative = EOSUncertainty(
        eos,
        parameter_errors=errors,
        correlation=[[1.0, -0.9], [-0.9, 1.0]],
    )

    assert positive.pressure(9.0).standard_error > negative.pressure(9.0).standard_error


def test_volume_measurement_error_is_propagated_to_pressure():
    eos = BM2(10.0, 120.0)
    uncertainty = EOSUncertainty(eos, covariance=[[0.0]], parameter_names=("K0",))
    volume = 9.0
    volume_sigma = 0.01

    prediction = uncertainty.pressure(volume, volume_sigma=volume_sigma)
    expected = eos.bulk_modulus(volume) / volume * volume_sigma

    assert np.isclose(prediction.standard_error, expected, rtol=1.0e-6)
    assert any("state-variable" in item for item in prediction.assumptions)


def test_monte_carlo_is_reproducible_and_agrees_for_linear_parameter():
    eos = BM2(10.0, 120.0)
    uncertainty = EOSUncertainty(eos, parameter_errors={"K0": 2.0})
    linear = uncertainty.pressure(9.0)
    first = uncertainty.pressure(
        9.0, method="monte_carlo", sample_count=4000, random_state=12
    )
    second = uncertainty.pressure(
        9.0, method="monte_carlo", sample_count=4000, random_state=12
    )

    assert first == second
    assert np.isclose(first.standard_error, linear.standard_error, rtol=0.05)
    assert first.lower < first.value < first.upper


def test_monte_carlo_rejects_invalid_eos_parameter_samples():
    eos = BM2(10.0, 1.0)
    uncertainty = EOSUncertainty(eos, parameter_errors={"K0": 2.0})

    prediction = uncertainty.pressure(
        9.0, method="monte_carlo", sample_count=100, random_state=3
    )

    assert prediction.rejected_fraction > 0.0


def test_monte_carlo_state_sampling_and_output_covariance():
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    prediction = uncertainty.pressure(
        np.array([8.5, 9.0]),
        volume_sigma=np.array([0.01, 0.02]),
        method="monte_carlo",
        full_covariance=True,
        sample_count=200,
        random_state=9,
    )

    assert prediction.covariance.shape == (2, 2)
    assert np.all(np.asarray(prediction.standard_error) > 0.0)
    assert any(
        "state-variable errors sampled" in item for item in prediction.assumptions
    )


@pytest.mark.parametrize("sample_count", [True, 2.5])
def test_monte_carlo_requires_integer_sample_count(sample_count):
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="integer"):
        uncertainty.pressure(9.0, method="monte_carlo", sample_count=sample_count)


def test_monte_carlo_rejects_too_few_or_unobtainable_samples():
    normal = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="at least two"):
        normal.pressure(9.0, method="monte_carlo", sample_count=1)

    locked = EOSUncertainty(
        LockedParameterEOS(1.0), parameter_errors={"coefficient": 0.1}
    )
    with pytest.raises(ArithmeticError, match="enough valid"):
        locked.pressure(1.0, method="monte_carlo", sample_count=2, random_state=2)


def test_one_sided_parameter_and_state_differences_are_supported():
    uncertainty = EOSUncertainty(
        BoundaryEOS(0.0), parameter_errors={"coefficient": 0.1}
    )

    parameter_boundary = uncertainty.pressure(2.0)
    state_boundary = uncertainty.pressure(
        5.0e-7,
        volume_sigma=1.0e-7,
        relative_step=1.0e-6,
    )

    assert np.isclose(parameter_boundary.standard_error, 0.2)
    assert state_boundary.standard_error > 0.0


def test_evaluation_errors_are_descriptive():
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="Unknown public"):
        uncertainty.evaluate("missing", 9.0)
    with pytest.raises(ValueError, match="Unknown public"):
        uncertainty.evaluate("_own_parameter_names")
    with pytest.raises(ArithmeticError, match="non-finite"):
        EOSUncertainty(
            NonFiniteEOS(1.0), parameter_errors={"coefficient": 0.1}
        ).pressure(1.0)
    with pytest.raises(ArithmeticError, match="output shape changed"):
        EOSUncertainty(
            ShapeChangingEOS(1.0), parameter_errors={"coefficient": 0.1}
        ).pressure(1.0)


@pytest.mark.parametrize(
    "options,match",
    [
        ({"confidence": 0.0}, "confidence"),
        ({"confidence": np.nan}, "confidence"),
        ({"relative_step": 0.0}, "relative_step"),
        ({"relative_step": np.inf}, "relative_step"),
        ({"method": "bootstrap"}, "method must be"),
        ({"argument_sigmas": {1: 0.1}}, "Invalid argument"),
        ({"argument_sigmas": {0: -0.1}}, "state uncertainties"),
        ({"argument_sigmas": {0: np.ones(2)}}, "broadcast"),
    ],
)
def test_evaluate_option_validation(options, match):
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match=match):
        uncertainty.evaluate("pressure", 9.0, **options)


def test_monte_carlo_state_uncertainty_validation():
    uncertainty = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="broadcast"):
        uncertainty.pressure(
            9.0,
            volume_sigma=np.ones(2),
            method="monte_carlo",
            sample_count=10,
        )
    with pytest.raises(ValueError, match="greater than zero"):
        uncertainty.pressure(
            9.0,
            volume_sigma=np.nan,
            method="monte_carlo",
            sample_count=10,
        )


def test_fit_result_creates_eos_uncertainty():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.5, 20)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.05,
        absolute_sigma=True,
    )

    uncertainty = result.eos_uncertainty()
    prediction = uncertainty.pressure(9.0)

    assert uncertainty.parameter_names == result.free_parameters
    assert np.allclose(uncertainty.covariance, result.covariance)
    assert prediction.standard_error > 0.0


def test_thermal_fit_can_add_independent_reference_eos_uncertainty():
    rt_eos = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(rt_eos, 300.0, 800.0, 1.6, 1.0, 2.0)
    volumes = np.repeat(np.array([0.8, 0.9, 1.0]), 5)
    temperatures = np.tile(np.linspace(500.0, 2000.0, 5), 3)
    pressures = expected.pressure(volumes, temperatures)
    result = fit_thermal_eos(
        MieGruneisenEinstein,
        rt_eos,
        volumes,
        temperatures,
        pressures,
        initial={"gamma0": 1.5, "q": 0.9},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
        pressure_sigma=0.01,
        absolute_sigma=True,
    )
    reference = EOSUncertainty(rt_eos, parameter_errors={"K0": 1.0})

    with pytest.raises(ValueError, match="assume_blocks_independent"):
        result.eos_uncertainty(additional=reference)
    combined = result.eos_uncertainty(
        additional=reference, assume_blocks_independent=True
    )

    assert "rt_eos.K0" in combined.parameter_names
    assert any(
        "covariance blocks" in item
        for item in combined.parameter_uncertainty.assumptions
    )
    assert combined.pressure(0.9, 1200.0).standard_error > 0.0

    wrong_reference = EOSUncertainty(BM3(1.0, 170.0, 4.0), parameter_errors={"K0": 1.0})
    with pytest.raises(ValueError, match="reference EOS"):
        result.eos_uncertainty(
            additional=wrong_reference, assume_blocks_independent=True
        )


def test_additional_uncertainty_requires_thermal_fit():
    volumes = np.linspace(8.0, 10.0, 8)
    model = BM2(10.0, 120.0)
    result = fit_rt_eos(
        BM2,
        volumes,
        model.pressure(volumes),
        initial={"K0": 110.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.1,
        absolute_sigma=True,
    )
    additional = EOSUncertainty(model, parameter_errors={"K0": 1.0})

    with pytest.raises(TypeError, match="thermal fitted EOS"):
        result.eos_uncertainty(additional=additional, assume_blocks_independent=True)


def test_uncertainty_convenience_methods_for_rt_and_thermal_eos():
    room = EOSUncertainty(BM2(10.0, 120.0), parameter_errors={"K0": 1.0})
    assert room.volume(20.0, pressure_sigma=0.1).standard_error > 0.0
    assert room.bulk_modulus(9.0, volume_sigma=0.01).standard_error > 0.0
    with pytest.raises(ValueError, match="only valid"):
        room.pressure(9.0, T=300.0)
    with pytest.raises(ValueError, match="only valid"):
        room.volume(20.0, temperature_sigma=1.0)
    with pytest.raises(ValueError, match="only valid"):
        room.bulk_modulus(9.0, T=300.0)

    thermal_eos = MieGruneisenEinstein(
        BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0
    )
    thermal = EOSUncertainty(thermal_eos, parameter_errors={"gamma0": 0.05})
    pressure = thermal_eos.pressure(0.9, 1200.0)
    assert thermal.pressure(0.9, 1200.0).standard_error > 0.0
    assert thermal.volume(pressure, 1200.0).standard_error > 0.0
    assert thermal.bulk_modulus(0.9, 1200.0).standard_error > 0.0
    for operation, argument in [
        (thermal.pressure, 0.9),
        (thermal.volume, pressure),
        (thermal.bulk_modulus, 0.9),
    ]:
        with pytest.raises(ValueError, match="Temperature is required"):
            operation(argument)
