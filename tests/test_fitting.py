"""Tests for P-V and P-V-T EOS fitting."""

import json

import numpy as np
import pytest

import peritheos.fitting as fitting_module
from peritheos import _rust
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
from peritheos.fitting import fit_joint_eos, fit_rt_eos, fit_thermal_eos


@pytest.mark.parametrize(
    "expected",
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
def test_every_builtin_rt_model_fits_end_to_end_natively(expected, monkeypatch):
    parameters = expected.parameter_values(include_reference=False)
    volumes = np.linspace(0.8 * expected.V0, expected.V0, 12)
    pressures = expected.pressure(volumes)

    def callback_solver_is_forbidden(*args, **kwargs):
        raise AssertionError("built-in fitting must not use a Python callback")

    monkeypatch.setattr(_rust, "fit_least_squares", callback_solver_is_forbidden)
    result = fit_rt_eos(
        type(expected),
        volumes,
        pressures,
        initial={"K0": 0.9 * parameters["K0"]},
        fixed={name: value for name, value in parameters.items() if name != "K0"},
    )

    assert result.success
    assert result.parameters["K0"] == pytest.approx(parameters["K0"], rel=1.0e-7)


@pytest.mark.parametrize(
    ("expected", "free_name", "initial_factor"),
    [
        (
            MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
            "gamma0",
            0.8,
        ),
        (
            MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
            "gamma0",
            0.8,
        ),
        (
            ThermalModifiedTait(
                ModifiedTait(1.0, 160.0, 4.0, -0.01),
                298.15,
                700.0,
                2.5e-5,
                2.0,
            ),
            "alpha0",
            0.8,
        ),
        (
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
            "delta",
            0.8,
        ),
    ],
)
def test_every_builtin_thermal_model_fits_end_to_end_natively(
    expected, free_name, initial_factor, monkeypatch
):
    parameters = expected.parameter_values(include_reference=False)
    volumes = np.repeat(expected.rt_eos.V0 * np.array([0.82, 0.9, 0.98]), 4)
    temperatures = np.tile(np.linspace(500.0, 2200.0, 4), 3)
    pressures = expected.pressure(volumes, temperatures)

    def callback_solver_is_forbidden(*args, **kwargs):
        raise AssertionError("built-in fitting must not use a Python callback")

    monkeypatch.setattr(_rust, "fit_least_squares", callback_solver_is_forbidden)
    result = fit_thermal_eos(
        type(expected),
        expected.rt_eos,
        volumes,
        temperatures,
        pressures,
        initial={free_name: initial_factor * parameters[free_name]},
        fixed={name: value for name, value in parameters.items() if name != free_name},
        max_nfev=500,
    )

    assert result.success
    assert result.parameters[free_name] == pytest.approx(
        parameters[free_name], rel=1.0e-6, abs=1.0e-10
    )


def test_custom_subclass_retains_python_residual_callback(monkeypatch):
    class ShiftedBM3(BM3):
        def pressure(self, V):
            return np.asarray(super().pressure(V)) + 2.0

    expected = ShiftedBM3(10.0, 120.0, 4.0)
    volumes = np.linspace(8.0, 10.0, 12)
    native_fast_path = _rust.fit_rt_eos_native

    def native_fast_path_is_forbidden(*args, **kwargs):
        raise AssertionError("custom subclasses must retain Python evaluation")

    monkeypatch.setattr(_rust, "fit_rt_eos_native", native_fast_path_is_forbidden)
    result = fit_rt_eos(
        ShiftedBM3,
        volumes,
        expected.pressure(volumes),
        initial={"K0": 100.0},
        fixed={"V0": 10.0, "K0_prime": 4.0},
    )
    monkeypatch.setattr(_rust, "fit_rt_eos_native", native_fast_path)

    assert result.success
    assert result.parameters["K0"] == pytest.approx(120.0)


def test_builtin_fit_never_calls_python_pressure_during_solve(monkeypatch):
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.0, 20)
    pressures = expected.pressure(volumes)
    python_evaluations = 0
    original_pressure = BM3.pressure

    def counting_pressure(self, V):
        nonlocal python_evaluations
        python_evaluations += 1
        return original_pressure(self, V)

    monkeypatch.setattr(BM3, "pressure", counting_pressure)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
    )

    assert result.success
    assert python_evaluations == 0


def test_native_fit_uses_rust_parameter_covariance(monkeypatch):
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.0, 20)

    def python_covariance_is_forbidden(*args, **kwargs):
        raise AssertionError("native fits must use the Rust covariance kernel")

    monkeypatch.setattr(
        fitting_module, "_parameter_covariance", python_covariance_is_forbidden
    )
    result = fit_rt_eos(
        BM3,
        volumes,
        expected.pressure(volumes),
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
    )

    assert result.success
    assert result.covariance.shape == (2, 2)
    assert np.all(np.isfinite(result.covariance))


def test_callable_loss_retains_scipy_compatibility_path(monkeypatch):
    expected = BM3(10.0, 120.0, 4.0)
    volumes = np.linspace(8.0, 10.0, 12)

    def linear_callable(squared_residuals):
        return np.vstack(
            [
                squared_residuals,
                np.ones_like(squared_residuals),
                np.zeros_like(squared_residuals),
            ]
        )

    def native_fast_path_is_forbidden(*args, **kwargs):
        raise AssertionError("callable losses must retain SciPy compatibility")

    monkeypatch.setattr(_rust, "fit_rt_eos_native", native_fast_path_is_forbidden)
    result = fit_rt_eos(
        BM3,
        volumes,
        expected.pressure(volumes),
        initial={"K0": 100.0},
        fixed={"V0": 10.0, "K0_prime": 4.0},
        loss=linear_callable,
    )

    assert result.success
    assert result.parameters["K0"] == pytest.approx(120.0)
    assert result.loss.endswith("linear_callable")


def test_fit_rt_eos_recovers_synthetic_parameters():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(7.5, 11.0, 30)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"V0": 9.8, "K0": 110.0, "K0_prime": 4.0},
        bounds={"V0": (8.0, 12.0), "K0": (50.0, 200.0)},
    )

    assert result.success
    assert np.isclose(result.parameters["V0"], expected.V0)
    assert np.isclose(result.parameters["K0"], expected.K0)
    assert np.isclose(result.parameters["K0_prime"], expected.K0_prime)
    assert result.covariance.shape == (3, 3)
    assert result.correlation.shape == (3, 3)
    assert result.degrees_of_freedom == volumes.size - 3


def test_fit_with_fixed_parameter_and_absolute_uncertainties():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.5, 20)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"K0": 100.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        sigma=0.05,
        absolute_sigma=True,
    )

    assert result.parameters["V0"] == 10.0
    assert result.standard_errors["V0"] == 0.0
    assert np.all(np.diag(result.covariance) > 0.0)


def test_fit_result_summary_and_json_export(tmp_path):
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.0, 12)
    result = fit_rt_eos(
        BM3,
        volumes,
        expected.pressure(volumes),
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.1,
        absolute_sigma=True,
    )

    summary = result.summary(precision=5)
    assert "FitResult (BM3)" in summary
    assert "K0_prime" in summary
    assert "fixed" in summary
    assert "chi-square" in summary
    assert "function evaluations" in summary

    payload = result.to_dict()
    assert payload["schema_version"] == 1
    assert payload["model"]["class"] == "BM3"
    assert payload["model"]["parameters"] == pytest.approx(result.parameters)
    assert payload["free_parameters"] == list(result.free_parameters)
    assert np.allclose(payload["covariance"], result.covariance)
    assert payload["adjusted_temperature"] is None
    assert payload["solver"]["loss"] == "linear"
    assert payload["solver"]["nfev"] > 0

    output_path = tmp_path / "fit-result.json"
    serialized = result.to_json(output_path, indent=2)
    decoded = json.loads(serialized)
    assert decoded == payload
    assert json.loads(output_path.read_text()) == payload

    with pytest.raises(ValueError, match="precision"):
        result.summary(precision=True)


def test_fit_result_json_replaces_nonfinite_diagnostics():
    expected = BM3(10.0, 120.0, 4.0)
    result = fit_rt_eos(
        BM3,
        volume=[9.0],
        pressure=[expected.pressure(9.0)],
        initial={"K0": 110.0},
        fixed={"V0": 10.0, "K0_prime": 4.0},
    )

    assert result.degrees_of_freedom == 0
    assert result.to_dict()["diagnostics"]["reduced_chi_square"] is None
    assert json.loads(result.to_json())["diagnostics"]["reduced_chi_square"] is None


def test_bm2_fit_matches_closed_form_weighted_least_squares():
    """Benchmark the optimizer against an independently solvable linear fit."""
    volumes = np.array([8.0, 8.5, 9.0, 9.5, 10.0])
    observed_pressure = np.array([35.2, 23.8, 14.5, 6.3, 0.2])
    pressure_sigma = np.array([0.4, 0.3, 0.25, 0.2, 0.2])
    unit_pressure = np.asarray(BM3(10.0, 1.0, 4.0).pressure(volumes))
    weights = pressure_sigma**-2
    expected_k0 = np.sum(weights * unit_pressure * observed_pressure) / np.sum(
        weights * unit_pressure**2
    )
    expected_variance = 1.0 / np.sum(weights * unit_pressure**2)

    result = fit_rt_eos(
        BM3,
        volumes,
        observed_pressure,
        initial={"K0": 90.0},
        fixed={"V0": 10.0, "K0_prime": 4.0},
        pressure_sigma=pressure_sigma,
        absolute_sigma=True,
    )

    assert np.isclose(result.parameters["K0"], expected_k0, rtol=1.0e-8)
    assert np.isclose(result.covariance[0, 0], expected_variance, rtol=1.0e-6)


def test_rank_deficient_fit_returns_finite_pseudoinverse_covariance():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.full(12, 9.0)
    pressures = expected.pressure(volumes)

    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"K0": 100.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.1,
        absolute_sigma=True,
    )

    assert result.success
    assert result.model.pressure(9.0) == pytest.approx(pressures[0])
    assert np.all(np.isfinite(result.covariance))
    assert np.linalg.eigvalsh(result.covariance).min() >= -1.0e-12


def test_fit_rt_eos_handles_pressure_and_volume_uncertainties():
    expected = BM3(10.0, 120.0, 4.3)
    true_volumes = np.linspace(8.0, 10.5, 20)
    measured_volumes = true_volumes + 0.008 * np.sin(np.arange(true_volumes.size))
    pressures = expected.pressure(true_volumes)

    result = fit_rt_eos(
        BM3,
        measured_volumes,
        pressures,
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.002,
        volume_sigma=0.01,
        absolute_sigma=True,
    )

    assert result.success
    assert np.isclose(result.parameters["K0"], expected.K0, rtol=0.01)
    assert np.isclose(result.parameters["K0_prime"], expected.K0_prime, rtol=0.02)
    assert np.any(np.abs(result.volume_corrections) > 0.0)
    assert result.adjusted_temperature is None
    assert result.temperature_corrections is None
    assert result.weighted_residuals.size == 2 * pressures.size
    assert result.degrees_of_freedom == pressures.size - 2


def test_large_latent_volume_fit_uses_structured_native_path():
    expected = BM3(10.0, 120.0, 4.3)
    true_volumes = np.linspace(8.0, 10.5, 1000)
    measured_volumes = true_volumes + 0.002 * np.sin(np.arange(true_volumes.size))
    pressures = expected.pressure(true_volumes)

    result = fit_rt_eos(
        BM3,
        measured_volumes,
        pressures,
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.01,
        volume_sigma=0.005,
        absolute_sigma=True,
        max_nfev=200,
    )

    assert result.success
    assert result.parameters["K0"] == pytest.approx(expected.K0, rel=1.0e-3)
    assert result.parameters["K0_prime"] == pytest.approx(expected.K0_prime, rel=1.0e-3)
    assert result.weighted_residuals.size == 2000
    assert result.covariance.shape == (2, 2)


def test_observation_covariance_matches_independent_rt_uncertainties():
    expected = BM3(10.0, 120.0, 4.3)
    true_volumes = np.linspace(8.0, 10.5, 20)
    measured_volumes = true_volumes + 0.008 * np.sin(np.arange(true_volumes.size))
    pressures = expected.pressure(true_volumes)
    common = {
        "initial": {"K0": 110.0, "K0_prime": 4.0},
        "fixed": {"V0": 10.0},
        "absolute_sigma": True,
    }

    independent = fit_rt_eos(
        BM3,
        measured_volumes,
        pressures,
        pressure_sigma=0.002,
        volume_sigma=0.01,
        **common,
    )
    covariance = np.diag([0.002**2, 0.01**2])
    correlated_api = fit_rt_eos(
        BM3,
        measured_volumes,
        pressures,
        observation_covariance=covariance,
        **common,
    )

    assert correlated_api.parameters == pytest.approx(independent.parameters)
    assert correlated_api.adjusted_volume == pytest.approx(independent.adjusted_volume)
    assert correlated_api.covariance == pytest.approx(independent.covariance)


def test_fit_thermal_eos_recovers_gamma():
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
        initial={"gamma0": 1.3, "q": 0.8},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
    )

    assert result.success
    assert np.isclose(result.parameters["gamma0"], 1.6)
    assert np.isclose(result.parameters["q"], 1.0)


def test_joint_fit_recovers_reference_and_thermal_parameters():
    reference = BM3(1.0, 160.0, 4.2)
    expected = MieGruneisenEinstein(reference, 300.0, 800.0, 1.6, 1.0, 2.0)
    volumes = np.tile(np.linspace(0.78, 1.02, 12), 4)
    temperatures = np.repeat(np.array([300.0, 700.0, 1300.0, 2000.0]), 12)
    pressures = expected.pressure(volumes, temperatures)

    result = fit_joint_eos(
        MieGruneisenEinstein,
        BM3,
        volumes,
        temperatures,
        pressures,
        initial={
            "rt_eos.V0": 0.99,
            "rt_eos.K0": 150.0,
            "rt_eos.K0_prime": 4.0,
            "gamma0": 1.4,
            "q": 0.8,
        },
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
        bounds={
            "rt_eos.V0": (0.9, 1.1),
            "rt_eos.K0": (100.0, 220.0),
            "gamma0": (0.5, 3.0),
        },
        pressure_sigma=0.01,
        absolute_sigma=True,
    )

    assert result.success
    assert result.parameters["rt_eos.V0"] == pytest.approx(1.0)
    assert result.parameters["rt_eos.K0"] == pytest.approx(160.0)
    assert result.parameters["rt_eos.K0_prime"] == pytest.approx(4.2)
    assert result.parameters["gamma0"] == pytest.approx(1.6)
    assert result.parameters["q"] == pytest.approx(1.0)
    assert result.model.rt_eos.K0 == pytest.approx(160.0)
    assert result.covariance.shape == (5, 5)
    assert result.eos_uncertainty().parameter_names == result.free_parameters


def test_fit_thermal_eos_handles_uncertainties_in_all_observables():
    rt_eos = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(rt_eos, 300.0, 800.0, 1.6, 1.0, 2.0)
    true_volumes = np.repeat(np.array([0.8, 0.9, 1.0]), 5)
    true_temperatures = np.tile(np.linspace(500.0, 2000.0, 5), 3)
    indices = np.arange(true_volumes.size)
    measured_volumes = true_volumes + 0.0003 * np.sin(indices)
    measured_temperatures = true_temperatures + 1.5 * np.cos(indices)
    pressures = expected.pressure(true_volumes, true_temperatures)

    result = fit_thermal_eos(
        MieGruneisenEinstein,
        rt_eos,
        measured_volumes,
        measured_temperatures,
        pressures,
        initial={"gamma0": 1.4, "q": 0.8},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
        pressure_sigma=0.001,
        volume_sigma=0.0005,
        temperature_sigma=2.0,
        absolute_sigma=True,
    )

    assert result.success
    assert np.isclose(result.parameters["gamma0"], 1.6, rtol=0.01)
    # The perturbed observables shift this weakly constrained parameter by
    # about 4%; allow for small optimizer differences across SciPy versions.
    assert np.isclose(result.parameters["q"], 1.0, rtol=0.05)
    assert np.any(np.abs(result.volume_corrections) > 0.0)
    assert np.any(np.abs(result.temperature_corrections) > 0.0)
    assert result.adjusted_temperature is not None
    assert result.weighted_residuals.size == 3 * pressures.size
    assert result.degrees_of_freedom == pressures.size - 2


def test_thermal_fit_accepts_per_observation_covariance():
    rt_eos = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(rt_eos, 300.0, 800.0, 1.6, 1.0, 2.0)
    volumes = np.repeat(np.array([0.82, 0.91, 1.0]), 4)
    temperatures = np.tile(np.linspace(500.0, 1800.0, 4), 3)
    pressures = expected.pressure(volumes, temperatures)
    standard_deviations = np.array([0.002, 0.0005, 2.0])
    correlation = np.array([[1.0, 0.3, -0.1], [0.3, 1.0, 0.2], [-0.1, 0.2, 1.0]])
    covariance = correlation * np.outer(standard_deviations, standard_deviations)
    per_observation_covariance = np.broadcast_to(
        covariance, (pressures.size, 3, 3)
    ).copy()

    result = fit_thermal_eos(
        MieGruneisenEinstein,
        rt_eos,
        volumes,
        temperatures,
        pressures,
        initial={"gamma0": 1.4, "q": 0.8},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
        observation_covariance=per_observation_covariance,
        absolute_sigma=True,
    )

    assert result.success
    assert result.parameters["gamma0"] == pytest.approx(1.6, rel=0.01)
    assert result.parameters["q"] == pytest.approx(1.0, rel=0.05)
    assert result.weighted_residuals.size == 3 * pressures.size
    assert result.degrees_of_freedom == pressures.size - 2
    raw_residuals = np.column_stack(
        [
            result.residuals,
            result.volume_corrections,
            result.temperature_corrections,
        ]
    )
    cholesky = np.linalg.cholesky(per_observation_covariance)
    expected_weighted = np.linalg.solve(cholesky, raw_residuals[..., None])[..., 0]
    assert result.weighted_residuals == pytest.approx(expected_weighted.T.ravel())


def test_robust_loss_reduces_outlier_bias():
    expected = BM3(10.0, 120.0, 4.0)
    volumes = np.linspace(8.0, 10.0, 25)
    pressures = np.asarray(expected.pressure(volumes)).copy()
    pressures[0] += 20.0
    common = {
        "initial": {"K0": 100.0},
        "fixed": {"V0": 10.0, "K0_prime": 4.0},
        "pressure_sigma": 0.1,
    }

    linear = fit_rt_eos(BM3, volumes, pressures, **common)
    robust = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        loss="soft_l1",
        f_scale=1.0,
        max_nfev=1000,
        **common,
    )

    assert robust.success
    assert robust.loss == "soft_l1"
    assert robust.f_scale == 1.0
    assert robust.max_nfev == 1000
    assert robust.nfev > 0
    assert np.isfinite(robust.cost)
    assert np.isfinite(robust.optimality)
    assert abs(robust.parameters["K0"] - expected.K0) < abs(
        linear.parameters["K0"] - expected.K0
    )


def test_fit_input_validation():
    volumes = np.linspace(8.0, 10.0, 5)
    pressures = BM3(10.0, 120.0, 4.0).pressure(volumes)
    with pytest.raises(ValueError, match="both initial and fixed"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"V0": 10.0},
            fixed={"V0": 10.0, "K0": 120.0, "K0_prime": 4.0},
        )
    with pytest.raises(ValueError, match="sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            sigma=0.0,
        )
    with pytest.raises(ValueError, match="either pressure_sigma or sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            pressure_sigma=0.1,
            sigma=0.1,
        )
    with pytest.raises(ValueError, match="volume_sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            volume_sigma=-0.1,
        )
    with pytest.raises(ValueError, match="cannot be combined"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            pressure_sigma=0.1,
            observation_covariance=np.eye(2),
        )
    with pytest.raises(ValueError, match="symmetric"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            observation_covariance=[[1.0, 0.2], [0.1, 1.0]],
        )
    with pytest.raises(ValueError, match="positive definite"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            observation_covariance=[[1.0, 2.0], [2.0, 1.0]],
        )
    with pytest.raises(ValueError, match="loss must be"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            loss="invalid",
        )
    with pytest.raises(ValueError, match="f_scale"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            f_scale=0.0,
        )
    with pytest.raises(ValueError, match="max_nfev"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            max_nfev=True,
        )


def test_joint_fit_requires_reference_parameters():
    reference = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(reference, 300.0, 800.0, 1.6, 1.0, 2.0)
    volumes = np.array([0.8, 0.9, 1.0])
    temperatures = np.array([500.0, 1000.0, 1500.0])

    with pytest.raises(ValueError, match=r"rt_eos\.\*"):
        fit_joint_eos(
            MieGruneisenEinstein,
            BM3,
            volumes,
            temperatures,
            expected.pressure(volumes, temperatures),
            initial={"gamma0": 1.4},
            fixed={"Tr": 300.0, "theta0": 800.0, "q": 1.0, "n": 2.0},
        )
