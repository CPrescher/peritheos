"""Keep ambient data, conditional refits and joint-fit projections distinct."""

import json

import numpy as np
import pytest

from scripts.brucite_ambient_comparison import (
    JOINT_REPORT,
    REPORT,
    curve_values,
    reproduce,
)
from scripts.brucite_mgd_comparison import data, pressure


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_joint_curves_reuse_previous_volume_objective_coefficients(report):
    joint = json.loads(JOINT_REPORT.read_text())
    assert (
        report["curves"]["empirical_joint"]["parameters"]
        == joint["fits"]["empirical_volume"]["parameters"]
    )
    assert (
        report["curves"]["mgd_joint"]["parameters"]
        == joint["fits"]["mgd_volume"]["parameters"]
    )
    assert report["ambient_observations"] == 11
    assert report["ambient_temperature_range_k"] == [299.3, 592.9]
    assert report["high_pressure_observations_already_in_joint_fits"] == {
        "300": 13,
        "473": 8,
        "673": 8,
        "873": 8,
    }


def test_conditional_mgd_fixes_mechanics_and_exposes_boundary(report):
    fit = report["conditional_ambient_mgd"]
    assert fit["solver_success"]
    assert fit["source_rows"] == list(range(1, 12))
    assert fit["free_parameters"] == ["V0_a3", "gamma0", "theta0_k"]
    for name, value in fit["fixed_parameters_from_joint_volume_fit"].items():
        assert value == report["curves"]["mgd_joint"]["parameters"][name]
        assert value == fit["parameters"][name]
    assert fit["active_bounds"] == {"theta0_k": -1}
    assert fit["parameters"]["theta0_k"] == pytest.approx(10)
    assert np.ptp(fit["multistart_sse"]) < 1e-11
    assert len(fit["formal_local_covariance"]["matrix"]) == 3


def test_ambient_residual_denominators_and_pressure_inversion(report):
    saved = json.loads(REPORT.read_text())
    dataset = data()
    temperature = dataset["temperature"][37:]
    observed = dataset["volume"][37:]
    for key, curve in report["curves"].items():
        predicted = curve_values(key, temperature, report)
        np.testing.assert_allclose(
            predicted - observed, curve["ambient_volume_residuals_a3"], atol=1e-12
        )
        assert np.sqrt(np.sum((predicted - observed) ** 2) / 11) == pytest.approx(
            curve["ambient_volume_rmse_a3"], abs=1e-12
        )
        assert curve["ambient_volume_rmse_a3"] == pytest.approx(
            saved["curves"][key]["ambient_volume_rmse_a3"], abs=1e-10
        )
        if "parameters" in curve:
            grid = np.array([299.3, 473, 592.9, 673, 873])
            predicted = curve_values(key, grid, report)
            model = "empirical" if key == "empirical_joint" else "mgd"
            np.testing.assert_allclose(
                pressure(model, predicted, grid, list(curve["parameters"].values())),
                0,
                atol=3e-11,
            )


def test_ambient_agreement_tradeoff_and_formal_high_temperature_reference(report):
    curves = report["curves"]
    assert curves["published"]["ambient_volume_rmse_a3"] == pytest.approx(
        0.0211674450, abs=1e-9
    )
    assert (
        curves["empirical_ambient"]["ambient_volume_rmse_a3"]
        < curves["empirical_joint"]["ambient_volume_rmse_a3"]
    )
    assert (
        curves["mgd_ambient_conditional"]["ambient_volume_rmse_a3"]
        < curves["mgd_joint"]["ambient_volume_rmse_a3"]
    )
    assert curves["published"]["formal_zero_pressure_volumes_a3"][
        "873"
    ] == pytest.approx(42.5159510616)
    assert curves["mgd_joint"]["formal_zero_pressure_volumes_a3"][
        "873"
    ] == pytest.approx(43.396349081, abs=1e-8)
