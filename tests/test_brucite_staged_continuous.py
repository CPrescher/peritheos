"""Verify fixed ambient expansion, linear modulus fits and continuous evaluation."""

import json

import numpy as np
import pytest

from scripts.brucite_replay_checks import check_report
from scripts.brucite_reproduce import bm3, observations
from scripts.brucite_staged_comparison import REPORT as STAGED_REPORT
from scripts.brucite_staged_continuous import (
    REPORT,
    moduli,
    reference_volume,
    reproduce,
    runtime_model,
)


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_fixed_ambient_law_and_regressed_intercepts(report):
    staged = json.loads(STAGED_REPORT.read_text(encoding="utf-8"))
    assert report["fixed_ambient_law"] == staged["stage1_ambient"]
    parameters = report["linear_moduli"]["parameters"]
    assert parameters == pytest.approx(
        {
            "K300_gpa": 40.411152632,
            "dK_dT_gpa_per_k": -0.02807347309,
            "Kprime300": 6.855997233,
            "dKprime_dT_per_k": 0.005172312941,
        },
        rel=1e-9,
    )
    # A fitted intercept is not silently forced through the 300 K staged value.
    assert parameters["K300_gpa"] != staged["isotherms"][0]["staged_refit"]["K0_gpa"]
    t = np.array([299.3, 300, 473, 573, 673, 873])
    eos = runtime_model(report)
    np.testing.assert_allclose(
        eos.pressure(reference_volume(t, report["fixed_ambient_law"]), t), 0, atol=1e-11
    )


def test_lines_solve_equal_weight_normal_equations(report):
    law = report["linear_moduli"]
    temperature = np.array(law["input_temperatures_k"])
    design = np.column_stack([np.ones(4), temperature - 300])
    residual = np.array(law["K_Kprime_residuals"])
    np.testing.assert_allclose(design.T @ residual, 0, atol=2e-10)
    for name in ["covariance_from_line_scatter", "covariance_propagated_from_stage2"]:
        covariance = np.array(law[name])
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-12)
        assert np.linalg.eigvalsh(covariance).min() > -1e-12
        assert covariance[0, 2] != 0  # Retain the K/Kprime cross-correlation.


def test_native_and_python_evaluation_between_measured_temperatures(report):
    eos = runtime_model(report)
    v, t = np.meshgrid([33, 35, 38, 41], [300, 400, 573, 773, 873])
    k0, kp = moduli(t, report["linear_moduli"]["parameters"])
    expected = bm3(v, reference_volume(t, report["fixed_ambient_law"]), k0, kp)
    np.testing.assert_allclose(eos.pressure(v, t), expected, atol=1e-11)
    native = eos._native
    del eos._native
    try:
        np.testing.assert_allclose(eos.pressure(v, t), expected, atol=1e-11)
    finally:
        eos._native = native
    for pressure, temperature in [(6, 300), (10, 400), (10, 573), (15, 773), (15, 873)]:
        volume = eos.volume(pressure, temperature)
        assert eos.pressure(volume, temperature) == pytest.approx(pressure, abs=1e-9)


def test_all37_residuals_preserve_rows_and_report_the_loss_of_slice_agreement(report):
    rows = observations("fukui_2003", "table4")
    eos = runtime_model(report)
    for result, count in zip(report["isotherms"], [13, 8, 8, 8]):
        selected = [
            row
            for row in rows
            if float(row["temperature_k"]) == result["temperature_k"]
        ]
        assert len(selected) == result["observations"] == count
        v = np.array([float(row["cell_volume_a3"]) for row in selected])
        p = np.array([float(row["pressure_gpa"]) for row in selected])
        np.testing.assert_allclose(
            eos.pressure(v, result["temperature_k"]) - p,
            result["continuous"]["pressure_residuals_gpa"],
            atol=1e-11,
        )
        np.testing.assert_allclose(
            eos.volume(p, result["temperature_k"]) - v,
            result["continuous"]["volume_residuals_a3"],
            atol=1e-9,
        )
        assert (
            result["continuous"]["volume_rmse_a3"]
            > result["staged_refit"]["volume_rmse_a3"]
        )
    assert report["all37_assessment"]["continuous"][
        "pressure_rmse_gpa"
    ] == pytest.approx(0.386160731, abs=1e-8)
    assert report["all37_assessment"]["continuous"]["volume_rmse_a3"] == pytest.approx(
        0.123518423, abs=1e-8
    )
    saved = json.loads(REPORT.read_text(encoding="utf-8"))
    check_report(report, saved)
