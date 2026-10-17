"""Check conditional V0 fits against the independent runtime equation/derivative."""

import json

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from scripts.brucite_replay_checks import check_report
from scripts.brucite_reproduce import observations, reference_volume
from scripts.brucite_v0_audit import MODES, REPORT, reproduce


@pytest.fixture(scope="module")
def report():
    return reproduce()


@pytest.mark.parametrize(
    "temperature,k0,kp,count",
    [
        (300, 41.8, 6.6, 13),
        (473, 32.6, 8.3, 8),
        (673, 28.8, 9.0, 8),
        (873, 23.9, 10.0, 8),
    ],
)
def test_fixed_moduli_v0_optimum_and_covariance_against_native_equation(
    report, temperature, k0, kp, count
):
    result = next(r for r in report["isotherms"] if r["temperature_k"] == temperature)
    assert result["fixed_published_moduli"] == {"K0_gpa": k0, "K0_prime": kp}
    rows = [
        r
        for r in observations("fukui_2003", "table4")
        if float(r["temperature_k"]) == temperature
    ]
    assert len(rows) == count == result["observations"]
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    fit = result["fits"][MODES[0]]
    model = BM3(fit["V0_a3"], k0, kp)
    residual = model.pressure(v) - p
    np.testing.assert_allclose(residual, fit["pressure_residuals_gpa"], atol=1e-12)
    # For P=f(V0/V), dP/dV0 = K(V)/V0. The least-squares normal equation
    # and residual-scaled covariance provide independent checks of the solver.
    derivative = model.bulk_modulus(v) / fit["V0_a3"]
    assert abs(residual @ derivative) < 1e-5
    covariance = (residual @ residual) / (count - 1) / (derivative @ derivative)
    assert fit["conditional_standard_error_a3"] == pytest.approx(
        np.sqrt(covariance), rel=1e-5
    )
    assert (
        fit["pressure_rmse_gpa"]
        <= result["printed_expansion_curve"]["pressure_rmse_gpa"]
    )
    assert all(r["solver_success"] for r in result["fits"].values())


def test_diagnostic_preserves_scope_and_weighting_sensitivity(report):
    check_report(report, json.loads(REPORT.read_text(encoding="utf-8")))
    assert not report["catalog_mutated"]
    all_rows = [row for result in report["isotherms"] for row in result["source_rows"]]
    assert sorted(all_rows) == list(range(1, 38))
    warm = report["isotherms"][1]
    assert warm["fits"]["effective_sigma"]["V0_a3"] == pytest.approx(
        reference_volume(473), abs=0.002
    )
    hot = report["isotherms"][-1]
    assert all(fit["V0_a3"] > 42.67 for fit in hot["fits"].values())
    assert (
        hot["printed_rounding_corner_sensitivity"]["unweighted_V0_range_a3"][0] > 42.68
    )
    assert hot["fits"][MODES[0]]["V0_a3"] == pytest.approx(42.714095, abs=1e-6)
    ambient = report["ambient_linear_alpha_diagnostics"][
        "unweighted_volume_free_300k_anchor"
    ]
    assert ambient["source_rows"] == list(range(1, 12))
    assert ambient["V0_at_isotherms_a3"][-1] == pytest.approx(42.68033, abs=1e-5)
    # The extrapolation is weakly constrained; a close match does not identify
    # the original authors' actual expansion law.
    assert ambient["conditional_V0_standard_errors_at_isotherms_a3"][-1] > 0.15
    assert ambient["parameter_standard_errors"]["dalpha_dT_per_k2"] > 0.9 * abs(
        ambient["dalpha_dT_per_k2"]
    )
