"""Independent convergence evidence for the author's confirmed theta integral."""

import json

import numpy as np
import pytest

from scripts.audit_speziale_2001_trapezoidal import (
    OUTPUT,
    audit,
    trapezoidal_theta,
    trial_pressure,
)
from scripts.reproduce_fei_2007_gold import check_saved
from scripts.reproduce_speziale_2001_mgo import bm3
from scripts.validate_speziale_2001_equations import closed_form_theta


@pytest.mark.parametrize("coordinate", ["log_volume", "volume"])
def test_theta_quadrature_converges_to_independent_analytic_solution(coordinate):
    volumes = 74.71 * np.array([0.6, 0.64, 0.8, 0.9, 1.0, 1.02])
    reference = closed_form_theta(volumes)
    errors = [
        np.max(np.abs(trapezoidal_theta(volumes, n, coordinate) - reference))
        for n in (16, 32, 64, 256, 1024)
    ]
    assert np.all(np.diff(errors) < 0)
    assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.02)
    assert errors[-1] < 2e-4
    # Every quadrature grid preserves the published reference boundary.
    assert trapezoidal_theta(74.71, 1, coordinate) == 773.0
    for debye in (False, True):
        assert trial_pressure(volumes, 300.0, 4, coordinate, debye) == pytest.approx(
            bm3(volumes), abs=1e-12
        )


def test_archived_sensitivity_and_comparison_scope():
    result = audit()
    check_saved(json.loads(OUTPUT.read_text()), result)
    assert len(result["trials"]) == 40
    assert result["reference"]["fei2004_rmse_gpa"] > 0.25
    assert result["reference"]["hirose2006_max_abs_residual_gpa"] > 1.9
    for row in result["trials"]:
        if row["debye_quadrature"] == "converged_gauss_legendre":
            # Changing only the confirmed theta quadrature does not recover
            # the later printed pressures or the digitized Figure 10 curve.
            assert row["fei2004_rmse_gpa"] > 0.25
            assert row["hirose2006_max_abs_residual_gpa"] > 1.9
            assert row["figure10_max_abs_residual_gpa"] > 0.29
            if row["intervals"] >= 16:
                assert row["figure10_max_pressure_change_gpa"] < 5e-6
            if row["intervals"] >= 64:
                assert row["later28_max_pressure_change_gpa"] < 1.4e-5
