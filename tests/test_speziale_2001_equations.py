"""Independent analytic, physical-limit and original-source equation checks."""

import hashlib
import json

import numpy as np
import pytest
from scipy.constants import R

from scripts.reproduce_speziale_2001_mgo import energy, gamma_theta, pressure
from scripts.validate_speziale_2001_equations import (
    CURVES,
    OUTPUT,
    PROVENANCE,
    closed_form_pressure,
    closed_form_theta,
    validate,
)


def test_independent_exponential_integral_solution_and_pressure():
    v = 74.71 * np.array([0.6, 0.64, 0.7, 0.8, 0.9, 1, 1.02])
    assert closed_form_theta(v) == pytest.approx(gamma_theta(v)[1], abs=1e-9)
    for volume in v:
        for t in (300, 1100, 3663):
            assert closed_form_pressure(volume, t) == pytest.approx(
                pressure(volume, t), abs=1e-9
            )


def test_debye_energy_has_mgo_dulong_petit_limit():
    # Two atoms per MgO formula unit imply 6R per mole MgO, independent of theta.
    theta = np.array([773, 1200, 1600])
    assert energy(theta, 1e8) / 1e8 == pytest.approx(np.full(3, 6 * R), rel=1e-5)
    # Eq. 4 with instantaneous q(V) fails the defining gamma/theta identity.
    v, h = 74.71 * 0.8, 1e-5
    _, lo = gamma_theta(v * np.exp(-h), theta_law="eq4_local_q")
    _, hi = gamma_theta(v * np.exp(h), theta_law="eq4_local_q")
    inferred_gamma = -(np.log(hi) - np.log(lo)) / (2 * h)
    assert abs(inferred_gamma - gamma_theta(v)[0]) > 1


def test_archived_equation_report_and_source_curve_provenance():
    from scripts.reproduce_fei_2007_gold import check_saved

    report = validate()
    check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
    source = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    assert source["csv_sha256"] == hashlib.sha256(CURVES.read_bytes()).hexdigest()
    assert source["kind"] == (
        "digitized published calculated isotherms, not experimental observations"
    )
    assert report["source_curve_checks"]["6"]["checkpoints"] == 11
    assert report["source_curve_checks"]["10"]["checkpoints"] == 4
    assert report["source_curve_checks"]["6"]["max_abs_volume_ratio_difference"] < 5e-4
    # Preserve the graphical discrepancy; it must not become exact author parity.
    assert report["source_curve_checks"]["10"]["max_abs_volume_ratio_difference"] > 1e-3
    assert report["exact_historical_variable_q_implementation_verified"] is False
    assert report["numerical"]["closed_form_max_difference_gpa"] < 1e-9
