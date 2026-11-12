"""Check preserved inputs and conditional fits, not experimental validation."""

import hashlib
import json

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA, observations, properties
from scripts.refit_pyrope_angel_2022_expanded import (
    EXPANDED_NAMES,
    expanded_rows,
    recalibrate,
    token,
)
from scripts.refit_pyrope_angel_2022_subset import (
    GAS_CONSTANT,
    LOWER,
    SCALING,
    UPPER,
    arrays,
    pressure_at_observations,
    solve,
)


def test_source_selection_tokens_and_missing_channels_are_preserved():
    assert token("16.74(50)") == (16.74, 0.5)
    assert token("1384(1.0)") == (1384, 1)
    np.testing.assert_allclose(token("1500.0(3)"), (1500, 0.3), rtol=1e-15)
    assert token("0.0") == (0, None)
    raw = observations("zhang-1998-table-1.csv")
    assert len(raw) == 39
    assert sum(r["source_excluded"] == "true" for r in raw) == 5
    rows = expanded_rows(True)
    assert len(rows) == 191
    counts = {
        name: sum(r["dataset"] == name for r in rows)
        for name in {r["dataset"] for r in rows}
    }
    assert counts["zhang_helium"] == 19
    assert counts["zhang_neon"] == 14
    assert counts["du_pv"] == 9
    assert counts["milani_pv"] == 11
    assert counts["chantel_pv"] == counts["chantel_ks"] == 10
    assert counts["sinogeikin_2000_ks"] == 17
    assert counts["sinogeikin_2002_ks"] == 9  # Do not synthesize three hidden pairs.
    assert next(r for r in rows if r["dataset"] == "chantel_ks")["P"] == 2.3
    assert [r for r in rows if r["dataset"] == "hartwig_vt"][10]["Y"] == 114.121
    assert all(r["sP"] is None for r in rows if r["dataset"].startswith("chantel"))


def test_quartz_recalibration_against_official_native_calculator():
    receipt = json.loads(
        (DATA / "expanded-eosfit/quartz-calibration/checks.json").read_text()
    )
    assert receipt["all_within_tolerance"]
    assert len(receipt["checks"]) == 22
    rows = [r for r in expanded_rows() if r["dataset"] == "milani_pv"]
    for row, checked in zip(rows, receipt["checks"][11:]):
        assert abs(row["P"] - checked["native_pressure_gpa"]) < 0.0006
        assert abs(row["P"] - checked["independent_pressure_gpa"]) < 1e-12
    p, error, v = recalibrate(8.456, 0.012)
    assert abs(p - 8.40913) < 0.0001
    assert error > 0 and 97 < v < 98


def test_modulus_scale_direction_against_independent_physical_state():
    p, _ = source("pyrope")
    c = np.array([p[i] for i in (2, 3, 11, 18)] + [1, 1, 0.97828])
    for pressure, temperature in ((0, 298), (14, 298), (12, 1000)):
        expected = properties(p, pressure, temperature, GAS_CONSTANT)
        data = arrays(
            [
                dict(
                    Y=expected["ks_gpa"] * c[6],
                    P=pressure,
                    T=temperature,
                    sY=2,
                    sP=None,
                    sT=None,
                    kind="KS",
                    scale_index=6,
                )
            ]
        )
        np.testing.assert_allclose(
            pressure_at_observations(c, data), pressure, atol=1e-9
        )


def test_expanded_fit_replay_and_nonparity_receipts():
    report = json.loads((DATA / "expanded-refits.json").read_text())
    assert report["missing_selected_channels"] == 3
    assert report["author_workflow_reproduced"] is False
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    assert (
        hashlib.sha256((DATA / "expanded-refit-plan.json").read_bytes()).hexdigest()
        == report["plan_sha256"]
    )
    for include_t in (False, True):
        rows = expanded_rows(include_t)
        saved = report["cases"][str(len(rows))]["refit"]
        initial = [report["published_parameters"][key] for key in EXPANDED_NAMES]
        fit = solve(
            rows,
            initial,
            list(range(7)),
            names=EXPANDED_NAMES,
            scaling=np.r_[SCALING, 1],
            lower=np.r_[LOWER, 0.9],
            upper=np.r_[UPPER, 1.1],
        )
        assert fit["numerically_converged"] and not fit["at_parameter_bound"]
        assert fit["scaled_jacobian_rank"] == 7
        np.testing.assert_allclose(
            list(fit["parameters"].values()),
            list(saved["parameters"].values()),
            rtol=1e-9,
        )
        native = json.loads(
            (DATA / "expanded-eosfit" / str(len(rows)) / "checks.json").read_text()
        )
        assert native["numerical_convergence_reported"]
        # Keep the failed, predeclared strict solver comparison visible.
        assert native["within_declared_tolerances"] is False
        assert saved["parameters"]["theta0_k"] < 550
        assert report["published_parameters"]["theta0_k"] > 770


def test_published_plot_crosschecks_preserve_unresolved_normalization():
    checked = json.loads((DATA / "expanded-figure-crosschecks.json").read_text())
    du = checked["points"]["du_pv"]
    assert max(abs(r["figure_minus_prepared_P_gpa"]) for r in du) < 0.02
    assert max(abs(r["figure_minus_prepared_V_cm3_mol"]) for r in du) < 0.004
    chantel = checked["points"]["chantel_pv"]
    assert max(abs(r["figure_minus_prepared_P_gpa"]) for r in chantel) < 0.02
    assert all(-0.06 < r["figure_minus_prepared_V_cm3_mol"] < -0.03 for r in chantel)
