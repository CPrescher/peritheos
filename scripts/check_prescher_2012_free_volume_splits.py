"""Repository-only free-V0 split fits, with the recovered ambient point.

Run ``python -m scripts.check_prescher_2012_free_volume_splits``.
These are independent diagnostic fits, not replacements for published records.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
import scipy

from scripts.recover_prescher_2012_fe3c import (
    ROOT,
    TABLE,
    bm3,
    fit,
    polynomial_pressure_fit,
)
from scripts.recover_prescher_2012_fe3c_splits import TARGETS, compare_fit

AMBIENT = (
    ROOT / "peritheos/data/research/cementite/prescher-2012-ambient-candidates.json"
)
OWN = ROOT / "peritheos/data/research/cementite/prescher-2012-own-shelx-cells.json"
REPORT = ROOT / "docs/data/prescher-2012-free-volume-split-check.json"


def replay():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    p, v = table["pressure_gpa"], table["volume_a3"]
    working = np.array([r["volume_a3"] for r in json.loads(OWN.read_text())["rows"]])
    definitions = [
        ("low_below_22_with_ambient", "low", p < 22, v, "primary"),
        ("high_22_up_final_table", "high", p >= 22, v, "primary"),
        ("low_through_22point6_with_ambient", "low", p < 24, v, "selection_control"),
        ("high_22_up_working_cells", "high", p >= 22, working, "selection_control"),
    ]
    cases = []
    for name, side, mask, volumes, role in definitions:
        ps, vs = p[mask], volumes[mask]
        rows = (np.flatnonzero(mask) + 1).tolist()
        if side == "low":
            ps, vs = np.r_[0.0, ps], np.r_[155.08, vs]
            rows = ["archived_ambient_candidate"] + rows
        result = fit(ps, vs, initial=[155.08, 160, 6])
        values = np.array(list(result["parameters"].values()))
        np.testing.assert_allclose(
            values, polynomial_pressure_fit(ps, vs), rtol=2e-5, atol=2e-5
        )
        # Independent complex-step derivatives check the covariance returned by
        # the nonlinear fit's finite-difference Jacobian.
        columns = []
        for index in range(3):
            perturbed = values.astype(complex)
            perturbed[index] += 1e-20j
            columns.append(bm3(vs, *perturbed).imag / 1e-20)
        jac = np.column_stack(columns)
        inverse = np.linalg.inv(jac.T @ jac)
        scaled = (
            inverse
            * result["sum_squared_objective_residuals"]
            / result["degrees_of_freedom"]
        )
        np.testing.assert_allclose(
            result["residual_scaled_covariance"], scaled, rtol=2e-4, atol=1e-7
        )
        sigma = {"low": 0.2, "high": 0.1}[side]
        cases.append(
            {
                "selection": name,
                "side": side,
                "role": role,
                "source_rows": rows,
                "pressure_gpa": ps.tolist(),
                "volume_a3": vs.tolist(),
                **result,
                "published_K0_K0prime_errorK0_errorK0prime": TARGETS[side],
                "residual_scaled_comparisons": compare_fit(
                    result, TARGETS[side], False
                ),
                "secondary_assumed_absolute_sigma_gpa": sigma,
                "secondary_absolute_errors_V0_K0_K0prime": (
                    sigma * np.sqrt(np.diag(inverse))
                ).tolist(),
                "verification": "Central parameters agree with independent pressure-polynomial solution; covariance agrees with complex-step derivatives",
            }
        )
    return {
        "format": "peritheos.prescher-2012-free-volume-split-check",
        "date": "2026-10-30",
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (TABLE, AMBIENT, OWN)
        },
        "ambient_input": {
            "pressure_gpa": 0.0,
            "volume_a3": 155.08,
            "stored_volume_error_a3": 0.19,
            "source": "May 4 CrysAlis log, 148-observation cell fit; see ambient-candidates evidence",
            "qualification": "P=0 assignment of the ambient candidate is the requested diagnostic assumption; the log does not explicitly attach pressure. The printed volume error is retained but not used as a weight or propagated in this pressure-residual fit.",
        },
        "protocol": {
            "equation": "BM3, conventional-cell A^3 and GPa",
            "free_parameters": ["V0", "K0", "K0_prime"],
            "objective": "Equal-weight pressure residuals; all input volumes treated as exact coordinates",
            "primary_errors": "Linearized conditional one-sigma errors: inverse(J^T J) times SSE/(N-3). Assumes independent homoscedastic pressure residuals; these are estimated from scatter, not propagated experimental coordinate uncertainties.",
            "secondary_errors": "Prior assumed absolute uniform pressure sigmas 0.2/0.1 GPa retained separately for comparison; not recovered author settings. Central values are unchanged by uniform weighting.",
            "source_error_confidence": "Unknown; numeric printed-error comparisons do not establish matching confidence conventions.",
            "primary_selection": "Final Table S1 split at 22 GPa; ambient added only to low. Actual low table pressures 1.8-20.7 GPa; high 22.6-47.4 GPa.",
            "controls": "Low includes 22.6 GPa; high uses working own/SHELX volumes, including non-final p16 volume. Kept separate from primary fits.",
        },
        "cases": cases,
        "author_workflow_reproduced": False,
        "catalog_action": "No source coefficients, adopted catalog V0, previous comparisons or Studio records changed.",
    }


if __name__ == "__main__":
    report = replay()
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(REPORT)
    for case in report["cases"]:
        print(
            case["selection"],
            case["parameters"],
            "residual-scaled errors:",
            case["residual_scaled_errors"],
        )
