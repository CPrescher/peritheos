"""Repository-only BM3 errors-in-both-coordinates diagnostics.

Run ``python -m scripts.check_prescher_2012_coordinate_errors``.
Published observations/coefficients and earlier diagnostic reports are preserved.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import scipy
from scipy import odr
from scipy.optimize import least_squares

from scripts.recover_prescher_2012_fe3c import (
    ROOT,
    SOURCE,
    TABLE,
    bm3,
    volume_at_pressure,
)
from scripts.recover_prescher_2012_fe3c_splits import matches_printed

AMBIENT = (
    ROOT / "peritheos/data/research/cementite/prescher-2012-ambient-candidates.json"
)
REPORT = ROOT / "docs/data/prescher-2012-coordinate-error-fits.json"


def coordinate_fit(p, v, sp, sv, ambient=None):
    """Jointly fit EOS coefficients and latent volumes; ambient P is fixed at 0.

    ``ambient=(volume, volume_error)`` adds a volume observation at P=0,
    rather than fixing V0. Its pressure error is unknown, not assigned zero.
    The exact-pressure treatment is an explicit diagnostic constraint.
    """
    p, v, sp, sv = (np.asarray(a, dtype=float) for a in (p, v, sp, sv))

    def residual(b):
        r = np.r_[(bm3(b[3:], *b[:3]) - p) / sp, (b[3:] - v) / sv]
        return np.r_[r, (b[0] - ambient[0]) / ambient[1]] if ambient else r

    opt = least_squares(
        residual,
        np.r_[[155.08, 160, 6], v],
        jac="cs",
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=5000,
    )
    if not opt.success:
        raise RuntimeError(opt.message)
    if np.linalg.matrix_rank(opt.jac) != len(opt.x):
        raise RuntimeError("Rank-deficient covariance")
    # The full inverse includes the fitted volume corrections; taking only
    # the coefficient Jacobian would silently condition on those corrections.
    cov = np.linalg.inv(opt.jac.T @ opt.jac)[:3, :3]
    errors = np.sqrt(np.diag(cov))
    n = len(p) + int(ambient is not None)
    dof = n - 3
    sse = float(opt.fun @ opt.fun)
    return {
        "parameters": dict(zip(("V0", "K0", "K0_prime"), opt.x[:3].tolist())),
        "absolute_errors_V0_K0_K0prime": errors.tolist(),
        "absolute_covariance": cov.tolist(),
        "residual_scaled_errors_V0_K0_K0prime": (errors * np.sqrt(sse / dof)).tolist(),
        "residual_scaled_covariance": (cov * sse / dof).tolist(),
        "correlation": (cov / np.outer(errors, errors)).tolist(),
        "observations": n,
        "degrees_of_freedom": dof,
        "weighted_sum_squares": sse,
        "reduced_chi_square": sse / dof,
        "latent_table_volumes_a3": opt.x[3:].tolist(),
        "normalized_residuals": opt.fun.tolist(),
        "success": bool(opt.success),
        "termination": opt.message,
    }


def odr_check(p, v, sp, sv, ambient):
    """Independent ODRPACK fit using P as input and V as output."""
    if ambient:
        p, v = np.r_[0, p], np.r_[ambient[0], v]
        # ODR fixes ambient input P. The placeholder weight at this fixed
        # coordinate is unused and is not an invented experimental sigma.
        sp, sv = np.r_[1, sp], np.r_[ambient[1], sv]
        fixed = np.r_[0, np.ones(len(p) - 1, dtype=int)]
    else:
        fixed = np.ones(len(p), dtype=int)
    data = odr.RealData(p, v, sx=sp, sy=sv, fix=fixed)
    result = odr.ODR(
        data,
        odr.Model(lambda b, pressure: volume_at_pressure(pressure, *b)),
        beta0=[155.08, 160, 6],
        maxit=300,
        sstol=1e-10,
        partol=1e-10,
    ).run()
    if result.info not in (1, 2, 3):
        raise RuntimeError(result.stopreason)
    return {
        "parameters_V0_K0_K0prime": result.beta.tolist(),
        "absolute_covariance": result.cov_beta.tolist(),
        "weighted_sum_squares": float(result.sum_square),
        "stop_reason": result.stopreason,
    }


def replay(full_without_ambient=False):
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    pressure = table["pressure_gpa"]
    definitions = [
        ("low_below_22_with_ambient", "below_22_gpa", pressure < 22, True, "primary"),
        ("high_22_up", "above_22_gpa", pressure >= 22, False, "primary"),
        ("full_with_ambient", "full_range", np.ones(len(table), bool), True, "primary"),
        (
            "low_through_22point6_with_ambient",
            "below_22_gpa",
            pressure < 24,
            True,
            "selection_control",
        ),
    ]
    if full_without_ambient:
        definitions = [
            (
                "full_without_ambient",
                "full_range",
                np.ones(len(table), bool),
                False,
                "primary",
            )
        ]
    cases = []
    for name, source_key, mask, include_ambient, role in definitions:
        p, v, sp, sv = [
            table[key][mask]
            for key in (
                "pressure_gpa",
                "volume_a3",
                "pressure_error_gpa",
                "volume_error_a3",
            )
        ]
        ambient = (155.08, 0.19) if include_ambient else None
        result = coordinate_fit(p, v, sp, sv, ambient)
        check = odr_check(p, v, sp, sv, ambient)
        np.testing.assert_allclose(
            list(result["parameters"].values()),
            check["parameters_V0_K0_K0prime"],
            rtol=1e-5,
            atol=5e-5,
        )
        np.testing.assert_allclose(
            result["absolute_covariance"],
            check["absolute_covariance"],
            rtol=5e-4,
            atol=1e-6,
        )
        np.testing.assert_allclose(
            result["weighted_sum_squares"],
            check["weighted_sum_squares"],
            rtol=1e-7,
            atol=1e-7,
        )
        published, errors = SOURCE[source_key]
        values = [result["parameters"][key] for key in ("K0", "K0_prime")]
        comparisons = {}
        for convention in ("absolute", "residual_scaled"):
            recovered = values + result[convention + "_errors_V0_K0_K0prime"][-2:]
            comparisons[convention] = {
                "K0_K0prime_errorK0_errorK0prime": recovered,
                "matches_printed_digits": matches_printed(
                    recovered, published + errors
                ),
            }
        cases.append(
            {
                "selection": name,
                "role": role,
                "table_row_indices_1based": (np.flatnonzero(mask) + 1).tolist(),
                "pressure_gpa": p.tolist(),
                "volume_a3": v.tolist(),
                "pressure_error_gpa": sp.tolist(),
                "volume_error_a3": sv.tolist(),
                "ambient_added": include_ambient,
                **result,
                "published_K0_K0prime_errorK0_errorK0prime": published + errors,
                "comparisons": comparisons,
                "independent_ODRPACK_check": check,
            }
        )
    return {
        "format": "peritheos.prescher-2012-coordinate-error-fits",
        "date": "2026-10-30",
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in ((TABLE,) if full_without_ambient else (TABLE, AMBIENT))
        },
        "protocol": {
            "equation": "BM3, GPa and conventional-cell A^3; V0, K0, K0_prime all free",
            "objective": "Sum of squared normalized pressure and volume corrections, using latent table volumes; plus ((V0-155.08)/0.19)^2 for the ambient point when included",
            "ambient_pressure": "Ambient observation omitted; no pressure coordinate is fixed. V0 remains free."
            if full_without_ambient
            else "Fixed at 0 GPa as an explicit diagnostic constraint; source pressure error and explicit pressure metadata are not recovered. V0 remains free.",
            "input_errors": "Table S1 pressure and volume errors, conditionally treated as independent absolute one-sigma widths. Ambient volume error 0.19 A^3 is additionally used only when the ambient point is included. Actual source confidence and coordinate/calibration covariance are unknown; pressure widths derive from before/after ruby differences.",
            "primary_errors": "Square roots of the coefficient block of the full inverse normal matrix, including latent-volume nuisance parameters; no residual-variance scaling",
            "secondary_errors": "Same covariance multiplied by weighted SSE/(N-3), retained separately",
            "full_fit_rows": "All 19 final Table S1 rows, P=1.8-47.4 GPa, with no ambient point"
            if full_without_ambient
            else "All 19 final Table S1 rows plus the ambient candidate; the other two figure-only markers have no recovered numeric coordinate errors and are excluded",
            "verification": "ODRPACK independently fits inverse V(P), fixing only the ambient pressure coordinate when included; parameters, covariance and objective agree with direct P(V) latent-volume least squares",
        },
        "cases": cases,
        "author_workflow_reproduced": False,
        "catalog_action": "No catalog, previous diagnostic or Studio values changed.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--without-ambient",
        action="store_true",
        help="Fit only the full 19-row table and save a separate report",
    )
    args = parser.parse_args()
    report = replay(full_without_ambient=args.without_ambient)
    output = (
        REPORT.with_name("prescher-2012-full-without-ambient.json")
        if args.without_ambient
        else REPORT
    )
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(output)
    for case in report["cases"]:
        print(
            case["selection"],
            case["parameters"],
            "absolute errors",
            case["absolute_errors_V0_K0_K0prime"],
            "residual-scaled errors",
            case["residual_scaled_errors_V0_K0_K0prime"],
        )
