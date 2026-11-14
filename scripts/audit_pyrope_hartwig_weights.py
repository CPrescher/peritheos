"""User-requested Hartwig weight sensitivity; diagnostic, not curated EOS fits.

Repository-only: python -m scripts.audit_pyrope_hartwig_weights.
All multiplier choices and controls are in hartwig-weight-plan.json.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.refit_pyrope_angel_2022_expanded import EXPANDED_NAMES, expanded_rows
from scripts.refit_pyrope_angel_2022_subset import (
    LOWER,
    SCALING,
    UPPER,
    arrays,
    effective_sigma,
    pressure_at_observations,
    solve,
)


def weighted_rows(rows, multiplier):
    """Change only copied Hartwig errors, giving weight multiplier w exactly."""
    if not np.isfinite(multiplier) or multiplier <= 0:
        raise ValueError("Weight multiplier must be finite and positive")
    root = np.sqrt(multiplier)
    return [
        {
            **r,
            "sY": r["sY"] / root,
            "sT": r["sT"] / root,
        }
        if r["dataset"] == "hartwig_vt"
        else dict(r)
        for r in rows
    ]


def fit_rows(rows, initial):
    return solve(
        rows,
        initial,
        list(range(7)),
        names=EXPANDED_NAMES,
        scaling=np.r_[SCALING, 1],
        lower=np.r_[LOWER, 0.9],
        upper=np.r_[UPPER, 1.1],
    )


def baseline_residual_summary(rows, parameters):
    """Comparable diagnostics under the original errors at every multiplier."""
    c = np.array([parameters[k] for k in EXPANDED_NAMES])
    d = arrays(rows)
    residual = pressure_at_observations(c, d) - d["P"]
    z = residual / effective_sigma(c, d)
    summary = {}
    for name in sorted({r["dataset"] for r in rows}):
        selected = np.array([r["dataset"] == name for r in rows])
        summary[name] = {
            "channels": int(sum(selected)),
            "baseline_error_sum_squares": float(z[selected] @ z[selected]),
            "pressure_residual_rms_gpa": float(
                np.sqrt(np.mean(residual[selected] ** 2))
            ),
            "max_abs_standardized_residual": float(max(abs(z[selected]))),
        }
    return summary


def run():
    plan_path = DATA / "hartwig-weight-plan.json"
    plan = json.loads(plan_path.read_text())
    baseline = json.loads((DATA / "expanded-refits.json").read_text())
    p, covariance = source("pyrope")
    initial = np.array([p[i] for i in (2, 3, 11, 18)] + [1.001, 1.0003, p[53]])
    cases = {}
    for include_t in (False, True):
        rows = expanded_rows(include_t)
        previous = initial.copy()
        results = []
        for w in plan["weight_multipliers"]:
            weighted = weighted_rows(rows, w)
            fit = fit_rows(weighted, initial)
            continuation = fit_rows(weighted, previous)
            previous = np.array([fit["parameters"][key] for key in EXPANDED_NAMES])
            summary = baseline_residual_summary(rows, fit["parameters"])
            # Large repeated row maps remain in the original expanded-fit file.
            # The unchanged inputs, w and complete conditional covariance suffice
            # to reproduce each sensitivity fit; keep all solver traces.
            fit.pop("rows")
            continuation.pop("rows")
            results.append(
                {
                    "weight_multiplier": w,
                    "effective_Hartwig_volume_error_cm3_mol": 0.03 / np.sqrt(w),
                    "effective_Hartwig_fractional_T_error": 0.01 / np.sqrt(w),
                    "fit_from_published_start": fit,
                    "fit_from_continuation_start": continuation,
                    "continuation_minus_published_start": {
                        key: continuation["parameters"][key] - fit["parameters"][key]
                        for key in EXPANDED_NAMES
                    },
                    "minus_published_coefficients": {
                        key: fit["parameters"][key] - initial[i]
                        for i, key in enumerate(EXPANDED_NAMES[:4])
                    },
                    "baseline_error_residuals": summary,
                    "baseline_error_total_sum_squares": sum(
                        v["baseline_error_sum_squares"] for v in summary.values()
                    ),
                }
            )
            print(
                len(rows),
                w,
                *(round(fit["parameters"][k], 6) for k in EXPANDED_NAMES[:4]),
                fit["numerically_converged"],
                flush=True,
            )
        cases[str(len(rows))] = results
    result = {
        "scope": "User-requested weight sensitivity, not full author-workflow reproduction or a recommended replacement",
        "author_workflow_reproduced": False,
        "plan_sha256": hashlib.sha256(plan_path.read_bytes()).hexdigest(),
        "baseline_fit_sha256": hashlib.sha256(
            (DATA / "expanded-refits.json").read_bytes()
        ).hexdigest(),
        "source_parameters": dict(zip(EXPANDED_NAMES[:4], initial[:4].tolist())),
        "source_covariance_marginal_widths": dict(
            zip(
                EXPANDED_NAMES[:4],
                np.sqrt(np.diag(covariance)[[1, 2, 10, 17]]).tolist(),
            )
        ),
        "input_sha256": baseline["input_sha256"],
        "cases": cases,
        "interpretation": "Changing w changes the objective; weighted reduced chi squares across w are not a goodness-of-fit ranking under one error model. Baseline-error residual summaries use unchanged coordinate errors for a matched descriptive comparison. Conditional covariance at any w excludes unknown systematic errors and does not recover source confidence.",
    }
    (DATA / "hartwig-weight-sensitivity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


if __name__ == "__main__":
    run()
