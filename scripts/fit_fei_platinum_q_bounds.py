"""Constrain q in the existing Fei Pt diagnostic fits; retain source EOS."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts.reproduce_fei_2004_pressure_scales import (
    PARAMETERS as PARAMETERS_2004,
)
from scripts.reproduce_fei_2004_pressure_scales import (
    load_rows,
)
from scripts.reproduce_fei_2004_pressure_scales import (
    source_pressure as pressure_2004,
)
from scripts.reproduce_fei_2007_platinum import (
    PARAMETERS as PARAMETERS_2007,
)
from scripts.reproduce_fei_2007_platinum import (
    observations,
)
from scripts.reproduce_fei_2007_platinum import (
    source_pressure as pressure_2007,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fei-platinum-q-bounded-diagnostics.json"


def fit_case(
    year, volume, temperature, target, *, q_min=None, q_fixed=None, start_q=0.5
):
    if year == 2007:
        base = np.array(PARAMETERS_2007["platinum"])
        names = ("K0_prime", "gamma0", "q")
        indices = [2, 3, 4]
    else:
        p = PARAMETERS_2004["platinum"]
        base = np.array([p[1], p[2], p[4], p[5]])
        names = ("K0", "K0_prime", "gamma0", "q")
        indices = [0, 1, 2, 3]

    def evaluate(c):
        return (
            pressure_2007(volume, temperature, coefficients=c)
            if year == 2007
            else pressure_2004(volume, temperature, "platinum", c)
        )

    published = evaluate(base)
    q_index = indices[-1]
    if q_fixed is not None:
        base[q_index] = q_fixed
        indices.pop()
    scales = np.array([100.0 if year == 2004 and i == 0 else 1.0 for i in indices])

    def unpack(x):
        c = base.copy()
        c[indices] = x * scales
        return c

    lower = np.full(len(indices), -np.inf)
    upper = np.full(len(indices), np.inf)
    initial = base[indices] / scales
    if q_fixed is None:
        initial[-1] = start_q
        if q_min is not None:
            lower[-1] = q_min
    result = least_squares(
        lambda x: evaluate(unpack(x)) - target,
        initial,
        bounds=(lower, upper),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=2000,
    )
    fitted = unpack(result.x)
    coefficients = fitted[[2, 3, 4]] if year == 2007 else fitted
    return {
        "parameters": dict(zip(names, coefficients.tolist())),
        "observations": len(volume),
        "q_lower_bound": q_min,
        "q_fixed": q_fixed,
        "q_start": start_q if q_fixed is None else None,
        "q_at_lower_bound": bool(
            q_min is not None and q_fixed is None and fitted[q_index] - q_min < 1e-7
        ),
        "solver_success": bool(result.success),
        "active_mask": result.active_mask.tolist(),
        "rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
        "published_rmse_gpa": float(np.sqrt(np.mean((published - target) ** 2))),
        "standard_errors": None,
        "uncertainty_note": "No symmetric regression errors reported: the q constraint changes inference at the boundary. Source covariance and calibration uncertainty remain unrecovered.",
    }


def main():
    volume2007, temperature2007, target2007, _, _ = observations()
    rows = load_rows("platinum")
    volume2004, temperature2004, target2004 = [
        np.array([float(r[k]) for r in rows])
        for k in ("volume_a3", "temperature_k", "pressure_gpa")
    ]
    results = {}
    for year, volume, temperature, target in (
        (2007, volume2007, temperature2007, target2007),
        (2004, volume2004, temperature2004, target2004),
    ):
        cases = {
            "unconstrained": fit_case(year, volume, temperature, target),
            "q_nonnegative": fit_case(year, volume, temperature, target, q_min=0.0),
            "q_positive_epsilon": fit_case(
                year, volume, temperature, target, q_min=1e-6
            ),
            "q_fixed_at_published": fit_case(
                year, volume, temperature, target, q_fixed=0.5
            ),
        }
        # Verify the boundary solution is stable to materially different starts.
        starts = [
            fit_case(year, volume, temperature, target, q_min=0.0, start_q=q)
            for q in (0.05, 2.0)
        ]
        assert all(case["solver_success"] for case in cases.values())
        assert all(case["solver_success"] for case in starts)
        assert all(
            abs(case["rmse_gpa"] - cases["q_nonnegative"]["rmse_gpa"]) < 1e-8
            for case in starts
        )
        assert cases["q_nonnegative"]["q_at_lower_bound"]
        assert cases["q_positive_epsilon"]["q_at_lower_bound"]
        cases["nonnegative_start_sensitivity"] = starts
        results[str(year)] = cases
    source_names = [
        "platinum-dewaele-2004-table1-compression.csv",
        "platinum-fei-2004-table2.csv",
    ]
    report = {
        "scope": "Separate bounded diagnostics. Published material coefficients and original observations are unchanged.",
        "objective": "Equal-weight least squares of pressure residuals; same datasets and pressure calibration as the preceding unconstrained diagnostics.",
        "constraints_2007": "78 observations; V0=60.38 A^3, K0=277 GPa and theta0=230 K fixed; K0-prime, gamma0 and q free except the stated q constraint.",
        "constraints_2004": "42 original Au-Pt rows and original Au-derived pressure targets; V0=60.38 A^3 and theta0=230 K fixed; K0, K0-prime, gamma0 and q free except the stated q constraint.",
        "interpretation": "Both nonnegative fits choose the lower q boundary. A strict-positive epsilon bound likewise reaches its boundary. Bounds remove increasing gamma under compression but do not recover a preferred positive q or establish the author's fit procedure.",
        "input_sha256": {
            name: hashlib.sha256(
                (ROOT / "peritheos/data/datasets" / name).read_bytes()
            ).hexdigest()
            for name in source_names
        },
        "results": results,
    }
    OUTPUT.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                year: {
                    name: {
                        key: value[key]
                        for key in ("parameters", "rmse_gpa", "q_at_lower_bound")
                    }
                    for name, value in cases.items()
                    if name != "nonnegative_start_sensitivity"
                }
                for year, cases in results.items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
