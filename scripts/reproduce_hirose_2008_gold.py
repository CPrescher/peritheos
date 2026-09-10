"""Audit Hirose (2008) Table 1 with independent BM3 equations and staged fits."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/gold-hirose-2008-table1.csv"


def source_pressure(volume, temperature, *, fit=2, coefficients=None):
    """Evaluate the printed equations independently of the Peritheos evaluator."""
    t = np.asarray(temperature)
    if fit == 0:
        v0, k0, kp = 67.85, 167.0, 5.58
        if coefficients is not None:
            kp = coefficients[0]
    elif fit == 1:
        dk, a0, a1 = (-0.028, 3.179e-5, 1.477e-8)
        if coefficients is not None:
            dk, a0, a1 = coefficients
        v0 = 67.85 * np.exp(a0 * (t - 300) + a1 * (t**2 - 300**2) / 2)
        k0, kp = 167 + dk * (t - 300), 5.58
    else:
        b1, b2, b3, a = (1.03e-6, 3.95e-10, 1.61e-13, 3.61e-4)
        if coefficients is not None:
            b1, b2, b3, a = coefficients
        v0 = 67.85 * np.exp(3.824e-5 * (t - 300) + 1.499e-8 * (t**2 - 300**2) / 2)
        k0 = 1 / (
            1 / 167 + b1 * (t - 300) + b2 * (t**2 - 300**2) + b3 * (t**3 - 300**3)
        )
        kp = 5.58 + a * (t - 300) * np.log(t / 300)
    x = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def reproduce():
    """Return output reproduction and explicitly diagnostic unweighted refits."""
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    v = np.array([float(r["volume_a3"]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    printed = np.array([float(r["gold_fit2_pressure_gpa"]) for r in rows])
    result = {}
    for fit, suffix in ((0, "300k"), (1, "fit1"), (2, "fit2")):
        mask = t == 300 if fit == 0 else t > 300
        predicted = source_pressure(v[mask], t[mask], fit=fit)
        record = get_eos_record(f"gold_hirose_2008_bm3_{suffix}")
        native = record.pressure(v[mask], temperature=t[mask])
        entry = {
            "rows": int(mask.sum()),
            "native_equation_max_difference_gpa": float(
                np.max(np.abs(native - predicted))
            ),
            "published_rmse_gpa": float(np.sqrt(np.mean((predicted - p[mask]) ** 2))),
            "published_mae_gpa": float(np.mean(np.abs(predicted - p[mask]))),
        }
        if fit < 2:
            initial = np.array([5.58] if fit == 0 else [-0.028, 3.179e-5, 1.477e-8])
            scales = np.array([1.0] if fit == 0 else [0.01, 1e-5, 1e-8])
            optimized = least_squares(
                lambda values: (
                    source_pressure(
                        v[mask], t[mask], fit=fit, coefficients=values * scales
                    )
                    - p[mask]
                ),
                initial / scales,
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
            )
            entry["refit_parameters"] = (optimized.x * scales).tolist()
            entry["refit_rmse_gpa"] = float(np.sqrt(np.mean(optimized.fun**2)))
            entry["solver_success"] = bool(optimized.success)
        else:
            entry["table1_output_max_difference_gpa"] = float(
                np.max(np.abs(source_pressure(v, t) - printed))
            )
            entry["table1_2070k_pressure_gpa"] = float(source_pressure(3.7152**3, 2070))
            entry["full_fit_reproduced"] = False
        result[suffix] = entry
    return result


def ledger_outcome(record):
    """Expose the source audit without treating a partial fit as complete."""
    suffix = record["identifier"].removeprefix("gold_hirose_2008_bm3_")
    metrics = reproduce()[suffix]
    common = {
        "dataset_identifiers": ["gold_hirose_2008_table1"],
        "observations": metrics["rows"],
        "published_rmse_gpa": metrics["published_rmse_gpa"],
    }
    if suffix == "fit2":
        return {
            **common,
            "status": "not_refittable",
            "reason": "All 21 Table 1 output pressures are reproduced, including 12 high-temperature observations. The full 38-row fit also used 26 Fei et al. (2004) observations not bundled by this initial Hirose audit. Their recovery and complete-fit reproduction are handled in the separate Fei follow-up.",
            "reproduction": metrics,
        }
    names = ["K0_prime"] if suffix == "300k" else ["dK_dT", "alpha0", "alpha1"]
    published = [5.58] if suffix == "300k" else [-0.028, 3.179e-5, 1.477e-8]
    errors = [0.02] if suffix == "300k" else [0.003, 1.39e-6, 3.10e-9]
    values = metrics["refit_parameters"]
    within = all(abs(a - b) <= e for a, b, e in zip(published, values, errors))
    return {
        **common,
        "status": "similar" if within else "parity_not_achieved",
        "fit_kind": "diagnostic_unweighted_pressure",
        "free_parameters": names,
        "parameters": [
            {
                "parameter": name,
                "published": p,
                "refit": v,
                "published_error": e,
                "refit_error": None,
                "relative_difference": abs(v - p) / abs(p),
                "within_combined_2sigma": None,
                "similar": abs(v - p) <= e,
            }
            for name, p, v, e in zip(names, published, values, errors)
        ],
        "rmse_gpa": metrics["refit_rmse_gpa"],
        "solver_success": metrics["solver_success"],
        "reason": "Published row selection and fixed reference coefficients retained. Unweighted pressure residuals are diagnostic because the source does not specify its objective, weights, or covariance. Published coefficients remain unchanged.",
        "reproduction": metrics,
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2))
