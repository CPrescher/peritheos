"""Independent BM fits to Nisr (2017) official Tables S1, S2 and S4.

The article specifies least squares but not its dependent variable or weights.
Unweighted pressure residuals give the principal diagnostic; volume residuals
and positive-sigma volume residuals quantify that unresolved choice. No result
from this script replaces the published parameterization.
"""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/nisr-2017-silica-tables-s1-s4.csv"
PUBLISHED = {
    "stishovite_nisr_2017_dry_bm3": ("S1", 46.569, 312.0, 4.59),
    "hydrous_stishovite_nisr_2017_bm3": ("S2", 47.191, 257.0, 4.59),
    "hydrous_silica_cacl2_nisr_2017_bm2": ("S4", 47.23, 286.0, 4.0),
}
ERRORS = {"S1": {"K0": 2.0}, "S2": {"K0": 9.0}, "S4": {"V0": 0.36, "K0": 18.0}}


def load_rows():
    with DATA.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def pressure(volume, v0, k0, kp):
    """Eulerian third-order Birch-Murnaghan; kp=4 is exactly BM2."""
    x = (v0 / np.asarray(volume)) ** (2.0 / 3.0)
    return 1.5 * k0 * x**2.5 * (x - 1.0) * (1.0 + 0.75 * (kp - 4.0) * (x - 1.0))


def volumes(pressures, v0, k0, kp):
    return np.array(
        [
            brentq(lambda v: pressure(v, v0, k0, kp) - p, 0.5 * v0, 1.1 * v0)
            for p in pressures
        ]
    )


def fit_branch(record_id, objective="pressure"):
    table, v0, k0, kp = PUBLISHED[record_id]
    selected = [r for r in load_rows() if r["table"] == table]
    p, v, sigma = np.array(
        [
            [float(r[n]) for n in ("pressure_gpa", "volume_a3", "volume_sigma_a3")]
            for r in selected
        ]
    ).T
    mask = sigma > 0 if objective == "volume_sigma" else np.ones(len(p), dtype=bool)
    names = ["V0", "K0"] if table == "S4" else ["K0"]

    def decode(values):
        return (values[0], values[1]) if table == "S4" else (v0, values[0])

    def residual(values):
        vr, kr = decode(values)
        if objective == "pressure":
            return pressure(v, vr, kr, kp) - p
        delta = volumes(p, vr, kr, kp) - v
        return delta[mask] / (sigma[mask] if objective == "volume_sigma" else 1.0)

    # Initial guesses deliberately do not equal the published free coefficients.
    initial = [47.0, 300.0] if table == "S4" else [300.0]
    fit = least_squares(
        residual, initial, x_scale="jac", ftol=1e-12, xtol=1e-12, gtol=1e-12
    )
    vr, kr = decode(fit.x)
    dof = int(mask.sum()) - len(names)
    # Residual-scaled diagnostic covariance, not the source covariance.
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * np.sum(fit.fun**2) / dof
    errors = np.sqrt(np.diag(covariance))
    curve_residual = pressure(v, vr, kr, kp) - p
    published_residual = pressure(v, v0, k0, kp) - p
    return {
        "table": table,
        "observations": int(mask.sum()),
        "source_rows": len(selected),
        "excluded_zero_sigma_pressures_gpa": p[~mask].tolist(),
        "objective": objective,
        "free_parameters": names,
        "parameters": dict(zip(names, map(float, fit.x))),
        "standard_errors": dict(zip(names, map(float, errors))),
        "covariance": covariance.tolist(),
        "degrees_of_freedom": dof,
        "solver_success": bool(fit.success),
        "pressure_range_gpa": [float(p.min()), float(p.max())],
        "rmse_gpa": float(np.sqrt(np.mean(curve_residual**2))),
        "published_rmse_gpa": float(np.sqrt(np.mean(published_residual**2))),
        "published_max_abs_residual_gpa": float(np.max(np.abs(published_residual))),
        "published_volume_rmse_a3": float(
            np.sqrt(np.mean((volumes(p, v0, k0, kp) - v) ** 2))
        ),
        "high_pressure_benchmark": {
            "source_pressure_gpa": float(p[-1]),
            "source_volume_a3": float(v[-1]),
            "calculated_pressure_gpa": float(pressure(v[-1], v0, k0, kp)),
            "calculated_volume_a3": float(volumes(p[-1:], v0, k0, kp)[0]),
        },
    }


@lru_cache(maxsize=1)
def reproduce():
    return {
        "doi": "10.1002/2017JB014055",
        "version": "Corrected 26 October 2017; official supporting information S1",
        "fits": {
            identifier: {
                mode: fit_branch(identifier, mode)
                for mode in ("pressure", "volume", "volume_sigma")
            }
            for identifier in PUBLISHED
        },
        "water_content_audit": {
            "source_formula": "Si0.954O2H0.184",
            "formula_implied_water_wt_percent": 100
            * (0.184 / 2 * 18.01528)
            / (0.954 * 28.0855 + 2 * 15.9994 + 0.184 * 1.00794),
            "volume_calibrated_water_wt_percent": 3.2,
            "reported_uncertainty_wt_percent": 0.5,
            "corrected_table2_volume_expansion_percent": 100 * (47.198 / 46.569 - 1),
        },
    }


def ledger_outcome(record):
    identifier = record["identifier"]
    result = reproduce()["fits"][identifier]["pressure"]
    published = record["eos"]["parameters"]
    comparisons = []
    for name, value in result["parameters"].items():
        error = ERRORS[result["table"]][name]
        fit_error = result["standard_errors"][name]
        comparisons.append(
            {
                "parameter": name,
                "published": published[name],
                "published_error": error,
                "refit": value,
                "refit_error": fit_error,
                "difference": value - published[name],
                "relative_difference": abs(value - published[name]) / published[name],
                "within_combined_2sigma": abs(value - published[name])
                <= 2 * np.hypot(error, fit_error),
                "similar": abs(value - published[name]) <= error,
            }
        )
    # Use Python bools so JSON serialization is independent of NumPy versions.
    for item in comparisons:
        item["within_combined_2sigma"] = bool(item["within_combined_2sigma"])
    return {
        "status": "parity"
        if all(p["similar"] for p in comparisons)
        else "parity_not_achieved",
        "dataset_identifiers": record["fit_datasets"],
        "observations": result["observations"],
        "selection": f"All Table {result['table']} rows, including its ambient row when present; Table S3 coexistence rows excluded.",
        "observed_pressure_range_gpa": result["pressure_range_gpa"],
        "fit_kind": "diagnostic_unweighted_pressure",
        "objective": "unweighted pressure residuals; source dependent variable and weights unspecified",
        "absolute_sigma": False,
        "free_parameters": result["free_parameters"],
        "parameters": comparisons,
        "rmse_gpa": result["rmse_gpa"],
        "published_rmse_gpa": result["published_rmse_gpa"],
        "degrees_of_freedom": result["degrees_of_freedom"],
        "solver_success": result["solver_success"],
        "qualification": "Coefficient parity within published 1-sigma widths, not recovery of an unpublished weighting protocol. Pressure errors, gold volumes and source covariance are absent. Printed zero volume errors are retained, never treated as infinite weights. Volume-objective sensitivity fits are in the paper-specific reproduction.",
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True, allow_nan=False))
