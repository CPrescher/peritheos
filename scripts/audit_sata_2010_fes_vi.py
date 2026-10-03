#!/usr/bin/env python3
"""Replay the cold FeS-VI mBM3 EOS from Sata (2010), Tables 1 and 6.

The source table re-reports Ohfuji (2007) and Sata (2008) observations. It
does not establish access to their complete original fit inputs. Refits are
conditional audit diagnostics, not replacement database coefficients.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data/fes-sata-2010-table1.csv"
OUTPUT = ROOT / "docs/data/sata-2010-fes-vi-audit.json"
VR = 98.96  # 12.37 A3/atom * 8 atoms/cell; selected reference, no quoted error.
PARAMETERS = np.array([0.0, 148.0, 4.53])  # Pr, Kr, K'r at fixed Vr.
SOURCE_ERRORS = np.array([4.2, 16.0, 0.34])


def parse_token(token):
    """Decode last-digit errors; absence is None, not zero or inferred sigma."""
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\((\d+)\))?", token)
    if match is None:
        raise ValueError(f"Invalid source number: {token!r}")
    value, error = match.groups()
    decimals = len(value.partition(".")[2])
    return float(value), None if error is None else int(error) * 10.0**-decimals


def observations():
    """Retain both phases and raw precision; fit only the phase-labelled VI."""
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for index, row in enumerate(rows, 1):
        row["table_fes_row"] = index
        row["z"] = int(row["z"])
        row["included"] = row["phase"] == "VI"
        row["temperature_k"] = 300
        for raw, name in (
            ("volume_raw", "volume_angstrom3"),
            ("pressure_raw", "pressure_gpa"),
            ("vs_vfe_raw", "vs_vfe"),
        ):
            row[name], row[f"error_{name}"] = parse_token(row[raw])
        row["molar_volume_cm3_per_mol_fes"] = (
            row["volume_angstrom3"] * Avogadro / (row["z"] * 1e24)
        )
    return rows


def pressure(volume, parameters=PARAMETERS, vr=VR):
    """Sata Eq. 1; finite-reference-pressure mBM3, pressure in GPa."""
    pr, kr, kp = parameters
    f = ((vr / np.asarray(volume)) ** (2 / 3) - 1) / 2
    h = pr + (3 * kr - 5 * pr) * f
    h += (4.5 * kr * (kp - 4) + 17.5 * pr) * f**2
    return (1 + 2 * f) ** 2.5 * h


def elastic_properties(volume, parameters=PARAMETERS, vr=VR):
    """Analytic K=-V*dP/dV and K'=dK/dP for equation/reference checks."""
    pr, kr, kp = parameters
    f = ((vr / np.asarray(volume)) ** (2 / 3) - 1) / 2
    s = 1 + 2 * f
    b, c = 3 * kr - 5 * pr, 4.5 * kr * (kp - 4) + 17.5 * pr
    h, dh = pr + b * f + c * f**2, b + 2 * c * f
    j = 5 * h + s * dh
    bulk = s**2.5 * j / 3
    derivative = (5 * j + s * (7 * dh + 2 * c * s)) / (3 * j)
    return bulk, derivative


def linear_design(volume, vr=VR):
    """Independent linear form in (Pr, Kr, Kr*K'r) at a fixed reference."""
    f = ((vr / volume) ** (2 / 3) - 1) / 2
    return (1 + 2 * f)[:, None] ** 2.5 * np.column_stack(
        [1 - 5 * f + 17.5 * f**2, 3 * f - 18 * f**2, 4.5 * f**2]
    )


def summary(residual):
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "mean_gpa": float(np.mean(residual)),
        "min_gpa": float(np.min(residual)),
        "max_gpa": float(np.max(residual)),
    }


def fit(volume, observed, dp, dv, objective):
    """Fit three coefficients; expose conditional covariance and objective."""
    design = linear_design(volume)
    scale = dp if objective == "quoted_pressure_errors" else np.ones_like(dp)
    weighted = design / scale[:, None]
    beta, _, rank, singular = np.linalg.lstsq(weighted, observed / scale, rcond=None)
    parameters = np.array([beta[0], beta[1], beta[2] / beta[1]])
    transform = np.array(
        [[1, 0, 0], [0, 1, 0], [0, -beta[2] / beta[1] ** 2, 1 / beta[1]]]
    )
    _, singular, vt = np.linalg.svd(weighted, full_matrices=False)
    if rank != 3:
        raise ValueError("Rank-deficient cold fit")
    covariance = transform @ ((vt.T / singular**2) @ vt) @ transform.T
    converged, at_bound, iterations = True, False, None
    if objective == "effective_variance":

        def residuals(values):
            k, _ = elastic_properties(volume, values)
            effective = np.sqrt(dp**2 + (k * dv / volume) ** 2)
            return (pressure(volume, values) - observed) / effective

        result = least_squares(
            residuals,
            parameters,
            bounds=([-100, 1, 0], [100, 1000, 20]),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=10000,
        )
        parameters = result.x
        k, _ = elastic_properties(volume, parameters)
        scale = np.sqrt(dp**2 + (k * dv / volume) ** 2)
        _, singular, vt = np.linalg.svd(result.jac, full_matrices=False)
        rank = int(np.sum(singular > singular[0] * 1e-12))
        if rank != 3:
            raise ValueError("Rank-deficient effective-variance fit")
        covariance = (vt.T / singular**2) @ vt
        converged = bool(result.success)
        at_bound = bool(np.any(result.active_mask))
        iterations = int(result.nfev)
    residual = pressure(volume, parameters) - observed
    dof = len(volume) - 3
    objective_sum = float(np.sum((residual / scale) ** 2))
    scaled_covariance = covariance * objective_sum / dof
    return {
        "objective": objective,
        "parameter_order": ["Pr_gpa", "Kr_gpa", "Kr_prime"],
        "parameters": parameters.tolist(),
        "difference_over_quoted_parameter_error": (
            (parameters - PARAMETERS) / SOURCE_ERRORS
        ).tolist(),
        "fixed_reference_volume_angstrom3": VR,
        "rank": int(rank),
        "singular_values": singular.tolist(),
        "degrees_of_freedom": dof,
        "objective_sum": objective_sum,
        "objective_sum_per_dof": objective_sum / dof,
        "covariance_unscaled": covariance.tolist(),
        "covariance_residual_scaled": scaled_covariance.tolist(),
        "standard_errors_residual_scaled": np.sqrt(np.diag(scaled_covariance)).tolist(),
        "covariance_qualification": (
            "Local conditional Jacobian covariance. Residual-scaled matrix uses "
            "objective_sum/(N-3); unscaled matrix uses the declared scales. "
            "Quoted source errors have unspecified confidence. Neither matrix "
            "includes calibration covariance, fixed-Vr uncertainty, rounding "
            "or cross-row systematics; neither is an author covariance."
        ),
        "converged": converged,
        "at_bound": at_bound,
        "function_evaluations": iterations,
        "algorithm": (
            "nonlinear least_squares, with parameter-dependent effective scale"
            if objective == "effective_variance"
            else "linear SVD in (Pr, Kr, Kr*Kprime), transformed to physical parameters"
        ),
        **summary(residual),
        "residuals_gpa": residual.tolist(),
    }


def reproduce():
    rows = observations()
    selected = [row for row in rows if row["included"]]
    v, p, dv, dp = (
        np.array([row[key] for row in selected])
        for key in (
            "volume_angstrom3",
            "pressure_gpa",
            "error_volume_angstrom3",
            "error_pressure_gpa",
        )
    )
    residual = pressure(v) - p
    k, _ = elastic_properties(v)
    fits = {
        objective: fit(v, p, dp, dv, objective)
        for objective in (
            "unweighted_pressure",
            "quoted_pressure_errors",
            "effective_variance",
        )
    }
    # Transform the source curve to a finite-P reference; same mBM3 curve.
    repaired_2007_vr = 12.615 * 4e24 / Avogadro
    kr, kp = elastic_properties(repaired_2007_vr)
    transformed = [float(pressure(repaired_2007_vr)), float(kr), float(kp)]
    finite_replay = pressure(v, [36.0, 306.0, 3.81], vr=repaired_2007_vr) - p
    return {
        "audit_id": "sata_2010_fes_vi_cold_mbm3",
        "reference_doi": "10.1029/2009JB006975",
        "dataset_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "disposition": "approximate_coefficient_reproduction_protocol_unresolved",
        "qualification": (
            "All 13 FeS-VI rows re-reported in Sata (2010) Table 1 are included; "
            "six FeS-VII rows are retained but excluded by phase. The original "
            "2007/2008 full fit inputs and author weighting are not available. "
            "Approximate coefficient parity is not exact protocol reproduction "
            "or independent physical/calibration validation. No refit is promoted."
        ),
        "source_rows": rows,
        "included_runs": [row["run"] for row in selected],
        "source_temperature_k": 300,
        "measurement_range": {
            "pressure_gpa": [float(p.min()), float(p.max())],
            "volume_angstrom3": [float(v.min()), float(v.max())],
        },
        "published_reference": {
            "volume_angstrom3_per_atom": 12.37,
            "volume_angstrom3_per_cell": VR,
            "molar_volume_cm3_per_mol_fes": VR * Avogadro / 4e24,
            "fixed_volume_error": None,
            "parameter_order": ["Pr_gpa", "Kr_gpa", "Kr_prime"],
            "parameters": PARAMETERS.tolist(),
            "quoted_errors": SOURCE_ERRORS.tolist(),
            "parameter_covariance": None,
            "weighting": "not recovered",
        },
        "published_replay": {
            **summary(residual),
            "residuals_gpa": residual.tolist(),
            "sum_squared_residual_over_pressure_error": float(
                np.sum((residual / dp) ** 2)
            ),
            "sum_squared_residual_over_effective_error": float(
                np.sum(residual**2 / (dp**2 + (k * dv / v) ** 2))
            ),
        },
        "conditional_refits": fits,
        "finite_reference_check": {
            "printed_2007_reference_volume_unit": "12.615 cm3/atom",
            "conditional_interpretation": "12.615 cm3/mol of FeS",
            "qualification": (
                "The 2010 HTML prints a dimensionally implausible unit. This "
                "diagnostic assumes cm3/mol FeS; the original Ohfuji paper "
                "has not been inspected, so no authoritative unit correction "
                "or original-source validation is claimed."
            ),
            "conditional_reference_volume_angstrom3": repaired_2007_vr,
            "quoted_parameters": [36.0, 306.0, 3.81],
            "quoted_errors": [1.7, 17.0, 0.28],
            "sata_curve_transformed_parameters": transformed,
            "max_reference_transform_pressure_difference_gpa": float(
                np.max(np.abs(pressure(v, transformed, repaired_2007_vr) - pressure(v)))
            ),
            "conditional_2007_replay": summary(finite_replay),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    if args.check:
        stored = json.loads(OUTPUT.read_text(encoding="utf-8"))

        # Numerical results may vary slightly between BLAS/SciPy platforms.
        def compare(a, b):
            if isinstance(a, dict):
                return a.keys() == b.keys() and all(compare(a[k], b[k]) for k in a)
            if isinstance(a, list):
                return len(a) == len(b) and all(compare(x, y) for x, y in zip(a, b))
            if isinstance(a, float):
                return bool(np.isclose(a, b, rtol=2e-6, atol=2e-8))
            return a == b

        if not compare(report, stored):
            raise SystemExit("Stored cold audit differs from replay")
        print("Sata FeS-VI cold audit verified")
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(OUTPUT)


if __name__ == "__main__":
    main()
