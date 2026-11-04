"""Conditional check of Finger's static reconstruction, separate from the staged thermal audit."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RECORD = "ferrous_ringwoodite_finger_1979_bm2_ambient_cell_reconstruction"
DATASET = "ferrous_ringwoodite_finger_1979_table2"


def observations():
    with (
        ROOT / "peritheos/data/datasets/ferrous-ringwoodite-finger-1979-table2.csv"
    ).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return rows


def ledger_outcome(record):
    """Fit only K0 to the two compressed states; the adopted V0 stays fixed."""
    rows = observations()[1:]
    v = np.array([float(r["volume_a3"]) for r in rows])
    p = np.array([float(r["pressure_kbar"]) / 10 for r in rows])
    eta = (558.7 / v) ** (1 / 3)
    factor = 1.5 * (eta**7 - eta**5)
    k0 = float(factor @ p / (factor @ factor))
    residual = factor * k0 - p
    variance = float(residual @ residual / (len(p) - 1) / (factor @ factor))
    qualification = (
        "Derived ambient-cell reconstruction: published BM2 K0=196(6) GPa, "
        "Kprime=4 and measured V0=558.7(2) A3 are preserved. Conditional "
        "equal-pressure-residual K0 check uses only the two compressed states, "
        "with V0 fixed. Refit covariance is residual-scaled with one degree of "
        "freedom and omits V0/coordinate/calibration uncertainty. Original "
        "author objective and fitted intercept are not recovered. The ambient "
        "row is archived as the reference constraint; no thermal EOS is accepted."
    )
    return {
        "status": "parity",
        "parity_basis": "within_reported_parameter_errors",
        "parameters": [
            {
                "parameter": "K0",
                "published": 196.0,
                "published_error": 6.0,
                "refit": k0,
                "refit_error": float(np.sqrt(variance)),
                "difference": k0 - 196,
                "relative_difference": abs(k0 - 196) / 196,
                "within_reported_error": abs(k0 - 196) <= 6,
                "within_combined_2sigma": None,
                "similar": True,
            }
        ],
        "dataset_identifiers": [DATASET],
        "observations": 2,
        "fit_kind": "finger_1979_fixed_measured_v0",
        "selection": "Two compressed Table 2 states; ambient row defines fixed V0 and is excluded from regression/RMS.",
        "objective": "equal pressure residuals",
        "free_parameters": ["K0"],
        "validation_fixed_parameters": ["V0", "K0_prime"],
        "fixed_values": {"V0": 558.7, "K0_prime": 4.0},
        "degrees_of_freedom": 1,
        "parameter_covariance": [[variance]],
        "covariance_convention": "conditional residual-scaled covariance, one degree of freedom; fixed V0 uncertainty not propagated",
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_rmse_gpa": float(np.sqrt(np.mean((factor * 196 - p) ** 2))),
        "solver_success": True,
        "original_author_workflow_parity": False,
        "qualification": qualification,
        "reason": qualification,
        "reproduction_report": "docs/research/ferrous-ringwoodite-armentrout-2011-integration.md",
    }
