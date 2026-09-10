#!/usr/bin/env python3
"""Audit the selected Ismailova et al. (2016) Fe-bridgmanite P-V rows.

The primary Figure 3 fit contains substantially more markers than the four
selected crystallographic states printed in supplementary Table S2.  This
script therefore reports checkpoint residuals and explicitly diagnostic fits;
it does not claim to reconstruct the source regression.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM2
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "fe088sio3-bridgmanite-ismailova-2016-table-s2-selected-crystallography.csv"
)
PUBLISHED = {"V0": 178.98, "K0": 190.0}
BOUNDS = {"V0": (100.0, 250.0), "K0": (1.0, 1000.0)}


def load_data() -> dict[str, np.ndarray]:
    """Load the lossless transcription of the four Table S2 states."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.asarray([float(row["pressure_gpa"]) for row in rows]),
        "pressure_sigma": np.asarray(
            [float(row["pressure_uncertainty_gpa"]) for row in rows]
        ),
        "volume": np.asarray(
            [float(row["volume_a3_conventional_cell"]) for row in rows]
        ),
        "volume_sigma": np.asarray(
            [float(row["volume_uncertainty_a3"]) for row in rows]
        ),
    }


def published_curve_diagnostic() -> dict[str, object]:
    """Evaluate the published BM2 at the four exact Table S2 states."""
    data = load_data()
    residuals = np.asarray(BM2(**PUBLISHED).pressure(data["volume"])) - data["pressure"]
    return {
        "observations": int(data["pressure"].size),
        "parameters": {**PUBLISHED, "K0_prime": 4.0},
        "pressure_residuals_gpa_in_source_order": residuals.tolist(),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "max_abs_pressure_residual_gpa": float(np.max(np.abs(residuals))),
    }


def diagnostic_refit(*, errors_in_variables: bool) -> dict[str, object]:
    """Fit the selected subset under an explicit, non-source audit objective."""
    data = load_data()
    fit = fit_rt_eos(
        BM2,
        volume=data["volume"],
        pressure=data["pressure"],
        initial=PUBLISHED,
        bounds=BOUNDS,
        pressure_sigma=data["pressure_sigma"] if errors_in_variables else None,
        volume_sigma=data["volume_sigma"] if errors_in_variables else None,
        absolute_sigma=errors_in_variables,
        max_nfev=5000,
    )
    residuals = np.asarray(fit.residuals, dtype=float)
    return {
        "objective": (
            "errors in pressure and volume using printed uncertainties"
            if errors_in_variables
            else "unweighted pressure residuals at printed volumes"
        ),
        "observations": int(data["pressure"].size),
        "free_parameters": list(fit.free_parameters),
        "fixed_parameters": ["K0_prime"],
        "parameters": {**fit.parameters, "K0_prime": 4.0},
        "standard_errors": fit.standard_errors,
        "pressure_residual_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "reduced_chi_square": float(fit.reduced_chi_square),
        "degrees_of_freedom": int(fit.degrees_of_freedom),
        "success": bool(fit.success),
    }


def reproduce() -> dict[str, object]:
    """Return deterministic diagnostics and the exact source-fit blocker."""
    return {
        "published_curve": published_curve_diagnostic(),
        "selected_table_s2_refits": {
            "unweighted_pressure": diagnostic_refit(errors_in_variables=False),
            "errors_in_variables": diagnostic_refit(errors_in_variables=True),
        },
        "source_fit": {
            "status": "not_refittable",
            "fixed_parameters": ["K0_prime"],
            "reported_exclusions": "none stated",
            "blockers": [
                "Figure 3A plots substantially more compression/decompression markers than the four selected states in Table S2",
                "the source supplies no row-level inclusion or exclusion flags",
                "the residual objective and numerical weights are not reported",
                "row-level Ne calibrant lattice parameters needed to recalculate pressure are absent",
            ],
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
