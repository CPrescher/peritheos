#!/usr/bin/env python3
"""Reproduce the Koemets et al. (2023) FA50 perovskite BM2.

Primary source: Frontiers in Chemistry 11, 1258389,
doi:10.3389/fchem.2023.1258389.  Supplementary Table S3 supplies the
single-crystal unit-cell measurements.  The two Pnma observations inside the
reported 8 +/- 2 GPa structural transition are retained in the transcription
but excluded from the phase-specific fit; the fitted branch therefore starts
with the 10.8 GPa observation.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM2
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "fa50-bridgmanite-koemets-2023-supplement-table3.csv"
)
PUBLISHED = {"V0": 171.2, "K0": 221.0}


def load_data() -> dict[str, np.ndarray]:
    """Load the direct Supplementary Table S3 transcription."""
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "pressure_sigma": np.array(
            [float(row["pressure_uncertainty_gpa"]) for row in rows]
        ),
        "volume": np.array([float(row["volume_a3_conventional_cell"]) for row in rows]),
        "volume_sigma": np.array([float(row["volume_uncertainty_a3"]) for row in rows]),
        "selected": np.array(
            [row["used_in_bm2_fit"] == "1" for row in rows], dtype=bool
        ),
    }


def reproduce() -> dict[str, object]:
    """Return published-curve and two independent-refit diagnostics."""
    data = load_data()
    selected = data["selected"]
    pressure = data["pressure"][selected]
    volume = data["volume"][selected]
    pressure_sigma = data["pressure_sigma"][selected]
    volume_sigma = data["volume_sigma"][selected]

    published_pressure = BM2(**PUBLISHED).pressure(volume)
    published_residual = published_pressure - pressure

    def pressure_residual(parameters: np.ndarray) -> np.ndarray:
        return BM2(V0=parameters[0], K0=parameters[1]).pressure(volume) - pressure

    unweighted = least_squares(
        pressure_residual,
        np.array([PUBLISHED["V0"], PUBLISHED["K0"]]),
        bounds=([150.0, 1.0], [200.0, 1000.0]),
        x_scale="jac",
    )
    if not unweighted.success:
        raise RuntimeError(unweighted.message)

    errors_in_variables = fit_rt_eos(
        BM2,
        volume,
        pressure,
        PUBLISHED,
        pressure_sigma=pressure_sigma,
        volume_sigma=volume_sigma,
        absolute_sigma=True,
    )

    unweighted_residual = pressure_residual(unweighted.x)
    return {
        "source_rows": int(len(data["pressure"])),
        "selected_rows": int(np.count_nonzero(selected)),
        "excluded_transition_rows": int(np.count_nonzero(~selected)),
        "selected_pressure_range_gpa": [
            float(np.min(pressure)),
            float(np.max(pressure)),
        ],
        "published_parameters": PUBLISHED,
        "published_curve_pressure_rmse_gpa": float(
            np.sqrt(np.mean(published_residual**2))
        ),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(published_residual))
        ),
        "unweighted_pressure_refit": {
            "parameters": {"V0": float(unweighted.x[0]), "K0": float(unweighted.x[1])},
            "pressure_rmse_gpa": float(np.sqrt(np.mean(unweighted_residual**2))),
        },
        "errors_in_variables_refit": {
            "parameters": errors_in_variables.parameters,
            "standard_errors": errors_in_variables.standard_errors,
            "reduced_chi_square": errors_in_variables.reduced_chi_square,
        },
    }


def main() -> None:
    """Print deterministic JSON diagnostics."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
