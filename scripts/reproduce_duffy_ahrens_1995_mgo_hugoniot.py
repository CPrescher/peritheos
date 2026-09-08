#!/usr/bin/env python3
"""Reproduce Duffy and Ahrens' (1995) Table 3 MgO Hugoniot fit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from peritheos.fitting import fit_linear_us_up

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mgo-duffy-ahrens-1995-table3-hugoniot.csv"
)
DATASET_SHA256 = "22a966924d96973ff66c8dcb37f51c7b51e426feb507dbe5b5c4665d075ee31c"
SOURCE_PDF_SHA256 = "daaf0a41440bb125dd0a062ca0582c358ceabf04f39c08c40cda2ea73dc67b2b"
SOURCE_URL = (
    "https://duffy.princeton.edu/sites/g/files/toruqf616/files/"
    "duffy_ahrens_seismo_1995.pdf"
)
PUBLISHED = {"c0": 6.87, "s": 1.24}
PUBLISHED_STANDARD_ERRORS = {"c0": 0.10, "s": 0.04}
V0 = 75.1559422385
RHO0 = 3.562


def sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest for *path*."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rows(path: Path = DATASET) -> list[dict[str, str]]:
    """Load the lossless Table 3 transcription."""
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 4:
        raise ValueError(f"Expected four Table 3 EOS rows, found {len(rows)}")
    return rows


def _fit(rows: list[dict[str, str]], *, mode: str) -> Any:
    selected = [row for row in rows if row["used_in_published_hugoniot_fit"] == "1"]
    up = np.asarray([float(row["particle_velocity_km_s"]) for row in selected])
    us = np.asarray([float(row["shock_velocity_km_s"]) for row in selected])
    kwargs: dict[str, Any] = {}
    if mode in {"shock_velocity_wls", "errors_in_variables"}:
        kwargs.update(
            shock_velocity_sigma=np.asarray(
                [
                    float(row["shock_velocity_standard_deviation_km_s"])
                    for row in selected
                ]
            ),
            absolute_sigma=True,
        )
    if mode == "errors_in_variables":
        kwargs["particle_velocity_sigma"] = np.asarray(
            [
                float(row["particle_velocity_standard_deviation_km_s"])
                for row in selected
            ]
        )
    return fit_linear_us_up(up, us, V0=V0, rho0=RHO0, P0=0.0, **kwargs)


def _fit_summary(result: Any) -> dict[str, Any]:
    covariance = np.asarray(result.covariance, dtype=float)
    standard_errors = np.sqrt(np.diag(covariance))
    correlation = covariance[0, 1] / (standard_errors[0] * standard_errors[1])
    return {
        "parameters": {name: float(result.parameters[name]) for name in ("c0", "s")},
        "standard_errors": {
            name: float(error) for name, error in zip(("c0", "s"), standard_errors)
        },
        "covariance": covariance.tolist(),
        "parameter_correlation_c0_s": float(correlation),
        "shock_velocity_rmse_km_s": float(
            np.sqrt(np.mean(np.square(result.residuals)))
        ),
        "chi_square": float(result.chi_square),
        "reduced_chi_square": float(result.reduced_chi_square),
        "degrees_of_freedom": int(result.degrees_of_freedom),
        "adjusted_particle_velocity_km_s": np.asarray(
            result.adjusted_particle_velocity, dtype=float
        ).tolist(),
    }


def reproduce(
    dataset: Path = DATASET, source_pdf: Path | None = None
) -> dict[str, Any]:
    """Return fit and Rankine-Hugoniot checks for the four Table 3 states."""
    dataset_digest = sha256(dataset)
    if dataset == DATASET and dataset_digest != DATASET_SHA256:
        raise ValueError("Bundled Table 3 CSV checksum does not match the audited file")
    if source_pdf is not None and sha256(source_pdf) != SOURCE_PDF_SHA256:
        raise ValueError("Source PDF checksum does not match the audited author copy")

    rows = load_rows(dataset)
    selected = [row for row in rows if row["used_in_published_hugoniot_fit"] == "1"]
    pressure_differences = []
    density_differences = []
    for row in selected:
        rho0 = float(row["initial_density_g_cm3"])
        up = float(row["particle_velocity_km_s"])
        us = float(row["shock_velocity_km_s"])
        pressure_differences.append(rho0 * us * up - float(row["axial_stress_gpa"]))
        density_differences.append(
            rho0 * us / (us - up) - float(row["shocked_density_g_cm3"])
        )

    fits = {
        mode: _fit_summary(_fit(rows, mode=mode))
        for mode in (
            "ordinary_least_squares",
            "shock_velocity_wls",
            "errors_in_variables",
        )
    }
    eiv = fits["errors_in_variables"]
    return {
        "source": {
            "doi": "10.1029/94JB02065",
            "url": SOURCE_URL,
            "source_pdf_sha256": SOURCE_PDF_SHA256,
            "location": "Table 3 and Equation (7), journal pages 533-534",
            "rights": (
                "The author-hosted article states copyright 1995 by the American "
                "Geophysical Union and gives no open data license. The article/PDF "
                "is not redistributed. The bundled CSV is a Peritheos-created factual "
                "transcription with an original machine-readable arrangement."
            ),
        },
        "dataset": {
            "path": str(dataset.relative_to(ROOT))
            if dataset.is_relative_to(ROOT)
            else str(dataset),
            "sha256": dataset_digest,
            "rows": len(rows),
        },
        "selection": {
            "included_rows": len(selected),
            "shot_numbers": [int(row["shot_number"]) for row in selected],
            "rule": (
                "Use all four final shock states in Table 3. Do not mix in the "
                "reverse/forward-impact sound-velocity shots from Tables 1-2, the "
                "single elastic-precursor value, single-crystal literature points, "
                "or the separately reduced 300 K hydrostat."
            ),
            "phase": (
                "Untransformed polycrystalline B1 MgO principal-Hugoniot branch; "
                "the source reports no transition for these 14-133 GPa states."
            ),
        },
        "published": {
            "parameters": PUBLISHED,
            "standard_errors": PUBLISHED_STANDARD_ERRORS,
            "uncertainty_convention": "one standard deviation in the last digit(s)",
        },
        "fits": fits,
        "published_rounding_reproduced": {
            name: round(eiv["parameters"][name], 2) == PUBLISHED[name]
            and round(eiv["standard_errors"][name], 2)
            == PUBLISHED_STANDARD_ERRORS[name]
            for name in ("c0", "s")
        },
        "rankine_hugoniot_checks": {
            "max_abs_pressure_difference_gpa": float(
                np.max(np.abs(pressure_differences))
            ),
            "max_abs_density_difference_g_cm3": float(
                np.max(np.abs(density_differences))
            ),
            "note": "Differences reflect only the precision printed in Table 3.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=DATASET)
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = (
        json.dumps(reproduce(args.dataset, args.source_pdf), indent=2, sort_keys=True)
        + "\n"
    )
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
