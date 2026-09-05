#!/usr/bin/env python3
"""Reproduce Kubo et al.'s (2006) MgGeO3 post-perovskite EOS fits.

Primary source: doi:10.1029/2006GL025686, paragraph 17, Table 2, and
Supporting Information Table S1.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mggeo3-post-perovskite-kubo-2006-table-s1-pv.csv"
)

PREFERRED = {"V0": 179.2, "K0": 207.0, "K0_prime": 4.4}
SENSITIVITY = {"V0": 175.9, "K0": 245.0, "K0_prime": 4.0}


def load_rows() -> list[dict[str, str]]:
    """Load all 25 source observations, including the three excluded rows."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def source_bm3(
    volume: np.ndarray | float,
    v0: float,
    k0: float,
    k0_prime: float,
) -> np.ndarray:
    """Evaluate the standard third-order Birch-Murnaghan pressure."""
    volume_array = np.asarray(volume, dtype=float)
    eta = (v0 / volume_array) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def selected_arrays() -> tuple[np.ndarray, np.ndarray]:
    """Return source pressures and derived volumes for the 22 selected rows."""
    selected = [row for row in load_rows() if row["fit_included"] == "1"]
    pressures = np.array([float(row["pressure_gpa"]) for row in selected])
    volumes = np.array([float(row["volume_a3"]) for row in selected])
    return pressures, volumes


def fit_with_fixed_derivative(k0_prime: float) -> np.ndarray:
    """Refit V0 and K0 with the selected K0' held fixed."""
    pressures, volumes = selected_arrays()
    result = least_squares(
        lambda beta: source_bm3(volumes, beta[0], beta[1], k0_prime) - pressures,
        x0=(179.0, 210.0),
        bounds=((160.0, 100.0), (200.0, 400.0)),
    )
    return result.x


def curve_metrics(parameters: dict[str, float]) -> tuple[float, float]:
    """Return pressure RMSE and maximum residual for one printed curve."""
    pressures, volumes = selected_arrays()
    calculated = source_bm3(
        volumes,
        parameters["V0"],
        parameters["K0"],
        parameters["K0_prime"],
    )
    residuals = calculated - pressures
    return (
        math.sqrt(float(np.mean(residuals**2))),
        float(np.max(np.abs(residuals))),
    )


def platinum_reduction_metrics() -> tuple[float, float]:
    """Recalculate source pressures from rounded Pt 111 d-spacings."""
    platinum_document = get_material_document("platinum")
    platinum = Material.from_eosmat(
        platinum_document,
        record_identifiers=["platinum_holmes_1989_vinet_1"],
    ).get_eos_record("platinum_holmes_1989_vinet_1")
    rows = load_rows()
    observed = np.array([float(row["pressure_gpa"]) for row in rows])
    d111 = np.array([float(row["d_pt111_angstrom"]) for row in rows])
    platinum_volumes = (np.sqrt(3.0) * d111) ** 3
    calculated = np.asarray(platinum.pressure(platinum_volumes), dtype=float)
    residuals = calculated - observed
    return (
        math.sqrt(float(np.mean(residuals**2))),
        float(np.max(np.abs(residuals))),
    )


def main() -> None:
    """Print source-curve, refit, selection, and pressure-scale checks."""
    rows = load_rows()
    included = [row for row in rows if row["fit_included"] == "1"]
    excluded = [row for row in rows if row["fit_included"] == "0"]
    print(f"Supporting Information Table S1 rows: {len(rows)}")
    print(
        f"Published fit selection: {len(included)} included, {len(excluded)} excluded"
    )
    print(
        "Excluded pressures: "
        + ", ".join(f"{float(row['pressure_gpa']):g}" for row in excluded)
        + " GPa"
    )

    for label, parameters in (
        ("preferred fixed-K0-prime=4.4 BM3", PREFERRED),
        ("nonpreferred fixed-K0-prime=4 BM2 sensitivity", SENSITIVITY),
    ):
        rmse, maximum = curve_metrics(parameters)
        refit = fit_with_fixed_derivative(parameters["K0_prime"])
        print(f"\n{label}:")
        print(f"  published-curve pressure RMSE: {rmse:.12f} GPa")
        print(f"  published-curve maximum residual: {maximum:.12f} GPa")
        print(f"  unweighted refit V0: {refit[0]:.12f} A^3")
        print(f"  unweighted refit K0: {refit[1]:.12f} GPa")

    pt_rmse, pt_maximum = platinum_reduction_metrics()
    print("\nHolmes et al. (1989) Pt-pressure recalculation from rounded d111:")
    print(f"  pressure RMSE: {pt_rmse:.12f} GPa")
    print(f"  maximum difference: {pt_maximum:.12f} GPa")


if __name__ == "__main__":
    main()
