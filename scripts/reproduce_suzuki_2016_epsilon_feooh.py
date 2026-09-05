#!/usr/bin/env python3
"""Reproduce the Suzuki (2016) epsilon-FeOOH P-V-T EOS.

Primary source: doi:10.2465/jmps.160719c, Equation 1 and Tables 1-2.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.fitting import fit_joint_eos

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "epsilon-feooh-suzuki-2016-table1-pvt.csv"
)

PUBLISHED = {
    "rt_eos.K0": (135.0, 3.0),
    "rt_eos.K0_prime": (6.1, 0.9),
    "alpha0": (2.6e-5, 0.7e-5),
    "alpha1": (1.0e-7, 0.3e-7),
    "dK_dT": (-0.05, 0.02),
}


def load_data() -> dict[str, np.ndarray]:
    """Load the 33 Table 1 rows at their printed precision."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: np.array(
            [float(row[name]) if row[name] else math.nan for row in rows],
            dtype=float,
        )
        for name in rows[0]
        if name != "source_row"
    }


def source_pressure(volume: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Evaluate Suzuki's printed equations independently of Peritheos."""
    tr = 300.0
    delta = temperature - tr
    v0_t = 66.278 * np.exp(2.6e-5 * delta + 0.5e-7 * delta**2)
    k0_t = 135.0 - 0.05 * delta
    eta = (v0_t / volume) ** (1.0 / 3.0)
    return 1.5 * k0_t * (eta**7 - eta**5) * (1.0 + 0.75 * (6.1 - 4.0) * (eta**2 - 1.0))


def published_model() -> ThermalReferenceStateEOS:
    """Construct the source-faithful Peritheos parameterization."""
    return ThermalReferenceStateEOS(
        BM3(V0=66.278, K0=135.0, K0_prime=6.1),
        Tr=300.0,
        alpha0=2.6e-5,
        alpha1=1.0e-7,
        dK_dT=-0.05,
        thermal_expansion_law="linear_reference_temperature",
        reference_volume_law="integrated_expansivity",
    )


def fit_table(data: dict[str, np.ndarray]):
    """Fit the 32 finite-pressure rows with V0 fixed by the ambient row."""
    selected = np.isfinite(data["pressure_sigma_gpa"])
    return fit_joint_eos(
        ThermalReferenceStateEOS,
        BM3,
        data["cell_volume_a3"][selected],
        data["temperature_k"][selected],
        data["pressure_gpa"][selected],
        initial={name: value for name, (value, _) in PUBLISHED.items()},
        fixed={"rt_eos.V0": 66.278, "Tr": 300.0},
        configuration={
            "thermal_expansion_law": "linear_reference_temperature",
            "reference_volume_law": "integrated_expansivity",
        },
        bounds={
            "rt_eos.K0": (1.0, 500.0),
            "rt_eos.K0_prime": (-5.0, 30.0),
            "alpha0": (-1.0e-4, 2.0e-4),
            "alpha1": (-1.0e-6, 1.0e-6),
            "dK_dT": (-0.5, 0.5),
        },
        pressure_sigma=data["pressure_sigma_gpa"][selected],
        volume_sigma=data["cell_volume_sigma_a3"][selected],
        absolute_sigma=True,
    )


def main() -> None:
    """Print curve reproduction and independent-refit diagnostics."""
    data = load_data()
    observed = data["pressure_gpa"]
    independent = source_pressure(data["cell_volume_a3"], data["temperature_k"])
    model = published_model()
    peritheos = np.asarray(
        model.pressure(data["cell_volume_a3"], data["temperature_k"]), dtype=float
    )
    residual = peritheos - observed

    print(f"Table 1 observations: {observed.size}")
    print(
        "Independent equation vs Peritheos maximum difference: "
        f"{np.max(np.abs(independent - peritheos)):.3e} GPa"
    )
    print(f"Published-curve pressure RMSE: {np.sqrt(np.mean(residual**2)):.12f} GPa")
    print(f"Published-curve maximum residual: {np.max(np.abs(residual)):.12f} GPa")
    print(f"700 K benchmark at 62.63 A^3: {model.pressure(62.63, 700.0):.12f} GPa")

    fit = fit_table(data)
    selected = np.isfinite(data["pressure_sigma_gpa"])
    raw_residual = (
        np.asarray(
            fit.model.pressure(
                data["cell_volume_a3"][selected],
                data["temperature_k"][selected],
            ),
            dtype=float,
        )
        - observed[selected]
    )
    print("\nPeritheos pressure-volume errors-in-variables refit:")
    for name, (published, sigma) in PUBLISHED.items():
        value = fit.parameters[name]
        print(
            f"  {name}: {value:.12g} "
            f"(published {published:.12g}; delta {((value - published) / sigma):+.3f} sigma)"
        )
    print(f"  reduced chi-square: {fit.reduced_chi_square:.12f}")
    print(f"  raw pressure RMSE: {np.sqrt(np.mean(raw_residual**2)):.12f} GPa")


if __name__ == "__main__":
    main()
