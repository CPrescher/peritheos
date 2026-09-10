#!/usr/bin/env python3
"""Reproduce the Gleason et al. (2008) epsilon-FeOOH BM2 record.

Primary source: doi:10.2138/am.2008.2942, manuscript page 7, Figure 3b,
and MSA depository item AM-08-056 Table 2.

The source does not print an epsilon-specific temperature-reduction equation.
This script makes the numerically recovered operation explicit and prints the
raw-hot-volume control that rules out treating the measurements as a 300 K
isotherm.
"""

from __future__ import annotations

import csv
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
    / "epsilon-feooh-gleason-2008-deposit-table2-pvt.csv"
)

TR = 300.0
ALPHA0 = 2.3e-5
V0 = 66.3
PUBLISHED_K0 = 158.0


def load_data() -> dict[str, np.ndarray]:
    """Load all 49 depository rows at their printed precision."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "volume": np.array([float(row["cell_volume_a3"]) for row in rows]),
        "temperature": np.array(
            [float(row["temperature_celsius"]) + 273.15 for row in rows]
        ),
    }


def reduce_volume_to_reference(
    volume: np.ndarray, temperature: np.ndarray
) -> np.ndarray:
    """Apply the recovered constant-expansivity reduction to 300 K."""
    return volume * np.exp(-ALPHA0 * (temperature - TR))


def fit_fixed_v0(volume: np.ndarray, pressure: np.ndarray):
    """Fit the sole free BM2 coefficient with source weights unavailable."""
    return fit_rt_eos(
        BM2,
        volume,
        pressure,
        initial={"K0": PUBLISHED_K0},
        fixed={"V0": V0},
        bounds={"K0": (1.0, 500.0)},
        absolute_sigma=False,
    )


def g_vs_g_volume(volume: np.ndarray, pressure: np.ndarray) -> float:
    """Return the zero-stress volume from a linear Jeanloz G-versus-g check."""
    g = 0.5 * ((V0 / volume) ** (2.0 / 3.0) - 1.0)
    normalized_stress = pressure / (3.0 * (1.0 + 2.0 * g) ** 2.5)
    slope, intercept = np.polyfit(g, normalized_stress, 1)
    zero_stress_g = -intercept / slope
    return float(V0 / (1.0 + 2.0 * zero_stress_g) ** 1.5)


def pressure_rmse(model: BM2, volume: np.ndarray, pressure: np.ndarray) -> float:
    """Return the ordinary pressure-residual RMSE in GPa."""
    residual = np.asarray(model.pressure(volume), dtype=float) - pressure
    return float(np.sqrt(np.mean(residual**2)))


def main() -> None:
    """Print the recovered fit and the deliberately incorrect raw-data control."""
    data = load_data()
    corrected_volume = reduce_volume_to_reference(data["volume"], data["temperature"])
    corrected_fit = fit_fixed_v0(corrected_volume, data["pressure"])
    raw_fit = fit_fixed_v0(data["volume"], data["pressure"])
    published = BM2(V0=V0, K0=PUBLISHED_K0)

    print(f"Table 2 epsilon-FeOOH observations: {data['pressure'].size}")
    print(
        "Measured temperature range: "
        f"{data['temperature'].min():.2f}-{data['temperature'].max():.2f} K"
    )
    print(
        f"Unweighted G-versus-g V0 check: {g_vs_g_volume(data['volume'], data['pressure']):.12f} A^3"
    )
    print("\nRecovered 300 K volume-reduction protocol:")
    print(f"  alpha0 fixed: {ALPHA0:.7g} K^-1")
    print(f"  V0 fixed: {V0:.6g} A^3")
    print(f"  K0 refit: {corrected_fit.parameters['K0']:.12f} GPa")
    print(f"  K0 standard error: {corrected_fit.standard_errors['K0']:.12f} GPa")
    print(
        "  published-curve pressure RMSE: "
        f"{pressure_rmse(published, corrected_volume, data['pressure']):.12f} GPa"
    )
    print("\nRaw-hot-volume control (not the source-faithful fit):")
    print(f"  K0 refit: {raw_fit.parameters['K0']:.12f} GPa")
    print(
        "  refit pressure RMSE: "
        f"{pressure_rmse(BM2(**raw_fit.parameters), data['volume'], data['pressure']):.12f} GPa"
    )
    print(
        "  published-curve pressure RMSE: "
        f"{pressure_rmse(published, data['volume'], data['pressure']):.12f} GPa"
    )


if __name__ == "__main__":
    main()
