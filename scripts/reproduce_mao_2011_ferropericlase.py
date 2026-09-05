#!/usr/bin/env python3
"""Reproduce the Mao et al. (2011) fp25 spin-branch audit.

The source does not publish HS/LS row flags.  This script reconstructs the
unclipped low-spin fraction from the two corrected published end-member curves,
uses a five-percent plateau tolerance, and refits only those pure-state inliers.
DOI: 10.1029/2011GL049915; correction: 10.1029/2011GL050814.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS

ROOT = Path(__file__).parents[1]
DATASET = ROOT / "peritheos/data/datasets/ferropericlase-fp25-mao-2011-table-s1-pvt.csv"
PLATEAU_TOLERANCE = 0.05
PUBLISHED = {
    "high_spin": {
        "V0": 76.34,
        "K0": 162.0,
        "K0_prime": 4.0,
        "alpha0": 3.76e-5,
        "dK_dT": -0.017,
    },
    "low_spin": {
        "V0": 74.4,
        "K0": 166.0,
        "K0_prime": 4.0,
        "alpha0": 3.37e-5,
        "dK_dT": -0.014,
    },
}


@dataclass(frozen=True)
class ThermalRefit:
    parameters: np.ndarray
    parameter_errors: np.ndarray
    reduced_chi_square: float
    pressure_rmse_gpa: float
    max_absolute_pressure_residual_gpa: float
    observations: int


def load_rows() -> np.ndarray:
    """Load the complete normalized Table S1 resource."""
    with DATASET.open(newline="", encoding="utf-8") as stream:
        source = list(csv.DictReader(stream))
    return np.array(
        [
            (
                int(row["source_order"]),
                float(row["temperature_k"]),
                float(row["pressure_gpa"]),
                float(row["pressure_standard_deviation_gpa"]),
                float(row["volume_a3_conventional_cell"]),
                float(row["volume_standard_deviation_a3"]),
                0.0
                if not row["temperature_standard_deviation_k"]
                else float(row["temperature_standard_deviation_k"]),
                float(row["reconstructed_low_spin_fraction_unclipped"]),
                int(row["used_in_diagnostic_high_spin_refit"]),
                int(row["used_in_diagnostic_low_spin_refit"]),
            )
            for row in source
        ],
        dtype=[
            ("source_order", int),
            ("temperature_k", float),
            ("pressure_gpa", float),
            ("pressure_sigma_gpa", float),
            ("volume_a3", float),
            ("volume_sigma_a3", float),
            ("temperature_sigma_k", float),
            ("stored_low_spin_fraction", float),
            ("high_spin_flag", int),
            ("low_spin_flag", int),
        ],
    )


def model(branch: str, alpha0: float, dK_dT: float) -> ThermalReferenceStateEOS:
    """Build the corrected published thermal BM3 branch in cell-volume units."""
    values = PUBLISHED[branch]
    return ThermalReferenceStateEOS(
        BM3(values["V0"], values["K0"], values["K0_prime"]),
        Tr=300.0,
        alpha0=alpha0,
        dK_dT=dK_dT,
        thermal_expansion_law="constant",
        reference_volume_law="integrated_expansivity",
    )


def published_model(branch: str) -> ThermalReferenceStateEOS:
    values = PUBLISHED[branch]
    return model(branch, values["alpha0"], values["dK_dT"])


def reconstructed_low_spin_fraction(rows: np.ndarray) -> np.ndarray:
    """Return Equation 6's volume-derived n_LS without clipping source scatter."""
    high_spin = published_model("high_spin")
    low_spin = published_model("low_spin")
    pressure = rows["pressure_gpa"]
    temperature = rows["temperature_k"]
    high_volume = np.asarray(high_spin.volume(pressure, temperature), dtype=float)
    low_volume = np.asarray(low_spin.volume(pressure, temperature), dtype=float)
    return (high_volume - rows["volume_a3"]) / (high_volume - low_volume)


def branch_mask(rows: np.ndarray, branch: str) -> np.ndarray:
    fractions = reconstructed_low_spin_fraction(rows)
    if branch == "high_spin":
        return fractions <= PLATEAU_TOLERANCE
    if branch == "low_spin":
        return fractions >= 1.0 - PLATEAU_TOLERANCE
    raise ValueError("branch must be 'high_spin' or 'low_spin'")


def effective_pressure_sigma(
    eos: ThermalReferenceStateEOS, rows: np.ndarray
) -> np.ndarray:
    """Propagate reported P/V errors and the stated 50 K hot-series error."""
    volume = rows["volume_a3"]
    temperature = rows["temperature_k"]
    volume_step = 1.0e-5
    temperature_step = 0.1
    dp_dv = (
        np.asarray(eos.pressure(volume + volume_step, temperature), dtype=float)
        - np.asarray(eos.pressure(volume - volume_step, temperature), dtype=float)
    ) / (2.0 * volume_step)
    dp_dt = (
        np.asarray(eos.pressure(volume, temperature + temperature_step), dtype=float)
        - np.asarray(eos.pressure(volume, temperature - temperature_step), dtype=float)
    ) / (2.0 * temperature_step)
    sigma = np.sqrt(
        rows["pressure_sigma_gpa"] ** 2
        + (rows["volume_sigma_a3"] * dp_dv) ** 2
        + (rows["temperature_sigma_k"] * dp_dt) ** 2
    )
    return np.where(sigma > 0.0, sigma, 1.0)


def fit_thermal_branch(rows: np.ndarray, branch: str) -> ThermalRefit:
    """Refit alpha0 and dK/dT with iterated effective-variance weights."""
    selected = rows[branch_mask(rows, branch)]
    published = PUBLISHED[branch]
    values = np.array([published["alpha0"], published["dK_dT"]])
    optimization = None
    for _ in range(100):
        sigma = effective_pressure_sigma(model(branch, *values), selected)

        def residual(candidate: np.ndarray) -> np.ndarray:
            calculated = np.asarray(
                model(branch, *candidate).pressure(
                    selected["volume_a3"], selected["temperature_k"]
                ),
                dtype=float,
            )
            return (calculated - selected["pressure_gpa"]) / sigma

        optimization = least_squares(
            residual,
            values,
            bounds=([-1.0e-4, -0.2], [1.0e-4, 0.2]),
            x_scale=[1.0e-5, 0.01],
            xtol=1.0e-14,
            ftol=1.0e-14,
            gtol=1.0e-14,
            max_nfev=5000,
        )
        if np.allclose(optimization.x, values, rtol=0.0, atol=1.0e-13):
            values = optimization.x
            break
        values = optimization.x
    assert optimization is not None
    final_sigma = effective_pressure_sigma(model(branch, *values), selected)
    calculated = np.asarray(
        model(branch, *values).pressure(
            selected["volume_a3"], selected["temperature_k"]
        ),
        dtype=float,
    )
    pressure_residual = calculated - selected["pressure_gpa"]
    weighted_residual = pressure_residual / final_sigma
    covariance = np.linalg.pinv(optimization.jac.T @ optimization.jac)
    return ThermalRefit(
        parameters=values,
        parameter_errors=np.sqrt(np.maximum(np.diag(covariance), 0.0)),
        reduced_chi_square=float(
            np.sum(weighted_residual**2) / (len(selected) - len(values))
        ),
        pressure_rmse_gpa=float(np.sqrt(np.mean(pressure_residual**2))),
        max_absolute_pressure_residual_gpa=float(np.max(np.abs(pressure_residual))),
        observations=len(selected),
    )


def fit_static_branch(rows: np.ndarray, branch: str) -> np.ndarray:
    """Refit the 300 K BM3 parameters with K0'=4 and source fixed values."""
    selected = rows[(rows["temperature_k"] == 300.0) & branch_mask(rows, branch)]
    if branch == "high_spin":
        values = np.array([PUBLISHED[branch]["K0"]])

        def static_model(candidate: np.ndarray) -> BM3:
            return BM3(76.34, candidate[0], 4.0)

    elif branch == "low_spin":
        values = np.array([PUBLISHED[branch]["V0"], PUBLISHED[branch]["K0"]])

        def static_model(candidate: np.ndarray) -> BM3:
            return BM3(candidate[0], candidate[1], 4.0)

    else:
        raise ValueError("branch must be 'high_spin' or 'low_spin'")
    for _ in range(100):
        eos = static_model(values)
        step = 1.0e-5
        dp_dv = (
            np.asarray(eos.pressure(selected["volume_a3"] + step), dtype=float)
            - np.asarray(eos.pressure(selected["volume_a3"] - step), dtype=float)
        ) / (2.0 * step)
        sigma = np.sqrt(
            selected["pressure_sigma_gpa"] ** 2
            + (selected["volume_sigma_a3"] * dp_dv) ** 2
        )
        sigma = np.where(sigma > 0.0, sigma, 1.0)
        result = least_squares(
            lambda candidate: (
                np.asarray(
                    static_model(candidate).pressure(selected["volume_a3"]),
                    dtype=float,
                )
                - selected["pressure_gpa"]
            )
            / sigma,
            values,
            xtol=1.0e-14,
            ftol=1.0e-14,
            gtol=1.0e-14,
            max_nfev=5000,
        )
        if np.allclose(result.x, values, rtol=0.0, atol=1.0e-13):
            return result.x
        values = result.x
    return values


def published_curve_residuals(rows: np.ndarray, branch: str) -> np.ndarray:
    selected = rows[branch_mask(rows, branch)]
    return (
        np.asarray(
            published_model(branch).pressure(
                selected["volume_a3"], selected["temperature_k"]
            ),
            dtype=float,
        )
        - selected["pressure_gpa"]
    )


def main() -> None:
    rows = load_rows()
    fractions = reconstructed_low_spin_fraction(rows)
    print(
        f"rows={len(rows)}; fraction range={fractions.min():.6f}..{fractions.max():.6f}"
    )
    for branch in ("high_spin", "low_spin"):
        static = fit_static_branch(rows, branch)
        thermal = fit_thermal_branch(rows, branch)
        published_residual = published_curve_residuals(rows, branch)
        print(
            f"{branch}: N={thermal.observations}; static={static}; "
            f"alpha0={thermal.parameters[0]:.10g}; "
            f"dK_dT={thermal.parameters[1]:.10g}; "
            f"published RMSE={np.sqrt(np.mean(published_residual**2)):.10g} GPa; "
            f"published max|residual|={np.max(np.abs(published_residual)):.10g} GPa"
        )


if __name__ == "__main__":
    main()
