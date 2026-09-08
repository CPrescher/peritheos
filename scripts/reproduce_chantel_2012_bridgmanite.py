#!/usr/bin/env python3
"""Audit Chantel et al. (2012) acoustic and thermoelastic claims."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from peritheos.acoustics import EulerianFiniteStrainAcoustic
from peritheos.fitting import fit_acoustic_finite_strain

ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "peritheos/data/datasets/bridgmanite-chantel-2012-table1-density-velocity.csv"
)

PUBLISHED_CURRENT_STUDY = {
    "K_S0": 247.0,
    "K_S0_prime": 4.5,
    "G0": 176.0,
    "G0_prime": 1.6,
}
PUBLISHED_CURRENT_STUDY_ERRORS = {
    "K_S0": 4.0,
    "K_S0_prime": 0.2,
    "G0": 2.0,
    "G0_prime": 0.1,
}
PUBLISHED_COMBINED = {
    "K_S0": 252.0,
    "K_S0_prime": 4.1,
    "G0": 175.0,
    "G0_prime": 1.7,
}
PUBLISHED_COMBINED_ERRORS = {
    "K_S0": 1.0,
    "K_S0_prime": 0.1,
    "G0": 1.0,
    "G0_prime": 0.1,
}


def bm3_pressure(
    rho: np.ndarray, rho0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate BM3 using density ratio rho/rho0 = V0/V."""
    eta = (rho / rho0) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _parameter_report(result, published, published_errors) -> dict[str, object]:
    return {
        name: {
            "published": published[name],
            "published_error": published_errors[name],
            "refit": result.parameters[name],
            "refit_standard_error": result.standard_errors[name],
            "difference": result.parameters[name] - published[name],
            "within_combined_two_sigma": bool(
                abs(result.parameters[name] - published[name])
                <= 2.0 * np.hypot(result.standard_errors[name], published_errors[name])
            ),
        }
        for name in published
    }


def reproduce() -> dict[str, object]:
    """Return source-equation, refit, and model-scope diagnostics."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    anchor = next(row for row in rows if row["reference_density_anchor"] == "1")
    acoustic = [row for row in rows if row["acoustic_velocity_fit_included"] == "1"]
    high_temperature = [
        row for row in rows if row["thermoelastic_validation_only"] == "1"
    ]

    rho0 = float(anchor["density_g_cm3"])
    rho = np.array([float(row["density_g_cm3"]) for row in acoustic])
    rho_sigma = np.array([float(row["density_sigma_g_cm3"]) for row in acoustic])
    vp = np.array([float(row["vp_km_s"]) for row in acoustic])
    vp_sigma = np.array([float(row["vp_sigma_km_s"]) for row in acoustic])
    vs = np.array([float(row["vs_km_s"]) for row in acoustic])
    vs_sigma = np.array([float(row["vs_sigma_km_s"]) for row in acoustic])

    unweighted = fit_acoustic_finite_strain(
        rho,
        vp,
        vs,
        rho0=rho0,
        initial=PUBLISHED_CURRENT_STUDY,
        max_nfev=5000,
    )
    diagonal_eiv = fit_acoustic_finite_strain(
        rho,
        vp,
        vs,
        rho0=rho0,
        initial=PUBLISHED_CURRENT_STUDY,
        density_sigma=rho_sigma,
        compressional_velocity_sigma=vp_sigma,
        shear_velocity_sigma=vs_sigma,
        absolute_sigma=True,
        max_nfev=5000,
    )

    combined_model = EulerianFiniteStrainAcoustic(rho0=rho0, **PUBLISHED_COMBINED)
    combined_vp, combined_vs = combined_model.velocities(rho)
    room = [row for row in rows if float(row["temperature_k"]) == 300.0]
    room_rho = np.array([float(row["density_g_cm3"]) for row in room])
    pressure = np.array([float(row["pressure_gpa"]) for row in room])
    predicted_pressure = bm3_pressure(room_rho, 100.387 / 24.42, 252.0, 4.1)
    pressure_residual = predicted_pressure - pressure

    return {
        "observations": len(rows),
        "reference_density_anchor": {"rows": 1, "rho0_g_cm3": rho0},
        "current_study_acoustic_fit": {
            "observations": len(acoustic),
            "objective": "simultaneous unweighted Vp and Vs residuals",
            "parameters": _parameter_report(
                unweighted,
                PUBLISHED_CURRENT_STUDY,
                PUBLISHED_CURRENT_STUDY_ERRORS,
            ),
            "vp_rmse_km_s": float(
                np.sqrt(np.mean(unweighted.compressional_velocity_residuals**2))
            ),
            "vs_rmse_km_s": float(
                np.sqrt(np.mean(unweighted.shear_velocity_residuals**2))
            ),
        },
        "diagonal_uncertainty_eiv_sensitivity": {
            "objective": (
                "Vp, Vs, and latent-density residuals divided by the printed "
                "marginal uncertainties; zero correlation assumed"
            ),
            "parameters": _parameter_report(
                diagonal_eiv,
                PUBLISHED_CURRENT_STUDY,
                PUBLISHED_CURRENT_STUDY_ERRORS,
            ),
            "reduced_chi_square": diagonal_eiv.reduced_chi_square,
            "maximum_absolute_density_adjustment_g_cm3": float(
                np.max(np.abs(diagonal_eiv.density_corrections))
            ),
        },
        "preferred_combined_coefficients_on_chantel_rows": {
            "scope": (
                "diagnostic only: the Li and Zhang (2005) rows included by the "
                "published combined fit are not republished by Chantel et al."
            ),
            "vp_rmse_km_s": float(np.sqrt(np.mean((combined_vp - vp) ** 2))),
            "vs_rmse_km_s": float(np.sqrt(np.mean((combined_vs - vs) ** 2))),
        },
        "table3_pressure_surface_at_300_k": {
            "observations": len(room),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(pressure_residual**2))),
            "maximum_absolute_pressure_residual_gpa": float(
                np.max(np.abs(pressure_residual))
            ),
        },
        "thermoelastic_fit_assessment": {
            "high_temperature_validation_rows": len(high_temperature),
            "fit_performed_by_source": False,
            "finding": (
                "The source states that the high-temperature data were insufficient "
                "to refine thermal parameters. Theta0, gamma0, q, and eta were "
                "adopted from Xu et al. (2008), and the 700/1200 K velocities only "
                "tested the resulting Stixrude-type thermoelastic model."
            ),
        },
        "unavailable_source_fit_information": [
            "numerical Li and Zhang (2005) observations used in the combined fit",
            "the exact residual definition and numerical least-squares weights",
            "within-row density/Vp/Vs error correlations",
            "the four-parameter acoustic covariance matrix",
            "a confidence convention for Chantel et al. parenthetical errors",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
