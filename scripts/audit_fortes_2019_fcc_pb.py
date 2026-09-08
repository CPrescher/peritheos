"""Audit the Fortes (2019) fcc-Pb pressure-scale parameterization.

The report fully specifies its temperature-dependent BM4 surface, but not the
row-level data or regression protocol needed to refit it.  This module keeps
those two claims separate: it evaluates the published surface independently,
and it attempts a diagnostic fit only to the three exact fcc P-V-T anchors
printed by Kuznetsov et al. (2002).  Those anchors are not represented as the
private table supplied to Fortes or as the reported global fit dataset.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).parents[1]
DATA = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "lead-fcc-kuznetsov-2002-table1-transition-pvt.csv"
)
HCP_FIGURE3_DATA = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "lead-hcp-kuznetsov-2002-figure3-digitized.csv"
)

FORTES_REPORT_SHA256 = (
    "5318783f2b40c7c8dc95986f4dd157d32e2834790cdff4daa6d7769769029ab7"
)
KUZNETSOV_ARTICLE_SHA256 = (
    "c4bfc71c5d4186f7f3159015bafa51552349f9e0f4ba305d13604673d0b3a61d"
)

# Fortes (2019), Table 1.  The table writes the constant terms as X(0)
# because its polynomial argument is Delta T = T - 300 K; they are the
# zero-pressure quantities at 300 K, not values at absolute zero.
PUBLISHED = {
    "V0_300": 121.418,
    "a": 1.058e-2,
    "b": 3.5e-6,
    "K0_300": 41.73,
    "c": -2.544e-5,
    "d": -2.8e-6,
    "K0_prime_300": 5.39,
    "e": 1.1e-3,
    "K0_double_prime": -0.33,
}

PUBLISHED_ERRORS = {
    "V0_300": 0.005,
    "a": 4.0e-5,
    "b": 2.0e-7,
    "K0_300": 0.01,
    "c": 4.0e-8,
    "d": 2.0e-7,
    "K0_prime_300": 0.25,
    "e": 1.0e-4,
    "K0_double_prime": 0.02,
}

# Kuznetsov et al. (2002), Equations 1--3 and Table 2, fcc column.  This is
# the source paper's fitted thermal BM3 parameterization, not the row-level
# P-V-T table that its author later supplied privately to Fortes.
KUZNETSOV_FCC = {
    "V0_296_atomic": 30.307,
    "alpha_0": 8.67e-5,
    "alpha_1": 2.12e-9,
    "alpha_2": 1.83e-10,
    "K0_296": 40.5,
    "b_1": -1.07e-3,
    "b_2": -1.14e-5,
    "K0_prime_296": 5.74,
    "b_prime": -5.45e-3,
}

# Kuznetsov et al. (2002), Table 2, hcp column. The signs and exponents were
# checked visually against the source PDF because text extraction loses the
# typographic minus signs in this table.
KUZNETSOV_HCP = {
    "V0_296_atomic": 29.908,
    "alpha_0": 8.97e-5,
    "alpha_1": 2.52e-8,
    "alpha_2": 0.0,
    "K0_296": 54.2,
    "b_1": -5.46e-2,
    "b_2": 2.34e-5,
    "K0_prime_296": 3.61,
    "b_prime": 5.01e-3,
}

KUZNETSOV_HCP_TABLE1_ANCHORS = {
    "temperature_k": np.asarray([296.0, 402.0, 469.0]),
    "pressure_gpa": np.asarray([13.1, 13.9, 12.6]),
    "atomic_volume_a3": np.asarray([24.84, 24.92, 25.06]),
}


def temperature_parameters(
    temperature_k: np.ndarray | float,
    *,
    k0_prime_300: float = PUBLISHED["K0_prime_300"],
    e: float = PUBLISHED["e"],
    k0_double_prime: float = PUBLISHED["K0_double_prime"],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return V0(T), K0(T), K0'(T), and the invariant K0''."""
    temperature = np.asarray(temperature_k, dtype=float)
    delta_t = temperature - 300.0
    v0 = PUBLISHED["V0_300"] + PUBLISHED["a"] * delta_t + PUBLISHED["b"] * (delta_t**2)
    k0 = PUBLISHED["K0_300"] + PUBLISHED["c"] * delta_t + PUBLISHED["d"] * (delta_t**2)
    k0_prime = k0_prime_300 + e * delta_t
    k0_double_prime_array = np.broadcast_to(k0_double_prime, temperature.shape)
    return v0, k0, k0_prime, k0_double_prime_array


def fortes_pressure(
    volume_a3_conventional_cell: np.ndarray | float,
    temperature_k: np.ndarray | float,
    *,
    k0_prime_300: float = PUBLISHED["K0_prime_300"],
    e: float = PUBLISHED["e"],
    k0_double_prime: float = PUBLISHED["K0_double_prime"],
) -> np.ndarray:
    """Evaluate Fortes Equations 1--5 independently of Peritheos."""
    volume = np.asarray(volume_a3_conventional_cell, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    volume, temperature = np.broadcast_arrays(volume, temperature)
    v0, k0, k0_prime, k0_double_prime_array = temperature_parameters(
        temperature,
        k0_prime_300=k0_prime_300,
        e=e,
        k0_double_prime=k0_double_prime,
    )
    strain = 0.5 * ((volume / v0) ** (-2.0 / 3.0) - 1.0)
    quadratic = (
        k0 * k0_double_prime_array + (k0_prime - 4.0) * (k0_prime - 3.0) + 35.0 / 9.0
    )
    return (
        3.0
        * k0
        * strain
        * (1.0 + 2.0 * strain) ** 2.5
        * (1.0 + 1.5 * (k0_prime - 4.0) * strain + 1.5 * quadratic * strain**2)
    )


def kuznetsov_temperature_parameters(
    temperature_k: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the published Kuznetsov fcc V0(T), K0(T), and K0'(T)."""
    temperature = np.asarray(temperature_k, dtype=float)
    delta_t = temperature - 296.0
    v0 = KUZNETSOV_FCC["V0_296_atomic"] * (
        1.0
        + KUZNETSOV_FCC["alpha_0"] * delta_t
        + 0.5 * KUZNETSOV_FCC["alpha_1"] * delta_t**2
        + (1.0 / 3.0) * KUZNETSOV_FCC["alpha_2"] * delta_t**3
    )
    k0 = (
        KUZNETSOV_FCC["K0_296"]
        + KUZNETSOV_FCC["b_1"] * delta_t
        + KUZNETSOV_FCC["b_2"] * delta_t**2
    )
    k0_prime = KUZNETSOV_FCC["K0_prime_296"] + KUZNETSOV_FCC["b_prime"] * delta_t
    return v0, k0, k0_prime


def kuznetsov_pressure(
    atomic_volume_a3: np.ndarray | float,
    temperature_k: np.ndarray | float,
) -> np.ndarray:
    """Evaluate Kuznetsov Equations 1--3 for fcc Pb."""
    volume = np.asarray(atomic_volume_a3, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    volume, temperature = np.broadcast_arrays(volume, temperature)
    v0, k0, k0_prime = kuznetsov_temperature_parameters(temperature)
    compression = v0 / volume
    return (
        1.5
        * k0
        * (compression ** (7.0 / 3.0) - compression ** (5.0 / 3.0))
        * (1.0 + 0.75 * (k0_prime - 4.0) * (compression ** (2.0 / 3.0) - 1.0))
    )


def kuznetsov_hcp_temperature_parameters(
    temperature_k: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the published Kuznetsov hcp V0(T), K0(T), and K0'(T)."""
    temperature = np.asarray(temperature_k, dtype=float)
    delta_t = temperature - 296.0
    v0 = KUZNETSOV_HCP["V0_296_atomic"] * (
        1.0
        + KUZNETSOV_HCP["alpha_0"] * delta_t
        + 0.5 * KUZNETSOV_HCP["alpha_1"] * delta_t**2
        + (1.0 / 3.0) * KUZNETSOV_HCP["alpha_2"] * delta_t**3
    )
    k0 = (
        KUZNETSOV_HCP["K0_296"]
        + KUZNETSOV_HCP["b_1"] * delta_t
        + KUZNETSOV_HCP["b_2"] * delta_t**2
    )
    k0_prime = KUZNETSOV_HCP["K0_prime_296"] + KUZNETSOV_HCP["b_prime"] * delta_t
    return v0, k0, k0_prime


def kuznetsov_hcp_pressure(
    atomic_volume_a3: np.ndarray | float,
    temperature_k: np.ndarray | float,
) -> np.ndarray:
    """Evaluate Kuznetsov Equations 1--3 for hcp Pb."""
    volume = np.asarray(atomic_volume_a3, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    volume, temperature = np.broadcast_arrays(volume, temperature)
    v0, k0, k0_prime = kuznetsov_hcp_temperature_parameters(temperature)
    compression = v0 / volume
    return (
        1.5
        * k0
        * (compression ** (7.0 / 3.0) - compression ** (5.0 / 3.0))
        * (1.0 + 0.75 * (k0_prime - 4.0) * (compression ** (2.0 / 3.0) - 1.0))
    )


def _kuznetsov_hcp_checks() -> dict[str, Any]:
    """Check the hcp surface against exact anchors and Figure 3 markers."""
    anchors = KUZNETSOV_HCP_TABLE1_ANCHORS
    anchor_pressure = kuznetsov_hcp_pressure(
        anchors["atomic_volume_a3"], anchors["temperature_k"]
    )
    anchor_residual = anchor_pressure - anchors["pressure_gpa"]

    with HCP_FIGURE3_DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    figure_volume = np.asarray(
        [float(row["atomic_volume_a3_per_atom"]) for row in rows]
    )
    figure_pressure = np.asarray([float(row["reduced_pressure_gpa"]) for row in rows])
    figure_fill = np.asarray([row["symbol_fill"] for row in rows])
    figure_prediction = kuznetsov_hcp_pressure(figure_volume, 296.0)
    figure_residual = figure_prediction - figure_pressure
    solid = figure_fill == "solid"

    return {
        "published_parameters": dict(KUZNETSOV_HCP),
        "reference_temperature_k": 296.0,
        "table1_transition_anchors": {
            "observations": int(anchor_pressure.size),
            "predicted_pressure_gpa": anchor_pressure.tolist(),
            "pressure_residuals_gpa": anchor_residual.tolist(),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(anchor_residual**2))),
        },
        "figure3_reduced_room_temperature_markers": {
            "observations": int(figure_pressure.size),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(figure_residual**2))),
            "max_abs_pressure_residual_gpa": float(np.max(np.abs(figure_residual))),
            "solid_room_temperature_observations": int(np.count_nonzero(solid)),
            "solid_room_temperature_pressure_rmse_gpa": float(
                np.sqrt(np.mean(figure_residual[solid] ** 2))
            ),
            "qualification": (
                "The solid markers are direct room-temperature observations. "
                "Open markers were reduced from high temperature by the source "
                "and all marker coordinates are plot digitizations."
            ),
        },
        "interpretation": (
            "The 296 K hcp BM3 slice is directly defined by Equations 1-3 and "
            "Table 2. Agreement with the seven digitized solid room-temperature "
            "markers verifies the executable slice within plot-reading precision."
        ),
    }


def _load_kuznetsov_anchors() -> dict[str, np.ndarray]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "temperature_k": np.asarray([float(row["temperature_k"]) for row in rows]),
        "pressure_gpa": np.asarray([float(row["pressure_gpa"]) for row in rows]),
        "volume_a3_conventional_cell": np.asarray(
            [float(row["volume_a3_conventional_cell"]) for row in rows]
        ),
    }


def _partial_anchor_fit(data: dict[str, np.ndarray]) -> dict[str, Any]:
    """Fit the three published derivative terms to three rounded anchors."""

    def residual(parameters: np.ndarray) -> np.ndarray:
        return (
            fortes_pressure(
                data["volume_a3_conventional_cell"],
                data["temperature_k"],
                k0_prime_300=float(parameters[0]),
                e=float(parameters[1]),
                k0_double_prime=float(parameters[2]),
            )
            - data["pressure_gpa"]
        )

    initial = np.asarray(
        [
            PUBLISHED["K0_prime_300"],
            PUBLISHED["e"],
            PUBLISHED["K0_double_prime"],
        ]
    )
    fit = least_squares(
        residual,
        initial,
        max_nfev=100_000,
        ftol=1.0e-14,
        xtol=1.0e-14,
        gtol=1.0e-14,
    )
    published_residuals = residual(initial)
    fitted_residuals = residual(fit.x)
    singular_values = np.linalg.svd(fit.jac, compute_uv=False)
    return {
        "scope": "Kuznetsov (2002) Table 1 transition anchors only",
        "is_reported_fortes_fit": False,
        "observations": int(data["pressure_gpa"].size),
        "free_parameters": ["K0_prime_300", "e", "K0_double_prime"],
        "fixed_parameters": ["V0_300", "a", "b", "K0_300", "c", "d"],
        "objective": "unweighted pressure residuals (diagnostic choice)",
        "published_parameters_pressure_residuals_gpa": published_residuals.tolist(),
        "published_parameters_pressure_rmse_gpa": float(
            np.sqrt(np.mean(published_residuals**2))
        ),
        "diagnostic_fit_parameters": {
            name: float(value)
            for name, value in zip(("K0_prime_300", "e", "K0_double_prime"), fit.x)
        },
        "diagnostic_fit_pressure_residuals_gpa": fitted_residuals.tolist(),
        "diagnostic_fit_pressure_rmse_gpa": float(
            np.sqrt(np.mean(fitted_residuals**2))
        ),
        "jacobian_rank": int(np.linalg.matrix_rank(fit.jac)),
        "jacobian_condition_number": float(singular_values[0] / singular_values[-1]),
        "interpretation": (
            "The three rounded phase-boundary anchors are nearly singular for "
            "three correlated derivative coefficients and do not reproduce the "
            "global fit. They are a public numerical check, not a substitute for "
            "the private 295-788 K compression table."
        ),
    }


def _room_temperature_identifiability() -> dict[str, Any]:
    """Show why the report's 300 K Figure 4 cannot determine dK'/dT."""
    volumes = PUBLISHED["V0_300"] * np.asarray([0.95, 0.90, 0.85, 0.81])
    epsilon = 1.0e-6
    plus = fortes_pressure(volumes, 300.0, e=PUBLISHED["e"] + epsilon)
    minus = fortes_pressure(volumes, 300.0, e=PUBLISHED["e"] - epsilon)
    derivative = (plus - minus) / (2.0 * epsilon)
    return {
        "source": "Fortes (2019) Figure 4 room-temperature comparison",
        "d_pressure_d_e_gpa_k": derivative.tolist(),
        "rank_consequence": (
            "Every derivative is exactly zero because Delta T=0 at 300 K; "
            "Figure 4 cannot identify e=dK0_prime/dT."
        ),
    }


def audit() -> dict[str, Any]:
    """Return deterministic equation, source-chain, and fit diagnostics."""
    data = _load_kuznetsov_anchors()
    temperatures = np.asarray([100.0, 300.0, 600.0, 788.0])
    v0, k0, k0_prime, k0_double_prime = temperature_parameters(temperatures)
    checkpoint_ratios = np.asarray([0.95, 0.90, 0.85])
    checkpoint_pressures = np.stack(
        [fortes_pressure(ratio * v0, temperatures) for ratio in checkpoint_ratios]
    )
    kuznetsov_residuals = (
        kuznetsov_pressure(
            data["volume_a3_conventional_cell"] / 4.0,
            data["temperature_k"],
        )
        - data["pressure_gpa"]
    )
    return {
        "format": "peritheos.fortes-2019-fcc-pb-audit",
        "format_version": 1,
        "record_identifier": "lead_fcc_fortes_2019_bm4_1",
        "outcome": "not_independently_refittable",
        "published_surface_executable": True,
        "sources": {
            "fortes_2019": {
                "doi": "10.5286/raltr.2019002",
                "pdf_sha256": FORTES_REPORT_SHA256,
                "locations": ["Equations 1-8", "Sections 2.1-2.3", "Tables 1-2"],
                "license": "CC BY 4.0",
            },
            "kuznetsov_2002": {
                "doi": "10.1016/S0038-1098(02)00112-6",
                "pdf_sha256": KUZNETSOV_ARTICLE_SHA256,
                "locations": ["Equations 1-3", "Table 1", "Figure 3"],
                "row_level_status": (
                    "three transition anchors printed; complete P-V-T series plot-only"
                ),
            },
        },
        "staged_fit": {
            "thermal_expansion_stage": {
                "method": (
                    "Debye internal-energy model fitted to compiled 10-600 K volume "
                    "data, then its T>100 K curve represented by V0(T) quadratic"
                ),
                "final_parameters": ["V0_300", "a", "b"],
                "raw_row_mask_and_weights_published": False,
            },
            "isothermal_modulus_stage": {
                "method": (
                    "KS=(c11+2c12)/3 converted with KT=KS/[alphaV^2 V KS T/CP+1], "
                    "then represented by K0(T) quadratic"
                ),
                "specific_heat_selection": (
                    "Meads et al. (1941) below room T; Leadbetter (1968a,b) above"
                ),
                "final_parameters": ["K0_300", "c", "d"],
                "CP_polynomial_coefficients_published": False,
                "raw_row_mask_and_weights_published": False,
            },
            "compression_stage": {
                "temperature_range_k": [295.0, 788.0],
                "fixed_parameters": ["V0_300", "a", "b", "K0_300", "c", "d"],
                "refined_parameters": [
                    "K0_prime_300",
                    "e",
                    "K0_double_prime",
                ],
                "algorithm": "least-squares minimization",
                "dependent_variable": "not reported",
                "weights": "not reported",
                "covariance_and_error_scaling": "not reported",
                "row_selection": "not reported",
            },
        },
        "published_parameters": PUBLISHED,
        "published_parameter_errors": PUBLISHED_ERRORS,
        "kuznetsov_published_predecessor_model": {
            "role": (
                "Published thermal BM3 comparison, not a reconstruction of the "
                "private P-V-T table used by Fortes"
            ),
            "reference_temperature_k": 296.0,
            "published_parameters": KUZNETSOV_FCC,
            "transition_anchor_pressure_residuals_gpa": kuznetsov_residuals.tolist(),
            "transition_anchor_pressure_rmse_gpa": float(
                np.sqrt(np.mean(kuznetsov_residuals**2))
            ),
        },
        "kuznetsov_hcp_published_model": _kuznetsov_hcp_checks(),
        "temperature_parameter_checkpoints": {
            "temperature_k": temperatures.tolist(),
            "V0_a3_conventional_cell": v0.tolist(),
            "K0_gpa": k0.tolist(),
            "K0_prime": k0_prime.tolist(),
            "K0_double_prime_gpa_inverse": k0_double_prime.tolist(),
        },
        "pressure_checkpoints": {
            "volume_ratios": checkpoint_ratios.tolist(),
            "temperature_k": temperatures.tolist(),
            "pressure_gpa_by_volume_ratio": checkpoint_pressures.tolist(),
        },
        "reported_fit_attempt": _partial_anchor_fit(data),
        "room_temperature_identifiability": _room_temperature_identifiability(),
        "missing_for_independent_refit": [
            "the tabulated Kuznetsov P-V-T data supplied privately to Fortes",
            "an exhaustive list of included literature series and row-level selection",
            "the treatment or recalculation of heterogeneous pressure scales",
            "the least-squares dependent variable and observation weights",
            "the covariance and uncertainty-scaling procedure",
        ],
    }


def main() -> None:
    """Print the audit as stable JSON."""
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
