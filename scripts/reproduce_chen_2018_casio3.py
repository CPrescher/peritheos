#!/usr/bin/env python3
"""Audit Chen et al.'s seven-point 300 K tetragonal CaSiO3 Vinet fit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos/data/datasets/ca-perovskite-tetragonal-chen-2018-table1-compression.csv"
)

V0_PER_FORMULA = 46.3
Z_I4MCM = 4.0
PUBLISHED_V0 = V0_PER_FORMULA * Z_I4MCM
PUBLISHED_K0 = 223.0
FIXED_K0_PRIME = 4.0


def vinet_pressure(volume, v0=PUBLISHED_V0, k0=PUBLISHED_K0, k0_prime=FIXED_K0_PRIME):
    """Return Vinet pressure for a scalar or array volume."""

    volume = np.asarray(volume, dtype=float)
    x = np.cbrt(volume / v0)
    return 3.0 * k0 * (1.0 - x) / x**2 * np.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))


def vinet_volume(pressure, v0, k0, k0_prime=FIXED_K0_PRIME):
    """Invert Vinet pressure on the physical compressed branch."""

    return np.array(
        [
            brentq(
                lambda volume: vinet_pressure(volume, v0, k0, k0_prime) - value,
                0.2 * v0,
                v0,
            )
            for value in np.asarray(pressure, dtype=float)
        ]
    )


def load_table():
    with DATASET.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    def values(name):
        return np.array([float(row[name]) for row in rows])

    pressure = values("pressure_pt_ye_gpa")
    i4_a = values("i4mcm_lattice_a_angstrom")
    i4_c = values("i4mcm_lattice_c_angstrom")
    i4_volume = values("i4mcm_volume_a3_conventional_cell")
    i4_volume_sigma = values("i4mcm_volume_1sigma_a3_conventional_cell")
    p4_a = values("p4mmm_lattice_a_angstrom")
    p4_c = values("p4mmm_lattice_c_angstrom")
    p4_volume_per_formula = values("p4mmm_volume_a3_conventional_cell")

    # Table 1 prints lattice parameters, not volumes. Check every derived value.
    np.testing.assert_allclose(i4_volume, i4_a**2 * i4_c, rtol=0.0, atol=5e-13)
    np.testing.assert_allclose(
        p4_volume_per_formula, p4_a**2 * p4_c, rtol=0.0, atol=5e-13
    )
    assert all(row["used_in_published_fit"] == "1" for row in rows)
    return rows, pressure, i4_volume, i4_volume_sigma, p4_volume_per_formula


def _fit_pressure(volume, pressure, initial=(PUBLISHED_V0, PUBLISHED_K0)):
    result = least_squares(
        lambda pars: vinet_pressure(volume, pars[0], pars[1]) - pressure,
        initial,
        bounds=([100.0, 1.0], [300.0, 1000.0]),
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=20_000,
    )
    return result.x


def _fit_volume(volume, pressure, sigma=None, initial=(PUBLISHED_V0, PUBLISHED_K0)):
    def residuals(pars):
        residual = vinet_volume(pressure, pars[0], pars[1]) - volume
        return residual if sigma is None else residual / sigma

    result = least_squares(
        residuals,
        initial,
        bounds=([100.0, 1.0], [300.0, 1000.0]),
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=20_000,
    )
    return result.x


def _fit_fixed_v0(volume, pressure, v0_per_formula=45.58):
    v0 = Z_I4MCM * v0_per_formula
    result = least_squares(
        lambda pars: vinet_pressure(volume, v0, pars[0], pars[1]) - pressure,
        [290.0, 2.3],
        bounds=([1.0, -10.0], [1000.0, 20.0]),
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=20_000,
    )
    return np.array([v0, *result.x])


def _summary(parameters, volume, pressure):
    v0, k0 = parameters[:2]
    k0_prime = parameters[2] if len(parameters) == 3 else FIXED_K0_PRIME
    residual = vinet_pressure(volume, v0, k0, k0_prime) - pressure
    return {
        "V0_a3_conventional_cell": float(v0),
        "V0_a3_per_formula_unit": float(v0 / Z_I4MCM),
        "K0_gpa": float(k0),
        "K0_prime": float(k0_prime),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "pressure_max_abs_residual_gpa": float(np.max(np.abs(residual))),
    }


def reproduce():
    rows, pressure, volume, volume_sigma, p4_volume_per_formula = load_table()
    published = _summary(np.array([PUBLISHED_V0, PUBLISHED_K0]), volume, pressure)
    published["pressure_residuals_gpa"] = [
        float(value) for value in vinet_pressure(volume) - pressure
    ]

    pressure_fit = _fit_pressure(volume, pressure)
    volume_fit = _fit_volume(volume, pressure)
    volume_weighted_fit = _fit_volume(volume, pressure, volume_sigma)
    fixed_v0_fit = _fit_fixed_v0(volume, pressure)

    # Diagnostic only: the paper explicitly identifies I4/mcm as the fit branch.
    both_volumes = np.concatenate([volume, Z_I4MCM * p4_volume_per_formula])
    both_pressures = np.concatenate([pressure, pressure])
    both_fit = _fit_pressure(both_volumes, both_pressures)

    return {
        "source_protocol": {
            "equation": "Vinet",
            "observations": len(rows),
            "included_runs": [row["run_id"] for row in rows],
            "excluded_runs": [],
            "structural_branch": "I4/mcm",
            "pressure_calibration": (
                "Pt EOS of Ye et al. (2017), DOI 10.1002/2016JB013811"
            ),
            "fixed_parameters": {"K0_prime": FIXED_K0_PRIME},
            "reported_weights": None,
            "reported_pressure_uncertainties": None,
            "reported_fit_software": None,
            "table_uncertainty_convention": "2 sigma",
        },
        "published_curve": published,
        "source_scope_refits": {
            "unweighted_pressure_residual": _summary(pressure_fit, volume, pressure),
            "unweighted_volume_residual": _summary(volume_fit, volume, pressure),
            "propagated_volume_1sigma_weighted": _summary(
                volume_weighted_fit, volume, pressure
            ),
        },
        "sensitivity_diagnostics": {
            "source_reported_fixed_V0_45_58_free_K0_K0_prime": _summary(
                fixed_v0_fit, volume, pressure
            ),
            "unsupported_both_refinement_branches_unweighted_pressure": _summary(
                both_fit, both_volumes, both_pressures
            ),
        },
        "conclusion": (
            "The complete rounded I4/mcm table supports a direct source-scope "
            "Vinet refit, but exact coefficient/covariance reproduction is not "
            "possible because the residual direction, pressure errors, weights, "
            "fit software, and unrounded inputs are not published."
        ),
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
