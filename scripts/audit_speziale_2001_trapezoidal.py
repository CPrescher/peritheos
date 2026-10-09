"""Sensitivity to trapezoidal quadrature, without refitting source coefficients.

python -m scripts.audit_speziale_2001_trapezoidal [--check]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import trapezoid

from scripts.reproduce_speziale_2001_mgo import PARAMETERS, RESIDUALS, bm3, energy
from scripts.validate_speziale_2001_equations import (
    closed_form_pressure,
    closed_form_theta,
    source_curve_rows,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/speziale-2001-trapezoidal-integration.json"
INTERVALS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 1024)


def trapezoidal_theta(volume, intervals, coordinate="log_volume"):
    """Integrate N5 with source Eq. 11 on two explicit uniform grids."""
    x = np.asarray(volume, dtype=float) / PARAMETERS["V0"]
    fractions = np.linspace(0.0, 1.0, intervals + 1)
    if coordinate == "log_volume":
        nodes = np.log(x[..., None]) * fractions
        gamma = 1.524 * np.exp(1.65 / 11.8 * np.expm1(11.8 * nodes))
        integral = trapezoid(gamma, nodes, axis=-1)
    elif coordinate == "volume":
        nodes = 1.0 + (x[..., None] - 1.0) * fractions
        gamma = 1.524 * np.exp(1.65 / 11.8 * (nodes**11.8 - 1.0))
        integral = trapezoid(gamma / nodes, nodes, axis=-1)
    else:
        raise ValueError(coordinate)
    return 773.0 * np.exp(-integral)


def trapezoidal_energy(theta, temperature, intervals):
    """Diagnostic Debye quadrature; not specified in the supplied author reply."""
    theta, temperature = np.broadcast_arrays(theta, temperature)
    y = theta / temperature
    nodes = y[..., None] * np.linspace(0.0, 1.0, intervals + 1)
    integrand = np.divide(
        nodes**3, np.expm1(nodes), out=np.zeros_like(nodes), where=nodes != 0.0
    )
    return 18.0 * R * temperature * trapezoid(integrand, nodes, axis=-1) / y**3


def trial_pressure(volume, temperature, intervals, theta_coordinate, debye_trapezoidal):
    volume, temperature = np.broadcast_arrays(volume, temperature)
    x = volume / 74.71
    gamma = 1.524 * np.exp(1.65 / 11.8 * (x**11.8 - 1.0))
    theta = (
        closed_form_theta(volume)
        if theta_coordinate is None
        else trapezoidal_theta(volume, intervals, theta_coordinate)
    )

    def energies(t):
        if debye_trapezoidal:
            return trapezoidal_energy(theta, t, intervals)
        return energy(theta, t)

    return (
        bm3(volume)
        + gamma
        * (energies(temperature) - energies(300.0))
        / (volume * N_A * 1e-30 / 4.0)
        / 1e9
    )


def audit():
    with RESIDUALS.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    later_v = np.array([float(r["mgo_volume_a3"]) for r in rows])
    later_t = np.array([float(r["temperature_k"]) for r in rows])
    targets = np.array([float(r["reported_pressure_gpa"]) for r in rows])
    fei = np.array([r["source"] == "Fei 2004 Table 1" for r in rows])
    curves = [r for r in source_curve_rows() if r["figure"] == "10"]
    curve_v = 74.71 * np.array([float(r["volume_ratio"]) for r in curves])
    curve_t = np.array([float(r["temperature_k"]) for r in curves])
    curve_targets = np.array([float(r["pressure_gpa"]) for r in curves])
    grid_v = np.repeat(74.71 * np.array([0.6, 0.64, 0.7, 0.8, 0.9, 1.0, 1.02]), 5)
    grid_t = np.tile(np.array([300.0, 500.0, 1100.0, 2000.0, 3663.0]), 7)

    def reference(v, t):
        return np.array([closed_form_pressure(a, b) for a, b in zip(v, t)])

    grid_ref = reference(grid_v, grid_t)
    later_ref = reference(later_v, later_t)
    curve_ref = reference(curve_v, curve_t)
    trials = []
    for coordinate, debye in [
        ("log_volume", False),
        ("volume", False),
        (None, True),
        ("log_volume", True),
    ]:
        for intervals in INTERVALS:
            grid = trial_pressure(grid_v, grid_t, intervals, coordinate, debye)
            later = trial_pressure(later_v, later_t, intervals, coordinate, debye)
            curve = trial_pressure(curve_v, curve_t, intervals, coordinate, debye)
            trials.append(
                {
                    "theta_grid": coordinate or "analytic",
                    "debye_quadrature": "trapezoidal"
                    if debye
                    else "converged_gauss_legendre",
                    "intervals": intervals,
                    "grid35_max_pressure_change_gpa": float(
                        np.max(np.abs(grid - grid_ref))
                    ),
                    "later28_max_pressure_change_gpa": float(
                        np.max(np.abs(later - later_ref))
                    ),
                    "later28_pressure_change_range_gpa": [
                        float(np.min(later - later_ref)),
                        float(np.max(later - later_ref)),
                    ],
                    "fei2004_rmse_gpa": float(
                        np.sqrt(np.mean((later[fei] - targets[fei]) ** 2))
                    ),
                    "hirose2006_max_abs_residual_gpa": float(
                        np.max(np.abs(later[~fei] - targets[~fei]))
                    ),
                    "figure10_max_pressure_change_gpa": float(
                        np.max(np.abs(curve - curve_ref))
                    ),
                    "figure10_max_abs_residual_gpa": float(
                        np.max(np.abs(curve - curve_targets))
                    ),
                }
            )
    return {
        "scope": "Quadrature sensitivity only; all published coefficients fixed. No original thermal fit is replayed.",
        "author_confirmation_scope": "Numerical integration of N5 using S11, generally trapezoidal; original grid unspecified. Debye-energy quadrature is not specified in the supplied reply.",
        "grid_definition": "Each state integrates theta from V0 to its volume using N equal intervals in log(V/V0) or V/V0. Debye diagnostic uses N equal intervals from 0 to theta/T.",
        "parameters": PARAMETERS,
        "reference": {
            "theta": "independent analytic exponential-integral solution",
            "pressure": "independent adaptive Debye integral and finite-strain BM3",
            "fei2004_rmse_gpa": float(
                np.sqrt(np.mean((later_ref[fei] - targets[fei]) ** 2))
            ),
            "hirose2006_max_abs_residual_gpa": float(
                np.max(np.abs(later_ref[~fei] - targets[~fei]))
            ),
            "figure10_max_abs_residual_gpa": float(
                np.max(np.abs(curve_ref - curve_targets))
            ),
        },
        "trials": trials,
        "limits": "Uniform grids tested here are sensitivity scenarios, not a recovered historical grid. Later Fei/Hirose pressures are not original Speziale fit targets. Digitized Figure 10 points have graphical extraction limits. Quadrature convergence does not establish original fitting parity.",
        "calibration_input_sha256": hashlib.sha256(RESIDUALS.read_bytes()).hexdigest(),
        "source_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.check:
        from scripts.reproduce_fei_2007_gold import check_saved

        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), result)
    else:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["reference"], indent=2))
