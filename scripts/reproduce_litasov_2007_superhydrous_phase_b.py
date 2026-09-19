"""Independently audit Litasov (2007) Tables 2, 4–7, without EOS-library calls.

The physical-cell Debye normalization (70 atoms) independently checks the
production formula-unit/molar adapter. Regressions are diagnostics, not new
published coefficients. Source least-squares weights and covariance are absent.
"""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Boltzmann, Planck
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/superhydrous-phase-b-litasov-2007-table2.csv"
OUTPUT = ROOT / "docs/data/litasov-2007-superhydrous-phase-b-reproduction.json"
PREFIX = "superhydrous_phase_b_lt_litasov_2007_"
DATASET = PREFIX + "table2"
# Explicit primary-table transcription, deliberately separate from eosmat.
CASES = {
    "bm3_a89": ([623.38, 138.7, 4.9], [0.39, 3.0, 0.3]),
    "bm3_t03": ([623.34, 134.7, 6.0], [0.04, 1.8, 0.3]),
    "htbm_a89": (
        [623.47, 135.8, 5.3, 3.21e-5, 1.19e-8, -0.026],
        [0.37, 2.6, 0.2, 0.14e-5, 0.37e-8, 0.003],
    ),
    "htbm_t03": (
        [623.34, 132.7, 6.2, 3.38e-5, 1.26e-8, -0.027],
        [0.04, 1.5, 0.2, 0.15e-5, 0.37e-8, 0.003],
    ),
    "thermal_pressure_a89": (
        [623.50, 135.3, 5.3, 3.78e-5, -0.002],
        [0.36, 2.3, 0.2, 0.19e-5, 0.002],
    ),
    "mgd_free_debye_a89": (
        [623.53, 135.5, 5.3, 1.18, 1.74, 552.0],
        [0.35, 2.3, 0.2, 0.06, 0.45, 130.0],
    ),
    "mgd_elastic_debye_a89": (
        [623.53, 135.5, 5.3, 1.33, 2.03, 860.0],
        [0.35, 2.3, 0.2, 0.05, 0.35, None],
    ),
}
NODES, WEIGHTS = leggauss(48)
NODES, WEIGHTS = (NODES + 1) / 2, WEIGHTS / 2


def rows():
    with DATA.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def observations(suffix):
    data = rows()
    if suffix.startswith("bm3"):
        data = [r for r in data if float(r["temperature_k"]) == 300]
    scale = "t03" if suffix.endswith("t03") else "a89"
    return data, tuple(
        np.array([float(r[key]) for r in data])
        for key in ("volume_a3", "temperature_k", f"pressure_{scale}_gpa")
    )


def bm3(v, v0, k0, kp):
    """Eq. (1), correcting its repeated 7/3 exponent to standard BM3."""
    eta = (v0 / np.asarray(v)) ** (1 / 3)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def pressure(suffix, v, t, parameters=None, *, literal_sign=False, atoms=70):
    x = CASES[suffix][0] if parameters is None else parameters
    v, t = np.broadcast_arrays(np.asarray(v, dtype=float), np.asarray(t, dtype=float))
    v0, k0, kp = x[:3]
    if suffix.startswith("bm3"):
        return bm3(v, v0, k0, kp)
    if suffix.startswith("htbm"):
        a0, a1, dk = x[3:]
        vt = v0 * np.exp(a0 * (t - 300) + 0.5 * a1 * (t**2 - 300**2))
        return bm3(v, vt, k0 + dk * (t - 300), kp)
    if suffix.startswith("thermal"):
        a0, dk = x[3:]
        return bm3(v, v0, k0, kp) + (a0 * k0 + dk * np.log(v0 / v)) * (t - 300)
    g0, q, theta0 = x[3:]
    gamma = g0 * (v / v0) ** q
    theta = theta0 * np.exp((gamma - g0) / q * (1 if literal_sign else -1))

    def energy(temperature):
        y = theta / temperature
        z = y[..., None] * NODES
        return (
            9
            * atoms
            * Boltzmann
            * temperature
            * np.sum(WEIGHTS * z**3 / np.expm1(z), axis=-1)
            / y**2
        )

    return bm3(v, v0, k0, kp) + gamma * (energy(t) - energy(300)) / (v * 1e-21)


def parameter_names(suffix):
    names = ["V0", "K0", "K0_prime"]
    if suffix.startswith("htbm"):
        names += ["alpha0", "alpha1", "dK_dT"]
    elif suffix.startswith("thermal"):
        names += ["alpha0", "dK_dT_V"]
    elif suffix.startswith("mgd"):
        names += ["gamma0", "q", "theta0"]
    return names


def fit(suffix, *, literal_sign=False, weighted=False):
    data, (v, t, p) = observations(suffix)
    published, errors = CASES[suffix]
    initial = np.array(published)
    free = np.ones(len(initial), dtype=bool)
    if "elastic" in suffix:
        free[-1] = False
    lower, upper = [600.0, 90.0, 1.0], [650.0, 200.0, 10.0]
    if suffix.startswith("htbm"):
        lower += [-1e-4, -1e-7, -0.09]
        upper += [1e-4, 1e-7, 0.0]
    elif suffix.startswith("thermal"):
        lower += [0.0, -0.1]
        upper += [1e-4, 0.1]
    elif suffix.startswith("mgd"):
        lower += [0.1, 0.01, 100.0]
        upper += [3.0, 10.0, 2000.0]

    def unpack(values):
        full = initial.copy()
        full[free] = values
        return full

    def calculate(values, volume=v):
        return pressure(suffix, volume, t, unpack(values), literal_sign=literal_sign)

    weights = np.ones(len(v))
    if weighted:
        # Fixed effective variances at the published curve. Ambient pressure is
        # assigned exactly for this sensitivity only, with measured sigma_V.
        derivative = (
            calculate(initial[free], v + 0.001) - calculate(initial[free], v - 0.001)
        ) / 0.002
        sp = np.array([float(r["pressure_sigma_gpa"] or 0) for r in data])
        sv = np.array([float(r["volume_sigma_a3"]) for r in data])
        weights = np.sqrt(sp**2 + (derivative * sv) ** 2)
    result = least_squares(
        lambda values: (calculate(values) - p) / weights,
        initial[free],
        x_scale=np.maximum(abs(initial[free]), 1e-12),
        bounds=(np.array(lower)[free], np.array(upper)[free]),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    residual = calculate(result.x) - p
    # Scaling avoids ill-conditioned covariance inversion from mixed units.
    scale = np.maximum(abs(initial[free]), 1e-12)
    jac = result.jac * scale
    covariance = (
        np.linalg.inv(jac.T @ jac) * np.sum(result.fun**2) / (len(v) - len(result.x))
    )
    standard_errors = np.sqrt(np.diag(covariance)) * scale
    comparisons = []
    for i, (name, value, old, error) in enumerate(
        zip(
            np.array(parameter_names(suffix))[free],
            result.x,
            initial[free],
            np.array(errors, dtype=object)[free],
        )
    ):
        within = abs(value - old) <= 2 * np.hypot(error, standard_errors[i])
        comparisons.append(
            {
                "parameter": str(name),
                "published": float(old),
                "published_error": error,
                "refit": float(value),
                "refit_error": float(standard_errors[i]),
                "difference": float(value - old),
                "relative_difference": float(abs((value - old) / old)),
                "within_combined_2sigma": bool(within),
                "similar": bool(within),
            }
        )
    published_p = calculate(initial[free])
    return {
        "parameters": comparisons,
        "observations": len(v),
        "solver_success": bool(result.success),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_rmse_gpa": float(np.sqrt(np.mean((published_p - p) ** 2))),
        "published_max_abs_residual_gpa": float(max(abs(published_p - p))),
        "published_pressure_gpa": published_p.tolist(),
        "pressure_residuals_gpa": (published_p - p).tolist(),
        "measurement_numbers": [int(r["measurement_number"]) for r in data],
        "observed_pressure_range_gpa": [float(min(p)), float(max(p))],
        "observed_temperature_range_k": [float(min(t)), float(max(t))],
    }


@lru_cache(maxsize=1)
def reproduce():
    results = {}
    for suffix in CASES:
        result = fit(suffix)
        result["effective_variance_sensitivity"] = fit(suffix, weighted=True)[
            "parameters"
        ]
        if suffix.startswith("mgd"):
            result["literal_equation9_control"] = fit(suffix, literal_sign=True)
            _, (v, t, p) = observations(suffix)
            result["wrong_atom_count_controls_rmse_gpa"] = {
                str(n): float(
                    np.sqrt(np.mean((pressure(suffix, v, t, atoms=2 * n) - p) ** 2))
                )
                for n in [26, 31]
            }
        results[suffix] = result
    data = [r for r in rows() if r["au_volume_a3"]]
    va, t, p = (
        np.array([float(r[k]) for r in data])
        for k in ["au_volume_a3", "temperature_k", "pressure_a89_gpa"]
    )
    av0 = 4.0786**3
    au_p = bm3(va, av0, 166.65, 5.4823) + (0.00714 - 0.0115 * np.log(av0 / va)) * (
        t - 300
    )
    results["gold_anderson_diagnostic"] = {
        "V0_a3": av0,
        "related_record": "gold_anderson_1989_bm3_1",
        "pressure_residuals_gpa": (au_p - p).tolist(),
        "rmse_gpa": float(np.sqrt(np.mean((au_p - p) ** 2))),
        "max_abs_residual_gpa": float(max(abs(au_p - p))),
    }
    vm = ((2 / 5400.0**3 + 1 / 9230.0**3) / 3) ** (-1 / 3)
    results["elastic_debye_temperature_k"] = float(
        Planck / Boltzmann * (3 * 70 / (4 * np.pi * 618.2e-30)) ** (1 / 3) * vm
    )
    # Figure 10 / p.151: independently printed density-difference bounds at
    # 35 GPa. Use the source's common rho0=3.30 normalization and mean cell.
    density_checks = []
    for temperature in [300.0, 1273.0]:
        densities = {
            suffix: 3.30
            * 623.34
            / brentq(
                lambda volume: pressure(suffix, volume, temperature) - 35.0,
                450.0,
                750.0,
            )
            for suffix in [
                "htbm_a89",
                "htbm_t03",
                "thermal_pressure_a89",
                "mgd_elastic_debye_a89",
            ]
        }
        models = [
            value for name, value in densities.items() if not name.endswith("t03")
        ]
        density_checks.append(
            {
                "temperature_k": temperature,
                "pressure_gpa": 35.0,
                "density_g_cm3": densities,
                "scale_difference_g_cm3": abs(
                    densities["htbm_a89"] - densities["htbm_t03"]
                ),
                "model_spread_g_cm3": max(models) - min(models),
            }
        )
    results["figure10_density_checks"] = density_checks
    results["ambient_mean_volume_a3"] = float(
        np.mean(
            [float(r["volume_a3"]) for r in rows() if int(r["measurement_number"]) < 0]
        )
    )
    return results


def ledger_outcome(record):
    suffix = record["identifier"].removeprefix(PREFIX)
    result = reproduce()[suffix]
    names = [p["parameter"] for p in result["parameters"]]
    return {
        "status": "similar"
        if all(p["similar"] for p in result["parameters"])
        else "parity_not_achieved",
        "reason": "Published coefficients compared with independent Table 2 unweighted pressure least squares; all reported-error comparisons use combined two-error diagnostic intervals (source confidence unspecified). Source weights/covariance unavailable. No diagnostic fit replaces a published record. Explicit BM3 and MGD sign interpretations, weighted/literal controls and held alternatives documented in literature-reproductions/litasov-2007-superhydrous-phase-b.md.",
        "dataset_identifiers": [DATASET],
        "observations": result["observations"],
        "fit_kind": "independent_nonlinear_pressure_least_squares",
        "objective": "unweighted pressure residuals; all 69 rows for thermal models, 20 at 300 K for BM3",
        "free_parameters": names,
        "parameters": result["parameters"],
        **{
            key: result[key]
            for key in [
                "rmse_gpa",
                "published_rmse_gpa",
                "observed_pressure_range_gpa",
                "observed_temperature_range_k",
                "solver_success",
            ]
        },
    }


if __name__ == "__main__":
    OUTPUT.write_text(json.dumps(reproduce(), indent=2, allow_nan=False) + "\n")
