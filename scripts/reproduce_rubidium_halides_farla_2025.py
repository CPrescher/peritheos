#!/usr/bin/env python3
"""Audit Farla 2025 preferred curves and conditional independent fits.

The reference calculation uses explicit BM3 and Gauss-Legendre integration of
Eqs. (6)-(9), without Peritheos evaluators. The independent fit is a diagnostic
effective-variance least-squares approximation, not an EoSFit7 workflow replay.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R, h, k
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
PREFIX = "rubidium-halides-farla-2025"
# Table I, case 1. Molar V0 [cm3/mol], theta0 [K], fitted [K0, K0', gamma0, q].
SOURCE = {
    "rbcl": (36.25, 177.0, [19.89, 5.00, 1.96, 1.05]),
    "rbbr": (41.74, 145.0, [16.28, 5.28, 2.18, 1.52]),
    "rbi": (51.11, 118.0, [13.69, 4.95, 2.21, 1.42]),
}
ERRORS = {
    "rbcl": [0.08, 0.02, 0.04, 0.09],
    "rbbr": [0.04, 0.02, 0.14, 0.24],
    "rbi": [0.04, 0.01, 0.07, 0.10],
}
NODES, WEIGHTS = np.polynomial.legendre.leggauss(48)


def debye_anchor_check(compound):
    """Check adopted theta using supplemental K/G/rho/M and physical SI laws.

    The supplement's printed negative exponents on sound speed/number density
    are dimensionally inconsistent. Its reported scalar speeds and theta agree
    with positive square/cube roots; this check does not alter the adopted theta.
    """
    bulk, shear, density, mass, source_speed = {
        "rbcl": (18.35, 9.18, 3.337, 120.92, 1849.6),
        "rbbr": (15.86, 8.04, 3.966, 165.37, 1586.9),
        "rbi": (12.57, 6.39, 4.157, 212.37, 1382.0),
    }[compound]
    vs = np.sqrt(shear * 1e9 / (density * 1000))
    vl = np.sqrt((bulk + 4 * shear / 3) * 1e9 / (density * 1000))
    vm = ((2 / vs**3 + 1 / vl**3) / 3) ** (-1 / 3)
    theta = (
        h
        / k
        * (3 * 2 * Avogadro * density * 1000 / (4 * np.pi * mass / 1000)) ** (1 / 3)
        * vm
    )
    return {
        "bulk_gpa": bulk,
        "shear_gpa": shear,
        "density_g_cm3": density,
        "molar_mass_g_mol": mass,
        "source_average_speed_m_s": source_speed,
        "calculated_average_speed_m_s": float(vm),
        "calculated_theta_k": float(theta),
        "source_fixed_theta_k": SOURCE[compound][1],
        "source": "Official supplement section 1, DFT-derived K and G",
    }


def read_rows(compound: str, kind: str) -> list[dict[str, str]]:
    with (DATA / f"{PREFIX}-{compound}-{kind}.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))


def reference_pressure(volume, temperature, molar_v0, theta0, parameters):
    """Independent pressure [GPa] for conventional Z=1 volume [A3]."""
    k0, kp, gamma0, q = np.asarray(parameters)
    v, t = np.broadcast_arrays(np.asarray(volume), np.asarray(temperature))
    v0 = molar_v0 * 1e24 / Avogadro
    ratio = v / v0
    eta = ratio ** (-1 / 3)
    cold = 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))
    gamma = gamma0 * ratio**q
    theta = theta0 * np.exp((gamma0 - gamma) / q)

    def energy(temp):
        upper = theta / temp
        xi = upper[..., None] * (NODES + 1) / 2
        integral = upper / 2 * np.sum(WEIGHTS * xi**3 / np.expm1(xi), axis=-1)
        return 9 * 2 * R * temp * integral / upper**3

    # Molar energy / molar cm3 becomes GPa by multiplying 1e-3.
    return (
        cold
        + gamma
        * (energy(t) - energy(np.full_like(t, 300.0)))
        / (v * Avogadro * 1e-24)
        * 1e-3
    )


def reference_ratio(pressure, temperature, molar_v0, theta0, parameters):
    v0 = molar_v0 * 1e24 / Avogadro
    return brentq(
        lambda ratio: (
            reference_pressure(ratio * v0, temperature, molar_v0, theta0, parameters)
            - pressure
        ),
        0.2,
        1.12,
        xtol=2e-14,
    )


def benchmark_diagnostics(compound):
    """Disclose last-digit differences; use a stated coefficient-rounding envelope."""
    from peritheos import get_material

    molar, theta, parameters = SOURCE[compound]
    parameters = np.array(parameters)
    record = get_material(compound).get_eos_record(
        f"{compound}_b2_farla_2025_bm3_mgd_cscl"
    )
    rows = read_rows(compound, "benchmarks")
    results = []
    for row in rows:
        p, t, expected = (
            float(row[key]) for key in ("pressure_gpa", "temperature_k", "volume_ratio")
        )
        calculated = reference_ratio(p, t, molar, theta, parameters)
        native = record.volume(p, t) / record.reference_volume
        # Linear sum of sensitivities to each half-last-digit coefficient
        # perturbation; not a confidence interval and not a fit tolerance.
        envelope = 0.00005  # four-decimal ratio table rounding
        for i in range(4):
            perturbation = np.eye(4)[i] * 0.005
            plus = reference_ratio(p, t, molar, theta, parameters + perturbation)
            minus = reference_ratio(p, t, molar, theta, parameters - perturbation)
            envelope += max(abs(plus - calculated), abs(minus - calculated))
        # theta is printed as an integer; molar V0 is printed to 0.01 cm3/mol.
        for dv, dt in ((0.005, 0.0), (0.0, 0.5)):
            plus = reference_ratio(p, t, molar + dv, theta + dt, parameters)
            minus = reference_ratio(p, t, molar - dv, theta - dt, parameters)
            envelope += max(abs(plus - calculated), abs(minus - calculated))
        results.append(
            {
                "pressure_gpa": p,
                "temperature_k": t,
                "source_ratio": expected,
                "reference_ratio": calculated,
                "difference": calculated - expected,
                "rounding_sensitivity_envelope": envelope,
                "within_rounding_envelope": abs(calculated - expected) <= envelope,
                "native_reference_difference": native - calculated,
                "molten_source_flag": bool(int(row["molten_source_flag"])),
            }
        )
    return {
        "rows": len(results),
        "max_absolute_ratio_difference": max(abs(r["difference"]) for r in results),
        "max_native_reference_ratio_difference": max(
            abs(r["native_reference_difference"]) for r in results
        ),
        "outside_last_digit_rounding": sum(
            abs(r["difference"]) > 5e-5 for r in results
        ),
        "outside_coefficient_rounding_envelope": sum(
            not r["within_rounding_envelope"] for r in results
        ),
        "states": results,
    }


def independent_fit(compound, *, include_cold_lvp=False, weighted=True):
    """Fit all DAC rows + available hot LVP rows, with V0/theta0 fixed.

    The effective variance is sP²+(dP/dV*sV)²+(dP/dT*sT)².
    It assumes independent errors because author correlation matrices are absent.
    A diagnostic 0.005 GPa floor is used for S1 printed-zero pressure errors.
    Original dataset zeros are unchanged; this floor is not a recovered error.
    """
    molar, theta, parameters = SOURCE[compound]
    rows = [
        r
        for r in read_rows(compound, "observations")
        if r["volume_a3"]
        and (
            include_cold_lvp
            or r["source_table"] == "S1"
            or float(r["temperature_k"]) > 300
        )
    ]
    v, t, p, sv, st, sp = np.array(
        [
            [float(r[key] or 0) for r in rows]
            for key in (
                "volume_a3",
                "temperature_k",
                "pressure_gpa",
                "volume_uncertainty_a3",
                "temperature_uncertainty_k",
                "pressure_uncertainty_gpa",
            )
        ]
    )
    sp = np.maximum(sp, 0.005)

    def residual(coefs):
        calculated = reference_pressure(v, t, molar, theta, coefs)
        if not weighted:
            return calculated - p
        dpdv = (
            reference_pressure(v + 1e-4, t, molar, theta, coefs)
            - reference_pressure(v - 1e-4, t, molar, theta, coefs)
        ) / 2e-4
        dpdt = (
            reference_pressure(v, t + 0.01, molar, theta, coefs)
            - reference_pressure(v, t - 0.01, molar, theta, coefs)
        ) / 0.02
        variance = sp**2 + (dpdv * sv) ** 2 + (dpdt * st) ** 2
        return (calculated - p) / np.sqrt(variance)

    fit = least_squares(
        residual,
        parameters,
        bounds=([1, 2, 0.1, 0.05], [50, 9, 5, 5]),
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-9,
    )
    published_residual = reference_pressure(v, t, molar, theta, parameters) - p
    fitted_residual = reference_pressure(v, t, molar, theta, fit.x) - p
    return {
        "rows": len(rows),
        "include_cold_lvp": include_cold_lvp,
        "objective": "independent_PVT_effective_variance"
        if weighted
        else "unweighted_P",
        "fixed": {"molar_V0_cm3_mol": molar, "theta0_k": theta, "n": 2, "Tr_k": 300},
        "source_parameters": parameters,
        "fit_parameters": fit.x.tolist(),
        "fit_minus_source": (fit.x - parameters).tolist(),
        "fit_minus_source_over_reported_esd": (
            (fit.x - parameters) / ERRORS[compound]
        ).tolist(),
        "source_pressure_rmse_gpa": float(np.sqrt(np.mean(published_residual**2))),
        "source_pressure_max_abs_gpa": float(np.max(abs(published_residual))),
        "fitted_pressure_rmse_gpa": float(np.sqrt(np.mean(fitted_residual**2))),
        "success": bool(fit.success),
        "author_fit_parity": False,
        "author_workflow_parity": False,
    }


def dac_refit(compound):
    """Independent unweighted BM3 fit of S1; not a source parameter replacement."""
    molar, theta, parameters = SOURCE[compound]
    rows = [r for r in read_rows(compound, "observations") if r["source_table"] == "S1"]
    v = np.array([float(r["volume_a3"]) for r in rows])
    p = np.array([float(r["pressure_gpa"]) for r in rows])

    def residual(coefs):
        return (
            reference_pressure(
                v,
                np.full_like(v, 300),
                coefs[0],
                theta,
                [coefs[1], coefs[2], *parameters[2:]],
            )
            - p
        )

    fit = least_squares(
        residual,
        [molar, *parameters[:2]],
        bounds=([molar * 0.8, 1, 2], [molar * 1.2, 50, 9]),
    )
    return {
        "rows": len(rows),
        "objective": "unweighted_pressure_residual",
        "parameter_order": ["molar_V0_cm3_mol", "K0_gpa", "K0_prime"],
        "fit_parameters": fit.x.tolist(),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "success": bool(fit.success),
        "author_fit_parity": False,
    }


def reproduce():
    return {
        "source_doi": "10.1063/5.0248905",
        "audit_date": "2026-10-15",
        "parameter_order": ["K0_gpa", "K0_prime", "gamma0", "q"],
        "limitations": [
            "Printed parameters and calculated ratios have finite precision.",
            "Conditional independent fits assume diagonal PVT errors and a 0.005 GPa "
            "diagnostic floor for printed-zero DAC errors; not EoSFit7 workflow parity.",
            "Cold LVP rows are separately included as a sensitivity case; source row "
            "flags and original EoSFit7 project files are unavailable.",
            "Benchmarks beyond data coverage and molten flags do not extend validity.",
        ],
        "compounds": {
            compound: {
                "debye_anchor": debye_anchor_check(compound),
                "benchmarks": benchmark_diagnostics(compound),
                "dac_refit": dac_refit(compound),
                "independent_fits": [
                    independent_fit(compound),
                    independent_fit(compound, weighted=False),
                    independent_fit(compound, include_cold_lvp=True),
                ],
            }
            for compound in SOURCE
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = reproduce()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
