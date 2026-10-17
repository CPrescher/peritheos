#!/usr/bin/env python3
"""Follow up Farla 2025: historical calibration, weighting and parity limits.

No network access is needed: public author-code input facts are pinned in the
family's author-inputs JSON. No published observations or coefficients change.
Reconstructed input choices are not the authors' recovered EoSFit7 files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R
from scipy.optimize import brentq, least_squares, minimize

from scripts.reproduce_rubidium_halides_farla_2025 import (
    DATA,
    ERRORS,
    PREFIX,
    SOURCE,
    read_rows,
    reference_pressure,
    reference_ratio,
)

INPUTS = DATA / f"{PREFIX}-author-inputs.json"
AUTHOR_THETA = {"rbcl": 177.04, "rbbr": 144.96, "rbi": 117.98}
OLD_NA = 6.0221367e23


def cscl_pressure(
    volume_a3, temperature, *, gui_override=True, legacy_na=True, reference_a0=4.1221
):
    """Independent Vinet+MGD replay of EosCross v1.3 module/GUI defaults.

    Source: commit 9518e361, Dewaele_2020.py, MBase_CsCl+Pt.py and
    single/CsCl1_Dewaele2020.py. The GUI overrides the module's V0 from a0.
    This is an identified historical default, not a recovered saved session.
    """
    facts = json.loads(INPUTS.read_text(encoding="utf-8"))
    source = next(s for s in facts["sources"] if s["id"] == "cscl_module")
    params = source["parameter_dictionaries_si"]["_params_Dewaele_CsCl"]
    na = OLD_NA if legacy_na else Avogadro
    molar_v = np.asarray(volume_a3) * na * 1e-24
    molar_v0 = reference_a0**3 * na * 1e-24 if gui_override else params["V_0"] * 1e6
    x = (molar_v / molar_v0) ** (1 / 3)
    cold = (
        3
        * params["K_0"]
        * 1e-9
        * (1 - x)
        / x**2
        * np.exp(1.5 * (params["Kprime_0"] - 1) * (1 - x))
    )
    # n=2 and Tr=300 in both reference_pressure and this source calibrant.
    thermal = reference_pressure(
        molar_v * 1e24 / Avogadro,
        temperature,
        molar_v0,
        params["Debye_0"],
        [0, 4, params["grueneisen_0"], params["q_0"]],
    )
    return cold + thermal


def calibration_replay():
    """Replay 109 distinct gauge readings, repeated in S2/S3/S4."""
    all_rows = {
        c: [r for r in read_rows(c, "observations") if r["source_table"] != "S1"]
        for c in SOURCE
    }
    # These are the same gauge observations printed with each compound.
    fields = [
        "run",
        "temperature_k",
        "cscl_volume_a3",
        "cscl_volume_uncertainty_a3",
        "cscl_pressure_gpa",
        "cscl_pressure_uncertainty_gpa",
    ]
    for compound in ("rbbr", "rbi"):
        assert [[r[k] for k in fields] for r in all_rows[compound]] == [
            [r[k] for k in fields] for r in all_rows["rbcl"]
        ]
    rows = all_rows["rbcl"]
    t, v, sv, p, sp = [np.array([float(r[key]) for r in rows]) for key in fields[1:]]
    result = {}
    for name, gui, legacy in (
        ("module_si", False, False),
        ("module_legacy_na", False, True),
        ("gui_defaults_legacy_na", True, True),
    ):

        def evaluate(volume, temp):
            return cscl_pressure(volume, temp, gui_override=gui, legacy_na=legacy)

        calculated = evaluate(v, t)
        rounding = 0.005 + np.maximum(
            abs(evaluate(v + 0.0005, t) - calculated),
            abs(evaluate(v - 0.0005, t) - calculated),
        )
        error = abs(calculated - p)
        dpdv = (evaluate(v + 1e-4, t) - evaluate(v - 1e-4, t)) / 2e-4
        dpdt = (evaluate(v, t + 0.01) - evaluate(v, t - 0.01)) / 0.02
        propagated = np.sqrt((dpdv * sv) ** 2 + (dpdt * 30) ** 2)
        outliers = []
        for i in np.flatnonzero(error > rounding):
            implied_t = brentq(lambda temp: evaluate(v[i], temp) - p[i], 100, 2200)
            outliers.append(
                {
                    "observation_index": int(i + 1),
                    "source_row_s2": int(rows[i]["source_row"]),
                    "run": rows[i]["run"],
                    "temperature_k": float(t[i]),
                    "cscl_volume_a3": float(v[i]),
                    "source_pressure_gpa": float(p[i]),
                    "calculated_pressure_gpa": float(calculated[i]),
                    "calculated_minus_source_gpa": float(calculated[i] - p[i]),
                    "rounding_allowance_gpa": float(rounding[i]),
                    "reported_pressure_error_gpa": float(sp[i]),
                    "implied_temperature_shift_k": float(implied_t - t[i]),
                }
            )
        result[name] = {
            "unique_rows": len(rows),
            "repeated_source_rows": 3 * len(rows),
            "within_pressure_rounding": int(sum(error <= 0.005)),
            "within_pressure_and_volume_rounding": int(sum(error <= rounding)),
            "within_reported_pressure_error_plus_rounding": int(
                sum(error <= sp + rounding)
            ),
            "rmse_gpa": float(np.sqrt(np.mean(error**2))),
            "max_abs_difference_gpa": float(max(error)),
            "propagated_sigma_matches_printed_rounding": int(
                sum(abs(propagated - sp) <= 0.005)
            ),
            "outliers_beyond_rounding": outliers,
        }
    result["qualification"] = (
        "GUI defaults and module inputs are recovered from public historical code. "
        "A saved runtime session is absent. Implied temperature shifts diagnose a "
        "pressure discrepancy only; they are not evidence of actual temperatures "
        "or author-applied corrections. Source observations remain unchanged."
    )
    return result


def weighted_refit(
    compound,
    *,
    pressure_floor=0.0,
    include_cold_lvp=False,
    author_theta=False,
    eosfit_r=False,
    recalibrate=False,
):
    """Iterate frozen effective-variance weights as described by Angel 2014.

    Unlike minimizing (residual/sigma(parameters)) directly, each inner solve
    holds sigma fixed. Neither variant recreates unknown author row flags or
    unrounded errors. Zero printed sP is allowed because finite sV supplies
    nonzero total variance; it is not treated as exactly known pressure.
    """
    molar, theta, source = SOURCE[compound]
    if author_theta:
        theta = AUTHOR_THETA[compound]
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
    v, t, p, sv, st, sp = [
        np.array([float(r[key] or 0) for r in rows])
        for key in (
            "volume_a3",
            "temperature_k",
            "pressure_gpa",
            "volume_uncertainty_a3",
            "temperature_uncertainty_k",
            "pressure_uncertainty_gpa",
        )
    ]
    sp = np.maximum(sp, pressure_floor)
    if recalibrate:
        for i, row in enumerate(rows):
            if row["source_table"] != "S1":
                p[i] = cscl_pressure(float(row["cscl_volume_a3"]), t[i])

    def evaluate(volume, temperature, coefs):
        value = reference_pressure(volume, temperature, molar, theta, coefs)
        if not eosfit_r:
            return value
        cold = reference_pressure(
            volume, np.full_like(temperature, 300), molar, theta, coefs
        )
        # CrysFML 2021 EthDebye uses R=8.314; test that convention explicitly.
        return cold + (value - cold) * 8.314 / R

    coefs = np.array(source)
    for cycle in range(50):
        dpdv = (evaluate(v + 1e-4, t, coefs) - evaluate(v - 1e-4, t, coefs)) / 2e-4
        dpdt = (evaluate(v, t + 0.01, coefs) - evaluate(v, t - 0.01, coefs)) / 0.02
        sigma = np.sqrt(sp**2 + (dpdv * sv) ** 2 + (dpdt * st) ** 2)
        fitted = least_squares(
            lambda x: (evaluate(v, t, x) - p) / sigma,
            coefs,
            bounds=([1, 2, 0.1, 0.05], [50, 9, 5, 5]),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-10,
        )
        change = max(abs(fitted.x - coefs))
        coefs = fitted.x
        if change < 1e-8:
            break
    delta = coefs - source
    return {
        "rows": len(rows),
        "cycles": cycle + 1,
        "converged": bool(change < 1e-8),
        "theta0_k": theta,
        "R_j_mol_k": 8.314 if eosfit_r else R,
        "pressure_floor_gpa": pressure_floor,
        "include_cold_lvp": include_cold_lvp,
        "pressure_targets": "reconstructed_CsCl_GUI_defaults_for_LVP"
        if recalibrate
        else "reported",
        "fit_parameters": coefs.tolist(),
        "fit_minus_source": delta.tolist(),
        "fit_minus_source_over_reported_esd": (delta / ERRORS[compound]).tolist(),
        "all_coefficients_round_to_source": bool(np.all(abs(delta) <= 0.005)),
        "max_coefficient_difference_in_reported_esd": float(
            max(abs(delta / ERRORS[compound]))
        ),
        "pressure_rmse_gpa": float(np.sqrt(np.mean((evaluate(v, t, coefs) - p) ** 2))),
        "author_input_files_recovered": False,
        "author_workflow_parity": False,
    }


def joint_benchmark_rounding(compound):
    """Seek one parameter vector consistent with all 54 rounded ratios.

    A per-state sensitivity envelope does not demonstrate this stronger claim.
    Fits here reconstruct calculated tables only, not measured-data EOS fits.
    A failed numerical search is not a proof of mathematical infeasibility.
    """
    molar, _, parameters = SOURCE[compound]
    theta = AUTHOR_THETA[compound]
    rows = read_rows(compound, "benchmarks")
    ratio, t, p = [
        np.array([float(r[k]) for r in rows])
        for k in ("volume_ratio", "temperature_k", "pressure_gpa")
    ]
    v0 = molar * 1e24 / Avogadro
    parameters = np.array(parameters)

    def pressure(z, r=ratio):
        return reference_pressure(r * v0, t, molar, theta, parameters + z * 0.005)

    scale = (
        pressure(np.zeros(4), ratio - 5e-5) - pressure(np.zeros(4), ratio + 5e-5)
    ) / 2
    fit = least_squares(
        lambda z: (pressure(z) - p) / scale,
        np.zeros(4),
        bounds=(-1, 1),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-10,
    )

    def constraints(z):
        return np.r_[
            (pressure(z, ratio - 5e-5) - p) / scale,
            (p - pressure(z, ratio + 5e-5)) / scale,
        ]

    optimum = minimize(
        lambda z: float(np.sum(((pressure(z) - p) / scale) ** 2)),
        fit.x,
        method="SLSQP",
        bounds=[(-1, 1)] * 4,
        constraints=[{"type": "ineq", "fun": constraints}],
        options={"ftol": 1e-11, "maxiter": 300},
    )
    reconstructed = parameters + 0.005 * optimum.x
    calculated = np.array(
        [reference_ratio(pi, ti, molar, theta, reconstructed) for pi, ti in zip(p, t)]
    )
    difference = calculated - ratio
    feasible = bool(max(abs(difference)) <= 5e-5 + 1e-9)
    return {
        "rows": len(rows),
        "fixed_molar_v0_cm3_mol": molar,
        "fixed_author_theta0_k": theta,
        "coefficient_half_rounding_bounds": [0.005] * 4,
        "ratio_half_rounding_interval": 0.00005,
        "joint_rounded_table_solution_found": feasible,
        "optimizer_success": bool(optimum.success),
        "optimizer_message": optimum.message,
        "candidate_coefficients": reconstructed.tolist(),
        "max_absolute_ratio_difference": float(max(abs(difference))),
        "outside_ratio_rounding": int(sum(abs(difference) > 5e-5 + 1e-9)),
        "outside_states": [
            {
                "pressure_gpa": float(p[i]),
                "temperature_k": float(t[i]),
                "ratio_difference": float(difference[i]),
            }
            for i in np.flatnonzero(abs(difference) > 5e-5 + 1e-9)
        ],
        "qualification": (
            "Calculated-table reconstruction, not original coefficients or a fit to observations. "
            "No candidate values are written to catalog records. Success establishes a compatible "
            "vector, not its uniqueness or author provenance. Failure is an unresolved diagnostic."
        ),
    }


def reproduce():
    return {
        "audit_date": "2026-10-15",
        "source_doi": "10.1063/5.0248905",
        "parameter_order": ["K0_gpa", "K0_prime", "gamma0", "q"],
        "author_inputs": str(INPUTS.relative_to(DATA.parent)),
        "original_eosfit_files_recovered": False,
        "author_contact": "none",
        "calibration": calibration_replay(),
        "compounds": {
            c: {
                "iterated_weight_refits": {
                    "literal_printed_errors": weighted_refit(c),
                    "diagnostic_floor_0p005": weighted_refit(c, pressure_floor=0.005),
                    "author_theta_and_eosfit_R": weighted_refit(
                        c, author_theta=True, eosfit_r=True
                    ),
                    "include_cold_lvp": weighted_refit(
                        c, author_theta=True, eosfit_r=True, include_cold_lvp=True
                    ),
                    "reconstructed_cscl_targets": weighted_refit(
                        c, author_theta=True, eosfit_r=True, recalibrate=True
                    ),
                },
                "joint_benchmark_rounding": joint_benchmark_rounding(c),
            }
            for c in SOURCE
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    serialized = json.dumps(reproduce(), indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized)


if __name__ == "__main__":
    main()
