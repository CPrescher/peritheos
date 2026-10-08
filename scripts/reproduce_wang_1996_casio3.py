#!/usr/bin/env python3
"""Audit Wang et al. (1996) compression and thermal fits, without replacing EOS.

The rounded experimental table is checksummed. Source equations are fitted to
pressure residuals; alternative masks and weights are explicit diagnostics.
The printed reciprocal-sum weight has a singularity at a zero-stress row, so
volume-only and quadrature weights must never be labelled original weights.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/wang-1996-casio3-refit.json"
DATASET_ID = "ca_perovskite_wang_1996_table1_pvt"


def load_table(identifier):
    document = get_material_document("ca_perovskite")
    dataset = next(d for d in document["datasets"] if d["identifier"] == identifier)
    path = ROOT / "peritheos/data" / dataset["resource"]["path"]
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != dataset["resource"]["sha256"]:
        raise ValueError(f"source checksum mismatch: {identifier}")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return dataset, rows


def solve(model, initial, observed, *, weight=None, names=None):
    weight = np.ones_like(observed) if weight is None else np.asarray(weight)
    result = least_squares(
        lambda x: (model(x) - observed) * weight,
        initial,
        x_scale="jac",
        max_nfev=10000,
        ftol=1e-11,
        xtol=1e-11,
        gtol=1e-11,
    )
    if not result.success:
        raise RuntimeError(result.message)
    dof = len(observed) - len(initial)
    covariance = np.linalg.inv(result.jac.T @ result.jac) * np.sum(result.fun**2) / dof
    names = names or [f"p{i}" for i in range(len(initial))]
    return {
        "parameters": dict(zip(names, map(float, result.x))),
        "standard_errors": dict(zip(names, map(float, np.sqrt(np.diag(covariance))))),
        "covariance": covariance.tolist(),
        "covariance_convention": "local Jacobian, scaled by residual variance; not source covariance",
        "pressure_rmse_gpa": float(np.sqrt(np.mean((model(result.x) - observed) ** 2))),
        "weighted_residual_sum_squares": float(np.sum(result.fun**2)),
        "degrees_of_freedom": dof,
        "solver_success": bool(result.success),
    }


def thermal_fit(rows, *, selection="all_66", weighting="unit", log_coefficient=0.0):
    selected = [
        r
        for r in rows
        if selection == "all_66"
        or (selection == "run13_34" and r["run"] == "13")
        or (selection == "exclude_two_low_64" and float(r["pressure_gpa"]) >= 2)
    ]
    v, t, p, dv, stress = [
        np.array([float(r[k]) for r in selected])
        for k in (
            "volume_a3",
            "temperature_k",
            "pressure_gpa",
            "volume_reported_uncertainty_a3",
            "differential_stress_gpa",
        )
    ]
    sigma_volume = 232 * dv / v
    sigma_stress = 2 * np.abs(stress)
    if weighting == "unit":
        weight = np.ones_like(v)
    elif weighting == "volume_only_proxy":
        weight = 1 / sigma_volume
    elif weighting == "quadrature_proxy":
        weight = 1 / np.sqrt(sigma_volume**2 + sigma_stress**2)
    elif weighting == "literal_reciprocal_sum":
        if np.any(sigma_stress == 0):
            raise ValueError(
                "literal source weight is singular at zero differential stress"
            )
        weight = np.sqrt(1 / sigma_volume**2 + 1 / sigma_stress**2)
    else:
        raise ValueError(weighting)
    free_log = log_coefficient is None

    def model(x):
        v0, k0, slope = x[:3]
        log_slope = x[3] if free_log else log_coefficient
        return BM3(V0=v0, K0=k0, K0_prime=4.8).pressure(v) + (
            slope - log_slope * np.log(v0 / v)
        ) * (t - 300)

    names = ["V0_a3", "K0_gpa", "dP_dT_v_gpa_per_k"]
    if free_log:
        names += ["dK_dT_v_gpa_per_k"]
    result = solve(
        model,
        [45.58, 233, 0.0071] + ([0.0] if free_log else []),
        p,
        weight=weight,
        names=names,
    )
    k0 = result["parameters"]["K0_gpa"]
    slope = result["parameters"]["dP_dT_v_gpa_per_k"]
    log_slope = result["parameters"].get("dK_dT_v_gpa_per_k", log_coefficient)
    result.update(
        {
            "selection": selection,
            "source_orders": [int(r["source_order"]) for r in selected],
            "observations": len(selected),
            "weighting": weighting,
            "objective": "pressure residuals",
            "fixed_parameters": {
                "K0_prime": 4.8,
                "T0_k": 300,
                **({} if free_log else {"dK_dT_v_gpa_per_k": log_coefficient}),
            },
            "derived": {
                "alpha0_per_k": slope / k0,
                "dK_dT_p_gpa_per_k": log_slope - 4.8 * slope,
            },
        }
    )
    return result


def high_temperature_fit(rows, *, weighting="unit", linear_alpha=False):
    v, t, p, dv = [
        np.array([float(r[k]) for r in rows])
        for k in (
            "volume_a3",
            "temperature_k",
            "pressure_gpa",
            "volume_reported_uncertainty_a3",
        )
    ]
    weight = np.ones_like(v) if weighting == "unit" else v / (232 * dv)

    def model(x):
        alpha, derivative = x[:2]
        alpha_slope = x[2] if linear_alpha else 0
        dt = t - 300
        reference_volume = 45.58 * np.exp(alpha * dt + alpha_slope * dt**2 / 2)
        bulk_modulus = 232 + derivative * dt
        eta = reference_volume / v
        return (
            1.5
            * bulk_modulus
            * (eta ** (7 / 3) - eta ** (5 / 3))
            * (1 + 0.75 * 0.8 * (eta ** (2 / 3) - 1))
        )

    names = ["alpha0_per_k", "dK_dT_p_gpa_per_k"] + (
        ["dalpha_dT_p_per_k2"] if linear_alpha else []
    )
    result = solve(
        model,
        [3.1e-5, -0.025] + ([0.0] if linear_alpha else []),
        p,
        weight=weight,
        names=names,
    )
    result.update(
        {
            "observations": len(rows),
            "source_orders": [int(r["source_order"]) for r in rows],
            "selection": "all_66",
            "weighting": weighting,
            "objective": "pressure residuals",
            "fixed_parameters": {
                "V0_a3": 45.58,
                "K0_gpa": 232,
                "K0_prime": 4.8,
                "T0_k": 300,
            },
            "thermal_law": "V0(T)=V0 exp(alpha0 dT + b dT^2/2); K0(T)=K0+k dT; c=0",
        }
    )
    return result


def reproduce():
    dataset, rows = load_table(DATASET_ID)
    rt_dataset, rt_rows = load_table("ca_perovskite_wang_1996_table1_room_temperature")
    rt_rows = [r for r in rt_rows if r["fit_included"] == "1"]
    v, p = [
        np.array([float(r[k]) for r in rt_rows]) for k in ("volume_a3", "pressure_gpa")
    ]
    room_temperature = solve(
        lambda x: BM3(V0=x[0], K0=x[1], K0_prime=4.8).pressure(v),
        [45.58, 232],
        p,
        names=["V0_a3", "K0_gpa"],
    )
    room_temperature.update(
        {
            "observations": 12,
            "objective": "unit-weight pressure residuals",
            "fixed_parameters": {"K0_prime": 4.8},
            "status": "parameters_within_published_uncertainties",
        }
    )
    mao_dataset, mao_rows = load_table("ca_perovskite_mao_1989_table1_compression")
    mao_above_one = [r for r in mao_rows if float(r["pressure_gpa"]) >= 1]
    mv, mp, ms = [
        np.array([float(r[k]) for r in mao_above_one])
        for k in (
            "unit_cell_volume_a3",
            "pressure_gpa",
            "pressure_standard_deviation_gpa",
        )
    ]

    def mao_model(x):
        return BM3(V0=x[0], K0=x[1], K0_prime=x[2]).pressure(mv)

    mao_equal = solve(
        mao_model, [45.71, 244, 4.8], mp, names=["V0_a3", "K0_gpa", "K0_prime"]
    )
    mao_equal.update(
        {
            "observations": 44,
            "objective": "unit-weight pressure residuals",
            "selection": "exclude both sub-1-GPa rows",
            "status": "rounded_parameters_and_errors_reproduced",
        }
    )
    mao_weighted = solve(
        mao_model,
        [45.47, 268, 4.3],
        mp,
        weight=1 / ms,
        names=["V0_a3", "K0_gpa", "K0_prime"],
    )
    mao_weighted.update(
        {
            "observations": 44,
            "objective": "pressure residuals / printed pressure errors",
            "selection": "exclude both sub-1-GPa rows",
            "status": "not_reproduced_at_published_precision",
        }
    )
    thermal = {
        "unweighted_all_66": thermal_fit(rows),
        "unweighted_exclude_low_64": thermal_fit(rows, selection="exclude_two_low_64"),
        "weighted1_volume_proxy_free_log": thermal_fit(
            rows, weighting="volume_only_proxy", log_coefficient=None
        ),
        "weighted2_run13_volume_proxy": thermal_fit(
            rows, selection="run13_34", weighting="volume_only_proxy"
        ),
        "run13_literal_reciprocal_sum": thermal_fit(
            rows, selection="run13_34", weighting="literal_reciprocal_sum"
        ),
        "all_quadrature_proxy": thermal_fit(rows, weighting="quadrature_proxy"),
        "text_three_term_all_volume_proxy": thermal_fit(
            rows, weighting="volume_only_proxy"
        ),
        "text_fixed_log_all_volume_proxy": thermal_fit(
            rows, weighting="volume_only_proxy", log_coefficient=-0.0036
        ),
        "text_fixed_log_run13_volume_proxy": thermal_fit(
            rows,
            selection="run13_34",
            weighting="volume_only_proxy",
            log_coefficient=-0.0036,
        ),
    }
    high_temperature = {
        "weighted_volume_proxy": high_temperature_fit(
            rows, weighting="volume_only_proxy"
        ),
        "unweighted1_constant_alpha": high_temperature_fit(rows),
        "unweighted2_linear_alpha": high_temperature_fit(rows, linear_alpha=True),
    }
    return {
        "format": "peritheos.wang-1996-casio3-reproduction-audit",
        "audit_date": "2026-10-08",
        "source_doi": "10.1029/95JB03254",
        "source_pdf_sha256": dataset["provenance"]["source_pdf_sha256"],
        "dataset_checksums": {
            d["identifier"]: d["resource"]["sha256"]
            for d in (dataset, rt_dataset, mao_dataset)
        },
        "published_eos_replaced": False,
        "overall_status": "partial_reproduction_with_explicit_diagnostics",
        "published_table2_targets": {
            "thermal_weighted1": {
                "V0_a3": [45.60, 0.05],
                "K0_gpa": [233, 7],
                "dP_dT_v_gpa_per_k": [0.0075, 0.0005],
                "dK_dT_v_gpa_per_k": [0.005, 0.022],
            },
            "thermal_weighted2": {
                "V0_a3": [45.61, 0.02],
                "K0_gpa": [229, 4],
                "dP_dT_v_gpa_per_k": [0.0069, 0.0001],
            },
            "thermal_unweighted": {
                "V0_a3": [45.58, 0.04],
                "K0_gpa": [233, 9],
                "dP_dT_v_gpa_per_k": [0.0071, 0.0001],
            },
            "high_temperature_weighted": {
                "alpha0_per_k": [3.55e-5, 0.18e-5],
                "dK_dT_p_gpa_per_k": [-0.036, 0.008],
            },
            "high_temperature_unweighted1": {
                "alpha0_per_k": [3.06e-5, 0.14e-5],
                "dK_dT_p_gpa_per_k": [-0.022, 0.010],
            },
            "high_temperature_unweighted2": {
                "alpha0_per_k": [3.19e-5, 0.48e-5],
                "dK_dT_p_gpa_per_k": [-0.035, 0.015],
                "dalpha_dT_p_per_k2": [4.5e-9, 4e-9],
            },
        },
        "room_temperature": room_temperature,
        "mao_equal_weight_reanalysis": mao_equal,
        "mao_pressure_weighted_excluding_sub1": mao_weighted,
        "mao_pressure_weighted_all": {
            "status": "not_exactly_reconstructible",
            "observations": 46,
            "finding": "The 0-GPa observation has no printed pressure error; its source weight is unknown. No error is invented.",
        },
        "thermal_pressure": thermal,
        "high_temperature_bm3": high_temperature,
        "unresolved": {
            "literal_source_weight": "w=1/(2|differential stress|)^2 + 1/(232*volume_error/volume)^2 is singular at Table 1 source_order=35, run 3, printed stress=0.000 GPa. No floor is silently introduced.",
            "fit_masks": "The two sub-2-GPa rows are explicitly excluded in the compression discussion; thermal figure captions say all Table 1 data are plotted. Thermal regression row masks are not explicitly given. Both 66- and 64-row unit-weight diagnostics are retained.",
            "isobaric_isochoric": "The numerical interpolated coordinates, interpolation algorithm, weights and covariance used for Figures 5-6 are not published. Their original interpolation regressions cannot be uniquely recovered from Table 1 alone.",
            "text_table_discrepancy": "The prose three-term thermal fit reports K0=231(7), V0=45.61(4), alpha0=3.33(6)e-5; Table 2 weighted-1 instead reports K0=233(7), V0=45.60(5), alpha0=3.23e-5 and a free logarithmic term. These are distinct source-reported fits, not interchangeable targets.",
            "confidence": "Agreement within source error bars or under proxy weights is a qualified numerical check, not reproduction of original weights, covariance or measurement accuracy.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = reproduce()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote Wang (1996) audit to {args.output}")


if __name__ == "__main__":
    main()
