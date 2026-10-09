"""Reconstruct Yokoo Pt pressure from Table V outputs, preserving provenance."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict, replace
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from scripts.fit_yokoo_2009_gold_thermal import (
    Parameters,
    check_reconstruction,
    metrics,
    quadrature,
)
from scripts.reproduce_yokoo_2009_gold import source_bm3
from scripts.reproduce_yokoo_2009_platinum import DATA, SOURCE, TABLE, V0, table_rows

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/yokoo-2009-platinum-thermal-reconstruction.json"
RESIDUALS = ROOT / "docs/data/yokoo-2009-platinum-thermal-residuals.csv"
ELECTRONIC = "tsuchiya-kawamura-2002-table1-electronic-pressure.csv"
PUBLISHED = Parameters(0.99382001, 288.4, 5.05, 2.63, 0.39, 5.2, 230)
FREE = ("vc_over_v0", "k0_prime", "gamma0", "a", "b")
BOUNDS = {
    "vc_over_v0": (0.95, 1),
    "k0_gpa": (250, 330),
    "k0_prime": (3, 8),
    "gamma0": (1, 4),
    "a": (0, 1),
    "b": (0.1, 12),
    "theta0_k": (100, 400),
}


@lru_cache(maxsize=1)
def electronic_nodes():
    with (DATA / ELECTRONIC).open() as stream:
        rows = list(csv.DictReader(stream))
    return (
        np.array([float(r["temperature_k"]) for r in rows]),
        np.array([float(r["platinum_electronic_pressure_gpa"]) for r in rows]),
    )


def pressure_gpa(
    ratio, temperature, parameters, *, ambient_cell_a3=V0, order=64, correction=None
):
    """Independent SI Gauss-Legendre pressure evaluator, absolute 0 K reference."""
    ratio, temperature = np.broadcast_arrays(
        np.asarray(ratio, dtype=float), np.asarray(temperature, dtype=float)
    )
    if (
        np.any(~np.isfinite(ratio))
        or np.any(~np.isfinite(temperature))
        or np.any((ratio < 0.6) | (ratio > 1))
        or np.any((temperature < 0) | (temperature > 3000))
    ):
        raise ValueError("Reconstruction requires 0.6<=V/V0<=1 and 0<=T<=3000 K.")
    p = parameters
    gamma = p.gamma0 * (1 + p.a * (ratio**p.b - 1))
    theta = (
        p.theta0_k
        * ratio ** (-(1 - p.a) * p.gamma0)
        * np.exp(-(gamma - p.gamma0) / p.b)
    )
    phonon = np.zeros(ratio.shape)
    positive = temperature > 0
    x = theta[positive] / temperature[positive]
    integral = np.full(x.shape, np.pi**4 / 15)
    moderate = x <= 80
    z, w = quadrature(order)
    argument = x[moderate, None] * z
    integral[moderate] = x[moderate] * np.sum(
        w * argument**3 / np.expm1(argument), axis=1
    )
    energy = 9 * R * temperature[positive] * integral / x**3
    vm = ambient_cell_a3 * ratio[positive] * N_A * 1e-30 / 4
    phonon[positive] = gamma[positive] * energy / vm / 1e9
    nodes, electronic = electronic_nodes()
    result = (
        source_bm3(ratio, p.vc_over_v0, p.k0_gpa, p.k0_prime)
        + phonon
        + np.interp(temperature, nodes, electronic)
    )
    if correction is not None:
        result = result + np.interp(temperature, *correction)
    return result


def independent_quad_pressure(ratio, temperature, p, *, correction=None):
    gamma = p.gamma0 * (1 + p.a * (ratio**p.b - 1))
    theta = (
        p.theta0_k
        * ratio ** (-(1 - p.a) * p.gamma0)
        * np.exp(-(gamma - p.gamma0) / p.b)
    )
    phonon = 0.0
    if temperature > 0:
        x = theta / temperature
        integral = quad(lambda z: z**3 / np.expm1(z), 0, x, epsabs=1e-11)[0]
        phonon = gamma * 9 * R * temperature * integral / x**3
        phonon /= V0 * ratio * N_A * 1e-30 / 4 * 1e9
    nodes, electronic = electronic_nodes()
    result = float(
        source_bm3(ratio, p.vc_over_v0, p.k0_gpa, p.k0_prime)
        + phonon
        + np.interp(temperature, nodes, electronic)
    )
    return (
        result
        if correction is None
        else result + float(np.interp(temperature, *correction))
    )


def fit_with_pressure_correction(
    ratio, temperature, pressure, selection, *, initial=PUBLISHED, ambient_cell_a3=V0
):
    """Cold-only Vc fit, then a/b and profiled pressure corrections at source T.

    The per-temperature correction is the mean residual over selected volumes.
    It is an empirical pressure reconstruction term, not recovered Pel or Eel.
    Its 0 K value is fixed to zero. Both the cold fit and the warm residual
    averages exclude all withheld volumes and source first-liquid markers.
    """
    cold = selection & (temperature == 0)
    cold_fit = least_squares(
        lambda x: source_bm3(ratio[cold], x[0], 288.4, 5.05) - pressure[cold],
        [PUBLISHED.vc_over_v0],
        bounds=([0.95], [1.0]),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-12,
    )
    nodes = np.unique(temperature)

    def evaluate(x):
        p = replace(PUBLISHED, vc_over_v0=float(cold_fit.x[0]), a=x[0], b=x[1])
        delta = (
            pressure_gpa(ratio, temperature, p, ambient_cell_a3=ambient_cell_a3)
            - pressure
        )
        correction = np.array(
            [
                0.0 if t == 0 else -np.mean(delta[(temperature == t) & selection])
                for t in nodes
            ]
        )
        return delta + np.interp(temperature, nodes, correction), p, correction

    result = least_squares(
        lambda x: evaluate(x)[0][selection],
        [initial.a, initial.b],
        bounds=([0, 0.1], [1, 12]),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-12,
        max_nfev=2000,
    )
    if not result.success or not cold_fit.success:
        raise RuntimeError("Pressure reconstruction did not converge")
    _, p, correction = evaluate(result.x)
    return p, (nodes, correction)


def fit(
    ratio,
    temperature,
    pressure,
    selection,
    free=FREE,
    initial=PUBLISHED,
    ambient_cell_a3=V0,
):
    values = asdict(initial)

    def decode(x):
        return Parameters(**{**values, **dict(zip(free, x))})

    def residual(x):
        return (
            pressure_gpa(ratio, temperature, decode(x), ambient_cell_a3=ambient_cell_a3)
            - pressure
        )[selection]

    result = least_squares(
        residual,
        [values[k] for k in free],
        bounds=([BOUNDS[k][0] for k in free], [BOUNDS[k][1] for k in free]),
        x_scale="jac",
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-12,
        max_nfev=2000,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return decode(result.x)


@lru_cache(maxsize=1)
def reconstruct():
    rows = table_rows()
    ratio, temperature, pressure = [
        np.array([float(row[key]) for row in rows])
        for key in ("volume_ratio", "temperature_k", "pressure_gpa")
    ]
    unmarked = np.array([row["source_phase_annotation"] == "unmarked" for row in rows])
    fits = {}
    for name, free, selected in [
        ("printed_coefficients_cold_volume_only", ("vc_over_v0",), unmarked),
        ("five_parameter_without_pressure_correction", FREE, unmarked),
        ("include_first_liquid_states", FREE, np.ones(len(rows), dtype=bool)),
        ("also_refine_cold_bulk_modulus", (*FREE, "k0_gpa"), unmarked),
        (
            "also_refine_bulk_modulus_and_debye_temperature",
            tuple(asdict(PUBLISHED)),
            unmarked,
        ),
    ]:
        p = fit(ratio, temperature, pressure, selected, free)
        delta = pressure_gpa(ratio, temperature, p) - pressure
        fits[name] = {
            "parameters": asdict(p),
            "refined_parameters": list(free),
            "fixed_parameters": [k for k in asdict(PUBLISHED) if k not in free],
            "fit_states": metrics(delta[selected]),
            "unmarked_states": metrics(delta[unmarked]),
            "first_liquid_states": metrics(delta[~unmarked]),
            "all_populated_states": metrics(delta),
        }
    p, correction = fit_with_pressure_correction(ratio, temperature, pressure, unmarked)
    delta = pressure_gpa(ratio, temperature, p, correction=correction) - pressure
    fits["primary"] = {
        "parameters": asdict(p),
        "refined_parameters": ["vc_over_v0", "a", "b"],
        "fixed_parameters": ["k0_gpa", "k0_prime", "gamma0", "theta0_k"],
        "residual_temperature_k": correction[0].tolist(),
        "residual_pressure_gpa": correction[1].tolist(),
        "fit_states": metrics(delta[unmarked]),
        "unmarked_states": metrics(delta[unmarked]),
        "first_liquid_states": metrics(delta[~unmarked]),
        "all_populated_states": metrics(delta),
    }
    distinct_starts = [
        PUBLISHED,
        replace(PUBLISHED, vc_over_v0=0.98, k0_prime=4.2, gamma0=2.1, a=0.2, b=2),
        replace(PUBLISHED, vc_over_v0=0.998, k0_prime=6, gamma0=3.3, a=0.6, b=8),
    ]
    starts = [
        {
            "initial": asdict(start),
            "fitted": asdict(
                fit_with_pressure_correction(
                    ratio, temperature, pressure, unmarked, initial=start
                )[0]
            ),
        }
        for start in distinct_starts
    ]
    ratios = np.sort(np.unique(ratio))
    held_out = unmarked & np.isin(ratio, ratios[1::2])
    holdout_p, holdout_correction = fit_with_pressure_correction(
        ratio, temperature, pressure, unmarked & ~held_out
    )
    holdout_delta = (
        pressure_gpa(ratio, temperature, holdout_p, correction=holdout_correction)
        - pressure
    )
    density_v0 = json.loads((DATA / SOURCE).read_text(encoding="utf-8"))[
        "volume_basis"
    ]["V0_cell_a3_derived_from_printed_density"]
    density_p, density_correction = fit_with_pressure_correction(
        ratio, temperature, pressure, unmarked, ambient_cell_a3=density_v0
    )
    return {
        "doi": "10.1103/PhysRevB.80.104114",
        "kind": "derived_table_output_pressure_reconstruction",
        "diagnostic_only": True,
        "published_analytical_pvt_reproduced": False,
        "author_fit_reproduced": False,
        "experimental_observations": 0,
        "qualification": "Fit to derived Table V outputs, not original experimental observations. Cold Vc, refined a/b and the separate empirical residual pressure table retain derived provenance; published coefficients and raw electronic nodes are unchanged. The residual correction is not an electronic-pressure observation. No caloric data or original fitting weights are needed for this bounded PVT reconstruction.",
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [TABLE, SOURCE, ELECTRONIC]
        },
        "conventions": {
            "ambient_cell_a3": V0,
            "atoms_per_cell": 4,
            "thermal_volume_normalization": "Ambient 300 K V0, distinct from cold Vc",
            "phonon_reference": "Absolute; no zero-point or 300 K subtraction",
            "electronic": "Raw 0 K-referenced Pt pressure nodes, volume independent, linear interpolation, no extrapolation",
            "electronic_reference_subtracted_at_300k": False,
            "electronic_interpolation_author_verified": False,
            "objective": "Fit Vc to selected cold rows with K0/K0-prime fixed; then a/b to selected rows with unit pressure weights, profiling the mean residual at each source temperature; residual at 0 K fixed to zero",
            "residual_pressure": "Separate empirical reconstruction correction, linear between source temperatures, not assigned to electronic or phonon physics",
            "optimizer": "scipy least_squares TRF, x_scale=jac, ftol=xtol=1e-13, gtol=1e-12",
            "uncertainty": "No experimental errors or covariance inferred from rounded model-output residuals",
            "bounds": {key: list(value) for key, value in BOUNDS.items()},
            "domain": {
                "volume_ratio": [0.6, 1],
                "temperature_k": [0, 3000],
                "phase_stability_verified": False,
            },
        },
        "parameter_comparison": [
            {
                "parameter": k,
                "published": None if k == "vc_over_v0" else getattr(PUBLISHED, k),
                "reconstructed": v,
                "relative_difference_percent": None
                if k == "vc_over_v0"
                else 100 * (v / getattr(PUBLISHED, k) - 1),
            }
            for k, v in asdict(p).items()
        ],
        "fits": fits,
        "multistart": starts,
        "volume_interpolation_holdout": {
            "held_out_volume_ratios": ratios[1::2].tolist(),
            "parameters": asdict(holdout_p),
            "residual_pressure_gpa": holdout_correction[1].tolist(),
            "training_states": metrics(holdout_delta[unmarked & ~held_out]),
            "held_out_states": metrics(holdout_delta[held_out]),
        },
        "unrounded_density_conversion": {
            "ambient_cell_a3": density_v0,
            "parameters": asdict(density_p),
            "residual_pressure_gpa": density_correction[1].tolist(),
            "unmarked_states": metrics(
                (
                    pressure_gpa(
                        ratio,
                        temperature,
                        density_p,
                        ambient_cell_a3=density_v0,
                        correction=density_correction,
                    )
                    - pressure
                )[unmarked]
            ),
        },
        "adaptive_quad_max_difference_gpa": max(
            abs(
                float(pressure_gpa(r, t, p, correction=correction))
                - independent_quad_pressure(r, t, p, correction=correction)
            )
            for r, t in zip(ratio, temperature)
        ),
        "cold_volume_cell_a3": p.vc_over_v0 * V0,
        "states": [
            {
                **row,
                "model_pressure_gpa": float(m),
                "model_minus_table_gpa": float(d),
                "used_in_primary_fit": bool(used),
            }
            for row, m, d, used in zip(rows, pressure + delta, delta, unmarked)
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reconstruct()
    if args.check:
        check_reconstruction(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
    else:
        OUTPUT.write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        with RESIDUALS.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(report["states"][0]))
            writer.writeheader()
            writer.writerows(report["states"])
    print(
        json.dumps(
            {
                "parameters": report["fits"]["primary"]["parameters"],
                "fit": report["fits"]["primary"]["fit_states"],
                "holdout": report["volume_interpolation_holdout"]["held_out_states"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
