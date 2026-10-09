"""Reconstruct Au thermal pressure from Yokoo Table III, not experimental data.

Run with ``python -m scripts.fit_yokoo_2009_gold_thermal`` to archive the fit,
``--check`` to reproduce it, or ``--pressure V_OVER_V0 T_K`` to evaluate it.
The published catalog records are never modified by this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from scripts.reproduce_yokoo_2009_gold import (
    DATA,
    ELECTRONIC,
    SOURCE,
    TABLE,
    V0,
    source_bm3,
    table_rows,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/yokoo-2009-gold-thermal-reconstruction.json"
RESIDUALS = ROOT / "docs/data/yokoo-2009-gold-thermal-residuals.csv"
FIGURE = ROOT / "docs/data/yokoo-2009-gold-thermal-reconstruction.png"


@dataclass(frozen=True)
class Parameters:
    """Cold volume is relative to ambient 300 K V0; thermal laws use V0."""

    vc_over_v0: float = 0.9901904131217347
    k0_gpa: float = 180.0
    k0_prime: float = 5.61
    gamma0: float = 2.96
    a: float = 0.45
    b: float = 4.2
    theta0_k: float = 170.0


PUBLISHED = Parameters()
NAMES = tuple(asdict(PUBLISHED))
PRIMARY_FREE = ("vc_over_v0", "k0_prime", "gamma0", "a", "b")
LOWER = Parameters(0.95, 150.0, 3.0, 1.0, 0.0, 0.1, 80.0)
UPPER = Parameters(1.0, 220.0, 8.0, 4.0, 1.0, 12.0, 300.0)


@lru_cache(maxsize=1)
def electronic_nodes():
    with (DATA / ELECTRONIC).open() as stream:
        rows = list(csv.DictReader(stream))
    return (
        np.array([float(row["temperature_k"]) for row in rows]),
        np.array(
            [float(row["gold_corrected_electronic_pressure_gpa"]) for row in rows]
        ),
    )


@lru_cache(maxsize=4)
def quadrature(order):
    nodes, weights = leggauss(order)
    return (nodes + 1) / 2, weights / 2


def pressure_gpa(ratio, temperature, parameters, *, ambient_cell_a3=V0, order=64):
    """Absolute cold + Debye phonon + interpolated electronic pressure.

    This pressure reconstruction is restricted to 0.6<=V/V0<=1, 0<=T<=3000 K.
    Those numerical bounds do not certify phase stability. The six marked
    liquid states are excluded from the primary fit and evaluated separately.
    No zero-point or 300 K subtraction term is inserted into the phonon energy.
    """
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
    # For low temperatures use the Debye integral's asymptotic value pi^4/15.
    # At x>80 the omitted integral tail is below double precision accuracy.
    integral = np.full(x.shape, np.pi**4 / 15)
    moderate = x <= 80
    z, w = quadrature(order)
    argument = x[moderate, None] * z
    integral[moderate] = x[moderate] * np.sum(
        w * argument**3 / np.expm1(argument), axis=1
    )
    energy = 9 * R * temperature[positive] * integral / x**3  # J/mol Au atoms
    molar_volume_si = ambient_cell_a3 * ratio[positive] * N_A * 1e-30 / 4
    phonon[positive] = gamma[positive] * energy / molar_volume_si / 1e9
    temperatures, electronic = electronic_nodes()
    return (
        source_bm3(ratio, p.vc_over_v0, p.k0_gpa, p.k0_prime)
        + phonon
        + np.interp(temperature, temperatures, electronic)
    )


def input_arrays():
    rows = table_rows()
    return (
        rows,
        np.array([float(row["volume_ratio"]) for row in rows]),
        np.array([float(row["temperature_k"]) for row in rows]),
        np.array([float(row["pressure_gpa"]) for row in rows]),
        np.array([row["source_phase_annotation"] == "unmarked" for row in rows]),
    )


def fit_parameters(
    ratio,
    temperature,
    pressure,
    selection,
    free=PRIMARY_FREE,
    initial=PUBLISHED,
    ambient_cell_a3=V0,
):
    """Equal pressure weights; no experimental or rounding-derived errors."""
    values = asdict(initial)

    def decode(x):
        return Parameters(**{**values, **dict(zip(free, x))})

    def residual(x):
        model = pressure_gpa(
            ratio, temperature, decode(x), ambient_cell_a3=ambient_cell_a3
        )
        return (model - pressure)[selection]

    result = least_squares(
        residual,
        [values[name] for name in free],
        bounds=(
            [asdict(LOWER)[name] for name in free],
            [asdict(UPPER)[name] for name in free],
        ),
        x_scale="jac",
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-12,
        max_nfev=2000,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return decode(result.x), result


def metrics(residual):
    return {
        "states": int(residual.size),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_gpa": float(np.max(abs(residual))),
        "mean_residual_gpa": float(np.mean(residual)),
    }


def fit_summary(parameters, residual, solid, free, selection):
    return {
        "parameters": asdict(parameters),
        "refined_parameters": list(free),
        "fixed_parameters": [name for name in NAMES if name not in free],
        "fit_states": metrics(residual[selection]),
        "unmarked_states": metrics(residual[solid]),
        "first_liquid_states": metrics(residual[~solid]),
        "all_populated_states": metrics(residual),
    }


def independent_quad_pressure(ratio, temperature, parameters):
    """Scalar SI implementation with adaptive quadrature for verification."""
    p = parameters
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
    return float(
        source_bm3(ratio, p.vc_over_v0, p.k0_gpa, p.k0_prime)
        + phonon
        + np.interp(temperature, nodes, electronic)
    )


@lru_cache(maxsize=1)
def reconstruct():
    rows, ratio, temperature, pressure, solid = input_arrays()
    all_states = np.ones(pressure.size, dtype=bool)
    fits = {}
    fitted = None
    for name, free, selection in [
        ("printed_coefficients_cold_volume_only", ("vc_over_v0",), solid),
        ("primary", PRIMARY_FREE, solid),
        ("include_first_liquid_states", PRIMARY_FREE, all_states),
        ("also_refine_cold_bulk_modulus", (*PRIMARY_FREE, "k0_gpa"), solid),
        ("also_refine_bulk_modulus_and_debye_temperature", NAMES, solid),
    ]:
        parameters, _ = fit_parameters(ratio, temperature, pressure, selection, free)
        residual = pressure_gpa(ratio, temperature, parameters) - pressure
        fits[name] = fit_summary(parameters, residual, solid, free, selection)
        if name == "primary":
            fitted = parameters
    assert fitted is not None
    residual = pressure_gpa(ratio, temperature, fitted) - pressure
    density_v0 = json.loads((DATA / SOURCE).read_text(encoding="utf-8"))[
        "volume_basis"
    ]["V0_cell_a3_derived_from_printed_density"]
    sensitivity, _ = fit_parameters(
        ratio, temperature, pressure, solid, ambient_cell_a3=density_v0
    )
    density_residual = (
        pressure_gpa(ratio, temperature, sensitivity, ambient_cell_a3=density_v0)
        - pressure
    )
    fits["unrounded_density_conversion"] = fit_summary(
        sensitivity, density_residual, solid, PRIMARY_FREE, solid
    )
    fits["unrounded_density_conversion"]["ambient_cell_a3"] = density_v0
    # Deterministic distinct starts, without randomness or altered data weights.
    starts = [
        PUBLISHED,
        replace(PUBLISHED, vc_over_v0=0.98, k0_prime=4.8, gamma0=2.5, a=0.2, b=2),
        replace(PUBLISHED, vc_over_v0=0.995, k0_prime=6.5, gamma0=3.4, a=0.7, b=7),
    ]
    multistart = []
    for initial in starts:
        parameters, _ = fit_parameters(
            ratio, temperature, pressure, solid, initial=initial
        )
        multistart.append(
            {
                "initial": asdict(initial),
                "fitted": asdict(parameters),
                "unmarked_states": metrics(
                    (pressure_gpa(ratio, temperature, parameters) - pressure)[solid]
                ),
            }
        )
    # Hold out alternate complete isochores to test interpolation in volume.
    ratios = np.sort(np.unique(ratio))
    held_out = solid & np.isin(ratio, ratios[1::2])
    train = solid & ~held_out
    cv_parameters, _ = fit_parameters(ratio, temperature, pressure, train)
    cv_residual = pressure_gpa(ratio, temperature, cv_parameters) - pressure
    comparison = []
    for name in NAMES:
        printed, reconstructed = getattr(PUBLISHED, name), getattr(fitted, name)
        comparison.append(
            {
                "parameter": name,
                "published": None if name == "vc_over_v0" else printed,
                "reconstructed": reconstructed,
                "relative_difference_percent": None
                if name == "vc_over_v0"
                else 100 * (reconstructed / printed - 1),
            }
        )
    quadrature_error = max(
        abs(float(pressure_gpa(r, t, fitted)) - independent_quad_pressure(r, t, fitted))
        for r, t in zip(ratio, temperature)
    )
    return {
        "doi": "10.1103/PhysRevB.80.104114",
        "kind": "derived_table_output_pressure_reconstruction",
        "scientific_validation": {
            "status": "not_reproduced",
            "audit_date": "2026-10-09",
            "scope": "Numerical verification of diagnostic Table III pressure reconstruction; published analytical PVT remains unreproduced",
            "verified_fields": [
                "equation",
                "units",
                "reference_state",
                "source_output_agreement",
                "reconstruction_parameters",
                "represented_range",
                "source_phase_annotations",
            ],
            "evidence": "docs/literature-reproductions/yokoo-2009-gold.md",
            "note": "Diagnostic pressure reconstruction: 156 unmarked source states, RMS 0.008923 GPa, maximum 0.021133 GPa; six first-liquid states held out. Altered phonon coefficients and inferred cold volume retain derived provenance and do not reproduce the published analytical parameter set. The pressure-convention audit identifies a thermal-increment discrepancy beyond the stated last-digit rounding test. Caloric data and original fitting weights are not PVT evaluation dependencies.",
        },
        "diagnostic_only": True,
        "published_analytical_pvt_reproduced": False,
        "author_fit_reproduced": False,
        "experimental_observations": 0,
        "qualification": "Fits derived Table III outputs, not original shock/ambient measurements. Published catalog coefficients are unchanged. This executable pressure reconstruction does not recover electronic energy, caloric properties or original shock-temperature reduction.",
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [TABLE, SOURCE, ELECTRONIC]
        },
        "equation_gpa": "Pc_BM3_GPa(V/Vc,0)+gamma(V)*9*R*T*integral_0^(theta/T)[z^3/(exp(z)-1)]/((theta/T)^3*Vm_m3_per_mol*1e9)+Pel_GPa(T)",
        "conventions": {
            "ambient_cell_a3": V0,
            "atoms_per_cell": 4,
            "thermal_volume_normalization": "Ambient 300 K V0, separate from cold Vc",
            "phonon_reference": "Absolute; no zero-point or 300 K subtraction term",
            "electronic": "Corrected Tsuchiya-Kawamura Table I, volume independent, linear interpolation, no extrapolation",
            "electronic_interpolation_author_verified": False,
            "fit_electronic_values": "All fit temperatures coincide with exact source nodes; interpolation has no effect on the fit",
            "objective": "Sum of squared pressure residuals in GPa with unit weights",
            "optimizer": "scipy.optimize.least_squares, trust-region reflective, x_scale=jac, ftol=xtol=1e-13, gtol=1e-12, max_nfev=2000",
            "debye_integration": "64-point Gauss-Legendre; asymptotic integral pi^4/15 for theta/T>80",
            "uncertainty": "No parameter errors or covariance assigned: residuals reflect rounded model outputs and reconstruction assumptions, not experimental errors",
            "bounds_lower": asdict(LOWER),
            "bounds_upper": asdict(UPPER),
            "domain": {
                "volume_ratio": [0.6, 1.0],
                "temperature_k": [0, 3000],
                "phase_stability_verified": False,
            },
        },
        "parameter_comparison": comparison,
        "fits": fits,
        "multistart": multistart,
        "volume_interpolation_holdout": {
            "held_out_volume_ratios": ratios[1::2].tolist(),
            "parameters": asdict(cv_parameters),
            "training_states": metrics(cv_residual[train]),
            "held_out_states": metrics(cv_residual[held_out]),
        },
        "adaptive_quad_max_difference_gpa": float(quadrature_error),
        "cold_volume_cell_a3": fitted.vc_over_v0 * V0,
        "states": [
            {
                **row,
                "model_pressure_gpa": float(model),
                "model_minus_table_gpa": float(delta),
                "used_in_primary_fit": bool(used),
            }
            for row, model, delta, used in zip(
                rows, pressure + residual, residual, solid
            )
        ],
    }


def check_reconstruction(saved, current):
    """Optimizer quantities allow tiny variation between numerical libraries."""
    if isinstance(current, dict):
        assert saved.keys() == current.keys()
        for name in current:
            check_reconstruction(saved[name], current[name])
    elif isinstance(current, list):
        assert len(saved) == len(current)
        for left, right in zip(saved, current):
            check_reconstruction(left, right)
    elif isinstance(current, float):
        np.testing.assert_allclose(saved, current, atol=2e-6, rtol=1e-5)
    else:
        assert saved == current


def plot_reconstruction(result):
    """Show residuals at actual source states; no manufactured observations."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, ratio, temperature, pressure, solid = input_arrays()
    del rows
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.8), sharex=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, 8))
    for ax, name, title in zip(
        axes,
        ["printed_coefficients_cold_volume_only", "primary"],
        [
            "Printed coefficients; cold volume reconstructed",
            "Reconstructed coefficients; 156 unmarked states fitted",
        ],
    ):
        parameters = Parameters(**result["fits"][name]["parameters"])
        residual = pressure_gpa(ratio, temperature, parameters) - pressure
        for t, color in zip(np.unique(temperature), colors):
            select = (temperature == t) & solid
            ax.plot(
                ratio[select],
                residual[select],
                ".-",
                color=color,
                markersize=4,
                linewidth=0.8,
                label=f"{t:g} K",
            )
            liquid = (temperature == t) & ~solid
            ax.plot(
                ratio[liquid],
                residual[liquid],
                "D",
                markerfacecolor="none",
                markeredgecolor=color,
                markersize=6,
            )
        ax.axhline(0, color="0.4", linewidth=0.7)
        ax.grid(alpha=0.2)
        ax.set_ylabel("Model − Table III (GPa)")
        ax.set_title(title, fontsize=11, loc="left")
    axes[0].legend(ncol=4, fontsize=8, loc="lower left")
    axes[1].set_xlabel("Volume / ambient 300 K volume")
    fig.suptitle("Au · Yokoo et al. (2009) · pressure reconstruction", fontsize=14)
    fig.text(
        0.5,
        0.01,
        "Table III contains model outputs. Open diamonds: first liquid states held out.\n"
        "Agreement measures reconstruction; original experimental fit is unverified.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.97))
    fig.savefig(FIGURE, dpi=180)
    fig.savefig(FIGURE.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--pressure", nargs=2, type=float, metavar=("V_OVER_V0", "T_K"))
    args = parser.parse_args()
    result = reconstruct()
    if args.pressure:
        parameters = Parameters(**result["fits"]["primary"]["parameters"])
        ratio, temperature = args.pressure
        print(
            json.dumps(
                {
                    "pressure_gpa": float(pressure_gpa(ratio, temperature, parameters)),
                    "volume_ratio": ratio,
                    "temperature_k": temperature,
                    "kind": result["kind"],
                    "phase_stability_verified": False,
                },
                indent=2,
            )
        )
        return
    if args.check:
        check_reconstruction(json.loads(OUTPUT.read_text(encoding="utf-8")), result)
        with RESIDUALS.open() as stream:
            archived_rows = list(csv.DictReader(stream))
        assert len(archived_rows) == len(result["states"])
        for archived, current in zip(archived_rows, result["states"]):
            assert archived.keys() == current.keys()
            for name, value in current.items():
                if isinstance(value, float):
                    check_reconstruction(float(archived[name]), value)
                else:
                    assert archived[name] == str(value)
    else:
        OUTPUT.write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        with RESIDUALS.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(result["states"][0]))
            writer.writeheader()
            writer.writerows(result["states"])
    if args.plot:
        plot_reconstruction(result)
    print(
        json.dumps(
            {
                "parameters": result["fits"]["primary"]["parameters"],
                "fit": result["fits"]["primary"]["unmarked_states"],
                "liquid_holdout": result["fits"]["primary"]["first_liquid_states"],
                "volume_holdout": result["volume_interpolation_holdout"][
                    "held_out_states"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
