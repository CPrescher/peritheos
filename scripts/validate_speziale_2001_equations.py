"""Equation-only checks against Speziale (2001), without fitting any data.

python -m scripts.validate_speziale_2001_equations [--check] [--plot]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expi

from scripts.reproduce_speziale_2001_mgo import (
    PARAMETERS,
    adaptive_pressure,
    gamma_theta,
    pressure,
)

ROOT = Path(__file__).resolve().parents[1]
CURVES = ROOT / "docs/data/speziale-2001-source-curve-checkpoints.csv"
PROVENANCE = ROOT / "docs/data/speziale-2001-source-curve-provenance.json"
OUTPUT = ROOT / "docs/data/speziale-2001-equation-validation.json"
PLOT = OUTPUT.with_suffix(".png")


def closed_form_theta(volume):
    """Integrate Eq. 11 analytically using the exponential integral, Ei.

    This independent expression applies to the fixed positive q0/q1 of this
    source, not arbitrary parameters near the constant-q or q0=0 limits.
    """
    x = np.asarray(volume, dtype=float) / 74.71
    a = 1.65 / 11.8
    log_theta_ratio = -1.524 * np.exp(-a) / 11.8 * (expi(a * x**11.8) - expi(a))
    return 773 * np.exp(log_theta_ratio)


def closed_form_pressure(volume, temperature):
    """Independent Ei theta, adaptive Debye integral, finite-strain BM3."""
    x = volume / 74.71
    gamma = 1.524 * np.exp(1.65 / 11.8 * (x**11.8 - 1))
    theta = float(closed_form_theta(volume))

    def e(t):
        y = theta / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-12, epsrel=1e-12)[
            0
        ]
        return 18 * R * t * integral / y**3

    f = (x ** (-2 / 3) - 1) / 2
    cold = 3 * 160.2 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (3.99 - 4) * f)
    molar_volume = volume * 1e-30 * N_A / 4
    return cold + gamma * (e(temperature) - e(300)) / molar_volume / 1e9


def source_curve_rows():
    with CURVES.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def validate():
    from scripts.reproduce_fei_2007_gold import mgo_pressure

    ratios = np.array([0.60, 0.64, 0.70, 0.80, 0.90, 1.0, 1.02])
    states = [(74.71 * x, t) for x in ratios for t in (300, 500, 1100, 2000, 3663)]
    gauss = np.array([pressure(v, t) for v, t in states])
    adaptive = np.array([adaptive_pressure(v, t) for v, t in states])
    analytical = np.array([closed_form_pressure(v, t) for v, t in states])
    existing = np.array([mgo_pressure(v, t) for v, t in states])
    _, theta = gamma_theta(74.71 * ratios)
    numerical = dict(
        states=len(states),
        volume_ratio_range=[float(min(ratios)), float(max(ratios))],
        temperature_range_k=[300, 3663],
        adaptive_max_difference_gpa=float(np.max(np.abs(gauss - adaptive))),
        closed_form_max_difference_gpa=float(np.max(np.abs(gauss - analytical))),
        existing_calibration_max_difference_gpa=float(np.max(np.abs(gauss - existing))),
        closed_form_theta_max_difference_k=float(
            np.max(np.abs(theta - closed_form_theta(74.71 * ratios)))
        ),
    )
    checkpoints = []
    for r in source_curve_rows():
        ratio, t, p = (
            float(r[k]) for k in ("volume_ratio", "temperature_k", "pressure_gpa")
        )
        kwargs = {"q1": 0} if r["gamma_law"] == "constant_q" else {}
        calculated_p = float(pressure(74.71 * ratio, t, **kwargs))
        calculated_ratio = brentq(
            lambda x: pressure(74.71 * x, t, **kwargs) - p, 0.7, 1.15
        )
        checkpoints.append(
            dict(
                figure=int(r["figure"]),
                gamma_law=r["gamma_law"],
                temperature_k=t,
                digitized_pressure_gpa=p,
                digitized_volume_ratio=ratio,
                calculated_pressure_gpa=calculated_p,
                pressure_difference_gpa=calculated_p - p,
                calculated_volume_ratio=calculated_ratio,
                volume_ratio_difference=calculated_ratio - ratio,
            )
        )
    figures = {}
    for fig in (6, 10):
        rows = [r for r in checkpoints if r["figure"] == fig]
        figures[str(fig)] = dict(
            checkpoints=len(rows),
            max_abs_pressure_difference_gpa=max(
                abs(r["pressure_difference_gpa"]) for r in rows
            ),
            max_abs_volume_ratio_difference=max(
                abs(r["volume_ratio_difference"]) for r in rows
            ),
            interpretation=(
                "Constant-q isotherms agree at graphical precision."
                if fig == 6
                else "Small plotted-curve differences remain; the inset does not "
                "identify the original variable-q theta prescription."
            ),
            exact_author_calculation_verified=False,
        )
    return dict(
        scope="Printed equations and numerical evaluation only; no fitting or later pressure-target matching.",
        parameters=PARAMETERS,
        normalization=dict(
            cell_volume_units="angstrom^3",
            formula_units_per_conventional_cell=4,
            atoms_per_formula_unit=2,
            molar_volume_expression="V_cell * 1e-30 * N_A / 4 (m^3/mol MgO)",
            energy_expression="9*n*R*T*(T/theta)^3*integral_0^(theta/T) z^3/(exp(z)-1) dz",
            thermal_pressure_expression="gamma*(E(T)-E(300))/V_m/1e9 (GPa)",
        ),
        theta_interpretation="Integrate d ln(theta) = -gamma d ln(V); printed Eq. 4 is the constant-q case.",
        numerical=numerical,
        gamma_infinite_compression_limit=float(1.524 * np.exp(-1.65 / 11.8)),
        source_curve_checks=figures,
        checkpoints=checkpoints,
        source_curve_csv_sha256=hashlib.sha256(CURVES.read_bytes()).hexdigest(),
        source_curve_provenance_sha256=hashlib.sha256(
            PROVENANCE.read_bytes()
        ).hexdigest(),
        printed_gamma_equation_verified=True,
        consistent_theta_integration_verified=True,
        energy_and_pressure_normalization_verified=True,
        exact_historical_variable_q_implementation_verified=False,
    )


def plot(report):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
    for ax, source_fig, temps, pmax in zip(
        axes, (6, 10), ((300, 1100, 3000), (1100,)), (50, 25)
    ):
        kwargs = {"q1": 0} if source_fig == 6 else {}
        for t, color in zip(temps, ("#377eb8", "#e07b22", "#28886b")):
            ratios = np.linspace(0.81 if source_fig == 6 else 0.90, 1.01, 300)
            ax.plot(
                pressure(74.71 * ratios, t, **kwargs),
                ratios,
                color=color,
                label=f"Calculated {t} K",
            )
            rows = [
                r
                for r in report["checkpoints"]
                if r["figure"] == source_fig and r["temperature_k"] == t
            ]
            ax.scatter(
                [r["digitized_pressure_gpa"] for r in rows],
                [r["digitized_volume_ratio"] for r in rows],
                color=color,
                marker="x",
                s=45,
            )
        ax.set(
            xlim=(0, pmax),
            ylim=(0.81 if source_fig == 6 else 0.9, 1.01),
            xlabel="Pressure (GPa)",
            ylabel="$V/V_0$",
        )
        ax.set_title(
            f"Figure {source_fig}: "
            + (
                "constant q = 1.65"
                if source_fig == 6
                else "Eq. 11, integrated Debye temperature"
            ),
            fontsize=10,
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle(
        "Speziale (2001): equation calculation vs. published calculated curves",
        fontsize=12,
    )
    fig.supxlabel(
        "Crosses: digitized source curves, not experimental data. Graphical precision only.",
        fontsize=9,
    )
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = validate()
    if args.check:
        from scripts.reproduce_fei_2007_gold import check_saved

        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
    else:
        OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.plot:
        plot(report)
    print(json.dumps(report["numerical"], indent=2))
    print(json.dumps(report["source_curve_checks"], indent=2))


if __name__ == "__main__":
    main()
