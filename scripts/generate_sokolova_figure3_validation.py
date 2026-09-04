#!/usr/bin/env python3
"""Generate the Sokolova et al. (2016) MgO Figure 3 equation audit plot."""

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import R
from scipy.optimize import brentq

from peritheos.eos.thermal.multi_oscillator import I_gamV
from peritheos.materials import MGO_SOKOLOVA_2013

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "docs" / "data" / "sokolova-2016-figure3-kt-digitized.csv"
OUTPUT_PATH = ROOT / "docs" / "images" / "sokolova-2016-figure3-validation.svg"


def paper_pressure(model, ratio, temperature):
    """Evaluate literal Sokolova et al. (2013) equations (6), (9), and (14)."""
    rt_eos = model.rt_eos
    volume = ratio * rt_eos.V0
    reference_pressure = rt_eos.pressure(volume)
    bulk_modulus = rt_eos.bulk_modulus(volume)
    bulk_modulus_derivative = rt_eos.bulk_modulus_derivative(volume)
    generalized_t = model.t - model.beta * np.cbrt(ratio)
    gamma = (
        -3.0 * bulk_modulus
        + 2.0 * reference_pressure * generalized_t
        + 9.0 * bulk_modulus * bulk_modulus_derivative
        - 6.0 * generalized_t * bulk_modulus
    ) / (6.0 * (3.0 * bulk_modulus - 2.0 * reference_pressure * generalized_t))
    gamma += model.delta
    volume_theta_factor = math.exp(
        I_gamV(ratio, model.delta, model.t, rt_eos, model.beta)
    )
    anharmonicity = model.a_0 * 1.0e-6 * ratio**model.m

    def absolute_thermal_pressure(current_temperature):
        pressure_factor = gamma - 0.5 * model.m * anharmonicity * current_temperature
        pressure_bar = 0.0
        for multiplicity, reference_theta in (
            (model.mE1, model.QE1o),
            (model.mE2, model.QE2o),
        ):
            theta = (
                reference_theta
                * volume_theta_factor
                * math.exp(0.5 * anharmonicity * current_temperature)
            )
            pressure_bar += (
                multiplicity
                * R
                * theta
                / math.expm1(theta / current_temperature)
                * pressure_factor
                / volume
            )
        pressure_bar += (
            1.5
            * model.n
            * R
            * model.e_0
            * 1.0e-6
            * ratio**model.g
            * model.g
            * current_temperature**2
            / volume
        )
        return pressure_bar / 10000.0

    return (
        reference_pressure
        + absolute_thermal_pressure(temperature)
        - absolute_thermal_pressure(model.Tr)
    )


def paper_zero_pressure_bulk_modulus(model, temperature):
    """Calculate K_T on the literal paper equation's zero-pressure branch."""
    ratio = brentq(
        lambda current_ratio: paper_pressure(model, current_ratio, temperature),
        1.0,
        1.2,
    )
    step = ratio * 1.0e-5
    pressure_derivative = (
        paper_pressure(model, ratio + step, temperature)
        - paper_pressure(model, ratio - step, temperature)
    ) / (2.0 * step)
    return -ratio * pressure_derivative


def load_digitized_rows():
    with DATA_PATH.open(newline="") as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def main():
    rows = load_digitized_rows()
    eos = MGO_SOKOLOVA_2013.eos
    temperatures = np.linspace(450.0, 2050.0, 161)
    workbook_curve = []
    paper_curve = []
    for temperature in temperatures:
        volume = eos.calculate_volume(0.0, temperature)
        workbook_curve.append(eos.bulk_modulus(volume, temperature))
        paper_curve.append(paper_zero_pressure_bulk_modulus(eos, temperature))

    observed_t = np.array([row["temperature_k"] for row in rows])
    observed_k = np.array([row["isothermal_bulk_modulus_gpa"] for row in rows])
    observed_t_error = np.array(
        [row["temperature_digitization_uncertainty_k"] for row in rows]
    )
    observed_k_error = np.array(
        [row["bulk_modulus_digitization_uncertainty_gpa"] for row in rows]
    )
    workbook_residual = np.array([row["workbook_residual_gpa"] for row in rows])
    paper_residual = np.array([row["printed_equation_residual_gpa"] for row in rows])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "svg.hashsalt": "sokolova-figure3-equation-audit",
        }
    )
    figure, (axis, residual_axis) = plt.subplots(
        2,
        1,
        figsize=(7.2, 5.8),
        sharex=True,
        gridspec_kw={"height_ratios": (3.0, 1.2), "hspace": 0.08},
    )

    workbook_color = "#0072B2"
    paper_color = "#D55E00"
    observed_color = "#202020"
    axis.plot(
        temperatures,
        workbook_curve,
        color=workbook_color,
        linewidth=2.2,
        label="VBA/workbook equation",
    )
    axis.plot(
        temperatures,
        paper_curve,
        color=paper_color,
        linewidth=2.2,
        linestyle="--",
        label="Literal printed equation",
    )
    axis.errorbar(
        observed_t,
        observed_k,
        xerr=observed_t_error,
        yerr=observed_k_error,
        color=observed_color,
        marker="o",
        markerfacecolor="white",
        markeredgewidth=1.2,
        markersize=5.0,
        linestyle="none",
        capsize=2.5,
        linewidth=0.9,
        label="Digitized published curve",
        zorder=3,
    )
    axis.set_title("Sokolova et al. (2016) MgO Figure 3 equation validation")
    axis.set_ylabel(r"Zero-pressure $K_T$ (GPa)")
    axis.grid(True, linewidth=0.5, alpha=0.25)
    axis.legend(frameon=False, loc="lower left")
    axis.text(
        0.98,
        0.97,
        "RMS residual\nVBA/workbook: 0.47 GPa\nPrinted: 9.77 GPa",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        bbox={"facecolor": "white", "edgecolor": "#B0B0B0", "pad": 4.0},
    )

    residual_axis.axhline(0.0, color="#707070", linewidth=0.8)
    residual_axis.errorbar(
        observed_t,
        workbook_residual,
        xerr=observed_t_error,
        yerr=observed_k_error,
        color=workbook_color,
        marker="o",
        markersize=4.5,
        linestyle="-",
        linewidth=1.1,
        capsize=2.0,
        label="VBA/workbook",
    )
    residual_axis.errorbar(
        observed_t,
        paper_residual,
        xerr=observed_t_error,
        yerr=observed_k_error,
        color=paper_color,
        marker="s",
        markersize=4.0,
        linestyle="--",
        linewidth=1.1,
        capsize=2.0,
        label="Printed",
    )
    residual_axis.set_xlim(400.0, 2100.0)
    residual_axis.set_xlabel("Temperature (K)")
    residual_axis.set_ylabel("Model - digitized\n(GPa)")
    residual_axis.grid(True, linewidth=0.5, alpha=0.25)
    residual_axis.legend(frameon=False, loc="lower left", ncol=2)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        OUTPUT_PATH,
        bbox_inches="tight",
        metadata={"Date": None, "Creator": "Peritheos equation audit"},
    )
    plt.close(figure)
    svg = OUTPUT_PATH.read_text(encoding="utf-8")
    OUTPUT_PATH.write_text(
        "\n".join(line.rstrip() for line in svg.splitlines()) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
