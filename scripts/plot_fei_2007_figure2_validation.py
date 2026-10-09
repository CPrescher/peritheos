"""Compare published Pt coefficients with calculated curves in Fei Figure 2.

Run with ``python -m scripts.plot_fei_2007_figure2_validation``. No coefficients
are fitted. Raster samples and calibration are preserved in the input JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_fei_2007_thermal_parts import PARAMETERS, pressure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"
INPUT = DATA / "fei-2007-platinum-source-curve-check.json"
OUTPUT = DATA / "fei-2007-platinum-figure2-validation.png"


def main():
    report = json.loads(INPUT.read_text())
    pc = np.polyfit(
        report["pressure_axis_ticks"]["x_pixel"],
        report["pressure_axis_ticks"]["pressure_gpa"],
        1,
    )
    vc = np.polyfit(
        report["volume_axis_ticks"]["y_pixel"],
        report["volume_axis_ticks"]["volume_a3"],
        1,
    )
    points = report["points"]
    # Recompute both calibration and model predictions rather than plot saved
    # pressures blindly. The pressure evaluator uses independent SI quadrature.
    for row in points:
        volume = float(np.polyval(vc, row["y_pixel"]))
        figure_pressure = float(np.polyval(pc, row["x_pixel"]))
        calculated = float(pressure(volume, row["temperature_k"], "platinum"))
        assert abs(volume - row["volume_a3"]) < 1e-10
        assert abs(calculated - row["published_coefficients_pressure_gpa"]) < 1e-9
        assert abs(figure_pressure - row["figure_pressure_gpa"]) < 1e-10
        # Readout sensitivity combines the colored stroke's observed horizontal
        # half-width with an illustrative one pixel per calibrated coordinate.
        dv = abs(vc[0])
        volume_sensitivity = max(
            abs(
                float(pressure(volume + sign * dv, row["temperature_k"], "platinum"))
                - calculated
            )
            for sign in (-1, 1)
        )
        lo, hi = row["x_pixel_span"]
        row["stroke_and_one_pixel_axis_sensitivity_gpa"] = float(
            ((hi - lo) / 2 + 1) * abs(pc[0]) + volume_sensitivity
        )

    colors = {300: "#2525cc", 1473: "#cc00cc", 1873: "#dd2222"}
    summary = {}
    for temperature in colors:
        selected = [r for r in points if r["temperature_k"] == temperature]
        difference = np.array([r["difference_gpa"] for r in selected])
        summary[str(temperature)] = {
            "samples": len(selected),
            "rmse_gpa": float(np.sqrt(np.mean(difference**2))),
            "mean_difference_gpa": float(difference.mean()),
            "max_absolute_difference_gpa": float(np.max(abs(difference))),
        }
    lower = {r["y_pixel"]: r for r in points if r["temperature_k"] == 1473}
    upper = {r["y_pixel"]: r for r in points if r["temperature_k"] == 1873}
    separation = [
        {
            "volume_a3": lower[y]["volume_a3"],
            "figure_separation_gpa": upper[y]["figure_pressure_gpa"]
            - lower[y]["figure_pressure_gpa"],
            "calculated_separation_gpa": upper[y]["published_coefficients_pressure_gpa"]
            - lower[y]["published_coefficients_pressure_gpa"],
            "difference_gpa": upper[y]["difference_gpa"] - lower[y]["difference_gpa"],
        }
        for y in lower
    ]
    differences = np.array([r["difference_gpa"] for r in separation])
    report["per_isotherm_summary"] = summary
    report["hot_curve_separation"] = {
        "description": "1873 K minus 1473 K pressure at the same volume; the cold pressure cancels. No coefficients refitted.",
        "rmse_gpa": float(np.sqrt(np.mean(differences**2))),
        "max_absolute_difference_gpa": float(np.max(abs(differences))),
        "points": separation,
    }
    report["readout_sensitivity_note"] = (
        "Error bars in the validation plot show observed horizontal half-stroke "
        "width plus an illustrative one pixel in each calibrated coordinate. "
        "These are deterministic readout sensitivities, not statistical errors. "
        "The earlier three-pixel criterion alone is not proof of exact parity."
    )
    report["published_coefficients"] = dict(
        zip(("V0", "K0", "K0_prime", "gamma0", "q", "theta0"), PARAMETERS["platinum"])
    )
    INPUT.write_text(json.dumps(report, indent=2) + "\n")

    plt.rcParams.update(
        {"font.size": 11, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig = plt.figure(figsize=(15, 8.5), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1])
    source_ax = fig.add_subplot(gs[:, 0])
    source_ax.imshow(plt.imread(ROOT / report["source_image"]))
    volume = np.linspace(47.05, 60.38, 800)
    for temperature in colors:
        calculated = pressure(volume, temperature, "platinum")
        x = (calculated - pc[1]) / pc[0]
        y = (volume - vc[1]) / vc[0]
        inside = (x >= 216) & (x <= 1780)
        source_ax.plot(x[inside], y[inside], "k--", linewidth=1.05, dashes=(5, 3))
    source_ax.set_xlim(0, 1800)
    source_ax.set_ylim(1400, 0)
    source_ax.axis("off")
    source_ax.set_title(
        "Original Figure 2 + published-coefficient calculation (black dashed)",
        fontsize=12,
    )

    residual_ax = fig.add_subplot(gs[0, 1])
    for temperature, color in colors.items():
        selected = [r for r in points if r["temperature_k"] == temperature]
        residual_ax.errorbar(
            [r["figure_pressure_gpa"] for r in selected],
            [r["difference_gpa"] for r in selected],
            yerr=[r["stroke_and_one_pixel_axis_sensitivity_gpa"] for r in selected],
            fmt="o",
            markersize=4,
            capsize=2,
            color=color,
            label=f"{temperature} K; RMS {summary[str(temperature)]['rmse_gpa']:.2f} GPa",
        )
    residual_ax.axhline(0, color="0.4", linewidth=0.8)
    residual_ax.set_xlabel("Pressure read from Figure 2 (GPa)")
    residual_ax.set_ylabel("Calculated - figure pressure (GPa)")
    residual_ax.set_title("Sampled curve residuals; no refitting")
    residual_ax.legend(fontsize=9, loc="lower right")
    residual_ax.grid(alpha=0.15)

    thermal_ax = fig.add_subplot(gs[1, 1])
    thermal_ax.plot(
        [r["volume_a3"] for r in separation],
        [r["figure_separation_gpa"] for r in separation],
        "o",
        markerfacecolor="white",
        markeredgecolor="black",
        label="Figure 2 curve separation",
    )
    thermal_ax.plot(
        volume,
        pressure(volume, 1873, "platinum") - pressure(volume, 1473, "platinum"),
        color="#008066",
        label="Published coefficients",
    )
    thermal_ax.set_xlim(47.5, 55.5)
    thermal_ax.set_xlabel("Pt cell volume ($\\AA^3$)")
    thermal_ax.set_ylabel("$P_{1873} - P_{1473}$ (GPa)")
    thermal_ax.set_title(
        "Thermal increment: cold pressure cancels\nRMS difference 0.032 GPa; maximum 0.073 GPa",
        fontsize=11,
    )
    thermal_ax.legend(fontsize=9)
    thermal_ax.grid(alpha=0.15)
    fig.suptitle(
        "Fei et al. (2007), Pt: Table 1 coefficients checked against Figure 2",
        fontsize=15,
    )
    fig.supxlabel(
        "Vinet + printed Debye law; V0=60.38, K0=277, K0'=5.08, gamma0=2.72, q=0.5, theta0=230 K.\n"
        "Error bars: line half-width + illustrative 1 pixel per axis; not statistical errors. Source: doi:10.1073/pnas.0609013104",
        fontsize=9,
    )
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)
    print(
        json.dumps(
            {
                "isotherms": summary,
                "thermal_separation": report["hot_curve_separation"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
