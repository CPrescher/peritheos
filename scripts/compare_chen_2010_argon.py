"""Compare the drawn Chen curve, reported-constant reconstruction and refit.

This is a diagnostic comparison, not a verification of the author's equation.
All differences compare volumes at the same pressure within the observed range.
"""

import argparse
import json

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from scripts.refit_chen_2010_argon import observations, pressure
from scripts.reproduce_chen_2010_argon import (
    ROOT,
    local_bm3_coefficients,
    polynomial_pressure,
    read_points,
    reproduce,
)

ORIGINAL = "argon_fcc_chen_2010_bm3_reported_constants"
REFIT = "argon_fcc_chen_2010_bm3_digitized_refit"
REPORT = ROOT / "docs/data/chen-2010-argon-comparison.json"
FIGURE = ROOT / "docs/data/chen-2010-argon-comparison.png"


def compare():
    refit = json.loads(
        (ROOT / "docs/data/chen-2010-argon-refit.json").read_text(encoding="utf-8")
    )
    source = reproduce()
    pars = refit["primary"]["density_parameters"]
    p, _, _, _ = observations()
    curve = read_points("figure5-curve")
    cp = np.array([float(row["pressure_gpa"]) for row in curve])
    cr = np.array([float(row["density_g_cm3"]) for row in curve])
    line = PchipInterpolator(cp, cr, extrapolate=False)
    coefficients = local_bm3_coefficients()

    def densities(pressures):
        return {
            "refit": [
                brentq(lambda r: pressure(r, pars) - x, pars[0], 5) for x in pressures
            ],
            "published_curve": line(pressures).tolist(),
            "reported_constants": [
                brentq(lambda r: polynomial_pressure(r, coefficients) - x, 1.68, 5)
                for x in pressures
            ],
            "pressure_offset_alternative": [
                brentq(lambda r: 2 + pressure(r, [2.18, 15.1, 5.4]) - x, 1.8, 5)
                for x in pressures
            ],
        }

    grid = np.linspace(min(p), max(p), 1000)
    rho = densities(grid)
    differences = {
        name: (100 * (np.array(rho["refit"]) / values - 1)).tolist()
        for name, values in rho.items()
        if name != "refit"
    }
    checkpoints = []
    for x in [2, 5, 10, 20, 26]:
        values = {key: value[0] for key, value in densities([x]).items()}
        checkpoints.append(
            {
                "pressure_gpa": x,
                "density_g_cm3": values,
                "volume_difference_from_refit_percent": {
                    key: 100 * (values["refit"] / value - 1)
                    for key, value in values.items()
                    if key != "refit"
                },
            }
        )
    return {
        "original_record_identifier": ORIGINAL,
        "refit_record_identifier": REFIT,
        "original_status": "not_reproduced",
        "original_usage": "DO NOT USE for scientific predictions or pressure calibration. Diagnostic reconstruction only; the author's exact equation is unverified.",
        "cause": "Unresolved. No error in the author's pressure calibration or adiabatic correction has been established.",
        "temperature_k": 290,
        "pressure_range_gpa": [float(min(p)), float(max(p))],
        "difference_definition": "100 * (V_comparison(P) / V_refit(P) - 1)",
        "curve_interpolation": "Shape-preserving PCHIP through digitized Figure 5 line vertices; no extrapolation.",
        "checkpoints": checkpoints,
        "max_abs_volume_difference_percent": {
            key: float(np.max(np.abs(value))) for key, value in differences.items()
        },
        "reported_constants_reconstruction": source[
            "standard_bm3_local_constraint_diagnostic"
        ],
        "refit_parameters": refit["primary"]["parameters"],
        "refit_pressure_rms_original_coordinates_gpa": refit["primary"][
            "pressure_rms_original_coordinates_gpa"
        ],
        "refit_at_2_gpa": refit["refit_at_2_gpa"],
        "published_curve_at_2_gpa": source["figure5_curve_diagnostic"],
        "plot": {
            "pressure_gpa": grid.tolist(),
            "density_g_cm3": rho,
            "volume_difference_from_refit_percent": differences,
        },
    }


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p, rho, sp, sr = observations()
    data = report["plot"]
    grid = data["pressure_gpa"]
    colors = {
        "refit": "#bd3b35",
        "published_curve": "#303946",
        "reported_constants": "#346ea0",
        "pressure_offset_alternative": "#ad6d15",
    }
    labels = {
        "refit": "Peritheos BM3 refit",
        "published_curve": "Published Figure 5 curve",
        "reported_constants": "Printed-constant reconstruction: DO NOT USE",
        "pressure_offset_alternative": "Alternative: 2 GPa pressure offset",
    }
    styles = {
        "refit": "-",
        "published_curve": "--",
        "reported_constants": "-",
        "pressure_offset_alternative": ":",
    }
    with plt.rc_context(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    ):
        fig = plt.figure(figsize=(13, 7.3), layout="constrained")
        gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1])
        ax = fig.add_subplot(gs[:, 0])
        zoom = fig.add_subplot(gs[0, 1])
        delta = fig.add_subplot(gs[1, 1])
        ax.errorbar(
            p,
            rho,
            xerr=sp,
            yerr=sr,
            fmt="o",
            ms=2.5,
            color="#8c98a4",
            ecolor="#d7dce1",
            elinewidth=0.6,
            label="80 digitized Brillouin density positions",
        )
        for name, values in data["density_g_cm3"].items():
            ax.plot(
                grid,
                values,
                color=colors[name],
                ls=styles[name],
                lw=2,
                label=labels[name],
            )
        ax.set(
            xlabel="Pressure (GPa)",
            ylabel="Density (g/cm³)",
            title="Density at the same pressure",
        )
        ax.legend(loc="lower right", fontsize=8)
        differences = data["volume_difference_from_refit_percent"]
        zoom.plot(
            grid, differences["published_curve"], color=colors["published_curve"], lw=2
        )
        zoom.set(
            title="Published curve vs refit — enlarged", ylabel="Volume difference (%)"
        )
        zoom.text(
            0.35,
            0.88,
            f"Maximum absolute difference: {report['max_abs_volume_difference_percent']['published_curve']:.2f}%",
            transform=zoom.transAxes,
        )
        for name in ["reported_constants", "pressure_offset_alternative"]:
            delta.plot(
                grid,
                differences[name],
                color=colors[name],
                ls=styles[name],
                lw=2,
                label=labels[name],
            )
        delta.set(
            title="Printed-constant reconstructions vs refit",
            xlabel="Pressure (GPa)",
            ylabel="Volume difference (%)",
        )
        delta.legend(fontsize=8, loc="upper left")
        for axis in [zoom, delta]:
            axis.axhline(0, color=colors["refit"], lw=1)
            axis.set_xlim(min(grid), max(grid))
        for axis in [ax, zoom, delta]:
            axis.grid(alpha=0.15)
        fig.suptitle(
            "Chen et al. (2010), argon: published curve, printed constants and independent refit",
            fontsize=15,
        )
        fig.supxlabel(
            "Digitized interval: 1.23–26.06 GPa, nominal 290 K. Difference = 100 × (Vcomparison / Vrefit − 1).\nReconstructions are diagnostics, not verified author equations. Error bars are graphical halfwidths; confidence levels are unknown.",
            fontsize=9,
        )
        fig.savefig(FIGURE, dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = compare()
    REPORT.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    if args.plot:
        plot(report)
