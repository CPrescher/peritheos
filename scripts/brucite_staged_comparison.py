"""Plot independent staged brucite refits, preserving the published EOS.

Stage 1 fits constant-alpha ambient expansion; stage 2 fixes its V0(T) and
fits each isotherm independently using unweighted pressure residuals.
No continuous modulus interpolation or catalog mutation is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
from scipy.optimize import least_squares

from scripts.brucite_reproduce import (
    DATA,
    ROOT,
    bm3,
    comparison_volume,
    metrics,
    observations,
    reference_volume,
    source,
)

REPORT = ROOT / "docs/data/brucite-staged-comparison.json"


def reproduce():
    ambient = observations("fukui_2003", "table1")
    temperature = np.array([float(row["temperature_k"]) for row in ambient])
    volume = np.array([float(row["cell_volume_a3"]) for row in ambient])
    fit = least_squares(
        lambda x: x[0] * np.exp(x[1] * 1e-5 * (temperature - 300)) - volume,
        [40.917, 6.69],
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
    )
    scale = np.array([1, 1e-5])
    covariance = (
        np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun @ fit.fun) / (len(volume) - 2)
    )
    covariance *= scale[:, None] * scale[None, :]
    anchor, alpha = float(fit.x[0]), float(fit.x[1] * 1e-5)
    stage1 = {
        "V300_a3": anchor,
        "alpha_per_k": alpha,
        "source_rows": list(range(1, 12)),
        "solver_success": bool(fit.success),
        "objective": "Unweighted sum over all 11 Table 1 rows of [V300*exp(alpha*(T-300))-Vobserved]^2, in A^6; V300 and constant alpha free.",
        "volume_residuals_a3": fit.fun.tolist(),
        "volume_rmse_a3": float(np.sqrt(np.mean(fit.fun**2))),
        "conditional_covariance": {
            "parameter_order": ["V300_a3", "alpha_per_k"],
            "matrix": covariance.tolist(),
        },
    }
    rows = observations("fukui_2003", "table4")
    isotherms = []
    for published in source()["sources"]["fukui_2003"]["parameters"]["isotherms"]:
        t = published["temperature_k"]
        group = [row for row in rows if float(row["temperature_k"]) == t]
        v = np.array([float(row["cell_volume_a3"]) for row in group])
        p = np.array([float(row["pressure_gpa"]) for row in group])
        v0 = anchor * np.exp(alpha * (t - 300))
        refit = least_squares(
            lambda x: bm3(v, v0, *x) - p,
            [published["K0_gpa"], published["K0_prime"]],
            bounds=([1, 1], [200, 20]),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        covariance = (
            np.linalg.inv(refit.jac.T @ refit.jac)
            * (refit.fun @ refit.fun)
            / (len(v) - 2)
        )
        variants = {}
        for key, parameters in {
            "published": [
                float(reference_volume(t)),
                published["K0_gpa"],
                published["K0_prime"],
            ],
            "staged_refit": [float(v0), float(refit.x[0]), float(refit.x[1])],
        }.items():
            pressure_residual = bm3(v, *parameters) - p
            predicted_volume = comparison_volume(p, *parameters)
            np.testing.assert_allclose(
                bm3(predicted_volume, *parameters), p, atol=1e-10
            )
            volume_residual = predicted_volume - v
            variants[key] = {
                "V0_a3": parameters[0],
                "K0_gpa": parameters[1],
                "K0_prime": parameters[2],
                "pressure_residuals_gpa": pressure_residual.tolist(),
                "volume_residuals_a3": volume_residual.tolist(),
                **metrics(pressure_residual),
                "volume_rmse_a3": float(np.sqrt(np.mean(volume_residual**2))),
            }
        isotherms.append(
            {
                "temperature_k": t,
                "observations": len(group),
                "source_rows": [int(row["source_row"]) for row in group],
                "pressure_range_gpa": [float(p.min()), float(p.max())],
                "solver_success": bool(refit.success),
                "conditional_covariance_fixed_V0": {
                    "parameter_order": ["K0_gpa", "K0_prime"],
                    "matrix": covariance.tolist(),
                },
                **variants,
            }
        )
    return {
        "format": "peritheos.brucite-staged-comparison",
        "format_version": 1,
        "primary_source_doi": source()["sources"]["fukui_2003"]["doi"],
        "source_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "brucite-source.json",
                "brucite-fukui-2003-table1.csv",
                "brucite-fukui-2003-table4.csv",
            ]
        },
        "stage1_ambient": stage1,
        "stage2_objective": "For each temperature independently, minimize sum over its Table 4 rows of [P_BM3(Vobserved;V0(T),K0,Kprime)-Pobserved]^2 [GPa^2]. V0(T) fixed from stage 1; K0 and Kprime free. No weighting, exclusions, pressure recalibration or modulus-temperature interpolation.",
        "covariance_qualification": "Local SSE/(N-p)*inverse(J.T J). Stage 2 covariance conditions on fixed V0(T); stage-1 uncertainty and shared cross-temperature correlations are not propagated. Coordinate errors, pressure-calibration systematics and phase/model uncertainty excluded. Covariance matrices are diagnostic, not author covariance.",
        "isotherms": isotherms,
        "author_workflow_reproduced": False,
        "qualifications": [
            "Replays the stated sequential fitting structure, not undisclosed author weights, objective, full-precision inputs or covariance. Unweighted pressure residuals are this explicit Peritheos baseline.",
            "Published curves use the original printed ambient law and tabulated individual-isotherm moduli. Refit curves use the freshly fitted ambient constant-alpha law, with no curved expansion or MGD replacement.",
            "All 11 Table 1 and 37 Table 4 observations retained. Printed errors shown in the plot but not used as regression weights.",
            "Four independent isotherm fits; no single continuous thermal EOS is inferred. Original phase and pressure-scale limitations persist. No catalog/default/source changes.",
        ],
    }


def write_plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rows = observations("fukui_2003", "table4")
    colors = {"published": "#c64c40", "staged_refit": "#187ca6"}
    fig = plt.figure(figsize=(13.2, 10.8), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.15, top=0.825, hspace=0.32, wspace=0.18
    )
    for i, result in enumerate(report["isotherms"]):
        t = result["temperature_k"]
        selected = [row for row in rows if float(row["temperature_k"]) == t]
        p, v, sp, sv = (
            np.array([float(row[name]) for row in selected])
            for name in [
                "pressure_gpa",
                "cell_volume_a3",
                "pressure_error_gpa",
                "cell_volume_error_a3",
            ]
        )
        inner = outer[i // 2, i % 2].subgridspec(
            2, 1, height_ratios=[2.5, 1], hspace=0.06
        )
        ax = fig.add_subplot(inner[0])
        residual_ax = fig.add_subplot(inner[1], sharex=ax)
        grid = np.linspace(p.min(), p.max(), 300)
        for key in ["published", "staged_refit"]:
            curve = result[key]
            parameters = [curve["V0_a3"], curve["K0_gpa"], curve["K0_prime"]]
            ax.plot(grid, comparison_volume(grid, *parameters), color=colors[key], lw=2)
            residual_ax.plot(
                p,
                curve["volume_residuals_a3"],
                "s" if key == "published" else "o",
                color=colors[key],
                ms=4,
                markerfacecolor="white" if key == "published" else colors[key],
            )
        ax.errorbar(
            p,
            v,
            xerr=sp,
            yerr=sv,
            fmt="o",
            color="#242424",
            ms=4,
            elinewidth=0.9,
            capsize=2,
            zorder=7,
        )
        residual_ax.errorbar(
            p,
            np.zeros_like(p),
            xerr=sp,
            yerr=sv,
            fmt="none",
            ecolor="#b8b8b8",
            elinewidth=0.8,
            capsize=2,
            zorder=0,
        )
        ax.set_title(f"{t} K · {len(p)} observations", loc="left", fontsize=12)
        ax.set_ylabel("Unit-cell volume (Å³)")
        ax.tick_params(labelbottom=False)
        ax.text(
            0.97,
            0.94,
            f"Volume RMS (Å³)\nPublished: {result['published']['volume_rmse_a3']:.3f}\nStaged refit: {result['staged_refit']['volume_rmse_a3']:.3f}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        residual_ax.axhline(0, color="#555555", lw=0.8)
        residual_ax.set_xlabel("Pressure (GPa)")
        residual_ax.set_ylabel("ΔV (Å³)")
        residual_ax.set_ylim(-0.55, 0.55)
        residual_ax.set_yticks([-0.4, 0, 0.4])
        for axis in [ax, residual_ax]:
            axis.grid(alpha=0.16)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Brucite: staged isotherm refits versus the published EOS", fontsize=17, y=0.975
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="o",
                color="#242424",
                ls="none",
                label="Published Table 4 observations + reported errors",
            ),
            Line2D(
                [],
                [],
                color=colors["published"],
                lw=2,
                label="Published expansion law + tabulated isotherm moduli",
            ),
            Line2D(
                [],
                [],
                color=colors["staged_refit"],
                lw=2,
                label="Staged refit: fitted ambient law fixed; K₀ and K′ fitted per isotherm",
            ),
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.945),
        frameon=False,
        fontsize=10,
    )
    ambient = report["stage1_ambient"]
    fig.text(
        0.075,
        0.067,
        f"Stage 1: unweighted fit to 11 ambient volumes gives V₃₀₀ = {ambient['V300_a3']:.6f} Å³ and constant α = {ambient['alpha_per_k']:.8g} K⁻¹.\nStage 2: hold V₀(T) fixed; independently fit K₀ and K′ using unweighted pressure residuals to all 13/8/8/8 Table 4 rows.\nΔV = calculated − observed volume at observed pressure. Grey bars show reported V/P errors about zero; model uncertainty is not shown.\nFour independent isotherms, not a continuous EOS. Original author objective and weights are unreported; published coefficients are preserved.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-staged-comparison"
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.check:
        from scripts.brucite_replay_checks import check_report

        check_report(report, json.loads(REPORT.read_text(encoding="utf-8")))
    else:
        REPORT.write_text(payload, encoding="utf-8")
    if args.plot:
        write_plot(report)
    print("Brucite staged comparison is current.")
    for result in report["isotherms"]:
        print(
            result["temperature_k"],
            json.dumps(
                {
                    key: {
                        field: result[key][field]
                        for field in [
                            "K0_gpa",
                            "K0_prime",
                            "pressure_rmse_gpa",
                            "volume_rmse_a3",
                        ]
                    }
                    for key in ["published", "staged_refit"]
                }
            ),
        )


if __name__ == "__main__":
    main()
