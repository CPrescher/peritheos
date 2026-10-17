"""Fit linear temperature laws to the staged brucite moduli, not raw PVT.

The fitted constant-alpha ambient law is held exactly fixed. Equal-weight
linear regressions have free intercepts and slopes. No catalog mutation.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
from scipy.linalg import block_diag

from scripts.brucite_reproduce import (
    ROOT,
    bm3,
    comparison_volume,
    metrics,
    observations,
)
from scripts.brucite_staged_comparison import REPORT as STAGED_REPORT
from scripts.brucite_staged_comparison import reproduce as staged_reproduce

REPORT = ROOT / "docs/data/brucite-staged-continuous.json"
ORDER = ["K300_gpa", "dK_dT_gpa_per_k", "Kprime300", "dKprime_dT_per_k"]


def moduli(temperature, parameters):
    dt = np.asarray(temperature) - 300
    return (
        parameters["K300_gpa"] + parameters["dK_dT_gpa_per_k"] * dt,
        parameters["Kprime300"] + parameters["dKprime_dT_per_k"] * dt,
    )


def reference_volume(temperature, ambient):
    return ambient["V300_a3"] * np.exp(
        ambient["alpha_per_k"] * (np.asarray(temperature) - 300)
    )


def runtime_model(report):
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import ThermalReferenceStateEOS

    ambient, parameters = (
        report["fixed_ambient_law"],
        report["linear_moduli"]["parameters"],
    )
    return ThermalReferenceStateEOS(
        BM3(ambient["V300_a3"], parameters["K300_gpa"], parameters["Kprime300"]),
        Tr=300,
        alpha0=ambient["alpha_per_k"],
        dK_dT=parameters["dK_dT_gpa_per_k"],
        dKprime_dT=parameters["dKprime_dT_per_k"],
    )


def reproduce():
    staged = staged_reproduce()
    if json.loads(STAGED_REPORT.read_text()) != staged:
        raise ValueError("Staged-input report is stale")
    temperatures = np.array([s["temperature_k"] for s in staged["isotherms"]])
    design = np.column_stack([np.ones(4), temperatures - 300])
    values = np.array(
        [
            [s["staged_refit"]["K0_gpa"], s["staged_refit"]["K0_prime"]]
            for s in staged["isotherms"]
        ]
    )
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    parameters = dict(zip(ORDER, coefficients.T.ravel().tolist()))
    residual = design @ coefficients - values
    residual_covariance = residual.T @ residual / (4 - 2)
    scatter_covariance = np.kron(residual_covariance, np.linalg.inv(design.T @ design))
    # Propagate the existing per-isotherm K/Kprime covariance separately from
    # the scatter-based regression convention. The two are not added together.
    inverse_design = np.linalg.solve(design.T @ design, design.T)
    transform = np.zeros((4, 8))
    transform[:2, 0::2] = inverse_design
    transform[2:, 1::2] = inverse_design
    input_covariance = block_diag(
        *[s["conditional_covariance_fixed_V0"]["matrix"] for s in staged["isotherms"]]
    )
    propagated_covariance = transform @ input_covariance @ transform.T
    ambient = staged["stage1_ambient"]
    rows = observations("fukui_2003", "table4")
    curves = []
    all_residuals = {
        key: {"pressure": [], "volume": []}
        for key in ["published", "staged_refit", "continuous"]
    }
    for s in staged["isotherms"]:
        t = s["temperature_k"]
        selected = [row for row in rows if float(row["temperature_k"]) == t]
        v = np.array([float(row["cell_volume_a3"]) for row in selected])
        p = np.array([float(row["pressure_gpa"]) for row in selected])
        k0, kp = moduli(t, parameters)
        v0 = float(reference_volume(t, ambient))
        dp = bm3(v, v0, k0, kp) - p
        dv = comparison_volume(p, v0, k0, kp) - v
        continuous = {
            "V0_a3": v0,
            "K0_gpa": float(k0),
            "K0_prime": float(kp),
            "pressure_residuals_gpa": dp.tolist(),
            "volume_residuals_a3": dv.tolist(),
            **metrics(dp),
            "volume_rmse_a3": float(np.sqrt(np.mean(dv**2))),
        }
        result = {
            "temperature_k": t,
            "observations": len(selected),
            "source_rows": s["source_rows"],
            "pressure_range_gpa": s["pressure_range_gpa"],
            "published": s["published"],
            "staged_refit": s["staged_refit"],
            "continuous": continuous,
        }
        curves.append(result)
        for key in all_residuals:
            all_residuals[key]["pressure"].extend(result[key]["pressure_residuals_gpa"])
            all_residuals[key]["volume"].extend(result[key]["volume_residuals_a3"])
    report = {
        "format": "peritheos.brucite-staged-continuous",
        "format_version": 1,
        "staged_report_sha256": hashlib.sha256(STAGED_REPORT.read_bytes()).hexdigest(),
        "source_sha256": staged["source_sha256"],
        "fixed_ambient_law": ambient,
        "linear_moduli": {
            "parameters": parameters,
            "reference_temperature_k": 300,
            "input_temperatures_k": temperatures.tolist(),
            "input_K_Kprime": values.tolist(),
            "predicted_K_Kprime": (design @ coefficients).tolist(),
            "parameter_order": ORDER,
            "objective": "Separate equal-weight least-squares lines K(T)=K300+dK_dT*(T-300) and Kprime(T)=Kprime300+dKprime_dT*(T-300) through the four staged modulus pairs. Both intercepts and slopes free; no 300 K modulus pinned and no optimization against raw PVT residuals in this stage.",
            "K_Kprime_residuals": residual.tolist(),
            "covariance_from_line_scatter": scatter_covariance.tolist(),
            "standard_errors_from_line_scatter": dict(
                zip(ORDER, np.sqrt(np.diag(scatter_covariance)).tolist())
            ),
            "covariance_propagated_from_stage2": propagated_covariance.tolist(),
            "covariance_qualification": "Line-scatter convention: residual.T@residual/(4-2) Kronecker inverse(X.T@X), order K300,dK_dT,Kprime300,dKprime_dT. Separate stage-input propagation: A*blockdiag(Cov(K,Kprime)_T)*A.T for the exact OLS estimator. These two covariance estimates describe different assumptions and are not added. Both condition on the fixed ambient law; shared ambient uncertainty, predictor errors, pressure calibration, phase/model discrepancy and inter-isotherm covariance excluded. No author covariance recovered.",
        },
        "isotherms": curves,
        "all37_assessment": {
            key: {
                **metrics(values["pressure"]),
                "volume_rmse_a3": float(
                    np.sqrt(np.mean(np.array(values["volume"]) ** 2))
                ),
            }
            for key, values in all_residuals.items()
        },
        "qualifications": [
            "Follows the requested stage: fit temperature dependence of the staged moduli, then evaluate the resulting continuous BM3. This is not a direct joint optimization of the 37 PVT rows.",
            "Fixed ambient constant-alpha law is copied exactly from the prior staged refit. All 11 ambient-volume predictions remain unchanged.",
            "Continuous temperature representation over 300-873 K, assessed at the four measured temperatures; no physical validation between them or outside experimental pressure coverage claimed.",
            "Linear trends need not pass through any individual staged modulus pair. The four slices are now linked, so their residuals can grow. Per-isotherm K/Kprime errors/correlation are preserved but not used as regression weights.",
            "The plotted published comparator uses original separate-isotherm parameters and printed expansion, not an inferred original continuous surface. Original source weights/author workflow remain unrecovered.",
            "Diagnostic only. No catalog, source, original record or default changed; phase and calibration qualifications persist.",
        ],
    }
    runtime = runtime_model(report)
    v = np.array([float(row["cell_volume_a3"]) for row in rows])
    t = np.array([float(row["temperature_k"]) for row in rows])
    k0, kp = moduli(t, parameters)
    expected = bm3(v, reference_volume(t, ambient), k0, kp)
    difference = float(np.max(abs(runtime.pressure(v, t) - expected)))
    if difference > 1e-10:
        raise ValueError("Continuous model disagrees with Peritheos runtime")
    report["runtime_max_pressure_difference_gpa"] = difference
    return report


def write_plots(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {
        "published": "#c64c40",
        "staged_refit": "#8a8f93",
        "continuous": "#187ca6",
    }
    styles = {"published": "-", "staged_refit": "--", "continuous": "-"}
    rows = observations("fukui_2003", "table4")
    fig = plt.figure(figsize=(13.2, 11.2), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.15, top=0.80, hspace=0.32, wspace=0.18
    )
    for i, s in enumerate(report["isotherms"]):
        t = s["temperature_k"]
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
        res = fig.add_subplot(inner[1], sharex=ax)
        grid = np.linspace(p.min(), p.max(), 250)
        for key in ["published", "staged_refit", "continuous"]:
            curve = s[key]
            ax.plot(
                grid,
                comparison_volume(
                    grid, curve["V0_a3"], curve["K0_gpa"], curve["K0_prime"]
                ),
                color=colors[key],
                ls=styles[key],
                lw=1.8,
            )
            res.plot(
                p,
                curve["volume_residuals_a3"],
                {"published": "s", "staged_refit": "^", "continuous": "o"}[key],
                color=colors[key],
                ms=4,
                markerfacecolor=colors[key] if key == "continuous" else "white",
            )
        ax.errorbar(
            p,
            v,
            xerr=sp,
            yerr=sv,
            fmt="o",
            color="#242424",
            ms=4,
            capsize=2,
            elinewidth=0.8,
            zorder=8,
        )
        res.errorbar(
            p,
            np.zeros_like(p),
            xerr=sp,
            yerr=sv,
            fmt="none",
            ecolor="#b8b8b8",
            capsize=2,
            elinewidth=0.8,
            zorder=0,
        )
        ax.set_title(f"{t} K · {len(p)} observations", loc="left", fontsize=12)
        ax.set_ylabel("Unit-cell volume (Å³)")
        ax.tick_params(labelbottom=False)
        ax.text(
            0.97,
            0.94,
            f"Volume RMS (Å³)\nPublished: {s['published']['volume_rmse_a3']:.3f}\nSeparate staged: {s['staged_refit']['volume_rmse_a3']:.3f}\nContinuous: {s['continuous']['volume_rmse_a3']:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        res.axhline(0, color="#555555", lw=0.8)
        res.set_xlabel("Pressure (GPa)")
        res.set_ylabel("ΔV (Å³)")
        res.set_ylim(-0.55, 0.55)
        res.set_yticks([-0.4, 0, 0.4])
        for axis in [ax, res]:
            axis.grid(alpha=0.16)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Brucite: continuous EOS from linear fits to the staged moduli",
        y=0.975,
        fontsize=17,
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
            *[
                Line2D([], [], color=colors[key], ls=styles[key], lw=2, label=label)
                for key, label in [
                    (
                        "continuous",
                        "Continuous refit: linear K₀(T) and K′(T); ambient expansion fixed",
                    ),
                    (
                        "staged_refit",
                        "Separate staged isotherm refits (before fitting temperature trends)",
                    ),
                    (
                        "published",
                        "Published individual isotherms: printed expansion + tabulated moduli",
                    ),
                ]
            ],
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.945),
        frameon=False,
        fontsize=10,
    )
    ambient = report["fixed_ambient_law"]
    pars = report["linear_moduli"]["parameters"]
    fig.text(
        0.075,
        0.065,
        f"Fixed V₀(T) = {ambient['V300_a3']:.6f} exp[{ambient['alpha_per_k']:.8g}(T − 300)] Å³; ambient-only fit retained exactly.\nK₀(T) = {pars['K300_gpa']:.5f} {pars['dK_dT_gpa_per_k']:+.8f}(T − 300) GPa; K′(T) = {pars['Kprime300']:.5f} {pars['dKprime_dT_per_k']:+.8f}(T − 300).\nBoth intercepts/slopes fitted with equal weights to the four staged modulus pairs; raw P–V–T residuals were not reoptimized.\nΔV = calculated − observed volume at observed pressure; grey bars show reported V/P errors about zero. Model uncertainty is not shown.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-staged-continuous"
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), layout="constrained")
    temperatures = np.array(report["linear_moduli"]["input_temperatures_k"])
    grid = np.linspace(300, 873, 200)
    line_values = moduli(grid, pars)
    staged = json.loads(STAGED_REPORT.read_text())
    for j, ax in enumerate(axes):
        values = np.array(report["linear_moduli"]["input_K_Kprime"])[:, j]
        sigma = [
            np.sqrt(s["conditional_covariance_fixed_V0"]["matrix"][j][j])
            for s in staged["isotherms"]
        ]
        ax.errorbar(
            temperatures,
            values,
            yerr=sigma,
            fmt="o",
            color="#242424",
            capsize=3,
            label="Staged moduli + conditional local errors",
        )
        ax.plot(
            grid,
            line_values[j],
            color=colors["continuous"],
            lw=2,
            label="Unweighted linear trend",
        )
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("K₀ (GPa)" if j == 0 else "K′ (dimensionless)")
        ax.grid(alpha=0.18)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8)
    fig.suptitle("Linear temperature laws fitted to the four staged modulus pairs")
    stem = ROOT / "docs/data/brucite-staged-moduli"
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
        if REPORT.read_text() != payload:
            raise SystemExit("Brucite staged continuous comparison is stale")
    else:
        REPORT.write_text(payload)
    if args.plot:
        write_plots(report)
    print("Brucite staged continuous comparison is current.")
    print(json.dumps(report["linear_moduli"]["parameters"], indent=2))
    print(json.dumps(report["all37_assessment"], indent=2))


if __name__ == "__main__":
    main()
