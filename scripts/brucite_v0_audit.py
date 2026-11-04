"""Infer conditional V0 for each Fukui isotherm with published moduli fixed.

Run with --plot to write the report and figures, or --check to check the
report. This diagnostic never modifies the material catalog or source data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product

import numpy as np
import scipy
from scipy.optimize import least_squares

from scripts.brucite_reproduce import (
    ROOT,
    bm3,
    bulk_modulus,
    comparison_volume,
    metrics,
    observations,
    reference_volume,
    source,
)

REPORT = ROOT / "docs/data/brucite-v0-audit.json"
TEMPERATURE_PLOT = ROOT / "docs/data/brucite-v0-temperature.png"
ISOTHERM_PLOT = ROOT / "docs/data/brucite-v0-isotherms.png"
MODES = ("unweighted_pressure", "pressure_sigma", "effective_sigma")


def fit_volume(rows, k0, kp, initial, mode):
    """Fit one reference volume; observational errors remain source metadata."""
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    sp = np.array([float(r["pressure_error_gpa"]) for r in rows])
    sv = np.array([float(r["cell_volume_error_a3"]) for r in rows])

    def residual(x):
        difference = bm3(v, x[0], k0, kp) - p
        if mode == "pressure_sigma":
            return difference / sp
        if mode == "effective_sigma":
            sigma = np.sqrt(sp**2 + (bulk_modulus(v, x[0], k0, kp) * sv / v) ** 2)
            return difference / sigma
        if mode != "unweighted_pressure":
            raise ValueError(mode)
        return difference

    fit = least_squares(
        residual,
        [initial],
        bounds=([35], [50]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
    )
    pressure_residual = bm3(v, fit.x[0], k0, kp) - p
    variance = float(
        (fit.fun @ fit.fun) / (len(rows) - 1) / (fit.jac.T @ fit.jac)[0, 0]
    )
    return {
        "V0_a3": float(fit.x[0]),
        "conditional_standard_error_a3": float(np.sqrt(variance)),
        "conditional_covariance_a6": [[variance]],
        "degrees_of_freedom": len(rows) - 1,
        "solver_success": bool(fit.success),
        "pressure_residuals_gpa": pressure_residual.tolist(),
        "mean_pressure_residual_gpa": float(np.mean(pressure_residual)),
        **metrics(pressure_residual),
    }


def expansion_volume(t, alpha0, alpha1=0, v300=40.917):
    """Integrate instantaneous alpha(T)=alpha0+alpha1*T from 300 K."""
    t = np.asarray(t)
    return v300 * np.exp(alpha0 * (t - 300) + 0.5 * alpha1 * (t**2 - 300**2))


def ambient_linear_alpha_fits():
    """Independent diagnostic using only Table 1, not inferred HP volumes."""
    rows = observations("fukui_2003", "table1")
    t = np.array([float(r["temperature_k"]) for r in rows])
    v = np.array([float(r["cell_volume_a3"]) for r in rows])
    sv = np.array([float(r["cell_volume_error_a3"]) for r in rows])
    dt = t - 300
    targets = np.array([300, 473, 673, 873]) - 300
    results = {}
    for mode in ("unweighted_volume", "volume_sigma"):
        weights = np.ones_like(sv) if mode == "unweighted_volume" else sv
        for free_anchor in (False, True):

            def volumes(x, delta):
                return (x[0] if free_anchor else 40.917) * np.exp(
                    x[-2] * 1e-5 * delta + 0.5 * x[-1] * 1e-8 * delta**2
                )

            fit = least_squares(
                lambda x: (volumes(x, dt) - v) / weights,
                [40.917, 6.69, 0] if free_anchor else [6.69, 0],
                xtol=1e-13,
                ftol=1e-13,
                gtol=1e-13,
            )
            scale = np.array([1, 1e-5, 1e-8] if free_anchor else [1e-5, 1e-8])
            names = (["V300_a3"] if free_anchor else []) + [
                "alpha_300_per_k",
                "dalpha_dT_per_k2",
            ]
            covariance = (
                np.linalg.inv(fit.jac.T @ fit.jac)
                * (fit.fun @ fit.fun)
                / (len(v) - len(fit.x))
            )
            covariance *= scale[:, None] * scale[None, :]
            predicted = volumes(fit.x, targets)
            derivatives = np.column_stack(
                [predicted * targets, 0.5 * predicted * targets**2]
            )
            if free_anchor:
                derivatives = np.column_stack([predicted / fit.x[0], derivatives])
            results[f"{mode}_{'free' if free_anchor else 'fixed'}_300k_anchor"] = {
                "V300_a3": float(fit.x[0]) if free_anchor else 40.917,
                "alpha_300_per_k": float(fit.x[-2] * 1e-5),
                "dalpha_dT_per_k2": float(fit.x[-1] * 1e-8),
                "alpha0_for_alpha0_plus_alpha1T_per_k": float(
                    fit.x[-2] * 1e-5 - 300 * fit.x[-1] * 1e-8
                ),
                "parameter_covariance": {
                    "parameter_order": names,
                    "matrix": covariance.tolist(),
                },
                "parameter_standard_errors": dict(
                    zip(names, np.sqrt(np.diag(covariance)).tolist())
                ),
                "source_rows": list(range(1, 12)),
                "volume_rmse_a3": float(
                    np.sqrt(np.mean((volumes(fit.x, dt) - v) ** 2))
                ),
                "V0_at_isotherms_a3": predicted.tolist(),
                "conditional_V0_standard_errors_at_isotherms_a3": np.sqrt(
                    np.einsum("ij,jk,ik->i", derivatives, covariance, derivatives)
                ).tolist(),
                "qualification": "New diagnostic alpha(T)=alpha(300)+dalpha/dT*(T-300), fitted only to the 11 ambient observations. Extrapolation to 873 K, not proof of the author's law. Residual-scaled local covariance omits systematics and model uncertainty.",
            }
    return results


def reproduce():
    src = source()["sources"]["fukui_2003"]
    rows = observations("fukui_2003", "table4")
    isotherms = []
    for parameters in src["parameters"]["isotherms"]:
        t = parameters["temperature_k"]
        group = [r for r in rows if float(r["temperature_k"]) == t]
        k0, kp = parameters["K0_gpa"], parameters["K0_prime"]
        fits = {
            mode: fit_volume(group, k0, kp, reference_volume(t), mode) for mode in MODES
        }
        v = np.array([float(r["cell_volume_a3"]) for r in group])
        p = np.array([float(r["pressure_gpa"]) for r in group])
        for result in fits.values():
            result["minus_printed_expansion_a3"] = float(
                result["V0_a3"] - reference_volume(t)
            )
            result["equivalent_constant_alpha_from_source_300k_anchor_per_k"] = (
                float(np.log(result["V0_a3"] / 40.917) / (t - 300))
                if t != 300
                else None
            )
        # The printed moduli have one decimal place. This checks rounding only,
        # separately from the much larger published coefficient uncertainties.
        corners = [
            fit_volume(group, k0 + dk, kp + dp, reference_volume(t), MODES[0])["V0_a3"]
            for dk, dp in product([-0.05, 0.05], repeat=2)
        ]
        isotherms.append(
            {
                "temperature_k": t,
                "observations": len(group),
                "source_rows": [int(r["source_row"]) for r in group],
                "fixed_published_moduli": {"K0_gpa": k0, "K0_prime": kp},
                "printed_expansion_V0_a3": float(reference_volume(t)),
                "printed_expansion_curve": metrics(
                    bm3(v, reference_volume(t), k0, kp) - p
                ),
                "fits": fits,
                "printed_rounding_corner_sensitivity": {
                    "K0_perturbation_gpa": [-0.05, 0.05],
                    "Kprime_perturbation": [-0.05, 0.05],
                    "unweighted_V0_range_a3": [min(corners), max(corners)],
                    "is_confidence_interval": False,
                },
            }
        )

    t = np.array([r["temperature_k"] for r in isotherms])
    comparisons = {}
    # Coefficients as printed in Fukui Table 3. A common V(300) isolates the
    # expansion laws; these are not reproductions of the other papers' EOSs.
    laws = {
        "fukui_2003_table3": (6.69e-5, 0),
        "fei_mao_1993_as_listed_by_fukui": (8.0e-5, 0),
        "xia_1998_as_listed_by_fukui": (7.3e-5, 3.6e-8),
        "redfern_wood_1992_as_listed_by_fukui": (0.23e-5, 2.37e-7),
    }
    for name, (alpha0, alpha1) in laws.items():
        predicted = expansion_volume(t, alpha0, alpha1)
        comparisons[name] = {
            "alpha0_per_k": alpha0,
            "alpha1_per_k2": alpha1,
            "common_V300_a3": 40.917,
            "V0_at_isotherms_a3": predicted.tolist(),
            "reference_volume_rmse_a3_by_mode": {
                mode: float(
                    np.sqrt(
                        np.mean(
                            (
                                predicted
                                - np.array(
                                    [r["fits"][mode]["V0_a3"] for r in isotherms]
                                )
                            )
                            ** 2
                        )
                    )
                )
                for mode in MODES
            },
        }

    trends = {}
    for mode in MODES:
        y = np.array([r["fits"][mode]["V0_a3"] for r in isotherms])
        for free_anchor in (False, True):

            def residual(x):
                return (
                    expansion_volume(
                        t, x[-1] * 1e-5, v300=x[0] if free_anchor else 40.917
                    )
                    - y
                )

            fit = least_squares(
                residual,
                [40.917, 7] if free_anchor else [7],
                xtol=1e-13,
                ftol=1e-13,
                gtol=1e-13,
            )
            trends[f"{mode}_{'free' if free_anchor else 'fixed'}_300k_anchor"] = {
                "V300_a3": float(fit.x[0]) if free_anchor else 40.917,
                "alpha_per_k": float(fit.x[-1] * 1e-5),
                "reference_volume_rmse_a3": float(np.sqrt(np.mean(fit.fun**2))),
                "reference_volume_residuals_a3": fit.fun.tolist(),
                "objective": "Unweighted volume residuals to four conditionally inferred V0 values, not a direct fit to ambient observations or original author workflow",
            }

    return {
        "audit_date": "2026-10-16",
        "software": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "optimizer": "scipy.optimize.least_squares; explicit bounds/tolerances in script",
        },
        "source": {
            "doi": src["doi"],
            "url": src["url"],
            "pdf_sha256": src["pdf_sha256"],
            "locations": "Table 4 observations, Tables 5-6 moduli, Eq. (1) and Tables 2-3 expansion",
        },
        "table4_csv_sha256": hashlib.sha256(
            (
                ROOT / "peritheos/data/datasets/brucite-fukui-2003-table4.csv"
            ).read_bytes()
        ).hexdigest(),
        "objectives": {
            "unweighted_pressure": "Minimize sum((P_BM3(V_observed;V0,K0_fixed,Kprime_fixed)-P_observed)^2), in GPa^2",
            "pressure_sigma": "Minimize sum((pressure residual / reported pressure standard deviation)^2)",
            "effective_sigma": "Minimize sum((pressure residual / sqrt(sigma_P^2+(K_T*sigma_V/V)^2))^2); parameter-dependent variance, diagnostic rather than confirmed author likelihood",
        },
        "covariance_convention": "SSE/(n-1) * inverse(J.T @ J), using the selected residual objective. Local linear conditional uncertainty only: published moduli treated as exact; their errors/covariance and pressure-scale systematics excluded. No author confidence level or workflow inferred.",
        "isotherms": isotherms,
        "temperature_law_comparisons": comparisons,
        "temperature_law_comparison_scope": "All Table 3 expansion laws compared with the same V300=40.917 A^3 anchor and exponential integral. Other authors' actual reference volumes/validity are not substituted; this is not validation of those EOSs. Extrapolation beyond ambient measurement intervals is diagnostic.",
        "constant_alpha_trend_fits": trends,
        "ambient_linear_alpha_diagnostics": ambient_linear_alpha_fits(),
        "conclusion": "The inferred reference volume near 42.7 A^3 at 873 K is robust across the three tested objectives and printed rounding corners. At 473 K the apparent unweighted offset is weighting-sensitive: effective-P/V weighting recovers the printed expansion volume closely. No single tested Table 3 expansion law consistently accounts for all inferred volumes. An independently fitted linear-temperature expansivity from the ambient observations extrapolates to about 42.68 A^3, a plausible explanation, but its curvature is weakly constrained and its 873 K prediction has about 0.15 A^3 conditional standard error. Actual author inputs/law cannot be uniquely recovered; no published parameter is changed.",
        "catalog_mutated": False,
    }


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = ["#176a9b", "#b7761c", "#22804a"]
    labels = [
        "Unweighted pressure fit",
        "Pressure-error weights",
        "Combined P/V-error weights",
    ]
    markers = ["o", "^", "s"]
    t = np.array([r["temperature_k"] for r in report["isotherms"]])
    grid_t = np.linspace(295, 890, 300)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.2, 7.0),
        sharex=True,
        height_ratios=[2, 1],
        layout="constrained",
    )
    fig.suptitle(
        "Brucite: reference volumes inferred with published moduli fixed", fontsize=12
    )
    ambient = observations("fukui_2003", "table1")
    axes[0].errorbar(
        [float(r["temperature_k"]) for r in ambient],
        [float(r["cell_volume_a3"]) for r in ambient],
        yerr=[float(r["cell_volume_error_a3"]) for r in ambient],
        fmt=".",
        color="0.65",
        label="Measured ambient volumes (Table 1)",
        capsize=2,
    )
    axes[0].plot(
        grid_t,
        reference_volume(grid_t),
        color="#b64b35",
        label="Printed law: α = 6.69 × 10⁻⁵ K⁻¹",
    )
    alpha = report["constant_alpha_trend_fits"][
        "unweighted_pressure_fixed_300k_anchor"
    ]["alpha_per_k"]
    axes[0].plot(
        grid_t,
        expansion_volume(grid_t, alpha),
        "--",
        color="0.35",
        label=f"Best constant α for blue points: {alpha / 1e-5:.3f} × 10⁻⁵ K⁻¹",
    )
    candidate = report["ambient_linear_alpha_diagnostics"][
        "unweighted_volume_free_300k_anchor"
    ]

    def ambient_candidate(temperature):
        dt = np.asarray(temperature) - 300
        return candidate["V300_a3"] * np.exp(
            candidate["alpha_300_per_k"] * dt
            + 0.5 * candidate["dalpha_dT_per_k2"] * dt**2
        )

    predicted = ambient_candidate(grid_t)
    covariance = np.array(candidate["parameter_covariance"]["matrix"])
    dt = grid_t - 300
    derivative = np.column_stack(
        [predicted / candidate["V300_a3"], predicted * dt, 0.5 * predicted * dt**2]
    )
    sigma = np.sqrt(np.einsum("ij,jk,ik->i", derivative, covariance, derivative))
    axes[0].fill_between(
        grid_t, predicted - sigma, predicted + sigma, color="#8156a8", alpha=0.10
    )
    axes[0].plot(
        grid_t,
        predicted,
        color="#8156a8",
        label="Diagnostic α(T) fit to ambient data (±1 SE)",
    )
    axes[1].plot(grid_t, predicted - reference_volume(grid_t), color="#8156a8", lw=1)
    for mode, color, label, marker in zip(MODES, colors, labels, markers):
        y = np.array([r["fits"][mode]["V0_a3"] for r in report["isotherms"]])
        errors = [
            r["fits"][mode]["conditional_standard_error_a3"]
            for r in report["isotherms"]
        ]
        axes[0].errorbar(
            t, y, yerr=errors, fmt=marker, color=color, capsize=3, label=label, zorder=3
        )
        axes[1].errorbar(
            t,
            y - reference_volume(t),
            yerr=errors,
            fmt=marker + "--",
            color=color,
            capsize=3,
        )
    axes[0].set(ylabel="Zero-pressure reference volume (Å³)")
    axes[0].legend(fontsize=8, loc="upper left")
    axes[1].axhline(0, color="#b64b35", lw=1)
    axes[1].set(
        xlabel="Temperature (K)", ylabel="Inferred V₀ − printed-law V₀\n(Å³)", xticks=t
    )
    for ax in axes:
        ax.grid(alpha=0.15)
    fig.savefig(TEMPERATURE_PLOT, dpi=170, bbox_inches="tight")
    plt.close(fig)

    rows = observations("fukui_2003", "table4")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), layout="constrained")
    for ax, result in zip(axes.flat, report["isotherms"]):
        temp = result["temperature_k"]
        group = [r for r in rows if float(r["temperature_k"]) == temp]
        v = np.array([float(r["cell_volume_a3"]) for r in group])
        p = np.array([float(r["pressure_gpa"]) for r in group])
        ax.errorbar(
            p,
            v,
            xerr=[float(r["pressure_error_gpa"]) for r in group],
            yerr=[float(r["cell_volume_error_a3"]) for r in group],
            fmt="o",
            markersize=3,
            color="black",
            capsize=2,
            label="Table 4 observations",
        )
        grid = np.linspace(p.min(), p.max(), 300)
        fixed = result["fixed_published_moduli"]
        ax.plot(
            grid,
            comparison_volume(
                grid,
                result["printed_expansion_V0_a3"],
                fixed["K0_gpa"],
                fixed["K0_prime"],
            ),
            color="#b64b35",
            label="Printed expansion V₀",
        )
        for mode, color, style, label in [
            (MODES[0], colors[0], "-", "Fitted V₀, unweighted"),
            (MODES[2], colors[2], "--", "Fitted V₀, combined errors"),
        ]:
            ax.plot(
                grid,
                comparison_volume(
                    grid,
                    result["fits"][mode]["V0_a3"],
                    fixed["K0_gpa"],
                    fixed["K0_prime"],
                ),
                style,
                color=color,
                label=label,
            )
        fit = result["fits"][MODES[0]]
        ax.set(
            title=f"{temp} K: K₀ = {fixed['K0_gpa']:g} GPa, K′ = {fixed['K0_prime']:g} fixed",
            xlabel="Pressure (GPa)",
            ylabel="Cell volume (Å³)",
        )
        ax.text(
            0.04,
            0.05,
            f"V₀: {result['printed_expansion_V0_a3']:.3f} → {fit['V0_a3']:.3f} Å³\nRMS: {result['printed_expansion_curve']['pressure_rmse_gpa']:.3f} → {fit['pressure_rmse_gpa']:.3f} GPa",
            transform=ax.transAxes,
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        ax.grid(alpha=0.15)
    axes[0, 0].legend(fontsize=7, loc="upper right")
    fig.suptitle(
        "All four brucite isotherms: fit only V₀, retain published moduli", fontsize=13
    )
    fig.savefig(ISOTHERM_PLOT, dpi=170, bbox_inches="tight")
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
        print("Brucite V0 audit report is current.")
    else:
        REPORT.write_text(payload, encoding="utf-8")
    if args.plot:
        plot(report)
    for r in report["isotherms"]:
        result = r["fits"][MODES[0]]
        print(
            f"{r['temperature_k']} K: printed V0 {r['printed_expansion_V0_a3']:.6f}; fitted V0 {result['V0_a3']:.6f} +/- {result['conditional_standard_error_a3']:.6f} A^3"
        )


if __name__ == "__main__":
    main()
