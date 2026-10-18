"""Independent diagnostic BM3-MGD fits to Fukui's brucite observations.

Run as ``python -m scripts.brucite_mgd_comparison``. Source data and catalog
records are never modified. Both pressure and volume objectives are unweighted.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.constants import Avogadro, R
from scipy.optimize import least_squares

from scripts.brucite_curved_expansion import REPORT as CURVED_REPORT
from scripts.brucite_reproduce import (
    ROOT,
    bm3,
    comparison_volume,
    observations,
    reference_volume,
    source,
)

REPORT = ROOT / "docs/data/brucite-mgd-comparison.json"
VOLUME_TO_J_PER_BAR_MOL = Avogadro * 1e-25  # Å³ per Z=1 cell -> J/bar/mol
N_ATOMS = 5  # Mg(OH)2: Mg + 2 O + 2 H; one formula unit per cell
NODES, WEIGHTS = leggauss(64)
NODES, WEIGHTS = (NODES + 1) / 2, WEIGHTS / 2
NAMES = {
    "mgd": ["V0_a3", "K0_gpa", "Kprime", "gamma0", "q", "theta0_k"],
    "empirical": [
        "V0_a3",
        "K0_gpa",
        "Kprime",
        "alpha300_per_k",
        "dalpha_dT_per_k2",
        "K873_gpa",
        "Kprime873",
    ],
}
SCALES = {
    "mgd": np.array([1, 1, 1, 1, 1, 1000]),
    "empirical": np.array([1, 1, 1, 1e-5, 1e-8, 1, 1]),
}
BOUNDS = {
    "mgd": ([40, 20, 3, 0.01, -5, 0.01], [42, 70, 12, 3, 10, 4]),
    "empirical": ([40, 20, 3, 0, -20, 5, 3], [42, 70, 12, 20, 30, 70, 15]),
}


def debye3(x):
    """Independent Gauss-Legendre integral of 3/x³ int t³/(exp(t)-1) dt."""
    x = np.asarray(x, dtype=float)
    argument = x[..., None] * NODES
    return 3 * np.sum(WEIGHTS * NODES**2 * argument / np.expm1(argument), axis=-1)


def pressure(model, volume, temperature, parameters):
    v, t = np.broadcast_arrays(np.asarray(volume), np.asarray(temperature))
    if model == "mgd":
        v0, k0, kp, gamma0, q, theta0 = parameters
        log_ratio = np.log(v / v0)
        gamma = gamma0 * np.exp(q * log_ratio)
        exponent = -gamma0 * (
            np.expm1(q * log_ratio) / q if abs(q) > 1e-10 else log_ratio
        )
        theta = theta0 * np.exp(exponent)
        energy_difference = (
            3 * N_ATOMS * R * (t * debye3(theta / t) - 300 * debye3(theta / 300))
        )
        return bm3(v, v0, k0, kp) + 1e-4 * gamma * energy_difference / (
            v * VOLUME_TO_J_PER_BAR_MOL
        )
    if model == "empirical":
        v0, k0, kp, alpha, beta, k873, kp873 = parameters
        dt = t - 300
        vref = v0 * np.exp(alpha * dt + 0.5 * beta * dt**2)
        return bm3(v, vref, k0 + (k873 - k0) * dt / 573, kp + (kp873 - kp) * dt / 573)
    raise ValueError(model)


def volumes(model, target_pressure, temperature, parameters, initial):
    """Vector Newton inversion near the observed, positive-modulus branch."""
    result = np.array(initial, dtype=float, copy=True)
    for _ in range(45):
        residual = pressure(model, result, temperature, parameters) - target_pressure
        if np.max(abs(residual)) < 2e-11:
            return result
        step = result * 1e-5
        derivative = (
            pressure(model, result + step, temperature, parameters)
            - pressure(model, result - step, temperature, parameters)
        ) / (2 * step)
        result -= np.clip(residual / derivative, -0.6, 0.6)
        result = np.clip(result, 26, 48)
    raise RuntimeError(f"{model}: volume inversion failed")


def data():
    table4 = observations("fukui_2003", "table4")
    table1 = observations("fukui_2003", "table1")
    # Table 1 pressures are ambient conditions, not measured pressure entries.
    rows = [("table4", row) for row in table4] + [("table1", row) for row in table1]
    return {
        "temperature": np.array([float(row["temperature_k"]) for _, row in rows]),
        "volume": np.array([float(row["cell_volume_a3"]) for _, row in rows]),
        "pressure": np.array(
            [
                float(row["pressure_gpa"]) if name == "table4" else 0.0
                for name, row in rows
            ]
        ),
        "rows": [
            {"table": name, "source_row": int(row["source_row"])} for name, row in rows
        ],
        "volume_error": np.array(
            [float(row["cell_volume_error_a3"]) for _, row in rows]
        ),
        "pressure_error": np.array(
            [
                float(row["pressure_error_gpa"]) if name == "table4" else np.nan
                for name, row in rows
            ]
        ),
    }


def initial_values(model):
    if model == "mgd":
        return [
            np.array([40.93, 41.8, 6.6, gamma, q, theta])
            for gamma, q, theta in [
                (0.6, 1, 0.3),
                (0.8, 3, 0.8),
                (1.2, 5, 1.6),
                (1.5, 0, 3.0),
            ]
        ]
    return [np.array([40.93, 41.8, 6.6, 6.1, 4.4, 23.5, 10])]


def fit_model(model, objective, dataset, starts=None, fixed_theta=None):
    scales = SCALES[model]
    lower, upper = (np.array(v, dtype=float) for v in BOUNDS[model])
    active = np.arange(len(scales))
    template = initial_values(model)[0].copy()
    if fixed_theta is not None:
        assert model == "mgd"
        active = active[:-1]
        template[-1] = fixed_theta / 1000

    def expand(x):
        result = template.copy()
        result[active] = x
        return result * scales

    def residual(x):
        parameters = expand(x)
        if objective == "pressure":
            return (
                pressure(model, dataset["volume"], dataset["temperature"], parameters)
                - dataset["pressure"]
            )
        if objective == "volume":
            return (
                volumes(
                    model,
                    dataset["pressure"],
                    dataset["temperature"],
                    parameters,
                    dataset["volume"],
                )
                - dataset["volume"]
            )
        raise ValueError(objective)

    results = []
    for start in initial_values(model) if starts is None else starts:
        result = least_squares(
            residual,
            start[active],
            bounds=(lower[active], upper[active]),
            x_scale="jac",
            ftol=2e-11,
            xtol=2e-11,
            gtol=2e-11,
            max_nfev=1200,
        )
        results.append(result)
    best = min(results, key=lambda result: result.fun @ result.fun)
    parameters = expand(best.x)
    rank = int(np.linalg.matrix_rank(best.jac))
    covariance = (
        np.linalg.pinv(best.jac.T @ best.jac)
        * (best.fun @ best.fun)
        / (len(best.fun) - len(active))
    )
    covariance *= scales[active, None] * scales[None, active]
    return {
        "model": model,
        "objective": objective,
        "parameters": dict(zip(NAMES[model], parameters.tolist())),
        "parameter_order": [NAMES[model][i] for i in active],
        "conditional_covariance": covariance.tolist(),
        "conditional_standard_errors": dict(
            zip(
                [NAMES[model][i] for i in active], np.sqrt(np.diag(covariance)).tolist()
            )
        ),
        "covariance_qualification": "Local SSE/(N-p) * inverse(J.T J); all row residuals unweighted. Ignores both-coordinate errors, pressure-scale systematics, model error and phase changes. Boundary optima invalidate symmetric confidence intervals.",
        "solver_success": bool(best.success),
        "solver_message": best.message,
        "jacobian_rank": rank,
        "jacobian_condition_number_scaled": float(np.linalg.cond(best.jac)),
        "active_bounds": {
            NAMES[model][i]: int(mask)
            for i, mask in zip(active, best.active_mask)
            if mask
        },
        "bounds": {
            name: [float(lo * scale), float(hi * scale)]
            for name, lo, hi, scale in zip(NAMES[model], lower, upper, scales)
        },
        "fixed_theta0_k": fixed_theta,
        "objective_sse": float(best.fun @ best.fun),
        "multistart_sse": [float(result.fun @ result.fun) for result in results],
        "observations": len(best.fun),
        "free_parameters": len(active),
        "row_pressure_residuals_gpa": (
            pressure(model, dataset["volume"], dataset["temperature"], parameters)
            - dataset["pressure"]
        ).tolist(),
        "row_volume_residuals_a3": (
            volumes(
                model,
                dataset["pressure"],
                dataset["temperature"],
                parameters,
                dataset["volume"],
            )
            - dataset["volume"]
        ).tolist(),
    }


def assessment(pressure_residual, volume_residual, dataset):
    """Use the same rows and denominators for every comparison."""
    groups = {
        "table4_all": np.arange(48) < 37,
        "table1_ambient": np.arange(48) >= 37,
        "all48": np.ones(48, dtype=bool),
    }
    groups.update(
        {
            f"table4_{t}k": (np.arange(48) < 37) & (dataset["temperature"] == t)
            for t in [300, 473, 673, 873]
        }
    )
    return {
        key: {
            "observations": int(mask.sum()),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(pressure_residual[mask] ** 2))),
            "volume_rmse_a3": float(np.sqrt(np.mean(volume_residual[mask] ** 2))),
            "mean_volume_residual_a3": float(np.mean(volume_residual[mask])),
            "maximum_absolute_volume_residual_a3": float(
                np.max(abs(volume_residual[mask]))
            ),
        }
        for key, mask in groups.items()
    }


def check_implementation(parameters):
    """Compare the independent quadrature with native and Python MGD paths."""
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import MieGruneisenDebye

    v0, k0, kp, gamma, q, theta = parameters
    runtime = MieGruneisenDebye(
        BM3(v0 * VOLUME_TO_J_PER_BAR_MOL, k0, kp),
        Tr=300,
        theta0=theta,
        gamma0=gamma,
        q=q,
        n=N_ATOMS,
    )
    v, t = np.meshgrid(np.linspace(32, 43, 9), [299.3, 300, 473, 592.9, 673, 873])
    expected = pressure("mgd", v, t, parameters)
    native_difference = float(
        np.max(abs(runtime.pressure(v * VOLUME_TO_J_PER_BAR_MOL, t) - expected))
    )
    del runtime._native
    python_difference = float(
        np.max(abs(runtime.pressure(v * VOLUME_TO_J_PER_BAR_MOL, t) - expected))
    )
    np.testing.assert_allclose(native_difference, 0, atol=1e-9)
    np.testing.assert_allclose(python_difference, 0, atol=1e-9)
    step = v * 1e-5
    modulus = (
        -v
        * (
            pressure("mgd", v + step, t, parameters)
            - pressure("mgd", v - step, t, parameters)
        )
        / (2 * step)
    )
    if np.any(modulus <= 0):
        raise ValueError("MGD comparison includes a nonpositive bulk modulus")
    return {
        "states": v.size,
        "native_max_difference_gpa": native_difference,
        "python_max_difference_gpa": python_difference,
        "minimum_bulk_modulus_gpa": float(modulus.min()),
    }


def reproduce():
    dataset = data()
    results = {}
    for model in ("mgd", "empirical"):
        pfit = fit_model(model, "pressure", dataset)
        start = np.array(list(pfit["parameters"].values())) / SCALES[model]
        vfit = fit_model(
            model, "volume", dataset, starts=[start, *initial_values(model)]
        )
        for fit in (pfit, vfit):
            key = f"{model}_{fit['objective']}"
            fit["assessment"] = assessment(
                np.array(fit["row_pressure_residuals_gpa"]),
                np.array(fit["row_volume_residuals_a3"]),
                dataset,
            )
            results[key] = fit
    curved = json.loads(CURVED_REPORT.read_text(encoding="utf-8"))["fit"]
    old_parameters = np.array(
        [
            curved["V300_a3"],
            41.8,
            6.6,
            curved["alpha_300_per_k"],
            curved["dalpha_dT_per_k2"],
            41.8 - 0.0319 * 573,
            6.6 + 0.00571 * 573,
        ]
    )
    old_residual_p = (
        pressure("empirical", dataset["volume"], dataset["temperature"], old_parameters)
        - dataset["pressure"]
    )
    old_residual_v = (
        volumes(
            "empirical",
            dataset["pressure"],
            dataset["temperature"],
            old_parameters,
            dataset["volume"],
        )
        - dataset["volume"]
    )
    profiles = {}
    implementation = {}
    for objective in ("pressure", "volume"):
        fit = results[f"mgd_{objective}"]
        initial = np.array(list(fit["parameters"].values())) / SCALES["mgd"]
        profiles[objective] = []
        for theta in [10, 300, 600, 900, 1200, 1800, 2400, 3000, 4000]:
            profile = fit_model(
                "mgd", objective, dataset, starts=[initial], fixed_theta=theta
            )
            profiles[objective].append(
                {
                    "theta0_k": theta,
                    "objective_sse": profile["objective_sse"],
                    "parameters": profile["parameters"],
                    "solver_success": profile["solver_success"],
                    "active_bounds": profile["active_bounds"],
                }
            )
        implementation[objective] = check_implementation(
            list(fit["parameters"].values())
        )
    return {
        "format": "peritheos.brucite-mgd-comparison",
        "format_version": 1,
        "primary_source": source()["sources"]["fukui_2003"]["doi"],
        "source_file_sha256": {
            name: hashlib.sha256(
                (ROOT / "peritheos/data/datasets" / name).read_bytes()
            ).hexdigest()
            for name in [
                "brucite-source.json",
                "brucite-fukui-2003-table1.csv",
                "brucite-fukui-2003-table4.csv",
            ]
        },
        "software": {"numpy": np.__version__, "scipy": scipy.__version__},
        "selection": "All 37 Fukui Table 4 P-V-T observations and all 11 Table 1 ambient T-V observations. Table 1 P=0 GPa is the ambient-condition assumption, not a printed pressure measurement. No rows excluded, no Xia dataset mixed in, no pressure recalibration.",
        "row_order": dataset["rows"],
        "primary_objective": "volume",
        "objectives": {
            "pressure": "Unweighted sum over all 48 rows of [Pmodel(Vobs,Tobs)-Pobs]^2 [GPa^2]. Explicit baseline; V and T treated as exact predictors.",
            "volume": "Unweighted sum over all 48 rows of [Vmodel(Pobs,Tobs)-Vobs]^2 [A^6]. P and T treated as exact predictors. Invert on the local branch near observations; all fitted MGD models checked for positive bulk modulus.",
        },
        "mgd_definition": {
            "reference_temperature_k": 300,
            "thermal_pressure_reference": "reference_temperature",
            "debye_temperature_law": "integrated_gruneisen",
            "n_atoms_per_formula_unit": N_ATOMS,
            "formula_units_per_cell": 1,
            "Cvmax_j_per_mol_k": 3 * N_ATOMS * R,
            "cell_a3_to_j_per_bar_mol": VOLUME_TO_J_PER_BAR_MOL,
            "free_parameters": NAMES["mgd"],
            "equation": "P=BM3(V;V0,K0,Kprime)+1e-4*gamma(V)*[ED(V,T)-ED(V,300)]/Vmolar; gamma=gamma0*(V/V0)^q; theta=theta0*exp[-gamma0*((V/V0)^q-1)/q], with continuous q=0 limit; ED=3*n*R*T*D3(theta/T). Effective single-Debye spectrum, all six parameters fitted, no adopted literature theta.",
        },
        "empirical_definition": "Same curved-expansion BM3 family as the current catalog variant, now jointly fitted to exactly the same 48 rows. Seven free parameters: V300,K300,Kprime300,alpha300,beta,K873,Kprime873. K and Kprime vary linearly between endpoints; V0(T)=V300*exp(alpha300*dT+0.5*beta*dT^2). Endpoint parameterization guarantees positive reference K over 300-873 K.",
        "current_continuous": {
            "identifier": "brucite_peritheos_fukui_2003_curved_expansion_thermal_refit",
            "parameters": dict(zip(NAMES["empirical"], old_parameters.tolist())),
            "fitted_observations": 11,
            "free_parameters": 3,
            "qualification": "Only ambient expansion was fitted; published moduli and slopes adopted. Comparison with joint fits includes the effect of refitting compression parameters, not just changing EOS family.",
            "assessment": assessment(old_residual_p, old_residual_v, dataset),
        },
        "fits": results,
        "theta_sensitivity_profiles": profiles,
        "independent_implementation_check": implementation,
        "interpretation_limits": [
            "These are in-sample residual comparisons; no external prediction or pressure-standard validation.",
            "Current model used only 11 ambient fit rows; both new families use the same 48 joint rows. Matched empirical control has seven free parameters; MGD has six.",
            "Single-Debye theta and gamma are effective fit coefficients, strongly objective-sensitive and correlated; no independent calorimetric or vibrational validation.",
            "Full conditional covariance retained. Unweighted objectives do not use printed coordinate errors; no error-weighted, orthogonal-distance or pressure-scale systematic validation claimed.",
            "Phase/H-order and dehydration limitations of the Fukui lattice-average observations persist. MGD fitting does not resolve phase equilibrium or justify extrapolation.",
            "Diagnostic only: no catalog record, default or source coefficient changed.",
        ],
    }


def write_plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    dataset = data()
    old = list(report["current_continuous"]["parameters"].values())
    mgd = list(report["fits"]["mgd_volume"]["parameters"].values())
    empirical = list(report["fits"]["empirical_volume"]["parameters"].values())
    colors = {
        "published": "#c64c40",
        "current": "#7639a5",
        "mgd": "#187ca6",
        "empirical": "#26875d",
    }
    fig = plt.figure(figsize=(13.2, 11.5), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.16, top=0.80, hspace=0.31, wspace=0.18
    )
    residual_axes = []
    for i, source_moduli in enumerate(
        source()["sources"]["fukui_2003"]["parameters"]["isotherms"]
    ):
        t = source_moduli["temperature_k"]
        mask = (np.arange(48) < 37) & (dataset["temperature"] == t)
        p, v, sp, sv = (
            dataset[key][mask]
            for key in ["pressure", "volume", "pressure_error", "volume_error"]
        )
        order = np.argsort(p)
        grid = np.linspace(p.min(), p.max(), 250)
        seed = np.interp(grid, p[order], v[order])
        curves = {
            "published": comparison_volume(
                grid,
                reference_volume(t),
                source_moduli["K0_gpa"],
                source_moduli["K0_prime"],
            ),
            "current": volumes("empirical", grid, t, old, seed),
            "mgd": volumes("mgd", grid, t, mgd, seed),
            "empirical": volumes("empirical", grid, t, empirical, seed),
        }
        inner = outer[i // 2, i % 2].subgridspec(
            2, 1, height_ratios=[2.4, 1], hspace=0.06
        )
        ax = fig.add_subplot(inner[0])
        res = fig.add_subplot(inner[1], sharex=ax)
        residual_axes.append(res)
        ax.errorbar(
            p,
            v,
            xerr=sp,
            yerr=sv,
            fmt="o",
            color="#242424",
            ms=4,
            capsize=2,
            elinewidth=0.9,
            zorder=8,
        )
        res.errorbar(
            p,
            np.zeros_like(p),
            xerr=sp,
            yerr=sv,
            fmt="none",
            ecolor="#bbbbbb",
            capsize=2,
            elinewidth=0.8,
            zorder=0,
        )
        for key, style in [
            ("published", "-"),
            ("current", ":"),
            ("empirical", "--"),
            ("mgd", "-"),
        ]:
            ax.plot(grid, curves[key], color=colors[key], ls=style, lw=1.8)
            if key == "published":
                residual = (
                    comparison_volume(
                        p,
                        reference_volume(t),
                        source_moduli["K0_gpa"],
                        source_moduli["K0_prime"],
                    )
                    - v
                )
            elif key == "current":
                residual = volumes("empirical", p, t, old, v) - v
            else:
                residual = np.array(
                    report["fits"][
                        "mgd_volume" if key == "mgd" else "empirical_volume"
                    ]["row_volume_residuals_a3"]
                )[mask]
            res.plot(
                p,
                residual,
                "o" if key != "empirical" else "s",
                color=colors[key],
                ms=3.5,
                markerfacecolor="white" if key == "empirical" else colors[key],
                alpha=0.9,
            )
        ax.set_title(f"{t} K · {len(p)} observations", loc="left", fontsize=12)
        ax.set_ylabel("Unit-cell volume (Å³)")
        ax.tick_params(labelbottom=False)
        res.axhline(0, color="#555555", lw=0.7)
        res.set_xlabel("Pressure (GPa)")
        res.set_ylabel("ΔV (Å³)")
        res.set_ylim(-0.56, 0.56)
        res.set_yticks([-0.4, 0, 0.4])
        for axis in [ax, res]:
            axis.grid(alpha=0.16)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Brucite: does a joint BM3–MGD fit improve the agreement?", y=0.97, fontsize=17
    )
    fig.legend(
        handles=[
            Line2D([], [], color=colors[k], ls=s, lw=2, label=label)
            for k, s, label in [
                ("mgd", "-", "Joint BM3–MGD refit · 6 free parameters"),
                (
                    "empirical",
                    "--",
                    "Joint empirical thermal refit · 7 free parameters",
                ),
                (
                    "current",
                    ":",
                    "Current continuous model · ambient expansion refit + published compression parameters",
                ),
                (
                    "published",
                    "-",
                    "Original published isotherms · printed expansion and tabulated moduli",
                ),
            ]
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.935),
        frameon=False,
        fontsize=10,
    )
    fig.text(
        0.075,
        0.065,
        "Both joint fits minimize unweighted volume residuals to the same 37 P–V–T + 11 ambient T–V rows.\nΔV = calculated − observed volume at observed pressure. Grey residual bars show reported V and P errors about zero.\nMGD uses Tr = 300 K, n = 5 atoms per Mg(OH)₂ formula unit, integrated γ–Θ law and reference-subtracted Debye energy.\nIn-sample diagnostic fits: no phase-equilibrium, calorimetric or pressure-scale validation; catalog records unchanged.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-mgd-comparison"
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
        if REPORT.read_text(encoding="utf-8") != payload:
            raise SystemExit("Brucite MGD comparison is stale")
    else:
        REPORT.write_text(payload, encoding="utf-8")
    if args.plot:
        write_plot(report)
    print("Brucite MGD comparison is current.")
    for key, fit in report["fits"].items():
        print(
            key,
            json.dumps(
                {
                    "parameters": fit["parameters"],
                    "table4": fit["assessment"]["table4_all"],
                    "table1": fit["assessment"]["table1_ambient"],
                }
            ),
        )


if __name__ == "__main__":
    main()
