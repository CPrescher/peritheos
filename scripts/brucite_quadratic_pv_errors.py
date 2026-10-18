"""Compare quadratic and linear brucite moduli with exact P/V error fitting.

All 37 PVT rows and the fixed staged ambient law are unchanged. Quadratics
are fitted in a positive Bernstein basis, then reported as ordinary powers
of T-300. The independent-isotherm fits are a diagnostic, not catalog changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json

# Retain the archived ODRPACK cross-check without changing its numerical method.
# SciPy 1.17+ warns on import; migration is required before SciPy 1.19 removes it.
import warnings

import numpy as np

from scripts.brucite_replay_checks import check_report

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore", message=r"`scipy\.odr` is deprecated.*", category=DeprecationWarning
    )
    from scipy import odr
from scipy.optimize import least_squares

from scripts.brucite_joint_pv_errors import REPORT as LINEAR_REPORT
from scripts.brucite_joint_pv_errors import assessment, read_data
from scripts.brucite_reproduce import ROOT, bm3, bulk_modulus, comparison_volume

REPORT = ROOT / "docs/data/brucite-quadratic-pv-errors.json"
SPAN = 573.0
ORDER = [
    "K300_gpa",
    "K_linear_gpa_per_k",
    "K_quadratic_gpa_per_k2",
    "Kprime300",
    "Kprime_linear_per_k",
    "Kprime_quadratic_per_k2",
]


def basis(temperature, degree):
    u = (np.asarray(temperature) - 300) / SPAN
    if degree == 2:
        return np.column_stack([(1 - u) ** 2, 2 * u * (1 - u), u**2])
    if degree == 0:
        return np.ones((u.size, 1))
    raise ValueError(degree)


def moduli(coefficients, temperature, degree=2):
    design = basis(temperature, degree)
    n = degree + 1
    return design @ coefficients[:n], design @ coefficients[n : 2 * n]


def polynomial_transform():
    block = np.array(
        [[1, 0, 0], [-2 / SPAN, 2 / SPAN, 0], [1 / SPAN**2, -2 / SPAN**2, 1 / SPAN**2]]
    )
    result = np.zeros((6, 6))
    result[:3, :3] = block
    result[3:, 3:] = block
    return result


def residual(x, data, degree, orthogonal):
    ng = 2 * (degree + 1)
    k, kp = moduli(x[:ng], data["t"], degree)
    volume = data["v"] + data["sv"] * x[ng:] if orthogonal else data["v"]
    delta = bm3(volume, data["v0"], k, kp) - data["p"]
    return np.r_[delta / data["sp"], x[ng:]] if orthogonal else delta


def jacobian(x, data, degree, orthogonal):
    ng = 2 * (degree + 1)
    n = len(data["v"])
    k, kp = moduli(x[:ng], data["t"], degree)
    volume = data["v"] + data["sv"] * x[ng:] if orthogonal else data["v"]
    eta = (data["v0"] / volume) ** (1 / 3)
    dk = bm3(volume, data["v0"], k, kp) / k
    dkp = 1.125 * k * (eta**7 - eta**5) * (eta**2 - 1)
    design = basis(data["t"], degree)
    global_jac = np.column_stack([dk[:, None] * design, dkp[:, None] * design])
    if not orthogonal:
        return global_jac
    result = np.zeros((2 * n, ng + n))
    result[:n, :ng] = global_jac / data["sp"][:, None]
    result[np.arange(n), ng + np.arange(n)] = (
        -bulk_modulus(volume, data["v0"], k, kp) / volume * data["sv"] / data["sp"]
    )
    result[n:, ng:] = np.eye(n)
    return result


def fit_data(data, degree, starts, orthogonal):
    ng = 2 * (degree + 1)
    n = len(data["v"])
    lower = np.r_[np.full(degree + 1, 5.0), np.full(degree + 1, 3.0)]
    upper = np.r_[np.full(degree + 1, 100.0), np.full(degree + 1, 20.0)]
    if orthogonal:
        lower = np.r_[lower, (1e-10 - data["v"]) / data["sv"]]
        upper = np.r_[upper, np.full(n, np.inf)]
    fits = [
        least_squares(
            residual,
            np.r_[start, np.zeros(n)] if orthogonal else start,
            jac=jacobian,
            args=(data, degree, orthogonal),
            bounds=(lower, upper),
            x_scale="jac",
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-13,
            max_nfev=1500,
        )
        for start in starts
    ]
    best = min(fits, key=lambda fit: fit.fun @ fit.fun)
    j = jacobian(best.x, data, degree, orthogonal)
    if orthogonal:
        a = j[np.arange(n), ng + np.arange(n)]
        info = np.einsum("ni,nj,n->ij", j[:n, :ng], j[:n, :ng], 1 / (1 + a**2))
        covariance = np.linalg.inv(info)
    else:
        covariance = (
            np.linalg.inv(np.einsum("ni,nj->ij", j, j))
            * (best.fun @ best.fun)
            / (n - ng)
        )
    k, kp = moduli(best.x[:ng], data["t"], degree)
    evaluated = assessment(data["v0"], k, kp, data)
    if orthogonal:
        np.testing.assert_allclose(
            evaluated["projection"]["chi_square"], best.fun @ best.fun, atol=1e-8
        )
    return {
        "bernstein_coefficients": best.x[:ng].tolist(),
        "degree": degree,
        "solver_success": bool(best.success),
        "solver_message": best.message,
        "global_parameters": ng,
        "latent_volumes": n if orthogonal else 0,
        "residual_components": 2 * n if orthogonal else n,
        "degrees_of_freedom": n - ng,
        "objective_value": float(best.fun @ best.fun),
        "objective_type": "Exact P/V errors-in-variables chi-square"
        if orthogonal
        else "Unweighted pressure-residual SSE [GPa^2]",
        "chi_square": float(best.fun @ best.fun)
        if orthogonal
        else evaluated["projection"]["chi_square"],
        "reduced_chi_square": float(best.fun @ best.fun) / (n - ng)
        if orthogonal
        else None,
        "active_global_bounds": best.active_mask[:ng].tolist(),
        "bernstein_covariance": covariance.tolist(),
        "covariance_convention": "Absolute supplied P/V errors; inverse profiled information, latent-volume block eliminated by Schur complement; not scaled by reduced chi-square."
        if orthogonal
        else "SSE/(N-ng)*inverse(J.T J); exact measured volumes and equal pressure weights.",
        "multistart": [
            {
                "success": bool(f.success),
                "objective": float(f.fun @ f.fun),
                "nfev": f.nfev,
                "optimality": float(f.optimality),
            }
            for f in fits
        ],
        **evaluated,
    }


def verify_odr(data, coefficients):
    """Independent ODRPACK solver; temperature arrays are fixed in the closure."""
    real = odr.RealData(data["v"], data["p"], sx=data["sv"], sy=data["sp"])
    model = odr.Model(
        lambda beta, volume: bm3(volume, data["v0"], *moduli(beta, data["t"]))
    )
    solver = odr.ODR(
        real, model, beta0=coefficients, maxit=500, sstol=1e-12, partol=1e-12
    )
    result = solver.run()
    return {
        "status": int(result.info),
        "stop_reasons": result.stopreason,
        "chi_square": float(result.sum_square),
        "degrees_of_freedom": len(data["v"]) - 6,
        "bernstein_coefficients": result.beta.tolist(),
        "scaled_covariance_not_used": True,
    }


def native_verification(ambient, coefficients):
    from peritheos.eos.rt import BM3

    max_p, max_v, min_bulk = 0.0, 0.0, np.inf
    for t in [300, 400, 473, 573, 673, 773, 873]:
        k, kp = (float(value[0]) for value in moduli(coefficients, np.array([t])))
        v0 = ambient["V300_a3"] * np.exp(ambient["alpha_per_k"] * (t - 300))
        native = BM3(v0, k, kp)
        pressures = np.array([0.0, 6.0, 10.0, 15.0])
        v = comparison_volume(pressures, v0, k, kp)
        max_p = max(max_p, float(np.max(abs(native.pressure(v) - bm3(v, v0, k, kp)))))
        max_v = max(max_v, float(np.max(abs(native.volume(pressures) - v))))
        min_bulk = min(min_bulk, float(np.min(native.bulk_modulus(v))))
    if max_p > 1e-10 or max_v > 1e-8 or min_bulk <= 0:
        raise ValueError("Quadratic model failed native evaluation/inversion check")
    return {
        "temperatures_k": [300, 400, 473, 573, 673, 773, 873],
        "max_pressure_difference_gpa": max_p,
        "max_volume_difference_a3": max_v,
        "minimum_sampled_bulk_modulus_gpa": min_bulk,
    }


def reproduce():
    staged, data = read_data()
    linear = json.loads(LINEAR_REPORT.read_text(encoding="utf-8"))
    if (
        linear["fixed_ambient_law"] != staged["stage1_ambient"]
        or linear["source_sha256"] != staged["source_sha256"]
    ):
        raise ValueError("Linear comparator does not use the same fixed expansion/data")
    pars = linear["fits"]["joint_pv_errors"]["parameters"]
    k_end = [pars["K300_gpa"], pars["K300_gpa"] + SPAN * pars["dK_dT_gpa_per_k"]]
    kp_end = [pars["Kprime300"], pars["Kprime300"] + SPAN * pars["dKprime_dT_per_k"]]
    embedded = np.array(
        [k_end[0], np.mean(k_end), k_end[1], kp_end[0], np.mean(kp_end), kp_end[1]]
    )
    starts = [
        embedded,
        np.array([42, 28, 25, 6.6, 8.3, 9.7]),
        np.array([41.8, 32, 23.9, 6.6, 8, 10]),
    ]
    baseline = fit_data(data, 2, starts, False)
    quadratic = fit_data(
        data, 2, [np.array(baseline["bernstein_coefficients"]), *starts], True
    )
    transform = polynomial_transform()
    for fit in [baseline, quadratic]:
        physical = transform @ fit["bernstein_coefficients"]
        cov = np.einsum(
            "ij,jk,lk->il", transform, np.array(fit["bernstein_covariance"]), transform
        )
        fit.update(
            parameters=dict(zip(ORDER, physical.tolist())),
            parameter_order=ORDER,
            covariance=cov.tolist(),
            standard_errors=dict(zip(ORDER, np.sqrt(np.diag(cov)).tolist())),
        )
    quadratic["covariance_scaled_by_observed_scatter"] = (
        np.array(quadratic["covariance"]) * quadratic["reduced_chi_square"]
    ).tolist()
    checked = verify_odr(data, quadratic["bernstein_coefficients"])
    difference = checked["chi_square"] - quadratic["chi_square"]
    if checked["status"] not in [1, 2, 3] or abs(difference) > 1e-5:
        raise ValueError("Independent ODRPACK solver disagrees with quadratic fit")
    checked["objective_difference"] = difference
    independent = []
    for t in [300, 473, 673, 873]:
        mask = data["t"] == t
        subset = {key: value[mask] for key, value in data.items()}
        fit = fit_data(
            subset,
            0,
            [np.array([40, 7]) if t == 300 else np.array([28, 9]), np.array([35, 8])],
            True,
        )
        independent.append(
            {
                "temperature_k": t,
                "source_rows": (np.flatnonzero(mask) + 1).tolist(),
                **fit,
            }
        )
    chi_independent = sum(fit["chi_square"] for fit in independent)
    comparisons = {
        "quadratic_pv_errors": quadratic,
        "linear_pv_errors": linear["fits"]["joint_pv_errors"],
        "published_isotherms": linear["fits"]["published_isotherms"],
    }
    per_temperature = []
    for t in [300, 473, 673, 873]:
        mask = data["t"] == t
        per_temperature.append(
            {
                "temperature_k": t,
                "observations": int(mask.sum()),
                "source_rows": (np.flatnonzero(mask) + 1).tolist(),
                "comparisons": {
                    key: {
                        "volume_rmse_a3": float(
                            np.sqrt(
                                np.mean(
                                    np.array(fit["row_volume_residuals_a3"])[mask] ** 2
                                )
                            )
                        ),
                        "pressure_rmse_gpa": float(
                            np.sqrt(
                                np.mean(
                                    np.array(fit["row_pressure_residuals_gpa"])[mask]
                                    ** 2
                                )
                            )
                        ),
                        "chi_square_contribution": float(
                            np.sum(
                                np.array(
                                    fit["projection"][
                                        "signed_orthogonal_distance_sigma"
                                    ]
                                )[mask]
                                ** 2
                            )
                        ),
                        "signed_orthogonal_distance_sigma": np.array(
                            fit["projection"]["signed_orthogonal_distance_sigma"]
                        )[mask].tolist(),
                    }
                    for key, fit in comparisons.items()
                },
            }
        )
    return {
        "format": "peritheos.brucite-quadratic-pv-errors",
        "format_version": 1,
        "source_sha256": staged["source_sha256"],
        "linear_report_sha256": hashlib.sha256(LINEAR_REPORT.read_bytes()).hexdigest(),
        "fixed_ambient_law": staged["stage1_ambient"],
        "source_rows": list(range(1, 38)),
        "objective": linear["objective"].replace(
            "Four global linear-modulus coefficients",
            "Six global quadratic-modulus coefficients",
        ),
        "temperature_law": "K(T)=K300+a*(T-300)+b*(T-300)^2 and Kprime(T)=Kprime300+c*(T-300)+d*(T-300)^2. b,d multiply the full squared term (no one-half); second derivatives are 2*b and 2*d. Six coefficients jointly fitted to all 37 raw PVT rows with both coordinate errors.",
        "numerical_parameterization": "For u=(T-300)/573, each modulus is c0*(1-u)^2+2*c1*u*(1-u)+c2*u^2. K coefficients bounded 5-100 GPa and Kprime coefficients 3-20, guaranteeing positive reference moduli over 300-873 K. All solutions checked for active bounds. Ordinary polynomial coefficients/covariance are transformed exactly.",
        "observational_assumptions": linear["observational_assumptions"],
        "ambient_data_treatment": linear["ambient_data_treatment"],
        "fits": {**comparisons, "quadratic_unweighted_pressure_baseline": baseline},
        "per_temperature": per_temperature,
        "independent_isotherm_pv_error_control": {
            "isotherms": independent,
            "chi_square": chi_independent,
            "global_parameters": 8,
            "degrees_of_freedom": 29,
            "reduced_chi_square": chi_independent / 29,
            "qualification": "Four independent BM3 fits with the same fixed V0(T), same P/V errors and all rows. Diagnostic flexibility control, not a continuous EOS and not the earlier unweighted staged moduli.",
        },
        "verification": {
            "independent_odrpack": checked,
            "native_bm3": native_verification(
                staged["stage1_ambient"], quadratic["bernstein_coefficients"]
            ),
        },
        "qualifications": [
            "The quadratic model adds two global parameters: N-p=31 degrees of freedom versus 33 for the linear model. Compare both raw chi-square and reduced chi-square; this is an in-sample comparison, not external prediction validation.",
            "Exact same fixed ambient expansion, source observations, P/V errors and independence assumptions as the linear fit. No outlier removal, pressure recalibration, robust clipping or ambient-law adjustment.",
            "Primary covariance is local and absolute-sigma, after profiling all 37 latent volumes. Separate dispersion-scaled covariance is diagnostic only. Both omit fixed ambient-law uncertainty, unreported error correlations, temperature errors and model/phase discrepancy.",
            "Quadratics are empirical interpolating forms over 300-873 K; only four distinct temperatures constrain their shape. Do not extrapolate or interpret coefficients as independently established thermodynamics.",
            "The independent-isotherm control can retain scatter even without a temperature law; that diagnostic does not identify a physical cause or prove the source observations wrong.",
            "No source, catalog, scientific validation status, original record or default changed.",
        ],
    }


def write_plots(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    staged, data = read_data()
    colors = {
        "quadratic_pv_errors": "#187ca6",
        "linear_pv_errors": "#8a8f93",
        "published_isotherms": "#c64c40",
    }
    fig = plt.figure(figsize=(13.2, 11.2), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.16, top=0.80, hspace=0.32, wspace=0.18
    )
    limit = max(
        6,
        np.ceil(
            max(
                np.max(
                    abs(
                        np.array(
                            report["fits"][key]["projection"][
                                "signed_orthogonal_distance_sigma"
                            ]
                        )
                    )
                )
                for key in colors
            )
        ),
    )
    for i, group in enumerate(report["per_temperature"]):
        t = group["temperature_k"]
        mask = data["t"] == t
        p, v, sp, sv = (data[key][mask] for key in ["p", "v", "sp", "sv"])
        grid = np.linspace(p.min(), p.max(), 250)
        inner = outer[i // 2, i % 2].subgridspec(
            2, 1, height_ratios=[2.5, 1], hspace=0.06
        )
        ax = fig.add_subplot(inner[0])
        res = fig.add_subplot(inner[1], sharex=ax)
        res.axhspan(-2, 2, color="#f3f3f3", zorder=-10)
        res.axhspan(-1, 1, color="#e5e5e5", zorder=-9)
        for key in ["published_isotherms", "linear_pv_errors", "quadratic_pv_errors"]:
            if key == "published_isotherms":
                curve = staged["isotherms"][i]["published"]
                params = [curve[name] for name in ["V0_a3", "K0_gpa", "K0_prime"]]
            else:
                if key == "quadratic_pv_errors":
                    k, kp = (
                        float(x[0])
                        for x in moduli(
                            report["fits"][key]["bernstein_coefficients"], np.array([t])
                        )
                    )
                else:
                    pars = report["fits"][key]["parameters"]
                    k = pars["K300_gpa"] + pars["dK_dT_gpa_per_k"] * (t - 300)
                    kp = pars["Kprime300"] + pars["dKprime_dT_per_k"] * (t - 300)
                params = [float(data["v0"][mask][0]), k, kp]
            ax.plot(
                grid,
                comparison_volume(grid, *params),
                color=colors[key],
                ls="--" if key == "linear_pv_errors" else "-",
                lw=1.8,
            )
            res.plot(
                p,
                group["comparisons"][key]["signed_orthogonal_distance_sigma"],
                {
                    "quadratic_pv_errors": "o",
                    "linear_pv_errors": "^",
                    "published_isotherms": "s",
                }[key],
                color=colors[key],
                ms=4,
                markerfacecolor=colors[key]
                if key == "quadratic_pv_errors"
                else "white",
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
        ax.set_title(f"{t} K · {len(p)} observations", loc="left", fontsize=12)
        ax.set_ylabel("Unit-cell volume (Å³)")
        ax.tick_params(labelbottom=False)
        ax.text(
            0.97,
            0.94,
            f"Volume RMS (Å³)\nPublished: {group['comparisons']['published_isotherms']['volume_rmse_a3']:.3f}\nLinear P/V fit: {group['comparisons']['linear_pv_errors']['volume_rmse_a3']:.3f}\nQuadratic P/V fit: {group['comparisons']['quadratic_pv_errors']['volume_rmse_a3']:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        res.axhline(0, color="#555555", lw=0.8)
        res.set_ylim(-limit - 0.25, limit + 0.25)
        res.set_yticks([-6, -3, 0, 3, 6])
        res.set_xlabel("Pressure (GPa)")
        res.set_ylabel("Signed distance\n(error units)")
        for axis in [ax, res]:
            axis.grid(alpha=0.15)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Brucite: quadratic versus linear temperature dependence of both moduli",
        y=0.975,
        fontsize=16,
    )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="o",
                color="#242424",
                ls="none",
                label="Published Table 4 observations + reported P/V errors",
            ),
            *[
                Line2D(
                    [],
                    [],
                    color=colors[key],
                    ls="--" if key == "linear_pv_errors" else "-",
                    lw=2,
                    label=label,
                )
                for key, label in [
                    (
                        "quadratic_pv_errors",
                        "Joint quadratic K₀(T), K′(T): both P/V errors, fixed ambient expansion",
                    ),
                    (
                        "linear_pv_errors",
                        "Joint linear K₀(T), K′(T): same P/V errors and fixed expansion",
                    ),
                    (
                        "published_isotherms",
                        "Published individual isotherms: original expansion + tabulated moduli",
                    ),
                ]
            ],
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.945),
        frameon=False,
        fontsize=10,
    )
    quad = report["fits"]["quadratic_pv_errors"]
    linear = report["fits"]["linear_pv_errors"]
    fig.text(
        0.075,
        0.069,
        f"Same 37 P–V–T rows, exact P/V errors-in-variables objective, and fixed ambient V₀(T); only modulus temperature dependence changes.\nLinear: χ² = {linear['chi_square']:.2f}, 4 global parameters, dof = 33, reduced χ² = {linear['reduced_chi_square']:.2f}. Quadratic: χ² = {quad['chi_square']:.2f}, 6 parameters, dof = 31, reduced χ² = {quad['reduced_chi_square']:.2f}.\nLower panels: signed orthogonal distance in reported P/V error units; shaded bands are ±1 and ±2. Model and ambient-law uncertainty not shown.\nOnly four temperatures constrain the quadratic form. In-sample diagnostic comparison; source phase/calibration limits persist, catalog unchanged.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-quadratic-pv-errors"
    fig.savefig(stem.with_suffix(".png"), dpi=200)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), layout="constrained")
    grid = np.linspace(300, 873, 250)
    qk, qkp = moduli(quad["bernstein_coefficients"], grid)
    lp = linear["parameters"]
    for i, ax in enumerate(axes):
        y = qk if i == 0 else qkp
        ly = (
            lp["K300_gpa"] + lp["dK_dT_gpa_per_k"] * (grid - 300)
            if i == 0
            else lp["Kprime300"] + lp["dKprime_dT_per_k"] * (grid - 300)
        )
        ax.plot(
            grid,
            y,
            color=colors["quadratic_pv_errors"],
            lw=2,
            label="Joint quadratic P/V-error fit",
        )
        ax.plot(
            grid,
            ly,
            color=colors["linear_pv_errors"],
            ls="--",
            lw=1.8,
            label="Joint linear P/V-error fit",
        )
        control = report["independent_isotherm_pv_error_control"]["isotherms"]
        ax.errorbar(
            [f["temperature_k"] for f in control],
            [f["bernstein_coefficients"][i] for f in control],
            yerr=[np.sqrt(f["bernstein_covariance"][i][i]) for f in control],
            fmt="o",
            color="#242424",
            capsize=3,
            label="Independent isotherms, same P/V errors",
        )
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("K₀ (GPa)" if i == 0 else "K′ (dimensionless)")
        ax.grid(alpha=0.16)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8)
    fig.suptitle(
        "Modulus temperature dependence: joint fits and independent-isotherm control"
    )
    stem = ROOT / "docs/data/brucite-quadratic-moduli"
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
        check_report(report, json.loads(REPORT.read_text(encoding="utf-8")))

    else:
        REPORT.write_text(payload, encoding="utf-8")
    if args.plot:
        write_plots(report)
    print("Brucite quadratic P/V-errors report is current.")
    for key in ["linear_pv_errors", "quadratic_pv_errors"]:
        fit = report["fits"][key]
        print(
            key,
            json.dumps(
                {
                    name: fit[name]
                    for name in [
                        "chi_square",
                        "degrees_of_freedom",
                        "reduced_chi_square",
                        "volume_rmse_a3",
                        "pressure_rmse_gpa",
                    ]
                }
            ),
        )
    print(
        "independent-isotherm control",
        report["independent_isotherm_pv_error_control"]["chi_square"],
        report["independent_isotherm_pv_error_control"]["reduced_chi_square"],
    )


if __name__ == "__main__":
    main()
