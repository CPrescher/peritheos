"""Joint fixed-expansion brucite fit with exact pressure/volume error residuals.

All 37 Table 4 points are fitted together. The 11 ambient Table 1 points
already set the staged constant-alpha law, which is copied without refitting.
This script produces diagnostic reports and plots, never catalog changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import scipy
from scipy.optimize import least_squares

from scripts.brucite_replay_checks import check_report
from scripts.brucite_reproduce import (
    ROOT,
    bm3,
    bulk_modulus,
    comparison_volume,
    metrics,
    observations,
)
from scripts.brucite_staged_comparison import REPORT as STAGED_REPORT
from scripts.brucite_staged_continuous import REPORT as LINE_REPORT
from scripts.brucite_staged_continuous import runtime_model

REPORT = ROOT / "docs/data/brucite-joint-pv-errors.json"
ORDER = ["K300_gpa", "Kprime300", "dK_dT_gpa_per_k", "dKprime_dT_per_k"]
SCALE = np.array([1, 1, 0.001, 0.001])
LOWER = np.array([25, 3, -40, -2])
UPPER = np.array([80, 15, 20, 20])


def read_data():
    staged = json.loads(STAGED_REPORT.read_text(encoding="utf-8"))
    for name, checksum in staged["source_sha256"].items():
        path = ROOT / "peritheos/data/datasets" / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
            raise ValueError("Staged source checksum is stale")
    rows = observations("fukui_2003", "table4")
    values = {
        key: np.array([float(row[column]) for row in rows])
        for key, column in [
            ("v", "cell_volume_a3"),
            ("p", "pressure_gpa"),
            ("t", "temperature_k"),
            ("sv", "cell_volume_error_a3"),
            ("sp", "pressure_error_gpa"),
        ]
    }
    if len(rows) != 37 or np.any(values["sv"] <= 0) or np.any(values["sp"] <= 0):
        raise ValueError(
            "Expected all 37 observations with positive reported P/V errors"
        )
    ambient = staged["stage1_ambient"]
    values["v0"] = ambient["V300_a3"] * np.exp(
        ambient["alpha_per_k"] * (values["t"] - 300)
    )
    values["u"] = (values["t"] - 300) / 1000
    return staged, values


def pressure_scaled(x, volume, data):
    return bm3(volume, data["v0"], x[0] + x[2] * data["u"], x[1] + x[3] * data["u"])


def parameter_derivatives(x, volume, data):
    """Independent analytic derivatives in the scaled four-parameter basis."""
    k = x[0] + x[2] * data["u"]
    eta = (data["v0"] / volume) ** (1 / 3)
    dk = pressure_scaled(x, volume, data) / k
    dkp = 1.125 * k * (eta**7 - eta**5) * (eta**2 - 1)
    return np.column_stack([dk, dkp, dk * data["u"], dkp * data["u"]])


def residual(x, data, orthogonal):
    if not orthogonal:
        return pressure_scaled(x, data["v"], data) - data["p"]
    volume = data["v"] + data["sv"] * x[4:]
    return np.r_[(pressure_scaled(x[:4], volume, data) - data["p"]) / data["sp"], x[4:]]


def jacobian(x, data, orthogonal):
    if not orthogonal:
        return parameter_derivatives(x, data["v"], data)
    n = len(data["v"])
    volume = data["v"] + data["sv"] * x[4:]
    k, kp = x[0] + x[2] * data["u"], x[1] + x[3] * data["u"]
    result = np.zeros((2 * n, 4 + n))
    result[:n, :4] = parameter_derivatives(x[:4], volume, data) / data["sp"][:, None]
    result[np.arange(n), 4 + np.arange(n)] = (
        -bulk_modulus(volume, data["v0"], k, kp) / volume * data["sv"] / data["sp"]
    )
    result[n:, 4:] = np.eye(n)
    return result


def fit_data(data, starts, orthogonal):
    n = len(data["v"])
    lower, upper = LOWER.copy(), UPPER.copy()
    if orthogonal:
        lower = np.r_[lower, (1e-10 - data["v"]) / data["sv"]]
        upper = np.r_[upper, np.full(n, np.inf)]
    fits = []
    for start in starts:
        initial = np.r_[start, np.zeros(n)] if orthogonal else np.asarray(start)
        fits.append(
            least_squares(
                residual,
                initial,
                jac=jacobian,
                args=(data, orthogonal),
                bounds=(lower, upper),
                x_scale="jac",
                ftol=1e-13,
                xtol=1e-13,
                gtol=1e-13,
                max_nfev=1200,
            )
        )
    best = min(fits, key=lambda fit: fit.fun @ fit.fun)
    j = jacobian(best.x, data, orthogonal)
    if orthogonal:
        # Schur complement profiles out the 37 independently fitted volumes.
        jtheta = j[:n, :4]
        latent_derivative = j[np.arange(n), 4 + np.arange(n)]
        # The diagonal latent block simplifies the Schur complement to this
        # sum, avoiding subtraction of nearly equal information matrices.
        information = np.einsum(
            "ni,nj,n->ij", jtheta, jtheta, 1 / (1 + latent_derivative**2)
        )
        covariance = np.linalg.inv(information)
    else:
        covariance = (
            np.linalg.inv(np.einsum("ni,nj->ij", j, j))
            * (best.fun @ best.fun)
            / (n - 4)
        )
    covariance *= SCALE[:, None] * SCALE[None, :]
    return (
        best,
        covariance,
        [
            {
                "success": bool(f.success),
                "objective": float(f.fun @ f.fun),
                "parameters": dict(zip(ORDER, (f.x[:4] * SCALE).tolist())),
                "nfev": f.nfev,
                "optimality": float(f.optimality),
            }
            for f in fits
        ],
    )


def project_curve(v0, k, kp, data):
    """Profile the same exact two-coordinate objective for a fixed comparator."""
    n = len(data["v"])

    def objective(z):
        v = data["v"] + data["sv"] * z
        return np.r_[(bm3(v, v0, k, kp) - data["p"]) / data["sp"], z]

    def derivative(z):
        v = data["v"] + data["sv"] * z
        a = -bulk_modulus(v, v0, k, kp) / v * data["sv"] / data["sp"]
        return np.vstack([np.diag(a), np.eye(n)])

    projection = least_squares(
        objective,
        np.zeros(n),
        jac=derivative,
        bounds=((1e-10 - data["v"]) / data["sv"], np.full(n, np.inf)),
        x_scale="jac",
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
    )
    z = projection.x
    vp = data["v"] + data["sv"] * z
    dp = (bm3(vp, v0, k, kp) - data["p"]) / data["sp"]
    sign = np.sign(bm3(data["v"], v0, k, kp) - data["p"])
    return {
        "adjusted_volume_a3": vp.tolist(),
        "adjusted_pressure_gpa": bm3(vp, v0, k, kp).tolist(),
        "volume_correction_sigma": z.tolist(),
        "pressure_residual_sigma": dp.tolist(),
        "signed_orthogonal_distance_sigma": (sign * np.hypot(dp, z)).tolist(),
        "chi_square": float(dp @ dp + z @ z),
        "projection_success": bool(projection.success),
    }


def assessment(v0, k, kp, data):
    v0, k, kp = np.broadcast_arrays(v0, k, kp, data["v"])[0:3]
    dp = bm3(data["v"], v0, k, kp) - data["p"]
    predicted_volume = np.array(
        [
            comparison_volume(p, vv, kk, kkprime)[0]
            for p, vv, kk, kkprime in zip(data["p"], v0, k, kp)
        ]
    )
    dv = predicted_volume - data["v"]
    return {
        "row_pressure_residuals_gpa": dp.tolist(),
        "row_volume_residuals_a3": dv.tolist(),
        **metrics(dp),
        "volume_rmse_a3": float(np.sqrt(np.mean(dv**2))),
        "projection": project_curve(v0, k, kp, data),
    }


def native_verification(data, ambient, best):
    """Use the public native fitting API as an independent solver check."""
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import ThermalReferenceStateEOS
    from peritheos.fitting import fit_joint_eos

    names = ["rt_eos.K0", "rt_eos.K0_prime", "dK_dT", "dKprime_dT"]
    fit = fit_joint_eos(
        ThermalReferenceStateEOS,
        BM3,
        data["v"],
        data["t"],
        data["p"],
        initial=dict(zip(names, (best.x[:4] * SCALE).tolist())),
        fixed={
            "rt_eos.V0": ambient["V300_a3"],
            "Tr": 300,
            "alpha0": ambient["alpha_per_k"],
        },
        bounds=dict(zip(names, zip(LOWER * SCALE, UPPER * SCALE))),
        pressure_sigma=data["sp"],
        volume_sigma=data["sv"],
        absolute_sigma=True,
        max_nfev=1500,
    )
    difference = fit.chi_square - best.fun @ best.fun
    if not fit.success or abs(difference) > 2e-5:
        raise ValueError(
            "Public native errors-in-variables fit disagrees with analytic replay"
        )
    return {
        "success": fit.success,
        "objective_difference": float(difference),
        "chi_square": fit.chi_square,
        "degrees_of_freedom": fit.degrees_of_freedom,
        "physical_parameters": {key: fit.parameters[key] for key in names},
        "maximum_parameter_difference_in_native_standard_errors": float(
            max(
                abs(fit.parameters[key] - value) / fit.standard_errors[key]
                for key, value in zip(names, best.x[:4] * SCALE)
            )
        ),
    }


def reproduce():
    staged, data = read_data()
    line = json.loads(LINE_REPORT.read_text(encoding="utf-8"))["linear_moduli"][
        "parameters"
    ]
    starts = [[41.8, 6.6, -31.9, 5.71], np.array([line[key] for key in ORDER]) / SCALE]
    baseline, baseline_covariance, baseline_starts = fit_data(data, starts, False)
    best, covariance, multistart = fit_data(data, [baseline.x, *starts], True)
    parameters = dict(zip(ORDER, (best.x[:4] * SCALE).tolist()))
    ambient = staged["stage1_ambient"]
    definitions = {}
    for key, fit, cov, runs in [
        ("unweighted_joint", baseline, baseline_covariance, baseline_starts),
        ("joint_pv_errors", best, covariance, multistart),
    ]:
        x = fit.x[:4]
        k, kp = x[0] + x[2] * data["u"], x[1] + x[3] * data["u"]
        definitions[key] = {
            "parameters": dict(zip(ORDER, (x * SCALE).tolist())),
            "solver_success": bool(fit.success),
            "solver_message": fit.message,
            "parameter_order": ORDER,
            "covariance": cov.tolist(),
            "standard_errors": dict(zip(ORDER, np.sqrt(np.diag(cov)).tolist())),
            "covariance_convention": "Absolute reported P/V errors; inverse profiled J.T J, 37 latent volumes eliminated by Schur complement. No chi-square scaling."
            if key == "joint_pv_errors"
            else "SSE/(37-4) times inverse(J.T J); exact measured volumes, equal pressure residual weights.",
            "active_global_bounds": fit.active_mask[:4].tolist(),
            "multistart": runs,
            **assessment(data["v0"], k, kp, data),
        }
    weighted = definitions["joint_pv_errors"]
    chi_square = float(best.fun @ best.fun)
    reduced = chi_square / (37 - 4)
    weighted.update(
        chi_square=chi_square,
        reduced_chi_square=reduced,
        degrees_of_freedom=33,
        global_parameters=4,
        latent_volumes=37,
        residual_components=74,
        covariance_scaled_by_observed_scatter=(covariance * reduced).tolist(),
        scaled_covariance_qualification="Alternative dispersion-scaled diagnostic only; inflation does not repair model inadequacy or unmodelled systematics.",
    )
    np.testing.assert_allclose(
        weighted["projection"]["chi_square"], chi_square, atol=1e-8
    )
    published_v0, published_k, published_kp = (
        np.concatenate(
            [
                np.full(s["observations"], s["published"][field])
                for s in staged["isotherms"]
            ]
        )
        for field in ["V0_a3", "K0_gpa", "K0_prime"]
    )
    # Table 4 is grouped by temperature; preserve its original row order.
    definitions["published_isotherms"] = assessment(
        published_v0, published_k, published_kp, data
    )
    per_temperature = []
    for t in [300, 473, 673, 873]:
        mask = data["t"] == t
        per_temperature.append(
            {
                "temperature_k": t,
                "observations": int(mask.sum()),
                "source_rows": (np.flatnonzero(mask) + 1).tolist(),
                "pressure_range_gpa": [
                    float(data["p"][mask].min()),
                    float(data["p"][mask].max()),
                ],
                "comparisons": {
                    key: {
                        "pressure_rmse_gpa": float(
                            np.sqrt(
                                np.mean(
                                    np.array(value["row_pressure_residuals_gpa"])[mask]
                                    ** 2
                                )
                            )
                        ),
                        "volume_rmse_a3": float(
                            np.sqrt(
                                np.mean(
                                    np.array(value["row_volume_residuals_a3"])[mask]
                                    ** 2
                                )
                            )
                        ),
                        "chi_square_contribution": float(
                            np.sum(
                                np.array(
                                    value["projection"][
                                        "signed_orthogonal_distance_sigma"
                                    ]
                                )[mask]
                                ** 2
                            )
                        ),
                        "signed_orthogonal_distance_sigma": np.array(
                            value["projection"]["signed_orthogonal_distance_sigma"]
                        )[mask].tolist(),
                    }
                    for key, value in definitions.items()
                },
            }
        )
    expected = pressure_scaled(best.x[:4], data["v"], data)
    runtime = runtime_model(
        {"fixed_ambient_law": ambient, "linear_moduli": {"parameters": parameters}}
    )
    difference = float(np.max(abs(runtime.pressure(data["v"], data["t"]) - expected)))
    if difference > 1e-10:
        raise ValueError("Joint fit does not match continuous runtime model")
    return {
        "format": "peritheos.brucite-joint-pv-errors",
        "format_version": 1,
        "source_sha256": staged["source_sha256"],
        "staged_report_sha256": hashlib.sha256(STAGED_REPORT.read_bytes()).hexdigest(),
        "software": {"numpy": np.__version__, "scipy": scipy.__version__},
        "fixed_ambient_law": ambient,
        "observations": 37,
        "source_rows": list(range(1, 38)),
        "objective": "Minimize sum_i [(P_BM3(Vtrue_i,T_i;theta)-Pobs_i)/sigmaP_i]^2 + [(Vtrue_i-Vobs_i)/sigmaV_i]^2. Four global linear-modulus coefficients and 37 latent volumes free. Exact errors-in-variables objective, not a parameter-dependent effective-error denominator approximation.",
        "model": "V0(T)=fixed V300*exp(fixed alpha*(T-300)); K0(T)=K300+dK_dT*(T-300); Kprime(T)=Kprime300+dKprime_dT*(T-300). Both 300 K compression intercepts and temperature slopes fitted jointly to all 37 rows.",
        "global_parameter_bounds": dict(
            zip(ORDER, (np.column_stack([LOWER * SCALE, UPPER * SCALE])).tolist())
        ),
        "observational_assumptions": "Printed P and V ESDs treated as independent Gaussian standard deviations within/between rows; no P-V or inter-row covariance supplied. Temperatures fixed because row-wise T errors are unreported. Au pressure-scale systematics are not included.",
        "ambient_data_treatment": "The 11 Table 1 ambient observations already determined the fixed expansion law. With V0(T) fixed, their predicted P=0 volumes do not depend on the four compression parameters; they add no information to this fit. No synthetic ambient pressure errors are assigned and the ambient law is not reoptimized.",
        "fits": definitions,
        "per_temperature": per_temperature,
        "verification": {
            "runtime_max_pressure_difference_gpa": difference,
            "native_public_fit": native_verification(data, ambient, best),
        },
        "qualifications": [
            "Lower chi-square under a weighting scheme does not imply lower unweighted volume/pressure RMS. The error-aware and explicit unweighted joint baselines use the same fixed ambient law and linear temperature model.",
            "Reduced chi-square is conditional on reported errors, the fixed ambient law and the linear-modulus model. A value far above one signals unresolved scatter relative to those assumptions; it is not a source-error or phase diagnosis.",
            "Absolute-sigma covariance and separately dispersion-scaled covariance are stored. Both exclude fixed ambient uncertainty, unknown pressure-scale correlations, temperature errors and phase/model discrepancy.",
            "Adjusted latent coordinates are diagnostic projections, not replacements for source observations. Published isotherms are assessed with their original coefficients and printed expansion.",
            "Fit coverage is the four measured temperature slices (300/473/673/873 K); interpolation is an empirical model, not independent phase/stability validation. Source phase and calibration caveats remain.",
            "No source, catalog record, scientific status or default changed.",
        ],
    }


def write_plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    staged, data = read_data()
    colors = {
        "joint_pv_errors": "#187ca6",
        "unweighted_joint": "#8a8f93",
        "published_isotherms": "#c64c40",
    }
    styles = {
        "joint_pv_errors": "-",
        "unweighted_joint": "--",
        "published_isotherms": "-",
    }
    fig = plt.figure(figsize=(13.2, 11.2), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.075, right=0.98, bottom=0.16, top=0.80, hspace=0.32, wspace=0.18
    )
    distance_limit = max(
        6,
        np.ceil(
            max(
                np.max(
                    abs(np.array(fit["projection"]["signed_orthogonal_distance_sigma"]))
                )
                for fit in report["fits"].values()
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
        for key in ["published_isotherms", "unweighted_joint", "joint_pv_errors"]:
            if key == "published_isotherms":
                source = staged["isotherms"][i]["published"]
                parameters = [
                    source[field] for field in ["V0_a3", "K0_gpa", "K0_prime"]
                ]
            else:
                pars = report["fits"][key]["parameters"]
                parameters = [
                    float(data["v0"][mask][0]),
                    pars["K300_gpa"] + pars["dK_dT_gpa_per_k"] * (t - 300),
                    pars["Kprime300"] + pars["dKprime_dT_per_k"] * (t - 300),
                ]
            ax.plot(
                grid,
                comparison_volume(grid, *parameters),
                color=colors[key],
                ls=styles[key],
                lw=1.8,
            )
            res.plot(
                p,
                group["comparisons"][key]["signed_orthogonal_distance_sigma"],
                {
                    "joint_pv_errors": "o",
                    "unweighted_joint": "^",
                    "published_isotherms": "s",
                }[key],
                color=colors[key],
                ms=4,
                markerfacecolor=colors[key] if key == "joint_pv_errors" else "white",
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
            f"Volume RMS (Å³)\nPublished: {group['comparisons']['published_isotherms']['volume_rmse_a3']:.3f}\nUnweighted joint: {group['comparisons']['unweighted_joint']['volume_rmse_a3']:.3f}\nP/V-error joint: {group['comparisons']['joint_pv_errors']['volume_rmse_a3']:.3f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )
        res.axhline(0, color="#555555", lw=0.8)
        res.set_ylim(-distance_limit - 0.25, distance_limit + 0.25)
        res.set_yticks([-6, -3, 0, 3, 6])
        res.set_xlabel("Pressure (GPa)")
        res.set_ylabel("Signed distance\n(error units)")
        for axis in [ax, res]:
            axis.grid(alpha=0.15)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Brucite: joint P–V–T fit with pressure and volume uncertainties",
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
                label="Published Table 4 observations + reported P/V errors",
            ),
            *[
                Line2D([], [], color=colors[key], ls=styles[key], lw=2, label=label)
                for key, label in [
                    (
                        "joint_pv_errors",
                        "Joint fit with both P and V errors; ambient expansion fixed",
                    ),
                    (
                        "unweighted_joint",
                        "Joint unweighted pressure-residual baseline; same fixed expansion",
                    ),
                    (
                        "published_isotherms",
                        "Published separate isotherms: original expansion + tabulated moduli",
                    ),
                ]
            ],
        ],
        loc="upper left",
        bbox_to_anchor=(0.075, 0.945),
        frameon=False,
        fontsize=10,
    )
    fit = report["fits"]["joint_pv_errors"]
    fig.text(
        0.075,
        0.069,
        f"All 37 P–V–T rows fitted simultaneously with 4 global linear-modulus parameters + 37 latent volumes; ambient V₀(T) retained exactly.\nMinimize Σ[(Pmodel(Vtrue,T)−Pobs)/σP]² + [(Vtrue−Vobs)/σV]²; χ² = {fit['chi_square']:.2f}, dof = 33, reduced χ² = {fit['reduced_chi_square']:.2f}.\nLower panels: signed orthogonal distance in reported P/V error units; shaded bands are ±1 and ±2. Temperature and ambient-law uncertainty excluded.\nIndependent Gaussian coordinate errors assumed; source phase/pressure-scale limitations remain. This diagnostic does not change the catalog.",
        fontsize=9,
        va="center",
        linespacing=1.45,
    )
    stem = ROOT / "docs/data/brucite-joint-pv-errors"
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
        write_plot(report)
    print("Brucite joint P/V-errors report is current.")
    for key, fit in report["fits"].items():
        print(
            key,
            json.dumps(
                {field: fit[field] for field in ["pressure_rmse_gpa", "volume_rmse_a3"]}
            ),
            "profiled chi-square",
            fit["projection"]["chi_square"],
        )
    print(json.dumps(report["fits"]["joint_pv_errors"]["parameters"], indent=2))


if __name__ == "__main__":
    main()
