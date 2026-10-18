"""Compare ambient-only expansion fits with the joint brucite EOS fits.

This reads the saved joint-fit report, checks its source hashes, and performs
only a conditional ambient MGD thermal refit. No source/catalog mutation.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
from scipy.optimize import least_squares

from scripts.brucite_curved_expansion import reproduce as curved_report
from scripts.brucite_mgd_comparison import (
    NAMES,
    data,
    volumes,
)
from scripts.brucite_mgd_comparison import (
    REPORT as JOINT_REPORT,
)
from scripts.brucite_reproduce import ROOT, reference_volume

REPORT = ROOT / "docs/data/brucite-ambient-comparison.json"
FREE_INDICES = np.array([0, 3, 5])
FREE_SCALES = np.array([1, 1, 1000])


def conditional_ambient_mgd(joint_parameters, dataset):
    """Optimize V300, gamma0, theta0 on Table 1; freeze joint K,Kprime,q."""
    template = np.array(joint_parameters, dtype=float)
    temperature = dataset["temperature"][37:]
    observed = dataset["volume"][37:]

    def expand(x):
        result = template.copy()
        result[FREE_INDICES] = x * FREE_SCALES
        return result

    def residual(x):
        return (
            volumes("mgd", np.zeros_like(temperature), temperature, expand(x), observed)
            - observed
        )

    results = [
        least_squares(
            residual,
            [40.93, 0.6, theta],
            bounds=([40, 0.01, 0.01], [42, 3, 4]),
            x_scale="jac",
            ftol=2e-12,
            xtol=2e-12,
            gtol=2e-12,
            max_nfev=1000,
        )
        for theta in [0.1, 0.6, 1.1, 2.4]
    ]
    fit = min(results, key=lambda result: result.fun @ result.fun)
    parameters = expand(fit.x)
    covariance = (
        np.linalg.pinv(fit.jac.T @ fit.jac)
        * (fit.fun @ fit.fun)
        / (len(observed) - len(FREE_INDICES))
    )
    covariance *= FREE_SCALES[:, None] * FREE_SCALES[None, :]
    return {
        "parameters": dict(zip(NAMES["mgd"], parameters.tolist())),
        "free_parameters": [NAMES["mgd"][i] for i in FREE_INDICES],
        "fixed_parameters_from_joint_volume_fit": {
            NAMES["mgd"][i]: float(parameters[i]) for i in [1, 2, 4]
        },
        "source_rows": list(range(1, 12)),
        "objective": "Unweighted sum of squared Vmodel(P=0,T)-Vobserved residuals, all 11 ambient Table 1 rows, in A^6; only V300,gamma0,theta0 optimized. K300,Kprime300,q fixed to the joint volume-objective MGD fit.",
        "objective_sse": float(fit.fun @ fit.fun),
        "solver_success": bool(fit.success),
        "multistart_sse": [float(result.fun @ result.fun) for result in results],
        "active_bounds": {
            NAMES["mgd"][i]: int(mask) if mask else -1
            for i, mask in zip(FREE_INDICES, fit.active_mask)
            if mask or (i == 5 and np.isclose(parameters[i], 10.0, rtol=0, atol=1e-5))
        },
        "bounds_for_free_parameters": {
            "V0_a3": [40, 42],
            "gamma0": [0.01, 3],
            "theta0_k": [10, 4000],
        },
        "formal_local_covariance": {
            "parameter_order": [NAMES["mgd"][i] for i in FREE_INDICES],
            "matrix": covariance.tolist(),
        },
        "uncertainty_qualification": "Formal local SSE/(11-3)*pinv(J.T J), scaled to physical parameters. The theta lower-bound optimum makes symmetric local intervals unreliable. Fixed mechanical/q coefficients come from the joint fit and are treated as exact; their uncertainty and all inter-fit covariance are omitted. This is a conditional thermal refit, not an independent six-parameter ambient determination or a physical theta measurement.",
    }


def curve_values(key, temperature, report):
    t = np.asarray(temperature, dtype=float)
    if key == "published":
        return reference_volume(t)
    if key == "empirical_ambient":
        fit = report["ambient_only_empirical"]
        dt = t - 300
        return fit["V300_a3"] * np.exp(
            fit["alpha_300_per_k"] * dt + 0.5 * fit["dalpha_dT_per_k2"] * dt**2
        )
    model = "empirical" if key == "empirical_joint" else "mgd"
    parameters = list(report["curves"][key]["parameters"].values())
    seed = reference_volume(t)
    return volumes(model, np.zeros_like(t), t, parameters, seed)


def reproduce():
    joint = json.loads(JOINT_REPORT.read_text(encoding="utf-8"))
    for name, checksum in joint["source_file_sha256"].items():
        if (
            hashlib.sha256(
                (ROOT / "peritheos/data/datasets" / name).read_bytes()
            ).hexdigest()
            != checksum
        ):
            raise ValueError("Joint-fit source checksum is stale")
    dataset = data()
    ambient_fit = curved_report()["fit"]
    mgd_joint = joint["fits"]["mgd_volume"]
    conditional = conditional_ambient_mgd(
        list(mgd_joint["parameters"].values()), dataset
    )
    report = {
        "format": "peritheos.brucite-ambient-comparison",
        "format_version": 1,
        "joint_report": "docs/data/brucite-mgd-comparison.json",
        "joint_report_sha256": hashlib.sha256(JOINT_REPORT.read_bytes()).hexdigest(),
        "source_file_sha256": joint["source_file_sha256"],
        "ambient_observations": 11,
        "ambient_temperature_range_k": [
            float(dataset["temperature"][37:].min()),
            float(dataset["temperature"][37:].max()),
        ],
        "high_pressure_observations_already_in_joint_fits": {
            str(t): int(np.sum(dataset["temperature"][:37] == t))
            for t in [300, 473, 673, 873]
        },
        "ambient_only_empirical": ambient_fit,
        "conditional_ambient_mgd": conditional,
        "curves": {
            "published": {
                "definition": "Printed constant-alpha ambient law: 40.917*exp(6.69e-5*(T-300))"
            },
            "empirical_ambient": {
                "definition": "Previous purple ambient-only curved expansion fit: free V300,alpha300,beta using all 11 Table 1 rows; no HP data in this regression."
            },
            "empirical_joint": {
                "parameters": joint["fits"]["empirical_volume"]["parameters"],
                "definition": "Green seven-parameter joint unweighted volume fit to all 48 rows.",
            },
            "mgd_joint": {
                "parameters": mgd_joint["parameters"],
                "definition": "Blue six-parameter joint unweighted volume fit to all 48 rows.",
            },
            "mgd_ambient_conditional": {
                "parameters": conditional["parameters"],
                "definition": "Ambient-only optimization of V300,gamma0,theta0 with K300,Kprime300,q held at the joint-fit values; conditional, not independent of HP-derived mechanical parameters.",
            },
        },
        "qualifications": [
            "Both joint models already use the higher-temperature pressure observations up to 873 K; they constrain thermal behaviour at their measured pressures.",
            "V(P=0,T) above 592.9 K is unmeasured. HP data at the same T constrain it indirectly via pressure extrapolation and the assumed EOS. It is not a measured high-temperature ambient volume.",
            "Fukui reports ambient dehydration at 600.6 K. Zero-pressure brucite curves at higher T are formal model reference volumes, not stable-phase expansion measurements.",
            "Conditional ambient-only MGD reaches the 10 K theta lower bound. This signals poor identifiability under its fixed compression/q choices; it does not establish a physical Debye temperature.",
            "Same volume-objective joint coefficients as the previous P-V figure; no new joint refit or silent coefficient change. Catalog records and source observations unchanged.",
        ],
    }
    t = dataset["temperature"][37:]
    v = dataset["volume"][37:]
    for key, curve in report["curves"].items():
        predicted = curve_values(key, t, report)
        curve.update(
            ambient_predicted_volumes_a3=predicted.tolist(),
            ambient_volume_residuals_a3=(predicted - v).tolist(),
            ambient_volume_rmse_a3=float(np.sqrt(np.mean((predicted - v) ** 2))),
            formal_zero_pressure_volumes_a3={
                str(temp): float(curve_values(key, np.array([temp]), report)[0])
                for temp in [300, 473, 673, 873]
            },
        )
    return report


def write_plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    dataset = data()
    temperature, observed, error = (
        dataset[key][37:] for key in ["temperature", "volume", "volume_error"]
    )
    measured_max = report["ambient_temperature_range_k"][1]
    styles = {
        "published": ("#c64c40", "-", "Published ambient expansion law"),
        "empirical_joint": ("#26875d", "-", "Green: joint empirical fit (48 rows)"),
        "empirical_ambient": (
            "#26875d",
            "--",
            "Empirical: ambient-only fit (11 rows; previous purple law)",
        ),
        "mgd_joint": ("#187ca6", "-", "MGD: joint fit (48 rows)"),
        "mgd_ambient_conditional": (
            "#187ca6",
            "--",
            "MGD: ambient-only thermal refit (fixed joint K, K′, q)",
        ),
    }
    fig = plt.figure(figsize=(14, 9.8), facecolor="white")
    outer = fig.add_gridspec(
        1, 2, left=0.075, right=0.98, bottom=0.18, top=0.79, wspace=0.22
    )
    left = outer[0].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.07)
    right = outer[1].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.07)
    ambient = fig.add_subplot(left[0])
    residual = fig.add_subplot(left[1], sharex=ambient)
    extension = fig.add_subplot(right[0])
    coverage = fig.add_subplot(right[1], sharex=extension)
    measured_grid = np.linspace(295, measured_max, 200)
    extrapolated_grid = np.linspace(measured_max, 873, 200)
    for key, (color, style, _) in styles.items():
        ambient.plot(
            measured_grid,
            curve_values(key, measured_grid, report),
            color=color,
            ls=style,
            lw=1.8,
        )
        extension.plot(
            measured_grid,
            curve_values(key, measured_grid, report),
            color=color,
            ls=style,
            lw=1.8,
        )
        extension.plot(
            extrapolated_grid,
            curve_values(key, extrapolated_grid, report),
            color=color,
            ls=style,
            lw=1.8,
            alpha=0.8,
        )
        residual.plot(
            temperature,
            report["curves"][key]["ambient_volume_residuals_a3"],
            "o" if style == "-" else "s",
            ms=4,
            color=color,
            markerfacecolor=color if style == "-" else "white",
        )
    for axis in [ambient, extension]:
        axis.errorbar(
            temperature,
            observed,
            yerr=error,
            fmt="o",
            color="#242424",
            capsize=2,
            ms=4,
            zorder=8,
        )
        axis.set_ylabel("Unit-cell volume (Å³)")
        axis.tick_params(labelbottom=False)
    residual.errorbar(
        temperature,
        np.zeros_like(temperature),
        yerr=error,
        fmt="none",
        ecolor="#b8b8b8",
        capsize=2,
        zorder=0,
    )
    residual.axhline(0, color="#555555", lw=0.8)
    residual.set_ylabel("ΔV (Å³)")
    residual.set_xlabel("Temperature (K)")
    residual.set_ylim(-0.105, 0.105)
    ambient.set_xlim(290, 605)
    ambient.set_title("Fit to the measured ambient volumes", loc="left", fontsize=12)
    extension.set_title(
        "Formal zero-pressure reference volumes to 873 K", loc="left", fontsize=12
    )
    extension.axvspan(measured_max, 885, color="#eeeeee", zorder=-5)
    extension.axvline(measured_max, color="#999999", ls=":", lw=1)
    extension.text(
        0.97,
        0.05,
        "Shaded: no ambient volume measurements\nCurves are model extrapolations",
        transform=extension.transAxes,
        ha="right",
        fontsize=9,
    )
    extension.set_xlim(290, 885)
    coverage.hlines(0.75, temperature.min(), temperature.max(), color="#242424", lw=1)
    coverage.plot(
        temperature, np.ones_like(temperature) * 0.75, "|", color="#242424", ms=11
    )
    hp_temperatures = [300, 473, 673, 873]
    coverage.plot(hp_temperatures, [0.25] * 4, "|", color="#c78321", ms=14)
    for temp, count in zip(hp_temperatures, [13, 8, 8, 8]):
        coverage.annotate(
            str(count),
            (temp, 0.25),
            xytext=(0, -15),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            color="#9b671e",
        )
    coverage.text(
        0.97,
        0.82,
        "Ambient observations: 11",
        transform=coverage.transAxes,
        ha="right",
        fontsize=9,
    )
    coverage.text(
        0.5,
        0.37,
        "Pressure observations: 37 (counts below)",
        transform=coverage.transAxes,
        ha="center",
        fontsize=9,
        color="#9b671e",
    )
    coverage.set_ylim(-0.1, 1)
    coverage.set_yticks([])
    coverage.set_xticks([300, 473, 593, 673, 873])
    coverage.set_xlabel("Temperature (K)")
    for axis in [ambient, residual, extension, coverage]:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.16)
    fig.suptitle(
        "Brucite: ambient-only expansion fits versus joint P–V–T fits",
        y=0.97,
        fontsize=17,
    )
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            color="#242424",
            ls="none",
            label="Ambient Table 1 observations + reported errors",
        )
    ]
    handles.extend(
        Line2D([], [], color=color, ls=style, lw=2, label=label)
        for color, style, label in styles.values()
    )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.075, 0.935),
        frameon=False,
        fontsize=9.5,
        ncol=2,
    )
    fig.text(
        0.075,
        0.075,
        "Solid green/blue: the same joint-fit parameters used in the previous P–V comparison; both include all 11 ambient + 37 pressure rows.\nDashed: optimize only ambient-volume residuals; MGD holds joint K300, K′300 and q fixed and reaches the theta lower bound (10 K).\nΔV = calculated − observed ambient volume; grey bars show reported errors. Our refits use unweighted volume residuals; red retains published coefficients.\nAbove 592.9 K ambient V is unmeasured; Fukui reports ambient dehydration at 600.6 K. High-pressure data do not establish stable ambient brucite there.",
        fontsize=9,
        va="center",
        linespacing=1.5,
    )
    stem = ROOT / "docs/data/brucite-ambient-comparison"
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
            raise SystemExit("Brucite ambient comparison is stale")
    else:
        REPORT.write_text(payload, encoding="utf-8")
    if args.plot:
        write_plot(report)
    print("Brucite ambient comparison is current.")
    for key, curve in report["curves"].items():
        print(
            key,
            curve["ambient_volume_rmse_a3"],
            curve["formal_zero_pressure_volumes_a3"],
        )


if __name__ == "__main__":
    main()
