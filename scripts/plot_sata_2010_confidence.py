"""Conditional K-Kprime ellipses and mean P(V) bands for both Sata diagnostics.

Use stored original/recalculated-pressure fit covariances; no refitting or
catalog writes. Joint pair regions use 2*F(2,nu); pointwise mean bands use t(nu).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import f as f_distribution
from scipy.stats import t as t_distribution

from scripts.audit_sata_2010_fe3c import ROOT, VR, csv_text, modified_pressure, pressure
from scripts.refit_sata_2010_recalculated_pressure import REPORT as FIT_REPORT
from scripts.report_plot_themes import dark_path, dark_svg

REPORT = ROOT / "docs/data/sata-2010-fe3c-confidence.json"
ELLIPSES = ROOT / "docs/data/sata-2010-fe3c-confidence-ellipses.csv"
BANDS = ROOT / "docs/data/sata-2010-fe3c-confidence-bands.csv"
ASSETS = ROOT / "docs/literature-reproductions/assets"
LEVELS = {"68": 0.6826894921370859, "95": 0.95}
VERSIONS = ("original_pressures", "recalculated_Sata_unchanged_Li")


def fit_state(report, version, model):
    if model == "fixed_V0_BM3":
        source = report[model][version]
        return (
            np.array(source["parameters"]),
            np.array(source["conditional_covariance"]),
            source["degrees_of_freedom"],
        )
    source = report["reference_mBM3_check"][version]
    return (
        np.array(source["parameters_Pr_Kr_Kprime"]),
        np.array(source["covariance"]),
        source["dof"],
    )


def mean_and_gradient(volume, parameters, model):
    """Analytic physical-parameter gradients, including all three mBM3 terms."""
    volume = np.atleast_1d(volume).astype(float)
    f = ((VR / volume) ** (2 / 3) - 1) / 2
    factor = (1 + 2 * f) ** 2.5
    if model == "fixed_V0_BM3":
        k, kp = parameters
        mean = pressure(volume, parameters)
        gradient = np.column_stack([mean / k, factor * 4.5 * k * f**2])
    else:
        _, k, kp = parameters
        mean = modified_pressure(volume, parameters)
        gradient = factor[:, None] * np.column_stack(
            [
                1 - 5 * f + 17.5 * f**2,
                3 * f + 4.5 * (kp - 4) * f**2,
                4.5 * k * f**2,
            ]
        )
    return mean, gradient


def ellipse_points(center, covariance, level, dof):
    theta = np.linspace(0, 2 * np.pi, 361)
    radius_squared = 2 * f_distribution.ppf(level, 2, dof)
    unit = np.column_stack([np.cos(theta), np.sin(theta)])
    points = (
        np.asarray(center)
        + np.sqrt(radius_squared) * unit @ np.linalg.cholesky(covariance).T
    )
    return theta, points, float(radius_squared)


def build():
    source = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    measured = np.array([r["volume_a3"] for r in source["inputs"]])
    volume = np.linspace(min(measured), max(measured), 241)
    summary, ellipse_rows, band_rows = {}, [], []
    for version in VERSIONS:
        summary[version] = {}
        for model in ("fixed_V0_BM3", "reference_mBM3"):
            parameters, covariance, dof = fit_state(source, version, model)
            pair = np.array([0, 1] if model == "fixed_V0_BM3" else [1, 2])
            pair_cov = covariance[np.ix_(pair, pair)]
            correlation = pair_cov[0, 1] / np.sqrt(pair_cov[0, 0] * pair_cov[1, 1])
            mean, gradient = mean_and_gradient(volume, parameters, model)
            variance = np.einsum("ij,jk,ik->i", gradient, covariance, gradient)
            if min(variance) < 0:
                raise ValueError("Negative propagated variance")
            se = np.sqrt(variance)
            entry = {
                "parameter_order": ["K0_gpa", "K0_prime"]
                if model == "fixed_V0_BM3"
                else ["Pr_gpa", "Kr_gpa", "Kr_prime"],
                "parameters": parameters.tolist(),
                "full_covariance": covariance.tolist(),
                "pair_covariance": pair_cov.tolist(),
                "K_Kprime_correlation": float(correlation),
                "residual_degrees_of_freedom": dof,
                "levels": {},
            }
            rows = [
                {
                    "pressure_version": version,
                    "model": model,
                    "volume_a3": float(v),
                    "mean_pressure_gpa": float(p),
                    "mean_standard_error_pressure_gpa": float(s),
                }
                for v, p, s in zip(volume, mean, se)
            ]
            for name, level in LEVELS.items():
                angles, points, threshold = ellipse_points(
                    parameters[pair], pair_cov, level, dof
                )
                t_factor = float(t_distribution.ppf((1 + level) / 2, dof))
                for angle, point in zip(angles, points):
                    ellipse_rows.append(
                        {
                            "pressure_version": version,
                            "model": model,
                            "level": level,
                            "angle_rad": float(angle),
                            "K_gpa": float(point[0]),
                            "K_prime": float(point[1]),
                        }
                    )
                for row, half_width in zip(rows, t_factor * se):
                    row[f"lower_{name}_pressure_gpa"] = row[
                        "mean_pressure_gpa"
                    ] - float(half_width)
                    row[f"upper_{name}_pressure_gpa"] = row[
                        "mean_pressure_gpa"
                    ] + float(half_width)
                entry["levels"][name] = {
                    "coverage": level,
                    "joint_pair_mahalanobis_squared": threshold,
                    "pointwise_mean_t_factor": t_factor,
                    "mean_pressure_half_width_range_gpa": [
                        float(min(t_factor * se)),
                        float(max(t_factor * se)),
                    ],
                }
            band_rows.extend(rows)
            summary[version][model] = entry
    report = {
        "kind": "conditional_diagnostic_confidence_visualization",
        "source_refit_report": FIT_REPORT.relative_to(ROOT).as_posix(),
        "source_sha256": hashlib.sha256(FIT_REPORT.read_bytes()).hexdigest(),
        "temperature_k": 300,
        "volume_range_a3": [float(min(measured)), float(max(measured))],
        "reference_volume_cell_a3": VR,
        "assumptions": "Existing equal-pressure-weight least squares, independent Gaussian pressure residuals of common unknown variance, RSS/(N-p) covariance scale, fixed measured volumes and fixed V0/Vr. These assumptions define conditional diagnostics; their experimental validity and the source's error convention are not established.",
        "ellipse_convention": "Local covariance ellipses for the K/Kprime pair: (theta_pair-center)^T C_pair^-1 (theta_pair-center) <= 2 F_(2,nu)(coverage). mBM3 uses the Kr/Krprime submatrix from the full fit with Pr free, not the covariance conditional on fixed Pr. Physical Kprime is a coefficient ratio, so these are locally linearized ellipses.",
        "band_convention": "Pointwise confidence intervals for the fitted mean P(V): P_hat(V) +/- t_nu((1+coverage)/2) sqrt(J C J^T). Displayed as horizontal P intervals with pressure on x and volume on y. mBM3 propagates full 3x3 covariance including Pr cross terms. Pair-joint ellipse coverage and pointwise mean coverage are different probability statements.",
        "excluded": "No residual variance term is added for future measurements: these are not prediction intervals or simultaneous curve bands. No V0, measured-volume, pressure-scale/coefficient or systematic uncertainty is added. No source-published covariance or one-/two-sigma convention is invented. Existing catalog, observations, refits, parity and recommendations are unchanged.",
        "method_reference": "https://itl.nist.gov/div898/handbook/pmd/section5/pmd511.htm",
        "linearized_region_reference": "https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir89-4102.pdf#page=9",
        "fits": summary,
    }
    return report, ellipse_rows, band_rows


def plot(report, ellipse_rows, band_rows, preview_dir=None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    source = json.loads(FIT_REPORT.read_text(encoding="utf-8"))
    volumes = np.array([r["volume_a3"] for r in source["inputs"]])
    sata = np.array([r["source"] == 0 for r in source["inputs"]])
    colors = ["#245b8a", "#ae4c14"]
    style = {
        "svg.hashsalt": "sata-confidence-v1",
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 10,
    }
    paths = []
    with plt.rc_context(style):
        for version in VERSIONS:
            original = version == "original_pressures"
            name = "original" if original else "recalculated"
            coordinate = np.array(
                [
                    r["published_pressure_gpa" if original else "fit_pressure_gpa"]
                    for r in source["inputs"]
                ]
            )
            fig = plt.figure(figsize=(14.6, 5.4), layout="constrained")
            grid = fig.add_gridspec(1, 3, width_ratios=[1.12, 1, 1])
            ellipse_ax = fig.add_subplot(grid[0, 0])
            pv_axes = [fig.add_subplot(grid[0, 1])]
            pv_axes.append(
                fig.add_subplot(grid[0, 2], sharex=pv_axes[0], sharey=pv_axes[0])
            )
            title = (
                "Original Sata + unchanged Li"
                if original
                else "Recalculated Sata + unchanged Li"
            )
            fig.suptitle(f"Fe₃C at 300 K — {title}", fontsize=14)
            for model, color, pv, label in zip(
                ("fixed_V0_BM3", "reference_mBM3"),
                colors,
                pv_axes,
                ("Fixed-V₀ BM3", "mBM3, Pᵣ fitted"),
            ):
                entry = report["fits"][version][model]
                for level_name, linestyle in (("95", "--"), ("68", "-")):
                    points = np.array(
                        [
                            [r["K_gpa"], r["K_prime"]]
                            for r in ellipse_rows
                            if r["pressure_version"] == version
                            and r["model"] == model
                            and r["level"] == LEVELS[level_name]
                        ]
                    )
                    ellipse_ax.plot(
                        points[:, 0],
                        points[:, 1],
                        color=color,
                        ls=linestyle,
                        label=f"{label}: {'68.3' if level_name == '68' else '95'}%",
                    )
                    if level_name == "68":
                        ellipse_ax.fill(
                            points[:, 0], points[:, 1], color=color, alpha=0.12
                        )
                pair = np.array(entry["parameters"])[
                    [0, 1] if model == "fixed_V0_BM3" else [1, 2]
                ]
                ellipse_ax.scatter(*pair, marker="+", s=70, color=color, zorder=5)
                curve = [
                    r
                    for r in band_rows
                    if r["pressure_version"] == version and r["model"] == model
                ]
                v = np.array([r["volume_a3"] for r in curve])
                mean = np.array([r["mean_pressure_gpa"] for r in curve])
                for level_name, alpha in (("95", 0.16), ("68", 0.36)):
                    pv.fill_betweenx(
                        v,
                        [r[f"lower_{level_name}_pressure_gpa"] for r in curve],
                        [r[f"upper_{level_name}_pressure_gpa"] for r in curve],
                        color=color,
                        alpha=alpha,
                        label=f"{'68.3' if level_name == '68' else '95'}% pointwise mean",
                    )
                pv.plot(mean, v, color=color, lw=1.4)
                pv.scatter(
                    coordinate[sata],
                    volumes[sata],
                    s=19,
                    color="#4e5660",
                    zorder=4,
                    label="Sata states",
                )
                pv.scatter(
                    coordinate[~sata],
                    volumes[~sata],
                    marker="D",
                    s=23,
                    facecolors="none",
                    edgecolors="#4e5660",
                    zorder=4,
                    label="Li states (unchanged)",
                )
                pv.set(title=label, xlabel="Pressure (GPa)")
                pv.legend(loc="upper right", fontsize=8, frameon=False)
                half_widths = entry["levels"]["95"][
                    "mean_pressure_half_width_range_gpa"
                ]
                pv.text(
                    0.04,
                    0.04,
                    f"95% P half-width: {half_widths[0]:.2f}–{half_widths[1]:.2f} GPa",
                    transform=pv.transAxes,
                    fontsize=8,
                )
            ellipse_ax.set(
                title="Joint K–K′ confidence ellipses",
                xlabel="K₀ (BM3) or Kᵣ (mBM3), GPa",
                ylabel="K₀′ (BM3) or Kᵣ′ (mBM3)",
            )
            ellipse_ax.legend(loc="upper right", fontsize=7.5, frameon=False)
            correlation = [
                report["fits"][version][m]["K_Kprime_correlation"]
                for m in ("fixed_V0_BM3", "reference_mBM3")
            ]
            ellipse_ax.text(
                0.04,
                0.04,
                f"K–K′ correlation: BM3 {correlation[0]:.3f}; mBM3 {correlation[1]:.3f}",
                transform=ellipse_ax.transAxes,
                fontsize=8,
            )
            pv_axes[0].set_ylabel("Volume (Å³ / conventional cell)")
            pv_axes[1].tick_params(labelleft=False)
            for ax in (ellipse_ax, *pv_axes):
                ax.grid(color="#e7eaee")
                ax.set_axisbelow(True)
                ax.spines[["top", "right"]].set_visible(False)
            fig.text(
                0.5,
                0.025,
                "Conditional diagnostic 68.3% / 95% regions: pair-joint F ellipses; pointwise t mean P(V) bands (horizontal). Not prediction limits.\nV₀/Vᵣ = 149.456 Å³/cell; mBM3 bands include full Pᵣ covariance. Calibration/systematic/input errors excluded; catalog findings unchanged.",
                ha="center",
                fontsize=8.5,
            )
            fig.get_layout_engine().set(rect=(0, 0.11, 1, 0.89))
            path = ASSETS / f"sata-2010-fe3c-confidence-{name}.svg"
            fig.savefig(
                path,
                metadata={
                    "Date": None,
                    "Description": report["band_convention"] + " " + report["excluded"],
                },
            )
            path.with_name(Path(dark_path(path.name)).name).write_bytes(
                dark_svg(path.read_bytes())
            )
            if preview_dir:
                preview_dir.mkdir(parents=True, exist_ok=True)
                fig.savefig(preview_dir / f"confidence-{name}.png", dpi=160)
            plt.close(fig)
            paths.append(path)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--preview-dir", type=Path)
    args = parser.parse_args()
    report, ellipses, bands = build()
    outputs = {
        REPORT: json.dumps(report, indent=2) + "\n",
        ELLIPSES: csv_text(ellipses),
        BANDS: csv_text(bands),
    }
    if args.write:
        for path, text in outputs.items():
            path.write_text(text)
    elif args.check:
        for path, text in outputs.items():
            if path.read_text(encoding="utf-8") != text:
                raise SystemExit(f"Confidence diagnostic needs regeneration: {path}")
        print("Sata conditional confidence report and plot-data checks passed")
    elif not args.plot:
        print(outputs[REPORT])
    if args.plot:
        for path in plot(report, ellipses, bands, args.preview_dir):
            print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
