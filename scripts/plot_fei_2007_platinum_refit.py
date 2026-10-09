"""Plot the published Fei Pt EOS and diagnostic refit against recovered data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from scripts.reproduce_fei_2007_platinum import (
    PARAMETERS,
    load_rows,
    observations,
    source_pressure,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fei-2007-platinum-fit-comparison"
PUBLISHED_COLOR = "#2166ac"
REFIT_COLOR = "#d95f02"


def main():
    result = json.loads(
        (ROOT / "docs/data/fei-2007-platinum-reproduction.json").read_text()
    )
    diagnostic = result["fixed_v0_joint_diagnostic"]
    coefficients = np.array(PARAMETERS["platinum"])
    for i, name in ((2, "K0_prime"), (3, "gamma0"), (4, "q")):
        coefficients[i] = diagnostic["parameters"][name]
    volume, temperature, observed, paired, _ = observations()
    published = source_pressure(volume, temperature)
    refitted = source_pressure(volume, temperature, coefficients=coefficients)
    for prediction, key in (
        (published, "published_rmse_gpa"),
        (refitted, "rmse_gpa"),
    ):
        assert (
            abs(np.sqrt(np.mean((prediction - observed) ** 2)) - diagnostic[key]) < 1e-8
        )
    cold_rows = load_rows("platinum_dewaele_2004_table1_compression")
    volume_error = np.r_[
        [4 * float(row["atomic_volume_uncertainty_a3"]) for row in cold_rows],
        [float(row["volume_uncertainty_a3"]) for row in paired],
    ]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 11,
            "savefig.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(12.5, 10.5), facecolor="white")
    outer = fig.add_gridspec(
        2, 2, left=0.085, right=0.97, bottom=0.14, top=0.83, hspace=0.38, wspace=0.26
    )
    fig.suptitle(
        "Platinum: published Fei EOS and our refit", x=0.5, y=0.985, fontsize=19
    )
    fig.text(
        0.5,
        0.948,
        "Published: K′₀ = 5.080, γ₀ = 2.720, q = +0.500   |   "
        "Refit: K′₀ = 5.08696, γ₀ = 2.61240, q = −1.66064",
        ha="center",
        fontsize=11,
    )
    handles = [
        Line2D([], [], color=PUBLISHED_COLOR, lw=2, label="Fei (2007), Table 1"),
        Line2D([], [], color=REFIT_COLOR, lw=2, ls="--", label="Our joint refit"),
        Line2D([], [], color="#333333", marker="x", ls="none", label="Dewaele RT data"),
        Line2D(
            [],
            [],
            color="#333333",
            marker="o",
            mfc="white",
            ls="none",
            label="Fei (2004) Au–Pt data",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.5,
        0.852,
        "Total pressure RMSE, all 78 observations: published 0.442 GPa → refit 0.347 GPa",
        ha="center",
        fontsize=11,
    )
    for slot, temp in enumerate((300.0, 1473.0, 1673.0, 1873.0)):
        mask = temperature == temp
        inner = outer[slot // 2, slot % 2].subgridspec(
            2, 1, height_ratios=(3.1, 1.3), hspace=0.08
        )
        ax = fig.add_subplot(inner[0])
        residual = fig.add_subplot(inner[1], sharex=ax)
        grid = np.linspace(volume[mask].min() - 0.15, volume[mask].max() + 0.15, 500)
        ax.plot(source_pressure(grid, temp), grid, color=PUBLISHED_COLOR, lw=2)
        ax.plot(
            source_pressure(grid, temp, coefficients=coefficients),
            grid,
            color=REFIT_COLOR,
            ls="--",
            lw=2,
        )
        for dewaele, marker in ((True, "x"), (False, "o")):
            selected = mask & ((np.arange(len(volume)) < 36) == dewaele)
            if selected.any():
                ax.errorbar(
                    observed[selected],
                    volume[selected],
                    yerr=volume_error[selected],
                    fmt=marker,
                    color="#333333",
                    markerfacecolor="white",
                    markersize=4.5,
                    elinewidth=0.7,
                    capsize=2,
                    zorder=3,
                )
        limits = (-1, 100) if temp == 300 else (10, 32)
        ax.set_xlim(limits)
        ax.set_ylim(volume[mask].min() - 0.3, volume[mask].max() + 0.3)
        ax.set_ylabel("Pt cell volume (Å³)")
        ax.set_title(
            f"{temp:,.0f} K  ·  {int(mask.sum())} observations", loc="left", fontsize=12
        )
        ax.tick_params(labelbottom=False)
        ax.grid(alpha=0.14)
        residual.axhline(0, color="#666666", lw=0.8)
        for prediction, color, marker in (
            (published, PUBLISHED_COLOR, "o"),
            (refitted, REFIT_COLOR, "D"),
        ):
            residual.scatter(
                observed[mask],
                prediction[mask] - observed[mask],
                s=22,
                marker=marker,
                edgecolors=color,
                facecolors="none",
                linewidths=1,
            )
        extent = max(
            abs(published[mask] - observed[mask]).max(),
            abs(refitted[mask] - observed[mask]).max(),
        )
        residual.set_ylim(-max(0.6, extent * 1.2), max(0.6, extent * 1.2))
        residual.set_ylabel("ΔP (GPa)", fontsize=10)
        residual.set_xlabel("Calibrated pressure (GPa)")
        residual.grid(alpha=0.14)
    fig.text(
        0.085,
        0.045,
        "ΔP = EOS pressure − calibrated data pressure. Fei data re-reduced with the published Fei (2007) Au scale.\n"
        "Volume error bars only. Equal pressure weights; V₀ = 60.38 Å³, K₀ = 277 GPa, θ₀ = 230 K fixed.\n"
        "Dewaele’s 298 K observations assigned to 300 K. Both curves use the printed Debye-temperature law.",
        fontsize=9,
        color="#444444",
        va="center",
    )
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=220)
    fig.savefig(OUTPUT.with_suffix(".pdf"))
    plt.close(fig)
    with OUTPUT.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "dataset",
                "source_row",
                "temperature_k",
                "volume_a3",
                "volume_uncertainty_a3",
                "calibrated_pressure_gpa",
                "published_pressure_gpa",
                "refit_pressure_gpa",
                "published_residual_gpa",
                "refit_residual_gpa",
            ]
        )
        for i in range(len(volume)):
            writer.writerow(
                [
                    "Dewaele 2004" if i < 36 else "Fei 2004",
                    i + 1 if i < 36 else paired[i - 36]["run"],
                    temperature[i],
                    volume[i],
                    volume_error[i],
                    observed[i],
                    published[i],
                    refitted[i],
                    published[i] - observed[i],
                    refitted[i] - observed[i],
                ]
            )
    print(OUTPUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
