"""Plot the retained Ono–Mibe observations and two existing BM3 curves.

Repository-only plotting utility; no fitting or catalog edits. Run with
``uv run python -m scripts.plot_ono_mibe_2010_fe3c``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from scripts.ono_mibe_2010_fe3c import (
    CARD,
    IDENTIFIER,
    NUMERICAL,
    ROOT,
    observations,
    volume,
)


def plot(output_dir: Path):
    result = json.loads((ROOT / NUMERICAL).read_text())
    diagnostic = result["independent_diagnostic_refit"]
    refit = diagnostic["parameters_K0_K0_prime"]
    source = next(
        r
        for r in json.loads(CARD.read_text())["eos_records"]
        if r["identifier"] == IDENTIFIER
    )
    assert {k: v["value"] for k, v in source["parameters"].items()} == {
        "V0": 155.4,
        "K0": 167,
        "K0_prime": 6.7,
    }
    assert diagnostic["fixed_V0_a3"] == 155.4
    rows = observations()[1]
    p = np.array([r["pressure_gpa"] for r in rows])
    v = np.array([r["volume_a3"] for r in rows])
    sv = np.array([r["volume_a3_error_width"] for r in rows])
    selected = np.array([r["selection_code"] == 1 for r in rows])
    outside = np.array([r["selection_code"] == 2 for r in rows])
    anchor = np.array([r["selection_code"] == 0 for r in rows])
    assert len(rows) == 14 and selected.sum() == 5 and outside.sum() == 8
    assert all(r["pressure_gpa_error_width"] is None for r in rows)
    grid = np.linspace(0, 34.3, 501)
    published_v, refit_v = volume(grid), volume(grid, refit)
    blue, orange, grey = "#2166ac", "#d95f02", "#7a828a"
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "svg.hashsalt": "ono-mibe-2010-pv-v1",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(12.4, 8.2), facecolor="white")
        gs = fig.add_gridspec(
            2,
            2,
            height_ratios=[2.3, 1],
            width_ratios=[1.3, 1],
            hspace=0.42,
            wspace=0.27,
        )
        close, full, residual = (
            fig.add_subplot(gs[0, 0]),
            fig.add_subplot(gs[0, 1]),
            fig.add_subplot(gs[1, :]),
        )
        fig.subplots_adjust(left=0.09, right=0.98, top=0.86, bottom=0.16)
        fig.suptitle(
            "Ono–Mibe (2010) · Fe₃C at 300 K",
            x=0.09,
            y=0.97,
            ha="left",
            fontsize=18,
            fontweight="bold",
        )
        fig.text(
            0.09,
            0.925,
            "Published BM3 and the separate fixed-V₀ diagnostic refit",
            fontsize=12,
            color="#444444",
        )
        fig.text(
            0.09,
            0.89,
            f"V₀ = 155.4 Å³ for both   |   Published: K₀ = 167 GPa, K′₀ = 6.7   |   Refit: K₀ = {refit[0]:.3f} GPa, K′₀ = {refit[1]:.5f}",
            fontsize=10,
        )
        for ax in (close, full):
            ax.plot(
                grid, published_v, color=blue, lw=2.3, label="Published BM3", zorder=2
            )
            ax.plot(
                grid,
                refit_v,
                color=orange,
                lw=2,
                ls="--",
                label="Diagnostic refit",
                zorder=3,
            )
            ax.errorbar(
                p[selected],
                v[selected],
                yerr=sv[selected],
                fmt="o",
                color="#222222",
                ms=5.5,
                capsize=4,
                elinewidth=1.2,
                label="Five nonambient fit points",
                zorder=5,
            )
            ax.errorbar(
                p[anchor],
                v[anchor],
                yerr=sv[anchor],
                fmt="D",
                color="#222222",
                mfc="white",
                ms=6,
                capsize=4,
                elinewidth=1.2,
                label="Measured ambient anchor",
                zorder=5,
            )
            ax.set(xlabel="Pressure (GPa)", ylabel="Unit-cell volume (Å³)")
            ax.grid(alpha=0.2, lw=0.7)
            ax.set_axisbelow(True)
        close.set(
            xlim=(-1, 36),
            ylim=(134.5, 156.7),
            title="Fit range · curves nearly overlap",
        )
        close.legend(loc="upper right", fontsize=9, frameon=False)
        full.errorbar(
            p[outside],
            v[outside],
            yerr=sv[outside],
            fmt="o",
            color=grey,
            mfc="white",
            ms=5.5,
            capsize=4,
            elinewidth=1.2,
            label="Above 35 GPa (outside fit)",
            zorder=5,
        )
        full.axvspan(34.3, 70, color="#e5e7eb", alpha=0.38, zorder=0)
        full.axvline(34.3, color=grey, lw=1, ls=":")
        full.set(
            xlim=(-1.5, 70), ylim=(123.5, 157), title="All 14 reported observations"
        )
        full.legend(
            handles=[full.get_legend_handles_labels()[0][-1]],
            labels=["Above 35 GPa (outside fit)"],
            loc="upper right",
            fontsize=9,
            frameon=False,
        )
        full.text(
            0.98,
            0.74,
            "Curves stop at 34.3 GPa",
            transform=full.transAxes,
            ha="right",
            fontsize=9,
            color="#666666",
        )
        low = selected | anchor
        residual.axhline(0, color=blue, lw=1.6, label="Published BM3 baseline")
        residual.plot(
            grid,
            refit_v - published_v,
            color=orange,
            lw=2,
            ls="--",
            label="Refit − published",
        )
        residual.errorbar(
            p[low],
            v[low] - volume(p[low]),
            yerr=sv[low],
            fmt="o",
            color="#222222",
            capsize=4,
            ms=5,
            elinewidth=1.2,
            label="Observed − published",
            zorder=4,
        )
        residual.set(
            xlim=(-1, 36),
            ylim=(-0.85, 0.85),
            xlabel="Pressure (GPa)",
            ylabel="V − Vpublished (Å³)",
            title="Expanded differences in the fit range",
        )
        residual.grid(alpha=0.2, lw=0.7)
        residual.set_axisbelow(True)
        residual.legend(loc="upper right", fontsize=9, frameon=False, ncol=3)
        fig.text(
            0.09,
            0.085,
            "Error bars: reported Table 1 volume-error widths, drawn literally; confidence level and covariance are unspecified.",
            fontsize=9.5,
            color="#444444",
        )
        fig.text(
            0.09,
            0.06,
            "Pressure errors are not reported. The diagnostic uses equal pressure weights on five nonambient points, with measured V₀ fixed.",
            fontsize=9.5,
            color="#444444",
        )
        fig.text(
            0.09,
            0.035,
            "Source: Ono & Mibe, PEPI 180 (2010), Tables 1 and 3 · doi:10.1016/j.pepi.2010.03.008 · Conventional Pnma cell, Z = 4",
            fontsize=9,
            color="#666666",
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for extension in ("png", "svg", "pdf"):
            path = output_dir / ("ono-mibe-2010-pv-comparison." + extension)
            fig.savefig(path, dpi=180, facecolor="white")
            paths.append(path)
        plt.close(fig)
    print("\n".join(str(p) for p in paths))
    print(
        f"Maximum curve separation within 0–34.3 GPa: {max(abs(refit_v - published_v)):.6f} Å³"
    )
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "output/ono-mibe-2010-fe3c-2026-11-01/pv-comparison",
    )
    plot(parser.parse_args().output_dir)
