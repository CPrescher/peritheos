"""Plot all six Hu shots against the two retained velocity parameterizations.

Run ``uv run python -m scripts.plot_hu_2019_fe3c``. No optimization, source
rewriting, porosity correction or catalog change is performed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from scripts.hu_2019_fe3c import (
    CARD,
    CAVEAT,
    IDENTIFIER,
    NUMERICAL,
    PUBLISHED,
    RHO0,
    ROOT,
    observations,
    shock_velocity,
)


def specific_volume(pressure, coefficients):
    """Full source-fit comparison on the explicitly adopted mean-density path."""
    c0, slope = coefficients
    p = np.asarray(pressure)
    up = 2 * p / (RHO0 * (np.sqrt(c0**2 + 4 * slope * p / RHO0) + c0))
    return (1 - up / shock_velocity(up, coefficients)) / RHO0


def plot(output_dir: Path):
    report = json.loads((ROOT / NUMERICAL).read_text())
    refit = report["independent_diagnostic_refit"]["parameters_c0_s"]
    source = next(
        r
        for r in json.loads(CARD.read_text())["eos_records"]
        if r["identifier"] == IDENTIFIER
    )
    assert [source["parameters"][k]["value"] for k in ("c0", "s")] == PUBLISHED.tolist()
    rows = observations()[1]
    assert len(rows) == 6

    def column(key):
        return np.array([r[key] for r in rows])

    up, us = column("particle_velocity_km_s"), column("shock_velocity_km_s")
    eup, eus = (
        column("particle_velocity_km_s_error_width"),
        column("shock_velocity_km_s_error_width"),
    )
    p, rho = column("pressure_gpa"), column("density_g_cm3")
    ep, erho = column("pressure_gpa_error_width"), column("density_g_cm3_error_width")
    v = 1 / rho
    # Exact transformed density intervals, rather than an invented symmetric V error.
    ev = np.vstack([v - 1 / (rho + erho), 1 / (rho - erho) - v])
    grid = np.linspace(up.min(), up.max(), 501)
    pgrid = np.linspace(p.min(), p.max(), 501)
    predictions = [shock_velocity(up, c) for c in (PUBLISHED, refit)]
    residuals = [predicted - us for predicted in predictions]
    np.testing.assert_allclose(
        residuals[0], report["source_equation_check"]["residuals_km_s"], atol=1e-12
    )
    np.testing.assert_allclose(
        residuals[1],
        report["independent_diagnostic_refit"]["residuals_km_s"],
        atol=1e-12,
    )
    blue, orange, grey = "#2166ac", "#d95f02", "#7a828a"
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "svg.hashsalt": "hu-six-shot-v1",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(13.6, 8.7), facecolor="white")
        gs = fig.add_gridspec(2, 2, height_ratios=[2.3, 1.1], hspace=0.44, wspace=0.28)
        ax, pv, residual, notes = [
            fig.add_subplot(gs[i, j]) for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]
        ]
        fig.subplots_adjust(left=0.08, right=0.98, top=0.84, bottom=0.19)
        fig.suptitle(
            "Hu et al. (2019) · All six Fe₃C shock states",
            x=0.08,
            y=0.97,
            ha="left",
            fontsize=18,
            fontweight="bold",
        )
        fig.text(
            0.08,
            0.925,
            "Published Hugoniot and the separate six-shot diagnostic refit",
            fontsize=12,
            color="#444444",
        )
        fig.text(
            0.08,
            0.887,
            "Both fits use all six shots through 245 GPa. Reported coordinate-error widths are shown without rescaling.",
            fontsize=10,
            color="#444444",
        )
        for coefficients, color, ls, label in [
            (PUBLISHED, blue, "-", "Published Hugoniot"),
            (refit, orange, "--", "Diagnostic refit"),
        ]:
            ax.plot(
                grid,
                shock_velocity(grid, coefficients),
                color=color,
                ls=ls,
                lw=2.2,
                label=label,
                zorder=2,
            )
            pv.plot(
                pgrid,
                specific_volume(pgrid, coefficients),
                color=color,
                ls=ls,
                lw=2.2,
                zorder=2,
            )
        for mask, marker, color, face, label in [
            (
                np.ones(len(rows), dtype=bool),
                "o",
                "#222222",
                "#222222",
                "All six reported shots",
            ),
        ]:
            ax.errorbar(
                up[mask],
                us[mask],
                xerr=eup[mask],
                yerr=eus[mask],
                fmt=marker,
                color=color,
                mfc=face,
                ms=6,
                capsize=4,
                elinewidth=1.1,
                label=label,
                zorder=5,
            )
            pv.errorbar(
                p[mask],
                v[mask],
                xerr=ep[mask],
                yerr=ev[:, mask],
                fmt=marker,
                color=color,
                mfc=face,
                ms=6,
                capsize=4,
                elinewidth=1.1,
                zorder=5,
            )
        for i, shot in enumerate(column("shot")):
            ax.annotate(
                str(int(shot)),
                (up[i], us[i]),
                xytext=(-9, 11),
                textcoords="offset points",
                fontsize=9,
                color="#444444",
            )
            pv.annotate(
                str(int(shot)),
                (p[i], v[i]),
                xytext=(6, 8),
                textcoords="offset points",
                fontsize=9,
                color="#444444",
            )
        ax.set(
            xlim=(1.44, 3.43),
            ylim=(6.95, 10.28),
            xlabel="Particle velocity, uₚ (km/s)",
            ylabel="Shock velocity, Uₛ (km/s)",
            title="Direct velocity-fit comparison",
        )
        ax.legend(loc="upper left", frameon=False, fontsize=9)
        pv.set(
            xlim=(75, 255),
            ylim=(0.0875, 0.1052),
            xlabel="Pressure (GPa)",
            ylabel="Specific volume (cm³/g)",
            title="Pressure–volume comparison",
        )
        for r, color, marker, label in [
            (residuals[0], blue, "o", "Published"),
            (residuals[1], orange, "s", "Diagnostic"),
        ]:
            residual.scatter(
                up, r, color=color, marker=marker, s=26, label=label, zorder=4
            )
        residual.axhline(0, color=grey, lw=1)
        residual.set(
            xlim=(1.44, 3.43),
            ylim=(-0.085, 0.065),
            xlabel="Particle velocity, uₚ (km/s)",
            ylabel="Uₛ(model) − Uₛ(data) (km/s)",
            title="Velocity residuals at the printed uₚ",
        )
        residual.legend(loc="lower left", frameon=False, fontsize=9, ncol=2)
        for pane in (ax, pv, residual):
            pane.grid(alpha=0.19, lw=0.7)
            pane.set_axisbelow(True)
        notes.axis("off")
        notes.text(
            0,
            1,
            "Velocity law: Uₛ = c₀ + s uₚ",
            va="top",
            fontsize=11,
            fontweight="bold",
        )
        table = notes.table(
            cellText=[
                ["c₀ (km/s)", "4.627", f"{refit[0]:.6f}"],
                ["s", "1.614", f"{refit[1]:.6f}"],
                ["Uₛ RMS (km/s)", "0.03468", "0.03138"],
            ],
            colLabels=["", "Published", "Diagnostic"],
            cellLoc="left",
            colLoc="left",
            bbox=[0, 0.2, 0.98, 0.62],
            colWidths=[0.43, 0.27, 0.3],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        for (row, _), cell in table.get_celld().items():
            cell.set_edgecolor("#e5e7eb")
            if row == 0:
                cell.set_facecolor("#f3f4f6")
                cell.set_text_props(fontweight="bold")
        notes.text(
            0,
            0.06,
            "Diagnostic: equal vertical Uₛ weights; no new fitting in this plot.",
            fontsize=9.5,
            va="bottom",
        )
        fig.text(
            0.08,
            0.125,
            "Bars: source Uₛ/uₚ and P/ρ uncertainty widths; confidence level and covariance are unspecified. Residual panel shows central values only.",
            fontsize=9.2,
            color="#444444",
        )
        fig.text(
            0.08,
            0.099,
            "P–V points use 1/ρ with exact transformed density-error intervals. Curves use adopted mean precursor ρ₀ = 7.525 g/cm³; porosity is retained.",
            fontsize=9.2,
            color="#444444",
        )
        fig.text(
            0.08,
            0.073,
            "Full experimental Hugoniot: 85.6–245 GPa. " + CAVEAT,
            fontsize=9.2,
            color="#444444",
        )
        fig.text(
            0.08,
            0.042,
            "Source: Hu et al., GRL 46 (2019), Table 1 · doi:10.1029/2019GL084545 · Data: Fei (2019), PANGAEA.905474, CC-BY-4.0",
            fontsize=9,
            color="#666666",
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for extension in ("png", "svg", "pdf"):
            path = output_dir / ("hu-2019-six-shot-comparison." + extension)
            fig.savefig(path, dpi=180, facecolor="white")
            paths.append(path)
        plt.close(fig)
    print("\n".join(str(p) for p in paths))
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "output/hu-2019-fe3c-full-range-2026-11-02/comparison",
    )
    plot(parser.parse_args().output_dir)
