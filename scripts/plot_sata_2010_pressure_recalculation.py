"""Show original and recalculated Sata Fe3C coordinates without refitting.

Repository-only figure generation. Recalculated P uses measured MgO volumes
and the printed Speziale BM3 coefficients; original observations remain intact.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from scripts.audit_sata_2010_fe3c import pressure
from scripts.audit_sata_2010_mgo_pressure import ROOT, recalculated_sample_points
from scripts.report_plot_themes import dark_path, dark_svg

FIGURE = (
    ROOT
    / "docs/literature-reproductions/assets/sata-2010-fe3c-pressure-recalculation.svg"
)


def plot(preview: Path | None = None):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    rows = recalculated_sample_points()
    arrays = {name: np.array([r[name] for r in rows]) for name in rows[0]}
    p = arrays["published_pressure_gpa"]
    rec = arrays["recalculated_pressure_gpa"]
    v = arrays["fe3c_volume_a3"]
    source_color, recalculated_color, curve_color = "#4e5660", "#ae4c14", "#245b8a"
    style = {
        "svg.hashsalt": "sata-pressure-recalculation-v1",
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 9,
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(11.8, 6.4), layout="constrained")
        grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1])
        pv = fig.add_subplot(grid[:, 0])
        delta = fig.add_subplot(grid[0, 1])
        residual = fig.add_subplot(grid[1, 1])
        fig.suptitle(
            "Sata (2010) Fe₃C: published and recalculated pressure coordinates",
            fontsize=13,
        )
        sample_v = np.linspace(min(v), max(v), 250)
        pv.plot(
            pressure(sample_v),
            sample_v,
            color=curve_color,
            label="Published Fe₃C BM3 curve",
        )
        pv.errorbar(
            p,
            v,
            xerr=arrays["published_pressure_reported_width_gpa"],
            yerr=arrays["fe3c_volume_reported_width_a3"],
            fmt="o",
            ms=6,
            mfc="none",
            mec=source_color,
            ecolor=source_color,
            elinewidth=0.7,
            capsize=2,
            label="Published P and V (reported error widths)",
        )
        pv.scatter(
            rec,
            v,
            marker="^",
            s=30,
            color=recalculated_color,
            label="Recalculated P; same measured Fe₃C V",
            zorder=4,
        )
        pv.set(
            xlabel="Pressure (GPa)",
            ylabel="Fe₃C volume (Å³ / conventional cell)",
            title="Nine Sata states at 300 K; no refit",
        )
        pv.legend(loc="upper right", fontsize=8, frameon=False)
        delta.scatter(p, rec - p, s=28, color=recalculated_color)
        delta.axhline(0, color=source_color, lw=0.8)
        delta.set(
            xlabel="Published pressure (GPa)",
            ylabel="Recalculated P − published P (GPa)",
            title="Pressure shifts from the printed MgO scale",
        )
        residual.scatter(
            p,
            arrays["published_pressure_residual_gpa"],
            marker="o",
            s=34,
            facecolors="none",
            edgecolors=source_color,
            label="Original pressure coordinates",
        )
        residual.scatter(
            rec,
            arrays["recalculated_pressure_residual_gpa"],
            marker="^",
            s=30,
            color=recalculated_color,
            label="Recalculated pressure coordinates",
        )
        for x1, x2, y1, y2 in zip(
            p,
            rec,
            arrays["published_pressure_residual_gpa"],
            arrays["recalculated_pressure_residual_gpa"],
        ):
            residual.plot([x1, x2], [y1, y2], color=source_color, lw=0.7, zorder=0)
        residual.axhline(0, color=source_color, lw=0.8)
        residual.set(
            xlabel="Pressure coordinate used (GPa)",
            ylabel="P(published Fe₃C curve) − P (GPa)",
            title="Residuals against the unchanged published curve",
        )
        for ax in (pv, delta, residual):
            ax.grid(color="#e7eaee", alpha=0.65)
            ax.set_axisbelow(True)
            ax.spines[["top", "right"]].set_visible(False)
        fig.text(
            0.5,
            0.025,
            "Recalculated P: measured MgO V + printed Speziale BM3 (V₀=74.71 Å³, K₀=160.2 GPa, K₀′=3.99).\nDerived coordinates, not established corrections. No new pressure-error bars or Fe₃C fit; Li points excluded from this comparison.",
            ha="center",
            fontsize=8,
        )
        fig.get_layout_engine().set(rect=(0, 0.09, 1, 0.91))
        FIGURE.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            FIGURE,
            metadata={
                "Date": None,
                "Title": "Sata 2010 Fe3C pressure-recalculation comparison",
                "Description": "Original and derived pressure coordinates from the same nine measured sample/calibrant pairs. Printed Speziale BM3 coefficients; no Fe3C refit, source replacement or uncertainty inference.",
            },
        )
        FIGURE.with_name(Path(dark_path(FIGURE.name)).name).write_bytes(
            dark_svg(FIGURE.read_bytes())
        )
        if preview:
            preview.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(preview, dpi=160)
        plt.close(fig)
    return FIGURE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path, help="optional PNG preview path")
    args = parser.parse_args()
    print(plot(args.preview).relative_to(ROOT))


if __name__ == "__main__":
    main()
