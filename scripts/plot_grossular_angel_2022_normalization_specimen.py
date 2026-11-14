"""Display paired normalization/specimen effects, not confidence regions."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scripts.audit_grossular_angel_2022_observations import DATA


def plot(dark=False, preview=None):
    report = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    colors = {
        "plotted": "#0072B2",
        "measured_room": "#E69F00",
        "linear_298": "#009E73",
        "free_milani_du": "#CC79A7",
        "free_diffraction": "#D55E00",
    }
    labels = {
        "plotted": "Fixed plotted normalization",
        "measured_room": "Tabulated room anchors",
        "linear_298": "Linear 298 K anchors",
        "free_milani_du": "Free Milani/Du offsets",
        "free_diffraction": "Free five diffraction offsets",
    }
    markers = {
        "as_plotted": "o",
        "reference_298": "D",
        "omit_expansion": "s",
        "calculated_endmember_alpha": "^",
    }
    specimen_labels = {
        "as_plotted": "GR-2B as plotted",
        "reference_298": "GR-2B referenced to 298 K",
        "omit_expansion": "Omit GR-2B expansion",
        "calculated_endmember_alpha": "Calculated endmember α",
    }
    bg, fg = ("#20242c", "#f1f3f5") if dark else ("white", "#20242c")
    with plt.rc_context(
        {
            "font.size": 11,
            "figure.facecolor": bg,
            "axes.facecolor": bg,
            "text.color": fg,
            "axes.labelcolor": fg,
            "axes.edgecolor": fg,
            "xtick.color": fg,
            "ytick.color": fg,
            "svg.fonttype": "none",
        }
    ):
        fig, ax = plt.subplots(figsize=(10.5, 5.6))
        fig.subplots_adjust(left=0.09, right=0.65, bottom=0.14, top=0.90)
        for case in report["cases"]:
            p = case["fits"][0]["parameters"]
            ax.scatter(
                p["K0_prime"],
                p["theta0_k"],
                marker=markers[case["isaak"]],
                s=75,
                color=colors[case["normalization"]],
                edgecolors=fg,
                linewidths=0.5,
                zorder=3,
            )
        source = report["cases"][0]["effects"]
        kp, theta = source["K0_prime"]["published"], source["theta0_k"]["published"]
        ax.axvline(kp, color=fg, linewidth=0.7, linestyle="--", alpha=0.4)
        ax.axhline(theta, color=fg, linewidth=0.7, linestyle="--", alpha=0.4)
        ax.scatter(kp, theta, marker="*", s=180, color=fg, zorder=4)
        ax.annotate(
            "Published",
            (kp, theta),
            xytext=(7, -16),
            textcoords="offset points",
            fontsize=10,
        )
        ax.set(
            xlabel="K′₀T",
            ylabel="θ₀ (K)",
            xlim=(5.045, 5.255),
            ylim=(730, 790),
            title="Normalization and Isaak expansion tests",
        )
        normalization = fig.legend(
            [
                Line2D([], [], marker="o", linestyle="none", color=c, markersize=7)
                for c in colors.values()
            ],
            list(labels.values()),
            title="Normalization",
            loc="upper left",
            bbox_to_anchor=(0.68, 0.91),
            frameon=False,
            fontsize=9.5,
        )
        fig.add_artist(normalization)
        fig.legend(
            [
                Line2D([], [], marker=m, linestyle="none", color=fg, markersize=7)
                for m in markers.values()
            ],
            list(specimen_labels.values()),
            title="Isaak expansion",
            loc="upper left",
            bbox_to_anchor=(0.68, 0.51),
            frameon=False,
            fontsize=9.5,
        )
        if preview:
            fig.savefig(preview, dpi=160)
        name = (
            "normalization-specimen-effects-dark.svg"
            if dark
            else "normalization-specimen-effects.svg"
        )
        fig.savefig(DATA / name)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
