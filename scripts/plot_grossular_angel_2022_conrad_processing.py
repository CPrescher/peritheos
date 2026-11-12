"""Show the declared Conrad treatments with every selected point retained."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA

LABELS = {
    "baseline": "Baseline formal errors",
    "cij_negative_correlation": "Formal Cᵢⱼ correlation −1",
    "cij_independent": "Formal Cᵢⱼ correlation 0",
    "extra_005": "Extra 0.5% Kₛ, independent",
    "extra_010": "Extra 1% Kₛ, independent",
    "extra_020": "Extra 2% Kₛ, independent",
    "shared_020": "Extra 2% Kₛ, fully shared",
    "pressure_02": "Assigned σ(P) = 0.2 GPa",
    "pressure_05": "Assigned σ(P) = 0.5 GPa",
    "round_01": "Round Kₛ to 0.1 GPa",
    "density_vp": "Density control from vₚ",
    "density_vs": "Density control from vₛ",
    "free_conrad_scale": "Free Conrad multiplier",
}


def plot(dark=False, preview=None):
    r = json.loads((DATA / "conrad-processing-refits.json").read_text())
    values = {c["id"]: c["fits"][0]["parameters"]["K0_prime"] for c in r["cases"]}
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
        fig, ax = plt.subplots(figsize=(10, 7.3))
        fig.subplots_adjust(left=0.35, right=0.96, bottom=0.12, top=0.84)
        for n, offset, marker, color, label in [
            ("plotted", -0.13, "o", "#0072B2", "Plotted normalization"),
            ("linear_298", 0.13, "s", "#E69F00", "Linear 298 K normalization"),
        ]:
            ax.scatter(
                [values[f"{n}__{v}"] for v in LABELS],
                np.arange(len(LABELS)) + offset,
                s=48,
                marker=marker,
                color=color,
                label=label,
                zorder=3,
            )
        ax.axvline(
            5.06534,
            linestyle="--",
            linewidth=1,
            color=fg,
            label="Published K′ = 5.06534",
        )
        ax.set(
            yticks=np.arange(len(LABELS)),
            yticklabels=list(LABELS.values()),
            xlabel="Fitted K′ (all 219 rows retained)",
            xlim=(5.045, 5.26),
            ylim=(len(LABELS) - 0.4, -0.6),
            title="Conrad preprocessing and assigned errors",
        )
        ax.grid(axis="x", color=fg, alpha=0.12)
        fig.legend(
            loc="upper center", bbox_to_anchor=(0.65, 0.98), frameon=False, fontsize=10
        )
        fig.savefig(
            DATA / ("conrad-processing-dark.svg" if dark else "conrad-processing.svg")
        )
        if preview:
            fig.savefig(preview, dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
