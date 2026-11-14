"""Display theta influence together with changes to Kprime."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA

LABELS = {
    "baseline": "All 219 rows",
    "du_vt": "Du V–T (13)",
    "bosenick_vt": "Bosenick V–T (19)",
    "milani_vt": "Milani V–T (35)",
    "rodehorst_vt": "Rodehorst V–T (9)",
    "isaak_vt": "Isaak GR-2B V–T (15)",
    "thieblot_vt_estimated": "Thiéblot V–T (7)",
    "milani_pv": "Milani P–V (7)",
    "kono_ks": "Kono Kₛ (36)",
    "gwanmesia_ks": "Gwanmesia Kₛ (52)",
    "isaak_ks": "Isaak GGD Kₛ (22)",
    "bass_ks": "Bass Kₛ (1)",
    "milani_all": "Milani, both series (42)",
    "isaak_all": "Isaak, both specimens (37)",
}


def plot(dark=False, preview=None):
    report = json.loads((DATA / "theta-influence-refits.json").read_text())
    parameters = {c["id"]: c["fits"][0]["parameters"] for c in report["cases"]}
    bg, fg = ("#20242c", "#f1f3f5") if dark else ("white", "#20242c")
    with plt.rc_context(
        {
            "font.size": 10.5,
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
        fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12.2, 7.6))
        fig.subplots_adjust(left=0.28, right=0.98, bottom=0.12, top=0.84, wspace=0.18)
        for ax, field, published, label in zip(
            axes, ["theta0_k", "K0_prime"], [749.825, 5.06534], ["θ (K)", "K′"]
        ):
            for norm, offset, marker, color, legend in [
                ("plotted", -0.13, "o", "#0072B2", "Plotted normalization"),
                ("linear_298", 0.13, "s", "#E69F00", "Linear 298 K normalization"),
            ]:
                x = [
                    parameters[f"{norm}__" + (g if g == "baseline" else "omit_" + g)][
                        field
                    ]
                    for g in LABELS
                ]
                ax.scatter(
                    x,
                    np.arange(len(LABELS)) + offset,
                    s=42,
                    marker=marker,
                    color=color,
                    label=legend,
                    zorder=3,
                )
            ax.axvline(
                published,
                color=fg,
                linestyle="--",
                linewidth=1,
                label="Published reference",
            )
            ax.axhline(0.6, color=fg, linewidth=0.6, alpha=0.35)
            ax.set(
                xlabel=label,
                yticks=np.arange(len(LABELS)),
                yticklabels=list(LABELS.values()),
                ylim=(len(LABELS) - 0.4, -0.7),
            )
            ax.grid(axis="x", color=fg, alpha=0.12)
        axes[0].set_xlim(700, 905)
        axes[1].set_xlim(4.84, 5.115)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.64, 0.955),
            frameon=False,
            fontsize=10,
        )
        fig.suptitle(
            "Dataset omissions with independent 2% Conrad diagnostic errors",
            y=0.995,
            fontsize=12,
        )
        fig.text(
            0.64,
            0.025,
            "Top row: baseline; other rows: omit named series. Conrad retained throughout.",
            ha="center",
            fontsize=9,
        )
        fig.savefig(
            DATA / ("theta-influence-dark.svg" if dark else "theta-influence.svg")
        )
        if preview:
            fig.savefig(preview, dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
