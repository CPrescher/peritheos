"""Plot the predeclared series deletions, not confidence intervals."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA

LABELS = {
    "baseline": "All 219 rows",
    "conrad_ks": "Conrad Kₛ (3)",
    "du_vt": "Du V–T (13)",
    "kono_ks": "Kono Kₛ (36)",
    "gwanmesia_ks": "Gwanmesia Kₛ (52)",
    "milani_pv": "Milani P–V (7)",
    "isaak_ks": "Isaak GGD Kₛ (22)",
    "bass_ks": "Bass Kₛ (1)",
    "milani_vt": "Milani V–T (35)",
    "bosenick_vt": "Bosenick V–T (19)",
    "rodehorst_vt": "Rodehorst V–T (9)",
    "isaak_vt": "Isaak GR-2B V–T (15)",
    "thieblot_vt_estimated": "Thiéblot V–T (7)",
    "milani_all": "Milani, both series (42)",
    "isaak_all": "Isaak, both specimens (37)",
}


def plot(dark=False, preview=None):
    report = json.loads((DATA / "dataset-influence-refits.json").read_text())
    by_id = {r["id"]: r["fits"][0]["parameters"]["K0_prime"] for r in report["cases"]}
    order = list(LABELS)
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
        fig, ax = plt.subplots(figsize=(10, 7.8))
        fig.subplots_adjust(left=0.32, right=0.96, bottom=0.12, top=0.85)
        for norm, offset, marker, color, label in [
            ("plotted", -0.13, "o", "#0072B2", "Plotted normalization"),
            ("linear_298", 0.13, "s", "#E69F00", "Linear 298 K normalization"),
        ]:
            x = [
                by_id[norm + "__" + (g if g == "baseline" else "omit_" + g)]
                for g in order
            ]
            ax.scatter(
                x,
                np.arange(len(order)) + offset,
                marker=marker,
                color=color,
                s=48,
                label=label,
                zorder=3,
            )
        ax.axvline(
            5.06534,
            color=fg,
            linestyle="--",
            linewidth=1.1,
            label="Published K′ = 5.06534",
        )
        ax.axhline(0.6, color=fg, linewidth=0.6, alpha=0.35)
        ax.set(
            yticks=np.arange(len(order)),
            yticklabels=[LABELS[g] for g in order],
            xlabel="K′ after omitting the named series",
            xlim=(5.035, 5.43),
            ylim=(len(order) - 0.4, -0.8),
            title="Grossular: dataset influence on K′",
        )
        ax.grid(axis="x", color=fg, alpha=0.12)
        fig.legend(
            loc="upper center", bbox_to_anchor=(0.64, 0.98), frameon=False, fontsize=10
        )
        name = "dataset-influence-dark.svg" if dark else "dataset-influence.svg"
        fig.savefig(DATA / name)
        if preview:
            fig.savefig(preview, dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
