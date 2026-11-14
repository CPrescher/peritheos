"""Compare declared covariance and preprocessing controls, not confidence bands."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.audit_grossular_angel_2022_observations import DATA


def plot(dark=False, preview=None):
    r = json.loads((DATA / "du-bosenick-refits.json").read_text())
    values = {c["id"]: c["fits"][0]["parameters"]["theta0_k"] for c in r["cases"]}
    bg, fg = ("#20242c", "#f1f3f5") if dark else ("white", "#20242c")
    with plt.rc_context(
        {
            "font.size": 10,
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
        fig, axes = plt.subplots(2, 2, sharey=True, figsize=(12.3, 8.0))
        fig.subplots_adjust(
            left=0.07, right=0.98, bottom=0.11, top=0.84, wspace=0.12, hspace=0.40
        )
        covariances = [
            "diagonal",
            "volume_anchor",
            "full_independent_T",
            "full_shared_T",
        ]
        treatments = [
            "table_nacl",
            "lattice",
            "tc_reconstructed",
            "double_corrected",
            "shared_T",
        ]
        for col, background in enumerate(["plotted", "linear_298"]):
            for anchor, color, marker, label in [
                ("measured_295", "#0072B2", "o", "Bosenick measured 295 K anchor"),
                ("linear_298", "#E69F00", "s", "Bosenick calculated 298 K anchor"),
            ]:
                top = [
                    values[f"{background}__{anchor}__{c}__table_nacl"]
                    for c in covariances
                ]
                bottom = [
                    values[f"{background}__{anchor}__full_independent_T__{t}"]
                    for t in treatments
                ]
                for ax, y in [(axes[0, col], top), (axes[1, col], bottom)]:
                    ax.plot(
                        range(len(y)),
                        y,
                        marker=marker,
                        color=color,
                        linestyle="--" if marker == "s" else "-",
                        markersize=6,
                        linewidth=1.0,
                        label=label,
                    )
            for ax in axes[:, col]:
                ax.axhline(
                    749.825,
                    color=fg,
                    linestyle=":",
                    linewidth=1.0,
                    label="Published θ = 749.825 K",
                )
                ax.set_ylim(740, 905)
                ax.grid(axis="y", color=fg, alpha=0.12)
            axes[0, col].set_title(
                "Milani/Du plotted background"
                if col == 0
                else "Milani/Du linear 298 K background"
            )
            axes[0, col].set(
                xticks=range(4),
                xticklabels=[
                    "Diagonal",
                    "Volume\nanchor",
                    "Full, independent\nraw T errors",
                    "Full, shared\nraw T errors",
                ],
            )
            axes[1, col].set(
                xticks=range(5),
                xticklabels=[
                    "Listed NaCl",
                    "a³ and σₐ",
                    "TC stage*",
                    "Double corr.*",
                    "Shared T",
                ],
            )
            axes[1, col].set_xlabel(
                "Du treatment, Bosenick full independent-T covariance"
            )
        axes[0, 0].set_ylabel("θ (K)")
        axes[1, 0].set_ylabel("θ (K)")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.53, 0.96),
            frameon=False,
            fontsize=10,
        )
        fig.suptitle(
            "Du processing and Bosenick anchor/covariance controls",
            y=0.995,
            fontsize=12,
        )
        fig.text(
            0.52,
            0.025,
            "All 219 rows retained; independent 2% Conrad diagnostic fixed.  *Unsupported calibration-stage error controls.",
            ha="center",
            fontsize=9,
        )
        fig.savefig(DATA / ("du-bosenick-dark.svg" if dark else "du-bosenick.svg"))
        if preview:
            fig.savefig(preview, dpi=145)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
