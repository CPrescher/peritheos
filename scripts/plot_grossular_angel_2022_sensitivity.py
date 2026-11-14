"""Plot declared weight sensitivity, not statistical confidence regions.

Repository-only: python -m scripts.plot_grossular_angel_2022_sensitivity.
"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from scripts.audit_grossular_angel_2022_observations import DATA


def plot(dark=False, preview=None):
    report = json.loads((DATA / "full-sensitivity-refits.json").read_text())
    cases = [r for r in report["cases"] if r["isaak_profile"] == "interpolated"]
    background, foreground = ("#20242c", "#f1f3f5") if dark else ("white", "#20242c")
    with plt.rc_context(
        {
            "font.size": 11,
            "figure.facecolor": background,
            "axes.facecolor": background,
            "text.color": foreground,
            "axes.labelcolor": foreground,
            "xtick.color": foreground,
            "ytick.color": foreground,
            "axes.edgecolor": foreground,
            "svg.fonttype": "none",
        }
    ):
        fig, axes = plt.subplots(2, 3, figsize=(11, 7))
        fig.subplots_adjust(
            left=0.08, right=0.91, bottom=0.10, top=0.94, wspace=0.24, hspace=0.30
        )
        for row, (parameter, label, cmap) in enumerate(
            [
                ("theta0_k", "θ₀ (K)", "Blues"),
                ("gamma0", "γ₀", "Oranges"),
            ]
        ):
            minimum = min(r["fits"][0]["parameters"][parameter] for r in cases)
            maximum = max(r["fits"][0]["parameters"][parameter] for r in cases)
            for col, isaak in enumerate([0.5, 1, 2]):
                values = np.empty((3, 3))
                for y, temperature in enumerate([0.5, 1, 2]):
                    for x, volume in enumerate([0.5, 1, 2]):
                        case = next(
                            r
                            for r in cases
                            if r["isaak_multiplier"] == isaak
                            and r["thieblot_volume_multiplier"] == volume
                            and r["thieblot_temperature_multiplier"] == temperature
                        )
                        values[y, x] = case["fits"][0]["parameters"][parameter]
                ax = axes[row, col]
                mesh = ax.pcolormesh(
                    np.arange(4),
                    np.arange(4),
                    values,
                    cmap=cmap,
                    vmin=minimum,
                    vmax=maximum,
                    edgecolors=background,
                    linewidth=2,
                )
                ax.set(
                    xticks=[0.5, 1.5, 2.5],
                    yticks=[0.5, 1.5, 2.5],
                    xticklabels=["0.5×", "1×", "2×"],
                    yticklabels=["0.5×", "1×", "2×"],
                    title=f"{label} · Isaak {isaak:g}×",
                    aspect="equal",
                )
                if row == 1:
                    ax.set_xlabel("Thiéblot σV multiplier")
                if col == 0:
                    ax.set_ylabel("Thiéblot σT multiplier")
                else:
                    ax.tick_params(labelleft=False)
                if isaak == 1:
                    ax.add_patch(
                        Rectangle(
                            (1, 1), 1, 1, fill=False, edgecolor="#d40050", linewidth=2.5
                        )
                    )
            colorbar = fig.colorbar(mesh, ax=axes[row, :], shrink=0.85, pad=0.02)
            colorbar.set_label(label)
        name = (
            "full-error-sensitivity-dark.svg" if dark else "full-error-sensitivity.svg"
        )
        if preview:
            fig.savefig(preview, dpi=160)
        fig.savefig(DATA / name)
        plt.close(fig)


if __name__ == "__main__":
    plot()
    plot(dark=True)
