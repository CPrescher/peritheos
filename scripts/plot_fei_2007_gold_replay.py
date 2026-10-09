"""Plot Au measured residuals and objective sensitivity; no EOS refitting."""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.reproduce_fei_2007_gold import (
    OUTPUT,
    PUBLISHED,
    ROOT,
    effective_sigma,
    observations,
    pressure,
)

PLOT = ROOT / "docs/data/fei-2007-gold-replay.png"


def main():
    report = json.loads(OUTPUT.read_text(encoding="utf-8"))
    _, d = observations()
    fits = report["fits"]
    colors = ["#374151", "#0072B2", "#D55E00"]
    fig, axes = plt.subplots(
        1, 3, figsize=(13, 4.9), gridspec_kw={"width_ratios": [1.15, 1, 1.1]}
    )
    for ax, mask, title in zip(
        axes[:2],
        [d["p"] < 50, d["p"] > 50],
        ["Fei (2004): 26 Au-MgO states", "Hirose (2006): two Au-MgO states"],
    ):
        published = pressure(d["v"], d["t"]) - d["p"]
        ax.errorbar(
            d["p"][mask],
            published[mask],
            yerr=effective_sigma(d)[mask],
            fmt="o",
            ms=4,
            color=colors[0],
            ecolor="#9CA3AF",
            capsize=2,
            label="Published q = 0.6",
        )
        for label, color, marker in [
            ("all28_equal", colors[1], "s"),
            ("all28_available_errors", colors[2], "^"),
        ]:
            ax.scatter(
                d["p"][mask],
                np.array(fits[label]["residuals_gpa"])[mask],
                color=color,
                marker=marker,
                s=22,
                label=f"{('Equal weights' if label.endswith('equal') else 'Available errors')}: q = {fits[label]['parameters']['q']:.3f}",
            )
        ax.axhline(0, lw=1, color="#6B7280")
        ax.set_title(title, fontsize=11, pad=13)
        ax.set_xlabel("Reported MgO pressure (GPa)")
        ax.set_ylabel("Au model - reported MgO pressure (GPa)")
        ax.grid(alpha=0.18)
    for i in (26, 27):
        axes[1].annotate(
            f"{d['t'][i]:.0f} K",
            (d["p"][i], pressure(d["v"][i], d["t"][i]) - d["p"][i]),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=9,
        )
    axes[1].set_xlim(102, 122)
    axes[0].legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.20), fontsize=9, frameon=False
    )
    q_values = np.linspace(0.05, 1.75, 180)
    for mode, color, name in [
        ("equal", colors[1], "Equal pressure weights"),
        ("available_errors", colors[2], "Available measurement errors"),
    ]:
        sigma = np.ones(28) if mode == "equal" else effective_sigma(d)
        costs = []
        for q in q_values:
            c = PUBLISHED.copy()
            c[4] = q
            costs.append(np.sum(((pressure(d["v"], d["t"], c) - d["p"]) / sigma) ** 2))
        fit = fits["all28_" + mode]
        axes[2].plot(
            q_values,
            np.array(costs) / fit["objective_sum"] - 1,
            color=color,
            label=name,
        )
        axes[2].axvline(fit["parameters"]["q"], color=color, lw=1, ls=":")
    axes[2].set(
        xlabel="q (other coefficients fixed)",
        ylabel="Objective / minimum - 1",
        ylim=(-0.01, 1.5),
    )
    axes[2].set_title("Objective depends on weighting", fontsize=11, pad=13)
    axes[2].grid(alpha=0.18)
    axes[2].legend(fontsize=9, frameon=False)
    fig.suptitle(
        "Fei (2007) Au thermal EOS: conditional measured-data replay",
        fontsize=15,
        y=0.985,
    )
    fig.text(
        0.035,
        0.035,
        "Bars show available pressure + Au-volume widths only; Hirose volume errors and temperature/calibration covariance are missing.\nAll diagnostics use reported MgO pressures. Objective curves are not confidence intervals; original author weights are unknown.",
        fontsize=9,
        color="#374151",
    )
    fig.subplots_adjust(top=0.82, bottom=0.32, wspace=0.38, left=0.06, right=0.98)
    fig.savefig(PLOT, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
