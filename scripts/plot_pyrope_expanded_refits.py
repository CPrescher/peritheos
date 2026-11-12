"""Figures for the expanded conditional pyrope audit, not recommended EOSs."""

from __future__ import annotations

import json

import matplotlib
import numpy as np

from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.report_plot_themes import dark_svg


def plot():
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = json.loads((DATA / "expanded-refits.json").read_text())
    old = json.loads((DATA / "subset-refits.json").read_text())
    fits = [old["cases"]["joint_recovered"]["refit"]] + [
        report["cases"][n]["refit"] for n in ("182", "191")
    ]
    with plt.rc_context({"svg.hashsalt": "pyrope-expanded", "font.size": 10}):
        fig, axes = plt.subplots(2, 2, figsize=(9, 6.5))
        for ax, key, label in zip(
            axes.flat,
            ["K0_gpa", "K0_prime", "theta0_k", "gamma0"],
            ["K0 (GPa)", "K0′", "θ0 (K)", "γ0"],
        ):
            ax.axhline(
                report["published_parameters"][key],
                color="#555555",
                ls="--",
                label="Published coefficient",
            )
            ax.errorbar(
                np.arange(3),
                [f["parameters"][key] for f in fits],
                yerr=[f["conditional_residual_scaled_widths"][key] for f in fits],
                fmt="o",
                color="#2463a2",
                capsize=4,
            )
            ax.set_xticks(np.arange(3), ["92 channels", "182 channels", "191 channels"])
            ax.set_ylabel(label)
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
        axes[0, 0].legend(loc="lower right", fontsize=8)
        fig.suptitle("Pyrope: expanded conditional source-model fits", fontsize=12)
        fig.text(
            0.5,
            0.015,
            "Bars: conditional residual-scaled fit widths; unknown weights and systematic covariance excluded.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.04, 1, 0.96))
        save(fig, "expanded-coefficients")

        fit = fits[-1]
        groups = [
            (["milani_vt", "hartwig_vt"], "T", "Expansion volumes"),
            (
                [
                    "zhang_helium",
                    "zhang_neon",
                    "du_pv",
                    "milani_pv",
                    "zou_pv",
                    "chantel_pv",
                ],
                "P",
                "Compression volumes",
            ),
            (
                ["zou_ks", "chantel_ks", "sinogeikin_2000_ks"],
                "P",
                "High-pressure moduli",
            ),
            (["sinogeikin_2002_ks"], "T", "Nine digitized ambient-pressure moduli"),
        ]
        fig, axes = plt.subplots(2, 2, figsize=(10, 7.6))
        for ax, (names, key, title) in zip(axes.flat, groups):
            for name in names:
                rows = [r for r in fit["rows"] if r["dataset"] == name]
                ax.errorbar(
                    [r[key] for r in rows],
                    [r["calculated_minus_observed_pressure_gpa"] for r in rows],
                    yerr=[r["effective_pressure_error_gpa"] for r in rows],
                    fmt="o",
                    ms=3,
                    capsize=2,
                    alpha=0.8,
                    label=name.replace("_", " "),
                )
            ax.axhline(0, color="#777777", lw=0.7)
            ax.set_xlabel("Temperature (K)" if key == "T" else "Pressure (GPa)")
            ax.set_ylabel("Calculated − observed P (GPa)")
            ax.set_title(title)
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
            ax.legend(
                fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2
            )
        fig.suptitle(
            "191 of 194 listed channels: conditional fit residuals", fontsize=12
        )
        fig.tight_layout(rect=(0, 0, 1, 0.96), h_pad=4)
        save(fig, "expanded-residuals")


def save(fig, name):
    import matplotlib.pyplot as plt

    path = DATA / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None})
    (DATA / f"{name}-dark.svg").write_bytes(dark_svg(path.read_bytes()))
    output = DATA.parents[3] / "output/pyrope-full-fit-2026-11-12"
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{name}.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    plot()
