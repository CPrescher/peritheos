"""Render the retained grossular comparisons, without re-fitting anything."""

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.audit_grossular_angel_2022_observations import DATA


def plot(preview=None):
    audit = json.loads((DATA / "observation-audit.json").read_text())
    colors = ["#337bb6", "#df8738", "#49a488", "#aa75bd"]
    for dark in (False, True):
        with plt.style.context("dark_background" if dark else "default"):
            plt.rcParams.update(
                {
                    "font.size": 10,
                    "svg.fonttype": "none",
                    "svg.hashsalt": "grossular-audit-20261112",
                }
            )
            fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), layout="constrained")
            for ax in axes:
                ax.axhline(0, color=".6", linewidth=0.8)
                ax.tick_params(direction="in", top=True, right=True)
                ax.grid(alpha=0.13)
            for (name, group), color, label in zip(
                audit["volume_shape_comparisons"].items(),
                colors,
                [
                    "Milani 2017",
                    "Du 2015 preprint",
                    "Bosenick–Geiger 1997*",
                    "Rodehorst 2002*",
                ],
            ):
                rows = group["rows"]
                axes[0].scatter(
                    [r["T"] for r in rows],
                    [r["observed_minus_calculated_cm3_mol"] for r in rows],
                    s=17,
                    color=color,
                    label=label,
                )
            axes[0].set(
                xlabel="Temperature (K)",
                ylabel="Observed − calculated V (cm³/mol)",
                title="Expansion shape; one volume scale per set",
            )
            axes[0].legend(
                loc="upper left", bbox_to_anchor=(0, -0.20), frameon=False, fontsize=8
            )
            for (name, group), color in zip(audit["moduli"].items(), colors):
                rows = [r for r in group["rows"] if r["selected"]]
                axes[1].scatter(
                    [r["temperature_k"] for r in rows],
                    [r["observed_minus_calculated_gpa"] for r in rows],
                    s=20,
                    color=color,
                    label=name.replace("_ggd", " GGD").title(),
                )
            axes[1].set(
                xlabel="Temperature (K)",
                ylabel="Observed − calculated Kₛ (GPa)",
                title="Moduli; source Kono scale = 0.99455",
            )
            axes[1].legend(
                loc="upper left", bbox_to_anchor=(0, -0.20), frameon=False, fontsize=8
            )
            rows = audit["milani_compression"]["rows"]
            axes[2].scatter(
                [r["recalibrated_pressure_gpa"] for r in rows],
                [r["observed_minus_calculated_cm3_mol"] for r in rows],
                s=25,
                color=colors[0],
            )
            axes[2].set(
                xlabel="Recalibrated pressure (GPa)",
                ylabel="Observed − calculated V (cm³/mol)",
                title="Milani compression; ambient anchor fixed",
                ylim=(-0.04, 0.04),
            )
            fig.suptitle("Grossular: unchanged Angel (2022) coefficients", fontsize=13)
            suffix = "-dark" if dark else ""
            fig.savefig(DATA / f"source-residuals{suffix}.svg", metadata={"Date": None})
            if preview is not None and not dark:
                fig.savefig(preview, dpi=160)
            plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview")
    args = parser.parse_args()
    plot(args.preview)
