"""Plot conditional pyrope coefficient fits and recovered-subset residuals."""

from __future__ import annotations

import json

import matplotlib
import numpy as np

from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.report_plot_themes import dark_svg


def plot():
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = json.loads((DATA / "subset-refits.json").read_text())
    cases = report["cases"]
    with plt.rc_context({"svg.hashsalt": "pyrope-subset-refits", "font.size": 10}):
        fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1))
        labels = ["Milani VT", "Hartwig VT", "Both VT", "Joint subset"]
        for ax, key, ylabel in zip(
            axes,
            ["theta0_k", "gamma0"],
            ["Debye temperature θ0 (K)", "Grüneisen parameter γ0"],
        ):
            source = report["source_parameters"][key]
            width = report["source_covariance_marginal_widths"][key]
            ax.axhspan(source - width, source + width, color="#777777", alpha=0.16)
            ax.axhline(
                source, color="#555555", ls="--", lw=1, label="Published coefficient"
            )
            ax.errorbar(
                np.arange(4),
                [case["refit"]["parameters"][key] for case in cases.values()],
                yerr=[
                    case["refit"]["conditional_residual_scaled_widths"][key]
                    for case in cases.values()
                ],
                fmt="o",
                capsize=4,
                color="#2463a2",
            )
            ax.set_xticks(np.arange(4), labels, rotation=15)
            ax.set_ylabel(ylabel)
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
        axes[0].legend(loc="lower right", fontsize=9)
        fig.suptitle(
            "Conditional source-model fits; original joint input not recovered",
            fontsize=11,
        )
        fig.text(
            0.5,
            0.015,
            "Bars: conditional fit widths. Shading: source covariance marginal width; confidence unspecified.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 0.94))
        save(fig, "subset-coefficients")

        fig, axes = plt.subplots(2, 2, figsize=(10, 6.8))
        joint = cases["joint_recovered"]
        for ax, dataset, label in zip(
            axes.flat,
            ["milani_vt", "hartwig_vt", "zou_pv", "zou_ks"],
            [
                "Milani VT: 33 volumes",
                "Hartwig VT: 12 volumes",
                "Zou PV: seven 300 K volumes",
                "Zou KS: 40 moduli",
            ],
        ):
            key = "P" if dataset == "zou_pv" else "T"
            before = [
                r
                for r in joint["source_eos_scales_only"]["rows"]
                if r["dataset"] == dataset
            ]
            after = [r for r in joint["refit"]["rows"] if r["dataset"] == dataset]
            ax.scatter(
                [r[key] for r in before],
                [r["calculated_minus_observed_pressure_gpa"] for r in before],
                marker="x",
                color="#999999",
                label="Published EOS; scales profiled",
                s=22,
            )
            ax.errorbar(
                [r[key] for r in after],
                [r["calculated_minus_observed_pressure_gpa"] for r in after],
                yerr=[r["effective_pressure_error_gpa"] for r in after],
                fmt="o",
                markersize=3,
                capsize=2,
                color="#2463a2",
                label="Conditional joint-subset refit",
                alpha=0.8,
            )
            ax.axhline(0, lw=0.7, color="#777777")
            ax.set_title(label)
            ax.set_xlabel("Pressure (GPa)" if key == "P" else "Temperature (K)")
            ax.set_ylabel("Calculated − observed P (GPa)")
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
        axes[0, 0].legend(fontsize=8, loc="upper left")
        fig.suptitle(
            "92 recovered data channels: partial source-model refit", fontsize=12
        )
        fig.text(
            0.5,
            0.015,
            "Bars: effective pressure errors from the stated initial coordinate-error protocol. Missing terms omitted.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.04, 1, 0.95))
        save(fig, "subset-residuals")


def save(fig, name):
    import matplotlib.pyplot as plt

    path = DATA / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None})
    (DATA / f"{name}-dark.svg").write_bytes(dark_svg(path.read_bytes()))
    plt.close(fig)


if __name__ == "__main__":
    plot()
