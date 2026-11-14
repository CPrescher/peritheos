"""Plot the explicitly diagnostic Hartwig weight sensitivity."""

from __future__ import annotations

import json

import matplotlib

from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.report_plot_themes import dark_svg


def plot():
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = json.loads((DATA / "hartwig-weight-sensitivity.json").read_text())
    with plt.rc_context({"svg.hashsalt": "pyrope-hartwig-weight", "font.size": 10}):
        fig, axes = plt.subplots(2, 2, figsize=(10, 7.6))
        for ax, key, label in zip(
            axes[0],
            ["theta0_k", "gamma0"],
            ["Debye temperature θ0 (K)", "Grüneisen parameter γ0"],
        ):
            center = report["source_parameters"][key]
            width = report["source_covariance_marginal_widths"][key]
            ax.axhspan(center - width, center + width, color="#777777", alpha=0.15)
            ax.axhline(
                center, ls="--", color="#555555", lw=1, label="Published coefficient"
            )
            for count, color in [("182", "#ae4c14"), ("191", "#245b8a")]:
                rows = report["cases"][count]
                ax.plot(
                    [r["weight_multiplier"] for r in rows],
                    [r["fit_from_published_start"]["parameters"][key] for r in rows],
                    "o-",
                    ms=4,
                    lw=1,
                    label=f"{count} channels",
                    color=color,
                )
            ax.set_xscale("log")
            ax.set_xlabel("Hartwig inverse-variance multiplier w")
            ax.set_ylabel(label)
            ax.legend(fontsize=8, loc="lower right")
        ax = axes[1, 0]
        rows = report["cases"]["191"]
        theta = [r["fit_from_published_start"]["parameters"]["theta0_k"] for r in rows]
        gamma = [r["fit_from_published_start"]["parameters"]["gamma0"] for r in rows]
        ax.plot(
            theta, gamma, "o-", ms=4, color="#245b8a", label="191-channel weight path"
        )
        ax.plot(
            report["source_parameters"]["theta0_k"],
            report["source_parameters"]["gamma0"],
            "D",
            color="#555555",
            label="Published coefficients",
        )
        for r, x, y in zip(rows, theta, gamma):
            if r["weight_multiplier"] in (1, 100, 500, 1000, 10000):
                offset = (
                    (5, 8)
                    if r["weight_multiplier"] == 1
                    else (-35, -14)
                    if r["weight_multiplier"] == 10000
                    else (4, -13)
                )
                ax.annotate(
                    str(r["weight_multiplier"]),
                    (x, y),
                    xytext=offset,
                    textcoords="offset points",
                    fontsize=8,
                )
        ax.set_xlabel("θ0 (K)")
        ax.set_ylabel("γ0")
        ax.legend(fontsize=8, loc="upper left")
        ax.set_title("Both thermal coefficients must be considered")
        ax = axes[1, 1]
        for name, label, color in [
            ("hartwig_vt", "Hartwig expansion", "#ae4c14"),
            ("milani_vt", "Milani expansion", "#245b8a"),
        ]:
            ax.plot(
                [r["weight_multiplier"] for r in rows],
                [
                    r["baseline_error_residuals"][name]["pressure_residual_rms_gpa"]
                    for r in rows
                ],
                "o-",
                color=color,
                ms=4,
                lw=1,
                label=label,
            )
        ax.set_xscale("log")
        ax.set_xlabel("Hartwig inverse-variance multiplier w")
        ax.set_ylabel("Pressure residual RMS (GPa)")
        ax.set_title("Tradeoff between expansion datasets")
        ax.legend(fontsize=8, loc="upper left")
        for ax in axes.flat:
            ax.tick_params(direction="in", top=True, right=True)
            ax.grid(alpha=0.15)
        fig.suptitle("Pyrope: sensitivity to increased Hartwig weights", fontsize=12)
        fig.text(
            0.5,
            0.015,
            "Diagnostic only. Shading: source marginal covariance widths; confidence unspecified. Author multiplier unknown.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.04, 1, 0.96))
        path = DATA / "hartwig-weight-sensitivity.svg"
        fig.savefig(path, metadata={"Date": None})
        (DATA / "hartwig-weight-sensitivity-dark.svg").write_bytes(
            dark_svg(path.read_bytes())
        )
        output = DATA.parents[3] / "output/pyrope-hartwig-weights-2026-11-12"
        output.mkdir(parents=True, exist_ok=True)
        fig.savefig(output / "sensitivity.png", dpi=140)
        plt.close(fig)


if __name__ == "__main__":
    plot()
