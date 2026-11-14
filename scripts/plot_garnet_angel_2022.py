"""Plot independent source-law curves and official calculated checkpoints.

Repository-only: ``python -m scripts.plot_garnet_angel_2022``. No observations
or optimizations enter this figure. Its residuals are calculator differences.
"""

from __future__ import annotations

import json

import matplotlib
import numpy as np

from scripts.audit_garnet_angel_2022 import FOLDER, independent_pressure, source
from scripts.report_plot_themes import dark_svg


def plot():
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    checks = json.loads((FOLDER / "eosfit/checks.json").read_text(encoding="utf-8"))
    with plt.rc_context({"svg.hashsalt": "garnet-angel-2022", "font.size": 10}):
        fig, axes = plt.subplots(
            2, 2, figsize=(10, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex="col"
        )
        for col, name in enumerate(["pyrope", "grossular"]):
            p, _ = source(name)
            volume = np.linspace(0.85, 1.04, 130) * p[1]
            for temperature, color in [(298.0, "#2463a2"), (1000.0, "#be4c20")]:
                pressure = [independent_pressure(p, v, temperature) for v in volume]
                axes[0, col].plot(
                    pressure,
                    volume / p[1],
                    color=color,
                    label=f"{temperature:g} K independent law",
                )
                rows = [
                    r
                    for r in checks["records"][name]
                    if r["temperature_k"] == temperature
                ]
                axes[0, col].scatter(
                    [r["eosfit_pressure_gpa"] for r in rows],
                    [r["volume_cm3_mol"] / p[1] for r in rows],
                    marker="x",
                    color=color,
                    s=35,
                )
                axes[1, col].scatter(
                    [r["eosfit_pressure_gpa"] for r in rows],
                    [r["source_gas_constant_difference_gpa"] for r in rows],
                    color=color,
                    s=23,
                )
            axes[0, col].set_title(name.capitalize() + " — calculated checks")
            axes[0, col].set_ylabel("Volume V/V0")
            axes[0, col].legend(fontsize=8)
            axes[1, col].axhline(0, color="#777777", lw=0.7)
            axes[1, col].set_xlabel("Pressure (GPa)")
            axes[1, col].set_ylabel("EosFit − source-R\ncalculation (GPa)")
            axes[1, col].ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
            for ax in axes[:, col]:
                ax.grid(alpha=0.15)
        fig.suptitle(
            "Source-law curves and EosFit calculator outputs; no experimental observations",
            fontsize=11,
        )
        fig.tight_layout()
        path = FOLDER / "calculated-comparison.svg"
        fig.savefig(path, metadata={"Date": None})
        (FOLDER / "calculated-comparison-dark.svg").write_bytes(
            dark_svg(path.read_bytes())
        )
        plt.close(fig)


if __name__ == "__main__":
    plot()
