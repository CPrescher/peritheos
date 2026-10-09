"""Published Matsui/Zhu isotherm offsets from Fei 2007 at 100-150 GPa."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.optimize import brentq

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.compare_pt_thermal_eos_with_fei import MODELS, predict

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/pt-isotherm-differences-100-150gpa"
TEMPERATURES = (300, 1000, 1873, 2000, 3000)
ALTERNATIVES = ("Matsui 2009", "Zhu v3")


def volume_at(name, pressure, temperature):
    volume = brentq(
        lambda v: float(predict(name, v, temperature)) - pressure,
        0.5 * 60.38,
        1.2 * 60.38,
        xtol=1e-12,
    )
    assert abs(float(predict(name, volume, temperature)) - pressure) < 1e-9
    return volume


def main():
    rows = []
    pressure_grid = np.linspace(100, 150, 101)
    curves = {}
    for temperature in TEMPERATURES:
        volume_grid = np.array(
            [volume_at("Fei 2007", p, temperature) for p in pressure_grid]
        )
        curves[str(temperature)] = {
            name: (predict(name, volume_grid, temperature) - pressure_grid).tolist()
            for name in ALTERNATIVES
        }
        for pressure in (100, 125, 150):
            fei_volume = volume_at("Fei 2007", pressure, temperature)
            for name in ALTERNATIVES:
                predicted = float(predict(name, fei_volume, temperature))
                own_volume = volume_at(name, pressure, temperature)
                rows.append(
                    {
                        "temperature_k": temperature,
                        "fei_pressure_gpa": pressure,
                        "fei_volume_a3": fei_volume,
                        "alternative": name,
                        "alternative_pressure_at_fei_volume_gpa": predicted,
                        "pressure_difference_at_same_volume_gpa": predicted - pressure,
                        "alternative_volume_at_same_pressure_a3": own_volume,
                        "volume_difference_at_same_pressure_percent": 100
                        * (own_volume / fei_volume - 1),
                    }
                )
    report = {
        "reference": MODELS["Fei 2007"],
        "alternatives": {n: MODELS[n] for n in ALTERNATIVES},
        "convention": "For each T and Fei reference pressure, invert Fei 2007 for V, then evaluate alternative pressure at that same V,T. Delta P = P_alternative - P_Fei. Volume differences are independently calculated at equal P,T.",
        "scope": "Published coefficients only, no refits. Zhu is the released-v3 optimizer/property set. High-temperature Fei predictions at 100-150 GPa extend beyond the recovered hot measurement pressures.",
        "checkpoints": rows,
        "fei_pressure_grid_gpa": pressure_grid.tolist(),
        "pressure_difference_curves_gpa": curves,
    }
    OUTPUT.with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    with OUTPUT.with_suffix(".csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(
        1, 2, figsize=(11.5, 4.8), layout="constrained", sharey=True
    )
    colors = ("#252525", "#007b98", "#a43685", "#d26917", "#bf2828")
    for ax, name in zip(axes, ALTERNATIVES):
        for t, color in zip(TEMPERATURES, colors):
            ax.plot(pressure_grid, curves[str(t)][name], color=color, label=f"{t} K")
        ax.axhline(0, color="0.5", linewidth=0.8)
        ax.set(
            title=f"{name} minus Fei 2007",
            xlabel="Fei reference pressure at this temperature (GPa)",
            xlim=(100, 150),
        )
        ax.grid(alpha=0.18)
    axes[0].set_ylabel("Pressure difference at the same V, T (GPa)")
    axes[1].legend(title="Isotherm", fontsize=9)
    fig.suptitle("Pt isotherm differences at 100-150 GPa", fontsize=14)
    fig.supxlabel(
        "Published coefficients; Fei sets the volume separately for each temperature. Zhu uses released-v3 optimizer constants.",
        fontsize=9,
    )
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=180)
    plt.close(fig)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
