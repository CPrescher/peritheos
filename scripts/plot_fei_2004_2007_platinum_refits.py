"""Compare published/diagnostic Pt curves with preserved numerical readouts."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_fei_2007_thermal_parts import PARAMETERS, pressure
from scripts.fei_curve_readouts import READOUTS_SHA256, curve_readouts
from scripts.reproduce_fei_2004_pressure_scales import source_pressure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"
SOURCE_2004 = "https://doi.org/10.1016/j.pepi.2003.09.018"
# Hash recorded when the original PDF was inspected; source images are not bundled.
SOURCE_2004_SHA256 = "b3fb00483e702e48d9d4cb5b4fa1dcefbbe176f94b4eca438d84b99a225342fd"


def main(*, write_report=False, output=None):
    figure2 = json.loads(
        (DATA / "fei-2007-platinum-source-curve-check.json").read_text(encoding="utf-8")
    )
    result2007 = json.loads(
        (DATA / "fei-2007-platinum-reproduction.json").read_text(encoding="utf-8")
    )
    fitted2007 = result2007["fixed_v0_joint_diagnostic"]["parameters"]
    coefficients2007 = list(PARAMETERS["platinum"])
    coefficients2007[2:5] = [fitted2007[k] for k in ("K0_prime", "gamma0", "q")]
    # Plot the preserved diagnostic, without running a new optimizer or changing
    # the historical fit just to replace its presentation.
    result2004 = json.loads(
        (DATA / "fei-2004-2007-platinum-refit-curve-comparison.json").read_text(
            encoding="utf-8"
        )
    )["fei2004_diagnostic"]
    coefficients2004 = result2004["refit_parameters"]
    source_image = "docs/data/fei-2007-thermal-source-figures/fei2004-figure4.jpg"
    readouts = curve_readouts("platinum_2004_figure4")
    pcal2004 = np.array(readouts["pressure_from_x_pixel"])
    vcal2004 = np.array(readouts["volume_from_y_pixel"])
    samples2004 = []
    for row in readouts["points"]:
        temperature, x, y = row["temperature_k"], row["x_pixel"], row["y_pixel"]
        volume = float(np.polyval(vcal2004, y))
        target = float(np.polyval(pcal2004, x))
        samples2004.append(
            {
                "x_pixel": x,
                "y_pixel": y,
                "temperature_k": temperature,
                "volume_a3": volume,
                "figure_pressure_gpa": target,
                "published_minus_figure_gpa": float(
                    source_pressure(volume, temperature, "platinum") - target
                ),
                "diagnostic_minus_figure_gpa": float(
                    source_pressure(volume, temperature, "platinum", coefficients2004)
                    - target
                ),
            }
        )
    curve_summary = {}
    for temperature in (300, 1473, 1873):
        selected = [r for r in figure2["points"] if r["temperature_k"] == temperature]
        volume = np.array([r["volume_a3"] for r in selected])
        target = np.array([r["figure_pressure_gpa"] for r in selected])
        difference = (
            pressure(volume, temperature, "platinum", coefficients2007) - target
        )
        curve_summary[str(temperature)] = {
            "samples": len(selected),
            "diagnostic_rmse_gpa": float(np.sqrt(np.mean(difference**2))),
            "diagnostic_max_absolute_difference_gpa": float(np.max(abs(difference))),
        }
    report = {
        "curve_readouts_sha256": READOUTS_SHA256,
        "replay_basis": "preserved_numerical_pixel_readouts_source_images_not_bundled",
        "fei2007_diagnostic_parameters": fitted2007,
        "fei2007_diagnostic_vs_figure2": curve_summary,
        "fei2004_diagnostic": result2004,
        "fei2004_figure4": {
            "source_pdf": SOURCE_2004,
            "source_pdf_sha256": SOURCE_2004_SHA256,
            "source_image": source_image,
            "source_image_sha256": readouts["source_image_sha256"],
            "pressure_from_x_pixel": pcal2004.tolist(),
            "volume_from_y_pixel": vcal2004.tolist(),
            "samples": samples2004,
            "scope": "Three manually reviewed raster rows per isotherm; readout comparison, not new observations or fit targets.",
        },
        "qualification": "Published coefficients unchanged. Both refits are unweighted pressure diagnostics with assumed constraints; neither is established as the original optimization. High-pressure thermal Figure 2 branches extend beyond the recovered thermal observations.",
    }
    if write_report:
        (DATA / "fei-2004-2007-platinum-refit-curve-comparison.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), layout="constrained")
    setups = [
        (
            axes[0],
            figure2["points"],
            np.linspace(47.05, 60.38, 800),
            (300, 1473, 1873),
            140,
            "2007: Figure 2 readouts",
            lambda v, t, fit: pressure(
                v, t, "platinum", coefficients2007 if fit else None
            ),
        ),
        (
            axes[1],
            samples2004,
            np.linspace(56.05, 60.38, 800),
            (300, 1473, 1673, 1873),
            30,
            "2004: Figure 4 readouts",
            lambda v, t, fit: source_pressure(
                v, t, "platinum", coefficients2004 if fit else None
            ),
        ),
    ]
    colors = {300: "#2525cc", 1473: "#cc00cc", 1673: "#b27300", 1873: "#dd2222"}
    for ax, points, volumes, temperatures, max_pressure, title, evaluate in setups:
        for temperature in temperatures:
            color = colors[temperature]
            for fit, style in ((False, "-"), (True, "--")):
                calculated = evaluate(volumes, temperature, fit)
                inside = (calculated >= 0) & (calculated <= max_pressure)
                ax.plot(
                    calculated[inside],
                    volumes[inside],
                    linestyle=style,
                    color=color,
                    linewidth=1.2,
                    label=f"{temperature} K {'diagnostic' if fit else 'published'}",
                )
            selected = [r for r in points if r["temperature_k"] == temperature]
            ax.scatter(
                [r["figure_pressure_gpa"] for r in selected],
                [r["volume_a3"] for r in selected],
                facecolors="none",
                edgecolors=color,
                s=26,
            )
        ax.set(
            xlabel="Pressure (GPa)",
            ylabel="Pt cell volume ($\\AA^3$)",
            xlim=(0, max_pressure),
            ylim=(volumes.min(), volumes.max()),
        )
        ax.grid(alpha=0.15)
        ax.legend(fontsize=8)
        ax.set_title(title, fontsize=13)
    fig.suptitle("Published Pt EOS versus our diagnostic refits", fontsize=17)
    fig.supxlabel(
        "Open circles: preserved curve readouts | Solid lines: published coefficients | Dashed lines: diagnostic refits\n"
        "2007: q = -1.661; gamma0 = 2.612 | 2004: q = -0.977; gamma0 = 2.559. The source figures have different pressure ranges.",
        fontsize=11,
    )
    fig.savefig(
        output or DATA / "fei-2004-2007-platinum-refit-curve-comparison.png", dpi=180
    )
    plt.close(fig)
    print(
        json.dumps(
            {
                "fei2007": curve_summary,
                "fei2004_table_output_max_difference_gpa": result2004[
                    "table_output_max_difference_gpa"
                ],
            },
            indent=2,
        )
    )

    return report


if __name__ == "__main__":
    main()
