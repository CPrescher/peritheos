"""Overlay published and diagnostic Pt curves on both Fei source figures."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts.audit_fei_2007_thermal_parts import PARAMETERS, pressure
from scripts.reproduce_fei_2004_pressure_scales import diagnostic, source_pressure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"
SOURCE_2004 = "https://doi.org/10.1016/j.pepi.2003.09.018"
# Hash recorded when the original PDF was inspected; plotting uses its bundled image.
SOURCE_2004_SHA256 = "b3fb00483e702e48d9d4cb5b4fa1dcefbbe176f94b4eca438d84b99a225342fd"


def main():
    figure2 = json.loads(
        (DATA / "fei-2007-platinum-source-curve-check.json").read_text(encoding="utf-8")
    )
    result2007 = json.loads(
        (DATA / "fei-2007-platinum-reproduction.json").read_text(encoding="utf-8")
    )
    fitted2007 = result2007["fixed_v0_joint_diagnostic"]["parameters"]
    coefficients2007 = list(PARAMETERS["platinum"])
    coefficients2007[2:5] = [fitted2007[k] for k in ("K0_prime", "gamma0", "q")]
    result2004 = diagnostic("platinum")
    coefficients2004 = result2004["refit_parameters"]
    source_image = DATA / "fei-2007-thermal-source-figures/fei2004-figure4.jpg"
    pcal2004 = np.polyfit(
        [107, 189.5, 273, 355, 437, 521, 603, 685], range(0, 29, 4), 1
    )
    vcal2004 = np.polyfit(
        [10, 62.5, 116.5, 169, 222, 274, 328, 381, 434, 488, 540.5],
        np.arange(61, 55.9, -0.5),
        1,
    )
    # Manually reviewed rows without overlapping data markers. Adjacent red
    # pixel components are joined across <=2-pixel JPEG gaps; ordered curves
    # are 300, 1473, 1673 and 1873 K, independently of model predictions.
    image2004 = plt.imread(source_image)
    samples2004 = []
    for y in (225, 350, 400):
        red, green, blue = image2004[y].astype(float).T
        xs = np.flatnonzero(
            (red > 150) & (green < 140) & (blue < 140) & (red - green > 60)
        )
        xs = xs[(xs > 120) & (xs < 730)]
        groups = np.split(xs, np.flatnonzero(np.diff(xs) > 2) + 1)
        assert len(groups) == 4
        for temperature, pixels in zip((300, 1473, 1673, 1873), groups):
            x = float(pixels.mean())
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
                        source_pressure(
                            volume, temperature, "platinum", coefficients2004
                        )
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
        "fei2007_diagnostic_parameters": fitted2007,
        "fei2007_diagnostic_vs_figure2": curve_summary,
        "fei2004_diagnostic": result2004,
        "fei2004_figure4": {
            "source_pdf": SOURCE_2004,
            "source_pdf_sha256": SOURCE_2004_SHA256,
            "source_image": str(source_image.relative_to(ROOT)),
            "source_image_sha256": hashlib.sha256(
                source_image.read_bytes()
            ).hexdigest(),
            "pressure_from_x_pixel": pcal2004.tolist(),
            "volume_from_y_pixel": vcal2004.tolist(),
            "samples": samples2004,
            "scope": "Three manually reviewed raster rows per isotherm; readout comparison, not new observations or fit targets.",
        },
        "qualification": "Published coefficients unchanged. Both refits are unweighted pressure diagnostics with assumed constraints; neither is established as the original optimization. High-pressure thermal Figure 2 branches extend beyond the recovered thermal observations.",
    }
    (DATA / "fei-2004-2007-platinum-refit-curve-comparison.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), layout="constrained")
    setups = [
        (
            axes[0],
            plt.imread(ROOT / figure2["source_image"]),
            figure2["linear_calibration"]["pressure_from_x"],
            figure2["linear_calibration"]["volume_from_y"],
            np.linspace(47.05, 60.38, 800),
            (300, 1473, 1873),
            140,
            "2007: Figure 2",
            lambda v, t, fit: pressure(
                v, t, "platinum", coefficients2007 if fit else None
            ),
        ),
        (
            axes[1],
            image2004,
            pcal2004,
            vcal2004,
            np.linspace(56.05, 60.38, 800),
            (300, 1473, 1673, 1873),
            30,
            "2004: Figure 4",
            lambda v, t, fit: source_pressure(
                v, t, "platinum", coefficients2004 if fit else None
            ),
        ),
    ]
    for (
        ax,
        original,
        pc,
        vc,
        volumes,
        temperatures,
        max_pressure,
        title,
        evaluate,
    ) in setups:
        ax.imshow(original)
        for fit, style, color in ((False, "--", "black"), (True, "--", "#00835e")):
            for temperature in temperatures:
                calculated = evaluate(volumes, temperature, fit)
                x = (calculated - pc[1]) / pc[0]
                y = (volumes - vc[1]) / vc[0]
                inside = (calculated >= 0) & (calculated <= max_pressure)
                ax.plot(
                    x[inside],
                    y[inside],
                    linestyle=style,
                    color=color,
                    linewidth=1.0 if not fit else 1.6,
                    dashes=(5, 3) if not fit else (2, 2),
                )
        ax.set_xlim(0, original.shape[1])
        ax.set_ylim(original.shape[0], 0)
        ax.axis("off")
        ax.set_title(title, fontsize=14)
    fig.suptitle("Published Pt EOS versus our diagnostic refits", fontsize=17)
    fig.supxlabel(
        "Original curves: source colors | Black dashed: published coefficients | Green dashed: our diagnostic refit\n"
        "2007: q = -1.661; gamma0 = 2.612 | 2004: q = -0.977; gamma0 = 2.559. The source figures have different pressure ranges.",
        fontsize=11,
    )
    fig.savefig(DATA / "fei-2004-2007-platinum-refit-curve-comparison.png", dpi=180)
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


if __name__ == "__main__":
    main()
