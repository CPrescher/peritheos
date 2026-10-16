"""Check Zhao Figure 2 against the retained table and published EOS.

Numerical digitization uses the 450 dpi rendering of page 26. The graphical
pressure is a separate approximate observation, never a correction to Table 1.
The optional comparison PNG includes an attributed source crop for inspection;
it is not part of the contributor CC0 numerical transcription.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from scripts.diopside_correct_dataset import corrected_rows
from scripts.diopside_reproduce import (
    ROOT,
    THERMAL,
    load_rows,
    thermal_pressure,
    values,
)

DEFAULT_REPORT = ROOT / "docs/data/diopside-plot-check.json"
CROP_BOX = (1610, 2279, 2968, 3335)
# Crop-relative pixel positions read from Figure 2. The lowest marker cluster
# near V=431.7 identifies source row 15 independently of the printed pressure.
X_TICKS = np.array([201, 371, 541, 711, 879, 1050, 1221], dtype=float)
P_TICKS = np.arange(0.0, 14.0, 2.0)
Y_TICKS = np.array([310, 451, 592, 733, 874], dtype=float)
V_TICKS = np.array([440, 430, 420, 410, 400], dtype=float)
ROW15_CENTRE = (590.0, 428.0)
AMBIENT_CENTRE = (201.0, 327.0)
CURVE_X = 1160.0
CURVE_Y = np.array([832.5, 815.5, 793.5, 766.5, 737.5, 704.5])
CURVE_T = np.array([300, 473, 673, 873, 1073, 1273])


def inspect_plots():
    p_map = np.polyfit(X_TICKS, P_TICKS, 1)
    v_map = np.polyfit(Y_TICKS, V_TICKS, 1)
    graph_p = float(np.polyval(p_map, ROW15_CENTRE[0]))
    graph_v = float(np.polyval(v_map, ROW15_CENTRE[1]))
    curve_p = float(np.polyval(p_map, CURVE_X))
    digitized_v = np.polyval(v_map, CURVE_Y)
    model_v = np.array(
        [brentq(lambda v: thermal_pressure(v, t) - curve_p, 350, 500) for t in CURVE_T]
    )
    rows = load_rows("diopside-zhao-1998-table1")
    v, t = (values(rows, column) for column in ("cell_volume_a3", "temperature_k"))
    # Use exactly the rounded pressure exported for subsequent dataset reuse.
    graphical_interpretation_p = values(corrected_rows(graph_p), "pressure_gpa")
    fit = least_squares(
        lambda x: thermal_pressure(v, t, x) - graphical_interpretation_p,
        THERMAL,
        x_scale=np.abs(THERMAL),
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    residual = thermal_pressure(v, t) - graphical_interpretation_p
    parameter_checks = {
        name: {
            "published_value": float(published),
            "reported_error": float(error),
            "independent_refit_value": float(fitted),
            "within_reported_error": bool(abs(fitted - published) <= error),
        }
        for name, published, error, fitted in zip(
            ["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"],
            THERMAL,
            [4.0, 0.6, 0.05e-5, 0.07e-8, 0.0004],
            fit.x,
        )
    }
    return {
        "source": "Zhao et al. (1998), DOI 10.4131/jshpreview.7.25, Figures 2-4, pages 26-27",
        "source_pdf_sha256": "946740290bc4533d4a09895d2f0f06e0869a3342ff468f47d3450925a67a4794",
        "method": "Manual tick/marker centre readings on the unmodified Figure 2 crop from a 450 dpi page-26 rendering; affine least-squares axis calibration. Coordinates are analysis pixels, not original experimental observations.",
        "render_dpi": 450,
        "render_size_pixels": [3157, 4513],
        "crop_box_pixels": list(CROP_BOX),
        "pressure_tick_pixels": X_TICKS.tolist(),
        "pressure_tick_values_gpa": P_TICKS.tolist(),
        "volume_tick_pixels": Y_TICKS.tolist(),
        "volume_tick_values_a3": V_TICKS.tolist(),
        "parity_assessment": {
            "assessment_date": "2026-10-16",
            "equation_mapping_validated": True,
            "published_curves_agree_at_graphical_resolution": bool(
                np.max(np.abs(digitized_v - model_v)) < 0.3
            ),
            "parameter_parity_with_documented_table_error": all(
                check["within_reported_error"] for check in parameter_checks.values()
            ),
            "parameter_checks": parameter_checks,
            "basis": "Five varied coefficients agree within their reported errors when the independently digitized Figure 2 pressure is used for row 15, with V0=438.67 A^3 and Tr=300 K fixed. This is published-model and parameter parity with a documented source-table error, rather than an exact recovery of the original numerical inputs and regression workflow.",
            "source_table_errors_block_published_eos_acceptance": False,
            "original_observations_preserved": True,
            "published_coefficients_preserved": True,
            "exact_original_author_workflow_recovered": False,
        },
        "row15": {
            "marker_centre_pixels": list(ROW15_CENTRE),
            "digitized_pressure_gpa": graph_p,
            "pressure_reading_bound_gpa": 0.05,
            "digitized_volume_a3": graph_v,
            "volume_reading_bound_a3": 0.2,
            "bound_convention": "Conservative graphical reading bounds, not statistical experimental errors or fitted uncertainty.",
            "table_pressure_gpa": 6.95,
            "table_volume_a3": 431.66,
            "table_temperature_k": 1072,
            "published_eos_pressure_gpa": float(thermal_pressure(431.66, 1072)),
            "finding": "Figure 2 places the matching marker near 4.6 GPa, while the Table 1 pressure would fall in an empty part of the plot. Supports a printed-table error; exact original pressure is not recovered.",
        },
        "ambient": {
            "marker_centre_pixels": list(AMBIENT_CENTRE),
            "digitized_pressure_gpa": float(np.polyval(p_map, AMBIENT_CENTRE[0])),
            "digitized_volume_a3": float(np.polyval(v_map, AMBIENT_CENTRE[1])),
            "volume_reading_bound_a3": 0.2,
            "table_volume_a3": 438.67,
            "table_lattice_metric_volume_a3": 443.5752871851424,
            "finding": "The plotted ambient point corroborates the printed volume column near 438.7 A^3; it does not recover the correct individual lattice dimensions.",
        },
        "figure2_curves": {
            "cross_section_pressure_gpa": curve_p,
            "cross_section_x_pixel": CURVE_X,
            "temperature_labels_k": CURVE_T.tolist(),
            "curve_centre_y_pixels": CURVE_Y.tolist(),
            "digitized_volumes_a3": digitized_v.tolist(),
            "independent_model_volumes_a3": model_v.tolist(),
            "digitized_minus_model_a3": (digitized_v - model_v).tolist(),
            "finding": "All six curves have the expected ordering and agree within 0.3 A^3 at this graphical cross section. This supports graphical consistency, not numerical identity with the original plotting code. The labelled temperatures are 300, 473, 673, 873, 1073 and 1273 K; the caption's statement of 200 K steps above 300 K is inexact.",
        },
        "figure3": {
            "finding": "Volume-temperature plot: raw points and pressure-corrected series are shown, but raw-point pressures are not encoded. It cannot independently determine row 15 pressure or recover the ambient lattice dimensions.",
        },
        "figure4": {
            "finding": "The plotted zero-pressure mean expansivity near 3.82e-5 K^-1 is consistent with the adopted alpha(T) coefficients. Figure slope label -1.54(2)e-6 and prose -1.52(2)e-6 differ; this remains a separate secondary-regression discrepancy.",
        },
        "graphical_interpretation_diagnostic": {
            "selection": "All 32 original rows; row 15 uses the Figure 2 pressure rounded to 4.58 GPa, exactly as in the separately bundled corrected P-V-T dataset. Original source tables and published coefficients are untouched.",
            "dataset_identifier": "diopside_zhao_1998_table1_corrected_pvt",
            "row15_analysis_pressure_gpa": float(graphical_interpretation_p[14]),
            "published_curve_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "refit_coefficients": dict(
                zip(["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"], map(float, fit.x))
            ),
            "refit_pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
            "fixed_V0_a3": 438.67,
            "fixed_Tr_k": 300,
            "objective": "Unweighted pressure residuals for sensitivity only; graphical bound is not an experimental row uncertainty.",
            "original_author_workflow_parity": False,
        },
    }


def comparison_image(page_render, output, result):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    page = Image.open(page_render)
    if page.size != (3157, 4513):
        raise ValueError("Expected page 26 at 450 dpi")
    crop = page.crop(CROP_BOX)
    # Pixel-to-axes transforms are fitted from observed ticks, rather than from
    # the scanned frame corners, which are slightly misaligned.
    p_map = np.polyfit(X_TICKS, P_TICKS, 1)
    v_map = np.polyfit(Y_TICKS, V_TICKS, 1)
    fig, (source_ax, model_ax) = plt.subplots(1, 2, figsize=(14, 5.9))
    source_ax.imshow(crop)
    source_ax.set_axis_off()
    source_ax.set_title("Published Figure 2", loc="left", fontweight="bold")
    xy_source = ROW15_CENTRE
    xy_table = (
        (6.95 - p_map[1]) / p_map[0],
        (431.66 - v_map[1]) / v_map[0],
    )
    source_ax.scatter(
        *xy_source, s=180, facecolors="none", edgecolors="#007a5e", linewidths=2
    )
    source_ax.scatter(*xy_table, s=110, marker="x", color="#b42232", linewidths=2)
    source_ax.annotate(
        "Plotted point: about 4.6 GPa",
        xy_source,
        xytext=(400, 110),
        arrowprops={"arrowstyle": "->", "color": "#007a5e"},
        fontsize=10,
        color="#00664d",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    source_ax.annotate(
        "Table says 6.95 GPa\nNo marker here",
        xy_table,
        xytext=(825, 320),
        arrowprops={"arrowstyle": "->", "color": "#b42232"},
        fontsize=10,
        color="#b42232",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.95},
    )
    rows = load_rows("diopside-zhao-1998-table1")
    p, v = (values(rows, c) for c in ("pressure_gpa", "cell_volume_a3"))
    curve_p = np.linspace(0, 12, 160)
    for t in CURVE_T:
        volumes = [
            brentq(lambda volume: thermal_pressure(volume, t) - pressure, 350, 500)
            for pressure in curve_p
        ]
        model_ax.plot(curve_p, volumes, color="0.55", lw=0.9)
    mask = np.arange(len(p)) != 14
    model_ax.scatter(
        p[mask],
        v[mask],
        s=24,
        facecolors="white",
        edgecolors="black",
        label="Unchanged Table 1 observations",
        zorder=3,
    )
    model_ax.scatter(
        6.95,
        431.66,
        s=80,
        marker="x",
        color="#b42232",
        label="Row 15 as printed",
        zorder=4,
    )
    graph_p = result["row15"]["digitized_pressure_gpa"]
    model_ax.errorbar(
        graph_p,
        431.66,
        xerr=0.05,
        marker="o",
        markersize=7,
        color="#007a5e",
        capsize=3,
        label="Row 15 graphical interpretation",
        zorder=4,
    )
    model_ax.annotate(
        "Published EOS gives 4.563 GPa\nat V = 431.66 A³, T = 1072 K",
        xy=(graph_p, 431.66),
        xytext=(5.9, 440),
        arrowprops={"arrowstyle": "->", "color": "#007a5e"},
        color="#00664d",
        fontsize=10,
    )
    model_ax.set(
        xlim=(-0.2, 12),
        ylim=(400, 451),
        xlabel="Pressure (GPa)",
        ylabel="Conventional cell volume (Å³)",
    )
    model_ax.set_title("Independent reconstruction", loc="left", fontweight="bold")
    model_ax.grid(alpha=0.16)
    model_ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    fig.suptitle(
        "Diopside: the plot supports a pressure-table error",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.02,
        0.025,
        "Source: Zhao et al. (1998), DOI 10.4131/jshpreview.7.25, Fig. 2. Source crop annotated for analysis. Original table values remain unchanged.",
        fontsize=8,
        color="0.3",
    )
    fig.subplots_adjust(left=0.025, right=0.985, bottom=0.14, top=0.87, wspace=0.13)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--page-render", type=Path)
    parser.add_argument("--image", type=Path)
    args = parser.parse_args()
    result = inspect_plots()
    if args.source_pdf:
        result["source_pdf_sha256"] = hashlib.sha256(
            args.source_pdf.read_bytes()
        ).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    if args.image:
        if args.page_render is None:
            parser.error("--image requires --page-render")
        comparison_image(args.page_render, args.image, result)
    print(
        json.dumps(
            {"row15": result["row15"], "curves": result["figure2_curves"]}, indent=2
        )
    )


if __name__ == "__main__":
    main()
