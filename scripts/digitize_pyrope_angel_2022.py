"""Extract factual Fig. 8 vectors from the original Angel PDF (PyMuPDF needed).

Repository maintenance only. Supply --pdf; no network access or source figures
are redistributed. These are digitized plotted observations, not author inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.audit_pyrope_angel_2022_observations import DATA


def extract(pdf: Path, output: Path):
    import fitz

    page = fitz.open(pdf)[10]
    drawings = page.get_drawings()
    # Axis major ticks: 20,15,10,5,0 GPa and 260..170 GPa, in PDF points.
    x = np.array([v[1].x for v in drawings[367]["items"][:5]])
    y = np.array([v[1].y for v in drawings[371]["items"][:10]])
    p_axis = np.polyfit(x, [20, 15, 10, 5, 0], 1)
    k_axis = np.polyfit(y, np.arange(260, 160, -10), 1)
    # Validate the identified paths, rather than accepting arbitrary PDFs.
    original = drawings[421:454:2]
    scaled = drawings[571:604:2]
    bars = drawings[454]["items"][1:]
    if len(original) != 17 or len(scaled) != 17 or len(bars) != 68:
        raise ValueError("Unexpected Fig. 8 paths")
    rows = []
    for i, (raw, adjusted) in enumerate(zip(original, scaled)):
        if raw["fill"] != (1, 1, 1) or adjusted["fill"] != (0, 0, 0):
            raise ValueError("Unexpected Fig. 8 marker colors")
        r, a = raw["rect"], adjusted["rect"]
        if abs(r.width - 3.081) > 0.002 or abs(r.height - 3.081) > 0.002:
            raise ValueError("Unexpected Fig. 8 marker shape")
        xp, yp = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        group = bars[4 * (16 - i) : 4 * (17 - i)]
        # Each error bar is four segments (two caps, two stems); paths reverse
        # the marker order. This pairing preserves overlapping measurements.
        if any(abs((g[1].x + g[2].x) / 2 - xp) > 0.003 for g in group):
            raise ValueError("Error-bar/marker association failed")
        cap1, cap2 = group[0][1].y, group[2][1].y
        rows.append(
            {
                "plot_row": i + 1,
                "pdf_marker_path": 421 + 2 * i,
                "pdf_x_pt": xp,
                "pdf_y_pt": yp,
                "pressure_gpa": float(np.polyval(p_axis, xp)),
                "ks_original_gpa": float(np.polyval(k_axis, yp)),
                "ks_error_bar_gpa": abs(cap1 - cap2) * abs(k_axis[0]) / 2,
                "ks_scaled_in_figure_gpa": float(np.polyval(k_axis, (a.y0 + a.y1) / 2)),
            }
        )
    output.mkdir(parents=True, exist_ok=True)
    with (output / "sinogeikin-2000-digitized-angel-figure-8.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "kind": "digitized_observations",
        "source": "Angel et al. 2022, Fig. 8; original open black squares attributed to Sinogeikin and Bass 2000",
        "source_doi": "10.1007/s00410-022-01918-5",
        "original_study_doi": "10.1016/S0031-9201(00)00143-6",
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "pdf_page_one_based": 11,
        "axis_P_from_pdf_x": p_axis.tolist(),
        "axis_KS_from_pdf_y": k_axis.tolist(),
        "max_axis_tick_residual_gpa": {
            "P": float(np.max(abs(np.polyval(p_axis, x) - [20, 15, 10, 5, 0]))),
            "KS": float(np.max(abs(np.polyval(k_axis, y) - np.arange(260, 160, -10)))),
        },
        "digitization_resolution_bounds_gpa": {"P": 0.01, "KS": 0.05},
        "resolution_qualification": "Conservative coordinate bounds from vector/tick quantization, not experimental standard deviations. Extra saved digits retain the deterministic extraction, not source precision.",
        "errors": "Vertical half error-bar spans digitized from original open-square series. No horizontal pressure bars identified. Confidence convention and original covariance not recovered.",
        "scaling_evidence": "All 17 markers, including both ambient-pressure points, also appear scaled in Fig. 8. This resolves the displayed mask only; S2 describes high-pressure scaling and original fit groups remain unavailable.",
    }
    (output / "sinogeikin-2000-digitization.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    compare_plotted_inputs(drawings, output)
    return report


def compare_plotted_inputs(drawings, output):
    """Use published plots to check column choices, not to optimize fit inputs."""
    from scripts.refit_pyrope_angel_2022_expanded import expanded_rows

    px = np.polyfit(
        [v[1].x for v in drawings[3]["items"][:8]], np.arange(35, -1, -5), 1
    )
    vy = np.polyfit(
        [v[1].y for v in drawings[7]["items"][:9]], np.arange(114, 97, -2), 1
    )
    rows = expanded_rows()
    checked = {}
    for name, indices in (
        ("du_pv", range(122, 131)),
        ("milani_pv", range(131, 142)),
        ("chantel_pv", range(142, 152)),
    ):
        data = [r for r in rows if r["dataset"] == name]
        points = []
        for r, index in zip(data, indices):
            box = drawings[index]["rect"]
            p = float(np.polyval(px, (box.x0 + box.x1) / 2))
            v = float(np.polyval(vy, (box.y0 + box.y1) / 2))
            points.append(
                {
                    "row": r["row"],
                    "figure_P_gpa": p,
                    "figure_V_cm3_mol": v,
                    "figure_minus_prepared_P_gpa": p - r["P"],
                    "figure_minus_prepared_V_cm3_mol": v - r["Y"],
                }
            )
        checked[name] = points
    result = {
        "source": "Angel 2022 Fig. 7 vector marker centers",
        "Du_preprint_check": "All nine preprint PV points match plotted pressures to 0.02 GPa and normalized volumes to 0.004 cm3/mol. This corroborates the recovered preprint table; it does not turn the preprint into the final primary version.",
        "Chantel_pressure_choice": "All ten plot abscissae match the calibrated (Au) pressure column to 0.02 GPa, not the internally derived pressure column. S2's 25.2 GPa endpoint alone would select the wrong column.",
        "normalization_limit": "Chantel common ambient-density normalization differs from plotted V by about 0.045 cm3/mol. Milani normalized volumes differ by up to 0.008 cm3/mol and recalibrated pressures by up to 0.042 GPa. Retained as unresolved original normalization/calibration details; no values/anchors tuned to fit the figure or the published coefficients.",
        "points": checked,
    }
    (output / "expanded-figure-crosschecks.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA)
    args = parser.parse_args()
    extract(args.pdf, args.output)
