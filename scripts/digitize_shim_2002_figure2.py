#!/usr/bin/env python3
"""Recover six Figure 2a observations from the author-hosted vector PDF.

Usage (PyMuPDF is an optional extraction dependency):
    uv run --with pymupdf python scripts/digitize_shim_2002_figure2.py PDF CSV

Figure 2 contains a white cubic marker over a black marker at the same state.
Only the final white marker and its replacement error bars are included.
Legend symbols, comparison curves and the c/a panel are excluded.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

SOURCE_SHA256 = "ad622245ea2b6b1171dcd533139be36362bf838fd7f9f3d1237b20f8c7739cdf"
PRESSURE_AXIS = ((328.29400634765625, 0.0), (537.4310302734375, 50.0))
RATIO_AXIS = ((481.09600830078125, 1.0), (574.0460205078125, 0.86))
ADOPTED_V0_A3 = 45.58
COORDINATE_BOUND_PT = 0.1
# Paragraph 7 explicitly prints all six pressures. Prefer these source values
# for the P-V observations; preserve the graph-derived pressures separately.
REPORTED_PRESSURES_GPA = (19.7, 24.2, 25.2, 26.6, 36.1, 45.8)


def linear(value, calibration):
    (position0, value0), (position1, value1) = calibration
    return value0 + (value - position0) * (value1 - value0) / (position1 - position0)


def extract_rows(source_pdf: Path) -> list[dict[str, object]]:
    import fitz

    digest = hashlib.sha256(source_pdf.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(f"unexpected source PDF SHA-256: {digest}")
    with fitz.open(source_pdf) as document:
        paths = document[1].get_drawings()
    if len(paths[416]["items"]) != 36 or len(paths[423]["items"]) != 6:
        raise ValueError("Figure 2 error-bar path inventory changed")

    pressure_scale = abs(50.0 / (PRESSURE_AXIS[1][0] - PRESSURE_AXIS[0][0]))
    volume_scale = abs(0.14 / (RATIO_AXIS[1][0] - RATIO_AXIS[0][0])) * ADOPTED_V0_A3
    # Marker path, six-line error-bar group, source refinement.
    selections = [(424, paths[423]["items"], "cubic")]
    selections += [
        (marker, paths[416]["items"][6 * group : 6 * (group + 1)], "tetragonal")
        for marker, group in [(419, 2), (418, 1), (420, 3), (421, 4), (422, 5)]
    ]
    rows = []
    for order, (marker, bars, refinement) in enumerate(selections, 1):
        path = paths[marker]
        rect = path["rect"]
        expected_fill = (1.0, 1.0, 1.0) if refinement == "cubic" else (0.0, 0.0, 0.0)
        if (
            path["fill"] != expected_fill
            or len(path["items"]) != 4
            or abs(rect.width - 3.873) > 0.002
            or abs(rect.height - 3.873) > 0.002
            or any(item[0] != "l" for item in bars)
        ):
            raise ValueError(f"Figure 2 marker/bar geometry changed at path {marker}")
        x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
        ratio = linear(y, RATIO_AXIS)
        # Last two lines are the vertical and horizontal error-bar stems.
        sigma_p = abs(bars[5][2].x - bars[5][1].x) * pressure_scale / 2
        sigma_v = abs(bars[4][2].y - bars[4][1].y) * volume_scale / 2
        rows.append(
            {
                "source_order": order,
                "pressure_gpa": f"{REPORTED_PRESSURES_GPA[order - 1]:.1f}",
                "plot_pressure_gpa": f"{linear(x, PRESSURE_AXIS):.6f}",
                "volume_ratio": f"{ratio:.9f}",
                "volume_a3": f"{ADOPTED_V0_A3 * ratio:.6f}",
                "pressure_plot_sigma_gpa": f"{sigma_p:.6f}",
                "volume_plot_sigma_a3": f"{sigma_v:.6f}",
                "pressure_digitization_bound_gpa": f"{COORDINATE_BOUND_PT * pressure_scale:.6f}",
                "volume_digitization_bound_a3": f"{COORDINATE_BOUND_PT * volume_scale:.6f}",
                "plot_x_pt": f"{x:.6f}",
                "plot_y_pt": f"{y:.6f}",
                "source_refinement": refinement,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_pdf", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    rows = extract_rows(args.source_pdf)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} Figure 2a observations to {args.output_csv}")


if __name__ == "__main__":
    main()
