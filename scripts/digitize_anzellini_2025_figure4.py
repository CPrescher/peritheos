#!/usr/bin/env python3
"""Extract the missing 300 K and 830 K Ir P-V markers from Figure 4.

The publisher PDF stores the marker circles as vector paths.  This script reads
their centers directly with pdfplumber and applies the linear axis calibration
encoded by the frame and major ticks.  It intentionally excludes fitted curves,
other temperature series, ambient-pressure literature triangles, and legend
symbols.

Usage:
    python scripts/digitize_anzellini_2025_figure4.py SOURCE.pdf OUTPUT.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import pdfplumber

SOURCE_SHA256 = "3c531753df70d77a32433c6663dc54b6796650d9412a2edbd981b130ea91aa01"
PAGE_INDEX = 2

# PDF-point axis calibration from the Figure 4 vector frame and major ticks.
PRESSURE_AXIS = ((346.745, 0.0), (552.845, 150.0))
VOLUME_AXIS = ((248.247, 44.0), (50.643, 59.0))

BLACK = (0.0, 0.0, 0.0)
BROWN_830_K = (0.4, 0.0, 0.0)


def _linear(value: float, calibration: tuple[tuple[float, float], ...]) -> float:
    (pixel0, value0), (pixel1, value1) = calibration
    return value0 + (value - pixel0) * (value1 - value0) / (pixel1 - pixel0)


def _marker_center(curve: dict) -> tuple[float, float]:
    return (curve["x0"] + curve["x1"]) / 2.0, (curve["top"] + curve["bottom"]) / 2.0


def _is_marker_circle(curve: dict) -> bool:
    width = curve["x1"] - curve["x0"]
    height = curve["bottom"] - curve["top"]
    return (
        curve.get("fill") is True
        and 1.80 < width < 1.90
        and 1.80 < height < 1.90
        and 346.0 < curve["x0"] < 553.0
        and 50.0 < curve["top"] < 249.0
    )


def extract_rows(source_pdf: Path) -> list[dict[str, object]]:
    digest = hashlib.sha256(source_pdf.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(
            f"unexpected source PDF SHA-256 {digest}; expected {SOURCE_SHA256}"
        )

    with pdfplumber.open(source_pdf) as document:
        curves = document.pages[PAGE_INDEX].curves

    black_markers = [
        curve
        for curve in curves
        if _is_marker_circle(curve)
        and curve.get("non_stroking_color") == BLACK
        and curve.get("stroking_color") == BLACK
    ]
    brown_markers = [
        curve
        for curve in curves
        if _is_marker_circle(curve)
        and curve.get("non_stroking_color") == BROWN_830_K
        and curve.get("stroking_color") in {BLACK, BROWN_830_K}
    ]

    if len(black_markers) != 91 or len(brown_markers) != 16:
        raise ValueError(
            "Figure 4 marker inventory changed: expected 91 black 300 K and "
            f"16 brown 830 K markers, found {len(black_markers)} and "
            f"{len(brown_markers)}"
        )

    rows: list[dict[str, object]] = []
    series = [
        (
            300,
            black_markers,
            "Monteseguro et al. (2019) 300 K observations replotted",
            "filled black circle",
        ),
        (
            830,
            brown_markers,
            "Anzellini et al. (2025) present-study 830 K isotherm",
            "filled dark-red circle",
        ),
    ]
    plot_order = 0
    for temperature, markers, source_series, source_marker in series:
        centers = sorted(
            (_marker_center(curve) for curve in markers), key=lambda x: x[0]
        )
        for series_order, (plot_x, plot_y) in enumerate(centers, start=1):
            plot_order += 1
            rows.append(
                {
                    "plot_order": plot_order,
                    "series_order": series_order,
                    "temperature_k": temperature,
                    "pressure_gpa": f"{_linear(plot_x, PRESSURE_AXIS):.9f}",
                    "volume_a3_conventional_cell": f"{_linear(plot_y, VOLUME_AXIS):.9f}",
                    "plot_x_pt": f"{plot_x:.6f}",
                    "plot_y_pt": f"{plot_y:.6f}",
                    "pressure_digitization_uncertainty_gpa": "0.75",
                    "volume_digitization_uncertainty_a3": "0.08",
                    "source_series": source_series,
                    "source_marker": source_marker,
                    "used_in_reconstructed_fit": 1,
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

    print(f"wrote {len(rows)} Figure 4 markers to {args.output_csv}")


if __name__ == "__main__":
    main()
