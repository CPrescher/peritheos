#!/usr/bin/env python3
"""Recover Figure 2 vector marker centers and curve vertices from Dewaele 2019.

Requires pdfplumber; the source PDF is supplied by the caller and is not bundled.
Coordinates remain plot-derived observations, not original experimental rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pdfplumber

PDF_SHA256 = "8ae1bb7611c6925a6a66ac39ed3b6c5bd435df304e7939f2de0362346cb9466d"


def linear_map(coordinates, values):
    xmean = sum(coordinates) / len(coordinates)
    ymean = sum(values) / len(values)
    slope = sum((x - xmean) * (y - ymean) for x, y in zip(coordinates, values)) / sum(
        (x - xmean) ** 2 for x in coordinates
    )
    return slope, ymean - slope * xmean


def extract(pdf_path: Path, output_dir: Path):
    checksum = hashlib.sha256(pdf_path.read_bytes()).hexdigest()
    if checksum != PDF_SHA256:
        raise ValueError(
            "Source PDF checksum differs; inspect geometry before reusing extraction"
        )
    with pdfplumber.open(pdf_path) as pdf:
        page = pdf.pages[3]
        # Labeled main-axis ticks, inspected on the complete rendered page.
        x_ticks = [page.lines[i]["x0"] for i in (92, 82, 72, 62)]
        y_ticks = [page.lines[i]["top"] for i in (37, 47)]
        pressure_map = linear_map(x_ticks, [0, 50, 100, 150])
        volume_map = linear_map(y_ticks, [25, 20])
        squares = [
            (i, r)
            for i, r in enumerate(page.rects)
            if r["fill"] and 215 < r["x0"] < 245 and 110 < r["top"] < 155
        ]
        circles = [
            (i, c)
            for i, c in enumerate(page.curves)
            if c["fill"]
            and 240 < c["x0"] < 360
            and c["top"] > 155
            and 2.9 < c["width"] < 3.0
        ]
        curves = [
            (i, c)
            for i, c in enumerate(page.curves)
            if len(c["pts"]) > 50 and c["x0"] < 240 and c["x1"] > 390
        ]
        if (len(squares), len(circles), len(curves)) != (9, 17, 1):
            raise ValueError("Unexpected Figure 2 object counts")
        observations = []
        for series, objects in (
            ("kuznetsov_squares", squares),
            ("dewaele_circles", circles),
        ):
            for index, obj in objects:
                observations.append(
                    (
                        series,
                        index,
                        (obj["x0"] + obj["x1"]) / 2,
                        (obj["top"] + obj["bottom"]) / 2,
                    )
                )
        for index, (x, y) in enumerate(curves[0][1]["pts"]):
            observations.append(("published_curve", index, x, y))
        output_dir.mkdir(parents=True, exist_ok=True)
        dataset = output_dir / "lead-dewaele-2019-figure2-vector.csv"
        with dataset.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(
                [
                    "series",
                    "pdf_object_or_vertex_index",
                    "pdf_x_pt",
                    "pdf_y_from_top_pt",
                    "plotted_pressure_gpa",
                    "atomic_volume_a3",
                ]
            )
            for series, index, x, y in observations:
                if series == "published_curve":
                    continue
                writer.writerow(
                    [
                        series,
                        index,
                        f"{x:.9f}",
                        f"{y:.9f}",
                        f"{pressure_map[0] * x + pressure_map[1]:.9f}",
                        f"{volume_map[0] * y + volume_map[1]:.9f}",
                    ]
                )
        curve_dataset = output_dir / "lead-dewaele-2019-figure2-curve.csv"
        with curve_dataset.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(
                [
                    "series",
                    "pdf_object_or_vertex_index",
                    "pdf_x_pt",
                    "pdf_y_from_top_pt",
                    "plotted_pressure_gpa",
                    "atomic_volume_a3",
                ]
            )
            for series, index, x, y in observations:
                if series == "published_curve":
                    writer.writerow(
                        [
                            series,
                            index,
                            f"{x:.9f}",
                            f"{y:.9f}",
                            f"{pressure_map[0] * x + pressure_map[1]:.9f}",
                            f"{volume_map[0] * y + volume_map[1]:.9f}",
                        ]
                    )
        metadata = {
            "source_doi": "10.3390/min9110684",
            "source_url": "https://www.mdpi.com/2075-163X/9/11/684",
            "source_pdf_sha256": checksum,
            "source_location": "Page 4, Figure 2 main panel",
            "method": "PDF vector filled-square/filled-circle bounding-box centers and fitted-curve vertices; duplicate stroke outlines and legend markers excluded",
            "pressure_ticks_pdf_x_pt": x_ticks,
            "pressure_tick_values_gpa": [0, 50, 100, 150],
            "volume_ticks_pdf_y_from_top_pt": y_ticks,
            "volume_tick_values_a3": [25, 20],
            "pressure_map_slope_intercept": pressure_map,
            "volume_map_slope_intercept": volume_map,
            "counts": {
                "kuznetsov_squares": len(squares),
                "dewaele_circles": len(circles),
                "published_curve": len(curves[0][1]["pts"]),
            },
            "dataset_sha256": hashlib.sha256(dataset.read_bytes()).hexdigest(),
            "curve_dataset_sha256": hashlib.sha256(
                curve_dataset.read_bytes()
            ).hexdigest(),
            "limitations": "Coordinates are quantized by the source plot, not exact observations. The nine-square selection is visibly attributed to Kuznetsov (legend says 2003; reference 18 is 2002). Curve vertices are a diagnostic and must never enter an observation fit. Pressure-scale identity is checked independently against Table 2 and Table 1.",
        }
        (output_dir / "lead-dewaele-2019-figure2-vector-source.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    extract(args.pdf, args.output_dir)
