"""Recover plotted grossular observations from the final Angel (2022) PDF.

Repository-only; requires optional PyMuPDF. Supply the privately obtained
publication with --pdf. This extracts numerical marker centers, not figures.
Digitization bounds are resolution estimates, never experimental errors.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "docs/data/garnet-angel-2022/grossular"
PDF_SHA256 = "eebef304cc8ef328f3d9219d9006f94791e4680e42d92655a396890eac81b692"


def calibration(coordinates, values):
    """Least-squares affine axis, with the largest tick residual retained."""
    xm, ym = (sum(a) / len(a) for a in (coordinates, values))
    slope = sum((x - xm) * (y - ym) for x, y in zip(coordinates, values)) / sum(
        (x - xm) ** 2 for x in coordinates
    )
    intercept = ym - slope * xm
    error = max(abs(slope * x + intercept - y) for x, y in zip(coordinates, values))
    return slope, intercept, error


def recover(pdf, output=DATA):
    import fitz  # Optional maintenance dependency, not part of the library.

    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError("This extraction is qualified only for the recorded final PDF")
    with fitz.open(pdf) as document:
        drawings = document[7].get_drawings()
    tx = [item[1].x for item in drawings[333]["items"][:6]]
    vy = [item[1].y for item in drawings[337]["items"][:9]]
    tcal = calibration(tx, [1500, 1200, 900, 600, 300, 0])
    vcal = calibration(vy, [129 - 0.5 * i for i in range(9)])
    if tcal[2] > 0.5 or vcal[2] > 0.0012:
        raise ValueError("Axis extraction exceeds declared tick residuals")
    groups = [
        ("isaak", range(419, 434), (1, 0, 0), True),
        ("thieblot_selected", range(453, 460), (1, 0, 1), True),
        ("thieblot_excluded_visible", range(518, 524), (1, 0, 1), False),
    ]
    output.mkdir(parents=True, exist_ok=True)
    receipt = {
        "recovered_on": "2026-11-12",
        "source": "Angel et al. (2022), final publication Figure 4, PDF page 8",
        "doi": "10.1007/s00410-022-01918-5",
        "pdf_sha256": PDF_SHA256,
        "pymupdf_version": fitz.VersionBind,
        "data_kind": "digitized_processed_plot_observations",
        "method": "Vector marker bounding-box centers; affine calibration to six temperature and nine volume major ticks. Filled/open markers checked against the rendered figure and legend.",
        "temperature_axis": dict(
            zip(["slope", "intercept", "max_tick_residual_k"], tcal)
        ),
        "volume_axis": dict(
            zip(["slope", "intercept", "max_tick_residual_cm3_mol"], vcal)
        ),
        "conservative_digitization_bounds": {
            "temperature_k": 2,
            "molar_volume_cm3_mol": 0.003,
            "meaning": "Resolution bounds, not experimental standard deviations or confidence intervals; axes, calibration and normalization induce shared errors.",
        },
        "qualification": "Thiéblot points are Angel's already processed plotted series, not recovered original Table values. Seven filled markers correspond to S1's selected subset; the article's exclusion prose conflicts with S1. Six open markers are retained separately as excluded visible evidence. Two of the 15 original observations are not identified here. No original uncertainties are recovered.",
        "groups": {},
    }
    for name, indices, color, filled in groups:
        rows = []
        for i in indices:
            drawing = drawings[i]
            rect = drawing["rect"]
            if (
                tuple(drawing["color"]) != color
                or tuple(drawing["fill"]) != (color if filled else (1, 1, 1))
                or len(drawing["items"]) != 4
                or not 3 < rect.width < 3.4
                or not 3 < rect.height < 3.4
            ):
                raise ValueError(f"Unexpected marker at path {i}")
            x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            rows.append(
                dict(
                    drawing_index=i,
                    x_pdf_point=x,
                    y_pdf_point=y,
                    temperature_k=tcal[0] * x + tcal[1],
                    molar_volume_cm3_mol=vcal[0] * y + vcal[1],
                )
            )
        filename = f"figure4-{name.replace('_', '-')}.csv"
        path = output / filename
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        receipt["groups"][name] = {
            "rows": len(rows),
            "transcription": filename,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    (output / "figure4-recovery.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA)
    args = parser.parse_args()
    recover(args.pdf, args.output)
