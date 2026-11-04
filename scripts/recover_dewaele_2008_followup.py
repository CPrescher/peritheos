"""Recover graph checkpoints for the scoped Dewaele 2008 follow-up audit.

Repository maintenance only; needs pdfplumber. Paths and symbol identities were
visually inspected against Figures 2 and 4 of the supplied publisher PDF.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

SOURCE = (
    Path(__file__).resolve().parents[1] / "docs/data/dewaele-2008-followup-source.json"
)
# Legend identities, not inferred from best numerical agreement.
CURVES = {
    "Al": 449,
    "Cu": 297,
    "Ta": 298,
    "Au": 300,
    "Ni": 366,
    "Ag": 330,
    "Pt": 299,
    "Mo": 352,
    "W": 263,
}


def recover(path):
    import pdfplumber

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != "a06fd5cb9433b99f35b20582be5214b1e805014c150a2907fdf8933776cf6daf":
        raise ValueError("This mapping requires the inspected publisher PDF")
    with pdfplumber.open(path) as pdf:
        page = pdf.pages[4]
        curves = {}
        for metal, index in CURVES.items():
            c = page.curves[index]
            curves[metal] = {
                "object_kind": "curve",
                "object_index": index,
                "pdf_points": c["pts"],
                "linewidth_pt": c["linewidth"],
                "stroking_color": c["stroking_color"],
            }
        lines = page.lines[718:755]
        assert all(c["stroking_color"] == (0.0, 0.52941, 0.79608) for c in lines)
        curves["Co"] = {
            "object_kind": "lines",
            "object_indices": list(range(718, 755)),
            "pdf_points": [
                [sum(x for x, _ in c["pts"]) / 2, c["pts"][0][1]] for c in lines
            ],
            "description": "midpoints of short horizontal dotted-curve segments",
        }
        for c in curves.values():
            c["pressure_difference_points_gpa"] = [
                [
                    (x - 362.154) / 164.720 * 160,
                    (215.718988 - y) / ((215.718988 - 168.424988) / 4),
                ]
                for x, y in c["pdf_points"]
            ]
        page = pdf.pages[8]
        fe = {}
        for name, i in [("thin_without_semicore", 24), ("thick_with_semicore", 25)]:
            c = page.curves[i]
            fe[name] = {
                "curve_index": i,
                "linewidth_pt": c["linewidth"],
                "pdf_points": c["pts"],
                "pressure_volume_points": [
                    [(x - 129.164) / 111.171 * 250, 12 - (y - 72.704988) / 98.498 * 5]
                    for x, y in c["pts"]
                ],
            }
    return {
        "source_pdf_sha256": digest,
        "source_doi": "10.1103/PhysRevB.78.104102",
        "observation_kind": "graphical model checkpoints; not measured or raw DFT data",
        "figure2b": {
            "page": 5,
            "x_axis": "2004 revised ruby pressure, GPa",
            "y_axis": "reduced-shock pressure minus 2004 ruby pressure, GPa",
            "axis_mapping": {
                "x0_pt": 362.154,
                "x160_pt": 526.874,
                "y0_pt": 215.718988,
                "y4_pt": 168.424988,
            },
            "reading_halfwidth_gpa": 0.03,
            "metals": curves,
        },
        "figure4_fe_lda": {
            "page": 9,
            "axis_mapping": {
                "x0_pt": 129.164,
                "x250_pt": 240.335,
                "y12_pt": 72.704988,
                "y7_pt": 171.202988,
            },
            "curves": fe,
            "all_electron_check": "Fe triangles in the inspected panel are filled GGA symbols; no independent LDA triangle series is identified.",
        },
    }


if __name__ == "__main__":
    SOURCE.write_text(json.dumps(recover(Path(sys.argv[1])), indent=2) + "\n")
    print(SOURCE)
