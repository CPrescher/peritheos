"""Recover Figure 10 marker centers from the supplied publisher PDF vectors.

Requires pdfplumber only for re-extraction; numerical refits use the bundled CSV.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/rhenium-sakai-2018-figure10.csv"
PROVENANCE = DATA.with_name("rhenium-sakai-2018-figure10-source.json")


def extract(pdf: Path) -> None:
    import pdfplumber

    record = next(
        r
        for r in json.loads(
            (ROOT / "peritheos/data/materials/rhenium.eosmat").read_text()
        )["eos_records"]
        if r["identifier"] == "rhenium_sakai_2018_yokoo_pt_vinet"
    )
    checksum = hashlib.sha256(pdf.read_bytes()).hexdigest()
    assert checksum == record["scientific_validation"]["primary_source_check"]["sha256"]
    with pdfplumber.open(pdf) as document:
        page = document.pages[10]
        # Axes identified from their labeled major ticks in Figure 10.
        xs = np.array([page.lines[i]["x0"] for i in range(39, 48)])
        ys = np.array([page.lines[i]["top"] for i in range(14, 21)])
        pressures = np.arange(0, 401, 50)
        volumes = np.arange(18, 31, 2)
        pressure_map = np.polyfit(xs, pressures, 1)
        volume_map = np.polyfit(ys, volumes, 1)
        rows = []
        groups = (
            ("rp01_dewaele_pt", range(3, 55, 2), True, (0.137, 0.122, 0.125)),
            ("rp01_yokoo_pt", range(55, 107, 2), True, (0.698, 0.706, 0.714)),
            ("micro17_yokoo_pt", range(107, 113), False, (0.698, 0.706, 0.714)),
        )
        for series, indices, filled, color in groups:
            for index in indices:
                marker = page.curves[index]
                assert marker["fill"] == filled
                assert tuple(marker["non_stroking_color"]) == color
                assert 2.3 < marker["width"] < 3.3
                assert 2.3 < marker["height"] < 3.3
                x = (marker["x0"] + marker["x1"]) / 2
                y = (marker["top"] + marker["bottom"]) / 2
                rows.append(
                    [
                        series,
                        index,
                        x,
                        y,
                        np.polyval(pressure_map, x),
                        np.polyval(volume_map, y),
                    ]
                )
        assert len(rows) == 58
        with DATA.open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(
                [
                    "series",
                    "pdf_curve_index",
                    "pdf_x_pt",
                    "pdf_top_pt",
                    "pressure_gpa",
                    "volume_a3",
                ]
            )
            writer.writerows(rows)
        provenance = {
            "doi": "10.1080/08957959.2018.1448082",
            "source_filename": pdf.name,
            "source_sha256": checksum,
            "pdf_page_one_based": 11,
            "article_page": 10,
            "figure": "10",
            "coordinate_basis": "PDF points, x from left and top from top of page",
            "pressure_ticks": list(zip(xs.tolist(), pressures.tolist())),
            "volume_ticks": list(zip(ys.tolist(), volumes.tolist())),
            "pressure_affine_map": pressure_map.tolist(),
            "volume_affine_map": volume_map.tolist(),
            "maximum_pressure_tick_residual_gpa": float(
                np.max(abs(np.polyval(pressure_map, xs) - pressures))
            ),
            "maximum_volume_tick_residual_a3": float(
                np.max(abs(np.polyval(volume_map, ys) - volumes))
            ),
            "series_counts": {name: len(indices) for name, indices, _, _ in groups},
            "method": (
                "Bounding-box centers of vector marker paths; one fill path per "
                "RP01 marker, excluding the duplicate stroke path. Micro17 uses "
                "six larger open diamond paths. Legend markers and three literature "
                "curves are excluded. Overlapping distinct markers are retained. "
                "Affine axes fitted to all labeled major ticks. No source experimental "
                "uncertainties are available; vector precision is not measurement precision."
            ),
            "csv_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        }
        PROVENANCE.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    extract(parser.parse_args().pdf)
