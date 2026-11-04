#!/usr/bin/env python3
"""Extract the Olivine sheet without recalculating the deposited workbook.

Source: Katsura, Zenodo 5644426 (CC BY 4.0). Every occupied cell is retained,
including formulas and their author-saved results. The CSV preserves all 41
source columns, with Excel addresses preventing ambiguous repeated labels.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "docs/data/forsterite-katsura-2022-sources"
DATA = ROOT / "peritheos/data/datasets"
PREFIX = "forsterite-katsura-2022"
NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
# Coordinates disambiguate repeated V/V0 and pressure headers. F:H are MgO;
# AN:AO are the olivine sample. The combined pressure is the saved AD/AE pair.
COLUMN_MAP = {
    "run_id": "A",
    "temperature_uncorrected_k": "B",
    "pressure_matsui_gpa": "D",
    "pressure_matsui_sigma_gpa": "E",
    "mgo_volume_ratio": "F",
    "mgo_volume_ratio_lower": "G",
    "mgo_volume_ratio_upper": "H",
    "temperature_corrected_upper_k": "S",
    "temperature_corrected_lower_k": "T",
    "temperature_corrected_k": "U",
    "pressure_tange_vinet_gpa": "X",
    "pressure_tange_vinet_sigma_gpa": "Y",
    "pressure_tange_bm3_gpa": "AB",
    "pressure_tange_bm3_sigma_gpa": "AC",
    "pressure_mean_gpa": "AD",
    "pressure_spread_gpa": "AE",
    "forsterite_volume_ratio": "AN",
    "forsterite_volume_ratio_sigma": "AO",
}


def extract() -> tuple[str, dict]:
    """Return a lossless numeric CSV and cell-addressed source evidence."""
    path = SOURCES / f"{PREFIX}-source.xlsx"
    with ZipFile(path) as archive:
        strings = [
            "".join(si.itertext())
            for si in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.find("x:sheets", NS)
        assert sheets is not None and sheets[0].attrib["name"] == "Olivine"
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        cells = {}
        for cell in sheet.findall(".//x:sheetData/x:row/x:c", NS):
            value = cell.find("x:v", NS)
            if value is None:
                continue
            raw = value.text
            if cell.attrib.get("t") == "s":
                raw = strings[int(raw)]
            formula = cell.find("x:f", NS)
            cells[cell.attrib["r"]] = {
                "value": raw,
                "kind": "text" if cell.attrib.get("t") == "s" else "number",
                **({"formula": formula.text} if formula is not None else {}),
            }
    columns = [chr(i) for i in range(65, 91)] + ["A" + chr(i) for i in range(65, 80)]
    out = io.StringIO(newline="")
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["source_row", *[f"col_{col.lower()}" for col in columns]])
    for row in range(4, 55):
        writer.writerow(
            [row, *[cells.get(f"{col}{row}", {}).get("value", "") for col in columns]]
        )
    evidence = {
        "source_doi": "10.5281/zenodo.5644426",
        "worksheet": "Olivine",
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "data_rows": [4, 54],
        "source_columns": columns,
        "cells": cells,
        "notes": "Cached author results retained verbatim; no workbook recalculation.",
    }
    return out.getvalue(), evidence


def main() -> None:
    csv_text, evidence = extract()
    DATA / f"{PREFIX}-olivine-workbook.csv".write_bytes((csv_text).encode("utf-8"))
    DATA / f"{PREFIX}-olivine-cells.json".write_bytes(
        (json.dumps(evidence, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    )
    out = io.StringIO(newline="")
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["source_row", *COLUMN_MAP])
    for row in range(4, 55):
        writer.writerow(
            [
                row,
                *[
                    evidence["cells"][f"{col}{row}"]["value"]
                    for col in COLUMN_MAP.values()
                ],
            ]
        )
    DATA / f"{PREFIX}-olivine-pvt.csv".write_bytes((out.getvalue()).encode("utf-8"))
    # Candidate six-column export in the author's xlsread order. Missing sT
    # stays empty: using the correction bracket as total sT is not established.
    out = io.StringIO(newline="")
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(
        [
            "pressure_gpa",
            "pressure_sigma_gpa",
            "volume_ratio",
            "volume_ratio_sigma",
            "temperature_k",
            "temperature_sigma_k",
        ]
    )
    for row in range(4, 55):
        writer.writerow(
            [
                *[
                    evidence["cells"][f"{col}{row}"]["value"]
                    for col in ["AD", "AE", "AN", "AO", "U"]
                ],
                "",
            ]
        )
    DATA / f"{PREFIX}-candidate-fitting-input.csv".write_bytes(
        (out.getvalue()).encode("utf-8")
    )
    DATA / f"{PREFIX}-candidate-fitting-input.json".write_bytes(
        (
            json.dumps(
                {
                    "status": "partial_candidate_reconstruction_not_original_author_export",
                    "source_doi": "10.5281/zenodo.5644426",
                    "source_sheet": "Olivine",
                    "source_rows": [4, 54],
                    "columns_in_author_script_order": [
                        "AD",
                        "AE",
                        "AN",
                        "AO",
                        "U",
                        None,
                    ],
                    "pressure_choice": "Author-saved combined AD/AE. Exact separate export selection cannot be confirmed; individual Vinet and BM3 coordinates remain in the complete source transcription.",
                    "missing": ["temperature_sigma_k"],
                    "notes": "All 51 rows retained. S/T are corrected temperatures evaluated at pressure brackets; their half-span is recoverable but not established as the total sigmaT supplied to the bootstrap. Blank sigmaT must not be silently read as zero. No claim of original bootstrap fit-input parity.",
                },
                indent=2,
            )
            + "\n"
        ).encode("utf-8")
    )


if __name__ == "__main__":
    main()
