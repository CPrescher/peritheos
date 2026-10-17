"""Recover Au calibration inputs without rewriting material documents.

The XLSX reader preserves the literal OOXML numeric tokens, cached formula
values, cell coordinates and blank errors. PDF extraction needs pdfplumber;
the numerical replay in reproduce_sakai_2025_gold.py has no PDF dependency.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
WORKBOOK_HASH = "c559b29fea9c793831e4bc6deef4d62712b5e78821fc420969f709cbd1060ff9"
DEWAELE_HASH = "38e8a84c98d7b72a12f092873b6e11e7a338fb05bd46f0a08eff5886da160515"
BASE_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs43246-025-00792-5/MediaObjects/"
)
WORKBOOK_URL = BASE_URL + "43246_2025_792_MOESM5_ESM.xlsx"
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
PAIRS = {
    1: ("Cu", "fcc", "Re", "hcp", "upstream"),
    2: ("Pt", "fcc", "Re", "hcp", "upstream"),
    3: ("Cu", "fcc", "W", "bcc", "upstream"),
    4: ("Au", "fcc", "W", "bcc", "direct_au_input"),
    7: ("Pt", "fcc", "MgO", "B1", "upstream"),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_workbook(path):
    """Read all source cells, including literal tokens and cached formulas."""
    if sha256(path) != WORKBOOK_HASH:
        raise ValueError("Unexpected Supplementary Data 1 workbook checksum")
    result = {}
    with ZipFile(path) as archive:
        strings = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = ["".join(node.itertext()) for node in strings.findall("s:si", NS)]
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        for sheet in book.findall("s:sheets/s:sheet", NS):
            rid = sheet.attrib[
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            ]
            target = targets[rid]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            content = archive.read(target)
            xml = ET.fromstring(content)
            cells = {}
            for cell in xml.findall("s:sheetData/s:row/s:c", NS):
                token = cell.findtext("s:v", default=None, namespaces=NS)
                formula = cell.findtext("s:f", default=None, namespaces=NS)
                kind = cell.attrib.get("t", "n")
                if kind == "s":
                    value = shared[int(token)] if token is not None else None
                elif kind == "inlineStr":
                    value = "".join(cell.find("s:is", NS).itertext())
                elif token is None:
                    value = None
                elif kind == "n":
                    value = float(token)
                else:
                    value = token
                cells[cell.attrib["r"]] = {
                    "value": value,
                    "source_token": token,
                    "type": kind,
                    "formula": formula,
                }
            result[sheet.attrib["name"]] = {
                "worksheet_xml_sha256": hashlib.sha256(content).hexdigest(),
                "cells": cells,
            }
    return result


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def paired_rows(sheets):
    """Normalize every row on the five observed Au ancestor/direct tables."""
    rows = []
    for number, (first, phase1, second, phase2, relevance) in PAIRS.items():
        name = next(n for n in sheets if n.startswith(f"Table S{number},"))
        cells = sheets[name]["cells"]
        indices = sorted({int(re.sub("[A-Z]+", "", c)) for c in cells})
        run = f"{number}"
        apparatus = "toroidal" if number != 7 else "conventional"
        for index in indices:
            if index < 5:
                continue
            marker = cells.get(f"A{index}", {}).get("value")
            if marker and "#" in str(marker):
                run = str(marker).replace("Run", "").replace("#", "").strip()
            if marker and str(marker).lower() in ("toroidal", "conventional"):
                apparatus = str(marker).lower()
            v1 = cells.get(f"B{index}", {}).get("value")
            v2 = cells.get(f"D{index}", {}).get("value")
            if not isinstance(v1, (int, float)) and not isinstance(v2, (int, float)):
                continue
            row = {
                "source_table": f"S{number}",
                "source_sheet": name,
                "source_row": index,
                "run": run,
                "apparatus": apparatus,
                "au_relevance": relevance,
                "source_marker": marker,
                "material_1": first,
                "phase_1": phase1,
                "material_2": second,
                "phase_2": phase2,
            }
            quantities = {
                "B": "volume_1_a3",
                "C": "volume_1_error_a3",
                "D": "volume_2_a3",
                "E": "volume_2_error_a3",
                "F": "secondary_a_angstrom",
                "G": "secondary_a_error_angstrom",
                "H": "secondary_c_angstrom",
                "I": "secondary_c_error_angstrom",
                "J": "secondary_c_over_a",
                "K": "secondary_c_over_a_error",
            }
            # S7's F/G are stress calculations, not lattice parameters.
            for col, quantity in quantities.items():
                cell = (
                    cells.get(f"{col}{index}", {}) if number != 7 or col < "F" else {}
                )
                row[quantity] = cell.get("value")
                row[quantity + "_source_token"] = cell.get("source_token")
            for col, quantity in (
                ("F", "pt_uniaxial_stress_gpa"),
                ("G", "mgo_uniaxial_stress_gpa"),
            ):
                cell = cells.get(f"{col}{index}", {}) if number == 7 else {}
                row[quantity] = cell.get("value")
                row[quantity + "_source_token"] = cell.get("source_token")
            rows.append(row)
    return rows


def dewaele_rows(path):
    """Extract Tables 4/5 by page and column position, preserving every token."""
    import pdfplumber

    if sha256(path) != DEWAELE_HASH:
        raise ValueError("Unexpected Dewaele (2018) supplementary PDF checksum")
    rows = []
    with pdfplumber.open(path) as pdf:
        configurations = [
            (8, 4, 1, [55, 73, 92, 127, 161, 184, 232, 278], 155, 627),
            (8, 4, 2, [281, 300, 318, 353, 387, 410, 458, 506], 155, 627),
            (9, 5, 3, [55, 87, 120, 171, 219, 259, 323, 391, 456, 506], 145, 366),
        ]
        for page, table, run, bounds, low, high in configurations:
            grouped = {}
            for word in pdf.pages[page].extract_words():
                if low < word["top"] < high and bounds[0] <= word["x0"] < bounds[-1]:
                    grouped.setdefault(round(word["top"], 1), []).append(word)
            for order, (_, words) in enumerate(sorted(grouped.items()), 1):
                tokens = []
                for left, right in zip(bounds[:-1], bounds[1:]):
                    tokens.append(
                        " ".join(w["text"] for w in words if left <= w["x0"] < right)
                    )
                if not tokens[1] or not re.fullmatch(r"\d+", tokens[1]):
                    raise ValueError(f"Unexpected PDF row tokens: {tokens}")
                fields = (
                    [
                        "membrane_pressure_bar",
                        "re_a_angstrom",
                        "re_c_angstrom",
                        "reported_re_pressure_gpa",
                        "au_a111_angstrom",
                        "au_a200_angstrom",
                    ]
                    if table == 4
                    else [
                        "membrane_pressure_bar",
                        "re_a_angstrom",
                        "re_c_angstrom",
                        "reported_re_pressure_gpa",
                        "au_a111_angstrom",
                        "al_fcc_a111_angstrom",
                        "al_hcp_a_angstrom",
                        "al_hcp_c_angstrom",
                    ]
                )
                values = dict(zip(fields, tokens[1:]))
                row = {
                    "source_table": table,
                    "source_page": page + 1,
                    "source_row_in_run": order,
                    "run": run,
                    "phase_au": "fcc",
                    "phase_re": "hcp",
                    "source_tokens_json": json.dumps(values, separators=(",", ":")),
                }
                for name in [
                    "membrane_pressure_bar",
                    "re_a_angstrom",
                    "re_c_angstrom",
                    "reported_re_pressure_gpa",
                    "au_a111_angstrom",
                    "au_a200_angstrom",
                    "al_fcc_a111_angstrom",
                    "al_hcp_a_angstrom",
                    "al_hcp_c_angstrom",
                ]:
                    row[name] = float(values[name]) if values.get(name) else None
                row["au_volume_a3"] = row["au_a111_angstrom"] ** 3
                row["re_volume_a3"] = (
                    math.sqrt(3) * row["re_a_angstrom"] ** 2 * row["re_c_angstrom"] / 2
                )
                row["volume_error_a3"] = None
                rows.append(row)
    counts = {run: sum(r["run"] == run for r in rows) for run in (1, 2, 3)}
    if counts != {1: 36, 2: 31, 3: 17}:
        raise ValueError(f"Unexpected recovered Dewaele run counts: {counts}")
    return rows


def recover(workbook, dewaele_pdf, output_dir):
    """Write source-specific normalized artifacts, never global metadata."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sheets = read_workbook(workbook)
    selected = {
        name: data
        for name, data in sheets.items()
        if name.startswith(tuple(f"Table S{i}," for i in range(1, 8)))
        or name in ("Table S8", "Table S14")
    }
    raw = {
        "workbook_sha256": WORKBOOK_HASH,
        "source_url": WORKBOOK_URL,
        "sheets": selected,
    }
    output_dir / "sakai-2025-au-workbook-cells.json".write_bytes(
        (json.dumps(raw, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    )
    rows = paired_rows(sheets)
    write_csv(output_dir / "sakai-2025-au-calibration-paired-volumes.csv", rows)
    write_csv(
        output_dir / "gold-sakai-2025-table-s4-au-w-volumes.csv",
        [row for row in rows if row["source_table"] == "S4"],
    )
    cells = sheets["Table S8"]["cells"]
    lineage = []
    for index in range(3, 48):
        if not cells.get(f"C{index}", {}).get("value"):
            continue
        lineage.append(
            {
                "source_row": index,
                **{
                    key: cells.get(f"{col}{index}", {}).get("value")
                    for col, key in zip(
                        "ABCDE",
                        [
                            "reference",
                            "reference_number",
                            "material",
                            "pressure_calibrant",
                            "note",
                        ],
                    )
                },
            }
        )
    write_csv(output_dir / "sakai-2025-table-s8-calibration-lineage.csv", lineage)
    anchor = []
    cells = sheets["Table S14"]["cells"]
    for index in range(6, 48):
        pressure = cells.get(f"F{index}", {}).get("value")
        if not isinstance(pressure, (int, float)):
            continue
        anchor.append(
            {
                "source_row": index,
                "volume_a3": cells[f"E{index}"]["value"],
                "fratanduono_2020_pressure_gpa": pressure,
                "sakai_rs_pressure_gpa": cells[f"G{index}"]["value"],
                "kind": "calculated_scale_conversion_input",
            }
        )
    write_csv(output_dir / "sakai-2025-table-s14-cu-anchor-calculated.csv", anchor)
    write_csv(
        output_dir / "gold-dewaele-2018-tables-4-5-au-re.csv", dewaele_rows(dewaele_pdf)
    )
    names = [
        "sakai-2025-au-calibration-paired-volumes.csv",
        "gold-sakai-2025-table-s4-au-w-volumes.csv",
        "sakai-2025-table-s8-calibration-lineage.csv",
        "sakai-2025-table-s14-cu-anchor-calculated.csv",
        "gold-dewaele-2018-tables-4-5-au-re.csv",
        "sakai-2025-au-workbook-cells.json",
    ]
    return {name: sha256(output_dir / name) for name in names}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--dewaele-pdf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DATA)
    args = parser.parse_args()
    print(
        json.dumps(recover(args.workbook, args.dewaele_pdf, args.output_dir), indent=2)
    )
