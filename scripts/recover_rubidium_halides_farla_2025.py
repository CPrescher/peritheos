#!/usr/bin/env python3
"""Normalize the official Farla 2025 supplement without repairing source tokens.

Usage: uv run python scripts/recover_rubidium_halides_farla_2025.py SOURCE.docx
The source is AIP Figshare file 52698776, doi:10.60893/figshare.jap.28515917.v1.
Only the standard library is required. No inferred observations are generated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
PREFIX = "rubidium-halides-farla-2025"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
COMPOUNDS = ("rbcl", "rbbr", "rbi")


def extract_tables(source: Path) -> list[list[list[str]]]:
    """Keep literal OOXML text, including malformed parentheses and empty cells."""
    with ZipFile(source) as archive:
        document = ET.fromstring(archive.read("word/document.xml"))
    return [
        [
            [
                "".join(t.text or "" for t in cell.findall(".//w:t", NS))
                for cell in row.findall("w:tc", NS)
            ]
            for row in table.findall("w:tr", NS)
        ]
        for table in document.findall(".//w:body/w:tbl", NS)
    ]


def number(token: str) -> tuple[str, str]:
    """Expand last-digit errors; retain a printed (0) as zero, never invent sigma.

    S1's '1.19(0' lacks a closing parenthesis; its literal token remains archived.
    The parsed uncertainty is the printed zero and its precision is unresolved.
    """
    if not token:
        return "", ""
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\((\d+)\)?)?", token)
    if match is None:
        raise ValueError(f"Unrecognized source token {token!r}")
    value, error = match.groups()
    if error is None:
        return value, ""
    places = len(value.partition(".")[2])
    return value, f"{int(error) * 10**-places:.{places}f}"


def observation_rows(tables: list[list[list[str]]], compound: str) -> list[dict]:
    """Retain all DAC/LVP rows, including blank sample or calibrant observations."""
    rows = []
    for index, cells in enumerate(tables[0], 1):
        if cells[0] != "300":
            continue
        pcol, vcol = {"rbcl": (1, 2), "rbi": (1, 3), "rbbr": (4, 5)}[compound]
        if not cells[vcol]:
            continue
        row = {
            "source_table": "S1",
            "source_row": index,
            "run": "DAC_RbBr" if compound == "rbbr" else "DAC_RbCl_RbI",
            "temperature_k": "300",
            "temperature_uncertainty_k": "",
            "pressure_gpa_source_token": cells[pcol],
            "volume_a3_source_token": cells[vcol],
        }
        row["pressure_gpa"], row["pressure_uncertainty_gpa"] = number(cells[pcol])
        row["volume_a3"], row["volume_uncertainty_a3"] = number(cells[vcol])
        rows.append(row)
    table_index = COMPOUNDS.index(compound) + 1
    run = ""
    fields = [
        "cscl_volume_a3",
        "cscl_pressure_gpa",
        "pt_volume_a3",
        "pt_pressure_gpa",
        "mo_volume_a3",
        "mo_pressure_gpa",
        "volume_a3",
        "density_g_cm3",
    ]
    for index, cells in enumerate(tables[table_index], 1):
        if cells[0].startswith("Run #"):
            run = re.search(r"BT\d+", cells[0])[0]
        if not cells[0].isdigit():
            continue
        row = {
            "source_table": f"S{table_index + 1}",
            "source_row": index,
            "run": run,
            "temperature_k": cells[0],
            # The 30 K is an author-assigned uncertainty, not an observed error.
            "temperature_uncertainty_k": "30",
        }
        for field, token in zip(fields, cells[1:]):
            row[field], error = number(token)
            row[field + "_source_token"] = token
            if field != "density_g_cm3":
                stem, unit = field.rsplit("_", 1)
                row[f"{stem}_uncertainty_{unit}"] = error
        row["pressure_gpa"] = row["cscl_pressure_gpa"]
        row["pressure_uncertainty_gpa"] = row["cscl_pressure_uncertainty_gpa"]
        row["pressure_gpa_source_token"] = row["cscl_pressure_gpa_source_token"]
        rows.append(row)
    return rows


def benchmark_rows(tables: list[list[list[str]]], compound: str) -> list[dict]:
    """Calculated tables are benchmarks, never experimental observations."""
    table_index = COMPOUNDS.index(compound) + 4
    result = []
    for index, cells in enumerate(tables[table_index], 1):
        if not cells[0].isdigit():
            continue
        for temperature, token in zip((300, 500, 1000, 1500, 2000, 2500), cells[1:]):
            result.append(
                {
                    "source_table": f"S{table_index + 1}",
                    "source_row": index,
                    "pressure_gpa": cells[0],
                    "temperature_k": temperature,
                    "volume_ratio": token.removesuffix("a"),
                    "volume_ratio_source_token": token,
                    "molten_source_flag": int(token.endswith("a")),
                }
            )
    return result


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write deterministic union columns, with missing observations left blank."""
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DATA)
    args = parser.parse_args()
    tables = extract_tables(args.source)
    if len(tables) != 7:
        raise ValueError("Expected the official seven-table supplement")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = {
        "source_file": args.source.name,
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "source_doi": "10.60893/figshare.jap.28515917.v1",
        "source_url": "https://ndownloader.figshare.com/files/52698776",
        "source_license": "CC BY 4.0 (AIP Figshare metadata, version 1)",
        "tables": tables,
    }
    (args.output_dir / f"{PREFIX}-source-tables.json").write_text(
        json.dumps(raw, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    for compound in COMPOUNDS:
        write_csv(
            args.output_dir / f"{PREFIX}-{compound}-observations.csv",
            observation_rows(tables, compound),
        )
        write_csv(
            args.output_dir / f"{PREFIX}-{compound}-benchmarks.csv",
            benchmark_rows(tables, compound),
        )


if __name__ == "__main__":
    main()
