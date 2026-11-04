#!/usr/bin/env python3
"""Normalize factual Fe-ringwoodite tables and check original PDF transcriptions.

No network access or source-file redistribution is required. The optional
--source-dir checks Table S1, Liu Table 1 and optional SLB Table A1 with Poppler.
Original observations, parenthesized errors, and the inconsistent printed
"density" column remain intact in the raw CSV.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
PREFIX = "ferrous-ringwoodite-"
AVOGADRO = 6.02214076e23
Z = 8
RAW = DATA / f"{PREFIX}armentrout-2011-table-s1-raw.csv"
NORMALIZED = DATA / f"{PREFIX}armentrout-2011-table-s1.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def normalized_csv() -> str:
    """Derive volumes and declared temperature hypotheses, retaining all 18 rows."""
    rows = []
    for raw in read_csv(RAW):
        nacl_a = float(raw["nacl_b1_a_angstrom"])
        sample_a = float(raw["ringwoodite_a_angstrom"])
        nacl_sd = int(raw["nacl_b1_a_parentheses"]) * 1e-4
        sample_sd = int(raw["ringwoodite_a_parentheses"]) * 1e-4
        upstream = float(raw["upstream_temperature_k"])
        downstream = float(raw["downstream_temperature_k"])
        mean_t = (upstream + downstream) / 2
        rows.append(
            {
                "measurement": raw["measurement"],
                "nacl_b1_a_angstrom": f"{nacl_a:.4f}",
                "nacl_b1_a_sd_angstrom": f"{nacl_sd:.4f}",
                "ringwoodite_a_angstrom": f"{sample_a:.4f}",
                "ringwoodite_a_sd_angstrom": f"{sample_sd:.4f}",
                "upstream_temperature_k": raw["upstream_temperature_k"],
                "downstream_temperature_k": raw["downstream_temperature_k"],
                "pressure_gpa": raw["pressure_gpa"],
                "pressure_reported_uncertainty_gpa": (
                    f"{int(raw['pressure_parentheses']) * 0.01:.2f}"
                ),
                "exclusion_footnote": raw["exclusion_footnote"],
                "included_by_table_s1": str(int(raw["exclusion_footnote"] == "0")),
                "derived_sample_volume_a3": f"{sample_a**3:.10f}",
                "derived_sample_volume_sd_a3": f"{3 * sample_a**2 * sample_sd:.10f}",
                "derived_sample_molar_volume_cm3_mol": (
                    f"{sample_a**3 * AVOGADRO * 1e-24 / Z:.10f}"
                ),
                "derived_nacl_b1_volume_a3": f"{nacl_a**3:.10f}",
                "derived_nacl_b1_volume_sd_a3": f"{3 * nacl_a**2 * nacl_sd:.10f}",
                "assumed_mean_sample_temperature_k": f"{mean_t:.2f}",
                "assumed_nacl_temperature_k": f"{(mean_t + 300) / 2:.2f}",
                "source_density_column": raw["source_density_column"],
                "source_density_parentheses": raw["source_density_parentheses"],
            }
        )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def check_source_pdfs(source_dir: Path) -> dict:
    """Compare original numeric tables and optional SLB italic/systematics flags."""
    table = source_dir / "armentrout-2011-table-s1.pdf"
    text = subprocess.check_output(["pdftotext", "-raw", str(table), "-"], text=True)
    segments = re.findall(
        r"Fay\d{3}.*?(?=Fay\d{3}|1\) Value excluded)", text, flags=re.S
    )
    pattern = re.compile(
        r"(Fay\d{3}(?:heat|quench|preheat))\s+"
        r"(\d\.\d{4})\s*\((\d+)\)\s*"
        r"(\d\.\d{4})([12]?)\s*\((\d+)\)\s*"
        r"(\d+)\s+(\d+)\s+(\d+\.\d+)\s*\((\d+)\)\s*"
        r"(\d+\.\d+)\s*\((\d+)\)"
    )
    extracted = []
    for segment in segments:
        segment = re.sub(r"quenc\s+h", "quench", segment)
        segment = re.sub(r"prehe\s+at", "preheat", segment)
        match = pattern.fullmatch(" ".join(segment.split()))
        if match is None:
            raise ValueError(f"Unparsed Table S1 row: {segment!r}")
        name, na, ne, fe, flag, ee, tu, td, density, de, p, pe = match.groups()
        if name == "Fay008heat" and tu == "3003":
            tu = "300"  # Printed superscript footnote 3, visually checked.
        extracted.append(
            [name, na, ne, fe, ee, tu, td, density, de, p, pe, flag or "0"]
        )
    expected = [list(row.values()) for row in read_csv(RAW)]
    if extracted != expected:
        raise ValueError("Official Table S1 does not match the raw CSV")

    liu = source_dir / "liu-2008a.pdf"
    text = subprocess.check_output(["pdftotext", "-layout", str(liu), "-"], text=True)
    extracted_liu = []
    for line in text.splitlines():
        values = line.split()
        if (
            len(values) == 10
            and re.fullmatch(r"\d\.\d\d", values[0])
            and values[1].isdigit()
        ):
            extracted_liu.append(values)
    expected_liu = [
        list(row.values()) for row in read_csv(DATA / f"{PREFIX}liu-2008-table1.csv")
    ]
    if extracted_liu != expected_liu:
        raise ValueError("Liu Table 1 does not match the CSV")
    result = {
        "table_s1_all_fields_match": True,
        "table_s1_rows": len(extracted),
        "liu_table1_all_cells_match": True,
        "liu_table1_rows": len(extracted_liu),
        "source_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (table, liu)
        },
        "finger_and_suzuki_tables": "Manually transcribed and visually checked; see source manifest.",
    }
    upstream = source_dir / "stixrude-2005.pdf"
    if upstream.is_file():
        # Table A1 is rotated. Poppler's XML preserves the italic/systematics
        # marks lost by ordinary text extraction; locate the Fe row itself.
        xml = subprocess.check_output(
            ["pdftohtml", "-f", "22", "-l", "22", "-xml", "-stdout", str(upstream)],
            text=True,
        )
        nodes = list(ET.fromstring(xml).iter("text"))
        # The initial capital F is a separate PDF text object.
        row_node = next(n for n in nodes if "e-ringw" in "".join(n.itertext()))
        numeric = [
            n
            for n in nodes
            if n.get("left") == row_node.get("left")
            and re.fullmatch(r"[0-9.]+|\([0-9]+\)", "".join(n.itertext()))
        ]
        expected = []
        flags = []
        for row in read_csv(DATA / f"{PREFIX}stixrude-2005-table-a1.csv"):
            expected.append(row["value"])
            flags.append(row["from_systematics"] == "1")
            if row["reported_error"]:
                decimals = len(row["value"].partition(".")[2])
                digits = round(float(row["reported_error"]) * 10**decimals)
                expected.append(f"({digits})")
        if ["".join(n.itertext()) for n in numeric] != expected:
            raise ValueError("SLB Table A1 Fe-row values/errors do not match CSV")
        if [
            n.find("i") is not None
            for n in numeric
            if not "".join(n.itertext()).startswith("(")
        ] != flags:
            raise ValueError("SLB Table A1 systematics flags do not match PDF italics")
        result["stixrude_table_a1_all_values_errors_and_systematics_match"] = True
        result["stixrude_table_a1_parameters"] = len(flags)
        result["source_sha256"][upstream.name] = hashlib.sha256(
            upstream.read_bytes()
        ).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--source-dir", type=Path)
    args = parser.parse_args()
    normalized = normalized_csv()
    if args.check:
        if NORMALIZED.read_text(encoding="utf-8") != normalized:
            raise SystemExit("Normalized Fe-ringwoodite observations are stale")
    else:
        NORMALIZED.write_text(normalized, encoding="utf-8")
    if args.source_dir:
        print(json.dumps(check_source_pdfs(args.source_dir), indent=2))


if __name__ == "__main__":
    main()
