#!/usr/bin/env python3
"""Recover all printed Table 1 rows and match them to the 2022 deposit.

Only numerical table content is packaged. The publisher PDF is not copied.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import tempfile
from pathlib import Path

from scripts.forsterite_recover_katsura_2022 import DATA

PREFIX = "forsterite-katsura-2009"
TOKEN = r"\d+\.\d+\(\d+\)"
ROW = re.compile(rf"^\s*({TOKEN})\s+(\d+)" + rf"\s+({TOKEN})" * 4)


def unpack(token: str) -> tuple[str, str]:
    value, digits = token[:-1].split("(")
    decimals = len(value.split(".")[1])
    return value, f"{int(digits) * 10 ** (-decimals):.{decimals}f}"


def recover(pdf: Path) -> tuple[list[dict], dict]:
    with tempfile.TemporaryDirectory(prefix="forsterite-2009-") as temp:
        text = Path(temp) / "article.txt"
        subprocess.run(["pdftotext", "-layout", str(pdf), str(text)], check=True)
        source = text.read_text()
    section = source.split(
        "Pressure, temperature, lattice parameter, and unit cell volume of forsterite"
    )[1].split("bulk modulus, Ks")[0]
    rows = []
    run = None
    names = ["pressure_gpa", "a_angstrom", "b_angstrom", "c_angstrom", "volume_a3"]
    for line in section.splitlines():
        if line.strip().startswith("Run M"):
            run = line.strip().split()[1]
        match = ROW.match(line)
        if match is None:
            continue
        tokens = match.groups()
        row = dict(
            table_row=str(len(rows) + 1),
            run_label=run,
            workbook_run_id="",
            temperature_k=tokens[1],
        )
        for name, token in zip(names, [tokens[0], *tokens[2:]]):
            value, error = unpack(token)
            row[name] = value
            row[name + "_reported_uncertainty"] = error
            row[name + "_source_token"] = token
        rows.append(row)
    assert len(rows) == 52
    with (DATA / "forsterite-katsura-2022-olivine-pvt.csv").open(newline="") as stream:
        workbook = list(csv.DictReader(stream))
    unlisted = []
    for observation in workbook:
        matches = [
            row
            for row in rows
            if row["run_label"] == observation["run_id"][:4]
            and float(row["temperature_k"])
            == float(observation["temperature_uncorrected_k"])
            and abs(
                float(row["pressure_gpa"]) - float(observation["pressure_matsui_gpa"])
            )
            < 0.02
        ]
        if not matches:
            unlisted.append(
                dict(
                    run_id=observation["run_id"],
                    workbook_row=int(observation["source_row"]),
                )
            )
            continue
        assert len(matches) == 1 and not matches[0]["workbook_run_id"]
        matches[0]["workbook_run_id"] = observation["run_id"]
    metadata = {
        "source_doi": "10.1016/j.pepi.2008.08.002",
        "source_version": "User-supplied final published PDF, Table 1, journal p.89 / PDF p.4",
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "row_count": len(rows),
        "workbook_row_count": len(workbook),
        "matched_workbook_rows": sum(bool(row["workbook_run_id"]) for row in rows),
        "workbook_rows_absent_from_printed_table": unlisted,
        "table_only_rows": [
            row["table_row"] for row in rows if not row["workbook_run_id"]
        ],
        "matching": "Same run and uncorrected T; original Matsui pressure within 0.02 GPa, accommodating printed/deposited rounding differences. Unique matches required.",
        "uncertainty": "Printed parenthetical errors retained literally. Section 3.1 describes artificially inflated profile errors; no confidence level or temperature sigma is inferred.",
        "normalization": "2009 Section 3.2 adopts V0,300=290.79 A^3 for its analysis; this differs from the 2022 geotherm script's V0=289.58 A^3. Raw printed volumes remain separate from deposited dimensionless ratios.",
        "disposition": "Absence of M375006 from Table 1 supports a separate table-matched sensitivity fit, not silent removal from the author-deposited 51-row input.",
    }
    return rows, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    args = parser.parse_args()
    rows, metadata = recover(args.pdf)
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    (DATA / f"{PREFIX}-table1.csv").write_text(out.getvalue())
    (DATA / f"{PREFIX}-source.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
