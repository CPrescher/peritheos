"""Export a reusable corrected P-V-T view without changing the source table.

Only row 15 pressure has an evidence-supported replacement: 4.58 GPa from
Figure 2, with a separate +/-0.05 GPa graphical reading bound. Original
experimental pressure errors are unreported and remain blank. The inconsistent
source lattice parameters are retained in the original transcription only.
"""

from __future__ import annotations

import argparse
import csv
import io
import json

from scripts.diopside_reproduce import DATA, ROOT, load_rows

OUTPUT = DATA / "diopside-zhao-1998-table1-corrected-pvt.csv"
IDENTIFIER = "diopside_zhao_1998_table1_corrected_pvt"
COLUMNS = [
    "source_row",
    "pressure_gpa",
    "pressure_error_gpa",
    "temperature_k",
    "temperature_error_k",
    "cell_volume_a3",
    "cell_volume_error_a3",
    "source_pressure_gpa",
    "pressure_correction_applied",
    "pressure_graphical_reading_bound_gpa",
    "pressure_source",
    "duplicate_of_source_row",
    "quality_notes",
    "correction_note",
]


def corrected_rows(digitized_pressure_gpa):
    output = []
    for original in load_rows("diopside-zhao-1998-table1"):
        row = {name: original[name] for name in COLUMNS[:7]}
        number = int(original["source_row"])
        row.update(
            source_pressure_gpa=original["pressure_gpa"],
            pressure_correction_applied="0",
            pressure_graphical_reading_bound_gpa="",
            pressure_source="table_1",
            duplicate_of_source_row="1" if number == 32 else "",
            quality_notes=(
                "source_lattice_volume_mismatch; use printed volume column"
                if number in (1, 6, 11, 16, 32)
                else ""
            ),
            correction_note="",
        )
        if number == 15:
            row.update(
                pressure_gpa=f"{digitized_pressure_gpa:.2f}",
                pressure_correction_applied="1",
                pressure_graphical_reading_bound_gpa="0.05",
                pressure_source="figure_2_digitized",
                correction_note=(
                    "Replaces printed 6.95 GPa with independently digitized "
                    "Figure 2 pressure rounded to 0.01 GPa. Reading bound "
                    "+/-0.05 GPa is not an experimental standard uncertainty."
                ),
            )
        output.append(row)
    return output


def corrected_csv(digitized_pressure_gpa):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(corrected_rows(digitized_pressure_gpa))
    return stream.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = json.loads((ROOT / "docs/data/diopside-plot-check.json").read_text())
    payload = corrected_csv(report["row15"]["digitized_pressure_gpa"])
    if args.check:
        if OUTPUT.read_bytes() != payload.encode("utf-8"):
            raise SystemExit("Corrected diopside P-V-T dataset is stale")
        print("Corrected diopside P-V-T dataset reproduces")
    else:
        OUTPUT.write_text(payload, encoding="utf-8")
        print(OUTPUT)


if __name__ == "__main__":
    main()
