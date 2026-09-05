"""Regenerate the repository-wide material/Hugoniot inventory CSV."""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos" / "data" / "materials"
OUTPUT = ROOT / "docs" / "data" / "hugoniot-material-inventory.csv"
FIELDS = (
    "material_identifier",
    "material_name",
    "formula",
    "represented_phase_or_polymorph",
    "phase_identity_source",
    "symmetry",
    "space_group",
    "formula_units_per_cell",
    "equation_kinds",
    "eos_record_count",
    "input_file",
)


def equation_kind(record: dict) -> str:
    """Return the explicit or structurally implied equation category."""
    if "equation_kind" in record:
        return str(record["equation_kind"])
    if record["eos"]["type"] == "LinearUsUpHugoniot":
        return "hugoniot"
    return "thermal" if record.get("thermal") is not None else "isothermal"


def main() -> None:
    rows = []
    for input_path in sorted(MATERIALS.glob("*.eosmat")):
        document = json.loads(input_path.read_text(encoding="utf-8"))
        phase = document.get("phase")
        kinds = tuple(dict.fromkeys(equation_kind(r) for r in document["eos_records"]))
        rows.append(
            {
                "material_identifier": document["identifier"],
                "material_name": document["name"],
                "formula": document["formula"],
                "represented_phase_or_polymorph": phase or document["name"],
                "phase_identity_source": (
                    "top-level phase field"
                    if phase
                    else "material name plus crystal structure"
                ),
                "symmetry": document.get("symmetry", ""),
                "space_group": document.get("space_group", ""),
                "formula_units_per_cell": document.get(
                    "formula_units_per_cell", ""
                ),
                "equation_kinds": ";".join(kinds),
                "eos_record_count": len(document["eos_records"]),
                "input_file": input_path.relative_to(ROOT).as_posix(),
            }
        )

    with OUTPUT.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
