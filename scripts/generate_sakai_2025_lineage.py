"""Build a separate calibration-evidence graph from primary numerical resources.

Repository maintenance tool. Does not read Table S8, raw workbook archives or
Sakai reconstruction reports. This is a broader, partial evidence graph, not
recovery of the authors' exact data selection, reduction or fitting sequence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
SOURCE = DATA / "sakai-2025-independent-lineage-sources.json"
OUTPUT = DATA / "sakai-2025-independent-calibration-lineage.csv"


def generate(source=SOURCE, data=DATA):
    configuration = json.loads(source.read_bytes())
    result = []

    def read(name):
        info = configuration["sources"][name]
        payload = (data / info["file"]).read_bytes()
        if hashlib.sha256(payload).hexdigest() != info["sha256"]:
            raise ValueError(f"Changed lineage input: {name}")
        return list(csv.DictReader(io.StringIO(payload.decode())))

    def emit(name, target, calibrant, relation, selected, fields, selection):
        if not selected:
            raise ValueError(f"Empty lineage edge: {name}/{target}/{calibrant}")
        info = configuration["sources"][name]
        fingerprint = hashlib.sha256(
            json.dumps(selected, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        result.append(
            {
                "target_material": target,
                "calibrant_or_partner": calibrant,
                "relation": relation,
                "source_doi": info["doi"],
                "source_location": info["source_location"],
                "source_resource": "datasets/" + info["file"],
                "source_sha256": info["sha256"],
                "observations": len(selected),
                "evidence_fields": ";".join(fields),
                "selection": selection,
                "selected_rows_sha256": fingerprint,
            }
        )

    # New Sakai experimental pairs: do not use the overlapping 326-row alias.
    for table, target, marker in [
        (1, "rhenium", "copper"),
        (2, "platinum", "rhenium"),
        (3, "tungsten", "copper"),
        (5, "molybdenum", "tungsten"),
        (6, "iron", "tungsten"),
        (7, "mgo", "platinum"),
    ]:
        name = f"sakai_s{table}"
        emit(
            name,
            target,
            marker,
            "paired_volumes",
            read(name),
            ["target_volume_source_cell", "calibrant_volume_source_cell"],
            "All present-study paired rows; not a recovered fit mask",
        )
    emit(
        "sakai_s4",
        "gold",
        "tungsten",
        "paired_volumes",
        read("sakai_s4"),
        ["volume_1_a3", "volume_2_a3"],
        "All present-study Au/W pairs",
    )

    # Ordinary reported gauge and co-measurement fields from original studies.
    for name, target, marker, fields in [
        (
            "dewaele2004_au",
            "gold",
            "ruby",
            [
                "atomic_volume_a3",
                "ruby_pressure_classical_gpa",
                "ruby_pressure_revised_gpa",
            ],
        ),
        (
            "dewaele2004_pt",
            "platinum",
            "ruby",
            [
                "atomic_volume_a3",
                "ruby_pressure_classical_gpa",
                "ruby_pressure_revised_gpa",
            ],
        ),
        (
            "dewaele2004_cu",
            "copper",
            "ruby",
            [
                "atomic_volume_a3",
                "ruby_pressure_classical_gpa",
                "ruby_pressure_revised_gpa",
            ],
        ),
        (
            "dewaele2008_mo",
            "molybdenum",
            "ruby",
            ["atomic_volume_a3", "ruby_pressure_mao_gpa"],
        ),
        (
            "takemura2008_au",
            "gold",
            "ruby",
            [
                "lattice_111_angstrom",
                "ruby_pressure_zha_gpa",
                "ruby_pressure_dorogokupets_gpa",
            ],
        ),
        ("sakai2011_nacl", "nacl_b2", "platinum", ["volume_a3", "pt_volume_a3"]),
        ("dewaele2018_au", "gold", "rhenium", ["au_volume_a3", "re_volume_a3"]),
        ("hirose2008_au", "gold", "mgo", ["volume_a3", "mgo_volume_a3"]),
        (
            "dewaele2000_mgo",
            "mgo",
            "platinum",
            ["mgo_unit_cell_volume_a3", "pt_lattice_parameter_a"],
        ),
        ("sata2002_nacl", "nacl_b2", "platinum", ["volume_a3", "pt_d111_angstrom"]),
        ("sata2002_nacl", "nacl_b2", "mgo", ["volume_a3", "mgo_d111_angstrom"]),
    ]:
        selected = [r for r in read(name) if all(r[k] for k in fields)]
        relation = "reported_ruby_pressure" if marker == "ruby" else "co_measured"
        emit(
            name,
            target,
            marker,
            relation,
            selected,
            fields,
            "All rows with the named fields; no Sakai exclusion or pressure recalibration applied",
        )

    copper, tungsten = read("dewaele2004_cu"), read("dewaele2004_w")
    if len(copper) != len(tungsten) or any(
        c[k] != w[k]
        for c, w in zip(copper, tungsten)
        for k in ("ruby_pressure_classical_gpa", "ruby_pressure_revised_gpa")
    ):
        raise ValueError("Dewaele Cu/W source coordinate differs")
    emit(
        "dewaele2004_w",
        "tungsten",
        "copper",
        "co_measured",
        [{"copper": c, "tungsten": w} for c, w in zip(copper, tungsten)],
        [
            "atomic_volume_a3",
            "ruby_pressure_classical_gpa",
            "ruby_pressure_revised_gpa",
        ],
        "Paired separate Cu/W resources with both reported ruby coordinates equal",
    )

    for gauge, marker in [
        ("ruby_fluorescence", "ruby"),
        ("helium_eos", "helium"),
        ("tungsten_eos", "tungsten"),
    ]:
        rows = [r for r in read("anzellini2014_re") if r["pressure_gauge"] == gauge]
        emit(
            "anzellini2014_re",
            "rhenium",
            marker,
            "reported_pressure_gauge",
            rows,
            ["run", "pressure_gauge", "pressure_gpa"],
            f"pressure_gauge={gauge}",
        )
    for gauge, marker in [("ruby", "ruby"), ("tungsten", "tungsten")]:
        rows = [
            r
            for r in read("dewaele2006_fe")
            if r["phase"] == "hcp" and r["pressure_gauge"] == gauge
        ]
        emit(
            "dewaele2006_fe",
            "iron",
            marker,
            "reported_pressure_gauge",
            rows,
            [
                "phase",
                "pressure_gauge",
                "ruby_wavelength_angstrom",
                "tungsten_lattice_parameter_angstrom",
            ],
            f"hcp phase; pressure_gauge={gauge}",
        )
    for gauge, marker in [("Ne", "neon"), ("MgO", "mgo"), ("NaCl-B2", "nacl_b2")]:
        rows = [r for r in read("fei2016_fe") if r["pressure_standard"] == gauge]
        emit(
            "fei2016_fe",
            "iron",
            marker,
            "reported_pressure_gauge",
            rows,
            ["pressure_standard", "marker_a_angstrom"],
            f"pressure_standard={gauge}",
        )
    for marker, field in [("nacl", "nacl_volume_a3"), ("mgo", "mgo_volume_a3")]:
        rows = [r for r in read("sakai2014_fe") if r["material"] == "iron" and r[field]]
        if not rows:
            continue  # Missing pure-Fe MgO pairs are not invented.
        emit(
            "sakai2014_fe",
            "iron",
            marker,
            "co_measured",
            rows,
            ["material", "temperature_k", field],
            "Pure iron rows; marker present; all temperatures; NaCl phase not inferred",
        )

    # Build every co-observed pair involving the nine target materials, rather
    # than importing a literature compiler's selected direction or phase mask.
    groups = defaultdict(dict)
    for row in read("dorfman2012"):
        group = groups[(row["run"], row["observation_index"])]
        if row["material"] in group:
            raise ValueError("Duplicate material in a source co-measurement group")
        group[row["material"]] = row
    names = {
        "Mo": "molybdenum",
        "MgO": "mgo",
        "Au": "gold",
        "NaCl B2": "nacl_b2",
        "Pt": "platinum",
        "Ne": "neon",
    }
    symbols = list(names)
    for i, a in enumerate(symbols):
        for b in symbols[i + 1 :]:
            pairs = [
                {a: group[a], b: group[b]}
                for group in groups.values()
                if a in group and b in group
            ]
            if pairs:
                emit(
                    "dorfman2012",
                    names[a],
                    names[b],
                    "co_measured",
                    pairs,
                    [
                        "run",
                        "observation_index",
                        "material",
                        "volume_a3_conventional_cell",
                    ],
                    "Same source run and observation index; direction is display only, not a fitted calibration",
                )
    result.sort(
        key=lambda r: (
            r["target_material"],
            r["source_doi"],
            r["calibrant_or_partner"],
            r["relation"],
        )
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(result[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(result)
    return stream.getvalue().encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = generate()
    if args.check:
        if OUTPUT.read_bytes() != payload:
            raise SystemExit("Independent Sakai lineage is stale")
    else:
        OUTPUT.write_bytes(payload)
    print(
        f"Independent lineage: {len(list(csv.DictReader(io.StringIO(payload.decode()))))} edges"
    )


if __name__ == "__main__":
    main()
