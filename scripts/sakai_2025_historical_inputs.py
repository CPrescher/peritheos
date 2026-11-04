"""Recover usable Table S8 constraints from existing source-owned numerical tables."""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from scripts.reproduce_sakai_2025_gold import DATA, read_rows, ruby2020_from_mao1986

FILES = (
    "copper-dewaele-2004-table1-compression.csv",
    "tungsten-dewaele-2004-table1-compression.csv",
    "platinum-dewaele-2004-table1-compression.csv",
    "rhenium-anzellini-2014-table3-compression.csv",
    "mgo-dewaele-2000-table2-pvt.csv",
    "molybdenum-dewaele-2008-table2-compression.csv",
    "dorfman-2012-tables-s1-s6-cocompression.csv",
    "iron-dewaele-2006-epaps-compression.csv",
    "sakai-2014-table2-pvt.csv",
)


def measured_pair(source, row, volume, calibrant, marker_volume):
    return {
        "source": source,
        "source_row": row,
        "target_volume_a3": float(volume),
        "calibrant": calibrant,
        "calibrant_volume_a3": float(marker_volume),
    }


def ruby_row(source, row, volume, reported_pressure, original_scale):
    if original_scale == "mao1986":
        pressure = ruby2020_from_mao1986(float(reported_pressure))
    elif original_scale == "dorogokupets2007":
        p = float(reported_pressure)
        shift = 2 * p / (1884 * (1 + np.sqrt(1 + 22 * p / 1884)))
        pressure = 1870 * shift * (1 + 5.63 * shift)
    else:
        raise ValueError(original_scale)
    return {
        "source": source,
        "source_row": row,
        "target_volume_a3": float(volume),
        "calibrant": "ruby2020",
        "calibrant_volume_a3": None,
        "reconstructed_pressure_gpa": float(pressure),
        "qualification": "Ruby coordinate inferred from rounded reported "
        + original_scale
        + " pressure, not raw fluorescence measurements.",
    }


def historical_rows(target):
    rows = []
    if target == "tungsten":
        cu = read_rows(FILES[0])
        w = read_rows(FILES[1])
        if len(cu) != len(w):
            raise ValueError("Dewaele2004 pair inventory differs")
        for i, (c, r) in enumerate(zip(cu, w), 1):
            # Table I's first two runs share the Cu/W/Al row coordinate.
            for key in ("ruby_pressure_classical_gpa", "ruby_pressure_revised_gpa"):
                if c[key] != r[key]:
                    raise ValueError("Dewaele2004 Cu/W source-row mismatch")
            rows.append(
                measured_pair(
                    "Dewaele2004 TableI Cu/W",
                    i,
                    2 * float(r["atomic_volume_a3"]),
                    "copper",
                    4 * float(c["atomic_volume_a3"]),
                )
            )
    elif target == "rhenium":
        for i, r in enumerate(read_rows(FILES[3]), 1):
            v = r["volume_a3_conventional_cell"]
            if r["run"] == "1" and r["pressure_gauge"] == "ruby_fluorescence":
                rows.append(
                    ruby_row(
                        "Anzellini2014 TableIII run1",
                        i,
                        v,
                        r["pressure_gpa"],
                        "dorogokupets2007",
                    )
                )
            elif r["run"] == "3" and r["pressure_gauge"] == "tungsten_eos":
                rows.append(
                    measured_pair(
                        "Anzellini2014 TableIII run3",
                        i,
                        v,
                        "tungsten",
                        float(r["tungsten_lattice_a_angstrom"]) ** 3,
                    )
                )
    elif target == "platinum":
        for i, r in enumerate(read_rows(FILES[2]), 1):
            rows.append(
                ruby_row(
                    "Dewaele2004 TableI Pt",
                    i,
                    4 * float(r["atomic_volume_a3"]),
                    r["ruby_pressure_classical_gpa"],
                    "mao1986",
                )
            )
    elif target == "mgo":
        for i, r in enumerate(read_rows(FILES[4]), 1):
            if float(r["temperature_k"]) == 300 and r["pt_lattice_parameter_a"]:
                rows.append(
                    measured_pair(
                        "Dewaele2000 TableII 300K",
                        i,
                        r["mgo_unit_cell_volume_a3"],
                        "platinum",
                        float(r["pt_lattice_parameter_a"]) ** 3,
                    )
                )
    elif target == "molybdenum":
        for i, r in enumerate(read_rows(FILES[5]), 1):
            rows.append(
                ruby_row(
                    "Dewaele2008 TableII Mo",
                    i,
                    2 * float(r["atomic_volume_a3"]),
                    r["ruby_pressure_mao_gpa"],
                    "mao1986",
                )
            )
    elif target == "gold":
        # Reuse the audited five-family input reconstruction, not fitted curve nodes.
        report = json.loads(
            (
                DATA.parents[2] / "docs/data/sakai-2025-gold-reconstruction.json"
            ).read_text(encoding="utf-8")
        )
        for r in report["pressure_reductions"]:
            if r["source"].startswith("Sakai") or not r["included_in_diagnostic"]:
                continue  # S4 enters separately; retain the established <200 GPa mask.
            row = (
                measured_pair(
                    r["source"],
                    r["source_row"],
                    r["volume_a3"],
                    r["calibrant"],
                    r["calibrant_volume_a3"],
                )
                if r["calibrant_volume_a3"] is not None
                else {
                    "source": r["source"],
                    "source_row": r["source_row"],
                    "target_volume_a3": r["volume_a3"],
                    "calibrant": "ruby2020",
                    "calibrant_volume_a3": None,
                    "reconstructed_pressure_gpa": r["calculated_pressure_gpa"],
                    "qualification": r["selection_note"],
                }
            )
            rows.append(row)
    elif target == "iron":
        for i, r in enumerate(read_rows(FILES[7]), 1):
            if r["phase"] == "hcp" and r["pressure_gauge"] == "tungsten":
                rows.append(
                    measured_pair(
                        "Dewaele2006 EPAPS hcp/W",
                        i,
                        2 * float(r["atomic_volume_a3_per_atom"]),
                        "tungsten",
                        float(r["tungsten_lattice_parameter_angstrom"]) ** 3,
                    )
                )
        for i, r in enumerate(read_rows(FILES[8]), 1):
            if (
                r["material"] == "iron"
                and float(r["temperature_k"]) == 300
                and r["nacl_a_angstrom"]
                and float(r["nacl_a_angstrom"]) < 4
            ):
                rows.append(
                    measured_pair(
                        "Sakai2014 TableII Fe/NaCl-B2",
                        i,
                        r["volume_a3_conventional_cell"],
                        "nacl_b2",
                        r["nacl_volume_a3"],
                    )
                )
    if target in ("mgo", "molybdenum"):
        groups = defaultdict(dict)
        for r in read_rows(FILES[6]):
            groups[(r["run"], r["observation_index"])][r["material"]] = r
        sample, marker = ("MgO", "Pt") if target == "mgo" else ("Mo", "MgO")
        calibrant = "platinum" if target == "mgo" else "mgo"
        for key, group in groups.items():
            if sample in group and marker in group:
                rows.append(
                    measured_pair(
                        "Dorfman2012 " + key[0],
                        key[1],
                        group[sample]["volume_a3_conventional_cell"],
                        calibrant,
                        group[marker]["volume_a3_conventional_cell"],
                    )
                )
    return rows
