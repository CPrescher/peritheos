#!/usr/bin/env python3
"""Bundle the attributed candidate observations for the Fu (2023) CaSiO3 audit.

The original Gréaux workbook is a checksummed input only; the output contains
numerical measurements, source row identifiers, and explicitly derived values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from audit_fu_2023_casio3_refit import (
    GRE_AUX_SHA256,
    REGISTERED_REFIT_RECORD,
    _arrays,
    parse_greaux_table,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data"
PUBLISHED_RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"
# Preserve the existing public identifier and fit lineage.
DATASET = "ca_perovskite_fu_2023_candidate_composite_external"
COMPOSITE_PATH = "datasets/ca-perovskite-fu-2023-candidate-observations.csv"
GREAUX_PATH = "datasets/ca-perovskite-greaux-2019-figure3b.csv"
RIGHTS = (
    "CC0-1.0 applies only to Peritheos contributors' transcription, column "
    "naming, normalization, and arrangement. Source measurements remain "
    "attributed to Sun et al. (2016) and Gréaux et al. (2019); no CC license "
    "is asserted for their articles or source files. See "
    "datasets/ca-perovskite-fu-2023-observations.LICENSE.md."
)
FINDING = (
    "The 174 reconstructed candidate observations are bundled with attribution: "
    "140 Sun Table 1 P-V-T rows at 1200-2200 K and 34 cubic Gréaux Figure 3b "
    "P-density-Vp-Vs rows. Gréaux volume and elastic moduli are explicitly "
    "derived from the retained measurements. Fu does not disclose the exact "
    "row manifest, objective, or weights; the candidate selection and the "
    "non-parity audit are not an exact reconstruction of the published fit."
)


def write_csv(path: Path, rows: list[dict]) -> dict:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return {
        "path": str(path.relative_to(DATA)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "media_type": "text/csv",
    }


def load_sun() -> tuple[list[dict], list[dict]]:
    with (DATA / "datasets/ca-perovskite-sun-2016-table1-pvt.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        original = list(csv.DictReader(stream))
    rows = [
        {
            "pressure_gpa": float(row["pressure_gpa"]),
            "pressure_sigma_gpa": float(row["pressure_uncertainty_gpa"]),
            "temperature_k": float(row["temperature_k"]),
            "volume_a3": float(row["volume_a3_per_formula_unit"]),
            "volume_sigma_a3": float(row["volume_uncertainty_a3"]),
        }
        for row in original
    ]
    assert len(rows) == 144
    return original, rows


def bundle(workbook: Path) -> None:
    if hashlib.sha256(workbook.read_bytes()).hexdigest() != GRE_AUX_SHA256:
        raise ValueError("Gréaux workbook does not match the audited source")
    greaux = parse_greaux_table(workbook)
    original_sun, sun = load_sun()
    arrays = _arrays(sun, greaux)
    raw_resource = write_csv(
        DATA / GREAUX_PATH,
        [{"source_order": i, **row} for i, row in enumerate(greaux, 1)],
    )
    columns = [
        {
            "name": "source_id",
            "quantity": "source_identifier",
            "unit": "1",
            "role": "flag",
        },
        {
            "name": "source_order",
            "quantity": "source_row_order",
            "unit": "1",
            "role": "flag",
        },
    ]
    for name, quantity, unit in (
        ("pressure_gpa", "pressure", "GPa"),
        ("temperature_k", "temperature", "K"),
        (
            "volume_a3",
            "conventional_unit_cell_volume",
            "angstrom^3/conventional_unit_cell",
        ),
        ("bulk_modulus_gpa", "adiabatic_bulk_modulus", "GPa"),
        ("shear_modulus_gpa", "shear_modulus", "GPa"),
    ):
        columns.append(
            {"name": name, "quantity": quantity, "unit": unit, "role": "value"}
        )
        if name != "temperature_k":
            prefix, suffix = name.rsplit("_", 1)
            columns.append(
                {
                    "name": f"{prefix}_standard_deviation_{suffix}",
                    "quantity": quantity,
                    "unit": unit,
                    "role": "standard_deviation",
                    "of": name,
                }
            )
    rows = []
    for prefix, source_id, source_orders in (
        (
            "sun",
            "sun_2016_table1",
            [
                int(r["source_order"])
                for r in original_sun
                if 1200 <= float(r["temperature_k"]) <= 2200
            ],
        ),
        ("g", "greaux_2019_figure3b", list(range(1, 35))),
    ):
        for i, source_order in enumerate(source_orders):
            row = dict.fromkeys((column["name"] for column in columns), "")
            row.update(source_id=source_id, source_order=source_order)
            for column, key in (
                ("pressure_gpa", "p"),
                ("pressure_standard_deviation_gpa", "ps"),
                ("temperature_k", "t"),
                ("volume_a3", "v"),
                ("volume_standard_deviation_a3", "vs"),
                ("bulk_modulus_gpa", "k"),
                ("bulk_modulus_standard_deviation_gpa", "ks"),
                ("shear_modulus_gpa", "mu"),
                ("shear_modulus_standard_deviation_gpa", "mus"),
            ):
                if f"{prefix}_{key}" in arrays:
                    row[column] = float(arrays[f"{prefix}_{key}"][i])
            rows.append(row)
    assert len(rows) == 174
    resource = write_csv(DATA / COMPOSITE_PATH, rows)
    path = DATA / "materials/ca_perovskite.eosmat"
    document = json.loads(path.read_text(encoding="utf-8"))
    dataset = next(d for d in document["datasets"] if d["identifier"] == DATASET)
    dataset.pop("rows", None)
    dataset.update(
        kind="pressure_volume_temperature",
        description="174 reconstructed candidate observations for Fu (2023): 140 Sun P-V-T rows and 34 cubic Gréaux elastic rows; exact published fit selection is undisclosed.",
        columns=columns,
        resource=resource,
        license=RIGHTS,
        used_by_eos_records=[PUBLISHED_RECORD, REGISTERED_REFIT_RECORD],
        source_datasets=[
            "ca_perovskite_sun_2016_table1_pvt",
            "ca_perovskite_greaux_2019_figure3b",
        ],
        uncertainty={
            "source": "Reported source uncertainties are retained; no row-wise temperature uncertainty is supplied.",
            "derived": "Gréaux V, KS, and shear errors use first-order propagation assuming independent density, Vp, and Vs errors; derived observables have shared inputs and are not statistically independent.",
        },
        notes="The legacy identifier ending in _external is retained for compatibility; rows are now bundled. Sun row order references the existing 144-row transcription; rows above 2200 K are excluded. Gréaux source_order indexes the 34 cubic Figure 3b observations, excluding the tetragonal Figure 3a sheet. Uses measured PNaCl, not derived PFS. V=116.162/(0.602214076*rho), KS=rho*(Vp^2-4*Vs^2/3), mu=rho*Vs^2. Cubic CaSiO3 has Z=1, so formula-unit and conventional-cell volumes coincide. Empty Sun elastic cells denote unavailable observations. These are experimental constraints, not model-generated checkpoints.",
    )
    raw_id = "ca_perovskite_greaux_2019_figure3b"
    raw_columns = [
        {
            "name": "source_order",
            "quantity": "source_row_order",
            "unit": "1",
            "role": "flag",
        }
    ]
    for name in greaux[0]:
        unit = (
            "GPa"
            if "pressure" in name
            else "K"
            if "temperature" in name
            else "g/cm^3"
            if "density" in name
            else "km/s"
        )
        quantity = (
            "temperature"
            if name == "temperature_k"
            else name.replace("_sigma", "").rsplit("_", 1)[0]
        )
        column = {
            "name": name,
            "quantity": quantity,
            "unit": unit,
            "role": "uncertainty" if "_sigma" in name else "value",
        }
        if "_sigma" in name:
            column["of"] = name.replace("_sigma", "")
        raw_columns.append(column)
    raw = {
        "identifier": raw_id,
        "kind": "pressure_density_sound_velocity_temperature",
        "description": "All 34 cubic CaSiO3 measurements in Gréaux et al. (2019) source worksheet Figure 3b, retaining both pressure columns and reported uncertainties.",
        "reference": "Gréaux et al. (2019), Nature 565, 218–221, DOI 10.1038/s41586-018-0816-5, source worksheet Figure 3b",
        "source_location": "Official source-data workbook, Figure 3b; source_order counts numeric observation rows",
        "source_url": "https://doi.org/10.1038/s41586-018-0816-5",
        "source_artifact": {
            "sha256": GRE_AUX_SHA256,
            "download_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0816-5/MediaObjects/41586_2018_816_MOESM1_ESM.xlsx",
        },
        "columns": raw_columns,
        "resource": raw_resource,
        "license": RIGHTS,
        "used_by_eos_records": [PUBLISHED_RECORD, REGISTERED_REFIT_RECORD],
        "notes": "Factual numerical transcription only; no publisher workbook layout, prose, or figures are bundled. PNaCl is measured pressure; PFS is a derived finite-strain diagnostic. The separate Figure 3a tetragonal rows are outside the cubic candidate selection.",
    }
    document["datasets"] = [
        d for d in document["datasets"] if d["identifier"] != raw_id
    ] + [raw]
    for record in document["eos_records"]:
        if record["identifier"] not in (PUBLISHED_RECORD, REGISTERED_REFIT_RECORD):
            continue
        check = record["scientific_validation"]["primary_data_check"]
        check.update(
            status="bundled",
            audit_date="2026-09-23",
            dataset_identifiers=[DATASET],
            finding=FINDING,
        )
        if record["identifier"] == PUBLISHED_RECORD:
            record["supporting_datasets"] = [DATASET]
    path.write_text(
        json.dumps(document, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--greaux-xlsx", type=Path, required=True)
    bundle(parser.parse_args().greaux_xlsx)
