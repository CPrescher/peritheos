"""Import the Dioptas material database as Peritheos-owned ``.eosmat`` data.

This is a provenance-preserving mechanical migration tool. It does not confer
primary-source validation on imported EOS parameters.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

DIOPTAS_VERSION = "0.10.0"
DIOPTAS_COMMIT = "5a8bfd81d10bfab3499039603380aae34576d60a"
DIOPTAS_URL = "https://github.com/Dioptas/Dioptas"
RT_MODELS = {
    "BM2": "birch_murnaghan_2",
    "BM3": "birch_murnaghan_3",
    "BM4": "birch_murnaghan_4",
    "Holzapfel": "holzapfel",
    "Murnaghan": "murnaghan",
    "Vinet": "vinet",
}
THERMAL_MODELS = {
    "AlphaKT": "thermal_reference_state",
    "MieGruneisenDebye": "mie_gruneisen_debye",
    "MieGruneisenEinstein": "mie_gruneisen_einstein",
    "Sokolova2016": "multi_oscillator_gruneisen_thermal_pressure",
}
FEI_2007_DOI = "10.1073/pnas.0609013104"


def _normalized_doi(reference: object) -> str | None:
    """Return a normalized DOI from structured reference metadata."""
    if not isinstance(reference, dict):
        return None
    doi = reference.get("doi")
    if not isinstance(doi, str):
        return None
    return doi.lower().removeprefix("https://doi.org/").removeprefix("doi:")


def _correct_debye_temperature_law(record: dict) -> dict | None:
    """Correct a source Debye-temperature law when primary evidence requires it.

    Dioptas 0.10.0 classified Fei et al. (2007) as the generic integrated
    constant-q Mie-Gruneisen-Debye model. Equation 3 and the definition
    immediately following it instead use a variable Gruneisen exponent in the
    Debye-temperature power law.
    """
    thermal = record.get("thermal")
    if not isinstance(thermal, dict):
        return None
    if (
        thermal.get("type") != "MieGruneisenDebye"
        or _normalized_doi(record.get("reference")) != FEI_2007_DOI
    ):
        return None

    thermal["debye_temperature_law"] = "variable_exponent"
    return {
        "path": "thermal.debye_temperature_law",
        "source_value": "integrated_gruneisen (implicit Dioptas default)",
        "value": "variable_exponent",
        "reason": (
            "Fei et al. (2007) use theta_D = theta0*(V/V0)^(-gamma(V)), "
            "not the integrated constant-q Debye-temperature relation."
        ),
        "primary_reference": {
            "doi": FEI_2007_DOI,
            "location": "Equation 3 and the definition immediately following it",
        },
    }


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_") or "record"


def _record_identifier(material_identifier: str, record: dict, index: int) -> str:
    reference = record.get("reference") or {}
    authors = reference.get("authors") or []
    author = authors[0] if authors else "unknown"
    year = reference.get("year", "undated")
    eos_type = (record.get("eos") or {}).get("type") or "eos"
    return _slug(f"{material_identifier}_{author}_{year}_{eos_type}_{index + 1}")


def migrate_document(path: Path) -> dict:
    """Return one Dioptas format-2 document with additive Peritheos fields."""
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("format_version") != 2:
        raise ValueError(f"{path}: expected Dioptas material format version 2")

    material_identifier = _slug(path.stem)
    document = dict(document)
    document["format"] = "peritheos.material"
    document["format_version"] = 3
    document["identifier"] = material_identifier
    document["units"] = {
        "pressure": "GPa",
        "temperature": "K",
        "volume": "angstrom^3/conventional_unit_cell",
    }

    identifiers: set[str] = set()
    migrated_records = []
    for index, source_record in enumerate(document.get("eos_records", [])):
        record = dict(source_record)
        record["eos"] = dict(record["eos"])
        record["eos"]["model"] = RT_MODELS[record["eos"]["type"]]
        if record.get("thermal") is not None:
            record["thermal"] = dict(record["thermal"])
            correction = _correct_debye_temperature_law(record)
            record["thermal"]["model"] = THERMAL_MODELS[record["thermal"]["type"]]
            if correction is not None:
                record.setdefault("migration_corrections", []).append(correction)
        identifier = _record_identifier(material_identifier, record, index)
        if identifier in identifiers:
            raise ValueError(f"{path}: duplicate generated record ID {identifier}")
        identifiers.add(identifier)
        record["identifier"] = identifier
        record["scientific_validation"] = {
            "status": "pending_primary_source_check",
            "note": (
                "Migrated for catalog and file-format compatibility. Peritheos "
                "has not independently revalidated this record against the cited "
                "primary publication unless a separately audited catalog record "
                "states otherwise."
            ),
            "migration_source": {
                "project": "Dioptas",
                "version": DIOPTAS_VERSION,
                "commit": DIOPTAS_COMMIT,
                "file": path.name,
            },
        }
        migrated_records.append(record)
    document["eos_records"] = migrated_records
    return document


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    source_files = sorted(args.source.glob("*.json"))
    if len(source_files) != 120:
        raise ValueError(f"Expected 120 Dioptas materials, found {len(source_files)}")
    args.destination.mkdir(parents=True, exist_ok=True)

    for stale_path in args.destination.glob("*.eosmat"):
        stale_path.unlink()

    material_count = 0
    record_count = 0
    for source_path in source_files:
        document = migrate_document(source_path)
        if not document["eos_records"]:
            continue
        material_count += 1
        record_count += len(document["eos_records"])
        destination = args.destination / f"{source_path.stem}.eosmat"
        destination.write_text(
            json.dumps(document, indent=1, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    if record_count != 147:
        raise ValueError(f"Expected 147 Dioptas EOS records, found {record_count}")
    if material_count != 116:
        raise ValueError(f"Expected 116 EOS materials, found {material_count}")

    manifest = {
        "format": "peritheos.material-library-manifest",
        "format_version": 1,
        "source": {
            "project": "Dioptas",
            "url": DIOPTAS_URL,
            "version": DIOPTAS_VERSION,
            "commit": DIOPTAS_COMMIT,
            "license": "MIT",
        },
        "materials": material_count,
        "eos_records": record_count,
        "scientific_validation": (
            "Migration preserves Dioptas data and attribution but does not replace "
            "Peritheos primary-source validation."
        ),
    }
    (args.destination / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
