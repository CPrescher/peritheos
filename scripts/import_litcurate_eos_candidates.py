"""Build a non-executable Peritheos candidate ledger from LitCurate output.

LitCurate records are discovery leads, not primary scientific authority.  This
script deliberately writes a review ledger rather than ``.eosmat`` files.  A
candidate can become executable only after the checks in
``docs/adding-materials-and-eos.md`` have been completed against the primary
publication and its official data.

The importer accepts LitCurate's merged ``database.json``, a flat CSV export,
or the public ``data.xlsx`` deposit used by Shakya et al. (2026).  XLSX reading
uses only the Python standard library so the curation tool does not add a
spreadsheet dependency to Peritheos.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import posixpath
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping
from xml.etree import ElementTree

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATERIALS = ROOT / "peritheos" / "data" / "materials"
DEFAULT_BACKLOG = ROOT / "docs" / "material-eos-candidates.md"

LEDGER_FORMAT = "peritheos.literature-eos-candidates"
LEDGER_VERSION = 1
LITCURATE_DATASET_DOI = "10.5281/zenodo.22118629"
LITCURATE_DATASET_URL = "https://doi.org/10.5281/zenodo.22118629"

ORIGIN_SOURCE = "Source-reported"
ORIGIN_CITATION = "Citation-reported"
UNSPECIFIED_VALUES = {"", "unknown", "unspecified", "none", "nan", "null"}
REJECTED_PUBLICATION_DOIS = {
    # Audited separately and deliberately excluded: the temperature-indexed BM3
    # rows are slices of a first-principles thermal calculation, not defensible
    # independent static EOS records for the current Peritheos model registry.
    "10.2138/am-2024-9562",
}

# These labels have an exact existing Peritheos equation family.  This is only
# a mechanical compatibility hint; the primary equation still has to be read.
SUPPORTED_MODELS = {
    "second-order Birch-Murnaghan": {
        "type": "BM2",
        "model": "birch_murnaghan_2",
    },
    "third-order Birch-Murnaghan": {
        "type": "BM3",
        "model": "birch_murnaghan_3",
    },
    "fourth-order Birch-Murnaghan": {
        "type": "BM4",
        "model": "birch_murnaghan_4",
    },
    "modified Tait": {
        "type": "ModifiedTait",
        "model": "modified_tait",
    },
    "Murnaghan": {"type": "Murnaghan", "model": "murnaghan"},
    "Vinet": {"type": "Vinet", "model": "vinet"},
}

AMBIGUOUS_MODELS = {
    "Birch-Murnaghan (order unspecified)",
    "finite strain (unspecified form)",
    "other named EOS",
    "unknown",
}

MODEL_WORK = {
    "AP2",
    "Keane",
    "Kunc",
    "Mie-Gr\u00fcneisen-Debye / thermal EOS",
    "spin-crossover EOS",
    "Stacey reciprocal K-primed",
}

NOT_EOS_MODELS = {"none (not an EOS fit)"}

FIELD_ALIASES = {
    "doi": ("doi",),
    "title": ("title",),
    "year": ("year",),
    "phase": ("phase",),
    "composition": ("composition",),
    "structure": ("structure",),
    "sample": ("sample",),
    "eos_model": ("eos_model",),
    "method": ("method",),
    "method_reported": ("method_reported",),
    "origin": ("origin",),
    "V0": ("V0", "Reported V0"),
    "V0_unit": ("V0_unit", "Reported V0_unit"),
    "V0_basis": ("V0_basis", "Reported V0_basis"),
    "V0_determination": ("V0_determination", "Reported V0_determination"),
    "V0_normalized_cm3_mol": ("V0(cm3/mol)",),
    "K0": ("K0", "Reported K0"),
    "K0_unit": ("K0_unit", "Reported K0_unit"),
    "K0_type": ("K0_type", "Reported K0_type"),
    "K0_determination": ("K0_determination", "Reported K0_determination"),
    "K0_normalized_gpa": ("K0(GPa)",),
    "Kp": ("Kp", "Reported Kp"),
    "Kp_unit": ("Kp_unit",),
    "Kp_determination": ("Kp_determination", "Reported Kp_determination"),
    "T_ref": ("T_ref",),
    "T_ref_unit": ("T_ref_unit",),
    "P_ref": ("P_ref",),
    "P_ref_unit": ("P_ref_unit",),
    "confidence": ("confidence",),
    "evidence": ("evidence",),
    "extra_info": ("extra_info",),
}


def normalize_doi(value: object) -> str:
    """Return a comparison-safe DOI without silently repairing its suffix."""
    if value is None:
        return ""
    doi = str(value).strip().lower()
    doi = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", doi)
    return doi.rstrip(".,; ")


def _is_present(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip().lower() not in UNSPECIFIED_VALUES
    return True


def _has_cell_value(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip().lower() not in {"", "nan", "null"}
    return True


def _first(record: Mapping[str, Any], names: Iterable[str]) -> Any:
    for name in names:
        if name in record and _has_cell_value(record[name]):
            return record[name]
    return None


def _canonical_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    expanded = dict(record)
    material = record.get("material")
    if isinstance(material, Mapping):
        for name in ("phase", "composition", "structure", "sample"):
            expanded.setdefault(name, material.get(name))

    equation = record.get("equation")
    if isinstance(equation, Mapping):
        expanded.setdefault(
            "eos_model", equation.get("reported_name") or equation.get("family")
        )

    parameters = record.get("parameters")
    if isinstance(parameters, list):
        parameter_aliases = {
            "v0": "V0",
            "vo": "V0",
            "k0": "K0",
            "ko": "K0",
            "b0": "K0",
            "kp": "Kp",
            "kprime": "Kp",
            "k0prime": "Kp",
            "koprime": "Kp",
        }
        for parameter in parameters:
            if not isinstance(parameter, Mapping):
                continue
            raw_name = (
                str(parameter.get("name", ""))
                .lower()
                .replace("′", "prime")
                .replace("'", "prime")
            )
            normalized_name = re.sub(r"[^a-z0-9]+", "", raw_name)
            target = parameter_aliases.get(normalized_name)
            if target is None or target in expanded:
                continue
            expanded[target] = parameter.get("reported_value")
            expanded[f"{target}_unit"] = parameter.get("unit")
            expanded[f"{target}_determination"] = parameter.get("determination")

    reference_state = record.get("reference_state")
    if isinstance(reference_state, Mapping):
        expanded.setdefault("T_ref", reference_state.get("temperature"))
        expanded.setdefault("T_ref_unit", reference_state.get("temperature_unit"))
        expanded.setdefault("P_ref", reference_state.get("pressure"))
        expanded.setdefault("P_ref_unit", reference_state.get("pressure_unit"))

    evidence_locations = record.get("evidence_locations")
    if isinstance(evidence_locations, list) and evidence_locations:
        expanded.setdefault("evidence", "structured_evidence_locations_present")

    result = {
        name: _first(expanded, aliases) for name, aliases in FIELD_ALIASES.items()
    }
    result["origin"] = {
        "this_study": ORIGIN_SOURCE,
        "source_reported": ORIGIN_SOURCE,
        "cited": ORIGIN_CITATION,
        "citation_reported": ORIGIN_CITATION,
        "unknown": "Unspecified",
    }.get(str(result["origin"] or "").lower(), result["origin"])
    result["all_parameters"] = parameters if isinstance(parameters, list) else None
    context_keys = (
        "record_scope",
        "equation",
        "reported_range",
        "pressure_calibration",
        "primary_data",
        "evidence_locations",
        "conflicts",
        "blockers",
    )
    result["candidate_context"] = {
        key: record[key] for key in context_keys if key in record
    }
    return result


def _column_index(cell_reference: str) -> int:
    letters = re.match(r"[A-Z]+", cell_reference.upper())
    if letters is None:
        raise ValueError(f"Invalid spreadsheet cell reference {cell_reference!r}")
    result = 0
    for character in letters.group(0):
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def _xlsx_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return ["".join(item.itertext()) for item in root]


def _xlsx_sheet_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    relationships = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    relationship_targets = {
        relation.attrib["Id"]: relation.attrib["Target"] for relation in relationships
    }
    relationship_id = None
    for sheet in workbook.iter():
        if sheet.tag.endswith("}sheet") and sheet.attrib.get("name") == sheet_name:
            relationship_id = next(
                value for key, value in sheet.attrib.items() if key.endswith("}id")
            )
            break
    if relationship_id is None:
        raise ValueError(f"Workbook has no sheet named {sheet_name!r}")
    target = relationship_targets[relationship_id].lstrip("/")
    return target if target.startswith("xl/") else posixpath.join("xl", target)


def _xlsx_value(cell: ElementTree.Element, shared_strings: list[str]) -> Any:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(
            child.text or "" for child in cell.iter() if child.tag.endswith("}t")
        )
    value_element = next((child for child in cell if child.tag.endswith("}v")), None)
    if value_element is None or value_element.text is None:
        return None
    value = value_element.text
    if cell_type == "s":
        return shared_strings[int(value)]
    if cell_type in {"str", "e"}:
        return value
    if cell_type == "b":
        return value == "1"
    try:
        number = float(value)
    except ValueError:
        return value
    return int(number) if number.is_integer() else number


def load_xlsx(path: Path, *, sheet_name: str = "database") -> list[dict[str, Any]]:
    """Read one flat XLSX worksheet using the Python standard library."""
    with zipfile.ZipFile(path) as archive:
        shared_strings = _xlsx_shared_strings(archive)
        root = ElementTree.fromstring(
            archive.read(_xlsx_sheet_path(archive, sheet_name))
        )
    rows: list[list[Any]] = []
    for row in root.iter():
        if not row.tag.endswith("}row"):
            continue
        values: list[Any] = []
        for cell in row:
            if not cell.tag.endswith("}c"):
                continue
            column = _column_index(cell.attrib["r"])
            if column >= len(values):
                values.extend([None] * (column + 1 - len(values)))
            values[column] = _xlsx_value(cell, shared_strings)
        rows.append(values)
    if not rows:
        return []
    headers = [str(value) if value is not None else "" for value in rows[0]]
    return [
        {
            header: values[index] if index < len(values) else None
            for index, header in enumerate(headers)
            if header
        }
        for values in rows[1:]
    ]


def load_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _envelope_payload(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        return {}
    payload = value.get("payload", value)
    return dict(payload) if isinstance(payload, Mapping) else {}


def load_json(path: Path) -> list[dict[str, Any]]:
    """Flatten LitCurate's merged envelope export into reported EOS rows."""
    document = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(document, list):
        return [dict(item) for item in document if isinstance(item, Mapping)]
    papers = document.get("papers") if isinstance(document, Mapping) else None
    if not isinstance(papers, list):
        raise ValueError("JSON must be a row list or LitCurate merged database export")

    flattened = []
    for paper in papers:
        if not isinstance(paper, Mapping):
            continue
        source = _envelope_payload(paper.get("source"))
        entries: list[Mapping[str, Any]] = []
        for value in paper.values():
            payload = _envelope_payload(value)
            candidate_entries = payload.get("eos_entries")
            if isinstance(candidate_entries, list):
                entries = [
                    entry for entry in candidate_entries if isinstance(entry, Mapping)
                ]
                break
        for entry in entries:
            flattened.append({**source, **entry})
    return flattened


def load_records(path: Path) -> list[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        return load_xlsx(path)
    if suffix == ".csv":
        return load_csv(path)
    if suffix == ".json":
        return load_json(path)
    raise ValueError("LitCurate source must be .xlsx, .csv, or .json")


def bundled_dois(materials_directory: Path) -> set[str]:
    dois = set()
    for path in sorted(materials_directory.glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        for record in document.get("eos_records", []):
            reference = record.get("reference")
            if isinstance(reference, Mapping):
                doi = normalize_doi(reference.get("doi"))
                if doi:
                    dois.add(doi)
    return dois


def backlog_dois(path: Path) -> set[str]:
    if not path.exists():
        return set()
    matches = re.findall(
        r"10\.\d{4,9}/[-._;()/:a-z0-9]+",
        path.read_text(encoding="utf-8").lower(),
    )
    return {match.rstrip(".,;)") for match in matches}


def equation_classification(model: object) -> tuple[str, dict[str, str] | None]:
    model_name = str(model).strip() if model is not None else ""
    if model_name in SUPPORTED_MODELS:
        return "existing_family_candidate", SUPPORTED_MODELS[model_name]
    if model_name in NOT_EOS_MODELS:
        return "not_an_eos_fit", None
    if model_name in MODEL_WORK:
        return "model_or_thermal_work", None
    if model_name in AMBIGUOUS_MODELS or not model_name:
        return "equation_identity_unresolved", None
    return "unrecognized_model_label", None


def _review_bucket(
    *,
    origin: object,
    doi: str,
    phase: object,
    equation_status: str,
    core_parameters_complete: bool,
    is_bundled: bool,
) -> str:
    if origin == ORIGIN_CITATION:
        return "citation_trace"
    if is_bundled:
        return "already_bundled_source"
    if equation_status == "not_an_eos_fit":
        return "not_an_eos_fit"
    if origin != ORIGIN_SOURCE or not doi:
        return "manual_triage"
    if equation_status == "model_or_thermal_work":
        return "model_work"
    if equation_status != "existing_family_candidate":
        return "equation_audit"
    if not core_parameters_complete or not _is_present(phase):
        return "primary_source_audit"
    return "review_first"


def _fingerprint(values: Iterable[object]) -> str:
    text = "\x1f".join("" if value is None else str(value).strip() for value in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _normalized_title(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def _parameter(record: Mapping[str, Any], name: str) -> dict[str, Any]:
    result = {"reported": record[name]}
    for suffix in ("unit", "basis", "type", "determination"):
        key = f"{name}_{suffix}"
        if record.get(key) is not None:
            result[suffix] = record[key]
    normalized_key = {
        "V0": "V0_normalized_cm3_mol",
        "K0": "K0_normalized_gpa",
    }.get(name)
    if normalized_key and record.get(normalized_key) is not None:
        result["litcurate_normalized"] = record[normalized_key]
        result["litcurate_normalized_unit"] = {
            "V0": "cm^3/mol",
            "K0": "GPa",
        }[name]
    return result


def build_ledger(
    source_records: Iterable[Mapping[str, Any]],
    *,
    source_path: Path,
    material_dois: set[str],
    candidate_backlog_dois: set[str],
) -> dict[str, Any]:
    records = []
    duplicate_keys: dict[tuple[object, ...], list[dict[str, Any]]] = defaultdict(list)

    for source_row, raw_record in enumerate(source_records, start=2):
        record = _canonical_fields(raw_record)
        doi = normalize_doi(record["doi"])
        if doi in REJECTED_PUBLICATION_DOIS:
            continue
        equation_status, suggested_mapping = equation_classification(
            record["eos_model"]
        )
        complete = all(_is_present(record[name]) for name in ("V0", "K0", "Kp"))
        is_bundled = bool(doi and doi in material_dois)
        is_backlog = bool(doi and doi in candidate_backlog_dois)
        bucket = _review_bucket(
            origin=record["origin"],
            doi=doi,
            phase=record["phase"],
            equation_status=equation_status,
            core_parameters_complete=complete,
            is_bundled=is_bundled,
        )
        identifier = "litcurate_" + _fingerprint(
            (
                doi,
                record["title"],
                record["year"],
                record["phase"],
                record["composition"],
                record["sample"],
                record["eos_model"],
                record["V0"],
                record["K0"],
                record["Kp"],
                record["T_ref"],
                source_row,
            )
        )
        flags = []
        if is_backlog:
            flags.append("already_in_documented_backlog")
        if record["eos_model"] == "fourth-order Birch-Murnaghan":
            flags.append("higher_order_coefficient_not_first_class_in_source_schema")
        if record["origin"] == ORIGIN_CITATION:
            flags.append("locate_and_audit_underlying_primary_publication")
        if not _is_present(record["evidence"]):
            flags.append("litcurate_evidence_missing")

        candidate = {
            "identifier": identifier,
            "source_row": source_row,
            "publication": {
                "doi": doi or None,
                "title": record["title"],
                "year": record["year"],
            },
            "material": {
                "phase": record["phase"],
                "composition": record["composition"],
                "structure": record["structure"],
                "sample": record["sample"],
            },
            "reported_eos": {
                "model_label": record["eos_model"],
                "method": record["method"],
                "method_reported": record["method_reported"],
                "origin": record["origin"],
                "parameters": {
                    name: _parameter(record, name) for name in ("V0", "K0", "Kp")
                },
                "all_parameters": record["all_parameters"],
                "reference_state": {
                    "temperature": record["T_ref"],
                    "temperature_unit": record["T_ref_unit"],
                    "pressure": record["P_ref"],
                    "pressure_unit": record["P_ref_unit"],
                },
            },
            "litcurate": {
                "confidence": record["confidence"],
                "evidence_present": _is_present(record["evidence"]),
                "extra_info_present": _is_present(record["extra_info"]),
                "candidate_context": record["candidate_context"] or None,
            },
            "classification": {
                "review_bucket": bucket,
                "equation_status": equation_status,
                "suggested_peritheos_mapping": suggested_mapping,
                "core_parameters_complete": complete,
                "bundled_source_doi": is_bundled,
                "documented_backlog_doi": is_backlog,
                "flags": flags,
            },
        }
        records.append(candidate)
        duplicate_key = (
            _normalized_title(record["title"]),
            record["year"],
            str(record["phase"] or "").lower(),
            str(record["composition"] or "").lower(),
            str(record["eos_model"] or "").lower(),
            str(record["V0"]),
            str(record["K0"]),
            str(record["Kp"]),
        )
        duplicate_keys[duplicate_key].append(candidate)

    for group in duplicate_keys.values():
        group_dois = {item["publication"]["doi"] for item in group}
        if len(group_dois - {None}) <= 1:
            continue
        duplicate_group = "possible_duplicate_" + _fingerprint(
            sorted(doi or "" for doi in group_dois)
        )
        for item in group:
            classification = item["classification"]
            classification["possible_duplicate_group"] = duplicate_group
            classification["flags"].append("possible_duplicate_publication_identifier")

    bucket_counts = Counter(
        record["classification"]["review_bucket"] for record in records
    )
    origin_counts = Counter(record["reported_eos"]["origin"] for record in records)
    source_dois = {
        record["publication"]["doi"]
        for record in records
        if record["reported_eos"]["origin"] == ORIGIN_SOURCE
        and record["publication"]["doi"]
    }
    new_source_dois = {doi for doi in source_dois if doi not in material_dois}
    return {
        "format": LEDGER_FORMAT,
        "format_version": LEDGER_VERSION,
        "source": {
            "name": "Lower-Mantle Equation-of-State Data",
            "doi": LITCURATE_DATASET_DOI,
            "url": LITCURATE_DATASET_URL,
            "license": "CC-BY-4.0",
            "input_file": source_path.name,
            "input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "worksheet": "database" if source_path.suffix.lower() == ".xlsx" else None,
        },
        "policy": (
            "Discovery-only ledger. No entry is executable authority. Audit the "
            "primary publication and satisfy docs/adding-materials-and-eos.md "
            "before creating or changing an eosmat record."
        ),
        "summary": {
            "records": len(records),
            "unique_publication_dois": len(
                {
                    record["publication"]["doi"]
                    for record in records
                    if record["publication"]["doi"]
                }
            ),
            "source_reported_publication_dois": len(source_dois),
            "source_reported_publication_dois_not_bundled": len(new_source_dois),
            "by_origin": dict(
                sorted(origin_counts.items(), key=lambda item: str(item[0]))
            ),
            "by_review_bucket": dict(sorted(bucket_counts.items())),
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="LitCurate .json, .csv, or .xlsx")
    parser.add_argument("destination", type=Path, help="Candidate ledger JSON")
    parser.add_argument("--materials", type=Path, default=DEFAULT_MATERIALS)
    parser.add_argument("--backlog", type=Path, default=DEFAULT_BACKLOG)
    args = parser.parse_args()

    source_records = load_records(args.source)
    ledger = build_ledger(
        source_records,
        source_path=args.source,
        material_dois=bundled_dois(args.materials),
        candidate_backlog_dois=backlog_dois(args.backlog),
    )
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text(
        json.dumps(ledger, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    summary = ledger["summary"]
    print(
        f"Wrote {summary['records']} candidates from "
        f"{summary['unique_publication_dois']} publications to {args.destination}"
    )


if __name__ == "__main__":
    main()
