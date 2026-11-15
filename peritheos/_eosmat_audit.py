"""EOSMAT v5 record ``check`` headlines and per-material audit sidecars.

A bundled material is stored as two files:

``peritheos/data/materials/<id>.eosmat``
    The readable record file: material header, EOS records (model,
    coefficients, reference state, ranges, notes, default flags and a short
    ``check`` block), then datasets.
``peritheos/data/audits/<id>.audit.json``
    The audit trail of every record, keyed by record identifier: ``evidence``,
    the detailed ``pressure_calibration`` and the other audit-only blocks
    listed in :data:`AUDIT_RECORD_KEYS`.

:func:`merge_audit` restores the complete inline v5 document. Standalone files
may keep all of this inline; a sidecar is never required for authoring.

Nothing here reclassifies a record. ``check.reproduction`` is a fixed mapping of
the stored assessment status and ``check.warning`` is built only from stored
coded restrictions (see :func:`usage_warning`).
"""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any

from peritheos.errors import EosmatError

AUDIT_PACKAGE = "peritheos.data.audits"
AUDIT_SUFFIX = ".audit.json"
AUDIT_FORMAT = "peritheos.eosmat-audit"

#: Record keys that describe how a record was checked, not what it is. They
#: move to the sidecar; every other record key stays in the record file.
AUDIT_RECORD_KEYS = (
    "evidence",
    "pressure_calibration",
    "fit_input_inventory",
    "network_reconstruction_audit",
    "reproduction_audit",
    "numerical_reproduction",
    "validation_selection",
    "migration_corrections",
    "scope_history",
    "author_input_followup",
)

#: Fixed mapping from ``evidence.reproduction.assessment.status``.
REPRODUCTION_HEADLINES = {
    "parity": "reproduced",
    "similar": "close",
    "not_reproduced": "disagrees",
    "not_refittable": "untested",
    "source_reconstruction": "untested",
    "not_assessed": "untested",
}
REPRODUCTION_VALUES = ("reproduced", "close", "disagrees", "untested")
SOURCE_STATUSES = (
    "primary_source_validated",
    "pending_primary_source_check",
    "not_reproduced",
    "deferred",
)
CHECK_KEYS = (
    "reproduction",
    "warning",
    "finding",
    "source",
    "source_checked_on",
    "audit",
)

#: Statements that apply to every assessment of their kind. They are defined
#: once in docs/eos-labels.md and are not repeated in bundled sidecars.
GENERIC_QUALIFICATIONS = (
    "This named numerical assessment alone does not establish author-workflow "
    "reproduction, physical accuracy or recommendation for use.",
    "No reconstructed coefficient table is recorded for this assessment; no "
    "missing values have been inferred.",
)
GENERIC_VALIDATION_NOTES = (
    "Independently checked against the cited primary publication (published "
    "article or official/author copy); no external software catalog was used as "
    "scientific authority.",
)

_REPORT_URL = (
    "https://peritheos.readthedocs.io/en/latest/reconstruction-reports/"
    "{material}/#record-{record}"
)
_REPORT_PATH = "docs/reconstruction-reports/{material}.md#record-{record}"
_REPORT_LEDGER = "docs/data/primary-eos-refits.json"

#: Short fixed wording for coded ``scientific_validation.usage_recommendation``
#: values. Free-text recommendations are never turned into a warning.
USAGE_RECOMMENDATION_WARNINGS = {
    "not_recommended": "Not recommended for use.",
    "not_recommended_for_quantitative_use": "Not recommended for quantitative use.",
    "consult_source_audit_before_quantitative_use": (
        "Consult the source audit before quantitative use."
    ),
    "qualified_use": "Qualified use only; see the audit.",
    "on_hold_not_recommended_pending_source_mapping_review": (
        "On hold; not recommended pending source mapping review."
    ),
}
_DO_NOT_USE_WARNING = "Do not use."
_HOLD_WARNING = "Default selection on hold pending review."
_EXPLICIT_SELECTION_WARNING = "Explicit selection only; not loaded by default."


def derived_report_links(material: str, record: str) -> dict[str, str]:
    """Return the reconstruction-report links implied by the identifiers."""
    return {
        "report_url": _REPORT_URL.format(material=material, record=record),
        "report_path": _REPORT_PATH.format(material=material, record=record),
        "report": _REPORT_LEDGER,
    }


def audit_file_name(material: str) -> str:
    return f"{material}{AUDIT_SUFFIX}"


def bundled_audit_reference(material: str, record: str) -> str:
    """Relative reference from a bundled record file to its sidecar entry."""
    return f"../audits/{audit_file_name(material)}#{record}"


def _evidence(record: Mapping[str, Any]) -> Mapping[str, Any]:
    """Evidence mapping of a native record, or the flattened v3 record itself."""
    evidence = record.get("evidence")
    return evidence if isinstance(evidence, Mapping) else record


def _section(value: Any, key: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        section = value.get(key)
        if isinstance(section, Mapping):
            return section
    return {}


def _mutable_section(value: Any, key: str) -> dict[str, Any]:
    if isinstance(value, dict):
        section = value.get(key)
        if isinstance(section, dict):
            return section
    return {}


def reproduction_headline(record: Mapping[str, Any]) -> str:
    """Map the stored assessment status; anything else is ``untested``."""
    assessment = _section(_section(_evidence(record), "reproduction"), "assessment")
    return REPRODUCTION_HEADLINES.get(str(assessment.get("status")), "untested")


def _hold_covers(hold: Any, identifier: str) -> bool:
    if not isinstance(hold, Mapping) or hold.get("status") != "on_hold":
        return False
    records = hold.get("records")
    if isinstance(records, list):
        return identifier in records
    return hold.get("record") == identifier


def usage_warning(
    record: Mapping[str, Any], material: Mapping[str, Any] | None = None
) -> str | None:
    """Join the fixed phrases for stored, coded usage restrictions.

    Sources, in order: ``evidence.reproduction.do_not_use`` (true), a coded
    ``evidence.scientific_validation.usage_recommendation``, a material-level
    ``recommendation_hold`` naming the record, and ``catalog_access:
    explicit_selection``. Prose fields are never interpreted.
    """
    evidence = _evidence(record)
    phrases: list[str] = []
    if _section(evidence, "reproduction").get("do_not_use") is True:
        phrases.append(_DO_NOT_USE_WARNING)
    stored = _section(evidence, "scientific_validation").get("usage_recommendation")
    recommendation = stored if isinstance(stored, str) else ""
    # Coded values are lower_snake_case words; prose is never interpreted.
    coded = recommendation.replace("_", "").isalpha() and recommendation.islower()
    if coded:
        phrases.append(
            USAGE_RECOMMENDATION_WARNINGS.get(
                recommendation, recommendation.replace("_", " ").capitalize() + "."
            )
        )
    on_hold_code = coded and recommendation.startswith("on_hold")
    if (
        material is not None
        and not on_hold_code
        and _hold_covers(material.get("recommendation_hold"), str(record["identifier"]))
    ):
        phrases.append(_HOLD_WARNING)
    if record.get("catalog_access") == "explicit_selection":
        phrases.append(_EXPLICIT_SELECTION_WARNING)
    return " ".join(phrases) or None


def check_block(
    record: Mapping[str, Any],
    material: Mapping[str, Any] | None = None,
    *,
    audit: str | None = None,
) -> dict[str, Any]:
    """Derive a record's ``check`` block from its inline evidence."""
    validation = _section(_evidence(record), "scientific_validation")
    result: dict[str, Any] = {
        "reproduction": reproduction_headline(record),
        "warning": usage_warning(record, material),
    }
    finding = _section(validation, "primary_source_check").get("finding")
    if isinstance(finding, str) and finding:
        result["finding"] = finding
    if "status" in validation:
        result["source"] = validation["status"]
    if "audit_date" in validation:
        result["source_checked_on"] = validation["audit_date"]
    if audit is not None:
        result["audit"] = audit
    return result


def headline(
    record: Mapping[str, Any], material: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Effective ``reproduction`` and ``warning`` of a record.

    A stored ``check`` value wins; otherwise the value is derived from inline
    evidence exactly as :func:`check_block` would. No assessment means
    ``untested``.
    """
    check = record.get("check")
    check = check if isinstance(check, Mapping) else {}
    return {
        "reproduction": check.get("reproduction") or reproduction_headline(record),
        "warning": (
            check["warning"] if "warning" in check else usage_warning(record, material)
        ),
    }


def validate_check(check: Any, location: str) -> None:
    """Validate the structure of an optional record ``check`` block."""
    if not isinstance(check, Mapping):
        raise EosmatError(f"{location}.check must be an object")
    unknown = set(check) - set(CHECK_KEYS)
    if unknown:
        raise EosmatError(f"{location}.check has unknown fields {sorted(unknown)}")
    if "reproduction" in check and check["reproduction"] not in REPRODUCTION_VALUES:
        raise EosmatError(
            f"{location}.check.reproduction must be one of "
            f"{', '.join(REPRODUCTION_VALUES)}"
        )
    if "source" in check and check["source"] not in SOURCE_STATUSES:
        raise EosmatError(f"{location}.check.source is invalid")
    if check.get("warning") is not None and (
        not isinstance(check["warning"], str) or not check["warning"].strip()
    ):
        raise EosmatError(f"{location}.check.warning must be null or text")
    for key in ("finding", "source_checked_on", "audit"):
        if key in check and (not isinstance(check[key], str) or not check[key]):
            raise EosmatError(f"{location}.check.{key} must be a non-empty string")


def validate_check_consistency(record: Mapping[str, Any], location: str) -> None:
    """A ``check`` headline must agree with inline evidence it summarizes."""
    check = record.get("check")
    if not isinstance(check, Mapping):
        return
    evidence = _evidence(record)
    assessment = _section(_section(evidence, "reproduction"), "assessment")
    if "reproduction" in check and "status" in assessment:
        expected = REPRODUCTION_HEADLINES.get(assessment["status"])
        if check["reproduction"] != expected:
            raise EosmatError(
                f"{location}.check.reproduction {check['reproduction']!r} "
                f"contradicts assessment status {assessment['status']!r}"
            )
    validation = _section(evidence, "scientific_validation")
    if "source" in check and "status" in validation:
        if check["source"] != validation["status"]:
            raise EosmatError(
                f"{location}.check.source contradicts scientific_validation.status"
            )


def audit_is_external(record: Mapping[str, Any]) -> bool:
    """Whether a record names a sidecar entry and carries no inline audit."""
    check = record.get("check")
    return (
        isinstance(check, Mapping)
        and isinstance(check.get("audit"), str)
        and not any(key in record for key in AUDIT_RECORD_KEYS)
    )


def _without(values: list[Any], removed: tuple[str, ...]) -> list[Any]:
    return [value for value in values if value not in removed]


def remove_generic_statements(record: Mapping[str, Any]) -> dict[str, Any]:
    """Return a record copy without the documented generic statements.

    A qualification list emptied by this removal and a removed note are
    omitted. Record-specific statements are kept verbatim and in order.
    """
    result = copy.deepcopy(dict(record))
    evidence = result.get("evidence")
    if not isinstance(evidence, dict):
        return result
    assessment = _mutable_section(evidence.get("reproduction"), "assessment")
    qualifications = assessment.get("qualifications")
    if isinstance(qualifications, list):
        kept = _without(qualifications, GENERIC_QUALIFICATIONS)
        if kept != qualifications:
            if kept:
                assessment["qualifications"] = kept
            else:
                del assessment["qualifications"]
    validation = evidence.get("scientific_validation")
    if isinstance(validation, dict) and validation.get("note") in (
        GENERIC_VALIDATION_NOTES
    ):
        del validation["note"]
    return result


def compact_entry(
    entry: Mapping[str, Any], material: str, record: str
) -> dict[str, Any]:
    """Drop report links that equal the identifier-derived values."""
    result = copy.deepcopy(dict(entry))
    assessment = _mutable_section(
        _mutable_section(result, "evidence").get("reproduction"), "assessment"
    )
    if assessment:
        for key, value in derived_report_links(material, record).items():
            if key not in assessment:
                raise EosmatError(
                    f"Record {record!r} has no assessment {key}; a bundled sidecar "
                    "regenerates that link and cannot represent its absence"
                )
            if assessment[key] == value:
                del assessment[key]
    return result


def expand_entry(
    entry: Mapping[str, Any], material: str, record: str
) -> dict[str, Any]:
    """Regenerate the identifier-derived report links of a bundled entry."""
    result = copy.deepcopy(dict(entry))
    assessment = _mutable_section(
        _mutable_section(result, "evidence").get("reproduction"), "assessment"
    )
    if assessment:
        for key, value in derived_report_links(material, record).items():
            assessment.setdefault(key, value)
    return result


def split_document(
    document: Mapping[str, Any], *, add_check: bool = True
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Split a complete native document into a record file and a sidecar.

    Generic statements are removed and derivable report links are omitted.
    Each record gets a ``check`` block derived from its evidence (an existing
    block is replaced only by the identical derivation plus the audit link).
    """
    if document.get("format_version") != 5:
        raise EosmatError("Audit sidecars are defined for EOSMAT v5 only")
    material = str(document["identifier"])
    records_document = copy.deepcopy(dict(document))
    sidecar: dict[str, Any] = {
        "format": AUDIT_FORMAT,
        "format_version": 5,
        "material": material,
        "records": {},
    }
    for record in records_document["eos_records"]:
        identifier = str(record["identifier"])
        if add_check:
            existing = record.get("check")
            derived = check_block(
                record,
                document,
                audit=bundled_audit_reference(material, identifier),
            )
            if existing is not None and {
                k: v for k, v in existing.items() if k != "audit"
            } != {k: v for k, v in derived.items() if k != "audit"}:
                raise EosmatError(
                    f"Record {identifier!r} check block differs from its evidence"
                )
            record["check"] = derived
        cleaned = remove_generic_statements(record)
        entry = {key: cleaned.pop(key) for key in AUDIT_RECORD_KEYS if key in cleaned}
        record.clear()
        record.update(cleaned)
        if entry:
            sidecar["records"][identifier] = compact_entry(entry, material, identifier)
    return records_document, sidecar


def validate_sidecar(sidecar: Any, material: str) -> None:
    if not isinstance(sidecar, Mapping):
        raise EosmatError("An EOSMAT audit sidecar must be a JSON object")
    if sidecar.get("format") != AUDIT_FORMAT or sidecar.get("format_version") != 5:
        raise EosmatError(
            f"Audit sidecar for {material!r} must declare format {AUDIT_FORMAT!r} "
            "and format_version 5"
        )
    if sidecar.get("material") != material:
        raise EosmatError(
            f"Audit sidecar names material {sidecar.get('material')!r}, "
            f"expected {material!r}"
        )
    records = sidecar.get("records")
    if not isinstance(records, Mapping):
        raise EosmatError("Audit sidecar records must be an object")
    for identifier, entry in records.items():
        if not isinstance(entry, Mapping) or set(entry) - set(AUDIT_RECORD_KEYS):
            raise EosmatError(
                f"Audit sidecar entry {identifier!r} may contain only "
                f"{', '.join(AUDIT_RECORD_KEYS)}"
            )


def merge_audit(
    document: Mapping[str, Any], sidecar: Mapping[str, Any]
) -> dict[str, Any]:
    """Return the complete inline document: record file plus sidecar entries."""
    material = str(document.get("identifier"))
    validate_sidecar(sidecar, material)
    entries = sidecar["records"]
    result = copy.deepcopy(dict(document))
    records = result.get("eos_records")
    if not isinstance(records, list):
        raise EosmatError("eos_records must be an array")
    seen = set()
    for record in records:
        identifier = str(record.get("identifier"))
        if identifier not in entries:
            continue
        seen.add(identifier)
        collision = set(entries[identifier]) & set(record)
        if collision:
            raise EosmatError(
                f"Record {identifier!r} holds {sorted(collision)} both inline and "
                "in its audit sidecar"
            )
        record.update(expand_entry(entries[identifier], material, identifier))
    unknown = set(entries) - seen
    if unknown:
        raise EosmatError(
            f"Audit sidecar for {material!r} names unknown records {sorted(unknown)}"
        )
    return result


def audit_entry(record: Mapping[str, Any]) -> dict[str, Any]:
    """The audit-only part of one complete native record."""
    return {
        key: copy.deepcopy(record[key]) for key in AUDIT_RECORD_KEYS if key in record
    }


@lru_cache(maxsize=None)
def _bundled_audit_text(material: str) -> str:
    resource = resources.files(AUDIT_PACKAGE).joinpath(audit_file_name(material))
    try:
        return resource.read_text(encoding="utf-8")
    except FileNotFoundError as error:
        raise EosmatError(
            f"The audit sidecar for bundled material {material!r} is missing "
            f"(peritheos/data/audits/{audit_file_name(material)}); the record file "
            "alone does not contain its evidence. Reinstall Peritheos.",
            code="eosmat.audit_missing",
            operation="load_audit",
            context={"material": material},
        ) from error


def load_bundled_audit(material: str) -> dict[str, Any]:
    """Load one bundled sidecar; raise ``EosmatError`` if it is missing."""
    sidecar = json.loads(_bundled_audit_text(material))
    validate_sidecar(sidecar, material)
    return sidecar


def load_referenced_audits(document: Mapping[str, Any], base: Path) -> dict[str, Any]:
    """Merge the sidecars named by ``check.audit`` references, relative to ``base``."""
    result = copy.deepcopy(dict(document))
    files: dict[str, Any] = {}
    for record in result.get("eos_records", []):
        if not isinstance(record, Mapping) or not audit_is_external(record):
            continue
        reference = str(record["check"]["audit"])
        relative, _, fragment = reference.partition("#")
        if fragment and fragment != record.get("identifier"):
            raise EosmatError(
                f"check.audit {reference!r} does not name record "
                f"{record.get('identifier')!r}"
            )
        if relative not in files:
            path = base / relative
            try:
                files[relative] = json.loads(path.read_text(encoding="utf-8"))
            except FileNotFoundError as error:
                raise EosmatError(
                    f"Audit sidecar {str(path)!r} referenced by record "
                    f"{record.get('identifier')!r} is missing",
                    code="eosmat.audit_missing",
                    operation="load_audit",
                    context={"path": str(path)},
                ) from error
    for sidecar in files.values():
        material = str(result.get("identifier"))
        validate_sidecar(sidecar, material)
        selected = {
            key: value
            for key, value in sidecar["records"].items()
            if any(
                isinstance(r, Mapping) and r.get("identifier") == key
                for r in result["eos_records"]
            )
        }
        result = merge_audit(result, {**sidecar, "records": selected})
    return result


__all__ = [
    "AUDIT_RECORD_KEYS",
    "GENERIC_QUALIFICATIONS",
    "GENERIC_VALIDATION_NOTES",
    "REPRODUCTION_HEADLINES",
    "audit_entry",
    "check_block",
    "derived_report_links",
    "load_bundled_audit",
    "merge_audit",
    "remove_generic_statements",
    "split_document",
    "usage_warning",
]
