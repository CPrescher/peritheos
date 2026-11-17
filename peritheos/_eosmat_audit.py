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
the stored assessment status, ``check.compared`` a fixed mapping of its status
and scope, ``check.why`` is composed only from the stored coefficient
comparison table, ``check.source_check`` is a fixed renaming of the stored
source-validation status and ``check.warning`` is built only from stored coded
restrictions (see :func:`usage_warning`). The audit sidecar keeps the original
codes. Two curated texts are copied, never composed: ``check.source_note`` is
the per-record ``primary_source_check.record_note`` when one is written and the
verbatim ``primary_source_check.finding`` otherwise, and
``check.warning_reason`` copies ``scientific_validation.warning_reason``.
"""

from __future__ import annotations

import copy
import json
from collections.abc import Callable, Mapping
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
    "superseded_record_text",
)

#: Free-text record fields whose replaced original text a record keeps in its
#: ``superseded_record_text`` audit block: ``label``, ``notes``,
#: ``validity_notes`` (``validity.notes``) and ``thermal_notes``
#: (``thermal.notes``). ``null`` means the field was absent before.
SUPERSEDED_TEXT_FIELDS = ("label", "notes", "validity_notes", "thermal_notes")

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

#: What ``reproduction`` refers to, from the stored assessment status and scope.
#: ``published_fit``: a refit compared with the published coefficients
#: (scope ``parameter_agreement`` or ``author_workflow``); ``peritheos_refit``:
#: a separately defined Peritheos refit rerun under its own inputs (scope
#: ``independent_refit``); ``source_equation``: the source equation or its
#: calculated checkpoints were rebuilt (status ``source_reconstruction``, or
#: scope ``source_equation`` with an assessed status); ``none``: nothing was
#: compared (status ``not_refittable`` or ``not_assessed``, or no assessment).
COMPARED_VALUES = ("published_fit", "peritheos_refit", "source_equation", "none")
_COMPARED_BY_SCOPE = {
    "parameter_agreement": "published_fit",
    "author_workflow": "published_fit",
    "independent_refit": "peritheos_refit",
    "source_equation": "source_equation",
}
_NOT_COMPARED_STATUSES = ("not_refittable", "not_assessed")

#: Fixed renaming of ``scientific_validation.status`` for ``check.source_check``.
#: The sidecar and compatibility exports keep the original codes.
SOURCE_CHECK_VALUES = {
    "primary_source_validated": "verified",
    "pending_primary_source_check": "pending",
    "not_reproduced": "unresolved",
    "deferred": "deferred",
}
SOURCE_STATUSES = tuple(SOURCE_CHECK_VALUES)
SOURCE_CHECK_CODES = {value: code for code, value in SOURCE_CHECK_VALUES.items()}
CHECK_KEYS = (
    "reproduction",
    "compared",
    "why",
    "warning",
    "warning_reason",
    "source_check",
    "source_checked_on",
    "source_note",
    "audit",
)
#: Earlier names of v5 ``check`` fields. They are rejected, never reinterpreted.
RENAMED_CHECK_KEYS = {"finding": "source_note", "source": "source_check"}
WHY_MAX_LENGTH = 160
#: Limits of the two curated sidecar texts copied into ``check``:
#: ``scientific_validation.primary_source_check.record_note`` (shown as
#: ``check.source_note``) and ``scientific_validation.warning_reason`` (shown as
#: ``check.warning_reason``).
RECORD_NOTE_MAX_LENGTH = 300
WARNING_REASON_MAX_LENGTH = 140
#: The Unicode ``White_Space`` characters. A curated text (``source_note``,
#: ``warning_reason`` and their sidecar fields) made only of these is empty;
#: the Rust reader trims the same set (``str::trim``).
WHITE_SPACE = (
    "\t\n\x0b\x0c\r \x85\xa0\u1680\u2000\u2001\u2002\u2003\u2004\u2005"
    "\u2006\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000"
)


def _blank(value: Any) -> bool:
    """A value that is not a string or holds only ``WHITE_SPACE``."""
    return not isinstance(value, str) or not value.strip(WHITE_SPACE)


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
_EXPLICIT_SELECTION_WARNING = "Never chosen automatically; select it explicitly."


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


def _assessment(record: Mapping[str, Any]) -> Mapping[str, Any]:
    return _section(_section(_evidence(record), "reproduction"), "assessment")


def reproduction_headline(record: Mapping[str, Any]) -> str:
    """Map the stored assessment status; anything else is ``untested``."""
    return REPRODUCTION_HEADLINES.get(
        str(_assessment(record).get("status")), "untested"
    )


def compared_scope(record: Mapping[str, Any]) -> str:
    """Map the stored assessment status and scope (see :data:`COMPARED_VALUES`)."""
    assessment = _assessment(record)
    status = assessment.get("status")
    if status is None or status in _NOT_COMPARED_STATUSES:
        return "none"
    if status == "source_reconstruction":
        return "source_equation"
    return _COMPARED_BY_SCOPE.get(str(assessment.get("scope")), "none")


_WHY_UNITS = {"K0": " GPa", "K0_double_prime": " GPa^-1", "theta0": " K"}


def _parameter_name(name: str) -> str:
    name = name.removeprefix("rt_eos.")
    if name.endswith("_double_prime"):
        return name.removesuffix("_double_prime") + "''"
    if name.endswith("_prime"):
        return name.removesuffix("_prime") + "'"
    return name


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _printed(value: Any) -> str:
    """The stored number as JSON prints it, without a spurious trailing ``.0``."""
    if isinstance(value, float) and value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return json.dumps(value)


def _rounded_refit(refit: float, published: str) -> str:
    """Round to at most one digit more than the printed published value."""
    mantissa, _, exponent = published.lower().partition("e")
    if exponent:
        digits = len(mantissa.replace("-", "").replace(".", "").lstrip("0"))
        return f"{refit:.{digits}e}"
    decimals = len(mantissa.partition(".")[2])
    text = f"{refit:.{decimals + 1}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def _join(items: list[str]) -> str:
    """``a``, ``a and b``, ``a, b and c``."""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _failed_criteria(row: Mapping[str, Any]) -> list[str]:
    """Phrases for the stored criterion flags that are ``false`` in one row.

    ``similar`` is the audit's own agreement screen (the generic limits or the
    criterion recorded in the assessment basis), so it is reported as "beyond
    the agreement limit" without naming a limit.
    """
    failed = []
    if row.get("within_combined_2sigma") is False:
        failed.append("outside the combined error")
    if row.get("similar") is False:
        failed.append("beyond the agreement limit")
    if row.get("within_reported_error") is False:
        failed.append("outside the published error")
    return failed


def _relative_change(row: Mapping[str, Any]) -> float | None:
    refit, base = _number(row.get("refit")), _number(row.get("published"))
    if refit is None or base is None or base == 0.0:
        return None
    return (refit - base) / abs(base)


def _why_row(row: Mapping[str, Any], *, numbers: bool) -> str:
    name = _parameter_name(str(row["parameter"]))
    published, refit = row.get("published"), _number(row.get("refit"))
    relative = _relative_change(row)
    change = None if relative is None else f"{relative * 100:+.1f}%"
    if not numbers:
        return name if change is None else f"{name} ({change})"
    if refit is None or _number(published) is None:
        return name
    unit = _WHY_UNITS.get(str(row["parameter"]).removeprefix("rt_eos."), "")
    printed = _printed(published)
    text = f"{name}: refit {_rounded_refit(refit, printed)} vs published "
    text += f"{printed}{unit}"
    return text + (f" ({change})" if change else "")


#: At most this many failing coefficients are named; the rest are counted.
_WHY_NAMED = 3


def _failing_lines(failing: list[tuple[Mapping[str, Any], list[str]]]) -> list[str]:
    """Candidate texts for the failing rows, most detailed first."""
    hidden = max(0, len(failing) - _WHY_NAMED)
    if hidden:
        # Name the largest relative differences; rows without one come last.
        def size(item: tuple[Mapping[str, Any], list[str]]) -> float:
            change = _relative_change(item[0])
            return -1.0 if change is None else abs(change)

        failing = sorted(failing, key=size, reverse=True)[:_WHY_NAMED]
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for row, failed in failing:
        groups.setdefault(_join(failed), []).append(row)
    more = f"{hidden} more"

    def grouped(name: Callable[[Mapping[str, Any]], str]) -> str:
        if hidden and len(groups) == 1:
            ((failed, group),) = groups.items()
            return f"{_join([name(r) for r in group] + [more])}: {failed}"
        text = "; ".join(
            f"{_join([name(r) for r in group])}: {failed}"
            for failed, group in groups.items()
        )
        return text + (f"; and {more}" if hidden else "")

    candidates = []
    if not hidden:
        candidates.append(
            "; ".join(
                f"{_why_row(row, numbers=True)}, {_join(failed)}"
                for row, failed in failing
            )
        )
    candidates.append(grouped(lambda r: _why_row(r, numbers=False)))
    candidates.append(grouped(lambda r: _parameter_name(str(r["parameter"]))))
    return candidates


def reproduction_why(record: Mapping[str, Any]) -> str | None:
    """One generated line saying which coefficients kept a refit from parity.

    Only for ``close`` and ``disagrees``. Built from the stored
    ``assessment.parameters`` rows and their criterion flags
    (``within_combined_2sigma``, ``similar``, ``within_reported_error``);
    no prose is interpreted. Published values print as stored (without a
    trailing ``.0``); the refit is rounded to at most one digit more. More
    than three failing coefficients: the three largest relative differences
    are named and the rest counted. Returns ``None`` when the table supports
    no statement within :data:`WHY_MAX_LENGTH` characters.
    """
    headline = reproduction_headline(record)
    if headline not in ("close", "disagrees"):
        return None
    rows = _assessment(record).get("parameters")
    if not isinstance(rows, list):
        return None
    rows = [r for r in rows if isinstance(r, Mapping) and "parameter" in r]
    failing = [(row, _failed_criteria(row)) for row in rows if _failed_criteria(row)]
    if failing:
        for text in _failing_lines(failing):
            if len(text) + 1 <= WHY_MAX_LENGTH:
                return text + "."
        return None
    if headline != "close":
        return None
    untested = [
        _parameter_name(str(row["parameter"]))
        for row in rows
        if "within_combined_2sigma" in row and row["within_combined_2sigma"] is None
    ]
    if not untested:
        return None
    verb = "has" if len(untested) == 1 else "have"
    text = f"{_join(untested)} {verb} no error to test against the combined-error rule"
    if rows and all(row.get("similar") is True for row in rows):
        text += "; all within the agreement limit"
    text += "."
    return text if len(text) <= WHY_MAX_LENGTH else None


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


def _text(value: Any) -> str | None:
    return value if isinstance(value, str) and value else None


def source_note(record: Mapping[str, Any]) -> str | None:
    """``check.source_note``: the per-record note, else the verbatim finding.

    ``scientific_validation.primary_source_check.record_note`` is a curated
    statement of what was checked for this record; without one the check's
    ``finding`` is copied unchanged. Neither is ever composed or shortened.
    """
    source = _section(
        _section(_evidence(record), "scientific_validation"), "primary_source_check"
    )
    return _text(source.get("record_note")) or _text(source.get("finding"))


def warning_reason(record: Mapping[str, Any]) -> str | None:
    """The stored ``scientific_validation.warning_reason``, or ``None``."""
    validation = _section(_evidence(record), "scientific_validation")
    return _text(validation.get("warning_reason"))


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
        "compared": compared_scope(record),
    }
    why = reproduction_why(record)
    if why is not None:
        result["why"] = why
    result["warning"] = usage_warning(record, material)
    # Copied whenever stored: a reason without a warning is invalid and is
    # reported by validation, never silently dropped.
    reason = warning_reason(record)
    if reason is not None:
        result["warning_reason"] = reason
    if validation.get("status") in SOURCE_CHECK_VALUES:
        result["source_check"] = SOURCE_CHECK_VALUES[validation["status"]]
    if "audit_date" in validation:
        result["source_checked_on"] = validation["audit_date"]
    note = source_note(record)
    if note is not None:
        result["source_note"] = note
    if audit is not None:
        result["audit"] = audit
    return result


def headline(
    record: Mapping[str, Any], material: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Effective ``check`` fields of a record, stored or derived.

    Returns ``reproduction``, ``compared``, ``why``, ``warning``,
    ``warning_reason``, ``source_check``, ``source_checked_on`` and
    ``source_note``. A stored ``check`` value wins;
    otherwise the value is derived from inline evidence exactly as
    :func:`check_block` would. No assessment means ``untested``, compared
    ``none``.
    """
    check = record.get("check")
    check = check if isinstance(check, Mapping) else {}
    validation = _section(_evidence(record), "scientific_validation")
    status = validation.get("status")
    warning = (
        check["warning"] if "warning" in check else usage_warning(record, material)
    )
    return {
        "reproduction": check.get("reproduction") or reproduction_headline(record),
        "compared": check.get("compared") or compared_scope(record),
        "why": check["why"] if "why" in check else reproduction_why(record),
        "warning": warning,
        "warning_reason": (
            None
            if warning is None
            else check.get("warning_reason") or warning_reason(record)
        ),
        "source_check": check.get("source_check")
        or SOURCE_CHECK_VALUES.get(str(status)),
        "source_checked_on": check.get("source_checked_on")
        or validation.get("audit_date"),
        "source_note": check.get("source_note") or source_note(record),
    }


def validate_check(check: Any, location: str) -> None:
    """Validate the structure of an optional record ``check`` block."""
    if not isinstance(check, Mapping):
        raise EosmatError(f"{location}.check must be an object")
    for old, new in RENAMED_CHECK_KEYS.items():
        if old in check:
            raise EosmatError(
                f"{location}.check.{old} is no longer read; EOSMAT v5 names this "
                f"field check.{new}"
                + (
                    f" ({', '.join(f'{c} -> {v}' for c, v in SOURCE_CHECK_VALUES.items())})"
                    if old == "source"
                    else ""
                )
            )
    unknown = set(check) - set(CHECK_KEYS)
    if unknown:
        raise EosmatError(f"{location}.check has unknown fields {sorted(unknown)}")
    if "reproduction" in check and check["reproduction"] not in REPRODUCTION_VALUES:
        raise EosmatError(
            f"{location}.check.reproduction must be one of "
            f"{', '.join(REPRODUCTION_VALUES)}"
        )
    if "compared" in check and check["compared"] not in COMPARED_VALUES:
        raise EosmatError(
            f"{location}.check.compared must be one of {', '.join(COMPARED_VALUES)}"
        )
    if "reproduction" in check and "compared" in check:
        # Only an assessed comparison can give a result; ``untested`` has none
        # or only a rebuilt source equation.
        allowed = (
            ("none", "source_equation")
            if check["reproduction"] == "untested"
            else ("published_fit", "peritheos_refit", "source_equation")
        )
        if check["compared"] not in allowed:
            raise EosmatError(
                f"{location}.check.compared {check['compared']!r} cannot accompany "
                f"reproduction {check['reproduction']!r}"
            )
    if "source_check" in check and check["source_check"] not in SOURCE_CHECK_CODES:
        raise EosmatError(
            f"{location}.check.source_check must be one of "
            f"{', '.join(SOURCE_CHECK_CODES)}"
        )
    if check.get("warning") is not None and (
        not isinstance(check["warning"], str) or not check["warning"].strip()
    ):
        raise EosmatError(f"{location}.check.warning must be null or text")
    for key in ("why", "warning_reason", "source_note", "source_checked_on", "audit"):
        if key in check and (not isinstance(check[key], str) or not check[key]):
            raise EosmatError(f"{location}.check.{key} must be a non-empty string")
    for key in ("warning_reason", "source_note"):
        if key in check and _blank(check[key]):
            raise EosmatError(f"{location}.check.{key} must be a non-empty string")
    if "warning_reason" in check:
        if len(check["warning_reason"]) > WARNING_REASON_MAX_LENGTH:
            raise EosmatError(
                f"{location}.check.warning_reason must be at most "
                f"{WARNING_REASON_MAX_LENGTH} characters"
            )
        if check.get("warning") is None:
            raise EosmatError(
                f"{location}.check.warning_reason is only given with a warning"
            )
    if "why" in check:
        if len(check["why"]) > WHY_MAX_LENGTH:
            raise EosmatError(
                f"{location}.check.why must be at most {WHY_MAX_LENGTH} characters"
            )
        if check.get("reproduction", "close") not in ("close", "disagrees"):
            raise EosmatError(
                f"{location}.check.why is only given for reproduction close or "
                "disagrees"
            )


def validate_superseded_record_text(value: Any, location: str) -> None:
    """Validate a record's ``superseded_record_text`` audit block.

    The block keeps the original value of every free-text field replaced by a
    text revision (see :data:`SUPERSEDED_TEXT_FIELDS`) and the ISO date
    ``superseded_on``. It never changes how the record is read.
    """
    import datetime
    import re

    where = f"{location}.superseded_record_text"
    if not isinstance(value, Mapping):
        raise EosmatError(f"{where} must be an object")
    unknown = set(value) - {*SUPERSEDED_TEXT_FIELDS, "superseded_on"}
    if unknown:
        raise EosmatError(f"{where} has unknown fields {sorted(unknown)}")
    if not any(key in value for key in SUPERSEDED_TEXT_FIELDS):
        raise EosmatError(f"{where} names no superseded field")
    date = value.get("superseded_on")
    try:
        valid_date = (
            isinstance(date, str)
            and re.fullmatch(r"\d{4}-\d{2}-\d{2}", date) is not None
            and datetime.date.fromisoformat(date) is not None
        )
    except ValueError:
        valid_date = False
    if not valid_date:
        raise EosmatError(f"{where}.superseded_on must be an ISO date (YYYY-MM-DD)")
    for key in ("label", "notes", "thermal_notes"):
        if key in value and value[key] is not None and not isinstance(value[key], str):
            raise EosmatError(f"{where}.{key} must be a string or null")
    notes = value.get("validity_notes")
    if notes is not None and (
        not isinstance(notes, list) or not all(isinstance(n, str) for n in notes)
    ):
        raise EosmatError(f"{where}.validity_notes must be an array of strings or null")


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
    if "compared" in check and "status" in assessment:
        expected = compared_scope(record)
        if check["compared"] != expected:
            raise EosmatError(
                f"{location}.check.compared {check['compared']!r} contradicts "
                f"assessment status {assessment['status']!r} and scope "
                f"{assessment.get('scope')!r} (expected {expected!r})"
            )
    validation = _section(evidence, "scientific_validation")
    if "source_check" in check and "status" in validation:
        if check["source_check"] != SOURCE_CHECK_VALUES.get(validation["status"]):
            raise EosmatError(
                f"{location}.check.source_check {check['source_check']!r} "
                f"contradicts scientific_validation.status {validation['status']!r}"
            )
    # The curated texts are copies: with the evidence inline (or merged from
    # the sidecar) a stored copy must equal its derivation.
    if "source_note" in check and _section(validation, "primary_source_check"):
        if check["source_note"] != source_note(record):
            raise EosmatError(
                f"{location}.check.source_note differs from the evidence "
                "(scientific_validation.primary_source_check.record_note, else "
                "its finding)"
            )
    if "warning_reason" in check and validation:
        if check["warning_reason"] != warning_reason(record):
            raise EosmatError(
                f"{location}.check.warning_reason differs from the evidence "
                "(scientific_validation.warning_reason)"
            )


def validate_check_text(
    record: Mapping[str, Any],
    material: Mapping[str, Any] | None,
    location: str,
) -> None:
    """Validate the curated texts ``check`` copies from a record's evidence.

    ``scientific_validation.primary_source_check.record_note`` and
    ``scientific_validation.warning_reason`` are optional non-empty strings
    of at most :data:`RECORD_NOTE_MAX_LENGTH` and
    :data:`WARNING_REASON_MAX_LENGTH` characters. A ``warning_reason`` needs
    a warning: it is an error on a record whose derived ``check.warning`` is
    ``null``. Records without inline evidence are not checked here.
    """
    evidence = record.get("evidence")
    if not isinstance(evidence, Mapping):
        return
    validation = _section(evidence, "scientific_validation")
    where = f"{location}.evidence.scientific_validation"
    source = _section(validation, "primary_source_check")
    for container, key, limit, place in (
        (source, "record_note", RECORD_NOTE_MAX_LENGTH, "primary_source_check."),
        (validation, "warning_reason", WARNING_REASON_MAX_LENGTH, ""),
    ):
        if key not in container:
            continue
        value = container[key]
        if _blank(value):
            raise EosmatError(f"{where}.{place}{key} must be a non-empty string")
        if len(value) > limit:
            raise EosmatError(
                f"{where}.{place}{key} must be at most {limit} characters"
            )
    if "warning_reason" in validation and usage_warning(record, material) is None:
        raise EosmatError(
            f"{where}.warning_reason is given but the record has no warning "
            "(no coded restriction)"
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
    Per-record default flags become the material ``default_records`` map.
    """
    from peritheos._eosmat_v4 import collect_default_records

    if document.get("format_version") != 5:
        raise EosmatError("Audit sidecars are defined for EOSMAT v5 only")
    material = str(document["identifier"])
    records_document = copy.deepcopy(dict(document))
    collect_default_records(records_document)
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
        if "superseded_record_text" in entry:
            validate_superseded_record_text(
                entry["superseded_record_text"], f"audit sidecar entry {identifier!r}"
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
    "COMPARED_VALUES",
    "SOURCE_CHECK_VALUES",
    "GENERIC_QUALIFICATIONS",
    "GENERIC_VALIDATION_NOTES",
    "RECORD_NOTE_MAX_LENGTH",
    "REPRODUCTION_HEADLINES",
    "SUPERSEDED_TEXT_FIELDS",
    "audit_entry",
    "check_block",
    "compared_scope",
    "derived_report_links",
    "load_bundled_audit",
    "merge_audit",
    "remove_generic_statements",
    "reproduction_why",
    "source_note",
    "split_document",
    "usage_warning",
    "validate_check_text",
    "validate_superseded_record_text",
    "WARNING_REASON_MAX_LENGTH",
    "warning_reason",
]
