"""EOSMAT v5 record ``check`` headlines and per-material audit sidecars.

A bundled material is stored as two files:

``peritheos/data/materials/<id>.eosmat``
    The readable record file: material header, EOS records (model,
    coefficients, reference state, ranges, notes, default flags and a short
    ``check`` block), then datasets.
``peritheos/data/audits/<id>.audit.json``
    The audit trail of every record, keyed by record identifier: ``evidence``,
    the detailed ``pressure_calibration`` and the other audit-only blocks
    listed in :data:`AUDIT_RECORD_KEYS`, and the material-level
    ``superseded_material_text`` block (:data:`MATERIAL_AUDIT_KEY`).

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
from collections.abc import Mapping
from decimal import Decimal
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

#: The material-level audit block: the replaced original of the material's
#: top-level ``notes`` (``null``: absent before) and ``superseded_on``. It sits
#: at the top level of a sidecar (after ``material``) and of a complete
#: inline document.
MATERIAL_AUDIT_KEY = "superseded_material_text"
SUPERSEDED_MATERIAL_TEXT_FIELDS = ("notes",)

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


#: Units ``why`` prints, as the record text standard writes them. EOSMAT stores
#: every volume per conventional cell in Å³ and every pressure in GPa
#: (docs/units.md); other coefficient units come from
#: :func:`peritheos._overview.parameter_unit`.
_WHY_VOLUME_UNIT = "Å³"
_WHY_PRESSURE_UNIT = "GPa"
#: The criterion flags of an assessment row and what a ``false`` one means:
#: the refit and published values differ by more than this.
_WHY_CRITERIA = (
    ("within_combined_2sigma", "twice the combined error"),
    ("similar", "the agreement limit"),
    ("within_reported_error", "the published error"),
)
#: From this relative difference (in percent) on, ``why`` gives the refit and
#: published values instead of a percentage.
_WHY_LARGE_PERCENT = 99.95
#: The largest coefficients named when a line counts the failing ones.
_WHY_LARGEST = 3
#: With more failing coefficients than this, the line counts them.
_WHY_NAMED_AT_MOST = 4
#: A published value printed with more significant digits than this (a
#: computed conversion) is not repeated in full; the percentage is given.
_WHY_PRINTED_DIGITS = 6


def _parameter_name(name: str) -> str:
    """The coefficient as the record summary shows it (``K0'``, ``dK/dT``)."""
    from peritheos._overview import PARAMETER_DISPLAY_NAMES

    name = name.removeprefix("rt_eos.")
    return PARAMETER_DISPLAY_NAMES.get(name, name)


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
    """Round to at most one digit more than the printed published value.

    A published value in e-notation gives the refit one more significant
    digit. Otherwise the refit gets one more decimal, without trailing zeros,
    and is written in e-notation where the record summary would write it
    (an exponent below -4 or from 6 on; ``9.224e-05``, not ``0.00009224``).
    """
    mantissa, _, exponent = published.lower().partition("e")
    if exponent:
        digits = len(mantissa.replace("-", "").replace(".", "").lstrip("0"))
        return f"{refit:.{digits}e}"
    decimals = len(mantissa.partition(".")[2])
    text = f"{refit:.{decimals + 1}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    if text in ("0", "-0"):
        return "0"
    power = Decimal(text).adjusted()
    if power < -4 or power >= 6:
        digits = len(text.replace("-", "").replace(".", "").lstrip("0").rstrip("0"))
        return f"{float(text):.{max(digits - 1, 0)}e}"
    return text


def _join(items: list[str], word: str = "and") -> str:
    """``a``, ``a and b``, ``a, b and c``."""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + f" {word} " + items[-1]


def _why_unit(record: Mapping[str, Any], parameter: str) -> str:
    """The coefficient's unit with a leading space; empty when dimensionless.

    The owning model is the reference (``rt_eos.``) model unless the record's
    thermal parameters hold the name. Native records name the models in
    ``model`` and ``thermal.model``, v3 exports in ``eos.type`` and
    ``thermal.type``.
    """
    from peritheos._overview import parameter_unit

    eos = record.get("eos")
    model = eos.get("type") if isinstance(eos, Mapping) else record.get("model")
    thermal = record.get("thermal")
    name = parameter.removeprefix("rt_eos.")
    if not parameter.startswith("rt_eos.") and isinstance(thermal, Mapping):
        parameters = thermal.get("parameters")
        if isinstance(parameters, Mapping) and name in parameters:
            model = thermal["type"] if "type" in thermal else thermal.get("model")
    unit = parameter_unit(
        model if isinstance(model, str) else "",
        name,
        volume_unit=_WHY_VOLUME_UNIT,
        pressure_unit=_WHY_PRESSURE_UNIT,
    )
    return "" if unit is None or unit == "1" else f" {unit}"


def _is_fixed(record: Mapping[str, Any], parameter: str) -> bool:
    """Whether the record holds the coefficient fixed.

    A native record marks it ``"fixed": true`` in ``parameters`` (or
    ``thermal.parameters``), a v3 export lists it in ``fixed_parameters`` (or
    ``thermal.fixed_parameters``). The reference (``rt_eos.``) model owns the
    name unless the record's thermal parameters hold it.
    """
    name = parameter.removeprefix("rt_eos.")
    owner: Any = record
    thermal = record.get("thermal")
    if not parameter.startswith("rt_eos.") and isinstance(thermal, Mapping):
        parameters = thermal.get("parameters")
        if isinstance(parameters, Mapping) and name in parameters:
            owner = thermal
    if not isinstance(owner, Mapping):
        return False
    listed = owner.get("fixed_parameters")
    if isinstance(listed, list) and name in listed:
        return True
    parameters = owner.get("parameters")
    if owner is record and isinstance(record.get("eos"), Mapping):
        parameters = record["eos"].get("parameters")
    entry = parameters.get(name) if isinstance(parameters, Mapping) else None
    return isinstance(entry, Mapping) and entry.get("fixed") is True


def _failed_criteria(row: Mapping[str, Any]) -> list[str]:
    """What the refit exceeds, for each stored criterion flag that is ``false``.

    ``similar`` is the audit's own agreement screen (the generic limits or the
    criterion recorded in the assessment basis), so it is reported as "the
    agreement limit" without naming a limit.
    """
    return [phrase for flag, phrase in _WHY_CRITERIA if row.get(flag) is False]


def _relative_change(row: Mapping[str, Any]) -> float | None:
    refit, base = _number(row.get("refit")), _number(row.get("published"))
    if refit is None or base is None or base == 0.0:
        return None
    return (refit - base) / abs(base)


def _percent(row: Mapping[str, Any]) -> str | None:
    """The relative difference, unless missing or too large to read as one."""
    relative = _relative_change(row)
    if relative is None or abs(relative) * 100.0 >= _WHY_LARGE_PERCENT:
        return None
    return f"{relative * 100:+.1f}%"


def _significant_digits(printed: str) -> int:
    mantissa = printed.lower().partition("e")[0]
    return len(mantissa.replace("-", "").replace(".", "").lstrip("0"))


def _values(row: Mapping[str, Any]) -> tuple[str, str] | None:
    """Rounded refit and published value as the summary prints it."""
    from peritheos._overview import format_number

    refit, published = _number(row.get("refit")), row.get("published")
    if refit is None or _number(published) is None:
        return None
    shown = format_number(published)
    return _rounded_refit(refit, shown), shown


class _Item:
    """One failing coefficient in the forms ``why`` can print it.

    ``full`` gives both values, the unit and the percentage; ``short`` the
    percentage alone, or both values when the difference is 100% or more;
    ``name`` the coefficient alone.
    """

    def __init__(self, record: Mapping[str, Any], row: Mapping[str, Any]) -> None:
        parameter = str(row["parameter"])
        self.name = _parameter_name(parameter)
        self.failed = _failed_criteria(row)
        relative = _relative_change(row)
        self.size = -1.0 if relative is None else abs(relative)
        percent, values = _percent(row), _values(row)
        if values is None:
            self.full = self.short = self.name
            return
        refit, published = values
        if percent is None:
            self.short = f"{self.name} ({refit} vs {published})"
        else:
            self.short = f"{self.name} ({percent})"
        stored = _significant_digits(_printed(row["published"]))
        if percent is not None and stored > _WHY_PRINTED_DIGITS:
            self.full = self.short
            return
        self.full = f"{self.name} {refit} vs {published}{_why_unit(record, parameter)}"
        if percent is not None:
            self.full += f" ({percent})"


def _union_criteria(items: list[_Item]) -> str:
    """Every criterion some item fails, once each, joined by "or"."""
    failed = [
        phrase
        for _, phrase in _WHY_CRITERIA
        if any(phrase in item.failed for item in items)
    ]
    return _join(failed, "or")


def _verb(names: list[str], one: str, several: str) -> str:
    return one if len(names) == 1 else several


def _grouped(items: list[_Item], texts: list[str], *, union: bool = False) -> str:
    """Join ``texts`` (one per item) by what each item's difference exceeds.

    Each criterion gets one clause naming the items that fail it; criteria
    failed by the same items share a clause. The clause with the most items
    comes first (ties in order of appearance, then of the criteria). An item
    named in an earlier clause is not repeated: "K0 also by more than ...".
    With ``union`` every item is listed in one clause that names each failed
    criterion once, joined by "or".
    """
    if union:
        clauses = [
            f"{_join(texts)} {_verb(texts, 'differs', 'differ')} by more than "
            f"{_union_criteria(items)}"
        ]
    else:
        groups: list[tuple[list[str], list[int]]] = []
        for _, phrase in _WHY_CRITERIA:
            members = [i for i, item in enumerate(items) if phrase in item.failed]
            if not members:
                continue
            for phrases, shared in groups:
                if shared == members:
                    phrases.append(phrase)
                    break
            else:
                groups.append(([phrase], members))
        groups.sort(key=lambda group: (-len(group[1]), group[1][0]))
        named: set[int] = set()
        clauses = []
        for phrases, members in groups:
            criterion = _join(phrases)
            fresh = [texts[i] for i in members if i not in named]
            seen = [items[i].name for i in members if i in named]
            if not fresh:
                clauses.append(f"{_join(seen)} also by more than {criterion}")
            else:
                clause = (
                    f"{_join(fresh)} {_verb(fresh, 'differs', 'differ')} by more "
                    f"than {criterion}"
                )
                if seen:
                    clause += f", as {_verb(seen, 'does', 'do')} {_join(seen)}"
                clauses.append(clause)
            named.update(members)
    text = "; ".join(clauses)
    if all(text == item.name for item, text in zip(items, texts)):
        return text
    return f"Refit vs published: {text}"


def _counted_lines(items: list[_Item], ordered: list[_Item]) -> list[str]:
    """The failing coefficients counted, with the largest named.

    The largest keep their values or percentage while the line fits, the
    smaller of them first giving way to their name.
    """
    count = f"{len(items)} coefficients differ by more than {_union_criteria(items)}"
    largest = ordered[:_WHY_LARGEST]
    candidates = []
    for shown in range(len(largest), -1, -1):
        texts = [
            item.short if position < shown else item.name
            for position, item in enumerate(largest)
        ]
        line = f"{count}; largest: {_join(texts)}"
        named = all(text == item.name for item, text in zip(largest, texts))
        candidates.append(line if named else f"Refit vs published: {line}")
    return candidates


def _failing_lines(items: list[_Item]) -> list[str]:
    """Candidate texts for the failing coefficients, most detailed first.

    Up to :data:`_WHY_NAMED_AT_MOST` coefficients: every coefficient with
    both values, its unit and its percentage; every one with its percentage;
    when the coefficients fail different criteria, every one with its
    percentage in one clause; the largest differences with their percentage
    and the others by name; names only. More coefficients, or when even the
    names do not fit, are counted and the :data:`_WHY_LARGEST` largest named.
    """
    # Largest first, so the coefficients that keep a percentage lead.
    ordered = sorted(items, key=lambda item: item.size, reverse=True)
    if len(items) > _WHY_NAMED_AT_MOST:
        return _counted_lines(items, ordered)
    mixed = len({_join(item.failed) for item in items}) > 1
    candidates = [
        _grouped(items, [item.full for item in items]),
        _grouped(items, [item.short for item in items]),
    ]
    if mixed:
        candidates.append(_grouped(items, [item.short for item in items], union=True))
    for shown in range(len(items) - 1, -1, -1):
        texts = [
            item.short if position < shown else item.name
            for position, item in enumerate(ordered)
        ]
        candidates.append(_grouped(ordered, texts))
        if mixed:
            candidates.append(_grouped(ordered, texts, union=True))
    if len(items) > _WHY_LARGEST:
        candidates.extend(_counted_lines(items, ordered))
    return candidates


def _missing_error(row: Mapping[str, Any]) -> str:
    """Which error the row lacks: ``published``, ``refit``, ``both`` or ``unknown``.

    An error counts as missing only when the row stores it as ``null``
    (``published_error``/``published_sigma``, ``refit_error``/``refit_sigma``).
    """

    def missing(side: str) -> bool:
        for key in (f"{side}_error", f"{side}_sigma"):
            if key in row:
                return row[key] is None
        return False

    published, refit = missing("published"), missing("refit")
    if published and refit:
        return "both"
    if published:
        return "published"
    return "refit" if refit else "unknown"


#: Why a coefficient was not tested against the combined error, by the missing
#: error: a sentence opening when one reason covers every such coefficient,
#: and what the coefficient has (or lacks) to list several reasons.
_UNTESTED_REASONS = {
    "published": ("the paper gives no error for", "no published error"),
    "refit": ("the refit gives no error for", "no refit error"),
    "both": (
        "neither the paper nor the refit gives an error for",
        "no published or refit error",
    ),
    "unknown": (
        "no error comparison is recorded for",
        "no recorded error comparison",
    ),
}


def _untested_lines(rows: list[Mapping[str, Any]], within: bool) -> list[str]:
    """Candidate texts for coefficients without a combined-error result.

    One reason: "the paper gives no error for V0 and K0, so they could not
    be tested ...", or shorter "no published error for V0 and K0, so ...".
    Several: "K0 has no published or refit error and alpha_KT no refit
    error, so neither could be tested ...".
    """
    reasons: dict[str, list[str]] = {}
    for row in rows:
        reasons.setdefault(_missing_error(row), []).append(
            _parameter_name(str(row["parameter"]))
        )
    head = "All coefficients are within the agreement limit; " if within else ""
    if len(reasons) == 1:
        ((reason, names),) = reasons.items()
        opening, lacking = _UNTESTED_REASONS[reason]
        pronoun = "it" if len(rows) == 1 else "they"
        result = f"{pronoun} could not be tested"
        bodies = [f"{opening} {_join(names)}", f"{lacking} for {_join(names)}"]
    else:
        # Starts with a coefficient name, which keeps its case.
        listed = [
            f"{_join(names)} "
            + (f"{_verb(names, 'has', 'have')} " if position == 0 else "")
            + _UNTESTED_REASONS[reason][1]
            for position, (reason, names) in enumerate(reasons.items())
        ]
        result = f"{'neither' if len(rows) == 2 else 'none'} could be tested"
        bodies = [_join(listed)]
    tail = f", so {result} against twice the combined error"
    texts = [head + body + tail for body in bodies]
    if head:
        texts.append(bodies[-1] + tail)
    if len(reasons) == 1:
        texts = [text[0].upper() + text[1:] for text in texts]
    return texts


def reproduction_why(record: Mapping[str, Any]) -> str | None:
    """One generated line saying why a refit missed parity, in plain words.

    Only for ``close`` and ``disagrees``. Built from the stored
    ``assessment.parameters`` rows and their criterion flags
    (``within_combined_2sigma``, ``similar``, ``within_reported_error``) and
    the stored errors; no prose is interpreted.

    With failing coefficients, each is named with its refit and published
    value, its unit and the relative difference, grouped by what the
    difference exceeds. A difference of 100% or more is shown by the two
    values alone. Coefficient names and published values print as the record
    summary prints them; the refit is rounded to at most one digit more.
    When that is too long the line falls back to percentages, then to names
    (the largest differences keep their percentage). More than four failing
    coefficients, or names that do not fit, are counted with the three
    largest named. A ``close`` record without failing coefficients names the
    coefficients that could not be tested against the combined error and
    which error is missing; a coefficient the record holds fixed is not a
    candidate for that test and is left out. Returns ``None`` when the table
    supports no statement within :data:`WHY_MAX_LENGTH` characters.
    """
    headline = reproduction_headline(record)
    if headline not in ("close", "disagrees"):
        return None
    rows = _assessment(record).get("parameters")
    if not isinstance(rows, list):
        return None
    rows = [r for r in rows if isinstance(r, Mapping) and "parameter" in r]
    items = [_Item(record, row) for row in rows if _failed_criteria(row)]
    if items:
        candidates = _failing_lines(items)
    elif headline == "close":
        untested = [
            row
            for row in rows
            if "within_combined_2sigma" in row
            and row["within_combined_2sigma"] is None
            and not _is_fixed(record, str(row["parameter"]))
        ]
        if not untested:
            return None
        within = bool(rows) and all(row.get("similar") is True for row in rows)
        candidates = _untested_lines(untested, within)
    else:
        return None
    for text in candidates:
        if len(text) + 1 <= WHY_MAX_LENGTH:
            return text + "."
    return None


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
    where = f"{location}.superseded_record_text"
    _validate_superseded_block(value, where, SUPERSEDED_TEXT_FIELDS)
    for key in ("label", "notes", "thermal_notes"):
        if key in value and value[key] is not None and not isinstance(value[key], str):
            raise EosmatError(f"{where}.{key} must be a string or null")
    notes = value.get("validity_notes")
    if notes is not None and (
        not isinstance(notes, list) or not all(isinstance(n, str) for n in notes)
    ):
        raise EosmatError(f"{where}.validity_notes must be an array of strings or null")


def validate_superseded_material_text(value: Any, location: str = "material") -> None:
    """Validate a material's ``superseded_material_text`` audit block.

    The block keeps the original top-level ``notes`` of a material replaced by
    a text revision (``null``: absent before) and the ISO date
    ``superseded_on``. It never changes how the material is read.
    """
    where = f"{location}.{MATERIAL_AUDIT_KEY}"
    _validate_superseded_block(value, where, SUPERSEDED_MATERIAL_TEXT_FIELDS)
    notes = value.get("notes")
    if notes is not None and not isinstance(notes, str):
        raise EosmatError(f"{where}.notes must be a string or null")


def _validate_superseded_block(value: Any, where: str, fields: tuple[str, ...]) -> None:
    """Keys and ``superseded_on`` of a ``superseded_*_text`` block."""
    import datetime
    import re

    if not isinstance(value, Mapping):
        raise EosmatError(f"{where} must be an object")
    unknown = set(value) - {*fields, "superseded_on"}
    if unknown:
        raise EosmatError(f"{where} has unknown fields {sorted(unknown)}")
    if not any(key in value for key in fields):
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
    A material-level ``superseded_material_text`` block moves to the sidecar.
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
    }
    if MATERIAL_AUDIT_KEY in records_document:
        sidecar[MATERIAL_AUDIT_KEY] = records_document.pop(MATERIAL_AUDIT_KEY)
    sidecar["records"] = {}
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
    if MATERIAL_AUDIT_KEY in sidecar:
        validate_superseded_material_text(
            sidecar[MATERIAL_AUDIT_KEY], f"audit sidecar for {material!r}"
        )


def merge_audit(
    document: Mapping[str, Any], sidecar: Mapping[str, Any]
) -> dict[str, Any]:
    """Return the complete inline document: record file plus sidecar entries.

    A sidecar's ``superseded_material_text`` block becomes a top-level key of
    the document.
    """
    material = str(document.get("identifier"))
    validate_sidecar(sidecar, material)
    entries = sidecar["records"]
    result = copy.deepcopy(dict(document))
    if MATERIAL_AUDIT_KEY in sidecar:
        if MATERIAL_AUDIT_KEY in result:
            raise EosmatError(
                f"Material {material!r} holds {MATERIAL_AUDIT_KEY} both inline and "
                "in its audit sidecar"
            )
        result[MATERIAL_AUDIT_KEY] = copy.deepcopy(sidecar[MATERIAL_AUDIT_KEY])
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
    "MATERIAL_AUDIT_KEY",
    "RECORD_NOTE_MAX_LENGTH",
    "REPRODUCTION_HEADLINES",
    "SUPERSEDED_MATERIAL_TEXT_FIELDS",
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
    "validate_superseded_material_text",
    "validate_superseded_record_text",
    "WARNING_REASON_MAX_LENGTH",
    "warning_reason",
]
