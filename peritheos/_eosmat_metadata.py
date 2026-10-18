"""Local metadata sharing and explicit reconstruction assessments for EOSMAT v4.

Only record evidence, citations and pressure-calibration metadata are resolved.
Coefficients, observations, array order and extension fields outside those
locations are never interpreted as references.
"""

from __future__ import annotations

import copy
import json
import math
import re
from collections import Counter
from collections.abc import Mapping
from typing import Any

from peritheos.errors import EosmatError

METADATA_FIELDS = ("reference", "evidence", "pressure_calibration")
STATUSES = (
    "parity",
    "similar",
    "not_reproduced",
    "not_refittable",
    "source_reconstruction",
    "not_assessed",
)
SCOPES = (
    "parameter_agreement",
    "author_workflow",
    "independent_refit",
    "source_equation",
)


def validate_assessment(evidence: Any) -> None:
    """Validate the opt-in standard; retain historical reproduction metadata."""
    if not isinstance(evidence, Mapping):
        return
    reproduction = evidence.get("reproduction")
    if not isinstance(reproduction, Mapping) or "assessment" not in reproduction:
        return
    assessment = reproduction["assessment"]
    where = "evidence.reproduction.assessment"
    if not isinstance(assessment, Mapping):
        raise EosmatError(f"{where} must be an object")
    if assessment.get("status") not in STATUSES:
        raise EosmatError(f"{where}.status must be one of {', '.join(STATUSES)}")
    if assessment.get("scope") not in SCOPES:
        raise EosmatError(f"{where}.scope must be one of {', '.join(SCOPES)}")
    for key in (
        "basis",
        "report",
        "report_url",
        "report_path",
        "source_status",
        "basis_code",
    ):
        if key in assessment and (
            not isinstance(assessment[key], str) or not assessment[key].strip()
        ):
            raise EosmatError(f"{where}.{key} must be a nonempty string")
    if assessment["status"] != "not_assessed" and not assessment.get("basis"):
        raise EosmatError(f"{where}.basis is required for an assessed result")
    qualifications = assessment.get("qualifications", [])
    if not isinstance(qualifications, list) or any(
        not isinstance(v, str) or not v.strip() for v in qualifications
    ):
        raise EosmatError(
            f"{where}.qualifications must be an array of nonempty strings"
        )
    if "fit" in assessment and not isinstance(assessment["fit"], Mapping):
        raise EosmatError(f"{where}.fit must be an object")
    if "parameters" in assessment:
        rows = assessment["parameters"]
        if not isinstance(rows, list):
            raise EosmatError(f"{where}.parameters must be an array")
        names: set[str] = set()
        for row in rows:
            if not isinstance(row, Mapping):
                raise EosmatError(f"{where}.parameters entries must be objects")
            name = row.get("parameter")
            if not isinstance(name, str) or not name.strip() or name in names:
                raise EosmatError(f"{where}.parameters requires unique parameter names")
            names.add(name)
            for key in ("published", "refit"):
                if key not in row:
                    raise EosmatError(f"{where}.parameters requires {key}")
            for key in (
                "published",
                "refit",
                "published_error",
                "refit_error",
                "published_sigma",
                "refit_sigma",
            ):
                if key not in row or row[key] is None:
                    continue
                value = row[key]
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or (isinstance(value, float) and not math.isfinite(value))
                ):
                    raise EosmatError(
                        f"{where}.parameters.{key} must be numeric or null"
                    )
                if key.endswith(("_error", "_sigma")) and value < 0:
                    raise EosmatError(f"{where}.parameters.{key} must be nonnegative")


def expand_metadata(document: Mapping[str, Any]) -> dict[str, Any]:
    """Return an inline v4 copy; resolve local metadata only, with no overrides."""
    result = copy.deepcopy(dict(document))
    datasets = result.get("datasets", [])
    for dataset in datasets if isinstance(datasets, list) else []:
        if isinstance(dataset, Mapping):
            citation = dataset.get("reference")
            if (
                isinstance(citation, Mapping)
                and "shared_ref" in citation
                and (document.get("shared_metadata") or len(citation) == 1)
            ):
                raise EosmatError(
                    "Dataset references must remain inline; shared_ref is metadata-only"
                )
    registry = result.pop("shared_metadata", {})
    had_registry = "shared_metadata" in document
    if not isinstance(registry, Mapping):
        raise EosmatError("shared_metadata must be an object")
    for name, value in registry.items():
        if (
            not isinstance(name, str)
            or not re.fullmatch(r"[a-z0-9]+(?:_[a-z0-9]+)*", name)
            or not isinstance(value, Mapping)
        ):
            raise EosmatError(
                "shared_metadata requires lower-snake-case names and object values"
            )

    nodes = 0

    def resolve(value: Any, stack: tuple[str, ...] = ()) -> Any:
        nonlocal nodes
        nodes += 1
        if registry and nodes > 100000:
            raise EosmatError("shared_metadata expansion exceeds 100000 metadata nodes")
        if len(stack) > 32:
            raise EosmatError("shared_metadata reference nesting exceeds 32")
        if isinstance(value, Mapping):
            if "shared_ref" in value and (registry or set(value) == {"shared_ref"}):
                if set(value) != {"shared_ref"}:
                    raise EosmatError(
                        "shared_ref cannot have sibling fields or overrides"
                    )
                name = value["shared_ref"]
                if not isinstance(name, str) or name not in registry:
                    raise EosmatError(f"Unknown shared_metadata reference {name!r}")
                if name in stack:
                    raise EosmatError(f"Cyclic shared_metadata reference {name!r}")
                return resolve(registry[name], (*stack, name))
            return {key: resolve(item, stack) for key, item in value.items()}
        if isinstance(value, list):
            return [resolve(item, stack) for item in value]
        return value

    # Check even unused definitions, so a saved pool cannot hide broken links.
    definitions = {name: resolve(value, (name,)) for name, value in registry.items()}
    if had_registry:
        # Even unused definitions are scientific information. Inline/legacy
        # exports retain the pool as inert extension metadata, never drop it.
        result["shared_metadata"] = definitions
    records = result.get("eos_records", [])
    for record in records if isinstance(records, list) else []:
        if not isinstance(record, dict):
            continue  # The ordinary structural validator diagnoses this.
        for key in METADATA_FIELDS:
            if key in record:
                record[key] = resolve(record[key])
        validate_assessment(record.get("evidence"))
    return result


def share_metadata(document: Mapping[str, Any]) -> dict[str, Any]:
    """Share only exactly repeated, substantial metadata objects in one file."""
    result = expand_metadata(document)
    counts: Counter[str] = Counter()

    def identity(value: Any) -> str:
        return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)

    def count(value: Any) -> None:
        if isinstance(value, dict):
            token = identity(value)
            if len(token) >= 160:
                counts[token] += 1
            for item in value.values():
                count(item)
        elif isinstance(value, list):
            for item in value:
                count(item)

    for record in result.get("eos_records", []):
        for key in METADATA_FIELDS:
            if key in record:
                count(record[key])
    registry: dict[str, Any] = copy.deepcopy(result.get("shared_metadata", {}))
    names: dict[str, str] = {}
    for name, value in registry.items():
        names.setdefault(identity(value), name)

    def replace(value: Any, context: str) -> Any:
        if isinstance(value, dict):
            token = identity(value)
            if counts[token] > 1:
                if token not in names:
                    authors = value.get("authors")
                    hint = (
                        authors[0]
                        if isinstance(authors, list)
                        and authors
                        and isinstance(authors[0], str)
                        else context
                    )
                    stem = (
                        re.sub(r"[^a-z0-9]+", "_", str(hint).lower()).strip("_")
                        or "metadata"
                    )
                    year = value.get("year")
                    if type(year) is int:
                        stem += f"_{abs(year)}" + ("_bce" if year < 0 else "")
                    name, index = stem, 2
                    while name in registry:
                        name, index = f"{stem}_{index}", index + 1
                    names[token] = name
                    registry[name] = copy.deepcopy(value)
                return {"shared_ref": names[token]}
            return {key: replace(item, key) for key, item in value.items()}
        if isinstance(value, list):
            return [replace(item, context) for item in value]
        return value

    for record in result.get("eos_records", []):
        for key in METADATA_FIELDS:
            if key in record:
                record[key] = replace(record[key], key)
    if registry:
        result["shared_metadata"] = registry
    return result
