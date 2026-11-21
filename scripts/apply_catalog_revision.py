"""Repository maintenance: change the scientific content of a bundled material.

One request file (JSON) describes one revision of one material, of one kind:

``add_record``
    append a new EOS record. The request supplies the record fields, its
    ``pressure_calibration`` and its audit ``evidence`` (source check,
    ``parameter_provenance`` and what the refit validator reads, such as
    ``primary_data_check`` and ``fit_datasets``) but no reproduction
    assessment and no ``check`` block: the tool computes the record's row in
    the registered refit ledger with ``scripts/validate_primary_eos_refits.py``,
    derives the assessment from it (``assessment_for_record``) and generates
    ``check``. The material's ``default_records`` is unchanged.
``set_fields``
    set or correct fields of existing records, each named by a path with its
    expected current value. Only the paths in :data:`ALLOW_LIST` are
    settable; a changed reproduction comparison re-derives the record's
    assessment and ``check``.

``remove_record`` and material-level fields (``default_records``, structure,
datasets) are not supported. Request format::

    {"id": "<date>-<material>-<what>", "material": "<material>",
     "kind": "set_fields", "date": "YYYY-MM-DD",
     "decided_by": "<name>", "reason": "<why>",
     "source": {"citation": "<reference>", "location": "<table, page>"},
     "rights": {"reviewed_by": "<name>", "date": "YYYY-MM-DD",
                "statement": "<what was reviewed>",
                "same_basis": true},             # add_record only
     "changes": [{"record": "<record id>", "path": ["parameter_error_confidence"],
                  "expected_absent": true,       # or "expected": <value>
                  "value": 0.6826894921370859}]} # or "remove": true
    # add_record: "record": {...} instead of "changes"

``rights`` is required for every revision because the rights review
(``docs/data/eos-structure-rights-policy.json``) pins the whole reviewed
material; the tool stores the statement, moves the pin and regenerates the
review outputs, and never writes the statement itself.

Everything is computed in a staging copy of the repository: the record file
and sidecar, the audit-split binding, the ledger, the rights policy and
review, the refit-ledger row and page, the primary-source audit entry and
manifest counts, catalog coverage, material summaries, reconstruction reports,
the paper investigation ledger and identical Rust fixture copies. The staged
catalog is then loaded through the library (every material, the
pressure-calibration graph and the revised records with primary validation),
and every recovery check runs on it, before any file is written. ``--undo ID``
reverses a material's most recent revision, also after later text revisions
of the same material; without them the files return byte for byte.

Usage::

    python scripts/apply_catalog_revision.py REQUEST.json --check  # dry run
    python scripts/apply_catalog_revision.py REQUEST.json          # apply
    python scripts/apply_catalog_revision.py --undo ID [--check]
    python scripts/apply_catalog_revision.py --verify

This is not an installed-package tool.
"""

from __future__ import annotations

import argparse
import copy
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from peritheos._eosmat_audit import (  # noqa: E402
    AUDIT_RECORD_KEYS,
    GENERIC_QUALIFICATIONS,
    GENERIC_VALIDATION_NOTES,
    bundled_audit_reference,
    check_block,
    compact_entry,
    expand_entry,
    merge_audit,
)
from peritheos._eosmat_format import format_audit, format_eosmat  # noqa: E402
from peritheos.eosmat import validate_eosmat_document  # noqa: E402
from scripts.catalog_revisions import (  # noqa: E402
    ABSENT,
    HEAD_KEY,
    LEDGER,
    SPLIT_RECEIPT_STATEMENT,
    apply_operations,
    content_sha256,
    fingerprint,
    ledger_text,
    prior_revision_snapshot,
    read_ledger,
    read_path,
    reverse_operations,
    verify_ledger,
)
from scripts.check_record_text import record_cases  # noqa: E402
from scripts.review_eos_structure_rights import (  # noqa: E402
    POLICY,
    build_review,
)
from scripts.scientific_source_snapshot import scientific_snapshot  # noqa: E402
from scripts.split_eosmat_audits import (  # noqa: E402
    LAYOUTS,
    RECEIPT,
    pre_split_bytes,
    receipt_text,
    recover_pre_split,
)

RECEIPT_PATH = RECEIPT.relative_to(ROOT).as_posix()
REFITS = "docs/data/primary-eos-refits.json"
REFITS_PAGE = "docs/primary-eos-refits.md"
SOURCE_AUDIT = "peritheos/data/primary-source-audit.json"
MANIFEST = "peritheos/data/materials/manifest.json"
FIXTURES = "crates/peritheos/tests"
#: Reconstruction inputs that pin whole source records by digest; the refit
#: validator refuses any change to them (scripts/reconstruct_similar_eos_uncertainties.py).
PINNED_RECORDS = "docs/data/uncertainty-reconstruction-inputs.json"
COMMON_KEYS = (
    "id",
    "material",
    "kind",
    "date",
    "decided_by",
    "reason",
    "source",
    "rights",
)
KIND_KEYS = {"add_record": "record", "set_fields": "changes"}
#: Optional request key: the exact changes of the ``check`` headline fields
#: (:data:`HEADLINE_FIELDS`) a revision causes on existing records,
#: ``{record: {field: {"from": old, "to": new}}}``.
ACKNOWLEDGE = "acknowledge_check_change"
HEADLINE_FIELDS = ("reproduction", "compared", "warning", "source_check")
CHANGE_KEYS = ("record", "path", "expected", "expected_absent", "value", "remove")
#: Settable fields of an existing reference-isotherm or thermal parameter. A
#: published coefficient ``value`` is not settable (see ``PUBLISHED_VALUE``).
PARAMETER_FIELDS = ("error", "fixed")
PUBLISHED_VALUE = (
    "correcting a published coefficient value is not supported by this tool "
    "yet; raise it with the maintainer"
)
RANGE_FIELDS = (
    "experimental_pressure_range_gpa",
    "experimental_temperature_range_k",
    "pressure_range_status",
)
VALIDITY_BOUNDS = ("pressure_gpa", "temperature_k", "volume_ratio")
#: Allowed pressure bounds of a range, GPa. The bundled catalog spans
#: -11.3 to 1374 GPa; the window leaves room without admitting typos.
PRESSURE_WINDOW_GPA = (-20.0, 2000.0)
#: Upper temperature bound of a range, K (the catalog reaches 8500 K).
TEMPERATURE_MAX_K = 20000.0
#: ``pressure_range_status`` values of the EOSMAT v5 schema.
PRESSURE_RANGE_STATUSES = (
    "reported_exactly",
    "reported_qualitatively",
    "theoretical",
    "reference_parameterization",
)
#: The only paths a ``set_fields`` revision may set, in words.
ALLOW_LIST = (
    "parameters.<name>.error|fixed and "
    "thermal.parameters.<name>.error|fixed (existing parameters); "
    "parameter_error_confidence; parameter_covariance; temperature_ref; "
    "experimental_pressure_range_gpa, experimental_temperature_range_k, "
    "pressure_range_status; validity.pressure_gpa|temperature_k|volume_ratio "
    "(existing bounds of an existing validity object)"
)
TEXT_FIELDS = ("label", "notes", "validity.notes", "thermal.notes")
#: Blocks a bundled record must carry (see docs/adding-materials-and-eos.md).
REQUIRED_RECORD_KEYS = (
    "identifier",
    "label",
    "model",
    "parameters",
    "reference",
    "determination_method",
    "record_kind",
    "experimental_pressure_range_gpa",
    "experimental_temperature_range_k",
    "pressure_calibration",
    "evidence",
)
VALIDATION_STATUSES = ("primary_source_validated", "not_reproduced", "deferred")
CALIBRATION_STATUSES = (
    "resolved",
    "partially_resolved",
    "not_applicable",
    "unresolved",
)
GENERATORS = (
    "scripts.generate_catalog_coverage",
    "scripts.generate_material_summaries",
    "scripts.generate_reconstruction_reports",
    "scripts.generate_paper_investigation_ledger",
)
_REVISION_ID = re.compile(r"[a-z0-9][a-z0-9._-]*\Z")
_RECORD_ID = re.compile(r"[a-z0-9][a-z0-9_]*\Z")
_STAGED_DIRECTORIES = ("peritheos", "scripts", "docs", "curation", "packages")
_STAGED_FILES = ("peritheos_build.py", "pyproject.toml", "README.md", "mkdocs.yml")
_IGNORED = shutil.ignore_patterns("__pycache__", "*.pyc")


class Refused(ValueError):
    """A request the catalog cannot take; nothing is written."""

    def __init__(self, reason: str, detail: str) -> None:
        super().__init__(f"{reason}: {detail}")
        self.reason = reason
        self.detail = detail


@dataclass
class Plan:
    """Files a revision (or its undo) writes, and what to show the operator."""

    revision: dict
    undo: bool
    files: dict[str, bytes | None] = field(default_factory=dict)
    shown: list[str] = field(default_factory=list)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _paths(material: str) -> tuple[str, str]:
    return (
        f"peritheos/data/materials/{material}.eosmat",
        f"peritheos/data/audits/{material}.audit.json",
    )


def _number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def _same(left: Any, right: Any) -> bool:
    """Equal values; numbers compare numerically (172 equals 172.0)."""
    if left is ABSENT or right is ABSENT:
        return left is right
    if _number(left) and _number(right):
        return float(left) == float(right)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(map(_same, left, right))
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(
            _same(left[k], right[k]) for k in left
        )
    return type(left) is type(right) and left == right


def _text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _date(value: Any, name: str) -> None:
    try:
        datetime.date.fromisoformat(value)
    except (TypeError, ValueError):
        raise Refused("request", f"{name} must be a YYYY-MM-DD date") from None


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


# Request --------------------------------------------------------------------


def check_request(request: Any) -> None:
    """Refuse a request whose shape is wrong, incomplete or carries extra keys."""
    if not isinstance(request, dict):
        raise Refused("request", "a request is one JSON object")
    kind = request.get("kind")
    if kind not in KIND_KEYS:
        raise Refused("request", f"kind must be one of {sorted(KIND_KEYS)}")
    allowed = {*COMMON_KEYS, KIND_KEYS[kind]}
    extra = sorted(set(request) - allowed - {ACKNOWLEDGE})
    if extra:
        raise Refused("undeclared_change", f"{kind} does not take {extra}")
    missing = sorted(allowed - set(request))
    if missing:
        raise Refused(
            "rights" if missing == ["rights"] else "request", f"missing {missing}"
        )
    if not isinstance(request["id"], str) or not _REVISION_ID.match(request["id"]):
        raise Refused("request", "id must be lower-case letters, digits, '.', '-', '_'")
    for key in ("material", "decided_by", "reason"):
        if not _text(request[key]):
            raise Refused("request", f"{key} must be non-empty text")
    _date(request["date"], "date")
    source = request["source"]
    if (
        not isinstance(source, dict)
        or set(source) != {"citation", "location"}
        or not all(_text(v) for v in source.values())
    ):
        raise Refused("request", "source needs exactly citation and location (text)")
    rights = request["rights"]
    if not isinstance(rights, dict):
        raise Refused("rights", "rights must be an object")
    required = {"reviewed_by", "date", "statement"}
    optional = {"same_basis", "differences"}
    if kind == "add_record":
        required.add("same_basis")
        if rights.get("same_basis") is False:
            required.add("differences")
    if not required <= set(rights) or set(rights) - required - optional:
        raise Refused(
            "rights", f"{kind} needs a rights statement with {sorted(required)}"
        )
    for key in ("reviewed_by", "statement", "differences"):
        if key in rights and not _text(rights[key]):
            raise Refused("rights", f"rights.{key} must be non-empty text")
    _date(rights["date"], "rights.date")
    if "same_basis" in rights and not isinstance(rights["same_basis"], bool):
        raise Refused("rights", "rights.same_basis must be true or false")
    acknowledged = request.get(ACKNOWLEDGE, {})
    if not isinstance(acknowledged, dict) or not all(
        isinstance(fields, dict)
        and set(fields) <= set(HEADLINE_FIELDS)
        and all(
            isinstance(pair, dict) and set(pair) == {"from", "to"}
            for pair in fields.values()
        )
        for fields in acknowledged.values()
    ):
        raise Refused(
            "request",
            f"{ACKNOWLEDGE} maps records to {{field: {{from, to}}}} for "
            f"{list(HEADLINE_FIELDS)}",
        )
    if kind == "set_fields":
        changes = request["changes"]
        if not isinstance(changes, list) or not changes:
            raise Refused("request", "changes must be a non-empty list")
        for change in changes:
            if not isinstance(change, dict) or set(change) - set(CHANGE_KEYS):
                raise Refused(
                    "undeclared_change", f"a change takes only {list(CHANGE_KEYS)}"
                )
    elif not isinstance(request["record"], dict):
        raise Refused("request", "record must be an object")


def read_request(path: Path) -> dict:
    """Read and check a request file; a missing or malformed file is refused."""
    try:
        request = json.loads(Path(path).read_text(encoding="utf-8"))
    except OSError as error:
        raise Refused("request", f"cannot read {path}: {error.strerror}") from None
    except ValueError as error:
        raise Refused("request", f"{path} is not valid JSON: {error}") from None
    check_request(request)
    return request


# Files and staging ----------------------------------------------------------


def _read(root: Path, material: str) -> tuple[dict, dict]:
    record_path, audit_path = (root / p for p in _paths(material))
    if not record_path.exists() or not audit_path.exists():
        raise Refused("unknown_material", f"{material} has no record file and sidecar")
    record_text = record_path.read_text(encoding="utf-8")
    audit_text = audit_path.read_text(encoding="utf-8")
    document, sidecar = json.loads(record_text), json.loads(audit_text)
    if format_eosmat(document) != record_text or format_audit(sidecar) != audit_text:
        raise Refused("not_canonical", f"{material}: files are not in canonical layout")
    return document, sidecar


def copy_repository(root: Path, stage: Path) -> None:
    """Copy the repository parts a revision reads or writes to ``stage``."""
    # Contents only: the copy stays writable whatever the tree's file modes.
    for relative in _STAGED_DIRECTORIES:
        if (root / relative).is_dir():
            shutil.copytree(
                root / relative,
                stage / relative,
                ignore=_IGNORED,
                copy_function=shutil.copyfile,
            )
    for relative in _STAGED_FILES:
        if (root / relative).is_file():
            shutil.copyfile(root / relative, stage / relative)
    for relative in (f"{FIXTURES}/data", f"{FIXTURES}/audits"):
        if (root / relative).is_dir():
            shutil.copytree(
                root / relative, stage / relative, copy_function=shutil.copyfile
            )
    for relative in (RECEIPT_PATH, LEDGER, POLICY):
        if (root / relative).exists() and not (stage / relative).exists():
            (stage / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / relative, stage / relative)


@contextmanager
def _stage(root: Path) -> Iterator[Path]:
    """A scratch copy of the repository parts a revision reads or writes."""
    with tempfile.TemporaryDirectory(prefix="catalog-revision-") as name:
        copy_repository(root, Path(name))
        yield Path(name)


def _tree(base: Path) -> dict[str, Path]:
    """Every staged file under ``base`` that a revision may change."""
    result = {}
    roots = [*_STAGED_DIRECTORIES, f"{FIXTURES}/data", f"{FIXTURES}/audits"]
    for relative in roots:
        directory = base / relative
        if not directory.is_dir():
            continue
        for path in directory.rglob("*"):
            if (
                path.is_file()
                and "__pycache__" not in path.parts
                and path.suffix not in (".pyc", ".so", ".pyd")
            ):
                result[path.relative_to(base).as_posix()] = path
    for relative in _STAGED_FILES:
        if (base / relative).is_file():
            result[relative] = base / relative
    return result


def _changed_files(stage: Path, root: Path) -> dict[str, bytes | None]:
    """Files whose staged bytes differ from the tree (``None``: removed)."""
    staged, current = _tree(stage), _tree(root)
    files: dict[str, bytes | None] = {}
    for relative, path in staged.items():
        content = path.read_bytes()
        old = current.get(relative)
        if old is None or old.read_bytes() != content:
            files[relative] = content
    for relative in sorted(set(current) - set(staged)):
        files[relative] = None
    return dict(sorted(files.items()))


def _run(stage: Path, arguments: list[str], what: str) -> str:
    """Run Python in the stage, importing the staged package and scripts."""
    environment = {
        **os.environ,
        "PYTHONPATH": str(stage),
        "PYTHONDONTWRITEBYTECODE": "1",
        "MPLBACKEND": "Agg",
    }
    result = subprocess.run(
        [sys.executable, *arguments],
        cwd=stage,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        tail = (result.stderr or result.stdout).strip().splitlines()[-3:]
        raise Refused(what, " | ".join(tail) or f"exit status {result.returncode}")
    return result.stdout


_SMOKE = """
import json, math, sys
from peritheos import get_material, list_eos_records, list_materials
from peritheos.eosmat import (
    get_material_document, list_eos_record_documents, list_material_documents,
    validate_pressure_calibration_references,
)
from peritheos.eos import ThermalEOS
from peritheos.eos.maltby2024 import Maltby2024
from peritheos.materials import Material

material, touched = sys.argv[1], json.loads(sys.argv[2])
validate_pressure_calibration_references()
list_eos_record_documents()
list_eos_records()
for item in list_materials():
    get_material(item.identifier)
for identifier in list_material_documents():
    get_material_document(identifier)
document = get_material_document(material)
records = {r["identifier"]: r for r in document["eos_records"]}
for identifier in touched:
    status = records[identifier]["scientific_validation"]["status"]
    if status == "deferred":
        continue
    loaded = Material.from_eosmat(
        document,
        record_identifiers=[identifier],
        require_primary_validation=status == "primary_source_validated",
    )
    record = loaded.get_eos_record(identifier)
    if record.is_hugoniot:
        record.pressure(record.reference_volume * 0.95, check_domain=False)
        continue
    # A state inside the enforced range, evaluated with uncertainty.
    validity = record.validity
    low, high = validity.pressure_gpa
    pressure = 0.5 * (max(low, -10.0) + min(high, 100.0))
    if not low <= pressure <= high:
        pressure = low if math.isfinite(low) else high
    thermal = isinstance(record.eos, (ThermalEOS, Maltby2024))
    temperature = None
    if thermal:
        t_low, t_high = validity.temperature_k
        temperature = record.reference_temperature
        if not t_low <= temperature <= t_high:
            temperature = 0.5 * (t_low + min(t_high, 4000.0))
    volume = record.volume(pressure, temperature)
    record.pressure(volume, temperature)
    result = record.pressure_with_uncertainty(volume, temperature)
    values = [result.value] if hasattr(result, "value") else []
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError(f"{identifier}: non-finite pressure with uncertainty")
print("ok")
"""

_ROWS = """
import json, sys
from scripts.validate_primary_eos_refits import validate_all
print(json.dumps(validate_all([sys.argv[1]])["records"], allow_nan=False))
"""


def _load_catalog(stage: Path, material: str, touched: list[str]) -> None:
    """Load the staged catalog through the library; refuse on any failure."""
    _run(stage, ["-c", _SMOKE, material, json.dumps(touched)], "not_loadable")


def _fresh_rows(stage: Path, material: str) -> dict[str, dict]:
    """The material's refit-ledger rows computed by the validator itself."""
    output = _run(stage, ["-c", _ROWS, material], "refit_ledger")
    rows = json.loads(output)
    # The saved ledger is written with sorted keys at every level.
    return {
        row["record_identifier"]: json.loads(json.dumps(row, sort_keys=True))
        for row in rows
    }


# Bookkeeping files ------------------------------------------------------------


def _refit_ledger(stage: Path) -> dict:
    from scripts.validate_primary_eos_refits import render_markdown

    path = stage / REFITS
    data = path.read_bytes()
    ledger = json.loads(data)
    if _json_bytes(ledger) != data or render_markdown(ledger) != (
        stage / REFITS_PAGE
    ).read_text(encoding="utf-8"):
        raise Refused(
            "refit_ledger", f"{REFITS} or {REFITS_PAGE} is not in its saved layout"
        )
    if ledger["summary"] != _refit_summary(ledger["records"]):
        raise Refused("refit_ledger", f"{REFITS} summary does not count its rows")
    return ledger


def _refit_summary(rows: list[dict]) -> dict:
    counts = Counter(row["status"] for row in rows)
    return {"total": len(rows), **dict(sorted(counts.items()))}


def _write_refit_ledger(stage: Path, ledger: dict) -> None:
    from scripts.validate_primary_eos_refits import render_markdown

    ledger["summary"] = _refit_summary(ledger["records"])
    (stage / REFITS).write_bytes(_json_bytes(ledger))
    (stage / REFITS_PAGE).write_text(render_markdown(ledger), encoding="utf-8")


def _insert_row(rows: list[dict], row: dict, material: str) -> None:
    """Insert after the material's last row (rows follow catalog order)."""
    positions = [i for i, r in enumerate(rows) if r["material_identifier"] == material]
    rows.insert(positions[-1] + 1 if positions else len(rows), row)


def _catalog_counts(stage: Path) -> tuple[list[dict], dict, dict]:
    """Source-audit entries, validation counts and calibration counts."""
    from scripts.apply_primary_source_audit import audit_entry
    from scripts.eosmat_io import read_material

    existing = {
        e["record"]: e
        for e in json.loads((stage / SOURCE_AUDIT).read_bytes())["records"]
    }
    entries, calibrations = [], []
    for path in sorted((stage / "peritheos/data/materials").glob("*.eosmat")):
        document = read_material(path)
        for record in document["eos_records"]:
            entries.append(
                audit_entry(
                    document, path.name, record, existing.get(record["identifier"])
                )
            )
            calibrations.append(record["pressure_calibration"])
    validation = dict(sorted(Counter(e["status"] for e in entries).items()))
    return entries, validation, _calibration_counts(calibrations)


def _calibration_counts(calibrations: list[dict]) -> dict:
    statuses = Counter(c["status"] for c in calibrations)
    return {
        "status_counts": {s: statuses[s] for s in CALIBRATION_STATUSES},
        "recalculation_counts": dict(
            sorted(Counter(c["recalculation"]["status"] for c in calibrations).items())
        ),
        "resolvable_reference_eos_uses": sum(
            "reference_eos_record" in m for c in calibrations for m in c["methods"]
        ),
        "resolvable_reference_calibration_uses": sum(
            "reference_calibration_record" in m
            for c in calibrations
            for m in c["methods"]
        ),
    }


def _manifest_with(
    manifest: dict, records: int, validation: dict, calibration: dict
) -> dict:
    """The manifest with recomputed counts; every other field and key order kept."""
    result = copy.deepcopy(manifest)
    result["materials"] = manifest["materials"]
    result["eos_records"] = records
    result["scientific_validation"]["counts"] = validation
    stored = result["pressure_calibration"]
    order = list(stored["status_counts"])
    stored["status_counts"] = {
        s: calibration["status_counts"][s]
        for s in order + [s for s in CALIBRATION_STATUSES if s not in order]
    }
    for key in (
        "recalculation_counts",
        "resolvable_reference_eos_uses",
        "resolvable_reference_calibration_uses",
    ):
        stored[key] = calibration[key]
    return result


def _bookkeeping(stage: Path, material: str) -> None:
    """Source-audit entries and manifest counts of the staged catalog.

    Entries of records still in the catalog are kept as they are; a record
    that is new gets the entry the audit builder gives it, after its
    material's last entry, and a record no longer in the catalog loses its
    entry. Only the counts of the manifest change.
    """
    report_path, manifest_path = stage / SOURCE_AUDIT, stage / MANIFEST
    report = json.loads(report_path.read_bytes())
    manifest = json.loads(manifest_path.read_bytes())
    entries, validation, calibration = _catalog_counts(stage)
    current = {e["record"] for e in entries}
    kept = [e for e in report["records"] if e["record"] in current]
    present = {e["record"] for e in kept}
    for entry in entries:
        if entry["record"] in present:
            continue
        positions = [i for i, e in enumerate(kept) if e["material"] == material]
        kept.insert(positions[-1] + 1 if positions else len(kept), entry)
    counts = Counter(e["status"] for e in kept)
    report["records"] = kept
    report["summary"] = {"records": len(kept), **dict(sorted(counts.items()))}
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps(
            _manifest_with(manifest, len(entries), validation, calibration),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )


def _check_bookkeeping(stage: Path) -> None:
    """Refuse when the source audit or manifest does not count the catalog."""
    report = json.loads((stage / SOURCE_AUDIT).read_bytes())
    manifest = json.loads((stage / MANIFEST).read_bytes())
    entries, validation, calibration = _catalog_counts(stage)
    if report["summary"] != {"records": len(entries), **validation} or (
        {e["record"] for e in report["records"]} != {e["record"] for e in entries}
    ):
        raise Refused("bookkeeping", f"{SOURCE_AUDIT} does not cover the catalog")
    if _manifest_with(manifest, len(entries), validation, calibration) != manifest:
        raise Refused("bookkeeping", f"{MANIFEST} does not count the catalog")


def _fixtures(stage: Path, root: Path, material: str) -> None:
    """Keep Rust fixture copies that equal the bundled files equal to them."""
    record_path, audit_path = _paths(material)
    for bundled, fixture in (
        (record_path, f"{FIXTURES}/data/{material}.eosmat"),
        (audit_path, f"{FIXTURES}/audits/{material}.audit.json"),
    ):
        if (root / fixture).is_file() and (root / fixture).read_bytes() == (
            root / bundled
        ).read_bytes():
            shutil.copy2(stage / bundled, stage / fixture)


def _generate(stage: Path) -> None:
    for module in GENERATORS:
        _run(stage, ["-m", module], "generator")


# Building a revision -----------------------------------------------------


def _catalog_names(root: Path) -> tuple[set[str], set[str]]:
    """Record identifiers and aliases, and every record key, in the catalog."""
    names: set[str] = set()
    keys: set[str] = set(AUDIT_RECORD_KEYS)
    for path in (root / "peritheos/data/materials").glob("*.eosmat"):
        for record in json.loads(path.read_bytes()).get("eos_records", []):
            names.add(record.get("identifier"))
            names.update(record.get("aliases") or [])
            keys.update(record)
    keys.discard("check")
    return names, keys


def _require(condition: bool, detail: str) -> None:
    if not condition:
        raise Refused("record", detail)


def _new_record(material: str, record: dict, document: dict, root: Path) -> dict:
    """Check an added record against what the catalog requires of it."""
    identifier = record.get("identifier")
    if not isinstance(identifier, str) or not _RECORD_ID.match(identifier):
        raise Refused("record", "identifier must be lower_snake_case")
    for key in ("check", "superseded_record_text", "default", "default_for"):
        if key in record:
            raise Refused(
                "undeclared_change" if key.startswith("default") else "record",
                f"{key} is not given in a request "
                + (
                    "(it would change the default selection)"
                    if key.startswith("default")
                    else "(it is generated)"
                ),
            )
    names, keys = _catalog_names(root)
    unknown = sorted(set(record) - keys)
    if unknown:
        raise Refused("record", f"unknown record keys {unknown}")
    missing = [key for key in REQUIRED_RECORD_KEYS if key not in record]
    if missing:
        raise Refused("record", f"a bundled record needs {missing}")
    for name in [identifier, *(record.get("aliases") or [])]:
        if name in names:
            raise Refused("identifier_taken", f"{name} is already used in the catalog")
    calibration = record["pressure_calibration"]
    _require(
        isinstance(calibration, dict)
        and calibration.get("status") in CALIBRATION_STATUSES
        and isinstance(calibration.get("methods"), list)
        and isinstance(calibration.get("recalculation"), dict)
        and _text(calibration["recalculation"].get("status"))
        and _text(calibration.get("audit_date")),
        "pressure_calibration needs status ("
        + ", ".join(CALIBRATION_STATUSES)
        + "), methods, recalculation.status and audit_date",
    )
    evidence = record["evidence"]
    _require(isinstance(evidence, dict), "evidence must be an object")
    validation = evidence.get("scientific_validation")
    _require(
        isinstance(validation, dict)
        and validation.get("status") in VALIDATION_STATUSES
        and _text(validation.get("audit_date"))
        and isinstance(validation.get("primary_source_check"), dict)
        and _text(validation["primary_source_check"].get("finding")),
        "evidence.scientific_validation needs status ("
        + ", ".join(VALIDATION_STATUSES)
        + "), audit_date and primary_source_check.finding",
    )
    _date(validation["audit_date"], "scientific_validation.audit_date")
    _require(
        isinstance(evidence.get("parameter_provenance"), dict)
        and bool(evidence["parameter_provenance"]),
        "evidence.parameter_provenance must say where each value comes from",
    )
    reproduction = evidence.get("reproduction")
    if isinstance(reproduction, dict) and "assessment" in reproduction:
        raise Refused(
            "derived_field",
            "evidence.reproduction.assessment is derived from the record's refit "
            "ledger row; leave it out",
        )
    if validation.get("note") in GENERIC_VALIDATION_NOTES:
        raise Refused(
            "generic_statement", "the generic note is stated once in the docs"
        )
    _check_record(record, identifier, complete=True)
    labels = [r.get("label") for r in document["eos_records"]]
    cases = record_cases({**record, "check": {}})
    if record["label"] in labels:
        cases.append("label_unique")
    if cases:
        raise Refused(
            "record_text",
            f"outside the record text standard: {cases} (docs/record-text-standard.md)",
        )
    return copy.deepcopy(record)


def _allowed(record: dict, path: list[str], change: Mapping[str, Any]) -> None:
    """Refuse a path outside :data:`ALLOW_LIST` or a value of the wrong kind."""
    joined = ".".join(path)
    if joined in TEXT_FIELDS or (path and path[0] in ("label", "notes")):
        raise Refused(
            "protected_field",
            f"{joined} is text; use scripts/apply_record_text.py (label, notes, "
            "validity.notes, thermal.notes) or scripts/apply_check_text.py",
        )
    if path[-1] == "value" and (
        (len(path) == 3 and path[0] == "parameters")
        or (len(path) == 4 and path[:2] == ["thermal", "parameters"])
    ):
        raise Refused("published_value", f"{joined}: {PUBLISHED_VALUE}")
    removal = bool(change.get("remove"))
    value = change.get("value")
    parameters = record.get("parameters")
    thermal = record.get("thermal")
    kind = None
    if (
        len(path) == 3
        and path[0] == "parameters"
        and isinstance(parameters, dict)
        and isinstance(parameters.get(path[1]), dict)
        and path[2] in PARAMETER_FIELDS
    ):
        kind = path[2]
    elif (
        len(path) == 4
        and path[:2] == ["thermal", "parameters"]
        and isinstance(thermal, dict)
        and isinstance(thermal.get("parameters"), dict)
        and isinstance(thermal["parameters"].get(path[2]), dict)
        and path[3] in PARAMETER_FIELDS
    ):
        kind = path[3]
    elif len(path) == 1 and path[0] in (
        "parameter_error_confidence",
        "parameter_covariance",
        "temperature_ref",
        *RANGE_FIELDS,
    ):
        kind = path[0]
    elif (
        len(path) == 2
        and path[0] == "validity"
        and isinstance(record.get("validity"), dict)
        and path[1] in VALIDITY_BOUNDS
        and path[1] in record["validity"]
    ):
        kind = "bound"
    if kind is None:
        raise Refused(
            "protected_field",
            f"{joined} is not settable by a catalog revision; settable: {ALLOW_LIST}",
        )
    if removal and kind not in ("error", "fixed", "parameter_covariance"):
        raise Refused("value", f"{joined} cannot be removed")
    if removal:
        return

    def pair(item: Any) -> bool:
        return (
            isinstance(item, list)
            and len(item) == 2
            and all(map(_number, item))
            and item[0] <= item[1]
        )

    if kind in RANGE_FIELDS or kind == "bound":
        _check_range(path[-1], value, joined)
        return
    if kind == "parameter_covariance":
        return  # checked on the record after all changes (_check_record)
    valid = {
        "error": value is None or (_number(value) and value >= 0),
        "fixed": isinstance(value, bool),
        "parameter_error_confidence": value is None
        or (_number(value) and 0 < value < 1),
        "parameter_covariance": value is None or isinstance(value, dict),
        "temperature_ref": _number(value) and value > 0,
        "experimental_pressure_range_gpa": pair(value),
        "experimental_temperature_range_k": pair(value),
        "pressure_range_status": _text(value),
        "bound": pair(value),
    }[kind]
    if not valid:
        raise Refused("value", f"{joined}: {json.dumps(value)} is not a valid value")


def _check_range(kind: str, value: Any, where: str) -> None:
    """A range or status value: two finite numbers inside the sane window."""
    if kind == "pressure_range_status":
        if value not in PRESSURE_RANGE_STATUSES:
            raise Refused(
                "value",
                f"{where}: {json.dumps(value)} is not one of {list(PRESSURE_RANGE_STATUSES)}",
            )
        return
    if not (
        isinstance(value, list)
        and len(value) == 2
        and all(map(_number, value))
        and value[0] <= value[1]
    ):
        raise Refused(
            "range",
            f"{where}: {json.dumps(value)} is not [lower, upper] finite numbers",
        )
    if "pressure" in kind:
        low, high = PRESSURE_WINDOW_GPA
        if value[0] < low or value[1] > high:
            raise Refused(
                "range", f"{where}: {value} GPa is outside [{low}, {high}] GPa"
            )
    elif "temperature" in kind:
        if value[0] <= 0 or value[1] > TEMPERATURE_MAX_K:
            raise Refused(
                "range", f"{where}: {value} K is outside (0, {TEMPERATURE_MAX_K}] K"
            )
    elif kind == "volume_ratio" and value[0] <= 0:
        raise Refused("range", f"{where}: {value} must be positive")


def _check_covariance(record: Mapping[str, Any], where: str) -> None:
    """Covariance of free parameters: square, symmetric, PSD, finite."""
    import numpy as np

    covariance = record.get("parameter_covariance")
    if covariance is None:
        return
    if not isinstance(covariance, dict):
        raise Refused("covariance", f"{where}: parameter_covariance must be an object")
    order = covariance.get("parameter_order")
    matrix = covariance.get("matrix")
    if (
        not isinstance(order, list)
        or not order
        or not all(isinstance(n, str) for n in order)
    ):
        raise Refused(
            "covariance", f"{where}: parameter_order must list parameter names"
        )
    if len(set(order)) != len(order):
        raise Refused("covariance", f"{where}: parameter_order repeats a name")
    items = _parameter_items(record)
    for name in order:
        item = items.get(name.split(".")[-1])
        if item is None:
            raise Refused(
                "covariance", f"{where}: {name} is not a parameter of the record"
            )
        if item.get("fixed"):
            raise Refused(
                "covariance", f"{where}: {name} is fixed, not a free parameter"
            )
    if (
        not isinstance(matrix, list)
        or len(matrix) != len(order)
        or not all(isinstance(r, list) and len(r) == len(order) for r in matrix)
        or not all(_number(v) for r in matrix for v in r)
    ):
        raise Refused(
            "covariance",
            f"{where}: matrix must be {len(order)}x{len(order)} finite numbers",
        )
    array = np.asarray(matrix, dtype=float)
    scale = max(float(np.max(np.abs(array))), 1e-300)
    if not np.allclose(array, array.T, rtol=1e-9, atol=1e-12 * scale):
        raise Refused("covariance", f"{where}: matrix is not symmetric")
    if float(np.min(np.linalg.eigvalsh((array + array.T) / 2))) < -1e-9 * scale:
        raise Refused("covariance", f"{where}: matrix is not positive semidefinite")


def _check_record(record: Mapping[str, Any], where: str, *, complete: bool) -> None:
    """Values of a record after the revision (``complete``: an added record)."""
    for name, item in _parameter_items(record).items():
        error = item.get("error")
        if error is not None and not (_number(error) and error >= 0):
            raise Refused(
                "value", f"{where}: {name}.error must be a finite number >= 0 or null"
            )
        if complete and not _number(item.get("value")):
            raise Refused("value", f"{where}: {name}.value must be a finite number")
        if "fixed" in item and not isinstance(item["fixed"], bool):
            raise Refused("value", f"{where}: {name}.fixed must be true or false")
    confidence = record.get("parameter_error_confidence")
    if confidence is not None and not (_number(confidence) and 0 < confidence < 1):
        raise Refused(
            "value", f"{where}: parameter_error_confidence must be in (0, 1) or null"
        )
    if complete:
        for key in RANGE_FIELDS:
            if record.get(key) is not None:
                _check_range(key, record[key], f"{where}/{key}")
        validity = record.get("validity")
        if isinstance(validity, dict):
            for key in VALIDITY_BOUNDS:
                if validity.get(key) is not None:
                    _check_range(key, validity[key], f"{where}/validity.{key}")
    _check_covariance(record, where)
    reference = record.get("temperature_ref")
    if reference is None:
        return
    if not (_number(reference) and reference > 0):
        raise Refused("value", f"{where}: temperature_ref must be a positive number")
    thermal = record.get("thermal") if isinstance(record.get("thermal"), dict) else {}
    tr = (thermal.get("parameters") or {}).get("Tr")
    if isinstance(tr, dict) and _number(tr.get("value")) and tr["value"] != reference:
        raise Refused(
            "value",
            f"{where}: temperature_ref {reference} differs from thermal Tr {tr['value']}",
        )
    span = record.get("experimental_temperature_range_k")
    if isinstance(span, list) and len(span) == 2 and all(map(_number, span)):
        if not span[0] <= reference <= span[1]:
            raise Refused(
                "value",
                f"{where}: temperature_ref {reference} K is outside the stated "
                f"temperature range {span}",
            )


def _check_pins(root: Path, identifiers: list[Any]) -> None:
    """Refuse records whose whole source record another receipt pins."""
    path = root / PINNED_RECORDS
    if not path.exists():
        return
    pinned = json.loads(path.read_bytes()).get("source_record_sha256", {})
    for identifier in identifiers:
        if identifier in pinned:
            raise Refused(
                "pinned_record",
                f"{identifier}: {PINNED_RECORDS} pins its source record by digest "
                "for scripts/reconstruct_similar_eos_uncertainties.py, which the "
                "refit validator runs; a revision of this record is not supported",
            )


def _largest_difference(saved: Any, fresh: Any) -> str:
    """The largest relative numerical difference between two rows, as text."""
    largest = [0.0, ""]

    def walk(a: Any, b: Any, path: str) -> None:
        if _number(a) and _number(b):
            scale = max(abs(a), abs(b), 1e-300)
            relative = abs(a - b) / scale
            if relative > largest[0]:
                largest[:] = [relative, path]
        elif isinstance(a, dict) and isinstance(b, dict):
            for key in a.keys() & b.keys():
                walk(a[key], b[key], f"{path}.{key}" if path else key)
        elif isinstance(a, list) and isinstance(b, list):
            for index, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{path}[{index}]")

    walk(saved, fresh, "")
    return f"{largest[0]:.2g} at {largest[1]}" if largest[1] else "not numerical"


def _set_fields(document: dict, changes: list[dict]) -> list[dict]:
    """Apply declared field changes to the record file; return the operations."""
    records = {r["identifier"]: r for r in document["eos_records"]}
    forward: list[dict] = []
    for change in changes:
        identifier = change.get("record")
        if identifier not in records:
            raise Refused("unknown_record", f"{identifier!r} is not in the material")
        path = change.get("path")
        if (
            not isinstance(path, list)
            or not path
            or not all(isinstance(k, str) and k for k in path)
        ):
            raise Refused("path", "path is a non-empty list of field names")
        if ("expected" in change) == bool(change.get("expected_absent")):
            raise Refused("expected", "give exactly one of expected or expected_absent")
        if ("value" in change) == bool(change.get("remove")):
            raise Refused("value", "give exactly one of value or remove")
        _allowed(records[identifier], path, change)
        for other in forward:
            shorter = min(len(other["path"]), len(path))
            if (
                other["record"] == identifier
                and other["path"][:shorter] == path[:shorter]
            ):
                raise Refused("overlapping_paths", f"{identifier}/{'.'.join(path)}")
        current = read_path(records[identifier], path)
        expected = ABSENT if change.get("expected_absent") else change["expected"]
        if not _same(current, expected):
            raise Refused(
                "expected_mismatch",
                f"{identifier}/{'.'.join(path)} is "
                f"{'absent' if current is ABSENT else json.dumps(current)}",
            )
        after_present = "value" in change
        if _same(current, change["value"] if after_present else ABSENT):
            raise Refused("no_change", f"{identifier}/{'.'.join(path)} already has it")
        operation = {
            "op": "set",
            "record": identifier,
            "path": list(path),
            "before_present": current is not ABSENT,
            "before": None if current is ABSENT else current,
            "after_present": after_present,
            "after": copy.deepcopy(change["value"]) if after_present else None,
        }
        apply_operations(document["eos_records"], [operation], "document")
        forward.append(operation)
    for identifier in dict.fromkeys(o["record"] for o in forward):
        _check_record(records[identifier], identifier, complete=False)
    return forward


def _apply_to_files(document: dict, sidecar: dict, operations: list[dict]) -> None:
    """Apply operations to a record file and sidecar (audit paths: sidecar)."""
    material = str(document["identifier"])
    for operation in operations:
        if operation["op"] == "add_record":
            record = copy.deepcopy(operation["record"])
            identifier = record["identifier"]
            audit = {k: record.pop(k) for k in AUDIT_RECORD_KEYS if k in record}
            document["eos_records"].append(record)
            sidecar["records"][identifier] = compact_entry(audit, material, identifier)
        elif operation["op"] == "remove_record":
            apply_operations(document["eos_records"], [operation], "document")
            sidecar["records"].pop(operation["record"], None)
        elif operation["path"][0] in AUDIT_RECORD_KEYS:
            identifier = operation["record"]
            entry = expand_entry(sidecar["records"][identifier], material, identifier)
            trial = [{"identifier": identifier, **entry}]
            apply_operations(trial, [operation], "document")
            trial[0].pop("identifier")
            sidecar["records"][identifier] = compact_entry(
                trial[0], material, identifier
            )
        else:
            apply_operations(document["eos_records"], [operation], "document")


def _refresh_checks(document: dict, sidecar: dict, identifiers: set[str]) -> None:
    merged = merge_audit(document, sidecar)
    complete = {r["identifier"]: r for r in merged["eos_records"]}
    material = str(document["identifier"])
    for record in document["eos_records"]:
        identifier = record["identifier"]
        if identifier in identifiers:
            audit = record.get("check", {}).get("audit") or bundled_audit_reference(
                material, identifier
            )
            check = check_block(
                {k: v for k, v in complete[identifier].items() if k != "check"},
                merged,
                audit=audit,
            )
            items = [(k, v) for k, v in record.items() if k != "check"]
            record.clear()
            record.update({**dict(items), "check": check})


def _files(document: dict, sidecar: dict) -> tuple[bytes, bytes]:
    order = [r["identifier"] for r in document["eos_records"]]
    entries = sidecar["records"]
    sidecar["records"] = {i: entries[i] for i in order if i in entries}
    return (
        format_eosmat(document).encode("utf-8"),
        format_audit(sidecar).encode("utf-8"),
    )


def _merged(files: tuple[bytes, bytes]) -> dict:
    try:
        record_file = json.loads(files[0])
        validate_eosmat_document(record_file)
        merged = merge_audit(record_file, json.loads(files[1]))
        validate_eosmat_document(merged)
    except Exception as error:  # noqa: BLE001 - reported as a refusal
        raise Refused("invalid", str(error)) from None
    return merged


def _validity(record: Mapping[str, Any]) -> Any:
    from peritheos.materials import _eosmat_validity

    return _eosmat_validity(record)


def _parameter_items(record: Mapping[str, Any]) -> dict[str, dict]:
    """A record's reference-isotherm and thermal parameters by name."""
    thermal = record.get("thermal") if isinstance(record.get("thermal"), dict) else {}
    return {
        **{
            k: v
            for k, v in (thermal.get("parameters") or {}).items()
            if isinstance(v, dict)
        },
        **{
            k: v
            for k, v in (record.get("parameters") or {}).items()
            if isinstance(v, dict)
        },
    }


def _row_disagreements(
    record: Mapping[str, Any], row: Mapping[str, Any], names: set[str] | None
) -> list[str]:
    """Published values, errors or fixed flags of ``row`` the record does not hold."""
    items = _parameter_items(record)
    found = []
    for entry in row.get("parameters") or []:
        name = str(entry.get("parameter", "")).split(".")[-1]
        item = items.get(name)
        if item is None or (names is not None and name not in names):
            continue
        if "published" in entry and not _same(entry["published"], item.get("value")):
            found.append(
                f"{name} published {entry['published']} vs {item.get('value')}"
            )
        if "published_error" in entry and not _same(
            entry["published_error"], item.get("error")
        ):
            found.append(
                f"{name} published_error {entry['published_error']} vs {item.get('error')}"
            )
    fixed = {str(n).split(".")[-1] for n in row.get("fixed_parameters") or []}
    for name, item in items.items():
        if names is not None and name not in names:
            continue
        compared = any(
            str(e.get("parameter", "")).split(".")[-1] == name
            for e in row.get("parameters") or []
        )
        if bool(item.get("fixed")) != (name in fixed) and (compared or name in fixed):
            found.append(f"{name} fixed {bool(item.get('fixed'))} vs row")
    return found


def _pinning_scripts(identifier: str, stage: Path) -> list[str]:
    """Repository scripts that name a record (where its row comes from)."""
    own = {"apply_catalog_revision.py", "catalog_revisions.py"}
    return sorted(
        f"scripts/{path.name}"
        for path in (stage / "scripts").glob("*.py")
        if path.name not in own and identifier in path.read_text(encoding="utf-8")
    )


def _check_rows(
    document: dict,
    forward: list[dict],
    rows: Mapping[str, dict],
    kind: str,
    stage: Path,
) -> None:
    """Refuse when a fresh refit-ledger row does not reflect the stored record."""
    records = {r["identifier"]: r for r in document["eos_records"]}
    for identifier, row in rows.items():
        names = None
        if kind == "set_fields":
            names = {
                o["path"][-2]
                for o in forward
                if o["op"] == "set"
                and o["record"] == identifier
                and o["path"][0] in ("parameters", "thermal")
            }
            if not names:
                continue
        found = _row_disagreements(records[identifier], row, names)
        if found:
            scripts = _pinning_scripts(identifier, stage) or [
                "scripts/validate_primary_eos_refits.py"
            ]
            raise Refused(
                "row_disagrees",
                f"{identifier}: the refit-ledger row does not take the stored values "
                f"({'; '.join(found)}); they are pinned in {', '.join(scripts)}",
            )


def _derive_assessments(
    stage: Path,
    material: str,
    document: dict,
    sidecar: dict,
    rows: Mapping[str, dict],
) -> list[dict]:
    """Set each named record's assessment to its derivation from its row.

    Returns the derived ``set`` operations (assessment path, complete form).
    """
    from scripts.standardize_reconstruction_metadata import (
        STATUS,
        assessment_for_record,
    )

    merged = merge_audit(document, sidecar)
    records = {r["identifier"]: r for r in merged["eos_records"]}
    operations = []
    for identifier, row in rows.items():
        record = records[identifier]
        current = record.get("evidence", {}).get("reproduction", {}).get("assessment")
        if row.get("status") not in STATUS:
            raise Refused(
                "refit_failed",
                f"{identifier}: the refit validator gives status "
                f"{row.get('status')!r} ({row.get('reason', 'no reason given')})",
            )
        try:
            derived = assessment_for_record(record, row)
        except (KeyError, ValueError, TypeError) as error:
            raise Refused(
                "refit_failed", f"{identifier}: no assessment from the row ({error!r})"
            ) from None
        if set(derived.get("qualifications") or []) & set(GENERIC_QUALIFICATIONS):
            raise Refused(
                "generic_statement", f"{identifier}: derived assessment is generic"
            )
        if current is not None and current == derived:
            continue
        path = ["evidence", "reproduction", "assessment"]
        if current is None and not isinstance(
            record.get("evidence", {}).get("reproduction"), dict
        ):
            path = ["evidence", "reproduction"]
            derived = {
                **record["evidence"].get("reproduction", {}),
                "assessment": derived,
            }
            current = ABSENT
        operation = {
            "op": "set",
            "record": identifier,
            "path": path,
            "derived": True,
            "before_present": current is not None and current is not ABSENT,
            "before": None if current is None or current is ABSENT else current,
            "after_present": True,
            "after": derived,
        }
        _apply_to_files(document, sidecar, [operation])
        operations.append(operation)
    return operations


def _revised_binding(
    binding: dict,
    document: dict,
    forward: list[dict],
    files: tuple[bytes, bytes],
    pins: tuple[str, str] | None,
) -> dict:
    """The audit-split binding of revised files (see scripts/catalog_revisions).

    ``forward`` are the operations applied since ``binding`` (for an undo:
    the reverse operations, with ``document`` the record file before them).
    ``pins`` are the sha256 of the revised reference and the content digest
    of the latest revision; ``None`` removes them.
    """
    result = copy.deepcopy(binding)
    result["after_record_sha256"] = _sha(files[0])
    result["after_audit_sha256"] = _sha(files[1])
    order = [r["identifier"] for r in document["eos_records"]]
    for operation in forward:
        if operation["op"] == "add_record":
            # Its check block is generated, not part of the reference.
            index = order.index(operation["record"]["identifier"])
            result["added_check"] = sorted({*result.get("added_check", []), index})
        elif operation["op"] == "remove_record":
            # Appended last, so no other record's index changes.
            index = order.index(operation["record"])
            for key in ("added_check", "superseded_record_text"):
                if key in result:
                    result[key] = [i for i in result[key] if i != index]
            if "check_text" in result:
                result["check_text"] = [
                    e for e in result["check_text"] if e[0] != index
                ]
                if not result["check_text"]:
                    del result["check_text"]
            if result.get("superseded_record_text") == []:
                del result["superseded_record_text"]
        elif operation["before_present"] != operation["after_present"]:
            # A recorded key order of the node must name the new key set.
            parent = ["eos_records", order.index(operation["record"])]
            parent += operation["path"][:-1]
            for ordering in result.get("object_orders", []):
                if ordering["path"] == parent:
                    key = operation["path"][-1]
                    if operation["after_present"]:
                        ordering["keys"].append(key)
                    else:
                        ordering["keys"].remove(key)
    ordered: dict[str, Any] = {}
    for key, value in result.items():
        if key in ("catalog_revision_sha256", HEAD_KEY):
            continue
        ordered[key] = value
        if key == "before_sha256" and pins is not None:
            ordered["catalog_revision_sha256"], ordered[HEAD_KEY] = pins
    return ordered


def _set_records(operations: list[dict]) -> list[str]:
    return list(dict.fromkeys(o["record"] for o in operations if o["op"] == "set"))


def _added(operations: list[dict]) -> list[str]:
    return [o["record"]["identifier"] for o in operations if o["op"] == "add_record"]


def _policy(
    stage: Path,
    material: str,
    reference: bytes,
    revision_id: str | None,
    remove: str | None,
) -> None:
    """Move the material's rights pin and the catalog hash, in the stage."""
    record_path = _paths(material)[0]
    policy = json.loads((stage / POLICY).read_bytes())
    policy["material_sha256"][record_path] = _sha(reference)
    policy["catalog_scientific_sha256"] = fingerprint(scientific_snapshot(stage))
    if revision_id is not None:
        pinned = policy.setdefault(
            "catalog_revisions", {"ledger": LEDGER, "materials": {}}
        )
        pinned["materials"][material] = [
            *pinned["materials"].get(material, []),
            revision_id,
        ]
        pinned["materials"] = dict(sorted(pinned["materials"].items()))
    if remove is not None:
        pinned = policy["catalog_revisions"]["materials"]
        pinned[material] = [i for i in pinned.get(material, []) if i != remove]
        if not pinned[material]:
            del pinned[material]
        if not pinned:
            del policy["catalog_revisions"]
    (stage / POLICY).write_bytes(_json_bytes(policy))
    try:
        outputs = build_review(stage)
    except (ValueError, KeyError, TypeError) as error:
        raise Refused("rights_review", f"the rights review cannot be built: {error!r}")
    for name, content in outputs.items():
        (stage / name).parent.mkdir(parents=True, exist_ok=True)
        (stage / name).write_bytes(content)


def _material_records(snapshot: dict, material: str) -> list:
    for document in snapshot["documents"]:
        if document["identifier"] == material:
            return document["eos_records"]
    return []


def plan_revision(request: dict, root: Path = ROOT, *, derived: bool = True) -> Plan:
    """Compute and verify one revision in a staging copy; write nothing.

    ``derived=False`` skips the documentation generators (tests only).
    """
    check_request(request)
    material = request["material"]
    kind = request["kind"]
    ledger = read_ledger(root)
    if any(r["id"] == request["id"] for r in ledger["revisions"]):
        raise Refused("duplicate_id", f"{request['id']} is already in {LEDGER}")
    document, sidecar = _read(root, material)
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    binding = receipt["materials"].get(material)
    if binding is None:
        raise Refused("not_split", f"{material} has no audit-split binding")
    merged_before = merge_audit(document, sidecar)
    record_path, audit_path = _paths(material)
    files_before = ((root / record_path).read_bytes(), (root / audit_path).read_bytes())
    reference_before = pre_split_bytes(material, root=root)
    validity_before = {
        r["identifier"]: _validity(r) for r in merged_before["eos_records"]
    }
    # Records another receipt pins cannot be revised; say so before any work.
    if kind == "set_fields":
        _check_pins(root, [c.get("record") for c in request["changes"]])
    # The request against the files, before anything is staged.
    if kind == "add_record":
        record = _new_record(material, request["record"], document, root)
        forward = [{"op": "add_record", "record": record}]
        _apply_to_files(document, sidecar, forward)
        touched = [record["identifier"]]
    else:
        forward = _set_fields(document, request["changes"])
        touched = _set_records(forward)
    _refresh_checks(document, sidecar, set(touched))
    _merged(_files(copy.deepcopy(document), copy.deepcopy(sidecar)))
    with _stage(root) as stage:
        snapshot_before = scientific_snapshot(stage)
        refits = _refit_ledger(stage)
        saved_rows = {r["record_identifier"]: r for r in refits["records"]}
        if kind == "set_fields":
            # The saved rows must be what the validator gives today.
            fresh = _fresh_rows(stage, material)
            from scripts.validate_primary_eos_refits import check_ledger_values

            for identifier in touched:
                try:
                    check_ledger_values(saved_rows[identifier], fresh[identifier])
                except (KeyError, ValueError) as error:
                    largest = _largest_difference(
                        saved_rows.get(identifier), fresh.get(identifier)
                    )
                    raise Refused(
                        "refit_ledger_stale",
                        f"{identifier}: the stored row of {REFITS} differs from a "
                        f"fresh computation on this machine before any change "
                        f"(largest relative difference {largest}; {error})",
                    ) from None
        _refresh_checks(document, sidecar, set(touched))
        files = _files(document, sidecar)
        _merged(files)
        (stage / record_path).write_bytes(files[0])
        (stage / audit_path).write_bytes(files[1])
        rows = _fresh_rows(stage, material)
        refit_rows: dict[str, dict] = {}
        derive: dict[str, dict] = {}
        for identifier in touched:
            row = rows[identifier]
            before = saved_rows.get(identifier)
            if before is not None:
                from scripts.validate_primary_eos_refits import check_ledger_values

                try:
                    check_ledger_values(before, row)
                    row = before  # unchanged up to solver roundoff
                except ValueError:
                    pass
            if before is None or row is not before:
                refit_rows[identifier] = {"before": before, "after": row}
            derive[identifier] = row
        _check_rows(document, forward, derive, kind, stage)
        derived_ops = _derive_assessments(stage, material, document, sidecar, derive)
        if kind == "add_record":
            # The added record carries its derived assessment itself.
            derived_ops = []
        _refresh_checks(document, sidecar, set(touched))
        files = _files(document, sidecar)
        merged_after = _merged(files)
        if kind == "add_record":
            forward[0]["record"] = next(
                r for r in merged_after["eos_records"] if r["identifier"] == touched[0]
            )
        forward += derived_ops
        reverse = reverse_operations(forward)
        _check_declared(merged_before, merged_after, forward, reverse, touched)
        headline_changes = _headline_changes(merged_before, merged_after, touched)
        if headline_changes != request.get(ACKNOWLEDGE, {}):
            raise Refused(
                "check_change",
                "the regenerated check block changes "
                + (json.dumps(headline_changes) if headline_changes else "nothing")
                + f"; the request must carry exactly this as {ACKNOWLEDGE}",
            )
        changed_ranges = {
            o["record"]
            for o in forward
            if o["op"] == "set" and o["path"][0] in (*RANGE_FIELDS, "validity")
        }
        validity: dict[str, tuple[Any, Any]] = {}
        for record in merged_after["eos_records"]:
            identifier = record["identifier"]
            if identifier not in touched:
                continue
            old, new = validity_before.get(identifier), _validity(record)
            if old is not None and old != new and identifier not in changed_ranges:
                raise Refused(
                    "validity_changed",
                    f"{identifier}: the enforced range would change from {old} to "
                    f"{new} although no range field is changed",
                )
            validity[identifier] = (old, new)
        # Bookkeeping in the stage.
        (stage / record_path).write_bytes(files[0])
        (stage / audit_path).write_bytes(files[1])
        for identifier, pair in refit_rows.items():
            records = refits["records"]
            if pair["before"] is None:
                _insert_row(records, pair["after"], material)
            else:
                records[records.index(pair["before"])] = pair["after"]
        _write_refit_ledger(stage, refits)
        if kind == "add_record":
            # Settable fields never change these counts.
            _bookkeeping(stage, material)
            _check_bookkeeping(stage)
        # A derived audit value reads differently in the pre-split reference,
        # where the split puts generic statements and bound originals back.
        _reference_values(forward, merged_after, binding, document, reference_before)
        reverse = reverse_operations(forward)
        # The pre-split reference: forward operations on the current reference
        # and recovery from the new files through the revised binding agree.
        reference_document = json.loads(reference_before)
        try:
            apply_operations(reference_document["eos_records"], forward, "reference")
        except ValueError as error:
            raise Refused(
                "not_bindable",
                f"the change cannot be expressed on the pre-split reference ({error}); "
                "the audit-split binding cannot express it",
            ) from None
        layout = LAYOUTS[binding.get("layout")]
        reference_after = format_eosmat(reference_document, layout=layout).encode(
            "utf-8"
        )
        revision: dict[str, Any] = {
            **{key: copy.deepcopy(request[key]) for key in COMMON_KEYS},
            ACKNOWLEDGE: copy.deepcopy(request.get(ACKNOWLEDGE, {})),
            "forward": forward,
            "reverse": reverse,
            "reference_sha256": {
                "before": _sha(reference_before),
                "after": _sha(reference_after),
            },
        }
        # Provisional: the complete revision's digest replaces it below.
        revision["content_sha256"] = content_sha256(revision)
        order_document = copy.deepcopy(document)
        new_binding = _revised_binding(
            binding,
            order_document,
            forward,
            files,
            (_sha(reference_after), revision["content_sha256"]),
        )
        try:
            recover_pre_split(material, *files, new_binding)
        except ValueError as error:
            raise Refused(
                "not_bindable",
                f"the audit-split binding cannot restore the revised reference ({error})",
            ) from None
        new_receipt = copy.deepcopy(receipt)
        new_receipt["materials"][material] = new_binding
        new_receipt.setdefault(*SPLIT_RECEIPT_STATEMENT)
        (stage / RECEIPT_PATH).write_text(receipt_text(new_receipt), encoding="utf-8")
        # A provisional ledger entry, so the staged binding's pins resolve
        # while the derived files are built; the complete entry follows.
        provisional = copy.deepcopy(ledger)
        provisional["revisions"].append(revision)
        (stage / LEDGER).parent.mkdir(parents=True, exist_ok=True)
        (stage / LEDGER).write_text(ledger_text(provisional), encoding="utf-8")
        snapshot_after = scientific_snapshot(stage)
        _check_snapshot(snapshot_before, snapshot_after, material, forward)
        policy = json.loads((stage / POLICY).read_bytes())
        if policy["material_sha256"].get(record_path) != _sha(reference_before):
            raise Refused(
                "rights_review_stale",
                f"{POLICY} does not pin the current reference of {material}; run "
                "python scripts/review_eos_structure_rights.py --check",
            )
        if policy["catalog_scientific_sha256"] != fingerprint(snapshot_before):
            raise Refused(
                "rights_review_stale",
                f"{POLICY} catalog_scientific_sha256 is not the current snapshot",
            )
        _policy(stage, material, reference_after, request["id"], None)
        _fixtures(stage, root, material)
        if derived:
            _generate(stage)
        _load_catalog(stage, material, touched)
        previous = [r for r in ledger["revisions"] if r["material"] == material]
        files_changed = _changed_files(stage, root)
        own = {record_path, audit_path, RECEIPT_PATH, LEDGER}
        revision.update(
            {
                "document_sha256": {
                    "before": _sha(format_eosmat(merged_before).encode("utf-8")),
                    "after": _sha(format_eosmat(merged_after).encode("utf-8")),
                },
                "reference_sha256": {
                    "before": _sha(reference_before),
                    "after": _sha(reference_after),
                },
                "scientific": {
                    "records_before": fingerprint(
                        _material_records(snapshot_before, material)
                    ),
                    "records_after": fingerprint(
                        _material_records(snapshot_after, material)
                    ),
                    "added_records": {
                        r["identifier"]: fingerprint(r)
                        for r in _material_records(snapshot_after, material)
                        if r["identifier"] in _added(forward)
                    },
                    "catalog_before": fingerprint(snapshot_before),
                    "catalog_after": fingerprint(snapshot_after),
                },
                "files": {
                    record_path: {
                        "before": _sha(files_before[0]),
                        "after": _sha(files[0]),
                    },
                    audit_path: {
                        "before": _sha(files_before[1]),
                        "after": _sha(files[1]),
                    },
                },
                "split_binding": {"before": copy.deepcopy(binding)},
                "policy": {
                    "material_sha256": {
                        "before": _sha(reference_before),
                        "after": _sha(reference_after),
                    },
                    "catalog_scientific_sha256": {
                        "before": fingerprint(snapshot_before),
                        "after": fingerprint(snapshot_after),
                    },
                },
                "refit_rows": refit_rows,
                # Generated pages are listed, not hashed: later revisions,
                # undos and generator changes legitimately rewrite them.
                "written_files": sorted(r for r in files_changed if r not in own),
                "created_files": sorted(
                    r
                    for r, content in files_changed.items()
                    if r not in own and content is not None and not (root / r).exists()
                ),
                "previous_sha256": fingerprint(previous[-1]) if previous else None,
            }
        )
        from scripts.catalog_revisions import REVISION_KEYS

        revision["content_sha256"] = None
        revision = {key: revision[key] for key in REVISION_KEYS}
        revision["content_sha256"] = content_sha256(revision)
        new_binding = _revised_binding(
            binding,
            order_document,
            forward,
            files,
            (_sha(reference_after), revision["content_sha256"]),
        )
        new_receipt["materials"][material] = new_binding
        (stage / RECEIPT_PATH).write_bytes(receipt_text(new_receipt).encode("utf-8"))
        new_ledger = copy.deepcopy(ledger)
        new_ledger["revisions"].append(revision)
        (stage / LEDGER).parent.mkdir(parents=True, exist_ok=True)
        (stage / LEDGER).write_text(ledger_text(new_ledger), encoding="utf-8")
        _verify_stage(stage, root, material, snapshot_before)
        plan = Plan(revision=revision, undo=False)
        plan.files = _changed_files(stage, root)
    plan.shown = _show(revision, merged_before, merged_after, touched, validity)
    return plan


def _reference_values(
    forward: list[dict],
    merged: dict,
    binding: dict,
    document: dict,
    reference: bytes,
) -> None:
    """Add the reference-view values of derived audit operations, in place."""
    from scripts.split_eosmat_audits import restore

    derived = [o for o in forward if o["op"] == "set" and o.get("derived")]
    if not derived:
        return
    trial = _revised_binding(binding, document, forward, (b"", b""), None)
    restored = {r["identifier"]: r for r in restore(merged, trial)["eos_records"]}
    current = {r["identifier"]: r for r in json.loads(reference)["eos_records"]}
    for operation in derived:
        path = operation["path"]
        before = read_path(current[operation["record"]], path)
        after = read_path(restored[operation["record"]], path)
        if (before, after) != (operation["before"], operation["after"]):
            operation["reference_before"] = None if before is ABSENT else before
            operation["reference_after"] = None if after is ABSENT else after


def _headline_changes(before: dict, after: dict, touched: list[str]) -> dict:
    """Changes of the headline fields of existing records' ``check`` blocks."""
    old = {r["identifier"]: r.get("check", {}) for r in before["eos_records"]}
    result: dict[str, dict] = {}
    for record in after["eos_records"]:
        identifier = record["identifier"]
        if identifier not in touched or identifier not in old:
            continue
        fields = {
            field: {
                "from": old[identifier].get(field),
                "to": record["check"].get(field),
            }
            for field in HEADLINE_FIELDS
            if old[identifier].get(field) != record["check"].get(field)
        }
        if fields:
            result[identifier] = fields
    return result


def _check_declared(
    before: dict,
    after: dict,
    forward: list[dict],
    reverse: list[dict],
    touched: list[str],
) -> None:
    """The operations turn ``before`` into ``after`` and back, ``check`` aside."""
    names = set(touched)

    def without_checks(document: dict) -> dict:
        result = copy.deepcopy(document)
        for record in result["eos_records"]:
            if record["identifier"] in names:
                record.pop("check", None)
        return result

    trial = copy.deepcopy(before)
    apply_operations(trial["eos_records"], forward, "document")
    if fingerprint(without_checks(trial)) != fingerprint(without_checks(after)):
        raise Refused("undeclared_change", "the files differ beyond the operations")
    back = copy.deepcopy(after)
    apply_operations(back["eos_records"], reverse, "document")
    if fingerprint(without_checks(back)) != fingerprint(without_checks(before)):
        raise Refused(
            "not_reversible", "the reverse operations do not restore the file"
        )
    for old, new in zip(before["eos_records"], after["eos_records"]):
        if old["identifier"] not in names and old != new:
            raise Refused("undeclared_change", f"{old['identifier']} would change")


def _check_snapshot(
    before: dict, after: dict, material: str, forward: list[dict]
) -> None:
    """The scientific snapshot changes by exactly the operations."""

    def rest(snapshot: dict) -> dict:
        result = copy.deepcopy(snapshot)
        for document in result["documents"]:
            if document["identifier"] == material:
                document.pop("eos_records")
        return result

    if fingerprint(rest(before)) != fingerprint(rest(after)):
        raise Refused("undeclared_change", "the scientific snapshot changed elsewhere")
    trial = copy.deepcopy(_material_records(before, material))
    apply_operations(trial, forward, "snapshot")
    if fingerprint(trial) != fingerprint(_material_records(after, material)):
        raise Refused(
            "undeclared_change", "the scientific snapshot differs from the operations"
        )


def _verify_stage(stage: Path, root: Path, material: str, before: dict) -> None:
    """Every recovery check on the staged files; raises on the first failure."""
    from scripts.apply_check_text import verify as verify_check_text
    from scripts.migrate_eosmat_v5 import prior_material_bytes
    from scripts.split_eosmat_audits import verify_catalog

    failures = verify_ledger(stage)
    if failures:
        raise Refused("lineage", "; ".join(failures))
    try:
        verify_catalog(stage)
    except ValueError as error:
        raise Refused("lineage", str(error)) from None
    if prior_material_bytes(material, root=stage) != prior_material_bytes(
        material, root=root
    ):
        raise Refused("lineage", "the v5-migration predecessor would change")
    now = scientific_snapshot(stage)
    base = prior_revision_snapshot(before, root) or before
    if fingerprint(prior_revision_snapshot(now, stage) or now) != fingerprint(base):
        raise Refused("lineage", "the reversed scientific snapshot would change")
    outputs = build_review(stage)
    if any((stage / name).read_bytes() != content for name, content in outputs.items()):
        raise Refused("lineage", "the rights review outputs are stale")
    failures = verify_check_text(stage)[1]
    if failures:
        raise Refused("lineage", "; ".join(failures))


def _values(record: Mapping[str, Any]) -> list[str]:
    """Coefficients, errors and fixed flags of a record, for display."""
    lines = []
    for prefix, parameters in (
        ("", record.get("parameters") or {}),
        ("thermal ", (record.get("thermal") or {}).get("parameters") or {}),
    ):
        for name, item in parameters.items():
            if not isinstance(item, dict):
                continue
            error = item.get("error")
            text = f"{prefix}{name} = {item.get('value')}"
            text += f" ± {error}" if error is not None else " (no error)"
            if item.get("fixed"):
                text += ", fixed"
            lines.append(f"    {text}")
    return lines


def _show(
    revision: dict,
    before: dict,
    after: dict,
    touched: list[str],
    validity: Mapping[str, tuple[Any, Any]],
) -> list[str]:
    """What a scientist checks: the values, the enforced range and ``check``."""
    old = {r["identifier"]: r for r in before["eos_records"]}
    new = {r["identifier"]: r for r in after["eos_records"]}
    lines = []
    for identifier in touched:
        record = new[identifier]
        if identifier in old:
            lines.append(f"Record {identifier} (changed):")
            for operation in revision["forward"]:
                if operation["op"] != "set" or operation["record"] != identifier:
                    continue
                if operation.get("derived"):
                    continue
                was = (
                    json.dumps(operation["before"])
                    if operation["before_present"]
                    else "absent"
                )
                now = (
                    json.dumps(operation["after"])
                    if operation["after_present"]
                    else "absent"
                )
                lines.append(f"  {'.'.join(operation['path'])}: {was} -> {now}")
            before_check = old[identifier].get("check", {})
            for key in (*HEADLINE_FIELDS, "why"):
                if before_check.get(key) != record["check"].get(key):
                    marker = "CHECK CHANGE" if key in HEADLINE_FIELDS else "check"
                    lines.append(
                        f"  {marker} {key}: {json.dumps(before_check.get(key))} -> "
                        f"{json.dumps(record['check'].get(key))}"
                    )
        else:
            lines.append(f"Record {identifier} (added): {record.get('label')}")
            lines.append(
                f"  model {record.get('model')}"
                + (
                    f" + {record['thermal'].get('model')}"
                    if isinstance(record.get("thermal"), dict)
                    else ""
                )
                + f", temperature_ref {record.get('temperature_ref')}"
            )
            lines += _values(record)
            lines.append(
                f"  pressure range {record.get('experimental_pressure_range_gpa')} GPa, "
                f"temperature range {record.get('experimental_temperature_range_k')} K"
            )
        was, now = validity.get(identifier, (None, None))
        if was is not None and was != now:
            lines.append(f"  enforced range: {was} -> {now}")
        elif now is not None:
            lines.append(
                f"  enforced range: {now.pressure_gpa} GPa, {now.temperature_k} K"
            )
        lines.append("  check block (record file):")
        lines += [
            "    " + line
            for line in json.dumps(
                record["check"], indent=2, ensure_ascii=False
            ).splitlines()
        ]
    return lines


# Undo ---------------------------------------------------------------------------


def plan_undo(revision_id: str, root: Path = ROOT, *, derived: bool = True) -> Plan:
    """Compute the undo of a material's most recent revision; write nothing.

    Later text revisions of the same material are kept: the reverse
    operations are applied to the current files and the binding entries of a
    removed record are dropped. Without later changes, every file returns to
    its bytes before the revision.
    """
    ledger = read_ledger(root)
    revision = next((r for r in ledger["revisions"] if r["id"] == revision_id), None)
    if revision is None:
        raise Refused("unknown_revision", f"{revision_id} is not in {LEDGER}")
    material = revision["material"]
    later = [
        r["id"]
        for r in ledger["revisions"][ledger["revisions"].index(revision) + 1 :]
        if r["material"] == material
    ]
    if later:
        raise Refused("out_of_order", f"undo the later revisions {later} first")
    record_path, audit_path = _paths(material)
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    binding = receipt["materials"][material]
    unchanged = all(
        _sha((root / p).read_bytes()) == revision["files"][p]["after"]
        for p in (record_path, audit_path)
    )
    document, sidecar = _read(root, material)
    order_before = copy.deepcopy(document)
    try:
        _apply_to_files(document, sidecar, revision["reverse"])
    except ValueError as error:
        raise Refused(
            "changed_since", f"cannot reverse {revision_id}: {error}"
        ) from None
    _refresh_checks(document, sidecar, set(_set_records(revision["reverse"])))
    files = _files(document, sidecar)
    _merged(files)
    previous = [r for r in ledger["revisions"] if r["material"] == material][:-1]
    head = previous[-1] if previous else None
    new_binding = _revised_binding(
        binding,
        order_before,
        revision["reverse"],
        files,
        None
        if head is None
        else (head["reference_sha256"]["after"], head["content_sha256"]),
    )
    if unchanged:
        exact = (
            _sha(files[0]) == revision["files"][record_path]["before"]
            and _sha(files[1]) == revision["files"][audit_path]["before"]
            and new_binding == revision["split_binding"]["before"]
        )
        if not exact:
            raise Refused(
                "not_exact", "the files would not return to their earlier bytes"
            )
    try:
        reference = recover_pre_split(material, *files, new_binding)
    except ValueError as error:
        raise Refused(
            "not_exact", f"the pre-split reference would not return ({error})"
        )
    if _sha(reference) != revision["reference_sha256"]["before"]:
        raise Refused("not_exact", "the pre-split reference would not return")
    receipt["materials"][material] = new_binding
    if not any("catalog_revision_sha256" in b for b in receipt["materials"].values()):
        receipt.pop(SPLIT_RECEIPT_STATEMENT[0], None)
    new_ledger = copy.deepcopy(ledger)
    new_ledger["revisions"] = [r for r in ledger["revisions"] if r["id"] != revision_id]
    with _stage(root) as stage:
        (stage / record_path).write_bytes(files[0])
        (stage / audit_path).write_bytes(files[1])
        (stage / RECEIPT_PATH).write_text(receipt_text(receipt), encoding="utf-8")
        (stage / LEDGER).write_text(ledger_text(new_ledger), encoding="utf-8")
        refits = _refit_ledger_raw(stage)
        for identifier, pair in revision["refit_rows"].items():
            rows = refits["records"]
            index = next(
                (i for i, r in enumerate(rows) if r["record_identifier"] == identifier),
                None,
            )
            if index is None:
                raise Refused(
                    "changed_since", f"the refit-ledger row of {identifier} is gone"
                )
            current = rows[index]
            if pair["before"] is None:
                del rows[index]  # the record leaves with its row
                continue
            if _without_text(current) != _without_text(pair["after"]):
                raise Refused(
                    "changed_since",
                    f"the refit-ledger row of {identifier} changed after "
                    f"{revision_id} beyond its text copies",
                )
            # Later text revisions keep their copies (label, notes, usage).
            restored = copy.deepcopy(pair["before"])
            for key in TEXT_ROW_FIELDS:
                if key in current:
                    restored[key] = copy.deepcopy(current[key])
                else:
                    restored.pop(key, None)
            if list(pair["before"]) == sorted(pair["before"]):
                restored = dict(sorted(restored.items()))
            rows[index] = restored
        _write_refit_ledger(stage, refits)
        added = _added(revision["forward"])
        if added:
            _bookkeeping(stage, material)
            _check_bookkeeping(stage)
        _policy(stage, material, reference, None, revision_id)
        _fixtures(stage, root, material)
        if derived:
            _generate(stage)
            # Generated files the revision created for its record and that no
            # generator removes (coefficient figures) leave with the record.
            for relative in revision["created_files"]:
                target = stage / relative
                if target.exists() and any(i in relative for i in added):
                    target.unlink()
        _load_catalog(stage, material, _set_records(revision["forward"]))
        records = _material_records(scientific_snapshot(stage), material)
        if fingerprint(records) != revision["scientific"]["records_before"]:
            raise Refused("not_exact", "the scientific snapshot would not return")
        failures = verify_ledger(stage)
        if failures:
            raise Refused("lineage", "; ".join(failures))
        plan = Plan(revision=revision, undo=True)
        plan.files = _changed_files(stage, root)
    plan.shown = [
        f"Undo {revision_id} ({revision['kind']}) of {material}"
        + ("" if unchanged else "; later text revisions of the material are kept"),
        *(f"  removes record {i}" for i in added),
        *(
            f"  {o['record']}/{'.'.join(o['path'])} returns to "
            + (json.dumps(o["after"]) if o["after_present"] else "absent")
            for o in revision["reverse"]
            if o["op"] == "set" and not o.get("derived")
        ),
    ]
    return plan


#: Refit-ledger row fields the text tools own (copies of record text).
TEXT_ROW_FIELDS = (
    "label",
    "source_notes",
    "usage_recommendation",
    "recommendation_reason",
)


def _without_text(row: Mapping[str, Any]) -> dict:
    return {k: v for k, v in row.items() if k not in TEXT_ROW_FIELDS}


def _refit_ledger_raw(stage: Path) -> dict:
    return json.loads((stage / REFITS).read_bytes())


# Writing --------------------------------------------------------------------------


def apply_plan(plan: Plan, root: Path = ROOT) -> None:
    """Write a verified plan; on any failure, put every file back."""
    originals = {
        relative: (root / relative).read_bytes() if (root / relative).exists() else None
        for relative in plan.files
    }
    try:
        for relative, content in plan.files.items():
            path = root / relative
            if content is None:
                path.unlink(missing_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
        failures = verify_ledger(root, [plan.revision["material"]])
        if failures:
            raise ValueError("; ".join(failures))
        for name, content in build_review(root).items():
            if (root / name).read_bytes() != content:
                raise ValueError(f"{name} is stale after writing")
    except Exception as error:
        for relative, content in originals.items():
            try:
                if content is None:
                    (root / relative).unlink(missing_ok=True)
                else:
                    (root / relative).write_bytes(content)
            except OSError:
                pass
        if isinstance(error, OSError):
            raise Refused(
                "not_writable", f"{error.filename}: {error.strerror}; nothing changed"
            ) from None
        raise


def verify(root: Path = ROOT) -> list[str]:
    """The whole ledger against the catalog (see ``verify_ledger``).

    With revisions, also that the manifest and the primary-source audit count
    the catalog.
    """
    failures = verify_ledger(root)
    if read_ledger(root)["revisions"]:
        try:
            _check_bookkeeping(root)
        except Refused as error:
            failures.append(str(error))
    return failures


def _report(plan: Plan, check: bool) -> None:
    revision = plan.revision
    verb = "Would" if check else "Did"
    action = "undo" if plan.undo else "apply"
    print(
        f"{verb} {action} {revision['id']} ({revision['kind']}) for "
        f"{revision['material']}."
    )
    for line in plan.shown:
        print(line)
    groups: dict[str, list[str]] = {}
    for relative, content in plan.files.items():
        if (
            relative.startswith(("peritheos/data/materials/", "peritheos/data/audits/"))
            and relative != MANIFEST
        ):
            group = "Record file and audit sidecar"
        elif relative in (RECEIPT_PATH, LEDGER, POLICY):
            group = "Receipts"
        else:
            group = "Derived files"
        groups.setdefault(group, []).append(
            relative + (" (removed)" if content is None else "")
        )
    for group, paths in groups.items():
        print(f"{group}:")
        for path in paths:
            print(f"  {path}")
    print(
        "No hand edits are needed except tests that pin a material's record "
        "count. Run the full checks in CONTRIBUTING.md."
    )


def main(argv: list[str] | None = None, root: Path = ROOT) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("request", nargs="?", type=Path, help="request JSON file")
    parser.add_argument("--check", action="store_true", help="dry run")
    parser.add_argument("--undo", metavar="ID", help="undo this revision")
    parser.add_argument("--verify", action="store_true", help="verify the ledger")
    arguments = parser.parse_args(argv)
    try:
        if arguments.verify:
            failures = verify(root)
            for failure in failures:
                print(f"FAILED {failure}")
            count = len(read_ledger(root)["revisions"])
            print(
                f"Verified {count} catalog revisions: "
                f"{'ok' if not failures else f'{len(failures)} failures'}."
            )
            return 1 if failures else 0
        if (arguments.request is None) == (arguments.undo is None):
            parser.error("give a request file or --undo ID")
        plan = (
            plan_undo(arguments.undo, root)
            if arguments.undo
            else plan_revision(read_request(arguments.request), root)
        )
        if not arguments.check:
            apply_plan(plan, root)
    except Refused as error:
        print(f"REFUSED {error}; nothing was written.")
        return 1
    except (OSError, ValueError) as error:
        print(f"REFUSED {type(error).__name__}: {error}; nothing was written.")
        return 1
    _report(plan, arguments.check)
    return 0


if __name__ == "__main__":
    sys.exit(main())


__all__ = [
    "ALLOW_LIST",
    "REQUIRED_RECORD_KEYS",
    "Plan",
    "Refused",
    "apply_plan",
    "check_request",
    "main",
    "plan_revision",
    "plan_undo",
    "read_request",
    "verify",
]
