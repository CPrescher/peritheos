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
    extra = sorted(set(request) - allowed)
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
import json, sys
from peritheos import get_material, list_eos_records, list_materials
from peritheos.eosmat import (
    get_material_document, list_eos_record_documents, list_material_documents,
    validate_pressure_calibration_references,
)
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
    loaded = Material.from_eosmat(document, record_identifiers=[identifier])
    record = loaded.get_eos_record(identifier)
    if not record.is_hugoniot:
        record.pressure(record.reference_volume * 0.95, record.reference_temperature)
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
    from scripts.standardize_reconstruction_metadata import assessment_for_record

    merged = merge_audit(document, sidecar)
    records = {r["identifier"]: r for r in merged["eos_records"]}
    operations = []
    for identifier, row in rows.items():
        record = records[identifier]
        current = record.get("evidence", {}).get("reproduction", {}).get("assessment")
        derived = assessment_for_record(record, row)
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
                    raise Refused(
                        "refit_ledger_stale",
                        f"{identifier}: the saved row differs from the validator "
                        f"before this change ({error})",
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
            "forward": forward,
            "reverse": reverse,
        }
        revision["content_sha256"] = content_sha256(revision)
        new_binding = _revised_binding(
            binding,
            document,
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
        revision["reference_sha256"] = {
            "before": _sha(reference_before),
            "after": _sha(reference_after),
        }
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
                "split_binding": {
                    "before": copy.deepcopy(binding),
                    "after_sha256": fingerprint(new_binding),
                },
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
                "derived_files": {
                    relative: {
                        "before": _sha((root / relative).read_bytes())
                        if (root / relative).exists()
                        else None,
                        "after": None if content is None else _sha(content),
                    }
                    for relative, content in files_changed.items()
                    if relative not in own
                },
                "previous_sha256": fingerprint(previous[-1]) if previous else None,
            }
        )
        from scripts.catalog_revisions import REVISION_KEYS

        revision = {key: revision[key] for key in REVISION_KEYS}
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
            if old[identifier].get("check") != record.get("check"):
                was = old[identifier].get("check", {}).get("reproduction")
                lines.append(
                    f"  reproduction: {was} -> {record['check']['reproduction']}"
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
    unchanged = (
        all(
            _sha((root / p).read_bytes()) == revision["files"][p]["after"]
            for p in (record_path, audit_path)
        )
        and fingerprint(binding) == revision["split_binding"]["after_sha256"]
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
            if index is None or rows[index] != pair["after"]:
                if pair["before"] is None and index is not None:
                    del rows[index]  # later text revisions may have edited it
                    continue
                raise Refused(
                    "changed_since",
                    f"the refit-ledger row of {identifier} changed after {revision_id}; "
                    "return its usage recommendation with scripts/apply_check_text.py "
                    "first",
                )
            if pair["before"] is None:
                del rows[index]
            else:
                rows[index] = pair["before"]
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
            for relative, pair in revision["derived_files"].items():
                target = stage / relative
                if (
                    pair["before"] is None
                    and target.exists()
                    and _sha(target.read_bytes()) == pair["after"]
                    and any(i in relative for i in added)
                ):
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
