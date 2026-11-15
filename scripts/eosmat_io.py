"""Version boundary for repository catalog maintenance.

Existing scientific audit algorithms use the historical dictionary layout.
Read either file version into that view, then write human-readable native v5.
This does not run an audit, refit data or infer scientific classifications.

A bundled material (a file in a ``materials`` directory beside an ``audits``
directory) is read with its audit sidecar merged and written back split:
record file plus ``audits/<id>.audit.json``, with the ``check`` headline
re-derived from the evidence. Other paths are read and written inline.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from peritheos._eosmat_audit import split_document
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos.eosmat import convert_eosmat, load_eosmat


def _audit_directory(path: Path) -> Path | None:
    path = Path(path)
    audits = path.parent.parent / "audits"
    if path.parent.name == "materials" and audits.is_dir():
        return audits
    return None


def read_native_material(path: Path) -> dict[str, Any]:
    """Read a complete native document, merging a bundled audit sidecar."""
    return load_eosmat(path)


def read_material(path: Path) -> dict[str, Any]:
    """Read scientific content in the stable v3-shaped maintenance view."""
    return convert_eosmat(
        read_native_material(path), version=3, legacy_uncertainty_units=False
    )


def native_material(document: Mapping[str, Any]) -> dict[str, Any]:
    """Convert to inline v5, retaining explicit material names and units."""
    legacy = convert_eosmat(document, version=3, legacy_uncertainty_units=False)
    native = convert_eosmat(document, version=5, inline=True)
    native.setdefault("name", legacy["name"])
    native.setdefault("units", legacy["units"])
    return native


def render_material(document: Mapping[str, Any]) -> str:
    """Render a complete native material without changing its scientific content."""
    return format_eosmat(native_material(document))


def write_material(path: Path, document: Mapping[str, Any]) -> None:
    """Write native v5 after conversion and validation have succeeded.

    Bundled materials are written as a record file plus audit sidecar.
    """
    audits = _audit_directory(path)
    if audits is None:
        path.write_bytes(render_material(document).encode("utf-8"))
        return
    native = native_material(document)
    for record in native.get("eos_records", []):
        record.pop("check", None)
    records, sidecar = split_document(native)
    record_bytes = format_eosmat(records).encode("utf-8")
    audit_bytes = format_audit(sidecar).encode("utf-8")
    path.write_bytes(record_bytes)
    (audits / f"{records['identifier']}.audit.json").write_bytes(audit_bytes)
