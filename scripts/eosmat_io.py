"""Version boundary for repository catalog maintenance.

Existing scientific audit algorithms use the historical dictionary layout.
Read either file version into that view, then write human-readable native v4.
This does not run an audit, refit data or infer scientific classifications.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from peritheos._eosmat_format import format_eosmat
from peritheos.eosmat import convert_eosmat, load_eosmat


def read_material(path: Path) -> dict[str, Any]:
    """Read scientific content in the stable v3-shaped maintenance view."""
    return convert_eosmat(load_eosmat(path), version=3)


def native_material(document: Mapping[str, Any]) -> dict[str, Any]:
    """Convert to inline v4, retaining explicit material names and units."""
    legacy = convert_eosmat(document, version=3)
    native = convert_eosmat(document, version=4, inline=True)
    native.setdefault("name", legacy["name"])
    native.setdefault("units", legacy["units"])
    return native


def render_material(document: Mapping[str, Any]) -> str:
    """Render a complete native material without changing its scientific content."""
    return format_eosmat(native_material(document))


def write_material(path: Path, document: Mapping[str, Any]) -> None:
    """Write native v4 after conversion and validation have succeeded."""
    path.write_text(render_material(document), encoding="utf-8")
