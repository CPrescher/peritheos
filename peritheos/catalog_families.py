"""Curated material families for catalog browsing, without numerical models."""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import cache
from importlib import resources
from types import MappingProxyType
from typing import TYPE_CHECKING

from peritheos.errors import MaterialError, MaterialLookupError

if TYPE_CHECKING:
    from peritheos.materials import Material


def _valid_family_id(value: object) -> bool:
    return (
        isinstance(value, str)
        and re.fullmatch(r"[a-z0-9]+(?:_[a-z0-9]+)*", value) is not None
    )


@dataclass(frozen=True)
class MaterialFamily:
    """A browsing category; its formula describes scope, not an executable sample."""

    identifier: str
    name: str
    description: str
    formula: str | None = None

    def __post_init__(self) -> None:
        if not _valid_family_id(self.identifier):
            raise MaterialError("Family identifier must be lower-snake-case")
        for key in ("name", "description"):
            value = getattr(self, key)
            if not isinstance(value, str) or not value.strip():
                raise MaterialError(f"Family {key} must be a non-empty string")
        if self.formula is not None and (
            not isinstance(self.formula, str) or not self.formula.strip()
        ):
            raise MaterialError("Family formula must be a non-empty string or None")


@dataclass(frozen=True)
class MaterialGroup:
    """Matching family members, or one standalone material when family is None.

    The members are the supplied materials, not the family's complete catalog
    membership. No member or EOS is selected as a default.
    """

    family: MaterialFamily | None
    materials: tuple[Material, ...]


def _parse_family_registry(document: object) -> Mapping[str, MaterialFamily]:
    if not isinstance(document, list):
        raise MaterialError("Material family registry must be a JSON array")
    families: dict[str, MaterialFamily] = {}
    for entry in document:
        if not isinstance(entry, dict):
            raise MaterialError("Material family entries must be objects")
        try:
            family = MaterialFamily(**entry)
        except TypeError as error:
            raise MaterialError(f"Invalid material family entry: {error}") from error
        if family.identifier in families:
            raise MaterialError(f"Duplicate material family {family.identifier!r}")
        families[family.identifier] = family
    return MappingProxyType(families)


@cache
def _family_registry() -> Mapping[str, MaterialFamily]:
    resource = resources.files("peritheos.data").joinpath("material-families.json")
    return _parse_family_registry(json.loads(resource.read_text(encoding="utf-8")))


def list_material_families() -> tuple[MaterialFamily, ...]:
    """Return bundled family definitions in stable identifier order."""
    registry = _family_registry()
    return tuple(registry[key] for key in sorted(registry))


def get_material_family(identifier: str) -> MaterialFamily:
    """Look up a family by exact ID, separately from material IDs and aliases."""
    if not _valid_family_id(identifier):
        raise MaterialError("Family identifier must be lower-snake-case")
    try:
        return _family_registry()[identifier]
    except KeyError:
        raise MaterialLookupError(
            f"Unknown material family {identifier!r}.",
            operation="lookup_material_family",
            field="identifier",
            context={"identifier": identifier},
        ) from None


def group_materials(materials: Iterable[Material]) -> tuple[MaterialGroup, ...]:
    """Group a material listing or search result without adding other members.

    Groups and members are sorted by stable identifiers. Unknown families fall
    back to standalone entries, preserving imported materials. Duplicate
    material identifiers raise MaterialError instead of silently losing data.
    """
    registry = _family_registry()
    grouped: dict[tuple[str, str], list[Material]] = {}
    seen: set[str] = set()
    for material in materials:
        if material.identifier in seen:
            raise MaterialError(
                f"Duplicate material identifier {material.identifier!r}"
            )
        seen.add(material.identifier)
        family = registry.get(material.family_id) if material.family_id else None
        key = (
            (family.identifier, "family")
            if family is not None
            else (material.identifier, "material")
        )
        grouped.setdefault(key, []).append(material)
    return tuple(
        MaterialGroup(
            family=registry[identifier] if kind == "family" else None,
            materials=tuple(
                sorted(grouped[(identifier, kind)], key=lambda m: m.identifier)
            ),
        )
        for identifier, kind in sorted(grouped)
    )


__all__ = [
    "MaterialFamily",
    "MaterialGroup",
    "get_material_family",
    "group_materials",
    "list_material_families",
]
