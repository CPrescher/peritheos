"""Canonical executable material catalog.

The bundled ``.eosmat`` resources are the only source used to construct the
canonical catalog. The smaller, pre-0.7 pressure-standard collection remains
available only as an isolated compatibility lookup in :mod:`peritheos.materials`.
Search implementation lives in :mod:`peritheos.catalog_search` so loading,
indexing, and discovery can evolve independently.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from difflib import get_close_matches
from functools import cache
from types import MappingProxyType
from typing import TypeVar

from peritheos.eosmat import get_material_document, list_material_documents
from peritheos.errors import MaterialError, MaterialLookupError
from peritheos.materials import EOSRecord, Material

_CatalogObject = TypeVar("_CatalogObject", Material, EOSRecord)


@dataclass(frozen=True)
class _CatalogIndex:
    materials: Mapping[str, Material]
    records: Mapping[str, EOSRecord]
    material_aliases: Mapping[str, Material]
    record_aliases: Mapping[str, EOSRecord]
    material_by_record: Mapping[str, Material]


def _index_aliases(
    *,
    kind: str,
    objects: Mapping[str, _CatalogObject],
) -> Mapping[str, _CatalogObject]:
    aliases: dict[str, _CatalogObject] = {}
    for item in objects.values():
        for alias in item.aliases:
            if alias in objects:
                raise MaterialError(
                    f"Bundled {kind} alias {alias!r} collides with an identifier"
                )
            previous = aliases.get(alias)
            if previous is not None:
                raise MaterialError(
                    f"Duplicate bundled {kind} alias {alias!r} for "
                    f"{previous.identifier!r} and {item.identifier!r}"
                )
            aliases[alias] = item
    return MappingProxyType(aliases)


@cache
def _catalog_index() -> _CatalogIndex:
    materials: dict[str, Material] = {}
    records: dict[str, EOSRecord] = {}
    material_by_record: dict[str, Material] = {}
    for document_identifier in list_material_documents():
        material = Material.from_eosmat(get_material_document(document_identifier))
        if material.identifier != document_identifier:
            raise MaterialError(
                f"Bundled material file {document_identifier!r} contains identifier "
                f"{material.identifier!r}"
            )
        if material.identifier in materials:
            raise MaterialError(
                f"Duplicate bundled material identifier {material.identifier!r}"
            )
        materials[material.identifier] = material
        for record in material.eos_records:
            if record.identifier in records:
                raise MaterialError(
                    f"Duplicate bundled EOS record identifier {record.identifier!r}"
                )
            records[record.identifier] = record
            material_by_record[record.identifier] = material

    material_mapping = MappingProxyType(materials)
    record_mapping = MappingProxyType(records)
    return _CatalogIndex(
        materials=material_mapping,
        records=record_mapping,
        material_aliases=_index_aliases(kind="material", objects=material_mapping),
        record_aliases=_index_aliases(kind="EOS record", objects=record_mapping),
        material_by_record=MappingProxyType(material_by_record),
    )


def _compatibility_materials() -> Mapping[str, Material]:
    # Imported lazily to keep canonical construction independent from the
    # historical convenience-record objects.
    from peritheos.materials import _MATERIAL_CATALOG

    return _MATERIAL_CATALOG


def _compatibility_records() -> Mapping[str, EOSRecord]:
    from peritheos.catalog_compat import LEGACY_RECORD_AUDIT
    from peritheos.materials import _EOS_RECORD_CATALOG

    if _EOS_RECORD_CATALOG.keys() != LEGACY_RECORD_AUDIT.keys():
        raise MaterialError("Historical EOS record audit manifest is out of sync")
    return _EOS_RECORD_CATALOG


def _suggestion(identifier: str, choices: Iterable[str]) -> str:
    matches = get_close_matches(identifier, sorted(set(choices)), n=3, cutoff=0.5)
    return f" Did you mean: {', '.join(matches)}?" if matches else ""


def get_material(identifier: str) -> Material:
    """Return a bundled executable material by stable identifier or exact alias.

    Canonical identifiers take precedence. Pre-0.7 generated identifiers are
    retained as compatibility lookups when they do not collide with a
    canonical identifier.
    """
    index = _catalog_index()
    material = index.materials.get(identifier) or index.material_aliases.get(identifier)
    if material is not None:
        return material
    compatibility = _compatibility_materials()
    if identifier in compatibility:
        return compatibility[identifier]
    choices = (*index.materials, *index.material_aliases, *compatibility)
    raise MaterialLookupError(
        f"Unknown material {identifier!r}." + _suggestion(identifier, choices),
        operation="lookup_material",
        field="identifier",
        context={"identifier": identifier},
    )


def get_eos_record(identifier: str) -> EOSRecord:
    """Return a bundled executable EOS record by identifier or exact alias."""
    index = _catalog_index()
    record = index.records.get(identifier) or index.record_aliases.get(identifier)
    if record is not None:
        return record
    compatibility = _compatibility_records()
    if identifier in compatibility:
        return compatibility[identifier]
    choices = (*index.records, *index.record_aliases, *compatibility)
    raise MaterialLookupError(
        f"Unknown EOS record {identifier!r}." + _suggestion(identifier, choices),
        operation="lookup_eos_record",
        field="identifier",
        context={"identifier": identifier},
    )


def list_materials(*, formula: str | None = None) -> tuple[Material, ...]:
    """List all 115 bundled executable materials in identifier order."""
    materials = tuple(
        _catalog_index().materials[key] for key in sorted(_catalog_index().materials)
    )
    if formula is None:
        return materials
    formula_key = formula.strip().casefold()
    return tuple(item for item in materials if item.formula.casefold() == formula_key)


def list_eos_records(*, formula: str | None = None) -> tuple[EOSRecord, ...]:
    """List all 150 bundled executable records in identifier order."""
    records = tuple(
        _catalog_index().records[key] for key in sorted(_catalog_index().records)
    )
    if formula is None:
        return records
    formula_key = formula.strip().casefold()
    return tuple(item for item in records if item.material.casefold() == formula_key)


# Re-export discovery here for users who naturally look for it beside lookup.
# The import is deliberately last: catalog_search imports the listing functions.
from peritheos.catalog_search import (  # noqa: E402, I001
    RangeQuery,
    RangeSemantics,
    ValidationStatus,
    search_eos_records,
    search_materials,
)


__all__ = [
    "RangeQuery",
    "RangeSemantics",
    "ValidationStatus",
    "get_eos_record",
    "get_material",
    "list_eos_records",
    "list_materials",
    "search_eos_records",
    "search_materials",
]
