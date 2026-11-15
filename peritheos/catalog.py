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
from difflib import SequenceMatcher, get_close_matches
from functools import cache
from types import MappingProxyType
from typing import TypeVar

from peritheos.catalog_families import (
    MaterialFamily,
    MaterialGroup,
    get_material_family,
    group_materials,
    list_material_families,
)
from peritheos.eosmat import (
    convert_eosmat,
    get_material_document,
    list_material_documents,
)
from peritheos.errors import MaterialError, MaterialLookupError
from peritheos.materials import EOSRecord, Material, _material_from_eosmat

_CatalogObject = TypeVar("_CatalogObject", Material, EOSRecord)


@dataclass(frozen=True)
class _CatalogIndex:
    materials: Mapping[str, Material]
    records: Mapping[str, EOSRecord]
    material_aliases: Mapping[str, Material]
    record_aliases: Mapping[str, EOSRecord]
    ambiguous_material_aliases: Mapping[str, tuple[str, ...]]
    ambiguous_record_aliases: Mapping[str, tuple[str, ...]]
    material_by_record: Mapping[str, Material]


def _index_aliases(
    *,
    kind: str,
    objects: Mapping[str, _CatalogObject],
) -> tuple[
    Mapping[str, _CatalogObject],
    Mapping[str, tuple[str, ...]],
]:
    alias_targets: dict[str, list[_CatalogObject]] = {}
    for item in objects.values():
        for alias in item.aliases:
            if alias in objects:
                raise MaterialError(
                    f"Bundled {kind} alias {alias!r} collides with an identifier"
                )
            alias_targets.setdefault(alias, []).append(item)
    aliases = {
        alias: targets[0]
        for alias, targets in alias_targets.items()
        if len(targets) == 1
    }
    ambiguous = {
        alias: tuple(sorted(item.identifier for item in targets))
        for alias, targets in alias_targets.items()
        if len(targets) > 1
    }
    return MappingProxyType(aliases), MappingProxyType(ambiguous)


@cache
def _catalog_index() -> _CatalogIndex:
    materials: dict[str, Material] = {}
    records: dict[str, EOSRecord] = {}
    material_by_record: dict[str, Material] = {}
    families = {family.identifier for family in list_material_families()}
    for document_identifier in list_material_documents():
        # The record file alone: audit sidecars load lazily per material.
        source_document = get_material_document(
            document_identifier, format_version=None, include_audit=False
        )
        if (
            source_document.get("family_id") is not None
            and source_document["family_id"] not in families
        ):
            raise MaterialError(
                f"Bundled material {document_identifier!r} references unknown "
                f"family {source_document['family_id']!r}"
            )
        document = convert_eosmat(
            source_document, version=3, legacy_uncertainty_units=False
        )
        # Deferred records remain document-only. An explicitly marked, nondefault
        # not-reproduced implementation can be selected for exploration; its
        # validation status is preserved. It cannot replace the default EOS.
        executable_ids = [
            record["identifier"]
            for record in document["eos_records"]
            if record["scientific_validation"]["status"] == "primary_source_validated"
            or (
                record["scientific_validation"]["status"] == "not_reproduced"
                and record.get("catalog_access") == "explicit_selection"
                and record.get("default") is not True
                and record.get("default_for") is None
            )
        ]
        if not executable_ids:
            continue
        material = _material_from_eosmat(
            source_document,
            record_identifiers=executable_ids,
            require_primary_validation=False,
            bundled_audit=True,
        )
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
    material_aliases, ambiguous_material_aliases = _index_aliases(
        kind="material", objects=material_mapping
    )
    record_aliases, ambiguous_record_aliases = _index_aliases(
        kind="EOS record", objects=record_mapping
    )
    return _CatalogIndex(
        materials=material_mapping,
        records=record_mapping,
        material_aliases=material_aliases,
        record_aliases=record_aliases,
        ambiguous_material_aliases=ambiguous_material_aliases,
        ambiguous_record_aliases=ambiguous_record_aliases,
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


_MAX_LISTED_CANDIDATES = 8

_MAX_SUGGESTIONS = 5


@cache
def _document_only_materials() -> frozenset[str]:
    """Bundled material cards without an executable EOS record."""
    return frozenset(list_material_documents()) - set(_catalog_index().materials)


def _material_label(material: Material) -> str:
    phase = " ".join(material.phase.split())
    if len(phase) > 40:
        phase = phase[:37].rstrip() + "..."
    return f"{material.identifier} ({phase})"


def _ambiguous_material(
    query: str, candidates: Iterable[Material], reason: str, subject: str = "material"
) -> MaterialLookupError:
    unique = {item.identifier: item for item in candidates}
    ordered = [unique[key] for key in sorted(unique)]
    listed = ", ".join(
        _material_label(item) for item in ordered[:_MAX_LISTED_CANDIDATES]
    )
    more = len(ordered) - _MAX_LISTED_CANDIDATES
    if more > 0:
        listed += f", and {more} more"
    return MaterialLookupError(
        f"Ambiguous {subject} {query!r}: {reason} {len(ordered)} materials: "
        f"{listed}. Use one of the identifiers.",
        operation="lookup_material",
        field="identifier",
        context={
            "identifier": query,
            "matches": tuple(item.identifier for item in ordered),
        },
    )


def _unique_or_ambiguous(
    query: str, candidates: list[Material], reason: str
) -> Material | None:
    unique = {item.identifier: item for item in candidates}
    if len(unique) == 1:
        return next(iter(unique.values()))
    if unique:
        raise _ambiguous_material(query, unique.values(), reason)
    return None


def _material_names(material: Material) -> tuple[str, ...]:
    return (*material.aliases, *((material.name,) if material.name else ()))


def _material_suggestions(query: str, index: _CatalogIndex) -> list[str]:
    """Rank identifier/name/alias/formula keys by prefix, containment, spelling."""
    folded = query.casefold()
    keys: dict[str, list[Material]] = {}
    for material in index.materials.values():
        for key in (material.identifier, *_material_names(material), material.formula):
            keys.setdefault(key.casefold(), []).append(material)
    scored: list[tuple[float, str]] = []
    for key in keys:
        ratio = SequenceMatcher(None, folded, key).ratio()
        contained = len(folded) >= 3 and folded in key
        if ratio < 0.6 and not contained:
            continue
        bonus = 0.5 if contained and key.startswith(folded) else 0.25 * contained
        scored.append((ratio + bonus, key))
    ranked = [
        material
        for _, key in sorted(scored, key=lambda item: (-item[0], item[1]))
        for material in keys[key]
    ]
    labels = dict.fromkeys(
        f"{item.identifier} ({item.name})" if item.name else item.identifier
        for item in ranked
    )
    return list(labels)[:_MAX_SUGGESTIONS]


def get_material(identifier: str) -> Material:
    """Return a bundled executable material by identifier, alias, name or formula.

    Resolution is ordered and never guesses between candidates:

    1. exact identifier or exact alias (an alias shared by several materials
       raises); then pre-0.7 compatibility identifiers;
    2. case-insensitive identifier, together with a case-sensitive formula
       (``"Au"``, ``"CO"`` versus ``"Co"``);
    3. case-insensitive alias or name (``"periclase"``, ``"GOLD"``);
    4. case-insensitive formula, only for an all-lowercase query (``"au"``).

    A step matching exactly one material returns it; a step matching several
    raises :class:`~peritheos.errors.MaterialLookupError` listing the candidate
    identifiers and phases. The identifier of a bundled card without an
    executable record is never redirected to another material. Unknown input
    raises with ranked suggestions.
    """
    index = _catalog_index()
    material = index.materials.get(identifier) or index.material_aliases.get(identifier)
    if material is not None:
        return material
    if identifier in index.ambiguous_material_aliases:
        matches = index.ambiguous_material_aliases[identifier]
        raise _ambiguous_material(
            identifier,
            (index.materials[item] for item in matches),
            "shared by",
            subject="material alias",
        )
    compatibility = _compatibility_materials()
    if identifier in compatibility:
        return compatibility[identifier]
    if not isinstance(identifier, str):
        raise MaterialLookupError(
            f"Material identifier must be a string, not {type(identifier).__name__}",
            operation="lookup_material",
            field="identifier",
            context={"identifier": repr(identifier)},
        )
    query = identifier.strip()
    folded = query.casefold()
    if folded in _document_only_materials():
        # Never redirect a material card's own identifier to another material.
        raise MaterialLookupError(
            f"Unknown material {identifier!r}: the bundled card {folded!r} has no "
            f"executable EOS record; get_material_document({folded!r}) returns it.",
            operation="lookup_material",
            field="identifier",
            context={"identifier": identifier},
        )
    materials = tuple(index.materials.values())
    for candidates, reason in (
        (
            [
                item
                for item in materials
                if item.identifier == folded or item.formula == query
            ],
            "the identifier or formula matches",
        ),
        (
            [
                item
                for item in materials
                if any(name.casefold() == folded for name in _material_names(item))
            ],
            "the name or alias matches",
        ),
        (
            [item for item in materials if item.formula.casefold() == folded]
            if query and query == query.lower()
            else [],
            "the formula matches",
        ),
    ):
        material = _unique_or_ambiguous(identifier, candidates, reason)
        if material is not None:
            return material
    suggestions = _material_suggestions(query, index)
    hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
    raise MaterialLookupError(
        f"Unknown material {identifier!r}.{hint} Browse with list_materials() "
        "or search_materials(text).",
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
    if identifier in index.ambiguous_record_aliases:
        matches = index.ambiguous_record_aliases[identifier]
        raise MaterialLookupError(
            f"Ambiguous EOS record alias {identifier!r}; matches: "
            f"{', '.join(matches)}.",
            operation="lookup_eos_record",
            field="identifier",
            context={"identifier": identifier, "matches": matches},
        )
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


def list_materials(
    *, formula: str | None = None, family_id: str | None = None
) -> tuple[Material, ...]:
    """List executable materials in identifier order, optionally within a family.

    Unknown family filters raise MaterialLookupError. None leaves membership
    unrestricted and retains the historical flat listing.
    """
    if family_id is not None:
        get_material_family(family_id)
    materials = tuple(
        _catalog_index().materials[key] for key in sorted(_catalog_index().materials)
    )
    if family_id is not None:
        materials = tuple(item for item in materials if item.family_id == family_id)
    if formula is None:
        return materials
    formula_key = formula.strip().casefold()
    return tuple(item for item in materials if item.formula.casefold() == formula_key)


def list_eos_records(
    *, formula: str | None = None, include_reference_variants: bool = False
) -> tuple[EOSRecord, ...]:
    """List source EOS choices; optionally include published scale variants.

    Reference variants remain executable by their original identifiers.
    """
    if not isinstance(include_reference_variants, bool):
        raise MaterialError("include_reference_variants must be a boolean")
    records = tuple(
        record
        for key in sorted(_catalog_index().records)
        if not (record := _catalog_index().records[key]).is_reference_variant
        or include_reference_variants
    )
    if formula is None:
        return records
    formula_key = formula.strip().casefold()
    return tuple(item for item in records if item.material.casefold() == formula_key)


def list_reference_variants(source_eos_record: str) -> tuple[EOSRecord, ...]:
    """List the original published pressure-scale alternatives of a source EOS."""
    source = get_eos_record(source_eos_record)
    return tuple(
        record
        for record in list_eos_records(include_reference_variants=True)
        if record.pressure_scale_variant_of == source.identifier
    )


# Re-export discovery here for users who naturally look for it beside lookup.
# The import is deliberately last: catalog_search imports the listing functions.
from peritheos.catalog_search import (  # noqa: E402, I001
    DeterminationMethod,
    RangeQuery,
    RangeSemantics,
    ValidationStatus,
    search_eos_records,
    search_materials,
)


__all__ = [
    "DeterminationMethod",
    "MaterialFamily",
    "MaterialGroup",
    "RangeQuery",
    "RangeSemantics",
    "ValidationStatus",
    "get_eos_record",
    "get_material",
    "get_material_family",
    "group_materials",
    "list_eos_records",
    "list_reference_variants",
    "list_material_families",
    "list_materials",
    "search_eos_records",
    "search_materials",
]
