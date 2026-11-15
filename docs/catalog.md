# Material catalog

Browse the [readable material summaries](material-summaries/index.md) for each
record's coefficients, units, reference state, source, limits, reproduction
headline and warning, and the [reconstruction reports](reconstruction-reports/index.md)
for how each record was checked. Both are generated from the bundled EOSMAT
files, and the report index counts the records by headline.

Peritheos exposes the curated library through the normal executable API.
`list_materials()` returns `Material` objects and `list_eos_records()` returns
`EOSRecord` choices; pass `include_reference_variants=True` to include published
pressure-scale reference variants. Both listings are ordered by stable
identifier and constructed from the same `.eosmat` files returned by the
advanced `get_material_document()` API. Deferred records are absent from
executable discovery; records available only by explicit selection carry that
restriction as their warning. See each material summary and the
[recommendation dispositions](research/eos-recommendation-dispositions-2026-10-17.md).

The [coverage overview](catalog-coverage.md) reports current material-level
static/thermal/Hugoniot coverage, calibration readiness and ledger differences.
The [material backlog](material-eos-candidates.md) separates open gaps,
completed additions and evidence-dependent holds.

## EOS record categories

Each record carries separate answers to three questions: **where the
coefficients come from** (`record_kind`), **whether the represented EOS was
checked against its source** (`scientific_validation.status`) and **whether a
refit reproduced it** (the `check.reproduction` headline and its detailed
code). They are independent fields, not a single quality ranking. All words
and codes are defined in [EOS labels](eos-labels.md). A new fit receives its
own record and never silently replaces the published coefficients. Older files
without `record_kind` mean `published`. Primary-source validation supports the
documented formulation and represented parameters; it does not establish the
accuracy of every prediction, reproduce the original regression, or validate
extrapolation. Author confirmation is evidence for a specified part of the
formulation; it is recorded with its scope and is not a separate global
validation status.

Other fields answer different questions: `equation_kind` distinguishes
isothermal, thermal and Hugoniot models; `determination_method` distinguishes
experimental, theoretical, hybrid and unknown origins of the underlying
determination. Pressure-calibration provenance and validity ranges have their
own assessments. None of these changes the classifications above.

### Fit reproduction

Besides the codes in [EOS labels](eos-labels.md#detailed-outcome-codes), the
[primary EOS refit ledger](primary-eos-refits.md) uses `refit_failed`:
extraction or numerical fitting failed before a usable comparison was produced.
This differs from a completed fit with discrepant coefficients. A paper-level
outcome can combine several record outcomes; see the
[paper investigation ledger](paper-investigation-ledger.md).

Individual source audits additionally use scoped `reproduction_status` values
such as `partially_reproduced` or `equation_reproduced`. These are additive audit
descriptions, rather than a closed schema enumeration. **Partially reproduced**
means that a specified component was verified or recovered while the complete
source procedure remains unreproduced; the accompanying note identifies which
component and which evidence remain missing. The
[Speziale (2001) MgO audit](literature-reproductions/speziale-2001-mgo.md) is a
worked example: a published, source-validated record with an author-confirmed
thermal formulation whose original fitting procedure is partially reproduced.

## Browse material families

Families group related composition variants for browsing. They contain no EOS,
default composition, or mixing model. Every individual material keeps its stable
identifier, exact formula, phase, datasets, and published records.

```python
from peritheos import (
    get_material_family,
    group_materials,
    list_material_families,
    list_materials,
    search_materials,
    search_eos_records,
)

families = list_material_families()
bridgmanite = get_material_family("bridgmanite")
members = list_materials(family_id="bridgmanite")
records = search_eos_records(family_id="bridgmanite", thermal=True)

# A consumer can use the same loop for families and standalone materials.
for group in group_materials(list_materials()):
    if group.family is None:
        material = group.materials[0]
        print(material.name, material.formula, material.phase)
    else:
        print(group.family.name, len(group.materials))
        for material in group.materials:
            print("  ", material.identifier, material.formula, material.phase)

# Group only the matching materials, without adding other family members.
results = search_materials(family_id="bridgmanite", formula="MgSiO3")
filtered_groups = group_materials(results)
```

`MaterialFamily` exposes `identifier`, `name`, `description`, and optional
`formula`. Its formula describes the family, not an exact sample.
`MaterialGroup` exposes `family` and a tuple of `materials`. An unassigned
material produces a standalone group with `family=None` and one material.
Family members still form a family group when only one matches a filter.
Consumers select a member before offering that material's EOS records.

The curated families are:

| Identifier | Scope |
|---|---|
| `bridgmanite` | Mg–Fe–Al orthorhombic silicate perovskites, including named samples, hydrous/vacancy-bearing variants and the hypothetical Al₂O₃ model endpoint. |
| `mg_fe_monoxide` | Solid Mg–Fe monoxides, including pure endmembers, substituted/deficient samples and their distinct structural phases. |
| `post_perovskite` | Mg–Fe–Al CaIrO₃-type Cmcm silicate post-perovskites, including the hypothetical Al₂O₃ model endpoint. |
| `phase_d` | Mg-rich, Al-bearing and Fe-bearing aluminous phase D, preserving each sample's reported composition and structural ordering. |
| `naalsio4_calcium_ferrite` | NaAlSiO₄ calcium-ferrite type and its measured Fe-free/Fe-bearing composition variants. |
| `magnesite` | MgCO₃ and natural Fe/Mn-substituted magnesite. |
| `olivine` | Forsterite and fayalite; wadsleyite and ringwoodite remain separate. |
| `garnet` | Pyrope, almandine and majorite, retaining cubic/tetragonal structural distinctions. |

CaSiO₃ perovskite entries remain standalone in this composition-focused curation.
The family mechanism can group phase variants, but does not automatically group
materials merely because they share a formula or structure type.

Membership follows the identity documented on each card. A shared structure
name alone is insufficient: CaSiO₃, CaIrO₃ and germanate perovskites are not
assigned to the bridgmanite family. Cmmm post-perovskite-II is outside the
Cmcm family. MgO liquid and higher iron oxides remain outside the solid
monoxide family. Hypothetical model endpoints retain their existing labels;
family membership does not assert stability or continuous miscibility.

Families have a separate lookup namespace: `get_material("bridgmanite")` still
returns the original MgSiO₃ material. Existing flat listings and searches retain
their return types and ordering. Family filters use exact IDs, combine with
existing filters, and raise `MaterialLookupError` for unknown IDs; `None` means
no family restriction. Free-text search retains its existing matching rules.

`group_materials()` accepts any iterable of materials and preserves the supplied
objects. Members are ordered by material ID; groups by family ID or standalone
material ID, with the family first on a collision. Duplicate material IDs raise
`MaterialError`. Imported files can carry an unknown `family_id`: they remain
loadable, preserve the ID on export, and appear as standalone entries until a
corresponding definition is available. Grouping an empty result returns `()`.

The bundled registry is `peritheos/data/material-families.json`. Membership is
stored only in each material's optional `family_id`; member lists are derived.
Bundled references must resolve, including on structure-only cards. The
executable catalog and grouping APIs continue to omit structure-only cards;
these remain accessible through the document API.

## Look up and execute a record

```python
from peritheos import get_eos_record, get_material

gold = get_material("gold")
scale = get_eos_record("gold_fei_2007_vinet_2")

pressure = scale.pressure(volume=55.0, temperature=2000.0)
assert scale in gold.eos_records
```

An unknown identifier raises `MaterialLookupError` and includes close matches
when useful. Material names and aliases such as `Periclase` are discovery
metadata; stable identifiers remain the safest choice for saved analyses.

## Search executable materials

```python
from peritheos import search_materials

carbonates = search_materials(text="carbonate")
post_aragonite = search_materials(
    formula="CaCO3",
    phase="Pmmn",
    model_family="BM3",
)
thermal_diamond = search_materials(
    formula="C",
    model_family="double Debye",
    thermal=True,
    caloric=True,
    pressure_gpa=500.0,
    temperature_k=5000.0,
)
```

Material identity filters apply to the material. Record-level filters are
combined on one record: a material is returned only if at least one of its
records satisfies every requested record criterion.

## Search executable records

```python
from peritheos import search_eos_records

fei_gold = search_eos_records(
    formula="Au",
    author="Fei",
    doi="10.1073/pnas.0609013104",
    thermal=True,
    uncertainty=True,
)

validated_vinet = search_eos_records(
    model_family="Vinet",
    validation_status="primary_source_validated",
)
```

Use `determination_method="experimental"`, `"theoretical"`, `"hybrid"`, or
`"unknown"` to select how the EOS parameters were established. For example,
`search_eos_records(formula="SrSiO3", determination_method="theoretical")`
selects the calculated Xiao records while excluding the experimental fit in
the same paper. Read the value from `record.determination_method`; the detailed
rationale is in `record.eosmat_metadata["parameter_provenance"]["determination_method"]`
for bundled records. This is independent of whether the record is published,
refitted, or source-validated. See the [schema definitions](eosmat-schema.md#eos-record).

Both search functions support these keyword filters:

| Filter | Meaning |
|---|---|
| `text`, `name`, `alias` | Case-insensitive word search over the corresponding identity fields. Free text also covers stable identifiers and references. |
| `formula` | Case-insensitive exact chemical-formula match. |
| `phase` | Phase, material name, symmetry, or phase alias text. |
| `model_family` | Public model discriminator or class name, including either component of a thermal composition. |
| `determination_method` | `experimental`, `theoretical`, `hybrid`, or `unknown`; `None` leaves parameter origin unrestricted. |
| `doi` | Exact DOI after removing an optional `doi:` or `https://doi.org/` prefix. |
| `author`, `reference` | Author text, or broader publication and primary-audit metadata. |
| `thermal`, `caloric` | Require or exclude temperature-dependent pressure and caloric capability. |
| `uncertainty` | Require or exclude published parameter error/covariance data. Measurement-only uncertainty remains available independently. |
| `validation_status` | One status or an iterable of statuses. |

Every returned item is directly executable; search does not return a separate
summary or result-wrapper type.

## Calibration-range semantics

`pressure_gpa` and `temperature_k` accept either a scalar or an ordered
two-value tuple. Bounds are closed. A scalar selects records containing that
point. For a range, the default `range_semantics="contains"` requires the
record's published calibration interval to contain the whole requested range:

```python
complete_coverage = search_eos_records(
    pressure_gpa=(50.0, 100.0),
    temperature_k=(300.0, 2000.0),
)
```

Use `range_semantics="overlaps"` when any closed-interval overlap is enough:

```python
partial_coverage = search_eos_records(
    pressure_gpa=(100.0, 150.0),
    range_semantics="overlaps",
)
```

For an isothermal record without a separately repeated temperature interval,
the searchable calibration temperature is its declared reference isotherm. A
missing pressure range, or a missing temperature range on a thermal record,
does not match a range query: unknown coverage is never treated as unbounded
coverage. Search ranges describe calibration/data coverage, never hard
evaluation limits.

## Raw documents and compatibility names

Use `get_material_document()` when editing, exchanging, or inspecting raw
format-3 data. Normal lookup and search already perform validated construction.

The pre-0.7 pressure-scale constants in `peritheos.materials` and their old
generated identifiers remain executable. They are an isolated compatibility
layer and do not add records to canonical listing or search results. This
distinction is necessary where primary-source audit corrections made a
historical convenience parameterization numerically different from its
canonical `.eosmat` counterpart.

Canonical identifiers win the single material-identifier collision:
`get_material("diamond")` now returns the seven-record document-built material
rather than the former five-record convenience grouping. Every historical
diamond record remains available through its module-level constant and record
identifier. The two colliding anchored record identifiers have equation-identical
canonical counterparts, so their numerical behavior is unchanged.

NaCl materials and EOS labels explicitly distinguish NaCl-B1 (rock-salt,
Fm-3m) from NaCl-B2 (CsCl-type, Pm-3m). The phase description is also
available through `Material.phase` and `EOSRecord.phase`.

The NaCl-B2 catalog retains eight [Sakai et al. (2011) published fits](literature-reproductions/sakai-2011-nacl-b2.md). Normal discovery shows the Matsui-Pt BM3 and Vinet sources; the other six fits and the later 2014 Yokoo-Pt recalibration are published pressure-scale reference variants. Use `source.on_pressure_scale(target_node)` to evaluate a source on another reachable scale. `list_reference_variants(source.identifier)` and `include_reference_variants=True` expose the original alternatives. See the [catalog-wide variant audit](research/pressure-scale-variant-audit.md) for other grouped materials and the [conversion interface](pressure-scale-normalization.md#one-source-eos-selectable-output-scale) for pressure and volume calculations.

Platinum and gold also have executable pressure-scale branches and explicit
diagnostic PVT reconstructions; see the
[Matsui (2009) Pt](literature-reproductions/matsui-2009-platinum.md),
[Yokoo (2009) Pt](literature-reproductions/yokoo-2009-platinum.md) and
[Yokoo PVT library](literature-reproductions/yokoo-2009-pvt-library.md) audits.

## Research collections

The [Dewaele (2008) PAW collection](research/dewaele-2008-paw.md) preserves the
paper's published computational parameterizations, provenance and graphical
checks outside the main database. Its records are absent from catalog listings,
searches and default selection; the paper's experimental Table IV EOS remain in
the catalog.
