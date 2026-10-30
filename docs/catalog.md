# Material catalog

Browse the [readable material summaries](material-summaries/index.md) for
parameter comparisons, uncertainties, calibration, selection holds and source
qualifications. These sheets are generated from the complete EOSMAT documents.

Peritheos exposes the main curated library through the normal executable
API. `list_materials()` returns 230 `Material` objects and
`list_eos_records()` returns 624 `EOSRecord` choices. Pass
`include_reference_variants=True` to include all 676 executable records. Both
material and EOS listings are ordered by
stable identifier and constructed from the same `.eosmat` files returned by
the advanced `get_material_document()` API.

Thorianite and cotunnite-type ThO2 have separate experimental thermal records
from Chidester et al. (2018); see the [source audit](literature-reproductions/chidester-2018-tho2.md)
for their phase-dependent data selections and reference-state limits.

The [coverage overview](catalog-coverage.md) reports current material-level
static/thermal/Hugoniot coverage, calibration readiness and ledger differences.
The [material backlog](material-eos-candidates.md) separates open gaps,
completed additions and evidence-dependent holds.

Fe-bearing phase E (`phase_e_fe0072`) preserves the Crichton–Ross (2000)
measured composition and published static BM3. Its independent refit achieves
parity within reported uncertainty; the [audit and residual plot](literature-reproductions/phase-e-crichton-ross-2000.md)
show the small coefficient and pressure-curve differences. Its phase-specific
diffraction fallback supplies indexed peak positions with uniform display
markers; it does not assert sample-specific atomic occupancies or intensities.

Brucite (`brucite`) retains the individually source-checked Fukui (2003)
300/473/673 K BM3 isotherms and a primary ambient neutron structure.
**The full published thermal EOS is not recommended: parity with the published
P–V–T data was not achieved.** The unresolved 873 K curve, continuous source
reconstruction and all Peritheos refits are withheld from executable selection
and preserved as audit evidence. These three isotherms do not provide a
continuous temperature model. The [audit](literature-reproductions/brucite.md)
documents the normalization discrepancy and diagnostic fits.

Seven additional published EOSs are explicitly not recommended: Fu's cubic
CaSiO3 thermal fit, Baty's experimental Pd BM3 and Sakai's five Fe-Ni thermal
variants. They remain available only for diagnostic selection, with
`not_reproduced` status. Wadsleyite's sole EOS is deferred and absent from
executable discovery; no replacement default is selected. See the
[recommendation dispositions](research/eos-recommendation-dispositions-2026-10-17.md).

Shen and Smith's MgO and NaCl-B2 records are also `not_reproduced` and not
recommended for quantitative use, with explicit diagnostic selection only.
Their published curves disagree with direct fixed-volume fits to the available
Table S1 observations. Exact author fit inputs remain unconfirmed; see the
[source-input audit](literature-reproductions/shen-smith-2026-cu-refits.md).

## EOS record categories

Peritheos records three separate assessments: **where the coefficients come
from**, **whether the represented EOS has been checked**, and **how much of
the fitting procedure has been reproduced**. Read them together with the
record's scope and qualifications. They are independent fields, not a single
quality ranking.

### Record origin

The `record_kind` field identifies the origin of the represented coefficients:

| Category | Meaning |
|---|---|
| Published (`published`) | The source publication's parameterization, including its fixed or adopted coefficients. |
| Refit (`refit`) | A separately attributed fit using a documented dataset, selection, equation, objective and weights. Its coefficients belong to that refit. |
| Derived (`derived`) | A model constructed through a documented transformation or composition of source results. |
| Diagnostic (`diagnostic`) | An investigation or comparison model whose qualifications limit its scientific use. |

A new fit receives its own record and provenance. It does not silently replace
the published coefficients. Older files without `record_kind` mean `published`.

### Scientific validation

The `scientific_validation.status` field describes the audit of the represented
EOS. Structural validity only establishes that the file has an accepted format.

| Category | Meaning |
|---|---|
| Primary-source validated (`primary_source_validated`) | Equation, coefficients, units, normalization and reference state have been checked against primary evidence. For a refit or derived record, source inputs and model choices are traced, while new coefficients or transformations are identified as Peritheos results. |
| Pending (`pending_primary_source_check`) | The primary-source audit is incomplete. |
| Deferred (`deferred`) | Required evidence is missing, ambiguous or inconsistent enough to withhold normal executable use. |
| Not reproduced (`not_reproduced`) | The record carries a documented unresolved reproduction problem and requires the stated qualifications and explicit acceptance. Some such records remain executable for investigation. |

Primary-source validation supports the documented formulation and represented
parameters. It does not establish the accuracy of every prediction, reproduce
the original regression, or validate extrapolation. Author confirmation is
evidence for a specified part of the formulation; it is recorded with its scope
and is not a separate global validation status.

### Fit reproduction

The [primary EOS refit ledger](primary-eos-refits.md) uses these outcomes:

| Category | Meaning |
|---|---|
| Parity (`parity`) | Recovered free coefficients meet the documented uncertainty and numerical comparison criteria. Conditional or source-specific criteria are stated explicitly. |
| Similar (`similar`) | The independent result meets documented numerical similarity criteria or a qualified uncertainty comparison. Similarity is not a statistical confidence claim. |
| Parity not achieved (`parity_not_achieved`) | A fit produced a comparison, but coefficient agreement or a dedicated published-curve reproduction check remains unresolved. |
| Source reconstruction (`source_reconstruction`) | Source equations, calculations or compositions are checked without independently refitting the complete EOS. |
| Direct refit unavailable (`not_refittable`) | Available observations or reduction details do not support an independent fit of the complete EOS. |
| Refit failed (`refit_failed`) | Extraction or numerical fitting failed before a usable comparison was produced. This differs from a completed fit with discrepant coefficients. |

These outcomes compare coefficients under the documented reconstruction. Even
parity does not by itself recover the author's row selection, objective,
weights, covariance, error estimates or calculator. For a Peritheos refit,
reproduction refers to that independent fit unless the original publication is
explicitly the comparison target. A paper-level outcome can also combine several
record outcomes; see the [paper investigation ledger](paper-investigation-ledger.md).

Individual source audits additionally use scoped `reproduction_status` values
such as `partially_reproduced` or `equation_reproduced`. These are additive audit
descriptions, rather than a closed schema enumeration. **Partially reproduced**
means that a specified component was verified or recovered while the complete
source procedure remains unreproduced. The accompanying note identifies which
component and which evidence remain missing.

### Example: Speziale (2001) MgO

The record `mgo_speziale_2001_bm3_2` is **published**, **primary-source validated**,
and has an **author-confirmed thermal formulation**. It uses the paper's final
BM3 plus variable-q Debye coefficients. Independent calculations verify the
pressure evaluation, and Sergio Speziale confirms integration of the defining
Debye-temperature relation using Equation (11).

Its complete original fitting procedure remains **partially reproduced**. The
cold coefficients are recovered within published error widths, while the full
thermal regression and original calculator settings are unavailable. A new
refit is not required to use this published EOS. A future independent refit
would be a separate record. See the [Speziale audit](literature-reproductions/speziale-2001-mgo.md).

Other fields answer different questions: `equation_kind` distinguishes
isothermal, thermal and Hugoniot models; `determination_method` distinguishes
experimental, theoretical, hybrid and unknown origins of the underlying
determination. Pressure-calibration provenance and validity ranges likewise
have their own assessments. None of these automatically changes the three
classifications above.

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

The [Matsui (2009) Pt audit](literature-reproductions/matsui-2009-platinum.md)
registers `platinum_matsui_2009_vinet_mgd_electronic` alongside the unchanged
300 K Vinet branch. It includes integrated-q Debye pressure and all 51
Tsuchiya-Kawamura electronic nodes, subtracting their 300 K value once.
Linear interpolation between 100 K nodes is explicitly a numerical choice;
no table extrapolation or electronic caloric properties are supplied.
The Au-only erratum does not change Pt. Source-grid reconstruction agrees
within 0.0071 GPa; the original shock/expansion fit remains unreproduced.

The [Yokoo (2009) Pt audit](literature-reproductions/yokoo-2009-platinum.md)
adds the published fitted 300 K Vinet branch, including its density-derived
V0=60.55 Å³ per four-atom cell. It resolves Sakai (2018)'s named calibration
reference while retaining the missing paired-observation boundary. The
derived Table V grid and original joint-fit status are documented separately.
Both `gold_yokoo_2009_pvt_reconstruction` and
`platinum_yokoo_2009_pvt_reconstruction` now provide bounded executable
[PVT reconstructions](literature-reproductions/yokoo-2009-pvt-library.md).
These explicit-selection diagnostic records have validation status
`not_reproduced`. They retain explicit interpolation and, for Pt, a separate
empirical residual-pressure correction; exact author-fit reproduction is
not claimed. The [pressure-convention audit](literature-reproductions/yokoo-2009-pressure-conventions.md)
tests the unreproduced analytical PVT independently of the fitted corrections.

## Research collections

The [Dewaele (2008) PAW collection](research/dewaele-2008-paw.md) preserves all
54 published computational parameterizations, provenance and graphical checks
outside the main database. Its records are absent from catalog listings, searches
and default selection. The paper's 13 experimental Table IV EOS remain in the catalog.
