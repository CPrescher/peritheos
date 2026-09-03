# Material catalog

Peritheos exposes the complete bundled library through the normal executable
API. `list_materials()` returns 115 `Material` objects and
`list_eos_records()` returns their 150 `EOSRecord` objects. Both are ordered by
stable identifier and constructed from the same `.eosmat` files returned by
the advanced `get_material_document()` API.

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

Both search functions support these keyword filters:

| Filter | Meaning |
|---|---|
| `text`, `name`, `alias` | Case-insensitive word search over the corresponding identity fields. Free text also covers stable identifiers and references. |
| `formula` | Case-insensitive exact chemical-formula match. |
| `phase` | Phase, material name, symmetry, or phase alias text. |
| `model_family` | Public model discriminator or class name, including either component of a thermal composition. |
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
