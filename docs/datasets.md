# Loading observation datasets

Material cards can include the observations used to fit or validate their EOS
records. Peritheos exposes those tables as typed, read-only NumPy columns. The
same API works for small tables embedded directly in an `.eosmat` document and
for CSV resources distributed with the package; pandas is not required.
See the [API reference](api.md#observation-datasets) for signatures, the
[EOSMAT dataset schema](eosmat-schema.md#primary-experimental-datasets) for the
storage contract, and [Error handling](error-handling.md) for exception types
and stable codes.

## Find and load a dataset

`Material.datasets` remains the serialized metadata representation. It is
useful for discovering stable identifiers without reading the table itself:

```python
from peritheos import get_material

material = get_material("coesite")
identifiers = tuple(item["identifier"] for item in material.datasets)
print(identifiers)
# ('coesite_levien_1981_table7_pv',)
```

Pass one of those identifiers to `get_dataset()` to load the typed table:

```python
dataset = material.get_dataset("coesite_levien_1981_table7_pv")

print(dataset.identifier)
print(dataset.kind)
print(dataset.source_location)
print(dataset.reference)
print(dataset.used_by_eos_records)
print(len(dataset))
```

The returned `Dataset` preserves its description, reference, source location,
notes, EOS-record links, resource metadata, and unrecognized extension
metadata. An unknown identifier raises `DatasetLookupError` and reports the
available identifiers, with close-match suggestions when possible.

### FeS reference data and author input variants

The `fes_vi` source document preserves Table S1 and its paired KCl markers,
plus these separately accessible source datasets:

| Dataset identifier | Rows | Source version |
|---|---:|---|
| `fes_morard_2026_author_without_additional_cold_pvt` | 167 | Author input excluding the 11 additional cold observations |
| `fes_morard_2026_author_literature_cold_pv` | 21 | Sata/Ohfuji reference observations as supplied in the author input |
| `fes_sata_2010_table1_vi_pv` | 13 | Sata (2010) numerical re-report with original table errors |

The EOS is deferred, but its observation tables can be loaded directly without
enabling the source equation:

```python
from peritheos import Dataset, get_material_document

document = get_material_document("fes_vi")
identifier = "fes_morard_2026_author_literature_cold_pv"
metadata = next(d for d in document["datasets"] if d["identifier"] == identifier)
dataset = Dataset.from_mapping(metadata)
pv = dataset.as_pressure_volume()
print(len(dataset), pv.pressure_unit, pv.volume_unit)  # 21 GPa cm^3/mol
```

The author datasets preserve cm³/mol, assigned fit errors and original input
row indices. Their shared `observation_id` identifies overlap between the input
and subsets; concatenating these tables would count measurements repeatedly.
Their thermal rows also overlap Table S1. Original Sata cell volumes and error
tokens remain separate, including 101.1 GPa where the author input uses 101.2.
Thirteen cold cell volumes match the Sata re-report; individual original-study
attribution for the remaining eight is unresolved. None of these uncertainty
columns is silently declared a standard deviation. The 167-row selection is confirmed
by user-reported personal communication from Guillaume Morard: the 11 unsuitable
quenched observations are excluded. Neither that subset nor the 178-row variant
containing it is exposed in the database. Original files and historical comparisons
remain in the audit archive outside the package. See the
[author-source audit](literature-reproductions/fes-eosfit-direct.md#author-input-follow-up-2026-10-03)
for EOS-file versus paper coefficients and replay results.

## Pressure-volume convenience view

`as_pressure_volume()` locates pressure, volume, and associated uncertainty
columns from their schema `quantity`, `role`, and `of` metadata. Units may be
requested as part of the conversion:

```python
pv = dataset.as_pressure_volume(pressure_unit="GPa")

print(pv.pressure)
print(pv.volume)
print(pv.pressure_unit)  # GPa
print(pv.volume_unit)  # angstrom^3
```

Levien and Prewitt report pressure in kbar, so this call converts it to GPa.
The conventional-cell volumes remain in the source's cubic-angstrom unit.
Missing numerical cells are represented by `numpy.nan`.

The coesite table declares both pressure and volume errors as standard
deviations, making the `*_sigma` aliases available:

```python
print(pv.pressure_sigma)
print(pv.volume_sigma)
```

The first pressure sigma is `nan`, because the ambient row does not report a
pressure error. The complete uncertainty arrays are also exposed as
`pressure_uncertainty` and `volume_uncertainty`. Their meanings are recorded in
`pressure_uncertainty_role` and `volume_uncertainty_role`, which can be
`standard_deviation`, `standard_error`, or `uncertainty`. A `*_sigma` property
is `None` unless the source column is explicitly a standard deviation; this
prevents a confidence interval or standard error from being silently treated
as one sigma. Bound columns remain available through the general column API and
are not automatically converted into symmetric uncertainties.

Some datasets contain multiple pressures or volumes, such as separate sample
and calibrant series. Peritheos does not guess between equally valid columns.
Select them explicitly:

```python
platinum = get_material("platinum")
two_scales = platinum.get_dataset("platinum_dewaele_2004_table1_compression")
pv = two_scales.as_pressure_volume(
    pressure_column="ruby_pressure_revised_gpa",
    volume_column="atomic_volume_a3",
    pressure_unit="GPa",
)
```

The column arguments use schema names, not necessarily the headings printed in
the archived CSV.

## General column access

Every schema column has typed metadata:

```python
for column in dataset.columns:
    print(column.name, column.quantity, column.unit, column.role, column.of)
```

Indexing a dataset by schema name returns a read-only NumPy array:

```python
pressure_kbar = dataset["pressure_kbar"]
volume = dataset["volume_a3"]
```

Use `get_column()` for one column's metadata and `find_columns()` for semantic
selection. This is the main interface for P-V-T, lattice-parameter, shock,
elasticity, density, energy-volume, and other dataset kinds:

```python
temperature_columns = dataset.find_columns(quantity="temperature")
value_columns = dataset.find_columns(role="value")

column = dataset.get_column("pressure_kbar")
pressure_gpa = dataset.values(column.name, unit="GPa")
```

Arrays returned directly or after conversion are read-only. Convert or copy
them explicitly before mutation:

```python
editable_pressure = dataset.values("pressure_kbar", unit="GPa").copy()
```

## Unit conversion

`values(name, unit=...)` and `as_pressure_volume()` perform explicit compatible
conversions. Supported families include:

| Family | Units |
| --- | --- |
| Pressure | `Pa`, `kPa`, `MPa`, `GPa`, `bar`, `kbar`, `Mbar` |
| Length | `m`, `mm`, `nm`, `pm`, `angstrom` |
| Volume | `m^3`, `cm^3`, `nm^3`, `angstrom^3` |
| Temperature | `K`, `degC` |
| Velocity | `m/s`, `km/s` |
| Density | `kg/m^3`, `Mg/m^3`, `g/cm^3` |
| Dimensionless | `1`, `dimensionless` |

Volume basis suffixes are significant. For example,
`angstrom^3/formula_unit` can be converted to another cubic-length unit with
the same `/formula_unit` basis, but not directly to `cm^3/mol`. Such a change
requires material-specific amount-of-substance information and is rejected
rather than inferred. Unsupported or dimensionally incompatible conversions
raise `DatasetError`.

Temperature uncertainty conversion applies only the scale, not the Celsius
offset. Thus a `2 degC` standard deviation is also `2 K`, while a value of
`26.85 degC` becomes `300 K`.

## Packaged resources and integrity

For a resource-backed dataset, loading follows this order:

1. Resolve the package-relative local path.
2. Read the resource bytes.
3. Calculate and compare the SHA-256 checksum.
4. Parse the CSV only after the checksum matches.
5. Map source columns to schema columns by their declared order.

The resulting object exposes `resource.path`, `resource.sha256`, and
`resource.media_type`. `checksum_verified` is `True` for a successfully loaded
resource and `None` for embedded rows. A checksum mismatch raises
`DatasetError` with code `dataset.checksum_mismatch`.

Archived CSV headings sometimes reproduce a paper or deposit verbatim, while
the `.eosmat` column definitions provide stable normalized names. Original CSV
headings are retained in `source_column_names`; indexing always uses the schema
names in `Dataset.columns`.

`Dataset.from_mapping()` also supports compatible external metadata. Supply a
directory containing its resource paths:

```python
from pathlib import Path

from peritheos import Dataset

dataset = Dataset.from_mapping(
    dataset_metadata,
    resource_root=Path("path/to/material-data"),
)
```

Resource paths must be relative and local. Only `text/csv` resources are
currently supported.

## Backward compatibility

The loader is additive. Existing code that reads `Material.datasets` as raw
mappings or serializes a material with `to_eosmat()` continues to use the
original representation. `get_dataset()` constructs a typed view without
rewriting the material card or resource.

The public classes and errors can be imported from the package root:

```python
from peritheos import (
    Dataset,
    DatasetColumn,
    DatasetError,
    DatasetLookupError,
    DatasetResource,
    PressureVolumeData,
)
```

## Resolve the pressure coordinate used by an EOS

For datasets with alternative pressure scales, use
[`resolve_dataset_pressure`](dataset-pressure-reductions.md) to select the
explicit per-EOS coordinate. `as_pressure_volume()` remains a view of reported
columns and never infers an EOS reduction or transforms its pressure scale.
