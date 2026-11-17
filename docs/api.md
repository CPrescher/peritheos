# API reference

## Errors

All deliberate Python errors derive from `PeritheosError`; domain classes such
as `EosError` and `FitError` support precise handling, while compatibility with
the previous built-in `ValueError`, `ArithmeticError`, `TypeError`, and related
categories is retained. Instances expose `code`, `operation`, `field`, and
read-only `context` attributes. The native extension raises the same classes.

```python
from peritheos import EosError, EosValidationError, PeritheosError
```

See [Error handling](error-handling.md) for the complete hierarchy, stable
codes, Rust error kinds, and source-chain examples.

## Material families

`MaterialFamily(identifier, name, description, formula=None)` describes a
browsing family. `Material.family_id` is optional; `MaterialGroup(family,
materials)` carries a family and its matching materials, or one standalone
material with `family=None`. These types and functions are exported from both
`peritheos` and `peritheos.catalog`:

- `list_material_families()` returns definitions in identifier order.
- `get_material_family(identifier)` resolves an exact family ID.
- `list_materials(family_id=...)` filters membership, optionally with `formula`.
- `search_materials(family_id=...)` and `search_eos_records(family_id=...)`
  combine membership with the existing search criteria.
- `group_materials(materials)` groups any listing or search result without
  expanding membership or selecting a default material.

See [family browsing](catalog.md#browse-material-families) for consumer examples,
ordering, error behavior, and imported-material fallback.

## Materials and EOS records

```python
from peritheos import (
    EOSRecord,
    HugoniotBranchDomain,
    HugoniotInitialState,
    HugoniotRecord,
    HugoniotVolumeBasis,
    Material,
    get_eos_record,
    get_material,
    list_eos_records,
    list_materials,
    search_eos_records,
    search_materials,
)
from peritheos.materials import DEFERRED_EOS_RECORDS
```

`get_material(identifier)` returns a material phase and
`list_materials(formula=...)` lists or filters all bundled materials.
`get_material` resolves, in order: an exact identifier or alias (then a pre-0.7
compatibility identifier); a case-insensitive identifier or a case-sensitive
formula (`"Au"`, `"CO"` is not `"Co"`); a case-insensitive alias or name
(`"periclase"`); and a case-insensitive formula for an all-lowercase query
(`"au"`). The first step with exactly one match returns it. A step matching
several materials raises `MaterialLookupError` listing the candidates and their
phases (`context["matches"]`), so `get_material("Fe")` never picks a phase. An
unknown name raises with ranked suggestions.

`print(material.summary())` gives a plain-text overview: identity, cell contents
and volume unit, one `Default` line (the designated record, that none is
designated, or the stored hold reason), then one row per record, default first,
with model, the P/T range `check_validity=True` enforces, reproduction headline
and short reference; `*` marks a designated default and `!` a warning, and each
distinct warning is listed once. `record.summary()` shows one record's
coefficients with errors, fixed flags and units, the check with what it compared
and why, any warning on one line with its reason, the labelled ranges, reference and DOI, and the source
check with its note. Both are rendered from the [overview](#material-and-record-overviews)
dictionaries and read only the record file, never the audit sidecar. `repr()`
of a `Material` or `EOSRecord` is one identifying line;
`Material.default_record()` is a method.

`EOSRecord.reproduction` is `"reproduced"`, `"close"`, `"disagrees"` or
`"untested"`; `EOSRecord.compared` is `"published_fit"`, `"peritheos_refit"`,
`"source_equation"` or `"none"`; `EOSRecord.why` names the failed coefficients
of a `close` or `disagrees` record, else `None`; `EOSRecord.warning` is the
stored usage restriction or `None`, and `EOSRecord.warning_reason` why it
applies (only with a warning, else `None`); `EOSRecord.source_check` is `"verified"`,
`"pending"`, `"unresolved"`, `"deferred"` or `None`, and
`EOSRecord.source_note` what the check found for this record (the curated
per-record note, else the verbatim finding). [EOS labels](eos-labels.md) defines the
words. `EOSRecord.audit` loads the record's audit trail on first access.
`Material.get_eos_record(identifier)` matches exact identifiers or aliases only;
a miss suggests close identifiers without selecting one.
`get_eos_record(identifier)` explicitly retrieves any executable record.
`list_eos_records(formula=...)` lists source choices; pass
`include_reference_variants=True` to include published pressure-scale alternatives,
or use `list_reference_variants(source_identifier)` to inspect one source's
alternatives. These canonical objects are constructed from the bundled `.eosmat`
documents. The pre-0.7 convenience constants remain in an isolated
compatibility layer and are excluded from canonical listing and search. Each
`Material` owns its `eos_records`, supports `get_eos_record(identifier)`, and
provides `to_dict()`/`from_dict()` and `to_eosmat()`/`from_eosmat()` for the
canonical executable material document, preserving the loaded v3, v4 or v5 version on export. Optional crystallographic
fields are preserved but not interpreted. Loading uses a fixed model registry
and never imports an implementation path. The legacy executable snapshot-v2
reader and `to_snapshot_dict()` remain compatibility-only APIs.

### Material and record overviews

`EOSRecord.overview(*, include_audit=False)` and
`Material.overview(*, include_audit=False)` return plain JSON-serialisable
dictionaries with everything an application, notebook or command-line tool
needs to present a record, so consumers need not re-implement mappings or parse
summary text. Values are stored record fields or facts derived in code; nothing
is inferred, and a missing value is `None`. By default no audit sidecar is read
and `pressure_scale` is `None`; `include_audit=True` adds it, loading a bundled
material's sidecar on first use. `material_overviews()` yields the overview of
every bundled material without reading any sidecar, and
`write_material_overviews(path)` writes them as one JSON index
(`{"overview_version", "peritheos_version", "materials"}`).

```python
import json

from peritheos import get_eos_record, get_material

record = get_eos_record("mgo_speziale_2001_bm3_2")
overview = record.overview()
print(overview["check"]["reproduction"], overview["check"]["display"]["compared"])
for row in overview["parameters"]:
    print(row["display_name"], row["value"], row["error"], row["fixed"], row["unit"])
assert overview["ranges"]["effective"]["pressure_gpa"] == [0.0, 203.0]

material = get_material("palladium").overview()
print(material["default_record"], material["default_hold"])
json.dumps(material)  # plain data
```

Record overview keys, in order:

| Key | Content |
|---|---|
| `overview_version` | Layout version, currently `1`; keys are only added within a version |
| `identifier`, `label`, `aliases` | Stored record identity |
| `material` | `identifier`, `name`, `formula`, `phase` (identifier and name are `None` when unknown) |
| `citation`, `short_citation`, `title`, `doi` | One-line citation (first author et al., year, journal, volume, locator, DOI), `Speziale et al. 2001`, the stored title and DOI |
| `model`, `thermal_model` | `identifier` and display `name` from a fixed table (the identifier when no documented name exists); `thermal_model` is `None` for isothermal and Hugoniot records |
| `category`, `equation_kind` | `equilibrium` or `hugoniot` (the `default_records` key); `isothermal`, `thermal` or `hugoniot` |
| `record_kind`, `determination_method` | As stored, else `None` |
| `reference_temperature_k` | Reference temperature in K |
| `is_default` | The record `default_equilibrium_record()` or `default_hugoniot_record()` returns because it is designated |
| `selected_when_no_default` | Returned only because the material designates no default in its category (the first record) |
| `catalog_access`, `pressure_scale_variant_of` | As stored, else `None` |
| `check` | `reproduction`, `compared`, `why`, `warning`, `warning_reason`, `source_check`, `source_checked_on`, `source_note` as stored, plus `display`: a short `compared` phrase and a one-sentence `explanation` of the headline from fixed tables consistent with [EOS labels](eos-labels.md) |
| `parameters` | In model order: `component` (`model` or `thermal_model`), `name`, `display_name` (`K0'`), stored `value` and `error` (parameter units), `fixed` (`None` for records built in Python), `unit` (see below) |
| `error_confidence`, `covariance_parameters` | Stated error confidence, else `None`; the covariance parameter order, else `None` |
| `ranges` | `validity` and `experimental` (with its `status`) as stored, and `effective`: what `check_validity=True` enforces (`source` names which stored range that is; open bounds are `None`), plus a Hugoniot's always-enforced `particle_velocity_km_s` branch domain |
| `pressure_scale` | With `include_audit=True`: the audit's `pressure_calibration` `status` and `methods` (`kind`, `material`, `reference_eos_record`, `reference_calibration_record`, `citation`, `scope`); otherwise `None` |
| `notes`, `thermal_notes`, `validity_notes` | Stored text, unmodified |
| `audit` | The stored `check.audit` reference |
| `report_url` | The reconstruction-report anchor for a bundled record; `None` otherwise |

`unit` is the file's volume unit for `V0`, GPa for `K0`, `P0` and (inverse)
`K0''`, and otherwise the unit the documentation states for that model's
coefficient; `"1"` marks a documented dimensionless coefficient and `None` a
coefficient whose unit no documentation states. Material overviews contain
`overview_version`, `identifier`, `name`, `formula`, `phase`, `cell_contents`,
`volume_unit`, `formula_units_per_cell`, `aliases`, `family_id`, `notes`,
`default_record` (the identifier `default_record()` returns), `default_records`
(designated `equilibrium` and `hugoniot` defaults with `identifier` and `label`,
else `None`), `default_hold` (the stored hold reason, else `None`) and
`records`, the record overviews with the default first. The JSON schema is
packaged as `peritheos/data/overview-v1.schema.json`.

The installed command shows the same information:

```bash
python -m peritheos list periclase            # one line per matching material
python -m peritheos show mgo                  # Material.summary()
python -m peritheos show mgo mgo_speziale_2001_bm3_2 --json  # EOSRecord.overview()
```

`show` resolves the material like `get_material` and the record exactly, with
suggestions on a miss; a record identifier alone is also accepted.
`--include-audit` adds the pressure scale.

`HugoniotRecord` specializes `EOSRecord` with mandatory loading-path,
branch-kind, initial-state, operational volume-basis, and branch-domain fields.
It stays in the common
`Material.eos_records` collection; `Material.hugoniot_records` and
`Material.equilibrium_records` return convenient filtered views. Hugoniot
records expose `shock_velocity()`, `particle_velocity()`, `density()`, and
`specific_internal_energy_change()` in addition to `pressure()` and `volume()`.

For an equilibrium source record, `record.on_pressure_scale(target_node)` returns
a `PressureScaledEOS` view. Its `pressure()` and `volume()` methods transform
through the same audited calibration route in opposite directions, preserving
the source volume convention. `recalculate()` returns the full intermediate-state
provenance. See [direct pressure-scale views](pressure-scale-normalization.md#one-source-eos-selectable-output-scale).

`search_materials()` and `search_eos_records()` accept typed filters for free
text, name, alias, formula, phase, model family, DOI, author/reference,
thermal or caloric capability, calibration pressure/temperature, uncertainty,
and validation status. Results are executable objects ordered by stable
identifier. A two-value range uses `range_semantics="contains"` by default;
`"overlaps"` requests any closed-interval overlap. Records with missing range
metadata do not satisfy a range filter. See [Material
catalog](catalog.md) for signatures and examples.

## Supporting-study datasets

`list_studies() -> tuple[str, ...]` lists supporting studies with typed dataset
declarations. This discovery is separate from material/EOS lookup and excludes
prose-only research records. `get_study_document(identifier) -> dict` returns an
independent metadata copy, retaining phase uncertainty and scientific limitations.
Core installations include public study context and declarations; private source
acquisition history and numerical CSV values are outside that export.

`get_study_dataset(study, identifier) -> Dataset` loads the named numerical table
through the same checksum, read-only-array and optional-package/version checks
as material datasets. No EOS or phase is inferred from a supporting study.
Unknown studies raise `DatasetLookupError` with `dataset.study_not_found`;
unknown dataset identifiers use `dataset.not_found`.

```python
from peritheos import get_study_dataset

liquid = get_study_dataset("argon_barker_1987", "argon_barker_1987_table1_liquid")
assert len(liquid) == 7
assert liquid.metadata["phase"] == "liquid"
print(liquid["pv_over_nkt"])
```

This installed example requires `peritheos[datasets]`. See the
[Barker source/rights review](research/barker-1987-redistribution.md) for the
liquid comparison and selected graphical dataset qualifications.

## Observation datasets

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

`resolve_dataset_pressure(document, eos_record, dataset_identifier, *,
resource_root=None, check_validity=False) -> DatasetPressure` resolves an explicit
per-EOS pressure coordinate. It retains original pressures, row identities,
calibration provenance and validity flags, and delegates ruby conversion to the
existing calibration API. Missing reductions remain unresolved. See
[Dataset pressure reductions](dataset-pressure-reductions.md) for the result
contract, uncertainty/temperature qualifications, and a consumption example.

`Material.get_dataset(identifier) -> Dataset` loads an observation table from
embedded rows or a checksummed packaged CSV. `Material.datasets` continues to
expose raw metadata mappings for discovery and serialization. For external
resources, use `Dataset.from_mapping(metadata, *, resource_root=None)` with a
local resource directory; the default root is the package's data directory.

| Interface | Result or metadata |
|---|---|
| `len(dataset)` | Number of observation rows |
| `dataset[name]` | Read-only NumPy array indexed by schema column name |
| `dataset.get_column(name)` | `DatasetColumn` with `name`, `quantity`, `unit`, `role`, optional `of`, and extension `metadata` |
| `dataset.find_columns(*, quantity=None, role=None)` | Tuple of matching `DatasetColumn` objects |
| `dataset.values(name, *, unit=None)` | Read-only array, optionally converted to compatible units |
| `dataset.as_pressure_volume(*, pressure_unit=None, volume_unit=None, pressure_column=None, volume_column=None)` | `PressureVolumeData` view with pressure/volume arrays, units, associated uncertainties, and uncertainty roles |
| `dataset.resource` | `DatasetResource` with `path`, `sha256`, `media_type`, and extension `metadata`, or `None` for embedded rows |

`Dataset` also preserves `identifier`, `kind`, `reference`, `source_location`,
`description`, `notes`, `used_by_eos_records`, and extension `metadata`.
`checksum_verified` is `True` after resource verification and `None` for
embedded rows; `source_column_names` retains archived CSV headings.
`PressureVolumeData.dataset` links back to the loaded table. Its
`pressure_sigma` and `volume_sigma` aliases are available only for explicitly
declared standard deviations; other uncertainty roles retain their original
meaning. Missing numeric cells become `numpy.nan`. Unit conversion preserves
volume basis and does not infer material-specific basis changes.

Unknown dataset or column names raise `DatasetLookupError`; invalid loading,
ambiguous selection, checksum failures, and unsupported conversions raise
`DatasetError`. See [Loading observation datasets](datasets.md) for examples,
[the dataset schema](eosmat-schema.md#primary-experimental-datasets) for the
storage contract, and [Error handling](error-handling.md) for exception bases
and stable codes.

## Standalone `.eosmat` files

`create_eosmat(identifier, formula, model, parameters, *, temperature_ref, ...)`
returns a small, structurally valid native EOSMAT v5 dictionary for a user's own material.
It accepts an exact model type such as `BM3` or its canonical identifier
`birch_murnaghan_3`, and fills in the matching exchange identifiers. Parameters
use GPa, kelvin and conventional-cell Å³. The reference temperature is required.

```python
from peritheos import create_eosmat, save_eosmat

document = create_eosmat(
    "my_material",
    "X",
    "BM3",
    {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0},
    temperature_ref=300,
)
save_eosmat("my_material.eosmat", document)
```

These coefficients are fictional. Optional keyword arguments are `name`,
`reference`, `parameter_errors`, `fixed_parameters`, `formula_units_per_cell`
and `thermal`, plus `format_version=3` for older consumers. A thermal input mapping can provide its `type` or `model`, parameters
and complete formulation-specific options. No cell contents are inferred.
The creator constructs the selected model to check its configuration without
adding unprovided metadata or changing the source inputs.

New records remain diagnostic and unreviewed. V5 stores detailed checks in
`evidence` and keeps errors/fixed flags beside parameter values. See the
[native format reference](eosmat-v5.md).
`Material.from_eosmat(document, require_primary_validation=False)` explicitly
accepts them for calculations. Creation and numerical construction do not
establish source validation or fit reproduction. All full-schema fields can
subsequently be added to the returned dictionary. See the
[standalone authoring workflow](eosmat-authoring.md) for the installed
`python -m peritheos create` command, uncertainty handling and database contribution.

### Readable files and local metadata sharing

New standalone files display their units explicitly. V5 formatting places EOS
records before long crystallographic detail, and a record label before its
coefficients. Ordinary files stay inline; sharing is optional for repeated
metadata in large database files.

`convert_eosmat(document, version=5, compact=True)` deduplicates exact repeated
metadata into a file-local `shared_metadata` table. `inline=True` expands local
references for first-generation v4 readers. The two options are mutually
exclusive; v3 exports are always expanded. Unused definitions and all original
scientific fields are retained. Load/save preserves compact files directly;
model-based reserialization may normalize metadata just as for other EOSMAT files.

The installed commands `python -m peritheos validate FILE` and
`python -m peritheos format FILE [--check]` validate declarations or format layout
without scientific inference. See [the authoring guide](eosmat-authoring.md) for
complete examples, and [native v5](eosmat-v5.md) for the optional standardized
`evidence.reproduction.assessment` and its explicit comparison scope/basis.

## Shared `.eosmat` material library

```python
from peritheos import (
    eosmat_schema,
    get_material_document,
    list_material_documents,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
```

`list_material_documents()` returns the identifiers of all bundled
material cards, including source-only cards with no executable EOS.
`get_material_document(identifier)` returns a defensive copy in the historical
v3 dictionary layout so existing integrations keep working. The bundled files
themselves use native v5. Pass `format_version=5` to get native fields, or
`format_version=None` to retain the stored version. The same options apply to
`get_eos_record_document(identifier)`, including document-only records.
`load_eosmat()` always retains the file's authored layout and also accepts native
Dioptas 0.10.0 format-2 files;
`save_eosmat()` validates and preserves optional fields and the input version.
`load_eosmat()` reads v2–v5. `convert_eosmat(document, version=5)`
migrates a complete v3/v4 document. `version=3` or `version=4` exports
model-unit uncertainties without new flags for older readers. V5 errors and
covariance use parameter units; see [the uniform rule](eosmat-v5.md).
`eosmat_schema(5)` returns the native schema; the historical no-argument schema
call retains v3 for existing callers.
The [`.eosmat` schema reference](eosmat-schema.md) documents the complete
field contract and consumer defaults.

## Pressure-calibration records

```python
from peritheos import (
    find_common_pressure_calibration_routes,
    find_pressure_calibration_path,
    get_cross_calibration_edge,
    get_pressure_calibration,
    list_cross_calibration_edges,
    list_diamond_raman_calibrations,
    list_pressure_calibrations,
    list_ruby_pressure_calibrations,
    list_ruby_xrd_bridges,
    list_xrd_pressure_standards,
    recalculate_diamond_raman_pressure,
    recalculate_eos_pressure_scale,
    recalculate_pressure_calibration_path,
    recalculate_ruby_pressure,
    recalculate_ruby_to_xrd_pressure,
    recalculate_xrd_pressure_scale,
    xrd_standard_pressure,
)
```

Seven published ruby R1 calibrations and two diamond-anvil Raman-edge
calibrations are bundled with their equations, coefficients, provenance, and
validity metadata. Ruby supports wavelength, wavelength shift, and normalized
wavelength ratio; diamond Raman supports wavenumber and normalized wavenumber
ratio. Both provide analytical inversion and scale-to-scale conversion. The
explicit cross-calibration edge
registry records simultaneous experiments and jointly optimized families.
XRD-to-XRD conversion inverts the
source EOS to a virtual volume of the common standard and evaluates the target
EOS at that coordinate. Ruby-to-XRD conversion uses a ruby-linked standard EOS
as the first edge. The high-level sample-EOS API recursively discovers a path
from recorded provenance and returns every node and edge used. An edge may
also transform its temperature coordinate; the Chidester KCl-to-Dorogokupets
Pt edge converts the lower KCl effective temperature to the nominal average
Pt-foil temperature and records both edge-end values in `intermediate_states`. See
[Pressure-scale normalization](pressure-scale-normalization.md) for the valid
graph edges, complete examples, cross-calibration literature, and the
difference between a derived transformation and observation-level
re-reduction.

`find_common_pressure_calibration_routes(source_nodes)` discovers XRD-standard
endpoints reachable from every supplied sample EOS and returns the shortest
route for each source. Results are ranked by maximum and total edge count, not
by scientific preference; callers should restrict `target_nodes` to a suitable
standard or internally consistent family and verify the common calibration/data
coverage.

Transferred Dioptas records have completed a primary-source classification,
and native primary-sourced records include aragonite, KCl, RbCl, diamond, MgO,
CaSiO3, stishovite, akimotoite, Phase Egg, and Rh2O3(II)-type alumina.
The library contains 663 source records in 232 material documents. Of these,
645 are `primary_source_validated`; fourteen qualified, nondefault
`not_reproduced` implementations require explicit selection. The executable
catalog contains 659 records across 228 materials. Four `deferred` records
preserve the Campbell fcc-Fe/FeO, Morard FeS and held wadsleyite findings and are excluded from
executable discovery. This status
records source validation, not numerical refit parity; see
[Validation](validation.md) for reproduction results and scientific limitations.

`Material.from_eosmat()` constructs validated records and refuses deferred ones
by default; callers can inspect legacy values with
`require_primary_validation=False` and select records with
`record_identifiers=(...)`. This opt-in never changes the stored status.

The Benedict et al. diamond branch is available directly from the convenience
catalog or from the bundled `diamond.eosmat` document. Its public volume is the
eight-atom conventional-cell volume:

```python
from peritheos.materials import DIAMOND_BENEDICT_2014

volume = 8.0 * 4.6542704116  # A^3/conventional cell
pressure = DIAMOND_BENEDICT_2014.pressure(volume, 3000.0)  # about 150 GPa

# Forward phase-line state from a cold pressure, hot temperature, and
# experimentally assumed retained thermal-pressure fraction.
heated_volume = DIAMOND_BENEDICT_2014.volume_with_dac_confinement(
    60.0,
    2000.0,
    f_dac=0.25,
)
thermal_increment = DIAMOND_BENEDICT_2014.thermal_pressure_increment(
    heated_volume, 2000.0
)
confinement_increment = DIAMOND_BENEDICT_2014.dac_thermal_pressure(
    heated_volume, 2000.0, f_dac=0.25
)
total_hot_pressure = 60.0 + confinement_increment

# Both inputs remain public eight-atom conventional-cell volumes.
temperature = DIAMOND_BENEDICT_2014.temperature_from_volumes(
    ambient_volume=39.61987225,
    heated_volume=40.0,
    f_dac=0.25,
)
```

The record's `reference_volume` is the Table I **0 K motionless-ion cold-curve**
volume, not a 300 K zero-total-pressure state. The record also carries the
source's narrower DFT-MD comparison domain and diamond phase-stability caveat.

The earlier Correa logarithmic-moment branch is a separate convenience record:

```python
from peritheos.materials import DIAMOND_CORREA_2008

pressure = DIAMOND_CORREA_2008.pressure(8.0 * 4.43, 5000.0)
# 202.628115 GPa
```

It preserves the published DFT-GGA cold-curve volume and records the authors'
approximately 3% ambient-volume caveat rather than applying an implicit shift.

For example, the independently reproducible staged aragonite BM2 record can be
used at its 298 K reference state or at the represented high temperatures:

```python
from peritheos import Material, get_material_document

document = get_material_document("aragonite")
aragonite = Material.from_eosmat(
    document,
    record_identifiers=["aragonite_martinez_1996_bm2_2"],
).eos_records[0]

pressure = aragonite.pressure(volume=215.0)  # GPa at 298 K
hot_pressure = aragonite.pressure(volume=215.0, temperature=873.0)
result = aragonite.pressure_with_uncertainty(
    volume=215.0,
    volume_sigma=0.05,
)
```

The separate Martinez global thermal BM3 reduction is intentionally absent:
its published fitted reference volume is missing, and its remaining coefficients
do not reproduce the printed dataset under the documented equations.

`get_material_document()` remains the advanced raw-document API for exchange,
inspection, and editing. Normal calculations and discovery do not require a
manual `Material.from_eosmat()` step. Each `EOSRecord` provides:

- `pressure(volume, temperature=None, check_validity=False)` in GPa;
- `volume(pressure, temperature=None, check_validity=False)` in
  angstrom^3/conventional unit cell;
- `pressure_with_uncertainty(...)` and `volume_with_uncertainty(...)`;
- `within_calibration_range(volume, temperature=None)` (preferred) and the
  compatibility alias `within_validity(...)`;
- `reference`, `validity`, `parameter_provenance`, `notes`, and volume metadata.

See [Pressure standards](pressure-standards.md) for the subset and workflows
commonly used for pressure calibration.

## Isothermal equations of state

```python
from peritheos.eos.rt import (
    Baonza,
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Morse3,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    RydbergStacey,
    SunMorse3,
    SunMorse4,
    Vinet,
    Vinet3,
)
```

Common methods:

- `pressure(V)`
- `bulk_modulus(V)`
- `bulk_modulus_derivative(V)` for $dK/dP$; the base implementation is numerical
- `volume(P)` and `calculate_volume(P)`

Constructor signatures and special requirements are:

| Class | Signature after class name | Special requirement |
|---|---|---|
| `BM2` | `(V0, K0)` | $K_0'=4$ is implied |
| `BM2Reference` | `(V0, K0, P0)` | $P(V0)=P0$; `V0`, `K0` are reference-pressure values; $K'=4$ at that reference |
| `BM3` | `(V0, K0, K0_prime)` | none |
| `BM4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `Murnaghan` | `(V0, K0, K0_prime)` | supports the `K0_prime=0` limit |
| `NaturalStrain2` | `(V0, K0)` | $K_0'=2$ is implied |
| `NaturalStrain3` | `(V0, K0, K0_prime)` | none |
| `NaturalStrain4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `ModifiedTait` | `(V0, K0, K0_prime, K0_double_prime)` | rejects singular coefficient sets and volumes outside its real domain |
| `Baonza` | `(V0, K0, K0_prime)` | uses the published fixed pseudospinodal exponent $\beta=0.85$ |
| `Vinet` | `(V0, K0, K0_prime)` | none |
| `Vinet3` | `(V0, K0, eta, beta, psi)` | none |
| `RydbergStacey` | `(V0, K0, K0_prime, K_infinity_prime)` | `K_infinity_prime` sets the limiting high-pressure derivative |
| `Morse3` | `(V0, K0, K0_prime)` | three-dimensional Morse-potential form |
| `SunMorse3` | `(V0, K0, K0_prime)` | Sun--Morse form with $n=3$ |
| `SunMorse4` | `(V0, K0, K0_prime)` | Sun--Morse form with $n=4$ |
| `Holzapfel` | `(V0, K0, K0_prime, n, Z)` | molar volume in `J bar^-1 mol^-1` |

`Holzapfel` overrides `bulk_modulus_derivative(V, eps=1e-6)` with its existing
specialized numerical form.
See the [equation reference](equation-reference.md#isothermal-equations) for
the mathematical definitions and coefficient domains.

The isothermal API also includes `DensityPolynomial3`,
`SecondOrderMurnaghan` and `OddInversePower`. Their signatures are respectively
`(V0, rho0, c0, c1, c2, c3)`, `(V0, K0, K0_prime, K0_double_prime, P0=0)`
and `(V0, C3, C5, C7, C9)`. Their normalization volumes need not
be zero-pressure roots. The odd inverse-volume family has independent NumPy
and Rust implementations; its Python class uses the NumPy backend.

## Thermal equations of state

```python
from peritheos.eos.thermal import (
    AsymptoticDebyeTabulatedPressure,
    AsymptoticPowerLawMieGruneisenDebyeExcess,
    Dewaele2006,
    DorogokupetsOganov2007,
    DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz,
    HollandPowell2011,
    HollandPowellThermalPressure,
    LinearThermalPressure,
    QCompromiseDebye,
    LogVolumeThermalPressure,
    DebyeQuadraticThermalPressure,
    DebyeTabulatedThermalPressure,
    DebyeAnharmonicHelmholtz,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    VariableQMieGruneisenDebye,
    MultiOscillatorGruneisenThermalEOS,
    Sokolova2016,
    SecondOrderTaylorThermalPressure,
    SoundVelocityDebyeHelmholtz,
    Tange2009Debye,
    ThermalModifiedTait,
    ThermalReferenceStateEOS,
)
```

Thermal constructor signatures are:

| Class | Parameters after `rt_eos` |
|---|---|
| `AsymptoticDebyeTabulatedPressure` | `Tr, theta0, gamma0, a, b, n, cold_volume_ratio`; bounded electronic/residual tables and volume/temperature ranges |
| `AsymptoticPowerLawMieGruneisenDebyeExcess` | `Tr, theta0, gamma0, a, b, n, beta0, m` |
| `DoubleDebyeHelmholtz` | `Vp, theta_a0, a_a, b_a, theta_b0, a_b, b_b, theta_1_0, a_1, b_1`, followed by optional `n, alpha0, Ve, kappa, phi0, Tr=None` |
| `DoubleDebyeLogMomentHelmholtz` | `Vp, theta_a0, a_a, b_a, theta_b0, a_b, b_b, theta_0_0, a_0, b_0`, followed by optional `n, anharmonic_a, phi0, Tr=None` |
| `Dewaele2006` | `Tr, theta0, gamma0, gamma_inf, beta, anharmonic_a, anharmonic_m, electronic_e, electronic_g, n` |
| `DorogokupetsOganov2007` | `Tr`, four oscillator-mode parameter groups, `gamma0, gamma_inf, beta`, anharmonic, electronic, defect, and atom-count terms |
| `LinearThermalPressure` | `Tr, alpha_KT` |
| `QCompromiseDebye` | `Tr, theta0, gamma0, n`; constant θ and γ/V; no q or caloric potential |
| `SecondOrderTaylorThermalPressure` | `Tr, eta0, c0, c1, c2, c3, c4, c5` |
| `SoundVelocityDebyeHelmholtz` | `Tr, molar_mass_g_mol, n, longitudinal_intercept, longitudinal_slope, shear_intercept, shear_slope` |
| `LogVolumeThermalPressure` | `Tr, alpha_KT_ref, dK_dT_V` |
| `DebyeTabulatedThermalPressure` | `Tr, theta0, gamma0, q, n`; electronic temperature/pressure tables and explicit linear interpolation |
| `DebyeQuadraticThermalPressure` | `Tr, theta0, gamma0, q, n, A, m` |
| `DebyeAnharmonicHelmholtz` | `Tr, theta0, gamma0, q, b1, b2, b3` (monatomic; absolute-zero thermal baseline; gas constant 8.31451 J/mol/K) |
| `ThermalReferenceStateEOS` | `Tr, alpha0, dK_dT, alpha1=0, thermal_expansion_law="constant", reference_volume_law="integrated_expansivity", bulk_modulus_law="linear_temperature", beta1=0, beta2=0, beta3=0, kprime_log_coefficient=0, dKprime_dT=0`; volume laws also include `linear_temperature` and `berman` |
| `MieGruneisenDebye` | `Tr, theta0, gamma0, q, n, debye_temperature_law="integrated_gruneisen", thermal_pressure_reference="reference_temperature", Cvmax=None` |
| `VariableQMieGruneisenDebye` | `Tr, theta0, gamma0, q0, q1, n`; volume-dependent q and the source-defined gamma law |
| `MieGruneisenEinstein` | `Tr, theta0, gamma0, q, n` |
| `HollandPowellThermalPressure` | `Tr, theta, alpha0, n` |
| `ThermalModifiedTait` | `Tr, theta, alpha0, n` |
| `MultiOscillatorGruneisenThermalEOS` | `Tr, QE1o, mE1, QE2o, mE2, delta, t, a_0, m, g, e_0`, followed by optional `beta, QBo, d, mb, QB1o, d1, mb1, n` |
| `Tange2009Debye` | `Tr, theta0, gamma0, a, b, n` |

`ThermalReferenceStateEOS` also supports `bulk_modulus_law="reciprocal_cubic"`
with `dK_dT=0` and the `beta1`, `beta2`, and `beta3` coefficients. The optional
`dKprime_dT` (K^-1) adds a linear temperature slope to `K0_prime`;
`kprime_log_coefficient` adds the logarithmic term for supported families.
Both are preserved in fitting, uncertainty propagation and interchange.
See [Cubic reference compressibility](eosmat-schema.md#cubic-reference-compressibility)
for the equations, supported reference models, and domain restrictions.

`Sokolova2016` is the compatibility alias for
`MultiOscillatorGruneisenThermalEOS`.

`Tange2009Debye.constrained_reference_parameters(...)` implements the Fit 3
ambient constraints. Given trial `gamma0` and fixed `T0`, `K_S0`, `alpha0`,
`C_P0`, and `n`, it returns the dependent `K0`, `theta0`, and `Cv0` values used
by the source's four-variable global objective.

The Mie-Gruneisen, multi-oscillator, and constant linear thermal-pressure
classes accept any `EosBase` reference. `SoundVelocityDebyeHelmholtz` also
requires the reference to implement `bulk_modulus`. `LogVolumeThermalPressure` and
`SecondOrderTaylorThermalPressure` require `V0`;
`HollandPowellThermalPressure` requires `K0`; thermal modified Tait requires `ModifiedTait`; and
`ThermalReferenceStateEOS` requires a reference that reconstructs through
`V0` and `K0`. Energy-based thermal classes require molar volume in
`J bar^-1 mol^-1`; `LinearThermalPressure`,
`SecondOrderTaylorThermalPressure`, and `ThermalReferenceStateEOS` inherit any
volume unit consistent with their reference EOS. The second-order Taylor model
adds absolute thermal pressure to a cold curve; its pressure need not vanish at
`Tr`, although `thermal_pressure_increment()` is zero there by definition.
The double-Debye Helmholtz classes accept `Vinet`, `BM2`, `BM3`, or `BM4`
references with matching reference energies. With the
constructor default `Tr=None`, it represents the classical 0 K cold curve and
`thermal_pressure()` is the absolute non-cold contribution, including
zero-point pressure. Supplying `Tr` instead treats the reference curve as an
isotherm and subtracts the non-cold Helmholtz contribution at that temperature.
Each class exposes `cold_energy()`, `zero_point_energy()`,
`ion_helmholtz_free_energy()`, `anharmonic_helmholtz_free_energy()`,
`helmholtz_free_energy()`, `ion_pressure()`, and `anharmonic_pressure()`.
The caloric interfaces `ion_internal_energy(V, T)`,
`anharmonic_internal_energy(V, T)`, and `internal_energy(V, T)` return J/mol;
`ion_entropy(V, T)`, `anharmonic_entropy(V, T)`, and `entropy(V, T)` return
J/mol/K. The reference-isotherm subtraction is preserved in total internal
energy, while fixed-volume energy increments and entropy retain the thermal
model. The absolute energy zero remains conventional. See the
[diamond composite audit](literature-reproductions/diamond-dewaele-thermal-composites.md)
for the distinction between source models and derived anchored composites.
Ordinary `temperature(P,V)` inversion and DAC `temperature_from_volumes()`
inversion are supported. DAC increments subtract the pressure at the chosen
reference temperature (300 K for the default absolute models).
`ThermalModifiedTait` is the modified-Tait compatibility subclass, and
`HollandPowell2011` remains its alias. Exact equations and parameter roles are documented
under [Thermal equations](equation-reference.md#thermal-equations).
`Sokolova2016` is a compatibility alias for
`MultiOscillatorGruneisenThermalEOS`.

`parameter_values()` contains only numeric parameters that fitting and
uncertainty propagation may vary. `configuration_values()` reports fixed
non-numeric constructor choices such as `debye_temperature_law`;
`with_parameters()` preserves those choices when reconstructing an EOS.

Common methods:

- `thermal_pressure(V, T)`
- `thermal_pressure_increment(V, T)`
- `dac_thermal_pressure(V, T, f_dac)`
- `pressure(V, T)`
- `volume(P, T)`
- `volume_with_dac_confinement(P_cold, T, f_dac=...)`
- `temperature(P, V)` and `calculate_temperature(P, V)`
- `temperature_from_volumes(V_ambient, V_heated, f_dac=...)`
- `bulk_modulus(V, T)`
- `isothermal_compressibility(V, T)`
- `thermal_expansivity(V, T)`
- `molar_heat_capacity_v(V, T)` when a caloric model exists
- `molar_heat_capacity_p(V, T)` when a caloric model exists
- `adiabatic_bulk_modulus(V, T)` when a caloric model exists
- `gruneisen_parameter(V, T)` when a caloric model exists

`thermal_pressure_increment()` is the heating pressure above the reference
isotherm. It equals `thermal_pressure()` for reference-relative models and
`pressure(V, T) - pressure(V, Tr)` for either absolute double-Debye Helmholtz
class. `dac_thermal_pressure()` returns only `f_dac` times that increment.
`volume_with_dac_confinement()` solves the forward hot-volume problem from a
cold pressure and temperature; `temperature_from_volumes()` solves its
two-volume inverse. Both follow the empirical confinement model described in
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md);
they require `0 <= f_dac < 1`. In this API, `f_dac` means
`(P_hot - P_ambient) / Delta_P_thermal(V_heated, T)`, where the denominator is
the pressure increment above the reference-temperature isotherm. It is not a
fraction of the cold pressure.

Mie-Gruneisen models additionally expose `gruneisen_parameter()`,
`characteristic_temperature()`, and the vibrational thermodynamic methods
documented under [Thermoelastic properties](thermoelastic-properties.md).

## Shock Hugoniot equations of state

```python
from peritheos import HugoniotBase, HugoniotState, LinearUsUpHugoniot
```

`LinearUsUpHugoniot(V0, rho0, c0, s, P0=0)` implements a phase-specific
`Us = c0 + s * up` relation and the Rankine--Hugoniot jump conditions. See
[Shock Hugoniot equations of state](hugoniots.md) for equations, units,
fitting, EOSMAT records, and phase-branch semantics.

## Fitting

```python
from peritheos import EulerianFiniteStrainAcoustic
from peritheos.fitting import (
    AcousticFitResult,
    FitResult,
    HugoniotFitResult,
    fit_acoustic_finite_strain,
    fit_joint_eos,
    fit_linear_us_up,
    fit_rt_eos,
    fit_thermal_eos,
)
```

`EulerianFiniteStrainAcoustic` evaluates the corresponding adiabatic bulk and
shear moduli and Vp/Vs curves. `fit_acoustic_finite_strain()` fits those
third-order Eulerian density--Vp--Vs relations and returns an
`AcousticFitResult`. Density uncertainties are handled as shared latent
coordinates; 2-by-2 or 3-by-3 observation covariance matrices can preserve
within-row velocity or velocity/density correlations. See the [acoustic
finite-strain guide](fitting.md#acoustic-finite-strain-fitting).

All fitting APIs reject `bounds` entries that name unknown or fixed parameters.
Each supplied bound must apply to a free parameter; unused bounds are errors.

`FitResult` contains the fitted `model`, parameter and uncertainty mappings,
covariance and correlation matrices, raw and weighted residuals, chi-square,
degrees of freedom, AIC, BIC, convergence status, and solver message. It also
reports `adjusted_volume`, `adjusted_temperature`, `volume_corrections`, and
`temperature_corrections` for errors-in-variables fits. See
[Fitting P-V and P-V-T data](fitting.md) for the pressure, volume, and
temperature uncertainty API.

`FitResult.summary()` returns a compact text report. `FitResult.to_dict()` and
`FitResult.to_json()` export a versioned, JSON-safe record containing the model
identity and reconstructable parameters, fit arrays, diagnostics, and solver
configuration. Passing a path to `to_json()` also writes the record to disk.

Named-loss fits of exact `BM3`, and `QCompromiseDebye` over exact `BM3`,
automatically use analytic pressure Jacobians in `fit_rt_eos()`,
`fit_thermal_eos()` and `fit_joint_eos()` as applicable. The weighting,
latent-coordinate treatment, robust losses and parameter covariance/scaling
conventions remain the same. Other models, custom subclasses and callable
losses retain their finite-difference paths; the public signatures are unchanged.

`peritheos.diagnostics.birch_murnaghan_finite_strain_diagnostic()` constructs
the conventional Birch-Murnaghan normalized-stress diagnostic from P-V
observations. The returned `BirchMurnaghanFiniteStrainDiagnostic` exposes the
transformed arrays and an optional numerical model curve. Plotting is left to
the caller; the fitting guide links to a complete executable example.

`fit_joint_eos()` estimates reference-isotherm and thermal parameters in one
regression. Reference parameters use dotted names such as `rt_eos.V0`; its
covariance includes reference/thermal cross-correlations and is directly
compatible with `FitResult.eos_uncertainty()`.

`fit_linear_us_up()` returns a `HugoniotFitResult`. With no uncertainties it
performs OLS; shock-velocity uncertainties select WLS, and particle-velocity
uncertainties add latent-coordinate errors-in-variables fitting.
`HugoniotFitResult.to_dict()` and `.to_json()` serialize the fitted model,
covariance, adjusted particle velocities, diagnostics, and solver status.

## Units

```python
from peritheos.units import (
    cell_volume_to_molar_volume,
    convert_density,
    convert_molar_volume,
    convert_pressure,
    convert_temperature,
    density_from_molar_volume,
    molar_volume_to_cell_volume,
    molar_volume_from_density,
)
```

Every conversion accepts scalars or NumPy arrays. Cell-volume conversion uses
angstrom cubed per conventional cell and requires the crystallographic number
of formula units per cell; its molar result is per mole of formula units.
Pressure and temperature conversion imports from `peritheos.utils` are
deprecated compatibility aliases.

## Uncertainty propagation

```python
from peritheos import EOSUncertainty, ParameterUncertainty, PredictionUncertainty
```

`EOSUncertainty` wraps a deterministic EOS and provides `pressure()`, `volume()`,
`bulk_modulus()`, and the generic `evaluate()` method. Each returns a
`PredictionUncertainty` containing the nominal value, standard error, confidence
limits, statistical assumptions, and an optional output covariance matrix.

`FitResult.eos_uncertainty()` constructs the wrapper directly from fitted
parameter covariance. See [Uncertainty in EOS calculations](uncertainty.md).
Partial parameter-error sets are supported. Parameters omitted from the error
mapping or covariance ordering are treated as exact, rather than as having an
unknown error that Peritheos will estimate.
