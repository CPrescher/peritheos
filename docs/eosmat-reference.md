# EOSMAT v5 record reference

This page defines every field of a current (v5) `.eosmat` file. The
[authoring guide](eosmat-authoring.md) shows how to write one, and
[EOS labels](eos-labels.md) defines the words used in `check` blocks. Older
files are described under [format versions](eosmat-v5.md#old-files-and-compatibility-exports).

## A minimal complete file

```json
{
  "format_version": 5,
  "identifier": "my_material",
  "formula": "MgO",
  "units": {"pressure": "GPa", "temperature": "K", "volume": "angstrom^3/conventional_unit_cell"},
  "eos_records": [
    {
      "identifier": "my_material_bm3",
      "model": "BM3",
      "temperature_ref": 300,
      "parameters": {
        "V0": {"value": 74.7, "error": 0.02},
        "K0": {"value": 160, "error": 2},
        "K0_prime": {"value": 4, "fixed": true}
      },
      "validity": {"pressure_gpa": [0, 30]},
      "reference": "Fictional example; no experimental validation."
    }
  ]
}
```

Save it as `my_material.eosmat`, then load and evaluate it:

```python
from peritheos import Material, load_eosmat

document = load_eosmat("my_material.eosmat")
material = Material.from_eosmat(document, require_primary_validation=False)
record = material.eos_records[0]
print(record.pressure(70.0))  # GPa at V = 70 Å³ per conventional cell
print(record.volume(10.0))    # Å³ per conventional cell at 10 GPa
```

`Material.from_eosmat` accepts only records whose source check is `verified`
unless you pass `require_primary_validation=False`; a file you write yourself
is `pending` until a curator checks it against its source.

## Units and volume

| Field | Meaning |
| --- | --- |
| `units` | Always `{"pressure": "GPa", "temperature": "K", "volume": "angstrom^3/conventional_unit_cell"}`. It may be omitted; any other value is an error. Dataset columns carry their own units. |
| `formula_units_per_cell` | Crystallographic Z: formula units per conventional cell. Required when a record uses a molar-volume thermal model (unless the record gives `volume.public_to_model_scale`) and when the material has a Hugoniot record. Never inferred from `formula`. |
| `cell_contents` | Text stating what one conventional cell holds, for example `"4 MgO formula units per conventional cubic cell"`. When absent, `Material.cell_contents` is built from Z. |
| `volume_basis` (record) | Amount of material that every `V` and `V0` of the record represents: `{"kind": "formula_units", "formula_units": Z, "molar_mass_g_mol": M}`. Required on Hugoniot records, where `formula_units` must equal `formula_units_per_cell` and V0, `rho0` and M must agree. Elsewhere it states the same basis and is not interpreted. |
| `volume` (record) | Explicit volume conversion: `public_to_model_scale` (multiplies the cell volume to give the model volume), `model_unit`, `reference_value` and `reference_meaning`. Only `public_to_model_scale` changes the calculation. |

Every volume in a record is Å³ per conventional cell. A molar-volume thermal
model receives `V × N_A × 10⁻²⁵ / Z` in J bar⁻¹ mol⁻¹; the conversion, the
meaning of `n` and worked examples are in
[Units and reference states](units.md#atom-count-and-molar-volume-basis).

## Material fields

| Field | Required | Type | Meaning |
| --- | --- | --- | --- |
| `format_version` | yes | `5` | Format version. A v5 file has no `format` key. |
| `identifier` | yes | lower snake case | Stable machine name of the material. |
| `name` | no | text | Display name; defaults to `identifier`. |
| `formula` | yes | text | Chemical formula. |
| `phase` | no | text | Phase name; `Material.phase` falls back to `symmetry`, then `"unspecified"`. |
| `aliases` | no | array of text | Earlier names and composition descriptions accepted by lookup. |
| `family_id` | no | lower snake case | Browsing family; see [Material catalog](catalog.md#browse-material-families). |
| `units`, `formula_units_per_cell`, `cell_contents` | see above | | Units and cell definition. |
| `symmetry`, `space_group`, `space_group_number` | no | text, text, integer | Crystal system (for example `CUBIC`) and space group. |
| `lattice` | no | object | `a`, `b`, `c` in Å and `alpha`, `beta`, `gamma` in degrees; `b` and `c` may be `null` when fixed by symmetry. Structural metadata; it need not reproduce any record's `V0`. |
| `atom_sites` | no | array | `element`, `wyckoff`, `x`, `y`, `z`, `occupancy` per site. |
| `peaks` | no | array | Fallback reflection list for diffraction applications, five numbers per reflection. |
| `source` | no | object | Provenance of the structure fields. |
| `notes` | no | text | Material-level remarks. |
| `recommendation_hold` | no | object | Review hold: `status`, `date`, `reason`, and `records` (array) or `record`. With `status: "on_hold"` and `scope: "equilibrium_default"`, no equilibrium default is chosen, not even the first record. A hold naming a record adds the hold sentence to its `warning`. |
| `default_records` | no | object | Default record per category: `equilibrium` (isothermal and thermal records) and `hugoniot`. See below. |
| `eos_records` | yes | array | The EOS records. An executable `Material` needs at least one. |
| `datasets` | no | array | Observation tables; see [Datasets](#datasets). |
| `shared_metadata` | no | object | Named citation, evidence or calibration objects used through `{"shared_ref": "<name>"}`; see [shared metadata](eosmat-v4.md#optional-shared-metadata-inside-the-same-file). |

`default_records` names each default once. Each value must name a record of
this file in that category. The per-record flags `default: true` and
`default_for` remain valid; when both forms are present they must agree.
Without either, `Material.default_record()` uses the first record.

## EOS record fields

Saved v5 records use this order. Fields not listed here are
[extension fields](#extension-fields), written after `aliases`; the audit
blocks `pressure_calibration` and `evidence` come last.

| Field | Required | Type and unit | Meaning |
| --- | --- | --- | --- |
| `identifier` | yes | lower snake case | Stable record name, unique within the material together with all `aliases`. |
| `label` | no | text | Display label; defaults to `identifier`. |
| `check` | no | object | Reproduction headline and source check; see [The check block](#the-check-block). |
| `model` | yes | model name | Reference isotherm or Hugoniot model (list below). |
| `temperature_ref` | no | K | Reference temperature. Absent: the Hugoniot initial temperature, else the thermal `Tr`, else 300 K. `0` only for a static isothermal record. |
| `reference_pressure_gpa` | no | GPa | Pressure at which the reference coefficients are stated. Informational; `BM2Reference` uses its parameter `P0`. |
| `parameters` | yes | object | Coefficients with errors and fixed flags; see [Parameters](#parameters). |
| `equation_options` | no | object | Extra settings of the reference component (the v3 `eos` object); may not contain `type`, `model` or `parameters`. |
| `thermal` | no | object | Thermal-pressure component; see [Thermal component](#thermal-component). |
| `parameter_error_confidence` | no | number in (0, 1) or `null` | Two-sided confidence of the reported errors when they are interval half-widths. `null` or absent: standard errors, or a convention stated in `notes`. |
| `parameter_covariance` | no | object or `null` | `parameter_order` and square `matrix`; see [Covariance](#covariance). |
| `validity` | no | object | The stated envelope of the record: `pressure_gpa` (GPa), `temperature_k` (K), `volume_ratio` (V/V0), each `[low, high]`, and `notes` (array of text). |
| `experimental_pressure_range_gpa` | no | `[low, high]` GPa | Pressure span of the data behind the fit. |
| `experimental_temperature_range_k` | no | `[low, high]` K | Temperature span of the data behind the fit. |
| `pressure_range_status` | no | enum | Where the pressure span comes from: `reported_exactly` (endpoints printed in the source), `reported_qualitatively` (source states the span without exact endpoints), `theoretical` (calculated states) or `reference_parameterization` (span of a reference or database parameterization, not of one experimental fit). |
| `range_provenance` | no | object or text | Source locations behind the ranges, usually keyed `pressure` and `temperature`. |
| `determination_method` | no | enum | How the coefficients were determined: `experimental`, `theoretical`, `hybrid` or `unknown` (the reader's value when absent). Defined in [Material catalog](catalog.md#eos-record-categories). |
| `record_kind` | no | enum | `published`, `refit`, `derived` or `diagnostic` ([record origin](eos-labels.md#record-origin)). The validator treats absence as `published`; `create_eosmat` writes `diagnostic`. A `refit` needs `fit_provenance` (in `evidence` or the sidecar); a `derived` record needs `derivation`. |
| `notes` | no | text | Record-specific qualifications. |
| `fit_datasets` | no | array of dataset identifiers | Datasets the coefficients were fitted to. |
| `supporting_datasets`, `validation_datasets`, `comparison_datasets`, `diagnostic_datasets`, `analysis_datasets`, `context_datasets` | no | arrays of dataset identifiers | Datasets that support, validate, are compared with, diagnose, were analysed for, or give context to the record. Not interpreted by the reader. |
| `reference` | no | text or object | Citation: `authors`, `year`, `title`, `source`, `volume`, `locator`, `doi`, optional `issue`, `details`, `locations`, `url`. Absent: `"User-supplied parameters"`. |
| `equation_kind` | no | enum | `isothermal`, `thermal` or `hugoniot`; when given it must match the model and the presence of `thermal`. |
| `default_for`, `default` | no | `equilibrium`/`hugoniot`; boolean | Per-record default flags; bundled files use `default_records` instead. |
| `catalog_access` | no | `explicit_selection` | Never chosen automatically: with the default `require_primary_validation=True`, `Material.from_eosmat` leaves out such a non-default record unless it is named in `record_identifiers`. Adds the explicit-selection sentence to `warning`. |
| `phase_scope`, `phase_selection` | no | text, array or object | Which phase(s) the fit covers and how the phase was assigned. Not interpreted. |
| `loading_path`, `branch_kind`, `initial_state`, `branch_domain` | Hugoniot only | | Shock path definition; see [Shock Hugoniots](hugoniots.md). |
| `volume_basis`, `volume` | see [Units and volume](#units-and-volume) | | |
| `sample_composition` | no | text or object | Composition of the measured or modelled sample. Not interpreted. |
| `experimental_configuration` | no | object | Sample, pressure cell, diffraction, calibration and similar setup notes. Not interpreted. |
| `computational_configuration` | no | object | Calculation method, code and settings of a theoretical record. Not interpreted. |
| `temperature_ref_provenance` | no | object | `{"kind": "assumed_room_temperature", "reported_condition": "room temperature", "note": "…"}`: marks `temperature_ref: 298` of an isothermal record as a nominal room temperature. |
| `pressure_scale_variant_of` | no | record identifier | Marks a published alternative pressure-scale fit of another record in this material; it shares the model, temperature and fixed parameters, is never a default and is left out of normal discovery. |
| `derived_from_record` | no | record identifier | Record that supplied model choices or fixed parameters to a refit. |
| `aliases` | no | array of text | Alternative record names accepted by `Material.get_eos_record`. |
| `derivation` | `derived` records | object | `source_kind` (`sesame_table`, `published_table` or `calculation`), `source_identifier`, `method`, optional `sampling_domain`, `software`, `source_version`, `access_and_licensing`. |
| `pressure_calibration` | no | object | How the pressures of the fit were obtained ([fields](eosmat-schema.md#pressure-calibration-and-reference-eos-records)). Bundled records keep it in the audit sidecar. |
| `evidence` | no | object | Audit trail; see [Audit sidecar](#audit-sidecar). |

v5 records must not contain the v3 names `eos`, `parameter_errors`,
`fixed_parameters` or the unit switch `parameter_error_volume_basis`, nor the
evidence blocks (`scientific_validation`, `parameter_provenance`,
`source_lineage`, `fit_provenance`, `fit_protocol`, `reproduction`,
`audit_corrections`) outside `evidence`.

**Models.** Reference and Hugoniot `model` values: `BM2`, `BM2Reference`,
`BM3`, `BM4`, `Vinet`, `Vinet3`, `Murnaghan`, `SecondOrderMurnaghan`,
`NaturalStrain2`, `NaturalStrain3`, `NaturalStrain4`, `Holzapfel`,
`ModifiedTait`, `Baonza`, `RydbergStacey`, `Morse3`, `SunMorse3`, `SunMorse4`,
`OddInversePower`, `DensityPolynomial3`, `Maltby2024` and `LinearUsUpHugoniot`.
Their parameters and equations are in [EOS models](models.md) and the
[equation reference](equation-reference.md). An unknown model is an error.

### Which range the reader uses

`EOSRecord.validity` is the envelope used by `within_calibration_range()` and
by evaluation with `check_validity=True`:

- When the record has a `validity` object, only that object is used. A range it
  omits is unbounded (temperature: any positive value), even if the record
  gives an experimental range for that quantity.
- Without `validity`, the envelope is `experimental_pressure_range_gpa` and
  `experimental_temperature_range_k`; a missing one is unbounded.

Catalog search (`search_eos_records`) instead takes each quantity separately:
`validity`, else the experimental range, else, for the temperature of an
isothermal record, `temperature_ref`. A quantity with no range never matches a
range query. Write `validity` with every quantity you want enforced.

## Parameters

**An error uses the same unit as its parameter value. A covariance entry uses
the product of the units of its two parameters.**

Each entry of `parameters` (and of `thermal.parameters`) is either a number or
an object with any of these keys:

| Key | Meaning |
| --- | --- |
| `value` | The coefficient: volumes in Å³ per conventional cell, pressures and moduli in GPa, temperatures in K, other coefficients in the units given on their model page. |
| `error` | Reported error, `≥ 0`. Absent or `null`: unknown, never zero. `0` is kept as reported and does not establish exactness. |
| `fixed` | `true` when the parameter was held fixed in the reported fit. Independent of `error`. |

An object with only `error` or `fixed: true` records information about a
parameter the equation does not store. `parameter_error_confidence` gives the
confidence level of all errors of the record.

### Covariance

```json
"parameter_covariance": {
  "parameter_order": ["V0", "K0"],
  "matrix": [[0.0004, 0.01], [0.01, 4.0]]
}
```

Here the diagonal units are Å⁶ and GPa², the off-diagonal entries Å³·GPa. In a
thermal record, `rt_eos.V0` and `rt_eos.K0` name reference-isotherm
parameters and keep the units of the record's `parameters`; unprefixed names
are thermal coefficients. When the model works in molar volume, Peritheos
scales `V0`, its error and its covariance rows and columns together: a factor
`s` on the volume gives `s` on its error and cross terms and `s²` on its
variance. Confidence, fixed flags and correlation are unchanged. A matrix says
nothing about how the source estimated it; state that in `notes`.

## Thermal component

```json
"thermal": {
  "model": "MieGruneisenDebye",
  "parameters": {"Tr": 300, "theta0": 760, "gamma0": 1.5, "q": 1, "n": 2},
  "debye_temperature_law": "integrated_gruneisen"
}
```

`model` and `parameters` are required. Thermal models: `AlphaKT`,
`LinearThermalPressure`, `LogVolumeThermalPressure`,
`SecondOrderTaylorThermalPressure`, `MieGruneisenDebye`,
`VariableQMieGruneisenDebye`, `MieGruneisenEinstein`, `QCompromiseDebye`,
`DebyeQuadraticThermalPressure`, `DebyeTabulatedThermalPressure`,
`AsymptoticDebyeTabulatedPressure`, `DebyeAnharmonicHelmholtz`,
`AsymptoticPowerLawMieGruneisenDebye`, `AsymptoticPowerLawMieGruneisenDebyeExcess`,
`Dewaele2006`, `DorogokupetsOganov2007`, `MultiOscillatorGruneisen`,
`Sokolova2016`, `HollandPowellThermalPressure`, `ThermalModifiedTait`,
`SoundVelocityDebyeHelmholtz`, `DoubleDebyeHelmholtz` and
`DoubleDebyeLogMomentHelmholtz`.

The optional law choices are `debye_temperature_law`,
`thermal_pressure_reference`, `thermal_expansion_law`, `reference_volume_law`,
`bulk_modulus_law` and `zero_point_pressure`; their allowed values and
defaults are given for each model in [Thermal models](models/thermal.md), the
linked equation pages and the [thermal equations](eosmat-schema.md#thermal-equations)
section of the v3 page, whose law definitions are unchanged in v5. A
`configuration` object passes further model settings. A value the model does
not support is an error. Other thermal keys found in bundled files
(`parameter_error_confidence`, `notes`, `gas_constant_j_mol_k`,
`parameter_constraints`) are kept and not interpreted.

## The check block

A curated record states its headline in `check`. Every field is optional, and
any other key is an error. The words are defined in [EOS labels](eos-labels.md).

| Field | Value |
| --- | --- |
| `reproduction` | `reproduced`, `close`, `disagrees` or `untested` |
| `compared` | `published_fit`, `peritheos_refit`, `source_equation` or `none`; `untested` allows only `none` or `source_equation`, other headlines exclude `none` |
| `why` | At most 160 characters, only with `close` or `disagrees`: the coefficients that failed |
| `warning` | `null` or text: an explicit usage restriction |
| `source_check` | `verified`, `pending`, `unresolved` or `deferred` |
| `source_checked_on` | Date of the source check (`YYYY-MM-DD` in bundled records) |
| `source_note` | The source check's finding |
| `audit` | `<path relative to this file>#<record identifier>`: the record's sidecar entry |

When evidence is present (inline or merged), `reproduction`, `compared` and
`source_check` must equal the values derived from it, or loading fails; the
[label page](eos-labels.md#detailed-outcome-codes) lists the stored codes they
rename. Without `check`, `EOSRecord.reproduction`, `.compared`, `.why`,
`.warning`, `.source_check` and `.source_note` are derived from the evidence;
a record without an assessment is `untested`, compared `none`. Without
evidence and without `check.source_check`, the record is `pending`.

## Datasets

Each entry of `datasets` is one observation table:

| Key | Required | Meaning |
| --- | --- | --- |
| `identifier` | yes | Unique dataset name (lower snake case). |
| `kind` | yes | Kind of table, for example `pressure_volume_temperature` or `calculated_checkpoint`. |
| `description` | no | One-line description. |
| `reference` | yes | Citation of the table's source (text or object, as for records). |
| `source_location` | yes | Exact table, figure or file in the source. |
| `source_url` | no | Where the source table can be obtained. |
| `columns` | yes | `name`, `quantity`, `unit`, `role` (`value`, `uncertainty`, `standard_deviation`, `standard_error`, `bound` or `flag`), `of` (the earlier value column an uncertainty belongs to), optional `pressure_scale`. |
| `rows` or `resource` | exactly one | Embedded rows (numbers or `null`, one per column), or a file `{path, sha256, media_type}` with a relative path. |
| `used_by_eos_records` | yes | Records in this material that use the table; `[]` for a supporting table with no EOS. |
| `uncertainty` | no | How errors and covariance were reported, for example `reported_as` and `covariance`. |
| `notes` | no | Remarks. |
| `license` | no | Redistribution terms of the table. |
| `pressure_reductions`, `observation_views` | no | Per-EOS pressure coordinates and selected observation views; see [Dataset pressure reductions](dataset-pressure-reductions.md). |

`Material.get_dataset(identifier)` returns typed columns
([Loading observation datasets](datasets.md)). It resolves `resource.path`
inside the installed package; for a table next to your own file use
`Dataset.from_mapping(mapping, resource_root=folder)`. Other dataset keys are
kept in `Dataset.metadata` and listed under [extension fields](#extension-fields).

## Audit sidecar

A curated material can keep its audit trail in a second file. This record file,
`materials/example_oxide.eosmat`, names its entry in `check.audit`:

<!-- example: record file -->
```json
{
  "format_version": 5,
  "identifier": "example_oxide",
  "formula": "MgO",
  "eos_records": [
    {
      "identifier": "example_oxide_bm3",
      "check": {"reproduction": "reproduced", "compared": "published_fit", "warning": null, "source_check": "pending",
                "audit": "../audits/example_oxide.audit.json#example_oxide_bm3"},
      "model": "BM3",
      "temperature_ref": 300,
      "parameters": {"V0": {"value": 74.7, "error": 0.02}, "K0": {"value": 160, "error": 2}, "K0_prime": {"value": 4, "fixed": true}},
      "validity": {"pressure_gpa": [0, 30]},
      "reference": "Fictional example; no experimental validation."
    }
  ]
}
```

The sidecar `audits/example_oxide.audit.json` holds the entry, keyed by record:

<!-- example: audit sidecar -->
```json
{
  "format": "peritheos.eosmat-audit",
  "format_version": 5,
  "material": "example_oxide",
  "records": {
    "example_oxide_bm3": {
      "evidence": {
        "scientific_validation": {"status": "pending_primary_source_check"},
        "reproduction": {"assessment": {
          "status": "parity",
          "scope": "parameter_agreement",
          "basis": "Fictional: a refit of example observations recovers V0 and K0 under the generic criterion.",
          "parameters": [
            {"parameter": "V0", "published": 74.7, "published_error": 0.02, "refit": 74.71, "refit_error": 0.02},
            {"parameter": "K0", "published": 160, "published_error": 2, "refit": 159.2, "refit_error": 1.8}
          ],
          "fit": {"free_parameters": ["V0", "K0"], "fixed_parameters": ["K0_prime"]}
        }}
      }
    }
  }
}
```

```python
from peritheos import Material, load_eosmat

document = load_eosmat("materials/example_oxide.eosmat")  # merges the sidecar
record = Material.from_eosmat(document, require_primary_validation=False).eos_records[0]
print(record.reproduction, record.warning)  # reproduced None
```

Rules:

- `material` equals the record file's `identifier`; every key of `records`
  names a record of that file.
- An entry may contain only `evidence`, `pressure_calibration`,
  `fit_input_inventory`, `network_reconstruction_audit`, `reproduction_audit`,
  `numerical_reproduction`, `validation_selection`, `migration_corrections`,
  `scope_history` and `author_input_followup`. Each keeps its inline meaning.
- A record uses its entry when `check.audit` names it and the record holds none
  of these blocks inline. The path is resolved relative to the record file; the
  fragment must equal the record identifier.
- Merging gives the complete inline record. A block present in both places is
  an error. Absent `report_url`, `report_path` and `report` of an assessment are
  filled in from the material and record identifiers.

The fields inside `evidence.reproduction.assessment` are described in
[EOSMAT v4](eosmat-v4.md#a-consistent-reconstruction-assessment) and the
outcome codes in [EOS labels](eos-labels.md#detailed-outcome-codes). How each
call treats a sidecar is listed under
[loading and sidecars](eosmat-v5.md#loading-missing-sidecars-and-exports).

## Extension fields

The schema leaves material, record and dataset objects open. Peritheos keeps
every key it does not interpret: `load_eosmat`, `save_eosmat` and `format`
preserve it, `Material.eosmat_metadata` and `EOSRecord.eosmat_metadata` expose
it, and `Dataset.metadata` holds extra dataset keys. Keys used in the bundled
catalog:

| Level | Keys | Meaning |
| --- | --- | --- |
| material | `structure_provenance`, `structure_reference` | Source of the structure fields. |
| material | `sample_composition`, `composition` | Reported sample composition and its analysis. |
| material | `crystal_system` | Crystal system name. |
| material | `supporting_studies`, `phase_references` | Further studies and phase-identification sources. |
| material | `phase_scope`, `integration_status`, `curation_disposition` | Which conditions and records the curation accepts, and why. |
| record | `spin_crossover_context`, `spin_path_context`, `iron_speciation`, `iron_valence` | Spin branch, transition data, compression path and iron state of the sample. |
| record | `constrained_fit_family` | Membership in a published family of constrained fits (family identifier, source row, preferred member). |
| record | `implicit_parameters`, `published_constraints` | Values the source held fixed without printing them as coefficients. |
| record | `parameter_uncertainty_convention`, `parameter_error_convention`, `uncertainties`, `uncertainty_provenance` | What the printed errors mean and where they come from. |
| record | `source_parameters`, `source_parameterization`, `source_compressibility_coefficients`, `source_quadratic_fit`, `source_coefficient_rounding`, `source_volume_normalization`, `source_heat_capacity` | The source's values in its own form or units, and how they map to `parameters`. |
| record | `published_fit_provenance`, `published_reference_fit`, `published_reference_state`, `reference_state` | The source's own fit description or reference state. |
| record | `source_fit_domain`, `source_fit_pressure_range_gpa`, `path_interpretation` | Domain and interpretation of the source's fit. |
| record | `fit_data_sources`, `recommended_refit_dataset` | Origin of the fitted data and the dataset preferred for a refit. |
| record | `author_fit_preference`, `supporting_references` | The authors' preferred fit; further citations. |
| dataset | `provenance`, `provenance_type`, `provenance_resource`, `provenance_only`, `artifact_provenance`, `transcription`, `digitization`, `processing`, `quality_control`, `generation_script` | How the table was obtained, transcribed or processed. |
| dataset | `source_doi`, `source_version`, `source_document`, `source_artifact`, `source_resource`, `original_resource`, `source_file_sha256`, `source_container_sha256`, `source_table_sha256`, `source_sha256`, `source_tokens`, `source_provenance`, `source_lineage`, `upstream_sources`, `study_identifier` | Identification and checksums of the source files and printed values. |
| dataset | `license_url`, `license_scope`, `license_notice`, `source_license`, `rights`, `source_rights`, `contributor_rights`, `source_attribution` | Redistribution terms. |
| dataset | `derivation`, `derived_from_dataset`, `derived_from_datasets`, `source_datasets`, `pressure_reconstruction`, `sample_length_corrections`, `quantity_origins` | How derived columns or tables were computed. |
| dataset | `related_eos_records`, `usage_recommendation`, `preferred_for_refitting`, `recommendation_reason`, `selection`, `supporting_references` | Intended use of the table and row selection. |
| dataset | `row_count`, `phase`, `determination_method`, `uncertainty_convention` | Row count, phase, origin and error convention of the table. |

None of these keys changes a calculation.
