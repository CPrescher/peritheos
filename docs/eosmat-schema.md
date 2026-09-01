# `.eosmat` schema reference

`.eosmat` is Peritheos's material exchange format. A single flat JSON document
can carry material identity, optional diffraction structure, and any number of
literature EOS records. Peritheos defines the scientific field semantics;
applications such as Dioptas may add or use the optional crystallographic
fields.

The canonical format is `peritheos.material` version 3. Its
[normative machine-readable definition](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/eosmat-v3.schema.json)
is distributed as `peritheos/data/eosmat-v3.schema.json` and uses JSON Schema
draft 2020-12. This page explains the fields and consumer behavior that cannot
be conveyed by shape validation alone.

## Top-level document

| Field | Required | Meaning |
|---|---:|---|
| `format` | yes | Must be `peritheos.material`. |
| `format_version` | yes | Must be `3`. |
| `identifier` | yes | Stable lower-snake-case material identifier. |
| `name` | yes | Human-readable material or phase name. |
| `formula` | yes | Chemical formula as text. |
| `units` | yes | Fixed public exchange units described below. |
| `eos_records` | yes | Array of zero or more EOS records. |
| `phase`, `aliases`, `notes` | no | Additional material identity and description. |
| `symmetry`, `lattice` | no | Unit-cell symmetry and reference lattice. |
| `formula_units_per_cell` | no | Crystallographic $Z$ for cell-to-molar conversion. |
| `space_group`, `space_group_number` | no | Optional space-group metadata. |
| `atom_sites`, `peaks`, `source` | no | Optional diffraction structure, fallback peaks, and source metadata. |

`eos_records` may be empty in a structure-only exchange file. Peritheos's
bundled catalog intentionally contains only materials with at least one EOS
record.

Unknown optional fields must be preserved by a round-trip-capable reader. They
must never be treated as Python import paths or executable code.

## Fixed exchange units

Every canonical document contains:

```json
"units": {
  "pressure": "GPa",
  "temperature": "K",
  "volume": "angstrom^3/conventional_unit_cell"
}
```

`eos.parameters.V0` therefore uses conventional-cell Å³ in the exchange file.
For a molar-energy thermal model, the consumer converts it using

\[
V_{\mathrm{model}}
=V_{\mathrm{cell}}\frac{N_A10^{-25}}{Z}
\quad\mathrm{J\,bar^{-1}\,mol^{-1}},
\]

where `formula_units_per_cell` is $Z$. A consumer must reject a thermal record
that requires molar volume when $Z$ is unavailable; guessing $Z$ from the
formula is not valid.

## EOS record

Each item in `eos_records` describes one source parameterization.

| Field | Required | Meaning |
|---|---:|---|
| `identifier` | yes | Stable lower-snake-case record identifier. |
| `label` | yes | Human-readable source/model label. |
| `reference` | yes | Citation string or structured citation with authors, year, source, and optional DOI. |
| `default` | no | Preferred record for this material; at most one may be true. |
| `eos` | yes | Reference-isotherm equation and parameters. |
| `thermal` | no | Thermal-pressure equation, parameters, and fixed equation choices. |
| `parameter_errors` | yes | Reported reference-isotherm parameter errors; JSON `null` means unavailable, not zero. |
| `parameter_error_confidence` | no | Two-sided confidence level when the reported errors are interval half-widths rather than standard errors. |
| `fixed_parameters` | yes | Reference-isotherm parameters fixed during the reported fit. |
| `experimental_pressure_range_gpa` | no | Two-element marginal pressure envelope. |
| `pressure_range_status` | no | Provenance of the pressure envelope: `reported_exactly`, `reported_qualitatively`, `theoretical`, or `reference_parameterization`. |
| `experimental_temperature_range_k` | no | Two-element marginal temperature envelope. |
| `temperature_ref` | no | Reference-isotherm temperature in K. |
| `parameter_provenance` | no | Field-level table, equation, page, or supplement provenance. |
| `source_lineage` | no | Ordered sources and their roles when a record combines an earlier fit, final parameter table, implementation, correction, or experimental context. |
| `parameter_covariance` | no | Published covariance metadata when available. |
| `scientific_validation` | yes | Validation boundary described below. |
| `notes` | no | Scientific qualifications that do not fit another field. |

An experimental range records data or fit coverage, not an assertion of phase
stability and not permission to extrapolate over every combination of its
marginal bounds. `reported_exactly` means the numerical endpoints occur in the
primary source; it does not imply a recommended extrapolation limit or exact
phase-stability boundary.

### Material identity versus EOS reference state

The top level describes the material and may contain crystallography used by a
diffraction application. Each EOS record independently defines a published
parameterization and its reference-volume convention. Consequently:

- `lattice` is structural metadata and need not numerically reproduce an EOS
  record's `V0`; the record-level `parameter_provenance` states which value is
  used by the equation;
- one material may have several literature records, and may also have several
  records that differ only in a scientifically meaningful reference state;
- a record spanning more than one phase is permitted only when the primary
  source itself publishes such a fit and the scope is explicit in its label,
  phase description, notes, and validity information;
- cell contents and volume convention belong to the material/record contract,
  not to the equation family name.

For example, the phase-D file has one crystallographic interchange card but
two BM2 records because Shieh et al. report distinct AntA and AntB ambient
volumes. The lithium record preserves Hanfland et al.'s single empirical Vinet
fit across bcc and fcc observations while recording `V0` in the conventional
two-atom bcc-cell convention. Neither case is flattened into a fictitious
universal reference volume.

## Reference-isotherm equations

Both `type` and `model` are required. `type` preserves the concise
application-facing name used by the original Dioptas format; `model` is the
stable mechanism-oriented Peritheos identifier. They must be paired exactly:

| `type` | `model` | Principal parameters |
|---|---|---|
| `BM2` | `birch_murnaghan_2` | `V0`, `K0` |
| `BM3` | `birch_murnaghan_3` | `V0`, `K0`, `K0_prime` |
| `BM4` | `birch_murnaghan_4` | `V0`, `K0`, `K0_prime`, `K0_prime_prime` |
| `Vinet` | `vinet` | `V0`, `K0`, `K0_prime` |
| `Murnaghan` | `murnaghan` | `V0`, `K0`, `K0_prime` |
| `Holzapfel` | `holzapfel` | `V0`, `K0`, `K0_prime`, `n`, `Z` |
| `ModifiedTait` | `modified_tait` | `V0`, `K0`, `K0_prime`, `K0_double_prime` |
| `NaturalStrain2` | `natural_strain_2` | `V0`, `K0` |
| `NaturalStrain3` | `natural_strain_3` | `V0`, `K0`, `K0_prime` |
| `NaturalStrain4` | `natural_strain_4` | `V0`, `K0`, `K0_prime`, `K0_double_prime` |

The complete definitions and domains are in the
[equation reference](equation-reference.md#isothermal-equations).

## Thermal equations

Thermal `type` and `model` must likewise match:

| `type` | `model` | Main parameters |
|---|---|---|
| `AlphaKT` | `thermal_reference_state` | `Tr`, `alpha0`, `dK_dT`; optional `alpha1` |
| `LinearThermalPressure` | `linear_thermal_pressure` | `Tr`, `alpha_KT` |
| `LogVolumeThermalPressure` | `log_volume_thermal_pressure` | `Tr`, `alpha_KT_ref`, `dK_dT_V` |
| `MieGruneisenDebye` | `mie_gruneisen_debye` | `Tr`, `theta0`, `gamma0`, `q`, `n` |
| `MieGruneisenEinstein` | `mie_gruneisen_einstein` | `Tr`, `theta0`, `gamma0`, `q`, `n` |
| `AsymptoticPowerLawMieGruneisenDebye` | `asymptotic_power_law_mie_gruneisen_debye` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n` |
| `MultiOscillatorGruneisen` | `multi_oscillator_gruneisen_thermal_pressure` | Oscillator, Grüneisen, anharmonic, and electronic parameters |
| `ThermalModifiedTait` | `thermal_modified_tait` | `Tr`, `theta`, `alpha0`, `n` |

`Sokolova2016` remains accepted as a legacy type alias for
`multi_oscillator_gruneisen_thermal_pressure`; new writers use the general
`MultiOscillatorGruneisen` type.

`AlphaKT` is retained as the Dioptas-facing interchange type. Its canonical
model identifier is the mechanism-oriented `thermal_reference_state`, and the
corresponding Peritheos class is `ThermalReferenceStateEOS`. It evaluates a
temperature-dependent reference volume and bulk modulus; it is not the
constant-`alpha_KT` pressure increment represented by `LinearThermalPressure`.

### Thermal-expansion law

An `AlphaKT` component may select the integrated volumetric expansion law:

```json
"thermal_expansion_law": "linear_temperature"
```

Allowed values are `constant` and `linear_temperature`. Omission means
`constant` for backward compatibility. The linear law requires `alpha1` and
represents $\alpha(T)=\alpha_0+\alpha_1T$; the implementation analytically
integrates this expression when constructing $V_0(T)$. Writers should store
the field explicitly for new or edited linear-temperature records. The field
is valid only for `AlphaKT` / `thermal_reference_state` components.

The independent `reference_volume_law` configuration controls how expansion
information constructs the reference volume. Omission means
`integrated_expansivity`, which uses the exponential integral above. The
alternative

```json
"reference_volume_law": "linear_temperature"
```

applies the direct relation
$V_0(T)=V_0(T_r)[1+\alpha_0(T-T_r)]$. Here `alpha0` is a mean expansion
coefficient over the represented interval, not a constant instantaneous
expansivity. This option requires `thermal_expansion_law="constant"` (or its
omission) and `alpha1=0`. It represents Martinez et al. (1996), Equation 3,
without replacing the paper's linear relation by an exponential approximation.

Thermal records may contain their own `parameter_errors` and
`fixed_parameters`. A missing error and an explicitly fixed parameter are
different facts and must remain distinguishable.

### Debye-temperature law

`MieGruneisenDebye` accepts an additional fixed field:

```json
"debye_temperature_law": "integrated_gruneisen"
```

Allowed values are:

| Value | Characteristic-temperature relation |
|---|---|
| `integrated_gruneisen` | $\Theta=\Theta_0\exp[-\gamma_0(x^q-1)/q]$, with the continuous $q=0$ limit |
| `variable_exponent` | $\Theta=\Theta_0x^{-\gamma(V)}$, where $\gamma(V)=\gamma_0x^q$ |

Here $x=V/V_0$. If `debye_temperature_law` is absent, it means
`integrated_gruneisen`. The JSON Schema `default` keyword documents this
semantic default but does not insert the field into a loaded JSON object.
Writers should store the field explicitly for new or edited MGD records.

The law is fixed equation configuration rather than a fitted numeric
parameter. Consumers must honor an explicit value or reject the record. They
must not silently replace `variable_exponent` with the integrated default.

The complete equations and the distinction between the two laws are in the
[thermal equation reference](equation-reference.md#mie-gruneisen-debye-and-einstein).

## Scientific validation status

Structural validity does not establish scientific validity. Every canonical
record declares one of:

| Status | Meaning |
|---|---|
| `primary_source_validated` | Equation, parameters, units, and reference state were checked against the primary publication or official supplement. |
| `pending_primary_source_check` | Structurally usable catalog data that has not completed that audit. |
| `deferred` | The record is intentionally unavailable or incomplete because its primary evidence is inaccessible, ambiguous, or inconsistent. |

Agreement with another software library is not sufficient for
`primary_source_validated`.

Bundled records additionally carry `audit_date`, a `primary_source_check`
object with DOI/URL and equation-table-page locations, and either
`verified_fields` or `unresolved`. These are additive extension fields. The
record-by-record package ledger is
`peritheos/data/primary-source-audit.json`. As of the 2026-09-01 audit, all 147
bundled records are validated, with no deferred or pending record.

## Complete EOS-only example

Optional diffraction fields are omitted here; adding them does not change the
EOS record semantics.

```json
{
  "format": "peritheos.material",
  "format_version": 3,
  "identifier": "gold",
  "name": "Gold",
  "formula": "Au",
  "phase": "fcc",
  "units": {
    "pressure": "GPa",
    "temperature": "K",
    "volume": "angstrom^3/conventional_unit_cell"
  },
  "formula_units_per_cell": 4,
  "eos_records": [
    {
      "identifier": "gold_fei_2007_vinet",
      "label": "Fei et al. (2007), Vinet",
      "reference": {
        "authors": ["Fei", "Ricolleau", "Frank", "Mibe", "Shen", "Prakapenka"],
        "year": 2007,
        "source": "Proc. Natl. Acad. Sci. USA",
        "doi": "10.1073/pnas.0609013104"
      },
      "default": true,
      "eos": {
        "type": "Vinet",
        "model": "vinet",
        "parameters": {
          "V0": 67.85,
          "K0": 167.0,
          "K0_prime": 6.0
        }
      },
      "parameter_errors": {
        "V0": 0.004,
        "K0": null,
        "K0_prime": 0.02
      },
      "fixed_parameters": ["K0"],
      "thermal": {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "variable_exponent",
        "parameters": {
          "Tr": 300.0,
          "theta0": 170.0,
          "gamma0": 2.97,
          "q": 0.6,
          "n": 1
        }
      },
      "experimental_pressure_range_gpa": [0.0, 89.0],
      "scientific_validation": {
        "status": "primary_source_validated",
        "note": "Equation 3 and Table 1 checked against the primary publication."
      }
    }
  ]
}
```

## Validation and I/O

Peritheos exposes the bundled schema and a stricter structural validator:

```python
from peritheos import Material, eosmat_schema, load_eosmat, validate_eosmat_document

schema = eosmat_schema()
document = load_eosmat("gold.eosmat")
validate_eosmat_document(document)
material = Material.from_eosmat(document)
```

`load_eosmat()` preserves an omitted default field rather than mutating the
document. Consumers constructing an EOS apply `integrated_gruneisen` when the
Debye-temperature field is absent.

The Python validator additionally checks ordered ranges, unique record
identifiers, at most one default record, exact `type`/`model` pairing, and that
`debye_temperature_law` appears only on `MieGruneisenDebye`.

## Compatibility rules

- Peritheos format 3 is the canonical write format.
- Native Dioptas 0.10 format-2 documents remain accepted as legacy input.
- Unknown optional fields are preserved, but unknown equation types, models,
  or law values must be rejected before calculation.
- Dioptas 0.10 does not interpret `debye_temperature_law`; the current Dioptas
  development implementation preserves format-3 fields and passes the law to
  Peritheos 0.5+.
- Changing field meaning, required units, or a stable identifier requires a
  documented migration; see [API stability](api-stability.md#eosmat-compatibility).

See [Dioptas integration](dioptas-integration.md) for application ownership,
the bundled-library migration, and current round-trip limitations.
