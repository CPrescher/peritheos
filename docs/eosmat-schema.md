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
| `datasets` | no | Primary experimental tables used by one or more EOS records. |

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

## Primary experimental datasets

`datasets` makes the observations behind a fitted EOS part of the material's
scientific provenance. Small tables should be embedded in `rows`, making a
standalone `.eosmat` self-contained. A large table may instead use a relative,
checksummed `resource`. Exactly one storage form is allowed.

Each dataset declares a stable `identifier`, `kind`, primary `reference`, exact
`source_location`, typed `columns`, and the `used_by_eos_records` identifiers.
Uncertainty columns explicitly identify the value column to which they apply.
The record should state whether the source reports standard deviations,
standard errors, bounds, and covariance information; absence of covariance
must not be interpreted as statistical independence.

```json
"datasets": [
  {
    "identifier": "example_2026_table2_pvt",
    "kind": "pressure_volume_temperature",
    "reference": {
      "authors": ["Example"],
      "year": 2026,
      "source": "Journal supplement",
      "doi": "10.example/data"
    },
    "source_location": "official supplement, Table 2",
    "columns": [
      {"name": "pressure_gpa", "quantity": "pressure", "unit": "GPa", "role": "value"},
      {"name": "pressure_sigma_gpa", "quantity": "pressure", "unit": "GPa", "role": "standard_deviation", "of": "pressure_gpa"},
      {"name": "volume_a3", "quantity": "conventional_unit_cell_volume", "unit": "angstrom^3", "role": "value"},
      {"name": "temperature_k", "quantity": "temperature", "unit": "K", "role": "value"}
    ],
    "rows": [[10.0, 0.2, 101.4, 300.0]],
    "used_by_eos_records": ["example_2026_bm3_1"],
    "uncertainty": {
      "reported_as": "standard_deviation",
      "covariance": "not_reported"
    }
  }
]
```

Dataset identifiers and column names use lower snake case. Every row must have
the declared number of columns, every value must be finite or JSON `null`, and
every `used_by_eos_records` entry must resolve within the same material.

## EOS record

Each item in `eos_records` describes one source parameterization.

| Field | Required | Meaning |
|---|---:|---|
| `identifier` | yes | Stable lower-snake-case record identifier. |
| `label` | yes | Human-readable source/model label. |
| `aliases` | no | Unique, non-empty historical or alternate record names accepted by executable lookup. An alias cannot collide with another record identifier or alias in the material. |
| `reference` | yes | Citation string or structured citation with authors, year, source, and optional DOI. |
| `default_for` | no | `equilibrium` or `hugoniot`; at most one default is allowed in each category. Legacy `default: true` is interpreted within the record's category. |
| `equation_kind` | required for Hugoniots | `isothermal`, `thermal`, or `hugoniot`; older equilibrium records may omit it and are inferred. |
| `loading_path` | for Hugoniots | `principal` or `precompressed`; reshocks require a future frame-aware model. |
| `branch_kind` | for Hugoniots | `untransformed` or `transformed`, independently of loading history. |
| `initial_state` | for Hugoniots | Required precursor material identifier, phase, temperature, pressure, and density; pressure and density must equal executable `P0` and `rho0`. |
| `volume_basis` | for Hugoniots | Required operational mass basis containing `formula_units` and `molar_mass_g_mol`; these values must be consistent with `V0` and `rho0`. |
| `branch_domain` | for Hugoniots | Required particle-velocity interval, scientific meaning, and boundary status; record-level evaluation enforces it. |
| `record_kind` | no | `published`, `refit`, `derived`, or `diagnostic`; omission means `published` for backward compatibility. |
| `derived_from_record` | for refits | Identifier of the published or prior record supplying the model choices and fixed parameters. |
| `fit_provenance` | for refits | Software/version, dataset, row selection, objective, weights, varied/fixed parameters, and fit statistics. |
| `derivation` | for derived records | Structured source kind and identifier, transformation method, sampling domain, software, and access/licensing information. |
| `eos` | yes | Primary equation and parameters: an equilibrium reference isotherm or a Hugoniot path model. |
| `thermal` | no | Thermal-pressure equation, parameters, and fixed equation choices. |
| `parameter_errors` | yes | Reported primary-equation parameter errors; JSON `null` means unavailable, not zero. |
| `parameter_error_confidence` | no | Two-sided confidence level when the reported errors are interval half-widths rather than standard errors. |
| `fixed_parameters` | yes | Primary-equation parameters fixed during the reported fit. |
| `experimental_pressure_range_gpa` | no | Two-element marginal pressure envelope. |
| `pressure_range_status` | no | Provenance of the pressure envelope: `reported_exactly`, `reported_qualitatively`, `theoretical`, or `reference_parameterization`. |
| `experimental_temperature_range_k` | no | Two-element marginal temperature envelope. |
| `temperature_ref` | no | Record reference/default temperature in K. It is normally the reference isotherm; an explicitly absolute-cold-curve thermal model documents its separate pressure baseline. Zero is allowed only for a static 0 K isothermal record; thermal and Hugoniot records require a positive temperature. |
| `parameter_provenance` | no | Field-level table, equation, page, or supplement provenance. |
| `source_lineage` | no | Ordered sources and their roles when a record combines an earlier fit, final parameter table, implementation, correction, or experimental context. |
| `pressure_calibration` | no | Audited pressure basis of the observations used for the fit, including resolvable links to reference EOS and optical-calibration records. |
| `parameter_covariance` | no | Published or reproducibly derived covariance metadata when available; its origin must be documented. |
| `scientific_validation` | yes | Validation boundary described below. |
| `notes` | no | Scientific qualifications that do not fit another field. |

Hugoniots are EOS records, not a parallel material collection. Their `eos`
component has type `LinearUsUpHugoniot` and parameters `V0`, `rho0`, `c0`, `s`,
and `P0`. They cannot contain a thermal component because temperature is not an
independent coordinate on the constrained shock path. See
[Shock Hugoniot equations of state](hugoniots.md).

An experimental range records data or fit coverage, not an assertion of phase
stability and not permission to extrapolate over every combination of its
marginal bounds. `reported_exactly` means the numerical endpoints occur in the
primary source; it does not imply a recommended extrapolation limit or exact
phase-stability boundary.

### Published records and refits

A refit is stored as a separate, opt-in EOS record. It never silently replaces
the source-reported parameters. Its identifier should use a `_refit` suffix,
`record_kind` must be `refit`, and `derived_from_record` must resolve within the
same material. `fit_provenance.dataset` must likewise resolve to an embedded or
checksummed primary dataset in that material. The fit metadata records enough
of the numerical experiment to distinguish software, objective, weighting,
row selection, varied parameters, and fixed assumptions.

For example, B4C retains the source-reported
`b4c_somayazulu_2023_berman_2` record and offers the independently reproduced
`b4c_somayazulu_2023_berman_refit` record explicitly:

```python
from peritheos import get_material

b4c = get_material("b4c")
refit = b4c.get_eos_record("b4c_somayazulu_2023_berman_refit")
pressure = refit.pressure(volume=289.1, temperature=2023.0)
```

### Pressure calibration and reference EOS records

`pressure_calibration` records how the pressures underlying a published fit
were obtained. This is separate from `reference`, which cites the publication
that reports the fitted material EOS, and from `eos`, which is the material EOS
itself. See [Pressure-scale normalization](pressure-scale-normalization.md) for
how these links may be composed into valid same-standard, ruby, and
cross-material conversion paths.

The object has a `status`, a list of `methods`, and a `recalculation` status.
Each method identifies its `kind`, the exact location where the material source
describes it, and the primary citation for the calibration. When the method is
an X-ray pressure standard and the exact calibration is bundled, it also has a
`reference_eos_record` containing a globally unique EOS-record identifier.
Ruby methods instead use `reference_calibration_record` to identify an
executable R1 wavelength calibration in the bundled calibration library.
Diamond-anvil Raman-edge calibrations and literature cross-calibration edges
are stored in the same registry and can participate in recursive conversion
paths:

```json
"pressure_calibration": {
  "status": "resolved",
  "methods": [
    {
      "kind": "equation_of_state",
      "reference_eos_record": "gold_fei_2007_vinet_2",
      "reference": {
        "authors": ["Fei", "Ricolleau", "Frank", "Mibe", "Shen", "Prakapenka"],
        "year": 2007,
        "source": "Proc. Natl. Acad. Sci. USA",
        "doi": "10.1073/pnas.0609013104"
      },
      "source_location": "Experimental methods, pressure determination"
    }
  ],
  "recalculation": {
    "status": "missing_calibrant_observations",
    "notes": "The reference EOS is executable, but the source does not tabulate the simultaneously measured calibrant volumes."
  },
  "audit_date": "2026-09-03"
}
```

Allowed method kinds distinguish `equation_of_state` from optical ruby or
other fluorescence gauges, shock-wave and ultrasonic measurements, ab initio
calculations, simultaneous self-consistent reductions, and ambient-pressure
anchors. A ruby scale is not mislabeled as a material EOS. For example,
`ruby_dorogokupets_oganov_2007` executes the quadratic normalized-shift form
published in that paper, whereas `ruby_holzapfel_2005` executes the modified
Freund--Ingalls form and its distinct three coefficients.

The recalculation status is deliberately independent of reference resolution.
It describes exact observation-level re-reduction and reproduction of the
published fit. That operation requires row-wise sample data and, for a measured
XRD marker, its paired volume and temperature. `ready` means those observables
are present in the material's `datasets`;
`missing_calibrant_observations` means the scale is known but the required
row-wise readings are unavailable. Other statuses identify an
unbundled reference, an unsupported reference model, a non-applicable primary
measurement, or a source from which recalculation is impossible.

A derived pressure-scale transformation has different requirements. The
published sample EOS can be evaluated on a volume grid, and its pressure axis
can be mapped through documented standard EOSs. XRD-to-XRD transformations use
the same virtual standard volume on source and target EOSs; ruby-to-XRD paths
use an XRD-standard EOS calibrated on the ruby scale as a bridge. This supports
comparable normalized EOS curves without pretending to reproduce the original
observation-level fit or its uncertainty.

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
| `SecondOrderTaylorThermalPressure` | `second_order_taylor_thermal_pressure` | `Tr`, `eta0`, `c0`, `c1`, `c2`, `c3`, `c4`, `c5` |
| `LogVolumeThermalPressure` | `log_volume_thermal_pressure` | `Tr`, `alpha_KT_ref`, `dK_dT_V` |
| `MieGruneisenDebye` | `mie_gruneisen_debye` | `Tr`, `theta0`, `gamma0`, `q`, `n` |
| `MieGruneisenEinstein` | `mie_gruneisen_einstein` | `Tr`, `theta0`, `gamma0`, `q`, `n` |
| `AsymptoticPowerLawMieGruneisenDebye` | `asymptotic_power_law_mie_gruneisen_debye` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n` |
| `DorogokupetsOganov2007` | `dorogokupets_oganov_2007` | Four oscillator modes, `gamma0`, `gamma_inf`, `beta`, anharmonic, electronic, and defect terms |
| `MultiOscillatorGruneisen` | `multi_oscillator_gruneisen_thermal_pressure` | Oscillator, Grüneisen, anharmonic, and electronic parameters |
| `ThermalModifiedTait` | `thermal_modified_tait` | `Tr`, `theta`, `alpha0`, `n` |
| `DoubleDebyeHelmholtz` | `double_debye_helmholtz` | Double-Debye coefficients; optional `Tr` |
| `DoubleDebyeLogMomentHelmholtz` | `double_debye_log_moment_helmholtz` | Logarithmic-moment double-Debye coefficients; optional `Tr` |

`Sokolova2016` remains accepted as a legacy type alias for
`multi_oscillator_gruneisen_thermal_pressure`; new writers use the general
`MultiOscillatorGruneisen` type.

`AlphaKT` is retained as the Dioptas-facing interchange type. Its canonical
model identifier is the mechanism-oriented `thermal_reference_state`, and the
corresponding Peritheos class is `ThermalReferenceStateEOS`. It evaluates a
temperature-dependent reference volume and bulk modulus; it is not the
constant-`alpha_KT` pressure increment represented by `LinearThermalPressure`.

`SecondOrderTaylorThermalPressure` represents an absolute polynomial thermal
pressure added to a cold curve. Its numeric `Tr` is the polynomial's
temperature coordinate, not a claim that the reference EOS is an isotherm at
that temperature.

For the two double-Debye Helmholtz types, `Tr` is always written explicitly.
`"Tr": null` means that the stored Vinet curve is the motionless-ion 0 K cold
curve and no non-cold contribution is subtracted. A positive numeric `Tr`
means that the Vinet curve is the complete isotherm at that temperature; the
ionic and anharmonic free energy and pressure at `Tr` are subtracted before the
thermal contribution is added. `null` is not accepted for any other thermal
parameter or model.

### Thermal-expansion law

An `AlphaKT` component may select the integrated volumetric expansion law:

```json
"thermal_expansion_law": "linear_temperature"
```

Allowed values are `constant`, `linear_temperature`, and
`linear_reference_temperature`. Omission means `constant` for backward
compatibility. Both linear laws require `alpha1`. `linear_temperature`
represents $\alpha(T)=\alpha_0+\alpha_1T$;
`linear_reference_temperature` represents
$\alpha(T)=\alpha_0+\alpha_1(T-T_r)$ and therefore defines `alpha0` at `Tr`.
The implementation analytically integrates the selected expression when
constructing $V_0(T)$. Writers should store the field explicitly for new or
edited linear-temperature records. The field is valid only for `AlphaKT` /
`thermal_reference_state` components.

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

The third option,

```json
"thermal_expansion_law": "linear_temperature",
"reference_volume_law": "berman"
```

applies EosFit7's Berman (1988) truncated polynomial
$V_0(T)=V_0(T_r)[1+\alpha_0(T-T_r)+\tfrac12\alpha_1(T-T_r)^2]$.
Here `alpha0` is the expansion coefficient at `Tr`; this is distinct from
EosFit's exact exponential-integral Fei form.

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

### MGD thermal-pressure baseline

`MieGruneisenDebye` also accepts the fixed field
`thermal_pressure_reference`. Its default, `reference_temperature`, adds the
Debye energy difference $E(V,T)-E(V,T_r)$ to a measured reference isotherm.
The alternative `absolute_zero` adds $E(V,T)$ to an explicitly 0 K cold curve.
These conventions are not interchangeable: a record using `absolute_zero`
must trace its cold coefficients and baseline directly to the source. `Tr`
remains positive and records the API's reference/default temperature even
when it is not subtracted from total thermal pressure.

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
`peritheos/data/primary-source-audit.json`. All 226 bundled records are
validated, with no deferred or pending record.

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

The Python API exposes the bundled schema and a stricter structural validator:

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

Rust validates the same document-level invariants and additionally constructs
every built-in model before accepting or writing a document:

```rust,no_run
use peritheos::{load_eosmat, load_eosmat_str};

let mut material = load_eosmat("gold.eosmat")?;
material.document["application_note"] = "calibration copy".into();
material.validate()?;

let json = material.to_json()?;
let round_tripped = load_eosmat_str(&json)?;
round_tripped.save("gold-copy.eosmat")?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

For a decoded `serde_json::Value`, the equivalent free functions are
`validate_eosmat_document`, `serialize_eosmat`, and `save_eosmat`. Both Rust
and Python check ordered ranges, unique record identifiers, at most one
default record, exact `type`/`model` pairing, and equation-specific choices.

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
