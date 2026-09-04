# Adding a material or equation of state

An addition to the Peritheos material library is a scientific claim, not only
a data-entry change. A reviewer must be able to identify the represented
material, trace every equation and parameter to its source, execute the same
formulation in Peritheos, and understand where the result is valid.

This guide defines the minimum evidence and verification needed for new
materials, literature EOS records, refits, and reusable equation families.
The [`.eosmat` schema reference](eosmat-schema.md) remains the authority for
individual fields.

## First classify the contribution

| Contribution | What changes | Main review burden |
|---|---|---|
| New material using an existing model | A new `.eosmat` document, primary data when available, and record tests | Identity, structure, provenance, units, and reproduction |
| New EOS for an existing material | A new record in an existing `.eosmat` document | Keeping distinct source parameterizations and reference states separate |
| Refit of a published EOS | A separate opt-in record with `record_kind: refit` | Reproducible data selection, objective, weights, fixed parameters, and statistics |
| New equation family | Native and Python implementations, interchange support where appropriate, tests, and equation documentation | Mathematical definition, domain, thermodynamic identities, numerical behavior, and API compatibility |

A material means a specific composition and phase. An EOS record means one
published, derived, or refitted parameterization for that material. An
equation family is the reusable mathematical model. Do not add a new equation
class when a publication can be represented exactly by an existing model and
an explicit configuration choice.

## Acceptance requirements

A bundled, executable EOS record must satisfy all of the following:

1. **Primary authority.** The publication of record or its official
   supplement defines the equation and parameters. A review article, database,
   another software package, or a migrated card may help locate the source but
   is not sufficient authority by itself.
2. **Unambiguous identity.** Composition, phase, structure, cell convention,
   and formula-unit count are explicit. Similar phases and different sample
   compositions remain separate materials or clearly separate records.
3. **Exact equation mapping.** The implemented expression, reference state,
   parameter definitions, fixed choices, and units match the source. A
   familiar label such as “BM3” or “MGD” is not enough when the paper uses a
   nonstandard convention.
4. **Traceable parameters.** Every stored value has a table, equation, page,
   supplement, or documented derivation. Reported errors, their confidence
   convention, fixed parameters, and covariance information are preserved;
   missing uncertainty is stored as missing, not zero.
5. **Numerical reproduction.** At least one independently specified source
   result is reproduced within a tolerance justified by reported uncertainty,
   printed precision, or digitization error. Reference-state identities alone
   are necessary tests but are not independent literature reproduction.
6. **Primary observations when available.** All recoverable rows relevant to
   the represented fit are bundled or checksummed and linked to the record,
   including reported coordinate uncertainties and calibrant observations.
7. **Honest scope.** Experimental pressure and temperature coverage, phase
   limitations, pressure calibration, data exclusions, and unresolved source
   contradictions are recorded. A marginal range is not presented as a
   rectangular phase-stability guarantee.
8. **Executable and tested.** The material loads through the public API, the
   record evaluates and inverts over its supported domain, and the complete
   project checks pass.

If any required equation coefficient, volume convention, or reference state
cannot be resolved, keep the candidate out of the executable bundle until it
can be resolved. Do not fill gaps from a different composition, phase,
pressure scale, or nominally similar parameterization.

## 1. Assemble the source record

Use the final publication and official supplementary files where possible.
Record:

- full citation, DOI or stable URL, and the exact source version;
- equation numbers and the definition of every symbol;
- parameter table and whether each value was fitted, fixed, adopted, or
  derived;
- reference temperature, reference pressure, and volume basis;
- pressure and temperature ranges actually represented by the experiments or
  fit;
- pressure calibration and the exact calibration version;
- reported uncertainty type, confidence level, weighting, covariance, and fit
  statistic; and
- inconsistencies between the abstract, prose, equations, tables,
  supplements, or source versions.

Prefer the authors' final parameterization. If two published alternatives are
scientifically meaningful, preserve them as separately named records rather
than selecting one silently. A correction or interpretation must be visible
in provenance metadata and explained in the reproduction note.

## 2. Establish material and volume identity

Before entering EOS coefficients, verify what one value of `V0` represents.
The canonical exchange volume is conventional-unit-cell volume in
angstrom cubed; molar, formula-unit, atomic, primitive-cell, density, and
lattice-parameter inputs must be converted explicitly.

A new bundled material must be diffraction-ready for PhaseSmith and
Dioptas. Include:

- a stable lower-snake-case identifier, name, formula, and phase;
- crystal system, conventional reference cell, space group, and
  `formula_units_per_cell`;
- occupied atomic sites with coordinates, multiplicities or Wyckoff labels,
  and occupancies; or a documented fallback peak list when no defensible
  atomic model is published; and
- separate provenance for the structure and for the EOS.

Check that site multiplicities and occupancies reproduce the stated cell
contents and that the lattice gives the stated cell volume. The structural
cell and an EOS record's reference cell may legitimately differ; document the
reason rather than forcing them to agree.

Extend an existing material when the composition and phase identity match.
Create a new material for a distinct phase, end-member, substitution level,
spin state, or other scientifically meaningful identity.

## 3. Capture primary observations

When a primary table or official machine-readable dataset is available, its
inclusion is required for a new bundled record. Preserve the reported data,
not a cleaned or pressure-recalculated replacement:

- retain every relevant row and identify the rows used by the published fit;
- retain original quantities and units when this prevents loss of provenance,
  while declaring their meaning in typed columns;
- include all reported pressure, volume, temperature, density, and lattice
  uncertainties and say whether they are standard deviations, standard
  errors, or interval bounds;
- retain paired pressure-marker observations needed to reproduce the pressure
  scale; and
- record the table, worksheet, figure, or supplement location and the data
  license when known.

Small tables may be embedded in `datasets[].rows`. Larger tables belong in
`peritheos/data/datasets/` and use a relative resource path, SHA-256 checksum,
and declared columns. A transcription test should check row count, columns,
representative rows, and the source extrema or grouping structure.

Digitized plots must be labeled as digitized or `plot_only`, include the
digitization method and uncertainty, and must not be described as the authors'
original fit table.

### When observations are unavailable

Absence of a recoverable table does not automatically exclude a published
parameterization. The contribution may still be accepted when the source
fully defines the equation, parameters, units, and reference state and a
published numerical value or curve provides an independent benchmark. In
that case:

- document which article, supplement, repository, and author-provided files
  were checked;
- state why a direct refit is impossible;
- classify the reproduction as parameterization-only or not directly
  refittable in the validation ledger; and
- do not synthesize observations or infer parameter uncertainties.

A parameter list with no auditable equation or volume convention is not
sufficient.

## 4. Reproduce the publication

Write a small, deterministic reproduction before treating the record as
validated. It should be independent of the test fixture: expected values must
come from the primary source, a transparent hand calculation, or a separately
documented reference implementation, not from Peritheos itself.

The reproduction should cover, as applicable:

- a reported state away from `V0` and, for a thermal model, away from `Tr`;
- the source reference-state convention and zero thermal increment at `Tr`;
- a state near the high-pressure or high-temperature end of the represented
  data;
- the stated equation variant, fixed coefficients, atom count, Gruneisen or
  characteristic-temperature law, and any staged fitting protocol; and
- residuals or fit statistics against the primary observations.

When the data used for the source fit are available, perform an independent
refit using the published row selection, objective, weights, fixed parameters,
and fitting order. Report coefficient differences and curve or residual
statistics even when the refit does not recover the published coefficients.
Failure to achieve parameter parity is evidence to investigate and disclose;
it is not permission to replace the publication silently.

A new Peritheos refit is always a separate record. Use `record_kind: refit`, a
`_refit` identifier suffix, `derived_from_record`, and complete
`fit_provenance`. The source-reported record remains unchanged and remains the
default unless there is a documented reason otherwise.

## 5. Create or update the `.eosmat` document

Follow the [schema reference](eosmat-schema.md) and an existing reviewed
material with the same model family. In addition to the equation parameters,
record:

- a stable record identifier and source/model label;
- structured citation and field-level parameter provenance;
- primary and thermal parameter errors and fixed-parameter lists;
- `temperature_ref`, experimental ranges, and range provenance;
- pressure-calibration methods and whether observation-level recalculation is
  possible;
- linked `fit_datasets` and `datasets[].used_by_eos_records`;
- `scientific_validation`, its audit date, source locations, verified fields,
  primary-data result, and unresolved issues; and
- notes for assumptions that affect interpretation or extrapolation.

Only set `primary_source_validated` after the source and numerical checks are
complete. Schema-valid JSON alone does not establish scientific validity.
Update the material manifest, primary-source audit ledger, refit ledger, and
documented catalog totals when the addition changes them.

## 6. Add a new equation family only when needed

Before implementation, write down the canonical equation, parameter domain,
units, reference state, supported branches, limiting cases, and whether its
thermal or energy terms are absolute or relative to `Tr`. Resolve notation
differences in prose instead of encoding publication-specific names into a
generic class.

Peritheos's built-in evaluator is native Rust with a public Python facade. A
new built-in family normally requires:

- the model and validation in `crates/peritheos/src/isothermal.rs`,
  `thermal.rs`, or `hugoniot.rs`;
- native batch, fitting, and uncertainty support where those operations are
  meaningful;
- PyO3 construction and dispatch plus the public Python class and exports;
- `.eosmat` type/model mapping, parsing, serialization, and JSON Schema
  changes when the model is exchangeable; and
- equation-reference, model-list, units, domain, and API documentation.

Do not expose a partially supported built-in silently through a Python-only
path. If a capability such as caloric properties, fitting, or serialization is
not mathematically defined, reject or omit it explicitly and document that
boundary.

## 7. Test the contribution

Tests should be proportional to the model but cover the following applicable
classes.

| Area | Required checks |
|---|---|
| Literature | Independent source value(s), justified tolerance, and primary-data/refit result |
| Reference state | `P(V0) = 0`, `K(V0) = K0`, and zero thermal increment at `Tr` |
| Derivatives | Analytic bulk modulus against `-V dP/dV`; stated `K0'`, `K0''`, Gruneisen, or caloric identities |
| Inversion | Pressure-volume and thermal round trips across the supported branch |
| Inputs | Invalid parameters, nonpositive volumes or temperatures, nonfinite values, and model-domain boundaries |
| Shapes | Scalars, arrays, broadcasting, stable return shapes, and ordered native batches |
| Interchange | Schema validation, Rust and Python loading, execution, serialization, and metadata-preserving round trip |
| Compatibility | Native/Python parity and existing public behavior unless a documented migration is intended |
| Data | Resource checksum, column count, row count, and representative source rows |
| Validity | Accepted in-range states and rejected out-of-range states when checking is enabled |

For a material-only addition, targeted `.eosmat`, transcription, literature,
and validity tests are usually sufficient. A new equation family needs the
full model matrix in both Rust and Python.

Run the normal contributor checks before review:

```bash
uv run ruff check .
uv run ruff format --check .
uv run --python 3.9 mypy
uv run pytest -q -W error --cov --cov-report=term-missing
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features --locked
uv run --group docs mkdocs build --strict
```

Also run `uv run python scripts/validate_primary_eos_refits.py --check` when
the change affects a bundled record or its primary dataset.

## Pull-request checklist

- [ ] The contribution is correctly classified as a material, record, refit,
      or equation-family change.
- [ ] Composition, phase, cell basis, `Z`, and reference state are explicit.
- [ ] Primary equation and every parameter have exact source locations.
- [ ] Fixed, fitted, adopted, and derived values are distinguished.
- [ ] Errors and covariance are preserved without inventing missing values.
- [ ] Recoverable primary observations and calibrant readings are included and
      tested, or their absence is documented.
- [ ] At least one independent literature value is reproduced with a justified
      tolerance.
- [ ] An available primary fit has been independently refitted and the outcome
      is recorded without overwriting the published record.
- [ ] Pressure calibration, experimental coverage, phase boundaries, and
      source contradictions are explicit.
- [ ] The `.eosmat` record validates, executes, inverts, and round-trips in
      Python and Rust.
- [ ] A new equation family includes native implementation, Python API,
      interchange support where applicable, equations, domains, and full
      numerical tests.
- [ ] Ledgers, manifests, catalog totals, references, user documentation, and
      `CHANGELOG.md` are updated.
- [ ] All relevant contributor checks pass.
