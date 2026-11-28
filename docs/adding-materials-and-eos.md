# Adding a material or equation of state

An addition to the Peritheos material library is a scientific claim, not only
a data-entry change. A reviewer must be able to identify the represented
material, trace every equation and parameter to its source, execute the same
formulation in Peritheos, and understand where the result is valid.

New and bundled library files use [native EOSMAT v5](eosmat-v5.md). The scientific requirements
below apply equally to v5 and historical v3/v4: native EOSMAT keeps parameter values/errors/fixed
flags together. A curated record carries a short `check` headline; its detailed
audit trail lives in a per-material sidecar (see [step 5](#5-create-or-update-the-eosmat-document)).
References to legacy field paths below retain the same scientific meaning.

This guide defines the minimum evidence and verification needed for new
materials, literature EOS records, refits, and reusable equation families.
The [EOSMAT v5 record reference](eosmat-reference.md) is the authority for
individual fields.

Start authoring with the small [EOSMAT template and editing guide](eosmat-authoring.md).
Use the generated [material summaries](material-summaries/index.md) to inspect
existing parameterizations and their qualifications. The template is fictional
and deliberately unaudited; the acceptance requirements below still apply.

### Reconstruction report presentation

Every [reconstruction report](reconstruction-reports/index.md) must expose the
original/source and reconstructed coefficients, available parameter errors, and
fit error for both curves. Show P–V curves (pressure on the horizontal axis),
or the appropriate fitted quantity for velocity, acoustic and joint objectives,
the selected input rows, and residuals for each curve on the same rows. Define
the residual sign, coordinate, units and RMS statistic. Distinguish parameter
uncertainty from descriptive residual error; retain confidence, weighting,
fixed-parameter and staged-fit qualifications. A stored Peritheos refit used as
the comparison target must be labelled as such.

Round displayed uncertainties to two significant digits when their leading digit
is 1 or 2, otherwise one, and round the corresponding value to the same decimal
place. For example, `25.447259 ± 0.091348` is displayed as `25.45 ± 0.09`, and
`24.9234 ± 0.1125` as `24.92 ± 0.11`. Preserve the selected decimal place on
carry (`0.099` becomes `0.10`). Values without a positive reported error use at
most six significant digits; this display limit does not estimate uncertainty.
Missing errors must be labelled; a reported zero does not mean an exact value.
Preserve full source and computed precision in EOSMAT and numerical evidence. A reconstruction
comparison may explicitly record `display_decimal_places` when an audited
printed-precision comparison requires the same decimal place on both source and
refit values/errors. This affects report display only: retain unrounded values,
error confidence and any uncertainty mismatch. The usual rounding convention
continues to apply when that optional field is absent.
The repository audit ledger may keep explicit `assessment_basis` and
`assessment_qualifications` text when a finding is more specific than the
default wording.

The repository-only maintenance command
`uv run --group docs python -m scripts.generate_reconstruction_reports` replays
stored coefficients without optimizing or changing scientific assessments.
`--check` checks generated pages and fingerprinted SVGs. Its P–V adapters reuse
the named audit's row selection, pressure reduction and volume normalization;
row-count and available source-RMS mismatches prevent an incorrect comparison.
Pressure RMS is evaluated at observed coordinates and may differ from an
errors-in-variables objective evaluated at fitted latent coordinates. Thermal
curves are separate labelled isotherms; residuals use every row's temperature.
Do not connect mixed-temperature observations as one fitted curve.

Curve SVGs retain an integrity-checked, full-precision input snapshot. Freshness
checks allow only cross-platform numerical roundoff (relative and absolute
tolerance `1e-9` in the displayed coordinates); labels, selection sizes, renderer
code and SVG body integrity must still match. This cache tolerance does not
round stored observations, coefficients, model evaluations or reported errors.

Specialized adapters replay the audit's equations in their actual coordinates:
Us–up, density–velocity, bulk/shear modulus, energy-balanced shock paths and
coupled mixed-spin trajectories. A reference HS/LS endmember must not be plotted
as though it were fitted directly to the mixed-spin observations. Co-compression
reports show the joint phase-pair equality residuals; their separate calculated
material-curve differences are explicitly distinguished from data residuals.
RMS values for different outputs retain their native units and separate panels.
When a statistic uses a different selection (for example, excluding an ambient
anchor), identify that difference beside the plot.

Available-subset plots must be prominently marked **partial replay**, both on
the report and in the index, with the missing inputs named there. Do not run a
new regression or silently substitute a different forward calculation to fill
such a gap.
These plots do not change the recorded reconstruction findings.

Generated curve/residual and coefficient plots have paired light and dark SVGs.
The documentation's theme toggle selects the appropriate image using Material's
`#only-light` / `#only-dark` image fragments, including when the site theme differs
from the operating-system preference. Dark variants change paint only: the
scientific geometry, text and precision are identical. Direct links to the
unsuffixed SVG remain light; `-dark.svg` is the standalone dark version.

Where the selected observations, coefficient mapping or pressure-reduction
replay are unavailable, state the specific gap beside the recorded error
statistics and link the detailed source audit. Do not invent observations,
uncertainties or a replacement model to fill a missing panel. New specialized
audits should supply a replayable row reader in
`scripts/reconstruction_plot_data.py`, with tests against their stored numerical
diagnostics, so both curves and residuals appear inline.

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

When a study publishes fits of the same observations on several pressure scales,
retain one source per equation family and data/constraint selection. Preserve
the other published fits with `pressure_scale_variant_of` pointing to that
source, rather than presenting them as independent EOS choices. Select output
pressure scales through an audited calibration route; the reference-variant
link itself does not authorize conversion. Different samples, temperatures,
pressure media, data subsets, and fitting constraints remain separate choices.
Document the source selection and unresolved calibration variants.

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
   missing uncertainty is stored as missing, not zero. When the source states
   the convention, set `parameter_error_confidence` (`0.6826894921370859` for
   1 sigma, `0.95` for 95%) and keep the printed errors unchanged.
5. **Numerical reproduction.** At least one independently specified source
   result is reproduced within a tolerance justified by reported uncertainty,
   printed precision, or digitization error. Reference-state identities alone
   are necessary tests but are not independent literature reproduction. A
   passing checkpoint does not override an unresolved discrepancy elsewhere
   within the record's claimed scope; apply the stopping rule below.
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
- record-level `determination_method` (`experimental`, `theoretical`, `hybrid`,
  or explicitly `unknown`) and its rationale in
  `parameter_provenance.determination_method`, considering fixed coefficients
  and thermal terms as well as the fitted data;
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

Do not treat an unspecified source-data license as automatic grounds for
omitting recoverable numerical observations, and do not label that source CC0.
Preserve explicit source terms and attribution. For factual transcriptions,
CC0 may cover only rights held by Peritheos contributors in transcription,
normalization, metadata and arrangement; state that scope separately and keep
unspecified source terms marked as unspecified. Article text, figures and
publisher files are not covered by this dedication.

The [licensing guide](licensing.md) explains redistribution scope. Add an
attributed dataset notice for explicit source terms and include it in the
installed licensing guide's index. Keep original author code under its own
terms. Inclusion for scientific reproducibility does not establish compliance
with a downstream distributor's licensing policy.

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
  refittable in the validation ledger (headline `untested`, not a failure); and
- do not synthesize observations or infer parameter uncertainties.

For a concrete unresolved data gap, add a reviewed entry to
[incomplete datasets and author requests](dataset-requests.md). State what is
missing, what it would enable, whether author outreach is worthwhile, and the
actual contact status. Link the audit and distinguish minimum inputs from
optional calibration or uncertainty information. Keep the request open until
received data have been validated and the scientific ledgers updated.

A parameter list with no auditable equation or volume convention is not
sufficient.

## 4. Reproduce the publication

An explicitly finite-pressure BM2 source can use `BM2Reference(V0, K0, P0)`
after establishing its equation mapping against independent source evidence.
V0/K0 then mean values at P0. Preserve original notation and values rather
than synthesizing ambient coefficients. A published-curve check establishes
equation mapping at graphical resolution; observation residuals and
author-workflow reproduction remain separate. The
[joint coesite-IV/V audit](literature-reproductions/coesite-iv-v-bykova-2018-experimental.md)
shows a reproduced source curve with substantial observation scatter.

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
it is not permission to develop a replacement model as part of curation.

### Stopping rule: unresolved reproduction discrepancies

If the published EOS or its reported agreement with the source data cannot be
reproduced within justified tolerances after checking transcription, equation
mapping, units, reference state, and the reported fitting protocol, document
the discrepancy and **do not recommend the affected EOS**. A documented
non-parity outcome completes the curation investigation; producing a better
fit is not required.

Preserve the published coefficients and observations unchanged. Record the
literal source calculation, quantitative differences, comparison tolerances,
affected pressure/temperature range, checks performed, and unresolved inputs
or contradictions in the audit sidecar. Set the headline to `disagrees` and the
source-validation status to `not_reproduced`; use `deferred` when the
source/model mapping is incomplete. Do not propose such a record as a default.
A maintainer may accept a scoped discrepancy for a specific record; the
acceptance and its reason belong in that record's audit.
Any retained executable diagnostic must be explicitly identified as such and
not recommended for scientific predictions or pressure calibration. Keep
catalog descriptions and scientific ledgers consistent with that finding.

Source-supported diagnostic checks may clarify the discrepancy. Stop before
trying successive alternative EOS families, thermal laws, extra parameters,
or new fitting protocols to create a replacement. Existing exploratory refits
remain diagnostic evidence; smaller residuals or reproduction of their own
coefficients do not establish published-EOS parity or a validated replacement.
Replacement-model development requires an explicit user request as a separate
research task and remains distinct from curation of the original publication.

Scope the finding to the affected parameterization and domain. Individually
verified isotherms may be retained without implying validation of an unresolved
continuous thermal EOS. Missing author-workflow details alone do not establish
a numerical reproduction failure when the published parameterization passes
independent source benchmarks. For unavailable observations, follow the
parameterization-only validation route above.

A new Peritheos refit is always a separate record. Use `record_kind: refit`, a
`_refit` identifier suffix, and complete `fit_provenance`. Add
`derived_from_record` when another record supplies model choices or fixed
parameters; a direct dataset refit does not invent a catalog parent. Any
source-reported record remains unchanged and remains the default unless there
is a documented reason otherwise.

## 5. Create or update the `.eosmat` document

Follow the [EOSMAT v5 record reference](eosmat-reference.md), which defines the
exact fields, and an existing reviewed material with the same model family.

The **record** holds what a reader needs to use the EOS: identifier and label,
citation, coefficients with errors and fixed flags, units, `temperature_ref`,
experimental ranges, linked datasets, and a `check`
block. `check.reproduction`, `compared`, `why` and `source_check` follow
mechanically from the audit's outcome code, scope, coefficient table and
source-validation status ([fields](eosmat-reference.md#the-check-block); words in
[EOS labels](eos-labels.md)); the split command below derives them.
Set `check.warning` only for an explicit usage restriction, such as a
reproduced fit that describes its data poorly; otherwise leave it null. It is
generated from a coded `scientific_validation.usage_recommendation`
([codes](eosmat-reference.md#giving-a-bundled-record-a-warning)); on a split
material set the code with `scripts/apply_check_text.py`, which keeps any
replaced free text as `recommendation_reason`.
The material's `default_records` names its default records. Propose default
changes to the maintainer, following the
[selection guideline](eos-labels.md#warning).

The **audit sidecar** `peritheos/data/audits/<id>.audit.json` holds the trail:
detailed outcome code, scope, basis, conditions, fixed parameters, refit
coefficients and statistics, pressure calibration, source-check locations,
verified fields and unresolved issues. Standalone user files may keep evidence inline and need
neither a sidecar nor a `check` block.

A **new** bundled material may be written as one complete inline-evidence v5
file in `peritheos/data/materials/`, without `check` blocks or a sidecar. The
repository maintenance command `python scripts/split_eosmat_audits.py --add`
(optionally naming materials) splits every bundled file that still holds
inline audit blocks or has no sidecar, derives each `check` block
mechanically from the stored evidence, removes only the generic statements
listed in the receipt `docs/data/eosmat-v5-audit-split.json`, and adds a
binding there that recovers your inline file byte for byte. It refuses a file
that is not canonically formatted or whose split would change anything else.
Run `python scripts/split_eosmat_audits.py --verify` afterwards.

### Changing a bundled material

`python scripts/split_eosmat_audits.py --add` refuses a material that is
already bundled: re-splitting it would replace its receipt binding, so its
original pre-split file, the stored text originals and the older receipts
that depend on them (v5 migration, rights review, scientific snapshots) could
no longer be verified. Use the repository tool for the kind of change:

| Change | Tool |
| --- | --- |
| `label`, `notes`, `validity.notes`, `thermal.notes`; the material's top-level `notes` | `python scripts/apply_record_text.py` ([record text standard](record-text-standard.md#6-applying-a-rewrite)) |
| `record_note`, `warning_reason`, coded `usage_recommendation` | `python scripts/apply_check_text.py` |
| Every generated `check.why` line, after the library changes its wording | `python scripts/apply_check_text.py --refresh-why` (only `why` may change) |
| Add an EOS record to an existing material; set or correct the fields listed below | `python scripts/apply_catalog_revision.py` |

`apply_catalog_revision.py` takes one request file per revision and records it
in the ledger `docs/data/catalog-revisions.json`: who decided it and why, the
source and location it rests on, a rights statement, explicit forward and
reverse operations and the hashes needed to verify and undo it. One request
changes one material; it may change several of its records. A field change
names each field and its current value (numbers compare numerically, so
`172` matches a stored `172.0`; setting the current value again is refused):

```json
{
  "id": "2026-11-20-qandilite-lv-2016-confidence",
  "material": "qandilite",
  "kind": "set_fields",
  "date": "2026-11-20",
  "decided_by": "<name>",
  "reason": "The source states its errors as one standard deviation.",
  "source": {"citation": "<reference>", "location": "<table, page>"},
  "rights": {"reviewed_by": "<name>", "date": "2026-11-20",
             "statement": "Sets the stated error confidence; no third-party content is added."},
  "changes": [
    {"record": "qandilite_lv_2016_bm2", "path": ["parameter_error_confidence"],
     "expected_absent": true, "value": 0.6826894921370859},
    {"record": "qandilite_lv_2016_bm3", "path": ["parameter_error_confidence"],
     "expected_absent": true, "value": 0.6826894921370859}
  ]
}
```

This change is already in the ledger as revision
`2026-11-21-qandilite-one-sigma`, so the request serves as a worked example
and would now be refused as written.

Use `"expected": <value>` when the field exists and `"remove": true` instead
of `"value"` to delete an `error`, a `fixed` flag or the covariance. Only these
fields are settable: a parameter's `error` and `fixed` flag (of an existing
reference-isotherm or thermal parameter), `parameter_error_confidence`,
`parameter_covariance`, `parameter_uncertainty_convention` (only on a record
that already carries it), `temperature_ref`, `experimental_pressure_range_gpa`,
`experimental_temperature_range_k`, `pressure_range_status`, the
`pressure_gpa`, `temperature_k` and `volume_ratio` bounds of a `validity`
object, and `volume_basis.molar_mass_g_mol` (only on a record that already
carries it). Everything else is refused, among it the identifier, all
other text (use the text tools), `check`, `model`, the thermal model and its laws,
`record_kind`, `determination_method`, `equation_kind`, `reference`,
`aliases`, default flags, the material's `default_records`, `catalog_access`
and dataset links.

A parameter stored as a bare number (`"K0": 444.5`) takes an `error` or
`fixed` flag like any other: the revision writes the object form with the same
value (`"K0": {"value": 444.5, "fixed": true}`), and removing its last field or
undoing the revision writes the bare number again. Its `value` stays
untouchable.

A `validity` object, when present, is the only enforced range (the
experimental ranges are then not enforced). A revision may add a missing bound
to an existing object (`"expected_absent": true`), for example to enforce a
stated pressure range on a record whose object holds only notes. It may also
create the object on a record that has none: the same request must then set
`pressure_gpa` if the record states `experimental_pressure_range_gpa` and
`temperature_k` if it states `experimental_temperature_range_k`, each
containing that range, so nothing silently becomes unbounded; `volume_ratio`
is optional. This is how a record gets an enforced range slightly wider than
its data (for example down to its reference temperature) while the
experimental ranges keep describing the data behind the fit. A bound cannot
be removed by a request, and an existing object cannot be removed or emptied;
undo removes what the revision added.

`volume_basis.molar_mass_g_mol` is a convenience value stating the amount of
material behind each volume, not a published coefficient: it is not
interpreted except on Hugoniot records, where V0, `rho0` and M must agree
([EOSMAT reference](eosmat-reference.md#units-and-volume)). A correction must
be a finite number above 0 and at most 10000 g/mol; there is no chemistry
check.

Values are checked before anything is staged: an error is a finite number
`>= 0` or `null`; `parameter_error_confidence` lies strictly between 0 and 1
or is `null`; `parameter_uncertainty_convention` is non-empty text; a range is
two finite numbers, lower first, with pressures
within -20 to 2000 GPa (the catalog spans -11.3 to 1374 GPa) and
temperatures above 0 K and at most 20000 K; `pressure_range_status` is one of
`reported_exactly`, `reported_qualitatively`, `theoretical` and
`reference_parameterization`; a changed `temperature_ref` is positive, equals
a thermal `Tr` and lies inside a stated temperature range; a covariance names
free (not fixed) parameters of the record after the change and is square,
finite, symmetric and positive semidefinite. Every touched record is then
loaded as the catalog loads it, evaluated at a state inside its enforced range
and its parameter uncertainty propagated.

The dry run shows each changed value, the enforced range and the `check`
block. A change of the enforced range is refused unless the request changes a
range field. A change of a headline field of `check` (`reproduction`,
`compared`, `warning`, `source_check`), which a changed error or fixed flag
can cause through the combined-error criterion, is refused unless the request
acknowledges exactly that change, for example
`"acknowledge_check_change": {"<record>": {"reproduction": {"from":
"disagrees", "to": "close"}}}`; the acknowledgement is stored in the
revision. A changed `why` line alone needs none.

Some records cannot be revised yet, and the tool says why before any other
work: records whose saved refit-ledger row differs from a fresh computation on
this machine (`refit_ledger_stale`, with the largest difference); and records
whose refit row takes its published values from a script rather than the
record (`row_disagrees`, naming the script). A record that
`docs/data/uncertainty-reconstruction-inputs.json` pins by digest can be
revised: the pin is checked on the record as it was before any catalog
revision, recovered through the ledger.

An added record uses `"kind": "add_record"` and `"record"` instead of
`"changes"`, with `"same_basis"` in `rights`. It is the complete record
without `check`, with these blocks (unknown keys are refused):

- `identifier` (unused anywhere in the catalog), `label`, `model`,
  `parameters`, `reference`, `determination_method`, `record_kind`,
  `experimental_pressure_range_gpa` and `experimental_temperature_range_k`;
- `pressure_calibration` with `status` (`resolved`, `partially_resolved`,
  `not_applicable` or `unresolved`), `methods`, `recalculation.status` and
  `audit_date`;
- `evidence.scientific_validation` with `status` (`primary_source_validated`,
  `not_reproduced` or `deferred`), `audit_date` and
  `primary_source_check.finding`, plus what the refit validator reads, such as
  `primary_data_check` with its `dataset_identifiers`;
- `evidence.parameter_provenance`, which the rights review binds.

The reproduction assessment is not written by hand: the tool computes the
record's row of the registered refit ledger with
`scripts/validate_primary_eos_refits.py`'s own function for that material,
derives `evidence.reproduction.assessment` (status and scope included) from
it and generates `check`. The material's `default_records` never changes.

The **rights statement** is required because the rights review
(`docs/data/eos-structure-rights-policy.json`) pins the whole reviewed content
of each material. It names who reviewed and when, and states what was
checked. For an added record it says whether the record is covered by the
same source and redistribution basis as the material's existing records
(`"same_basis": true`) or, with `"same_basis": false`, what differs
(`"differences"`). For a field change it says that no third-party content is
added. The tool stores the statement in the revision, moves the material's
pin in the policy, lists the revision there and regenerates the review
outputs; it never writes the statement itself.

```bash
python scripts/apply_catalog_revision.py request.json --check   # show the result, write nothing
python scripts/apply_catalog_revision.py request.json           # apply
python scripts/apply_catalog_revision.py --verify
```

Everything is computed in a scratch copy first. The tool writes the record
file, the sidecar, the audit-split receipt, the ledger, the rights policy and
review, the refit-ledger row and `docs/primary-eos-refits.md`, the
primary-source audit entry and manifest counts, catalog coverage, the material
summaries, the reconstruction reports, the paper investigation ledger and Rust
fixture copies that equal the bundled files. Before writing, it loads the
staged catalog through the library (every material, the pressure-calibration
graph and the revised records with primary validation) and runs every
recovery check; a failure writes nothing. Do not run the full
`python scripts/validate_primary_eos_refits.py` for a revision: it refits every
record and rewrites the whole ledger, and its `--check` currently fails: the
saved ledger lacks one policy sentence the script generates, and many saved
rows differ from a fresh computation in copied record text and in the last
digits of fitted numbers. Nothing needs a hand edit, except
a test that pins a material's record count.

`python scripts/apply_catalog_revision.py --undo <id>` reverses a material's
most recent revision and refuses an older one. Later text and check-text
revisions of the same material are kept (those of a removed record go with
it); without them every file returns byte for byte.

`--verify` checks each revision's fields and its content digest, which covers
every stored field and which the audit-split binding pins for a material's
latest revision; the chain of a material's revisions; the exact walk back to
the original pre-split file; the scientific snapshot; the rights pins; the
refit-ledger rows; and the manifest and source-audit counts. Generated pages a
revision wrote are listed by path only: later changes rewrite them, and their
generators' `--check` keeps them current.

The ledger catches accidental and unexplained changes: every check of a
bundled material fails when a stored value, hash, order or decision differs
from the files. It does not prevent a deliberate, fully consistent rewrite by
someone with commit access; that is visible only in the git history.

Not supported yet: correcting a published coefficient `value` (raise it with
the maintainer), removing a record, material-level fields
(`default_records`, structure, datasets), changes to the audit evidence of an
existing record other than through the two text tools, and new materials
(use `--add`). Propose default changes to the maintainer.

Write both briefly. State each finding once, not again in basis, notes and
qualifications. Do not restate a label's generic caveat; [EOS labels](eos-labels.md)
holds it. Keep solver diagnostics and full-precision refit numbers out of the
record. Do not hard-code catalog counts in prose or add sections that only say
nothing is recorded. Prefer a short specific statement over defensive
qualification. Labels, notes, validity notes and thermal notes follow the
[record text standard](record-text-standard.md); rewrite existing record text
and material notes only with `scripts/apply_record_text.py`, which keeps every
replaced original in the sidecar.

Write the source check's complete result in
`scientific_validation.primary_source_check.finding`. The record shows it as
`check.source_note`; when it is long, restates coefficients or is shared with
other records, also write a per-record
`primary_source_check.record_note` (at most 300 characters: what was checked
against the paper for this record and any conversion or correction that
matters to a user), which `check.source_note` shows instead. A record with a
`warning` may get a one-clause `scientific_validation.warning_reason` (at most
140 characters) saying why, taken from the audit. In an inline file these are
ordinary evidence fields; for an already split material use the repository
tool `scripts/apply_check_text.py` (`--check`, apply, `--verify`; format in
the [record text standard](record-text-standard.md#6-applying-a-rewrite)).

| | Record `check` | Audit entry |
|---|---|---|
| Bad | "Parity achieved. This does not establish experimental accuracy, phase stability or author-workflow reproduction …", repeated in notes, plus refit V0=263.929881998 | The same disclaimer again under basis and qualifications |
| Good | `reproduction: reproduced`, `compared: published_fit`, `warning: null` | Scope `parameter_agreement`; basis "within printed source widths, K0 fixed at 210 GPa"; refit coefficients with errors; source "Table 2, pp. 64-65" |

Only set `primary_source_validated` after the source and numerical checks are
complete. Schema-valid JSON alone does not establish scientific validity.
Update the material manifest and ledgers, and regenerate generated totals, when
the addition changes them.

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
| Reference state | `P(V0) = 0` (or explicit `P0` for `BM2Reference`), `K(V0) = K0`, and zero thermal increment at `Tr` |
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
the change affects a bundled record or its primary dataset (repository-only
tool). It recomputes every row and matches saved and fresh rows by
`record_identifier`, so the order of rows in the file does not matter, and it
lists every differing row with its first differing leaf. Numbers are compared
with a 1e-5 relative tolerance, except the source facts named in the script's
`SOURCE_FACT_KEYS` (published, printed, fixed and observed values, their
errors and ranges), which must match exactly, as must text, counts and
classifications. A few named calculated keys have their own rule, documented
next to their key lists in the script: self-consistency maxima of an
equation against an independent implementation (1e-9 absolute), covariance
off-diagonals (on the correlation scale), differences from a published value
(on the coefficient scale), and optimizer evaluation and iteration counts
(not compared when both fits report success). It currently fails, independent
of any record, on rows that need a maintainer decision: the two rank-deficient
diamond double-Debye fits, whose coefficients drift on this build; the
NaCl-B2 and neon Fei (2007) rows, whose calculated cold-stage K0 and K0'
are stored under `fixed_parameters` and so are compared exactly; and five
small drifts (1e-5 to 1e-4 relative) in weakly determined quantities (Morard
2026 FeS, two Hirose 2008 gold rows, Luo 2023 MgO, Litasov 2013 siderite).
A revision of a bundled material instead gets its rows from
`python scripts/apply_catalog_revision.py`
([changing a bundled material](#changing-a-bundled-material)).

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
- [ ] Unresolved reproduction discrepancies are documented and the record's
      headline is `disagrees`; replacement-model development has not been
      substituted for the curation outcome.
- [ ] The record's `check` headline and any warning are set; the audit trail is
      in the sidecar and no finding is repeated.
- [ ] Pressure calibration, experimental coverage, phase boundaries, and
      source contradictions are explicit.
- [ ] The `.eosmat` record validates, executes, inverts, and round-trips in
      Python and Rust.
- [ ] A new equation family includes native implementation, Python API,
      interchange support where applicable, equations, domains, and full
      numerical tests.
- [ ] Ledgers, manifests, generated totals, references, user documentation, and
      `CHANGELOG.md` are updated.
- [ ] All relevant contributor checks pass.

## Assign an optional browsing family

Set `family_id` on a material card only when its documented composition and
phase fall within a definition in `peritheos/data/material-families.json`.
Family membership is curated, not inferred from names or formulas. Preserve
exact composition, phase, structure and provenance on the individual card.

For a new family, add a unique stable `identifier`, `name`, `description` that
states inclusion/exclusion boundaries, and optionally a descriptive `formula`.
Do not add a member list: membership is derived from material cards. Leave
unrelated or uncertain materials unassigned; no single-member family is needed.
Run the catalog family and `.eosmat` tests after editing membership. Families
are discovery metadata and do not define mixing or interpolation rules.

For independently calculated reference peak arrays, keep the structural input
license/assessment separate from the generated-output dedication. Preserve the
source cell versus reference-cell distinction, occupational disorder, proxy
substitutions, scattering convention and omitted intensity corrections. See
[the ten-pattern replacement](research/calculated-reference-patterns.md) for a
repository-only reproducible calculation and private preservation receipt.
Installed core records expose `source.pattern_kind`, `source.pattern_inputs` and
`source.pattern_input_source_terms`; no observation extra is needed to read them:

```python
from peritheos import get_material_document

record = get_material_document("feh2")
assert record["source"]["pattern_kind"] == "independently_calculated_xray"
print(record["source"]["pattern_input_source_terms"])
```
