# EOSMAT v4

`BM2Reference` requires V0, K0 and P0, with V0/K0 defined at P0; ordinary BM2
remains referenced to zero pressure. Native v4 and explicit v3 exports retain
this model identifier and its parameters. A reader without the model must
reject it. See the [reference-state convention](equation-reference.md#finite-pressure-bm2-reference)
and [complete standalone example](examples/bm2-reference.eosmat).

EOSMAT v4 is a simpler **native material format**. The `.eosmat` file itself is
the complete editable record. Readers load it directly, and saving it preserves
v4. There is no separate short authoring file or generated companion.

Source verification and coefficient parity are separate from recommendation for
use. A nondefault record marked `"catalog_access": "explicit_selection"` is
omitted by ordinary `Material.from_eosmat(document)` loading, regardless of its
parity or source-verification status. An explicit `record_identifiers` selection
retains access for diagnostic work. The record, its warning and published
coefficients remain in the source document and round-trip exports.

## Small material files

This complete file contains fictional teaching coefficients:

```json
{
  "format_version": 4,
  "identifier": "my_material",
  "formula": "X",
  "eos_records": [
    {
      "identifier": "my_material_bm3",
      "model": "BM3",
      "temperature_ref": 300,
      "parameters": {"V0": 10, "K0": 160, "K0_prime": 4},
      "record_kind": "diagnostic"
    }
  ]
}
```

The model is declared once. There is no `eos` wrapper, redundant `type`/`model`
pair, separate parameter-error map, empty fixed-parameter list or required audit
boilerplate. Display names and labels are optional and fall back to identifiers.
Material and record identifiers remain explicit and stable.

**Version 4 defines exchange units:** GPa, kelvin and Å³ per conventional unit
cell. Omitted `units` means these exact conventions, not arbitrary or inferred
units. An explicit `units` object can repeat the same conventions. Source-table
columns retain their original units independently. Molar/atomic volume cannot
be entered as conventional-cell volume without a documented conversion.

## Reading order in saved files

The formatter keeps each record's model, reference temperature and coefficients
together. Thermal components consistently start with `model` and `parameters`,
followed by their explicit reference convention and model-law choices.

The record then shows the reported error confidence, experimental pressure and
temperature coverage, model validity, origin and selection restrictions. Sample,
branch and volume conventions and source-specific qualifications appear before
the longer covariance, citation, pressure-calibration and `evidence` sections.
Unrecognized record extensions also stay ahead of those detailed sections.
Material-level recommendation holds appear before `eos_records`.

This is a presentation convention, not a required JSON key order. Readers accept
any key order. Formatting moves whole fields within their existing objects;
it does not rename fields, change their nesting, omit values or infer scientific
information. Parameter
order, covariance parameter names/matrix rows, source arrays, record order and
unknown nested metadata remain intact. Historical v2/v3 formatting is retained.
Use `python -m peritheos format FILE` to apply this layout to an existing file.

## Parameters and uncertainties together

A parameter without reported error or fixed status is a number. Annotate it
when more information exists:

```json
"parameters": {
  "V0": {"value": 10, "error": null},
  "K0": {"value": 160, "error": 2},
  "K0_prime": {"value": 4, "fixed": true}
}
```

Missing or `null` error means unknown, never zero. Zero remains an explicitly
supplied zero. `fixed` describes fitting constraints independently of reported
errors. Annotation-only entries can retain uncertainty or fixed-state metadata
for implicit/adopted parameters that were not stored as equation coefficients.
The confidence convention, covariance matrix and ordered covariance parameter
names remain explicit record metadata under their established field names.

## Thermal and advanced models

`thermal` uses the same single model name and annotated parameters. All existing
law choices remain explicit, including `debye_temperature_law`,
`thermal_pressure_reference` and the expansivity/reference-volume laws.
Unknown model names or unsupported laws are rejected, never substituted.

Use the existing exact type names as native model identifiers, such as `BM3`,
`Vinet`, `MieGruneisenDebye` or `VariableQMieGruneisenDebye`. In particular,
`Sokolova2016` and `MultiOscillatorGruneisen` retain their distinct source type
identities even though the historical runtime identifier is shared.

Supply `formula_units_per_cell` when molar-volume conversion is required.
Reference temperatures follow the existing model conventions; new standalone
files should supply `temperature_ref` explicitly, as the creator requires.
A calculation reference temperature is not evidence of a measured temperature.
Hugoniots retain their complete initial state, mass basis, loading path, branch
domain and validity restrictions. Primary-component configuration/extensions
that previously lived inside `eos` are retained in `equation_options`.

## Complete database evidence

Detailed scientific evidence is optional for a standalone file and retained
for database records in one `evidence` object:

```json
"evidence": {
  "scientific_validation": {"status": "pending_primary_source_check"},
  "parameter_provenance": {"K0": "Source table and fit selection"}
}
```

This object also holds `source_lineage`, `fit_provenance`, `fit_protocol`,
`reproduction` and `audit_corrections` when present. Absence of scientific review
means **unreviewed**; creating or converting a file cannot promote it to
primary-source-validated. A refit still requires its complete fit provenance.
Scientific review, numerical reproduction and recommendation for use remain
separate claims.

References, pressure calibration, experimental coverage, operational validity,
default selection, recommendation holds and source/reference-variant links
retain their meanings. Structure and original `datasets` stay attached to the
material. Rows, column definitions, errors, hashes and explicit observation
views are retained. No observation is refitted or recalibrated during conversion.
Other material/record extension fields survive round trips. A collision with a
new reserved field is reported rather than silently discarding an extension.

## Optional shared metadata, inside the same file

Small files should stay inline. Large files may declare `shared_metadata`, an
object mapping readable lower-snake-case names to metadata objects:

```json
"shared_metadata": {
  "source_example": {
    "authors": ["Teaching example"],
    "year": 2026,
    "source": "Fictional teaching source"
  }
}
```

Inside a record's `reference`, `evidence` or `pressure_calibration`, a singleton
`{"shared_ref": "source_example"}` substitutes that exact object. References can
be nested inside those metadata locations. They do not apply to models,
parameters, structure, datasets or unrelated extensions. No remote reference,
inheritance, merge or override is supported. Dangling names, cycles, sibling
fields on a reference and invalid definitions are errors, including in unused
definitions. Expansion is bounded to 32 reference levels and 100,000 metadata
nodes when a sharing table is present.

`convert_eosmat(document, compact=True)` shares only exactly repeated objects
large enough to benefit from sharing; it never merges merely similar sources.
The result is still a complete native EOSMAT file. Record order, source arrays,
coefficients, unknown errors and fixed flags are unchanged. Normal load/save
preserves the authored document; converting a `Material` back to a document can
normalize model metadata, as for inline files.

`convert_eosmat(document, inline=True)` expands local references. It retains all
shared definitions, including unused supporting evidence, as inert metadata.
V3 export also expands them. A first-generation v4 reader needs an inline export;
local sharing is an opt-in capability, not a requirement for old or new files.
The updated JSON Schema describes the reference shapes; use the Python or Rust
semantic validator as well to check resolution, cycles and expanded model data.

## A consistent reconstruction assessment

See [EOS labels and reconstruction status](eos-labels.md) for the complete
label definitions, default comparison criteria and scope qualifications.

All bundled EOS records now carry a standard `assessment` first inside
`evidence.reproduction`. Its outcome, scope, comparison basis and limitations
are copied from the existing numerical reconstruction ledger and record evidence;
no new fit or scientific classification is performed. Historical reproduction
fields remain intact alongside it. Standalone user files may omit assessments:

```json
"evidence": {
  "reproduction": {
    "assessment": {
      "status": "not_assessed",
      "scope": "parameter_agreement",
      "qualifications": ["No independent refit has been performed."]
    }
  }
}
```

`status` is one of `parity`, `similar`, `not_reproduced`, `not_refittable`,
`source_reconstruction` or `not_assessed`. `scope` is `parameter_agreement`,
`author_workflow`, `independent_refit` or `source_equation`. An assessed result
also requires a nonempty `basis` describing the exact comparison criterion.
Optional `qualifications` is an array of nonempty strings and `report` is a
nonempty report identifier/path. Unknown additional assessment fields are retained.
Python, Rust and JSON Schema validate this standard block. Historical free-form
reproduction findings are not rewritten or used to invent an assessment.

`source_status` preserves the original ledger code. In particular,
`parity_not_achieved` maps to native `not_reproduced` **for the stated comparison
scope**, without changing the record's separate scientific-validation status.
Stored independent Peritheos refits use `scope: independent_refit`; reproducing
their coefficients does not reproduce the original publication. Source-calculation
reconstructions use `source_equation`. Parameter agreement never automatically
becomes an `author_workflow` claim.

### Actual reconstructed coefficients

When available, `assessment.parameters` stores the original coefficient
comparison rows verbatim, including unsuccessful and partial reconstructions:

```json
"parameters": [
  {
    "parameter": "K0",
    "published": 215.3,
    "published_error": null,
    "refit": 254.92768974032103,
    "refit_error": null,
    "difference": 39.62768974032102,
    "relative_difference": 0.18405801087004653,
    "similar": false,
    "within_combined_2sigma": null
  }
]
```

This row is from the partial 298 K reconstruction of
`akimotoite_dorogokupets_2015_298k_rydberg_stacey`; it is not a replacement EOS.
The main record's `parameters` remain the published source parameterization.
`published` denotes the comparison target; for an independent-refit assessment
it denotes the stored Peritheos refit. `refit` is the actual reconstructed value.
Error fields retain the audit's convention, including historical `*_sigma`
spellings; missing or null is unknown, never zero. Ordered comparisons and all
additional flags remain intact. Parameter names, component prefixes and units
follow the named audit, without an inferred conversion.

`assessment.fit` records the available objective, selected observations, weights,
free/fixed parameters, covariance and fit statistics. Fixed/adopted coefficients
are listed separately and must not be called independently reconstructed.
Missing tables or fit details remain missing. `report` identifies the full ledger;
look up its `record_identifier` using the EOS record's identifier. Dedicated audit
links and detailed historical evidence remain available in the same record.

### Open the human-readable reconstruction report

Bundled assessments also contain two direct destinations:

```json
"report_url": "https://peritheos.readthedocs.io/en/latest/reconstruction-reports/akimotoite/#record-akimotoite_dorogokupets_2015_298k_rydberg_stacey",
"report_path": "docs/reconstruction-reports/akimotoite.md#record-akimotoite_dorogokupets_2015_298k_rydberg_stacey"
```

`report_url` opens the specific EOS section on **Read the Docs latest**.
`report_path` is relative to the repository and identifies the same section in
offline or version-specific documentation builds. Material summaries use relative
links, so they remain inside the documentation version being read. The hosted
URL becomes available only after the corresponding docs are published. `latest`
can change after subsequent audits; use the repository path at the relevant
commit or documentation version when an exact historical assessment is needed.
These optional fields do not require standalone user files to have online reports.
The existing `report` field still identifies the full machine-readable ledger.

The [reconstruction reports](reconstruction-reports/index.md) cover every main-catalog
EOS. They explain scope, comparison basis, recorded data/method, fixed constraints
and limitations. Where numeric comparison pairs exist, a graph shows source/target
and reconstructed **coefficient central values**, with full values and reported
errors in a table. These graphs are not observation/residual plots, uncertainty
bands or new fits. Missing values are never inferred. Supporting source reports
and their existing observation/residual plots remain linked where explicitly
identified by the evidence. Illustrated discussions that explicitly name an EOS
are also linked as related reports; mentioning an EOS does not establish that
the discussion reproduces its fit.

Repository maintainers generate these pages and graphs with
`python -m scripts.generate_reconstruction_reports`; `--check` detects drift.
It reads stored assessments without rerunning scientific fits. Plot checks use
input and SVG-content fingerprints, avoiding byte differences caused solely by
different Matplotlib versions. Plot generation needs the documentation dependencies.
After a reviewed audit change, update the assessment and regenerate reports and
material summaries together.

Repository maintainers can run
`python -m scripts.standardize_reconstruction_metadata --check` to detect drift.
Without `--check`, the tool adds missing assessments from the existing ledger.
It refuses to replace a differing assessment: review existing findings explicitly.
Update the ledger and affected assessment together when a scientific audit changes.
This repository tool does not add assessments to arbitrary user files or infer
parity from successful parsing, source verification or recommended/default status.

A numerical parameter match is not automatically an author-workflow reproduction,
a reproduced independent refit does not establish source-author parity, and a
source-equation checkpoint is not an independent data fit. Neither creation,
conversion nor successful parsing sets a parity status or recommends an EOS.

## Migration and compatibility

Python and Rust read v2/v3/v4, including native local-sharing files with updated
readers. Standalone files and all 232 bundled library files use native v4. The
library keeps units visible in the header and evidence inline; migration changes
the organization, not scientific values, observations, source evidence or usage
restrictions. Older
applications whose pinned reader supports only v3 need an explicit v3 export:

```bash
python -m peritheos convert old.eosmat new.eosmat
python -m peritheos convert new.eosmat compatible.eosmat --format-version 3
```

The destination must be new. These commands migrate complete files, not a
separate authoring representation. No second file is needed during normal
editing. Original versions can be kept as backups.

```python
from peritheos import convert_eosmat, eosmat_schema, load_eosmat, save_eosmat

native = convert_eosmat(load_eosmat("old.eosmat"), version=4)
save_eosmat("new.eosmat", native)
legacy = convert_eosmat(native, version=3)
schema = eosmat_schema(4)
```

`Material.from_eosmat()` supports v4 directly. `material.to_eosmat()` retains the
loaded format; use `format_version=3` or `4` for an explicit target. `from_dict()`
also recognizes native v4. Python's historical `eosmat_schema()` default remains
v3 for callers validating existing archives; request `eosmat_schema(4)` for
native files. Rust provides `export_eosmat_v3()` alongside direct v4
load/validate/save. Unsupported models, conflicting layouts and malformed
parameter annotations fail validation.

Conversion preserves scientific values, fixed-state/error semantics, source
arrays and unknown extensions. Optional empty thermal error/fixed containers
can be omitted because their meaning is unchanged. V4 metadata remains
unreviewed unless the original record already contains explicit scientific
review. See [legacy v3 field definitions](eosmat-schema.md) for unchanged
advanced semantics, and [standalone authoring](eosmat-authoring.md) to create a
new file.
