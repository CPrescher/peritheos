# EOSMAT v4

EOSMAT v4 is a simpler **native material format**. The `.eosmat` file itself is
the complete editable record. Readers load it directly, and saving it preserves
v4. There is no separate short authoring file or generated companion.

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

Legacy `reproduction` metadata remains unchanged. New v4 records may place a
standard `assessment` inside `evidence.reproduction`:

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
