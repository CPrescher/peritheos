# EOSMAT v5

EOSMAT v5 is the native format for new files and the bundled material library.
It keeps the readable v4 layout and adds one uniform uncertainty rule:

**An error uses the same unit as its parameter value. A covariance entry uses
the product of the units of its two parameters.**

A user file is one self-contained `.eosmat` file. Curated records may add two
optional pieces: a short [`check` block](#the-check-block) with the record's
reproduction headline, and an [audit sidecar](#audit-sidecars) that holds the
audit trail in a second file. Neither is needed to write or use your own EOS.
The [authoring guide](eosmat-authoring.md) has complete editable examples.

## Small standalone file

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
      "reference": "Fictional example; no experimental validation."
    }
  ]
}
```

Here both `V0.value` and `V0.error` use Å³ per conventional cell. Both `K0.value`
and `K0.error` use GPa. Adding a molar-volume thermal model does not change those
file units. Peritheos converts the reference volume, its error and the relevant
covariance rows and columns together when constructing the numerical model.
Thermal coefficients retain their model's documented coefficient units.

## Covariance

For the example above, a covariance may be written as:

```json
"parameter_covariance": {
  "parameter_order": ["V0", "K0"],
  "matrix": [[0.0004, 0.01], [0.01, 4.0]]
}
```

The diagonal units are Å⁶ and GPa². The off-diagonal entries use Å³·GPa.
In a thermal record, `rt_eos.V0` and `rt_eos.K0` identify reference-isotherm
parameters; their covariance still uses the units of the corresponding values
in the file's reference `parameters` block. Unprefixed names identify the
thermal model's own coefficients.

If internal volume is `s * V0`, its error becomes `s * error`. A covariance
row or column associated with that volume gets one factor of `s`; its variance
gets two. Confidence levels, fixed flags, correlation, and the source meaning of
the reported errors remain separate and unchanged. Missing or `null` errors
remain unknown. A reported zero is preserved and does not establish exactness.
Negative errors are invalid.

<a id="record-check-headline-and-audit-sidecars"></a>

## The `check` block

A record may end with a `check` block. Every field is optional; unknown fields
are rejected.

| Field | Allowed values |
| --- | --- |
| `reproduction` | `reproduced`, `close`, `disagrees` or `untested` |
| `warning` | `null` or non-empty text: an explicit usage restriction |
| `finding` | Non-empty text: the source check's finding |
| `source` | `primary_source_validated`, `pending_primary_source_check`, `not_reproduced` or `deferred` |
| `source_checked_on` | Non-empty text, the date of the source check (`YYYY-MM-DD` in bundled records) |
| `audit` | Non-empty text: the record's entry in an audit sidecar, `<path relative to this file>#<record identifier>` |

The words are defined in [EOS labels](eos-labels.md). `reproduction` is a fixed
mapping of `evidence.reproduction.assessment.status`:

| Assessment status | `reproduction` |
| --- | --- |
| `parity` | `reproduced` |
| `similar` | `close` |
| `not_reproduced` | `disagrees` |
| `not_refittable`, `source_reconstruction`, `not_assessed`, or no assessment | `untested` |

When the evidence is present (inline or merged from a sidecar), `reproduction`
and `source` must agree with it, or loading fails. `finding`, `source` and
`source_checked_on` copy `scientific_validation.primary_source_check.finding`,
`scientific_validation.status` and `scientific_validation.audit_date`.
In bundled records, `warning` is composed only from coded fields, in this order,
with fixed wording: `evidence.reproduction.do_not_use: true` ("Do not use."), a
coded `scientific_validation.usage_recommendation` such as
`not_recommended_for_quantitative_use` ("Not recommended for quantitative use."),
a material `recommendation_hold` naming the record ("Default selection on hold
pending review.") and `catalog_access: "explicit_selection"` ("Explicit
selection only; not loaded by default."). Free-text recommendations are never
turned into a warning.

Without a `check` block, Python derives the same `reproduction` and `warning`
from the evidence (`EOSRecord.reproduction`, `EOSRecord.warning`); a record with
no assessment is `untested`. Loading never adds a positive headline.

## Audit sidecars

A sidecar holds the audit trail of one material's records, keyed by record
identifier:

```json
{
  "format": "peritheos.eosmat-audit",
  "format_version": 5,
  "material": "<material identifier>",
  "records": {
    "<record identifier>": {"evidence": {}, "pressure_calibration": {}}
  }
}
```

`material` must equal the record file's `identifier`, and every key of
`records` must name a record in that file. An entry may contain only these
blocks, which keep the meaning they have inline: `evidence`,
`pressure_calibration`, `fit_input_inventory`, `network_reconstruction_audit`,
`reproduction_audit`, `numerical_reproduction`, `validation_selection`,
`migration_corrections`, `scope_history` and `author_input_followup`. Every
other record field stays in the record file. Merging an entry into its record
gives the complete inline record; a block present in both places is an error.

A record uses its sidecar entry when `check.audit` names it and the record holds
none of those blocks inline. The path before `#` is resolved relative to the
record file; the fragment, when present, must equal the record identifier.

On merge, an assessment's `report_url`, `report_path` and `report` are filled
in when absent, with the values implied by the material and record identifiers
(the bundled [reconstruction report](reconstruction-reports/index.md) anchor and
`docs/data/primary-eos-refits.json`). Bundled sidecars therefore omit them;
stored values that differ are kept.

The bundled library stores each material as
`peritheos/data/materials/<id>.eosmat` (header, EOS records with their `check`
blocks, then datasets) and `peritheos/data/audits/<id>.audit.json`. Generic
caveats are stated once in [EOS labels](eos-labels.md) rather than per record;
the receipt `docs/data/eosmat-v5-audit-split.json` lists the three removed
statements and their original positions. A bundled material written later as
a complete inline file is split by the same procedure with the repository
command `python scripts/split_eosmat_audits.py --add`
([adding materials](adding-materials-and-eos.md#5-create-or-update-the-eosmat-document)).

## Loading, missing sidecars and exports

| Call | Sidecar handling |
| --- | --- |
| `load_eosmat(path)` | Merges every sidecar referenced by `check.audit`; returns the complete document. |
| `load_eosmat(path, include_audit=False)` | Returns the file exactly as written. |
| `get_material_document(id)` | Bundled material, merged; `include_audit=False` returns the record file alone. |
| `python -m peritheos validate` | Requires referenced sidecars to exist and merge cleanly. |
| `python -m peritheos format` | Formats the record file only. |
| `python -m peritheos convert` | Writes one complete file with all evidence inline. |

A referenced sidecar that does not exist raises `EosmatError` with
`code="eosmat.audit_missing"`. A record loaded without its sidecar still has its
coefficients and `check` headline and can be evaluated; reading
`EOSRecord.audit` then raises `EosmatError` with `code="eosmat.audit_not_loaded"`.
Bundled catalog records evaluate and report their headline without reading the
sidecar, and load it on first access to evidence-backed attributes.

Standalone files need no sidecar: inline `evidence` remains fully supported,
`check` is optional, and `EOSRecord.audit` returns the inline blocks.

## Complete example: record file and audit entry

This fictional record file is `materials/example_oxide.eosmat`:

<!-- example: record file -->
```json
{
  "format_version": 5,
  "identifier": "example_oxide",
  "formula": "MgO",
  "units": {"pressure": "GPa", "temperature": "K", "volume": "angstrom^3/conventional_unit_cell"},
  "eos_records": [
    {
      "identifier": "example_oxide_bm3",
      "model": "BM3",
      "temperature_ref": 300,
      "parameters": {
        "V0": {"value": 74.7, "error": 0.02},
        "K0": {"value": 160, "error": 2},
        "K0_prime": {"value": 4, "fixed": true}
      },
      "validity": {"pressure_gpa": [0, 30]},
      "reference": "Fictional example; no experimental validation.",
      "check": {
        "reproduction": "reproduced",
        "compared": "published_fit",
        "warning": null,
        "source_check": "pending",
        "audit": "../audits/example_oxide.audit.json#example_oxide_bm3"
      }
    }
  ]
}
```

Its sidecar is `audits/example_oxide.audit.json`:

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
        "reproduction": {
          "assessment": {
            "status": "parity",
            "scope": "parameter_agreement",
            "basis": "Fictional: a refit of example observations recovers V0 and K0 under the generic criterion.",
            "parameters": [
              {"parameter": "V0", "published": 74.7, "published_error": 0.02, "refit": 74.71, "refit_error": 0.02},
              {"parameter": "K0", "published": 160, "published_error": 2, "refit": 159.2, "refit_error": 1.8}
            ],
            "fit": {"free_parameters": ["V0", "K0"], "fixed_parameters": ["K0_prime"]}
          }
        }
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

`parity` maps to `reproduced`, and `source` repeats the evidence status, so the
two files are consistent. The record remains an unreviewed teaching example.
Writing the same blocks inline under the record, without `check.audit`, gives
an equivalent single file.

## Old files and compatibility exports

Python and Rust continue to read v2, v3 and v4. Their historical uncertainty
interpretations are retained, including the optional v3/v4 public/model volume
annotations. Loading and saving a file preserves its version; it does not
silently migrate the file. New standalone creation defaults to v5. Newly
constructed `Material` objects also default to v5; loaded v3/v4/v5 objects retain
their input version. The older v2 loader supports export through a `Material`:
`Material.from_eosmat(document, require_primary_validation=False).to_eosmat(format_version=5)`.
Raw `load_eosmat()` / `save_eosmat()` also preserve v2 without migration.

Migrate a v3/v4 file with the installed package:

```bash
python -m peritheos convert old.eosmat current.eosmat
```

The destination is protected against overwriting. Conversion changes only the
required unit representation and format layout; it performs no fit and changes
no coefficients, observations, source qualifications or confidence semantics.

For older software, including Dioptas versions with a v3 reader, export the
version its reader supports:

```bash
python -m peritheos convert current.eosmat legacy-v4.eosmat --format-version 4 --inline
python -m peritheos convert current.eosmat legacy-v3.eosmat --format-version 3
```

`convert_eosmat(document, version=3)` and `version=4`, and the corresponding
`Material.to_eosmat(format_version=...)` calls, convert uncertainties back to
historical model units and omit the newer unit-switch fields. Old readers
therefore receive the numerical convention they expect. Readers predating v5
must reject v5 rather than guess its uncertainty semantics. Compatibility does
not add support for a model that an older consumer never implemented.

v3 and v4 exports of a complete document (any sidecar merged, as `load_eosmat`,
`get_material_document` and the `convert` command do by default) carry the full
evidence inline and omit `check`. A record whose evidence is still in an
unmerged sidecar keeps its `check` block in the export, so the audit reference
is not lost.

The catalog document API keeps its historical v3-shaped dictionary view:
`get_material_document(identifier)` preserves source coefficient units with
explicit legacy annotations and includes the merged evidence. Request
`format_version=5` or `None` for the native file. For maintenance views only,
`convert_eosmat(..., version=3, legacy_uncertainty_units=False)` retains source
uncertainty units and annotations; that view requires a current reader.

## Other fields and reading order

All other [v4 field definitions](eosmat-v4.md), model identifiers, optional
metadata sharing, source evidence and reconstruction assessment semantics remain
supported. Saved v5 files use the same
[human-readable order](eosmat-v4.md#reading-order-in-saved-files), placing model,
reference state and coefficients before detailed evidence.

`parameter_error_volume_basis` and `parameter_covariance.volume_basis` are
prohibited in v5 records: the version already specifies the units. The schema
is included in the installed package as `eosmat_schema(5)` and
`peritheos/data/eosmat-v5.schema.json`.
