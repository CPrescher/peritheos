# EOSMAT v5 format overview

An `.eosmat` file is one JSON document describing one material: its identity,
fixed units, one or more EOS records with their coefficients, errors,
reference state, source and limits, and optionally the observation tables
behind them. EOSMAT v5 is the format of new files and of the bundled library.
Every field is defined in the [record reference](eosmat-reference.md); the
[authoring guide](eosmat-authoring.md) shows how to write your own file.

## A small standalone file

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

```python
from peritheos import Material, load_eosmat

document = load_eosmat("my_material.eosmat")
material = Material.from_eosmat(document, require_primary_validation=False)
print(material.eos_records[0].pressure(70.0))
```

The flag accepts records not yet checked against their source, as every file
you write yourself is. `V0` and its error are in Å³ per conventional cell.

<a id="record-check-headline-and-audit-sidecars"></a>

## The `check` block

A curated record may carry a short `check` block: its reproduction headline
(`reproduction`, `compared`, `why`), any usage `warning`, the source check
(`source_check`, `source_checked_on`, `source_note`) and an `audit` link to its
sidecar entry. When evidence is present, the headline must agree with it.
Field values are in the [record reference](eosmat-reference.md#the-check-block);
the words are defined in [EOS labels](eos-labels.md). Your own files need no
`check` block, and loading never gives a record a positive headline.

## Default records

A material names its default records once, directly before `eos_records`,
for example `"default_records": {"equilibrium": "my_material_bm3"}`; the keys
are `equilibrium` and `hugoniot`. Per-record `default` and `default_for` flags
remain valid and must agree with the map. With neither,
`Material.default_record()` uses the first record.

## What changed from v4

| Topic | v5 |
| --- | --- |
| Uncertainty units | An error uses the unit of its parameter value; a covariance entry the product of its two parameters' units ([details](eosmat-reference.md#covariance)). `parameter_error_volume_basis` and `parameter_covariance.volume_basis` are not allowed. |
| Reproduction headline | Optional record `check` block. |
| Audit trail | May live in a per-material [audit sidecar](eosmat-reference.md#audit-sidecar). |
| Default selection | Material `default_records` map. |
| Reading order | `check` follows the label; the [field table](eosmat-reference.md#eos-record-fields) gives the order. |

<a id="audit-sidecars"></a>
<a id="complete-example-record-file-and-audit-entry"></a>

## Loading, missing sidecars and exports

The bundled library stores each material in `peritheos/data/materials/<id>.eosmat`
and its audit trail in `peritheos/data/audits/<id>.audit.json`; the
[sidecar rules and a complete example](eosmat-reference.md#audit-sidecar) are in
the reference.

| Call | Sidecar handling |
| --- | --- |
| `load_eosmat(path)` | Merges every sidecar referenced by `check.audit`; returns the complete document. |
| `load_eosmat(path, include_audit=False)` | Returns the file exactly as written. |
| `get_material_document(id)` | Bundled material, merged; `include_audit=False` returns the record file alone. |
| `python -m peritheos validate` | Requires referenced sidecars to exist and merge cleanly. |
| `python -m peritheos format` | Formats the record file only. |
| `python -m peritheos convert` | Writes one complete file with all evidence inline. |

A referenced sidecar that does not exist raises `EosmatError` with
`code="eosmat.audit_missing"`. A record loaded without its sidecar keeps its
coefficients and `check` headline and can be evaluated; reading
`EOSRecord.audit` then raises `code="eosmat.audit_not_loaded"`. Bundled records
load their sidecar on first access to evidence. Maintainers split a bundled
file with the repository command `python scripts/split_eosmat_audits.py --add`
([adding materials](adding-materials-and-eos.md#5-create-or-update-the-eosmat-document)).

## Old files and compatibility exports

Python and Rust read v2, v3, v4 and v5 with their historical uncertainty
interpretations. Loading and saving keeps a file's version. Migrate a v3/v4
file, or export for an older reader such as Dioptas with a v3 loader:

```bash
python -m peritheos convert old.eosmat current.eosmat
python -m peritheos convert current.eosmat legacy-v4.eosmat --format-version 4 --inline
python -m peritheos convert current.eosmat legacy-v3.eosmat --format-version 3
```

Destinations are never overwritten. Conversion changes only uncertainty units
and layout, never a coefficient, observation or qualification. v3/v4 exports
(also `convert_eosmat(document, version=3)` or `Material.to_eosmat(format_version=4)`)
use historical model units for errors, write `default` on every record and
carry merged evidence inline without `check`; a record whose sidecar was not
merged keeps its `check`. Readers older than v5 must reject v5.
`get_material_document(id)` returns a v3-shaped view unless given
`format_version=5`. Older layouts: [v4](eosmat-v4.md), [v3](eosmat-schema.md).
The schema ships as `eosmat_schema(5)` and `peritheos/data/eosmat-v5.schema.json`.
