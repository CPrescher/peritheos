# Creating and editing EOSMAT

EOSMAT v5 is the native format for new material files. Your own EOS is one
self-contained `.eosmat` file: identity, units, model, coefficients, source and
limits. There is no separate authoring format, database registration, `check`
block or audit sidecar to write. Curated records may add those optional pieces;
see the [v5 format reference](eosmat-v5.md).

## Read a file before editing it

Read the file from the top, in this order. You do not need database audit
metadata to create an EOS for your own calculations.

| Question | Field to read |
|---|---|
| Which material and phase? | `name`, `formula`, optional `phase`; `identifier` is its stable machine name |
| What are the units and cell contents? | `units`, and `formula_units_per_cell` when a model needs it |
| Which source or alternative am I looking at? | The record's `label`, `identifier` and `reference` |
| Which equation and coefficients? | `model`, `temperature_ref`, `parameters`; optional `thermal` |
| What do the errors mean? | Each parameter's `error`/`fixed`, then confidence convention and covariance if available |
| Where may I use it? | `validity`; measured coverage is separate |
| Was it checked? (curated records) | `check`: the reproduction headline and what it compared, any warning and the source check |
| Which record is the default? | The material's `default_records` (or a record's `default` flag) |

For catalog display names, put the mineral family first and keep useful sample,
structure, spin or model qualifiers after it. Store the exact composition in
`formula`; applications may append it to `name` for display. For example,
`name: Bridgmanite natural` and `formula: Mg0.90Fe0.09Al0.005Ca0.005SiO3`
describe the existing natural sample without embedding its formula twice.
Retain previous names and useful composition descriptions in material-level
`aliases`. See the [catalog naming conventions](catalog.md) for mixed Mg–Fe
oxide synonyms. Display changes must preserve stable identifiers and scientific
content.

Start with one of these complete, editable files:

- [Basic material](examples/simple-material.eosmat): one EOS, explicit units,
  a fixed coefficient and a clearly marked fictional citation.
- [Thermal material](examples/thermal-material.eosmat): a cold EOS plus Debye
  thermal pressure, explicit cell contents and explicit thermal laws.
- [Database-style material](examples/curated-material.eosmat): two alternatives,
  uncertainty/provenance, validity and a small calculated-checkpoint table.

Every example uses fictional teaching coefficients and remains unreviewed.
Calculated checkpoints are labeled as model output; they are not experimental
observations or independent reproduction evidence. These files are examples to
edit, not recommended pressure standards.

## Create a material

```bash
python -m peritheos create my_material.eosmat \
  --formula X --model BM3 --temperature 300 \
  --parameters V0=10 K0=160 K0_prime=4
```

The coefficients above are fictional. Replace them with your material's values.
EOSMAT fixes the exchange units to GPa, kelvin and **conventional-cell Å³**.
The reference temperature is supplied explicitly. The command creates a native
v5 file with its units visible in the header and refuses to overwrite an existing destination.
Use `--id` if the filename stem is not lower snake case, `--name` for a display
name and `--reference` to record a source citation.

Alternatively, edit the [small native template](examples/simple-material.eosmat)
directly. It has a fictional fixed pressure derivative and remains unreviewed.

Python provides the same native file structure:

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

The creator checks that the selected model can be constructed. It does not fit
observations or establish scientific validation. New records remain diagnostic
and unreviewed, even when you supply a citation. To deliberately use your own
unaudited parameters:

```python
from peritheos import Material, load_eosmat

document = load_eosmat("my_material.eosmat")
material = Material.from_eosmat(document, require_primary_validation=False)
volume = material.eos_records[0].volume(10.0)
```

## Add scientific information when available

Errors and fixed flags are adjacent to parameter values in the actual file:

```json
"parameters": {
  "V0": 10,
  "K0": {"value": 160, "error": 2},
  "K0_prime": {"value": 4, "fixed": true}
}
```

Missing and explicit `null` errors remain unknown; a supplied zero stays zero.
The creation API accepts `parameter_errors={"K0": 2}` and
`fixed_parameters=["K0_prime"]`; the command accepts `--errors K0=2 --fixed K0_prime`.
These are input conveniences; the saved file uses the native annotations above.

In v5, every error uses the same unit as its parameter value. The reference
`V0` and its error both use conventional-cell Å³, including for molar-volume
thermal models. Internal conversions are handled together by the reader.
No uncertainty unit-switch fields are needed or allowed in v5 records.

This complete example uses fictional coefficients and a cell-volume error:

```python
from peritheos import Material, create_eosmat

document = create_eosmat(
    "my_thermal",
    "MgO",
    "BM3",
    {"V0": 74.7, "K0": 160, "K0_prime": 4},
    temperature_ref=300,
    parameter_errors={"V0": 0.02},  # Å³ per conventional cell
    formula_units_per_cell=4,
    thermal={
        "model": "MieGruneisenDebye",
        "parameters": {"Tr": 300, "theta0": 760, "gamma0": 1.5, "q": 1, "n": 2},
    },
)
record = Material.from_eosmat(document, require_primary_validation=False).eos_records[0]
prediction = record.pressure_with_uncertainty(70, 1000, volume_sigma=0.01)
print(prediction.value, prediction.standard_error)
```

Each covariance entry uses the product of the corresponding parameter units:
V0 variance uses Å⁶, K0 variance uses GPa², and a V0–K0 cross term uses Å³·GPa.
Use `rt_eos.V0` for reference-isotherm volume in a thermal record. State the
matrix's source, parameter order and confidence assumptions separately.
Missing errors remain unknown; zero is preserved with a precision qualification,
and negative errors are invalid. See [the v5 uncertainty rule](eosmat-v5.md).

Old v2/v3/v4 files keep their original interpretations. Explicit v3/v4 exports
convert errors and covariance back to historical model units for older readers.
Use `convert_eosmat(document, version=3)` or `version=4`, or the corresponding
`material.to_eosmat(format_version=...)` call. Same-version v5 conversion
retains literal source values. No fit or scientific reclassification is run.

A thermal component declares its model once and retains all formulation
options. In Python, pass `thermal` and supply `formula_units_per_cell` for
molar-volume models. Command-line equivalents are `--thermal-model`,
`--thermal-parameters` and `--z`. Cell contents are never inferred from the
formula. Model-specific law choices use the Python component mapping or direct
file editing. Advanced Hugoniot states use the complete schema.

For the published garnet q-compromise pressure approximation, use the named
thermal model `QCompromiseDebye` with `Tr`, `theta0`, `gamma0`, `n` and a known
crystallographic Z. Do not represent it as full MGD with a chosen q. See the
[source audit](literature-reproductions/garnet-angel-2022.md) and
[executable example](examples/garnet-q-compromise.py).

An explicitly supported finite-pressure BM2 law uses `BM2Reference` with V0,
K0 and P0. V0 and K0 are values at P0; retain the source's actual reference
notation in its provenance. Do not put high-pressure reference values in an
ordinary `BM2` record. A complete, unaudited demonstration is
[bm2-reference.eosmat](examples/bm2-reference.eosmat):

```python
from peritheos import Material, load_eosmat

document = load_eosmat("bm2-reference.eosmat")
material = Material.from_eosmat(document, require_primary_validation=False)
record = material.get_eos_record("finite_reference_demo")
assert record.pressure(375.04) == 32.7
```

The published Bykova record has separate source review, joint-phase scope and
use limitations.

Add crystallography, more EOS records, source observations, covariance,
calibration and detailed `evidence` to the same file. Keep experimental coverage
distinct from permitted model evaluation, and preserve all known scientific
qualifications. Database contribution requires the additional review described
in [Adding materials and EOS](adding-materials-and-eos.md).

## Optional: reproduction headline and audit trail

A file needs neither of these to load or calculate.

- A record's `check` block states its reproduction headline (`reproduced`,
  `close`, `disagrees` or `untested`), what it `compared`, a generated `why`
  for `close` or `disagrees`, any usage `warning` and the `source_check`.
  When `evidence` is present, these must agree with it. Without `check`, the
  `EOSRecord` properties of the same names are derived from the evidence;
  without evidence the record is `untested` and compared `none`.
- `evidence.reproduction.assessment` holds the reconstruction outcome, basis and
  limitations; its `parameters` table holds the source-versus-reconstructed
  values, including failed reconstructions, and `fit` the selection, objective,
  weights and fixed constraints. The EOS's own `parameters` are never replaced.
  See the [reconstruction field reference](eosmat-v4.md#actual-reconstructed-coefficients).
- The bundled library keeps that audit trail in a per-material sidecar,
  `peritheos/data/audits/<id>.audit.json`, referenced by `check.audit`. Your own
  files can keep it inline or use the same layout.

The [v5 reference](eosmat-v5.md#complete-example-record-file-and-audit-entry)
has a complete record file with its audit entry. Words and caveats are defined
once in [EOS labels](eos-labels.md).

## Save, format and exchange

`load_eosmat()` reads v2/v3/v4/v5. `save_eosmat()` retains the input version, all
source values, record order and unknown material/record extensions. It writes
existing destinations so deliberate edits can be saved. The creation and
migration commands use new destinations to protect existing files.
Saved files use UTF-8 and LF line endings on every platform, preserving accented
names and scientific symbols as readable text.

```bash
python -m peritheos validate my_material.eosmat
python -m peritheos format my_material.eosmat
python -m peritheos format --check my_material.eosmat
```

These commands are included in the installed package; a repository checkout is
not needed. `validate` checks structure, model declarations and local references.
It does not establish experimental accuracy or fit reproduction. `format` changes
layout only, retaining every JSON value and array order; `--check` is read-only.

In v5, saved records put model, reference temperature, coefficients and thermal
terms first, followed by uncertainty confidence, experimental coverage, validity
and usage qualifications. Thermal blocks use `model`, then `parameters`, then
the explicit law/reference choices. Full covariance, citations, calibration and
audit evidence follow the practical information. Material-level recommendation
holds appear before the EOS list. These are formatting choices; field names and
scientific meanings are unchanged, and readers accept any JSON key order.
See the [native reading order](eosmat-v4.md#reading-order-in-saved-files).

For an application whose pinned reader supports only v3:

```bash
python -m peritheos convert my_material.eosmat compatible.eosmat --format-version 3
```

Use `--format-version 3` with `create` if that older format is required directly.
Full [v3 compatibility fields](eosmat-schema.md) remain documented. New readers
accept old files, while old readers need the explicit compatibility export.


## Repeated metadata in large database files

Keep ordinary user files inline. For a large file with identical citations or
source-check objects repeated across records, optional local sharing can reduce
repetition. It never needs a second file:

```bash
python -m peritheos convert material.eosmat compact.eosmat --compact
python -m peritheos convert compact.eosmat inline.eosmat --inline
```

The new destination is protected from overwriting. In a compact file,
`{"shared_ref": "source_name"}` means "use exactly this definition from the
file's `shared_metadata` table". There are no overrides, inherited parameters,
external lookups or implicit scientific decisions. Only citations, evidence and
pressure-calibration metadata can use this notation. Coefficients, observations
and unknown extensions outside those locations are untouched.

The [shared-metadata example](examples/shared-material.eosmat) contains the same
scientific content as the inline database-style example. Prefer the inline file
when that is easier to read. The first-generation v4 reader did not implement
local sharing: use `--format-version 4 --inline` for it, or `--format-version 3` for v3 consumers.
Both compatibility exports retain unused definitions as inert supporting metadata.

## Inspect database records

Use the [readable material summaries](material-summaries/index.md) to compare
existing records and qualifications. They remain generated views of the
authoritative catalog files. Database-only summary/coverage regeneration is
independent of standalone authoring:

```bash
python -m scripts.generate_material_summaries --check
python -m scripts.generate_catalog_coverage --check
```

The bundled library uses the same native v5 structure, with each material's
audit trail in its sidecar. `python -m peritheos convert INPUT OUTPUT` migrates a
standalone file and writes explicit compatibility exports with the installed
package. Catalog migration has repository-only, checksum-bound receipts in
`docs/data/eosmat-v5-migration.json` and `docs/data/eosmat-v5-audit-split.json`;
they are not required for user files.

Three historical refit scripts are fingerprinted by their scientific audit
reports. Their original code and reports remain unchanged. Use the repository
compatibility runner for their file access:

```bash
python -m scripts.run_eosmat_audit refit_hirose_fei_gold_speziale --check
python -m scripts.run_eosmat_audit refit_miozzi_2020_speziale --check
python -m scripts.run_eosmat_audit refit_miozzi_2020_tange --check
```

The runner gives the original script a temporary v3 export. A successful
`--register` saves any material changes back as v5, after checking for concurrent
edits. For its historical calibration fingerprints, its catalog lookup excludes
only the later standardized ledger assessment, retaining all original evidence,
coefficients and datasets. The archived hashes and audit code remain unchanged;
ordinary library lookups and compatibility exports retain the full assessment.
This compatibility step is only for these historical repository scripts;
normal file reading, authoring and calculation use v5 directly.

## Reading reconstruction reports

Each bundled record has a section in the
[reconstruction reports](reconstruction-reports/index.md). It opens with the
headline and any warning, then the scope and basis of the check, and shows only
what the audit records: method, curve and residual plots, and the coefficient
comparison. The assessment's `report_url` opens it on Read the Docs **latest**;
`report_path` names the same page and anchor in the repository. Bundled sidecars
omit both because they follow from the identifiers (see
[audit sidecars](eosmat-v5.md#audit-sidecars)). Standalone files can omit them or
link to their own reports.
