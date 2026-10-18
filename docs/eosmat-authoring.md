# Creating and editing EOSMAT

EOSMAT v4 is the native format for new material files. It supports small
standalone user definitions and complete curated database records in the same
`.eosmat` file. There is no separate authoring format, database registration or
required companion file. See the [v4 format reference](eosmat-v4.md).

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
| Where may I use it, and was it reconstructed? | `validity` and `evidence.reproduction`; measured coverage is separate |

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
Version 4 fixes the exchange units to GPa, kelvin and **conventional-cell Å³**.
The reference temperature is supplied explicitly. The command creates a native
v4 file with its units visible in the header and refuses to overwrite an existing destination.
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

A thermal component declares its model once and retains all formulation
options. In Python, pass `thermal` and supply `formula_units_per_cell` for
molar-volume models. Command-line equivalents are `--thermal-model`,
`--thermal-parameters` and `--z`. Cell contents are never inferred from the
formula. Model-specific law choices use the Python component mapping or direct
file editing. Advanced Hugoniot states use the complete schema.

Add crystallography, more EOS records, source observations, covariance,
calibration and detailed `evidence` to the same file. Keep experimental coverage
distinct from permitted model evaluation, and preserve all known scientific
qualifications. Database contribution requires the additional review described
in [Adding materials and EOS](adding-materials-and-eos.md).

## Save, format and exchange

`load_eosmat()` reads v2/v3/v4. `save_eosmat()` retains the input version, all
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

In v4, saved records put model, reference temperature, coefficients and thermal
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

For the meaning of parity, similarity, source checks, refits and diagnostics,
use [the central EOS label guide](eos-labels.md).

Bundled records put reconstruction outcome, basis and limitations at
`evidence.reproduction.assessment`. Its `parameters` table contains the actual
source-versus-reconstructed values, including failed reconstructions; its `fit`
object records the available selection, objective, weights and fixed constraints.
The EOS's main `parameters` are preserved. The readable material summaries show
the scoped assessment and a separate reconstructed-coefficient table. See the
[reconstruction field reference](eosmat-v4.md#actual-reconstructed-coefficients).
Standalone user files do not require an audit block.

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
local sharing: use `--inline` for it, or `--format-version 3` for v3 consumers.
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

The bundled library uses the same native v4 structure, with explicit units and
inline evidence. For repository maintenance, migrate older files in place with
`python -m scripts.migrate_eosmat_v4 [paths ...]`; without paths it selects the
bundled library. `--check` reports files needing migration without writing them.
The migration compares the complete scientific content before writing and is
safe to repeat. Catalog maintenance tools read either layout and write v4.
Use the installed `peritheos convert` command above for standalone files and
explicit compatibility exports.

Three historical refit scripts are fingerprinted by their scientific audit
reports. Their original code and reports remain unchanged. Use the repository
compatibility runner for their file access:

```bash
python -m scripts.run_eosmat_audit refit_hirose_fei_gold_speziale --check
python -m scripts.run_eosmat_audit refit_miozzi_2020_speziale --check
python -m scripts.run_eosmat_audit refit_miozzi_2020_tange --check
```

The runner gives the original script a temporary v3 export. A successful
`--register` saves any material changes back as v4, after checking for concurrent
edits. For its historical calibration fingerprints, its catalog lookup excludes
only the later standardized ledger assessment, retaining all original evidence,
coefficients and datasets. The archived hashes and audit code remain unchanged;
ordinary library lookups and compatibility exports retain the full assessment.
This compatibility step is only for these historical repository scripts;
normal file reading, authoring and calculation use v4 directly.

## Reading reconstruction reports

Bundled records link directly to their [human-readable reconstruction report](reconstruction-reports/index.md)
through `evidence.reproduction.assessment.report_url`. The link opens the specific
EOS on Read the Docs **latest** once the documentation is published. The companion
`report_path` identifies the same Markdown page and anchor relative to the repository;
local and version-specific documentation builds use relative links.

The reports show the recorded outcome, basis, limitations, method and coefficient
comparisons, including unsuccessful reconstructions. Available comparison pairs
have a graph of stored central values; detailed source reports supply observation
and residual plots where available. A missing reconstruction is explained rather
than illustrated with invented observations. Source coefficients remain unchanged.
Standalone files can omit these fields or link to their own reports.
