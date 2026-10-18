# Creating and editing EOSMAT

EOSMAT v4 is the native format for new material files. It supports small
standalone user definitions and complete curated database records in the same
`.eosmat` file. There is no separate authoring format, database registration or
required companion file. See the [v4 format reference](eosmat-v4.md).

## Create a material

```bash
python -m peritheos create my_material.eosmat \
  --formula X --model BM3 --temperature 300 \
  --parameters V0=10 K0=160 K0_prime=4
```

The coefficients above are fictional. Replace them with your material's values.
Version 4 fixes the exchange units to GPa, kelvin and **conventional-cell Å³**.
The reference temperature is supplied explicitly. The command creates a native
v4 file of about 14 lines and refuses to overwrite an existing destination.
Use `--id` if the filename stem is not lower snake case, `--name` for a display
name and `--reference` to record a source citation.

Alternatively, edit the [small native template](examples/simple-material.eosmat)
directly. It has a fictional fixed pressure derivative and remains unreviewed.

Python provides the same native file structure:

```python
from peritheos import create_eosmat, save_eosmat

document = create_eosmat(
    "my_material", "X", "BM3",
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

```bash
python -m scripts.format_eosmat my_material.eosmat
python -m scripts.format_eosmat --check my_material.eosmat
```

For an application whose pinned reader supports only v3:

```bash
python -m peritheos convert my_material.eosmat compatible.eosmat --format-version 3
```

Use `--format-version 3` with `create` if that older format is required directly.
Full [v3 compatibility fields](eosmat-schema.md) remain documented. New readers
accept old files, while old readers need the explicit compatibility export.

## Inspect database records

Use the [readable material summaries](material-summaries/index.md) to compare
existing records and qualifications. They remain generated views of the
authoritative catalog files. Database-only summary/coverage regeneration is
independent of standalone authoring:

```bash
python -m scripts.generate_material_summaries --check
python -m scripts.generate_catalog_coverage --check
```
