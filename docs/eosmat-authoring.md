# Reading and editing EOSMAT

Start with the [readable material summaries](material-summaries/index.md) to
compare EOS choices, parameters, reported errors and scientific qualifications.
Each sheet links to the complete material document. The EOSMAT file remains
the authoritative record; the summaries are generated views.

EOSMAT is still ordinary JSON, in `peritheos.material` version 3. Formatting
groups short parameter dictionaries and arrays together and puts material
identity and units first. Every scientific field, numerical value, array order
and extension remains in the file. Detailed audit history stays available.

## Start with one isothermal record

Copy [simple-material.eosmat](examples/simple-material.eosmat) to start a
material document. It contains a fictional BM3 isotherm with three coefficients,
missing errors represented by `null`, and a fixed pressure derivative. Its
validation status is deliberately `pending_primary_source_check`.

Replace the example identity, composition, cell contents, reference temperature,
coefficients and citation with source-supported values. Change the record kind
and determination method to describe the represented parameterization. Keep
the validation status pending until the primary source has actually been
checked. The template is an authoring starting point; it does not meet the
scientific acceptance requirements for a new bundled material by itself.

The main sections to edit are:

| Section | What belongs here |
| --- | --- |
| Material identity | Stable identifier, name, composition and phase |
| Units and cell contents | GPa, kelvin, conventional-cell Å³ and the correct formula-unit count |
| Reference | The actual primary citation and exact source location |
| `eos` | Both existing model identifiers and source coefficients |
| `parameter_errors` | Reported errors; `null` means unavailable, never zero |
| `fixed_parameters` | Parameters held fixed, independently of reported errors |
| `scientific_validation` | What was checked and any limitations on use |

Add measured coverage, calibration, parameter provenance, structure and
observations when supported by the source. Missing coverage does not authorize
extrapolation. A room-temperature isotherm does not establish a thermal model.
See [Adding materials and EOS models](adding-materials-and-eos.md) before
contributing a catalog record, and the [schema reference](eosmat-schema.md)
for thermal components, datasets, Hugoniots and advanced metadata.

## Save and format a document

The public writer validates the existing format and applies the readable layout:

```python
from peritheos import load_eosmat, save_eosmat

document = load_eosmat("my_material.eosmat")
save_eosmat("my_material.eosmat", document)
```

This does not change scientific validation status or select a different EOS.
An unaudited template is structurally valid but is rejected by the default
executable material loader. For deliberate inspection of fictional teaching
coefficients, the opt-out is explicit:

```python
from peritheos import Material

example = Material.from_eosmat(document, require_primary_validation=False)
```

Do not use that opt-out as a substitute for scientific review.

From the repository, format individual files or a directory with:

```bash
python -m scripts.format_eosmat my_material.eosmat
python -m scripts.format_eosmat --check my_material.eosmat
```

With no paths, the command covers the bundled materials and the simple template.
`--check` reports formatting drift without modifying files. Long source notes
remain intact even when their JSON string exceeds the usual line width. The
Rust serializer and other JSON writers can use different whitespace while
remaining compatible with the same v3 document.

## Regenerate the material summaries

```bash
python -m scripts.generate_material_summaries
python -m scripts.generate_material_summaries --check
```

The generator includes all stored records, including deferred evidence and
explicit diagnostic selections. It uses the executable catalog to report actual
automatic choices, including recommendation holds and pressure-scale reference
variants. It does not infer fit reproduction from a source-validation label.
Expand the scientific qualifications, pressure calibration and parameter
provenance sections for the underlying findings.

After rewriting bundled files, also regenerate the existing coverage snapshot:

```bash
python -m scripts.generate_catalog_coverage
python -m scripts.generate_validation_summary
```

That snapshot fingerprints the material file bytes, so whitespace changes its
fingerprint. Dataset checksums and observation-view bindings use the underlying
resources or normalized content and remain unchanged by formatting.
