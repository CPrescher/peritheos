# EOSMAT v5: uncertainties in parameter units

EOSMAT v5 is the native format for new files and the bundled material library.
It keeps the readable v4 layout and introduces one uniform uncertainty rule:

**An error uses the same unit as its parameter value. A covariance entry uses
the product of the units of its two parameters.**

There are no error-unit switches, registration requirements or companion files.
The [authoring guide](eosmat-authoring.md) contains complete working examples.
Identity, models, coefficients, units, sources and limits remain in one file.

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

For older software, export the version its reader supports:

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

The catalog document API keeps its historical v3-shaped dictionary view:
`get_material_document(identifier)` preserves source coefficient units with
explicit legacy annotations. Request `format_version=5` or `None` for the native
file. For maintenance views only,
`convert_eosmat(..., version=3, legacy_uncertainty_units=False)` retains source
uncertainty units and annotations; that view requires a current reader.

## Other fields and reading order

All other [v4 field definitions](eosmat-v4.md), model identifiers, optional
metadata sharing, source evidence and reconstruction assessment semantics remain
supported. Ordinary files stay inline. Saved v5 files use the same
[human-readable order](eosmat-v4.md#reading-order-in-saved-files), placing model,
reference state and coefficients before detailed evidence.

`parameter_error_volume_basis` and `parameter_covariance.volume_basis` are
prohibited in v5 records: the version already specifies the units. The schema
is included in the installed package as `eosmat_schema(5)` and
`peritheos/data/eosmat-v5.schema.json`.
