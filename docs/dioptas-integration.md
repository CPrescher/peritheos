# Dioptas and `.eosmat`

Peritheos owns the canonical `.eosmat` material schema and bundled EOS
material library. Dioptas consumes that library and uses the optional
crystallographic fields to calculate and display diffraction reflections.
Peritheos preserves those fields but does not need them for EOS calculations.

## One flat material document

There is no nested "Dioptas material containing a Peritheos material". A
single document stores the shared identity, optional structure, and zero or
more EOS records:

```json
{
  "format": "peritheos.material",
  "format_version": 3,
  "identifier": "gold",
  "name": "Gold",
  "formula": "Au",
  "phase": "fcc",
  "symmetry": "CUBIC",
  "lattice": {
    "a": 4.0862,
    "b": null,
    "c": null,
    "alpha": 90.0,
    "beta": 90.0,
    "gamma": 90.0
  },
  "formula_units_per_cell": 4,
  "space_group": "Fm-3m",
  "space_group_number": 225,
  "atom_sites": [],
  "source": {},
  "peaks": [],
  "eos_records": []
}
```

Only `format`, `format_version`, `identifier`, `name`, `formula`, `units`, and
`eos_records` are required by the canonical format. `symmetry`, lattice, space
group, formula units per cell, atom sites, lossless CIF source, and fallback
peaks are optional. This permits both an EOS-only Peritheos material and a
complete diffraction-ready Dioptas material without introducing two kinds of
material.

The normative machine-readable definition is distributed as
`peritheos/data/eosmat-v3.schema.json`. Structural validation never executes a
Python implementation path or treats a missing reported error as zero.
See the dedicated [`.eosmat` schema reference](eosmat-schema.md) for the
field-by-field contract, equation tables, defaults, and a complete document.

## Bundled library

The package contains 116 `.eosmat` materials and 147 EOS records migrated from
[Dioptas 0.10.0](https://github.com/Dioptas/Dioptas/releases/tag/0.10.0), commit
`5a8bfd81d10bfab3499039603380aae34576d60a`. The four Dioptas structure-only
entries are not bundled because Peritheos's catalog is EOS-focused. The format
still permits `eos_records` to be empty for Dioptas-created structure-only
materials.

```python
from peritheos import get_material_document, list_material_documents

print(len(list_material_documents()))  # 116
gold = get_material_document("gold")
print(len(gold["eos_records"]))
```

The migration preserves Dioptas's structures, references, parameters, errors,
fixed-fit flags, ranges, notes, and record order. It adds stable material and
record identifiers plus exact migration provenance.

One primary-source correction is recorded rather than preserved incorrectly:
the [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) Au and Ne thermal records use
`MieGruneisenDebye` with `debye_temperature_law: variable_exponent`.
Dioptas 0.10.0 implicitly used the default `integrated_gruneisen` relation, but
Fei equation 3 and the definition immediately following it specify
`theta_D = theta0 * (V/V0)**(-gamma(V))`. Each changed record contains a
`migration_corrections` entry with the source value, corrected value, DOI, and
equation location. The independent audit subsequently validates the complete
Fei records; the correction alone would not have been sufficient.

For `MieGruneisenDebye`, an absent `debye_temperature_law` means
`integrated_gruneisen`. Writers should nevertheless store an explicit law for
new or edited records so the scientific choice is visible without consulting
the schema default.

Migration is not scientific validation. A separate, reproducible audit now
marks 116 transferred records `primary_source_validated` and 31 `deferred`; no
record remains pending. A record was promoted only after its equation, every
parameter, units, reference state, phase, uncertainty convention, and data
range were traced to the cited primary publication or official supplement.
The complete ledger is bundled as
`peritheos/data/primary-source-audit.json`; deferred records remain readable
but are not executable by default.

The audit restores model inputs omitted by the flat migration when primary
evidence supplies them (`n`/`Z` for Sokolova, formula-unit `n` for two silica
Debye records, and the ice reference temperature). It also normalizes the
Dioptas type `AlphaKT` to the mechanism-oriented Peritheos model identifier
`thermal_reference_state`; the application-facing type is retained for
cross-compatible round trips.

## Reading and writing

```python
from peritheos import load_eosmat, save_eosmat

material = load_eosmat("sample.eosmat")
save_eosmat("copy.eosmat", material)
```

The reader accepts:

- canonical Peritheos `.eosmat` format 3; and
- unmodified Dioptas 0.10.0 material format 2, which has the same flat shape
  but lacks the Peritheos format discriminator and stable identifiers.

The writer emits the document supplied to it after structural validation and
does not discard unknown optional fields. This preservation rule is important
for applications that use only one part of the document.

## Dioptas compatibility

[Dioptas 0.10.0's material loader](https://github.com/Dioptas/Dioptas/blob/0.10.0/dioptas/model/eos/material.py)
already accepts the canonical format-3 files because the material and record
shape is an additive extension of its format 2. Its
[`EosPhase`](https://github.com/Dioptas/Dioptas/blob/0.10.0/dioptas/model/util/eos_phase.py)
uses Peritheos as the numerical engine and performs the conventional-cell to
molar-volume conversion required by thermal and Holzapfel models.

All 116 bundled Peritheos documents and all 147 records were loaded with the
Dioptas 0.10.0 `Material` implementation as an interoperability regression.
Dioptas 0.10.0 preserves complete EOS record dictionaries but rewrites the
top-level version as 2 and does not dispatch `debye_temperature_law`. The
coordinated Dioptas development update now preserves the format discriminator,
identifier, units, phase, cell definition, crystallography, and unknown
format-3 extensions. It also passes `debye_temperature_law` to Peritheos,
rejects unknown values, and records `variable_exponent` on its Fei Au and Ne
entries. That numerical update requires Peritheos 0.5+; the Dioptas lockfile
must be refreshed after the Peritheos release is available. Other records that
omit the field intentionally retain the integrated default.

Runtime ownership remains separate. Whether a record is bundled, loaded from
a file, read-only, duplicated, or locally edited is Dioptas application state
and is not stored as a scientific material role.

## Scientific and unit boundary

- Pressure is GPa and temperature is K.
- Dioptas-facing reference volumes are conventional-cell Å³.
- `formula_units_per_cell` links that cell to formula-molar thermal models.
- Structure lattice parameters and a record's EOS reference volume need not be
  numerically identical unless their stated reference conditions are also
  identical.
- Experimental fit ranges are source-data coverage, not phase-stability
  ranges.
- JSON `null` uncertainty means not reported or not verified, never zero.
- Unknown model types remain data; a reader must not coerce them to a different
  equation.

Peritheos remains responsible for the EOS model registry and scientific record
semantics. Dioptas remains responsible for interpreting structures,
calculating reflections, and storing GUI/project state.
