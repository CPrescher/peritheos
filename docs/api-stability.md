# API stability

Peritheos follows [Semantic Versioning](https://semver.org/). Before 1.0, a
minor release may include a necessary breaking correction, particularly when a
published equation or unit convention was implemented incorrectly. Such a
change must be prominent in the changelog and release notes.

## Public API

The supported public API consists of:

- names documented in the [API reference](api.md);
- names exported by a module's `__all__`;
- constructor parameters and documented public methods of exported EOS classes;
- the documented GPa, kelvin, and molar-volume conventions; and
- scalar and NumPy broadcasting behavior described in the documentation.

Names beginning with an underscore, implementation details not documented in
the API reference, and exact optimizer messages are private.

The exception hierarchy and documented machine-readable error codes are
public. Existing built-in exception compatibility is retained; exact message
text is diagnostic and is not a compatibility surface. See
[Error handling](error-handling.md) for the hierarchy and source-chain rules.

The Python distribution includes a `py.typed` marker. Public inline type
annotations are therefore part of the supported developer experience, although
improvements that make an annotation more precise are not runtime API changes.

## `.eosmat` compatibility

The canonical `.eosmat` format version 3 is a public exchange contract.
Readers must preserve unknown optional fields, and additive optional fields do
not require a format-version increment. Removing a field, changing its unit or
meaning, or making an optional field mandatory requires a new format version
and a documented migration. Native Dioptas 0.10.0 format-2 material documents
remain supported as legacy input.

Stable material and record identifiers are not reused for a different phase or
parameterization. A scientific correction creates an auditable data change;
saved project reproducibility must not depend on silently resolving an old
identifier to altered parameters.

The bundled `.eosmat` identifiers are canonical for the executable catalog.
Historical pressure-scale identifiers and module-level convenience constants
remain supported as explicit compatibility lookups, preserving their prior
parameters and numerical results when a primary-source audit found that the
new canonical record is not identical. They are intentionally omitted from
`list_materials()`, `list_eos_records()`, and search results, whose deterministic
contents are exactly the bundled 145-material/233-record collection.

The migration audit covered all 37 historical convenience records: 9 have an
equation-identical canonical counterpart, 17 have a counterpart whose audited
parameters, reference state, or equation composition differs, and 11 have no
record in the bundled 233-record collection. Consequently, the old names are not
blindly redirected to scientifically different records. Their compatibility
objects remain numerically stable; all new discovery and canonical-identifier
lookups use document-built objects. An exhaustive internal manifest is tested
against both catalogs so no historical record can be added, removed, or
redirected without an explicit audited disposition.

Canonical identifiers take precedence on collision. The only material-level
collision is `diamond`, which therefore resolves to the complete seven-record
document-built material instead of the former five-record convenience grouping.
The historical diamond records remain available individually through their
constants and record identifiers. The two record-level identifier collisions
are equation-identical anchored variants and retain the same numerical results.

An equation discriminator is part of a record's scientific meaning. Readers
must reject an unknown discriminator rather than substituting a related model.
In particular, `MieGruneisenDebye` with `debye_temperature_law:
variable_exponent` cannot fall back to its `integrated_gruneisen` default: the
characteristic-temperature laws differ for nonzero `q`.

## Compatibility commitments

`volume()` is the preferred pressure-to-volume method. `calculate_volume()` is
a supported compatibility alias and will remain available throughout the 1.x
series.

`temperature()` is the preferred pressure-and-volume-to-temperature method.
`calculate_temperature()` is its supported compatibility alias and will remain
available throughout the 1.x series.

After 1.0, planned public API removals will normally emit
`DeprecationWarning` for at least two minor releases. Scientific corrections
that cannot preserve old behavior will include a reproducible before-and-after
case and a migration note.

Patch releases may improve numerical precision, validation, warning text, or
solver robustness without treating the last floating-point bit or exact error
message as stable. Published calculations should record the Peritheos version,
model order, complete parameters, units, and uncertainty convention.
