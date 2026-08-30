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
