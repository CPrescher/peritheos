# Bundled material library

This directory contains all 116 materials with EOS records from the
120-material, 147-EOS-record Dioptas 0.10.0 database, tag commit
`5a8bfd81d10bfab3499039603380aae34576d60a`. Dioptas is distributed under the
MIT License. Its project source is <https://github.com/Dioptas/Dioptas>.

The migration preserves the Dioptas crystallographic and EOS data and adds
stable identifiers plus explicit migration provenance. It does **not** make
Dioptas the scientific authority for an EOS record. The primary-source audit
dated 2026-08-31 classifies all 147 records: 116 are
`primary_source_validated`, while 31 are `deferred` with a concrete missing or
ambiguous evidence reason. No bundled record remains pending. The complete
machine-readable ledger is `../primary-source-audit.json`.

Only validated records are executable by default. Deferred records remain in
the files so Dioptas and other consumers can preserve the catalog without
silently converting inherited values into Peritheos-endorsed pressure scales.

The Fei et al. (2007) Au and Ne Debye-temperature laws are a documented
exception to byte-for-byte preservation. They are corrected from generic
implicit integrated behavior to `MieGruneisenDebye` with
`debye_temperature_law: variable_exponent`, based on equation 3 and the
definition immediately following it. Each affected record retains the source
value and primary-source location under `migration_corrections`. The audit also
restores the published `V0 = 35.12(2) angstrom^3` uncertainty for the Hanfland
graphite record and records, without hiding, the conflicting BM2/BM3 wording
in the Somayazulu B4C article. The Hazen--Finger zircon record is corrected
from BM2 to its published BM3 parameters with fixed `K0' = 6.5`. It also
restores model-required `n`/`Z` values
for Sokolova records, `n` for the Sun silica Debye records, and `Tr` for the
Bezacier ice records, each with field-level `audit_corrections` provenance.

The Dioptas-facing thermal type `AlphaKT` is preserved for interchange, while
its canonical mechanism identifier is `thermal_reference_state`. Peritheos
evaluates validated instances with `ThermalReferenceStateEOS` rather than
making the source/application label part of the public equation name.

Four Dioptas structure-only entries (`fe_fcc`, `fes_iii`, `nitrogen_epsilon`,
and `o8`) are intentionally not bundled because they contain no EOS record.
The `.eosmat` schema still permits an empty `eos_records` array for
application-created structure-only documents.

The `.eosmat` documents use the Peritheos-owned format version 3. Their
material, structure, and record layout is an additive evolution of the Dioptas
0.10.0 format-2 layout, so Dioptas 0.10.0 can read and preserve them. A future
Dioptas writer must preserve the version-3 top-level fields for a lossless
round trip and honor the variable-exponent law before numerically
evaluating the corrected Fei records. In the shared schema, an omitted
`debye_temperature_law` means `integrated_gruneisen`.
