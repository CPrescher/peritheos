# Bundled material library

This directory contains 141 curated materials with 223 EOS records. The
collection began with the 120-material, 147-EOS-record Dioptas 0.10.0 database,
tag commit
`5a8bfd81d10bfab3499039603380aae34576d60a`. Its project source is
<https://github.com/Dioptas/Dioptas>.

Dioptas is migration provenance, not a licensor of the underlying scientific
data. Many entries were assembled earlier in legacy JCPDS collections used in
GSECARS and related beamline workflows. Dioptas provenance does not establish
copyright ownership of, or relicense, numerical observations and EOS parameters
taken from beamline files and cited publications. Dataset-specific licenses are
recorded where known. Otherwise, a source citation records scientific provenance
rather than a copyright license.

The migration preserves supported Dioptas crystallographic and EOS data and adds
stable identifiers plus explicit migration provenance. It does **not** make
Dioptas the scientific authority for an EOS record. The primary-source audit
dated 2026-09-05 classifies all 223 bundled records as
`primary_source_validated`. No bundled record remains pending or deferred. The complete
machine-readable ledger is `../primary-source-audit.json`.

Additional records native to Peritheos have no invented Dioptas migration
provenance. They include the primary-sourced staged aragonite BM2 P-V-T
parameterization from Martinez et al. (1996), the Dewaele et al. (2012) B2-KCl
P-V-T pressure calibration, the complete Correa and Benedict diamond Helmholtz
models and experimental anchors, and the Tange (2009), Dewaele (2004), and
Takemura-Dewaele (2008) reference standards added for pressure-calibration
lineage, the Hemley neon refit, Campbell and Heinz's RbCl-B2 record, the audited
MgO, CaSiO3, stishovite, akimotoite, and Phase Egg records added through the
LitCurate intake, the selected Shi et al. Rh2O3(II)-type alumina thermal EOS,
the Fu et al. Fe-Al bridgmanite records, the phase-separated Sun et al. cubic
and tetragonal CaSiO3 records, and phase-restricted MgO and NiO Hugoniots. The
KCl record is the preferred
`kcl.eosmat` record and
keeps its measured 298 K range distinct from its molecular-dynamics thermal
extension. The Martinez staged result uses its exact Equation (3) direct-linear
reference-volume law. The paper's separate global thermal BM3 entry is excluded
because its fitted reference volume is omitted and the remaining coefficients
do not reproduce the printed dataset under the documented equations.

## Primary observation tables

When a reviewed primary paper prints a recoverable observation table, the
material document links that table through its top-level `datasets` array.
The current bundle contains 162 distinct primary datasets with 43,034 observation
rows, represented by 181 material-document links to 205 EOS records.
The Ono et al. cubic-SnO2 table is linked from both legacy diffraction-pattern
entries. The Shen--Smith Table S1 workbook is linked across its ten calibrant
and phase entries while retaining the simultaneously measured Cu reference
volumes even though the library has no separate Shen--Smith Cu EOS record.
Small tables may remain inline; larger tables are stored as
SHA-256-checksummed CSV resources under `../datasets/`. Quantities and
uncertainties retain the paper's reported conventions and units, including
formula-unit volumes, molar volumes, densities, lattice parameters, and
pressure-calibrant readings.

Dataset links describe the observations associated with a record; they do not
assert that every row entered the published regression. Row-selection details
are kept in dataset notes, for example for the stressed niobium run and the
full thermal ice-VII table. Conversely, the absence of a dataset is not a claim
that no underlying measurements ever existed: plots, unavailable supplements,
and papers that report only fitted parameters are not reconstructed into
invented point tables.

The sibling `../pressure-calibrations.json` registry contains executable,
versioned ruby R1 calibration records. Ruby-based EOS entries link to these
with `reference_calibration_record`; XRD-based entries use
`reference_eos_record` to link to the exact material EOS.

The final full-text audit resolved the earlier CsCl, magnetite, Li, majorite,
MW60, NiS, phase-D, cubic-SnO2, and SrO blockers. Phase D is intentionally two
EOS records because the primary paper reports distinct AntA and AntB ambient
volumes. The lithium fit is explicitly labeled as one empirical Vinet curve
spanning bcc and fcc observations. The cubic SnO2 records expose only the
published 300 K reference isotherm; a separate single-pressure expansivity is
not silently promoted to a complete thermal EOS.

Primary review also consolidated duplicate majorite cards and removed two
records whose cited sources do not define the migrated EOS: the Fei-labeled
FeO static BM3 and the Hixson--Fritz tungsten BM3 reduction. The corrected
InN entry follows Muñoz and Kunc's theoretical Murnaghan fit. The Campbell
B2-KCl entry is explicitly labeled as a Campbell-ratio/Dewaele-B1-volume
composite, including the limited uncertainty propagation that combination
permits.

All eleven Sokolova marker records distinguish scientific-fit provenance from
software lineage. Sokolova et al. (2013), Tables 1 and 4, supply the reference
inputs and final cross-calibrated coefficients. Dorogokupets et al. (2012) is
the preceding fit source for diamond and the nine metals; Dorogokupets (2010)
is the earlier MgO source. Sokolova et al. (2016) supplies the Excel/VBA
implementation, conventions, corrected equations, and the implemented MgO
anharmonic-coefficient correction. Each `.eosmat` record stores these roles in
`source_lineage`. Its identifier uses `_sokolova_2013` for the scientific fit
year; the former workbook-year `_sokolova_2016` identifier is not retained.

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
The ten Shen--Smith (2026) Cu-anchored 300 K Vinet fits are validated directly
against Equation (4), Tables I--II, and the phase-specific range discussion in
the supplied version of record. Their printed `K0` and `K0'` errors are stored,
but no confidence level or covariance is inferred because the paper states
neither.

The Dioptas-facing thermal type `AlphaKT` is preserved for interchange, while
its canonical mechanism identifier is `thermal_reference_state`. Peritheos
evaluates validated instances with `ThermalReferenceStateEOS` rather than
making the source/application label part of the public equation name.
The corrected Anderson Au record instead uses the format-3 extension
`LogVolumeThermalPressure` / `log_volume_thermal_pressure`. Older consumers
must preserve this unknown component rather than evaluating it as `AlphaKT`;
the latter is a different equation.

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
An omitted `thermal_expansion_law` means `constant`; the explicit
`linear_temperature` value requires `alpha1` and integrates
`alpha(T)=alpha0+alpha1*T`. An omitted `reference_volume_law` means
`integrated_expansivity`; `linear_temperature` instead applies the direct
relation `V0(T)=V0(Tr)*[1+alpha0*(T-Tr)]` used by the staged aragonite record.
`berman` applies EosFit7's truncated quadratic
`V0(T)=V0(Tr)*[1+alpha0*(T-Tr)+0.5*alpha1*(T-Tr)^2]`.
