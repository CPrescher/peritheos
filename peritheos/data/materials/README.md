# Bundled material library

This directory contains 225 curated material documents with 614 EOS records. The
collection began with the 120-material, 147-EOS-record Dioptas 0.10.0 database,
tag commit
`5a8bfd81d10bfab3499039603380aae34576d60a`. Its project source is
<https://github.com/Dioptas/Dioptas>.

Every EOS record has an explicit `determination_method` and a rationale under
`parameter_provenance.determination_method`. The 2026-09-23 review classifies
473 records as experimental, 114 as theoretical, 20 as hybrid, and seven as
unknown. These counts include the two deferred Campbell source records.
The seven unknowns are the Holland et al. (2013) THERMOCALC modified-Tait
endmembers: the existing source audit verifies the coefficients but does not
establish each endmember's measurement-versus-calculation constraints.
They remain explicitly unknown rather than inferring their origin from the
paper title or `self_consistent` calibration method.

Classification applies to the represented parameter set, including fixed
coefficients and thermal terms. Thus Ono's experimental CaSiO3 reference plus
AIMD thermal term, Oganov's experimentally shifted MgO curves, and Holmes's
calculated Pt curve with a measured reference volume are hybrid. The
unshifted calculated MgO curves remain theoretical. A measured compression
fit retains experimental classification when using a theoretical pressure
calibrant; calibrant ancestry is recorded separately. See the
[schema definition](../../../docs/eosmat-schema.md#eos-record).

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
dated 2026-09-20 classifies 612 of 614 bundled records as
`primary_source_validated`. The Campbell fcc-Fe and FeO records are retained
as `deferred` source evidence because of unresolved coefficient-refit
discrepancies; they are not exposed by the executable catalog. The complete
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
library also includes the Datchi et al. absolute-zero c-BN MGD record, the
Suzuki epsilon-FeOOH reference-temperature thermal EOS, and the Noguchi et al.
700 K CaSiO3 BM2-MGD record, plus Nisr et al.'s dry stishovite, hydrous
stishovite, and hydrous CaCl2-type silica reference-temperature EOS records.
The Katsura et al. ringwoodite thermal EOS remains
excluded because its published atom-count normalization could not be
reproduced consistently. The KCl record is the preferred
`kcl.eosmat` record and
keeps its measured 298 K range distinct from its molecular-dynamics thermal
extension. The Martinez staged result uses its exact Equation (3) direct-linear
reference-volume law. The paper's separate global thermal BM3 entry is excluded
because its fitted reference volume is omitted and the remaining coefficients
do not reproduce the printed dataset under the documented equations.

The Chidester et al. (2018) thorianite and cotunnite-type ThO2 cards retain
separate experimental thermal fits and separately sourced experimental
structures. All 100 official supplemental rows are bundled, including the
12 room-temperature fluorite observations excluded near the volume anomaly.
Cotunnite's 300 K reference curve is extrapolated from high-temperature data.

## Space-group completeness

The 2026-09-22 review checked all 225 material documents, including names,
phase labels, and existing structure notes. Three missing assignments were
recoverable from an explicit elemental hcp phase: `beryllium_hcp`, `lead_hcp`,
and `zinc_hcp`. These now carry `symmetry: HEXAGONAL`, `space_group: P63/mmc`,
and `space_group_number: 194`, derived from the
[A3 hcp prototype](https://aflow.org/p/A_hP2_194_c-001/).
Each card records the phase citation and prototype mapping under
`source.space_group_assignment`. This does not supply an experimental lattice
or refined atomic coordinates.

The remaining nine omissions are intentional:

| Material identifier | Reason no unique space group is assigned |
| --- | --- |
| `fesio3_liquid`, `mgo_liquid`, `mgsio3_liquid` | Liquids have no crystallographic space group. |
| `ca0988mg0918fe0078mn0016c2o6_dolomite_iii` | The source supplies a possible monoclinic indexing cell, with unresolved structure and alternative indexings. |
| `casio3_perovskite_tetragonal` | The source's P4/mmm pseudocell is a volume/indexing convention; the diffraction data leave the exact tetragonal space group open. |
| `cu_handbook_1972_legacy` | The historical derived isotherm has unresolved crystallographic phase identity. |
| `fe093o_rhombohedral`, `mgfe94o_rhombohedral` | Rhombohedral distortion is identified, but the sources do not assign a unique space group. |
| `fesio3_bridgmanite` | The source identifies orthorhombic Pv and its cell contents without specifying a space-group setting. |

Space-group metadata is curated in the material documents; the loader does
not guess from names. A crystal-system label alone, an ambiguous indexing cell,
or a space group mentioned as an alternative is insufficient for an assignment.

The [tetragonal CaSiO3 source review](../../../docs/literature-reproductions/casio3-tetragonal-space-group.md)
distinguishes the Shim (2002) P4/mmm indexing model, now explicit in
`source.indexing_model`, from the later I4/mcm material card. It records the
original papers' evidence and qualifications for both interpretations.

## Primary observation tables

When a reviewed primary paper prints a recoverable observation table, the
material document links that table through its top-level `datasets` array.
Litasov et al. (2007) adds seven experimental LT superhydrous-phase-B records,
69 complete P–V–T rows with both Au pressure scales and marker observations,
and 30 diffraction peaks. The independent Pnn2 structure applies the authors'
O6 correction and official 2014 H2 erratum; the original archive CIF is retained.
The source audit documents equation typos, calibration qualifications and held
branches, separately from the Al-bearing Xu (2024) investigation.

The current bundle contains 312 distinct primary datasets with
22,186 observation rows, represented by
351 material-document links to
472 EOS records.
This inventory includes source-reported derived grids and diagnostics; dataset
metadata distinguishes them from independent observations. For Zha et al.
(2004), eight paired Au/Re observations constrain the continuous thermal EOS,
while Tables IV/V contain derived isotherms and Table VII diagnoses strain.
The Ono et al. cubic-SnO2 table is linked from both legacy diffraction-pattern
entries. The Shen--Smith Table S1 workbook is linked across its ten calibrant
and phase entries while retaining the simultaneously measured Cu reference
volumes even though the library has no separate Shen--Smith Cu EOS record.
Small tables may remain inline; larger tables are stored as
SHA-256-checksummed CSV resources under `../datasets/`. Quantities and
uncertainties retain the paper's reported conventions and units, including
formula-unit volumes, molar volumes, densities, lattice parameters, and
pressure-calibrant readings.

Dataset-specific reuse terms are recorded in metadata. An unspecified source
data license is recorded as unspecified, rather than used as an automatic
exclusion or relabeled CC0. Scoped CC0 applies only to contributors' own rights
in factual transcriptions, normalization, metadata and arrangement; explicit
source terms and attribution remain intact. Noguchi (2013) Table 1 is bundled
with all 54 rows, both pressure scales, errors and the three fit exclusions. For the Dorfman et al.
co-compression table, CC0 applies only to rights held by Peritheos contributors
in the factual CSV transcription and arrangement; it does not relicense the
publisher's article, auxiliary PDF, or third-party rights.

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

Only validated records are executable by default. Deferred evidence and
source-only material cards remain in the files so Dioptas and other consumers
can preserve the catalog without silently converting inherited values into
Peritheos-endorsed pressure scales.

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

Three Dioptas structure-only entries (`fes_iii`, `nitrogen_epsilon`,
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
