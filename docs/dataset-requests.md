# Incomplete datasets and author requests

This page tracks whether obtaining additional source data would resolve a
documented reproduction gap, exactly what to request, and the outreach status.
The [paper investigation ledger](paper-investigation-ledger.md#papers-with-unavailable-direct-refits)
provides the complete inventory of catalog papers with unavailable direct
refits and their record-level blockers. Its
[withheld/deferred papers](paper-investigation-ledger.md#withheld-or-deferred-papers)
also include candidates awaiting primary evidence.

The list below contains reviewed requests. Papers absent from it have **not
yet been assessed for author outreach**; absence does not mean their data are
complete or contacting the authors would be unhelpful. Some ledger gaps concern
computed grids, adopted coefficients, or model definitions rather than missing
experimental observations.

## Reviewed requests

| Paper / material | Recommendation | What additional data would enable | Outreach status | Last reviewed |
|---|---|---|---|---|
| [Shen and Smith (2026), MgO and NaCl-B2](#shen-and-smith-2026-mgo-and-nacl-b2) | Exact author-used inputs needed to resolve both fit holds | Establish whether additional observations or unpublished processing entered Table II; replay Cu reduction and regression | Not contacted | 2026-10-17 |
| [Frost et al. (2023), Pd](#frost-et-al-2023-palladium) | Recover original regression inputs/settings to resolve the recommendation hold | Reconcile both published fits with all 93 verified observations and the stated coordinate errors | Not contacted | 2026-10-17 |
| [Morard et al. (2026), FeS](#morard-et-al-2026-fes) | Accepted on parameter agreement; no author email needed | GUI and console agree; residual-envelope discrepancy retained as informational | User supplied author communication, supplement, EOS and selected input | 2026-10-03 |
| [Miozzi et al. (2020), hcp iron](#miozzi-et-al-2020-hcp-iron) | Author clarification would resolve a substantial source-fit discrepancy | Replay the actual pressure reduction and EosFit-7c objective; verify thermal normalization | Not contacted | 2026-10-02 |
| [Chen et al. (2010), fcc argon](literature-reproductions/chen-2010-argon.md) | Original coefficients and data needed to resolve source inconsistency | Reconcile the 2 GPa elastic constants with Figure 5 and the extrapolated density; reproduce integration and fit | Not contacted; user requested no outreach | 2026-09-24 |
| [Crichton et al. (2016), bcc vanadium](#crichton-et-al-2016-bcc-vanadium) | Contact authors; high priority for original-fit verification | Refit the complete 62-state experiment; resolve the two incomplete fixed-K′ alternatives | Not contacted; historical contact route available | 2026-09-19 |
| [Nisr et al. (2017), hydrous silica](#nisr-et-al-2017-hydrous-silica) | Optional follow-up for exact regression and calibration replay; published coefficients already reproduced | Resolve the original fit objective, weights and Au pressure reduction | Not contacted | 2026-09-19 |
| [Maltby et al. (2024), fcc argon](#maltby-et-al-2024-fcc-argon) | Contact authors for sample-state and neighbour-table clarification | Reproduce Table 8 and the original multiproperty fit before promoting a replacement | Not contacted; precise request drafted | 2026-09-24 |

### October 2026 additions: optional workflow recovery

The following requests do not block use of the accepted equations. They concern
exact original fit or calibration replay. No new author outreach is authorized.

| Source | Evidence already recovered | Remaining optional evidence | Contact status |
|---|---|---|---|
| Katsura (2022) forsterite | Author workbook and MATLAB fitting/geotherm code; both published coefficient vectors; 51 observations | Original six-column `Olivine_PVT2021.xlsx`, exact pressure-branch selection, temperature errors and bootstrap settings; script-specific covariance/errors if available | Not contacted |
| Katsura (2009) forsterite | Final article and all 52 Table 1 rows | Complete earlier thermal-expansion inputs and combined-fit selection/normalization; this is distinct from the accepted 2022 records | Not contacted |
| Zhao (1998) diopside | Complete source table, published equation/coefficients, plot-corroborated row-15 error and separately corrected reusable dataset | Exact original selected fit input, unrounded row-15 pressure, V0 treatment, weighting/covariance and paired NaCl readings | Not contacted |
| [Crichton–Ross (2000) phase E](literature-reproductions/phase-e-crichton-ross-2000.md) | All eleven primary P,a,c,V rows/errors, author-thesis comparison, published BM3 and parity within reported uncertainty, coefficient differences and residual plot | Optional original unrounded fit table/settings and covariance; paired quartz volumes/errors for independent pressure recalculation. These inputs do not block use or parity of the accepted EOS; no further outreach is needed for integration | Not contacted |
| Farla (2025) rubidium halides | Primary numerical supplement, historical CsCl implementation, independent weighted fits and direct EosFit7c replay | Original saved EosFit7/GSAS-II projects, exact row flags, unrounded input/errors/correlations, raw ruby wavelengths and saved CsCl runtime settings | Not contacted; user excludes author contact |

See the [forsterite audit](literature-reproductions/forsterite-katsura-2009-2022.md),
[diopside audit](literature-reproductions/diopside.md) and
[rubidium-halide audit](literature-reproductions/rubidium-halides-farla-2025.md).
The diopside parameter comparisons have conditional parity within reported
error widths; rubidium fits are similar and do not reproduce every coefficient
at publication rounding. Running EosFit independently is not recovery of the
original author workflow or original saved projects.

### Shen and Smith (2026): MgO and NaCl-B2

Paper: [Simultaneous x-ray diffraction measurements of nine pressure calibrants to 140 GPa](https://doi.org/10.1103/fxgq-96sg).
Records: `mgo_shen_2026_vinet_3`, `nacl_b2_shen_2026_vinet_2`.
Evidence: [source-input and fixed-volume fit audit](literature-reproductions/shen-smith-2026-cu-refits.md).

The bundled Table S1 transcription accounts for all 194 DAC-2 MgO steps
and 131 NaCl-B2 observations, with simultaneous Cu volumes/errors. Both
materials were absent from DAC-1. The described Table II fit uses this
study's direct Cu-referenced P-V data and presents literature for comparison;
no additional experimental fit series is identified. Nevertheless, the data
availability statement mentions other author-held data, so completeness of
the actual author inputs cannot be claimed. The exact NaCl coexistence mask
also remains ambiguous.

Both records are `not_reproduced`, not recommended for quantitative use and
available only by explicit diagnostic selection. Minimum useful author inputs:
the exact unrounded MgO/NaCl-B2 fit tables with run IDs and exclusions, any
additional historical or unpublished experimental rows used, Cu pressure
reduction code/constants, fitting coordinate/weights and confirmation that
Table II reports that fit. Also request the EOS-fitting software name/version
(EosFit or otherwise), original saved input/project files and commands, or the
actual fitting script. Neither the article nor the inspected figures/stress
supplement identifies the EOS regression software. Covariance and uncertainty conventions would
additionally permit a justified statistical comparison. No outreach has been
made or authorized.

### Frost et al. (2023): Palladium

Paper: [The equations of state of statically compressed palladium and rhodium](https://doi.org/10.1063/5.0161038).
Records: `palladium_frost_2023_vinet_1`, `palladium_frost_2023_bm3_2`.
Evidence: [source replay and acceptance decision](literature-reproductions/palladium-frost-2023-audit.md).

All 93 official supplementary Pd observations and their printed lattice
errors are already recovered and verified. Independent orthogonal and
effective-variance interpretations of the stated uncertainty-aware fit give
different central coefficients. The combined-two-width screen fails with
absolute-error covariance but passes with residual-scaled covariance. Source
error confidence and covariance treatment are unspecified; both fits remain
numerically similar. Both records are `not_reproduced` and not recommended,
and automatic Pd equilibrium default selection is on hold.

Minimum useful author inputs are the original unrounded fit table, row/run
manifest and ambient-point treatment, exact residual objective and error
transformation, fixed parameters or priors, uncertainty confidence/rescaling,
and saved fit project or script. These would distinguish an undocumented
workflow detail from a source-side inconsistency. Paired W/Au volumes and
their calibration settings are additionally useful for independent pressure
reduction, but are a separate verification layer. No outreach has been made
or authorized by this audit.

### Morard et al. (2026): FeS

Paper: [FeS phase diagram and thermal Equation of State under high pressure](https://doi.org/10.1103/h4pj-rvxx).
Record: `fes_vi_morard_2026_bm3_mgd` (accepted on parameter agreement, 2026-10-19).
Evidence: [source audit](literature-reproductions/morard-2026-fes.md).

All 146 supplied Table S1 observations are bundled. Standard integrated MGD
replay gives pressure RMS 1.911 GPa and a maximum absolute residual 6.149 GPa;
15 rows exceed the article's ±3 GPa bound, including 13 above 40 GPa.
The variable-exponent alternative also fails. The public Figshare PT workbook
and supplied Table S1 have different pressures and remain separate sources.
Exact GSAS lattice matching recovers all final-table KCl markers; an independent
calibration replay reproduces every supplied pressure within 1.3×10⁻¹³ GPa.
Table S1 C is reciprocal to the printed definition, and its dV uses linear
lattice-error addition. Those values are preserved unchanged.

On 2026-10-03 the user supplied the figures supplement (including S5 but no
Table S2), `FeS6.eos`, and both combined author input variants. These recover
21 literature cold rows. User-reported author personal communication confirms
that the 11 unsuitable quenched rows were excluded. Three additional typed datasets
preserve the selected 167-row input, its literature cold subset and the original
Sata table errors; the 11 rows and the 178-row variant are not database datasets. The saved March 2024 EOS confirms full MGD; its
reference temperature and several coefficients differ from the article.
Actual EosFit replay has thermal RMS 1.632 GPa, maximum 5.058 GPa and 14 rows
outside ±3 GPa. See the [author-input follow-up](literature-reproductions/fes-eosfit-direct.md#author-input-follow-up-2026-10-03).

The user subsequently confirmed that the supplied files are all inputs used by
Guillaume Morard. They are therefore treated as the complete author-used inputs;
a new request for the same data/EOS files is not needed. The distributed
20211014 GUI comparison is complete and agrees closely with the console.
On 2026-10-19 the maintainer accepted staged parameter agreement within reported
uncertainties as sufficient for catalog use. The ±3 GPa envelope discrepancy
is informational; exact original workflow remains unverified. No email is
needed on this criterion and no author outreach was performed by this audit.

### Crichton et al. (2016): bcc vanadium

Paper: *High-temperature equation of state of vanadium*,
[DOI 10.1080/08957959.2015.1123256](https://doi.org/10.1080/08957959.2015.1123256).
Record: `vanadium_bcc_crichton_2016_bm3_thermal`.
Evidence: [source audit and reproduction](literature-reproductions/crichton-2016-vanadium.md).

**Why contact is worthwhile.** The published equation reproduces the plotted
curves, but only 29 distinguishable experimental symbols could be digitized,
and their exact temperatures are missing. The nominal-temperature proxy does
not recover all published thermal coefficients. Obtaining the complete input
table would make a source-faithful refit possible; it would not guarantee
coefficient parity. The existing published-coefficient record remains usable
within its documented range.

**Minimum request for the main fit:**

- The original 62 P–V–T rows in their original units and volume normalization,
  with exact row-wise temperatures, run identifiers, and the selected/excluded
  points, including treatment of the initial stressed-foil cycle.
- The EOS-FIT5.2 input/output files or equivalent fit settings: residual
  definition, weighting, fixed/free parameters, and reference-state convention.
- Measurement uncertainties and parameter covariance, if retained. These would
  support statistical checks; their absence should not prevent sharing the
  central observations.

**Useful additional requests:**

- Paired NaCl/Au lattice measurements, the precise calibrant EOS references and
  coefficients, and PTX-Cal settings, to replay the original P–T reduction.
- Complete fitted parameter sets for the K′ = 3.5 and K′ = 4 alternatives,
  especially V0 and the thermal coefficients, to assess the withheld records.
- Clarification of the Conclusion's thermal-expansion percentages, which differ
  from calculations using the printed coefficients.
- A public repository deposit or permission to redistribute the supplied data
  with attribution, so the reproduction can be bundled and independently run.

**Access checked.** The UCL accepted manuscript and publisher abstract were
inspected. Publisher full text/supplement routes, institutional repositories,
correction searches and the local Zotero EOS collection did not yield the
complete numerical dataset. The final typeset full text was inaccessible.
See the audit for the exact routes and limitations.

**Contact route and history.** The accepted manuscript lists W. A. Crichton as
the corresponding author with `crichton@esrf.fr`. This is the paper's historical
address; current deliverability has not been verified. No request has been
sent and no response or permission has been received.

### Nisr et al. (2017): hydrous silica

Paper: *Phase transition and equation of state of dense hydrous silica up to
63 GPa*, [DOI 10.1002/2017JB014055](https://doi.org/10.1002/2017JB014055).
Records: `stishovite_nisr_2017_dry_bm3`,
`hydrous_stishovite_nisr_2017_bm3`, and
`hydrous_silica_cacl2_nisr_2017_bm2`.
Evidence: [source audit and reproduction](literature-reproductions/nisr-2017-hydrous-silica.md).

**Why follow-up could help.** All 43 official supplement observations are
bundled, and independent fits recover the published free coefficients within
one reported sigma. Exact replay remains limited by unspecified regression
weights and missing calibrant measurements. These gaps do not withhold the
three published EOS records.

**Minimum request for exact regression replay:** the original fit input/output
files or equivalent residual definition, weights, row selections, treatment
of printed zero volume errors, and fixed/free parameter settings. Parameter
covariance and clarification of the conflicting ambient volumes would improve
the uncertainty and reference-state audit.

**Optional calibration inputs:** paired Au lattice measurements and the exact
Fei et al. (2007) equation variant and coefficients; numerical water-volume
calibration coefficients and their uncertainty propagation would separately
enable replay of the inferred water content. These requests concern different
reductions and are not prerequisites for the existing coefficient comparison.

**Access and contact status.** The corrected publication of record and official
supporting information were inspected; their versions and checksums are in the
audit. No author has been contacted, no current contact route has been verified,
and no additional data or sharing permission has been received. This optional
request remains open.

### Maltby et al. (2024): fcc argon

Paper: *Equation of State for Solid Argon Valid for Temperatures up to 300 K
and Pressures up to 16 GPa*, [DOI 10.1063/5.0237497](https://doi.org/10.1063/5.0237497).
Study: `argon_fcc_maltby_2024`; no validated executable EOSMAT record.
Evidence: [source audit](literature-reproductions/argon-maltby-2024.md) and
[rounding/constrained-refit investigation](literature-reproductions/argon-maltby-2024-refit.md).

**Why contact is worthwhile.** The printed equations and geometric fcc sum do
not reproduce Table 8's sample volume. Tested coefficient-rounding intervals,
cutoffs and derivative alternatives do not resolve the difference. Official
Table SI.1 additionally includes an impossible fcc shell at squared distance
14. A constrained pressure refit improves grouped errors on a screened subset
but worsens low-pressure and external checks; it is not a replacement for
understanding the published implementation.

**Minimum request:** full-precision final coefficients; numerical cutoff and
neighbour list; executable Table 8 calculation, including explicit potential
density dependence and tail differentiation; clarification of SI.1 shell 14.

**For the complete refit:** exact selected primary data, especially the 38
Dewaele rows, plus caloric, expansivity and compressibility observations;
coordinate uncertainties, weights, fixed/free coefficients, staging and
covariance where available. Request attribution/redistribution terms or a
public deposit for supplied material. The unreported complete objective
cannot be reconstructed from the currently bundled compression data alone.

**Access and contact status.** The NVA accepted manuscript and official AIP
supplement were checked, including rendered equations and the neighbour table.
Title/DOI/correction searches, GitHub repository searches and accessible
ThermoPack main source did not yield an author implementation or correction.
The publisher page lists corresponding author Øivind Wilhelmsen at
`oivind.wilhelmsen@ntnu.no`. A [specific draft request](author-requests/maltby-2024.md)
is prepared. No message has been sent, and no additional data or permissions
have been received. This request remains open.

## Brucite: unresolved source questions (no active research follow-up)

**Status: not contacted; no request or further fitting planned in this task.**
The [audit](literature-reproductions/brucite.md) records failure to establish
parity for the full published thermal EOS, which is not recommended. Only
three individually source-checked isotherms remain executable; all refits are
historical diagnostics. The following questions preserve the evidence needed
if an authoritative correction or author clarification becomes available.
They are not a requirement to conduct new research within Peritheos.

- **Fukui et al. (2003), [DOI](https://doi.org/10.1007/s00269-003-0353-z):**
  full-precision V0 used for each isotherm; the static fit’s fixed/free V0 choice;
  original objective, weights and covariance; independently regressed thermal
  intercepts; clarification of the 873 K normalization mismatch (graphical V0 about 42.70 A^3
  versus printed-expansion 42.516 A^3); paired Au lattice
  readings and row temperature information. Final primary article and all
  48 printed observations are already recovered. No missing marker reading is
  replaced by an inverse-EOS calculation.
- **Xia et al. (1998), [DOI](https://doi.org/10.2138/am-1998-0107):** final
  combined-fit V0 and dKprime/dT, omitted from Table 7, and the exact adopted
  Fei–Mao ambient data. All 54 Table 6 rows are already preserved.
- **Fei and Mao (1993), [DOI](https://doi.org/10.1029/93JB00701):** accessible
  final full article and authoritative reference volume/data before any
  executable alternative can be accepted. The publisher abstract alone is
  insufficient; publisher and institutional/public-repository searches have
  not recovered the full paper.
- **Wu et al. (2026), [DOI](https://doi.org/10.3390/molecules31101631):**
  clarification of Section 4’s 1.8/3.2 GPa and 573/823 K assignments versus its
  negative-slope fitted boundary. Printed values are retained without swapping.

The primary ambient structural gap is resolved through He (2023), Table 4.3;
the old Zigan/Rothbauer article is no longer required to establish the catalog
structure. None of these requests has been sent.

## Maintaining the list

Add one entry per paper when an audit identifies a concrete information gap.
Link the affected records and audit, distinguish essential inputs from optional
improvements, record the access routes already checked, and explain what a
response would unlock. Recommend outreach when the requested material could
resolve that gap. A coefficient mismatch alone is not evidence of missing data;
first inspect the model, units, fit selection and objective.

Use explicit outreach states: **not contacted**, **contacted** (with date),
**response received**, **data received—validation pending**, **resolved**, or
**closed without data** (with reason). Record follow-up dates and sharing terms
when known. Receipt of a file does not establish parity: preserve its provenance,
rerun the reproduction, and update the scientific ledgers before marking the
gap resolved. Keep unresolved secondary requests visible if only part of the
requested data arrives.

### Miozzi et al. (2020): hcp iron

Paper: *A New Reference for the Thermal Equation of State of Iron*,
[DOI 10.3390/min10020100](https://doi.org/10.3390/min10020100).
Evidence: [source audit and reproduction](literature-reproductions/miozzi-2020-iron.md).

All 131 published supplementary rows, reported coordinate errors and paired MgO
volumes are bundled. The preferred thermal coefficients give 6.98 GPa RMS
pressure residual under standard MGD normalization, exceeding the stated
−3 to +3 GPa interval; both printed room-temperature fits also miss the supplied
MgO-series pressures. Three records retain explicit-selection not-reproduced
status. This finding does not establish a mistake in the authors' measurements.

**Minimum clarification:** the actual EosFit-7c input files with original row
selection, the exact Speziale MgO pressure coefficients/implementation, the
intended energy/molar-volume normalization, and the preferred thermal reference
volume (22.80 versus 22.81 A³ versus 6.87 cm³/mol).

**Additional useful inputs:** original regression weights and error handling,
parameter covariance/confidence definitions, and paired ruby/Sr gauge readings
for the helium series. The nonpreferred joint fits additionally need their
Dewaele (2006) row selection. No author contact has been made.

### Morard GUI follow-up, 2026-10-19

The [actual distributed-GUI comparison](literature-reproductions/fes-gui-replay.md) has been completed for saved-cold and cold-then-hot controls on the complete confirmed selected input. The GUI agrees closely with the console and retains residuals above ±3 GPa. The exact 20210609 executable/workflow remains unverified. The author-used input files are complete; do not request them again. The maintainer accepts the parameter agreement for catalog use; residual-envelope non-parity is informational and no email is needed on this criterion.
