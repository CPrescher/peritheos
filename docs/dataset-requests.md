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
| [Miozzi et al. (2020), hcp iron](#miozzi-et-al-2020-hcp-iron) | Author clarification would resolve a substantial source-fit discrepancy | Replay the actual pressure reduction and EosFit-7c objective; verify thermal normalization | Not contacted | 2026-10-02 |
| [Chen et al. (2010), fcc argon](literature-reproductions/chen-2010-argon.md) | Original coefficients and data needed to resolve source inconsistency | Reconcile the 2 GPa elastic constants with Figure 5 and the extrapolated density; reproduce integration and fit | Not contacted; user requested no outreach | 2026-09-24 |
| [Crichton et al. (2016), bcc vanadium](#crichton-et-al-2016-bcc-vanadium) | Contact authors; high priority for original-fit verification | Refit the complete 62-state experiment; resolve the two incomplete fixed-K′ alternatives | Not contacted; historical contact route available | 2026-09-19 |
| [Nisr et al. (2017), hydrous silica](#nisr-et-al-2017-hydrous-silica) | Optional follow-up for exact regression and calibration replay; published coefficients already reproduced | Resolve the original fit objective, weights and Au pressure reduction | Not contacted | 2026-09-19 |
| [Maltby et al. (2024), fcc argon](#maltby-et-al-2024-fcc-argon) | Contact authors for sample-state and neighbour-table clarification | Reproduce Table 8 and the original multiproperty fit before promoting a replacement | Not contacted; precise request drafted | 2026-09-24 |

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
