# Candidate material equations of state

Work-status review against the local catalog and source audits: 2026-11-01.
Original candidate discovery pass: 2026-09-02.

This is a review backlog, not a set of validated Peritheos records. A candidate
is included when a primary paper reports a pressure-volume or
pressure-volume-temperature relation for a homogeneous solid phase, the phase
is absent from the bundled library or materially extends an existing record,
and the published formulation is plausibly representable in Peritheos.

Every accepted candidate must also be diffraction-ready for PhaseSmith and
Dioptas. Its material file must contain the phase-specific crystal system,
conventional unit cell, formula-unit count, space group, and occupied atomic
sites (or a documented fallback peak list when a defensible atomic model is
not published). Crystal structure and EOS provenance are audited separately;
an EOS is not accepted into the library with only a phase name or nominal
formula.

The [generated coverage overview](catalog-coverage.md) is the current inventory
and records which materials have static, thermal, and Hugoniot coverage. Its
[machine-readable snapshot](data/catalog-coverage.json) also identifies differences
between the live material files and the scientific ledgers. It includes local
work and must not be read as the contents of a published release.

Current material/record counts and scientific-ledger reconciliation are maintained
in the generated coverage overview rather than duplicated in this backlog.
Forsterite, diopside, the rubidium halides, phase E, and the scoped brucite
addition are integrated. Finger (1979) Fe-ringwoodite static is integrated at
296.15 K and 0–4 GPa; its separate 2011 thermal candidate remains staged.
Brucite does not yet provide an accepted broad high-pressure thermal or
disordered-phase EOS; its precise limits are in the source audit.

## How work status is assigned

An integrated record or a documented non-reproduction finding can complete a
curation investigation. Under the [curation stopping rule](adding-materials-and-eos.md#stopping-rule-unresolved-reproduction-discrepancies),
`not_reproduced`, `deferred`, `similar`, or `not_refittable` is not by itself an
instruction to keep fitting. Work status is separate from scientific validation,
fit parity, recommendation, and catalog availability; this review changes none
of those record-level findings.

The active queues below contain additions, integration decisions, or concrete
checks still to perform. Completed investigations retain their limitations and
specific reopening conditions separately. Missing exact author-workflow details
are evidence gaps, not automatically active tasks or author requests. New model
development and outreach require their own scope; a backlog row does not
authorize either.

This review covers the listed backlog and the source audits for the catalog's
currently deferred or `not_reproduced` records. Other numerical-ledger outcomes
retain their individual comparison scopes; their labels alone do not add them
to an active queue.

## Active material-gap queue

### Fe3C expansion, reviewed 1 November 2026

The [cementite investigation](literature-reproductions/cementite.md) reviews
Prescher (2012), Litasov (2013) and Takahashi (2019). The original
[Sata (2010) static BM3](literature-reproductions/sata-2010-fe3c.md) is now
retained as a qualified nondefault source EOS, with conditional parameter parity
and the simple caveat that some quoted point errors are smaller than the overall
experimental scatter. The recalculated-pressure refit remains separate.
[Ono–Mibe (2010)](literature-reproductions/ono-mibe-2010-fe3c.md) is now integrated
as the original experimental BM3 at 300 K and 0–34.3 GPa, preserving all 14
observations and keeping the 0 K calculations and later Litasov comparison separate.
[Hu (2019)](literature-reproductions/hu-2019-fe3c.md) is now integrated as a
full experimental Hugoniot at 85.6–245 GPa (scope updated 2 November), retaining
the porous precursor and all six source shots. Possible melting at higher
pressures remains an explicit phase caveat; the refit and thermal reduction stay separate.
The existing Scott coefficients remain unchanged; their digitized reconstruction
is now qualified similarity. All 81 Litasov Table 1 rows, 60 Takahashi Table S1 rows and
19 user-supplied Prescher Table S1 states are recovered. Both Takahashi thermal
models are now catalog records of `cementite_nonmagnetic`, with the Prescher
47.4 GPa structure and distinct pressure-scale provenance. Prescher's published
full-range EOS is the sole public Prescher choice, with unchanged K coefficients
and errors and adopted Figure 1 full-curve V0=154.5678863704812 Å³. With V0
fixed, all 22 digitized graph markers reproduce both K coefficients and formal
errors at published precision. The complete 19-row table and 22 graph markers
remain separate public datasets. Split fits are audit-only; original confidence,
V0 constraints and author workflow remain unresolved. See the
[integration and comparison](literature-reproductions/prescher-2012-recovery.md#database-integration).
The separate [Litasov (2013) investigation](literature-reproductions/litasov-2013-fe3c.md)
is now complete: KE source calculations are reproduced within printed precision,
the prose HT fit has conditional central-parameter agreement, and body-MGD has
partial agreement; the other alternatives retain quantified limitations or
non-parity. No Litasov record is integrated. Selecting an eligible source variant
and establishing an exact supported model mapping are separate integration
decisions, not a reason to repeat the completed investigation.
[Li (2002)](literature-reproductions/li-2002-fe3c.md) is now integrated with
its original BM3 coefficients, fixed measured V0 and all 14 Table 1 rows.
Source-curve checks and conditional parameter parity pass; original weights,
confidence and ruby/neon pressure-reduction replay remain unresolved.

### Open additions and integration decisions

These priorities supersede the historical A/B/C ordering below. Candidate status
does not establish an executable EOS; the [acceptance requirements](adding-materials-and-eos.md)
still apply. Work is grouped by material or related source family so unrelated
sources, compositions, and pressure standards remain distinguishable.

| Priority | Work item | Current coverage / disposition | Next concrete action |
|---|---|---|---|
| Fe3C | Sata (2010), Ono–Mibe (2010) experiment and Hu (2019) full experimental shock addition completed | Original static records and all observations preserved. Hu’s published six-shot Us–up relation retains the full 85.6–245 GPa experimental path and porous precursor, with a simple caveat about possible melting at higher pressures. The refit and thermal reduction remain separate. | Completed source investigations are not an active refitting queue. The distinct Ono 0 K calculations, Hu thermal reduction and eligible Litasov thermal mapping require separately scoped integration work. Proceed to the next material gap. [Hu audit](literature-reproductions/hu-2019-fe3c.md); [source inventory](literature-reproductions/cementite.md#additional-primary-sources). |
| G05 / C10 | Fe2SiO4 ringwoodite | Finger (1979) derived static ambient-cell reconstruction integrated 1 November; one installed static record, valid only at 296.15 K / 0–4 GPa | Static addition complete. The Armentrout–Kavner (2011) thermal candidate remains separate in the unchanged staged card, not installed; reopen only for its specific source-equation/calibration evidence. [Integration and thermal acceptance conditions](research/ferrous-ringwoodite-armentrout-2011-integration.md). |
| G07 | Garnet thermal coverage | Pyrope, almandine and majorite each have one static record; queued literature search | Find source-owned, composition-specific P–V–T parameterizations and primary observations; do not treat this family as a mixing model. |
| G08 | Akimotoite continuous thermal EOS | Eleven static/temperature-slice records; queued | Review a source-faithful continuous model; separate the Karki temperature slices and partial Dorogokupets synthesis from a complete thermal EOS. |
| G09 / C16 | Phase-D thermal compositions | Existing Mg-rich, Al-bearing and Fe-bearing cards are static; queued | Review Litasov et al. (2008), DOI [10.1029/2007JB004937](https://doi.org/10.1029/2007JB004937); add its measured composition separately if it differs from existing cards. |
| G10 | Feldspars | Albite, anorthite and K-feldspar absent; queued literature search | Start from phase-specific static data and authoritative structures, including Benusa et al. (2005) albite; preserve anomalous compression and polymorph limits. |
| G11 / C18 | Hcp-Zn thermal EOS | Two static Dewaele (2019) records; thermal extension remains open | Audit Errandonea et al. (2018), DOI [10.1088/1361-648X/aacac0](https://doi.org/10.1088/1361-648X/aacac0), including its exact thermal formulation and axial behavior. |
| G12 / C19 | Independent experimental Mo thermal scale | Existing thermal coverage is the Sokolova synthesis; alternate scale open | Audit Huang et al. (2016), DOI [10.1038/srep19923](https://doi.org/10.1038/srep19923), as a separate experimental record. |

## Existing-record verification and maintenance queue

| Work item | Current disposition | Next action / reopening condition |
|---|---|---|
| Coverage and integration maintenance | Routine checks after accepted additions; earlier catalog-count expectations and Rust floating-point JSON round-trip issue are resolved | Run repository maintenance command `python -m scripts.generate_catalog_coverage --check`; regenerate coverage after additions and reconcile both scientific ledgers. Use the generated overview for current totals. |
| Calibration dependencies | See the generated coverage overview for current `reference_eos_not_bundled` and `reference_model_not_supported` totals | Check for stale metadata/links first, then audit genuinely missing standards (including Tsuchiya Au and historical NaCl/Pt formulations). These status counts do not establish that every cited scale is absent. |
| Sakai (2025) nine-material network | All nine published RT records accepted; Au conditional central-value parity integrated in the common ledger | Extend source reconstruction to the remaining sequential calibration stages using measured paired volumes. Keep calculated Table S9 outputs separate. [Audit](literature-reproductions/sakai-2025-nine-materials.md). |
| Morard (2026) FeS-VI | Reviewed 1 November: numerical/source-input audit complete; live conditional parameter assessment is `parity`, while the residual-envelope finding and deferred availability remain. Original 20210609 GUI comparison awaits a located executable; no newer completed comparison found | If the original executable becomes available, compare the confirmed complete selected inputs and saved model, including weighting, groups, cold/thermal selections and refinement controls. Current replacement-release downloads do not establish historical GUI parity. Do not repeat fits or input requests. The subsequently authorized G05/C10 Finger static addition is now integrated separately. [Availability review](literature-reproductions/fes-eosfit-direct.md#current-status-and-original-gui-availability-review-2026-11-01). |
| Uncertainty coverage | See generated coverage for covariance and error-confidence availability; source confidence is often unspecified | Preserve missing errors as missing; distinguish conditional refit covariance, fixed-input uncertainty, and shared calibrant systematics. |
| Continuous spin crossover and liquids | Bounded branch/reference curves exist; general coverage remains limited | Treat continuous crossover and broad liquid thermal models as separately designed model work. |

## Completed investigations with documented limitations

These are completed reviews of the available evidence, not claims of full
author-workflow reproduction or unrestricted scientific use. Their retained
source discrepancies are findings. Reopen only for the specific evidence or
separately scoped integration/model task below; no new outreach is planned by
this review. The linked audits retain the quantitative comparisons and source
provenance.

| Investigation | Completed finding / current catalog disposition | Reopening condition or separate task |
|---|---|---|
| Miozzi (2020) hcp-Fe | Dominant near-factor-of-two thermal-gamma discrepancy explained by author-reported n=2 at per-mole-Fe volume. Physically normalized n=1 Speziale and Tange independent refits are integrated; the original three published records remain `not_reproduced`. [Audit and registered refit](literature-reproductions/miozzi-2020-iron.md#registered-published-n2-record-and-speziale-n1-refit-2026-10-03). | Exact publication replay would need the actual pressure reducer, selected input, weights/error handling and original session. Those residual gaps do not keep the completed curation/refit task active. |
| Wadsleyite, Katsura (2009) | Calibration/model-mapping investigation completed 30 October. All 85 source rows verified; cold-stage K′=4.108 agrees with 4.1(1). Standard n=7 thermal calculation misses both hot observations and source-calculated curves; pressure recalibration alone cannot reconcile this internal source discrepancy. Native/Python and consistent cell conversions agree. The n=5 diagnostic is unphysical and not adopted; existing deferred hold remains. User reports on 1 November that the email to Tomoo Katsura was sent; reply awaited. [Audit](literature-reproductions/wadsleyite-katsura-2009.md). | Original thermal/pressure calculator or saved project resolving normalization, volume arguments, reference subtraction and Figure 2; Inoue ambient joint-fit inputs/weights; measured paired MgO cells for an independent calibration test. No replacement EOS or new outreach is implied. |
| G04 / C12, brucite | Three source-checked isotherms integrated; published 873 K/continuous thermal EOS not reproduced or recommended. Historical refits remain diagnostics. The audit explicitly ends further fitting under the curation scope. [Disposition](literature-reproductions/brucite.md). | Authoritative source correction or original implementation resolving the discrepancy; a replacement thermal model would be a separate research request. |
| Litasov (2013) Fe3C | Independent equation/fit investigation complete on 30 October: KE calculations reproduced within printed precision, prose HT conditional parameter agreement, body-MGD partial agreement, and scoped non-parity/limitations for other alternatives. No catalog addition. [Complete verdicts](literature-reproductions/litasov-2013-fe3c.md). | An explicit integration decision and exact supported model mapping for an eligible variant; full workflow recovery additionally needs original masks, constraints, objective and thermochemical inputs. Do not repeat the completed source replay or integrate all alternatives as a batch. |
| Campbell (2009) FeO / fcc-Fe | Combined-data and thermal-correction diagnostics completed; hcp-calibrated Seagle rows dominate the FeO mismatch. Both affected published records remain deferred. [Audit](literature-reproductions/campbell-2009-mismatch.md). | Original combined/recalibrated input table and reduction/fit controls; fcc recalibration alone is insufficient. |
| Fei (2007) NaCl-B2 | Published equation and figure-data replay completed; original coefficient parity unresolved. Separate staged digitized-data refit is integrated with its uncertainty limitations. [Audit and refit](literature-reproductions/fei-2007-nacl-b2.md#available-independent-thermal-refit). | Original unrounded measurements, measured Pt cells, errors, weights and calibration/covariance procedure for exact source replay. This does not resume the paused Pt investigation. |
| Dorfman (2012) Au–Mo–Pt | Coupled-fit, constraint-mode and input-recovery audits completed; six published records retain their non-parity qualifications. AN012's role in the mismatch is unproven. [Audit](literature-reproductions/dorfman-2012-cocompression.md). | Additional primary input/code identifying the actual pairs, masks and constraints. The inaccessible dissertation is a possible source-recovery lead, not established missing fit evidence or an automatic new task. |
| Katsura (2004) Mg-ringwoodite thermal / C05 | Thermal normalization discrepancy quantified; source thermal addition withheld, existing static record preserved. The n=5 diagnostic is not a permissible correction. [Audit](literature-reproductions.md#ringwoodite-katsura-2004). | Authoritative numeric atom count, Debye-temperature law and energy/volume normalization or original code. |
| Yokoo (2009) Au/Pt PVT | Qualified diagnostic reconstructions retained; analytical PVT unresolved. The convention investigation is closed pending author evidence; Au follow-up found no evidence satisfying its reopening conditions. [Evidence checklist](literature-reproductions/yokoo-2009-au-follow-up.md). | Provenance-confirmed pressure code, complete specification or correction addressing the established mismatch. More fitting to calculated outputs is insufficient. |
| Chen (2010) argon | Original-constant reconstruction and separate digitized-data refit completed; source parameterization discrepancy documented. Diagnostic original-constant record remains not recommended. [Audit](literature-reproductions/chen-2010-argon.md). | Original source equation, numerical coefficients/inputs or correction; preserve the user's no-outreach instruction. |
| Maltby (2024) argon | Equation, neighbour-table, cutoff, rounding and constrained-refit investigations completed. Original sample-state mismatch persists; exploratory replacement remains research-only. [Source audit](literature-reproductions/argon-maltby-2024.md), [follow-up](literature-reproductions/argon-maltby-2024-refit.md). | Author implementation and sample-state/neighbour-table clarification, or a separate explicit model-development request. |
| Frost (2023) Pd | Source/protocol audit complete; all 93 observations verified. Coefficient comparison depends on an unspecified covariance convention, so a failed replay of the full author protocol is not established. Recommendation remains on hold. [Audit](literature-reproductions/palladium-frost-2023-audit.md). | Original objective, row/ambient-point treatment, constraints, uncertainty convention and saved fit project; no repeated request for the already recovered table. |
| Shen and Smith (2026) MgO / NaCl-B2 | Source transcription, Cu pressure reduction, fixed-volume fits and source-input/software searches completed. Two published fits remain `not_reproduced` and not recommended; other phase fits retain their scoped conditional parity. [Audit](literature-reproductions/shen-smith-2026-cu-refits.md). | Exact unrounded author inputs, NaCl coexistence mask, any additional processed observations, Cu reducer and actual fitting code/project. Available-input non-parity does not establish failure of the unrecovered author workflow. |
| Baty (2024) experimental Pd | Source equation, all 78 observations and source-fit discrepancy audited; affected experimental record remains `not_reproduced` and not recommended. Separate DFT records are unaffected. [Audit](literature-reproductions.md#palladium-baty-2024). | Author-supported input/constraint/weighting or reduction evidence addressing the experimental-fit mismatch. |
| Fu (2023) CaSiO3 thermal fit | Source/input-study reconstruction, joint-fit attempt and normalization audit completed; published thermal discrepancy retained, with a separately identified independent refit available. [Audit](literature-reproductions/fu-2023-eos.md). | Actual source selection, objective/weights and thermal implementation or correction; no further normalization tuning as routine curation. |
| Sakai (2014) Fe0.9Ni0.1 thermal variants | Five source equations and extrapolation checkpoints independently evaluated; the roughly 7.1 GPa hot-data RMS versus the reported roughly 1.1 GPa is documented. Affected variants remain explicit-selection and not recommended. [Audit](literature-reproductions/iron-fei-input-papers.md). | Original pressure reduction, selected regression inputs and implementation reconciling the hot-data residual claim. This finding is scoped to these five thermal variants. |

## Paused by the user

Fei (2007) Pt remains paused pending an author reply and explicit user resumption.
Its coefficients and hold are preserved. This backlog review does not resume
that investigation. [Disposition](literature-reproductions/fei-2007-platinum.md).

## Completed or partly completed candidates

These items are no longer wholly open material gaps. Their historical source
details remain below, with current status taking precedence.

| Historical item | Reconciled status |
|---|---|
| G01, forsterite thermal EOS | Katsura (2022) integration and equation/data/refit investigations complete for both source variants; full-precision author-script values at theta0=760 K recommended, static default preserved. Exact bootstrap recovery is optional. Original Katsura (2009) remains a distinct source hold needing earlier combined-fit inputs/normalization, not a reason to reopen the accepted 2022 addition. [Audit](literature-reproductions/forsterite-katsura-2009-2022.md). |
| G02, diopside | Thompson (2008) static and Zhao (1998) thermal records integrated with conditional parameter parity. Original tables and the separate corrected row-15 analysis dataset retained. Exact author-workflow recovery is optional, not an acceptance blocker. [Audit](literature-reproductions/diopside.md). |
| G03 / C07, B2 RbCl/RbBr/RbI | Integration and follow-up completed; all three now meet their scoped conditional parameter-parity criteria. Exact publication digits, historical CsCl replay and benchmark qualifications remain. The user deferred author contact; missing original projects are evidence gaps rather than pending requests. [Decision](literature-reproductions/rubidium-halides-farla-2025.md#assessment-and-follow-up-decision-2026-10-16). |
| G06 / C11, dense hydrous phase E | Measured Fe-bearing static BM3 integrated with all 11 rows/errors, conditional parameter parity and phase-specific peak-position fallback. Exact original weighting/covariance and paired quartz readings are optional recovery, not unfinished integration. [Audit](literature-reproductions/phase-e-crichton-ross-2000.md). |
| Sueda (2009), Dewaele (2000) MgO, Wang (2012) stishovite, C01 B4C, C02 c-BN, C03 epsilon-FeOOH, C04 CaSiO3 | Implemented; individual reproduction qualifications remain attached to their records. |
| C06 vanadium | Main thermal record implemented; incomplete fixed-K′ alternatives withheld. |
| C08 siderite | Published static record and separate joint thermal refit implemented; the published thermal alternatives remain withheld in its [audit](literature-reproductions/litasov-2013-siderite.md). |
| C09 superhydrous phase B | Seven records implemented, including five thermal alternatives. [Audit](literature-reproductions/litasov-2007-superhydrous-phase-b.md). |
| C13 delta-AlOOH | Both Vanpeteghem (2002) static alternatives implemented; no thermal record. |
| C14 Fe–Si | Source-specific accepted records and two separate physically normalized B2 refits implemented; published B2 thermal normalization remains unresolved. [Audit](literature-reproductions/fischer-2014-fesi.md). |
| C15 ferropericlase Mg0.75Fe0.25O | High-/low-spin reference branches and a Matsui thermal record exist; a continuous crossover model remains open. |
| C17 dolomite/ankerite | Mao (2011) measured-composition dolomite and dolomite-III records implemented; stoichiometric low-pressure dolomite/ankerite coverage remains open. |
| Ono (2006) CaTi2O4-type MgAl2O4 | Implemented as `mgal2o4_cati2o4_ono_2006_bm2_2`; retained below with provenance. |
| C20 / C21 | Continuous spin-crossover phase H remains model work; K-hollandite remains deferred, not an accepted static-slice addition. |

## Historical candidate notes

The A/B/C labels and source coefficient summaries below preserve the original
review context. They are not the active execution order and do not supersede
subsequent source audits or the current coverage overview.

## Ranking key

- **A — review first:** strong scientific or calibration value, explicit
  parameters, and likely implementable with existing classes.
- **B — useful next:** worthwhile coverage, but composition, crystallography,
  reference-state conversion, or source reconciliation needs more work.
- **C — model work:** scientifically attractive, but faithful implementation
  likely needs a new phase-transition or spin-crossover formulation.

"Low effort" means no new public EOS class is expected. It does not waive the
normal primary-source, units, uncertainty, phase, and numerical reproduction
audit.

## A — review first

### Sueda et al. (2009) — CaFe2O4-type MgAl2O4

**Status: accepted and implemented (2026-09-08), two thermal records.**
DOI [10.1016/j.pepi.2008.07.046](https://doi.org/10.1016/j.pepi.2008.07.046).
The HTBM and BM3-MGD fits extend `mgal2o4_cafe2o4`, sharing the published
300 K static triplet. All 46 Table 1 observations and calibrant data are
bundled; independent thermal refits recover every free coefficient within
its printed error. The [source audit](literature-reproductions/sueda-2009-mgal2o4-cafe2o4.md)
documents equation printing errors, uncertainty limitations, and staged fitting.
This manually investigated DOI has no row in the imported LitCurate ledger.

### Audited discovery lead — MgO/periclase, Dewaele et al. (2000)

**Status: accepted and implemented (2026-09-04) as one record.**

- **Identity correction:** DOI
  [10.1029/1999JB900364](https://doi.org/10.1029/1999JB900364) is Dewaele,
  Fiquet, Andrault, and Hausermann (2000), not Fei et al. Fei (1999) is one of
  the comparison datasets used in the paper's combined room-temperature
  analysis.
- **Authoritative model:** the recommended full P-V-T result is a 300 K BM3
  reference isotherm (`V0 = 74.71 A^3`, `K0 = 161 GPa`, `K0' = 3.94`) plus
  integrated-Gruneisen MGD thermal pressure (`theta0 = 800 K`,
  `gamma0 = 1.45`, `q = 0.8`, `n = 2`).
- **Why one record:** Table 3's Murnaghan, Vinet, logarithmic, and alternate
  constrained BM3 rows are room-temperature extrapolation comparisons. They
  are not separate preferred P-V-T parameterizations and are not promoted to
  duplicate executable records.
- **Primary data and calibration:** all 61 Table 2 rows are bundled, including
  the simultaneous Pt lattice parameters. Pressures use the Jamieson et al.
  (1982) Pt shock-Hugoniot/Debye scale. That exact Pt implementation is not yet
  a bundled reference record, so the calibration is explicit but not silently
  recalculated.
- **Reproduction:** the stored BM3 gives 144.947 GPa at `V/V0 = 0.667`, versus
  the paper's 145 GPa. An errors-in-variables refit to the 41 new heated rows,
  with the source's staged/adopted quantities including `gamma0` fixed, gives
  `q = 0.847(122)`, within combined two-sigma uncertainty of the reported
  value. See the
  [full reproduction](literature-reproductions.md#mgo-dewaele-2000).
### Stishovite: lower-mantle Vinet-MGD

**Status: accepted and implemented (2026-09-04).**

- **Source:** Wang et al. (2012), [P-V-T equation of state of stishovite up to
  mid-lower mantle conditions](https://doi.org/10.1029/2011JB009100).
- **Implemented scope:** one Vinet-MGD record and all 56 Table 1 P-V-T
  observations, including simultaneous Au volumes and printed uncertainties.
  The equally successful BM3-MGD coefficients are preserved as the documented
  alternative; the source does not explicitly prefer either cold curve.
- **Audit result:** the existing asymptotic-power-law MGD model exactly matches
  equations (4)--(8), reproduces two independent Table 3 values, and recovers
  all fitted coefficients within combined two-standard-deviation intervals.
  Exact pressure re-reduction awaits an executable Tsuchiya (2003) Au EOS.

### C01 — Boron carbide B4C: add the thermal branch

**Status: accepted and implemented (2026-09-03).**

- **Change:** extend the existing `b4c` material rather than add a duplicate;
  add both the MGD and Berman thermal parameterizations.
- **Source:** Somayazulu et al. (2023), [P-V-T equation of state of boron
  carbide](https://doi.org/10.1098/rsta.2022.0331).
- **Published scope:** 8--50 GPa and 300--2500 K. The implemented record uses
  the quench-data BM3 (`V0 = 328.4 A^3` fixed, `K0 = 221(2) GPa`,
  `K0' = 3.3(1)`) plus the final article's MGD pair (`gamma0 = 0.8` fixed,
  fitted `q = 2.1`). `Tr = 298 K` and `theta0 = 1425 K` follow the explicit
  official supplement.
- **Mapping:** `BM3` + `MieGruneisenDebye` with the explicit
  `integrated_gruneisen` Debye-temperature law. The paper's generalized Tange
  relation reduces exactly to this model for its fixed `a = 1`.
- **Reproduction:** the official supplement includes 51 complete P-V-T rows
  with uncertainties. At their printed V and T, the final parameters give a
  1.39 GPa pressure RMSE and chi-squared per heated observation of 0.961. A
  Peritheos errors-in-variables refit of the 41 heated rows with `gamma0 = 0.8`
  fixed gives `q = 1.050(519)` and reduced chi-squared 0.767; it does not
  independently recover the article's `q = 2.1`. See the
  [full reproduction record](literature-reproductions.md#c01-boron-carbide-somayazulu-et-al-2023).
  The 289.1 A^3, 2023 K row is reproduced as 40.45 GPa versus 40.4(2.3) GPa
  observed.
- **Berman alternative:** the second record maps the paper's EosFit7c fit to
  `BM3 + ThermalReferenceStateEOS` with EosFit's truncated-quadratic Berman
  volume law. The published `alpha0 = 1.94(16)e-5 K^-1` and
  `dK0/dT = -0.008(3) GPa/K` are retained. Matching EosFit's effective-variance
  objective on the 41 printed heated rows gives `alpha0 = 1.811e-5 K^-1` and
  `dK0/dT = -0.01312 GPa/K`; Peritheos's latent-coordinate fit gives
  `1.781e-5 K^-1` and `-0.01232 GPa/K`. The remaining discrepancy may reflect
  a different internal data revision, unavailable unrounded inputs, or an
  undocumented fit selection. The direct EosFit7c result is available as the
  opt-in `b4c_somayazulu_2023_berman_refit` record; the published record is
  unchanged.
- **Primary data:** all 51 supplementary rows are embedded in `b4c.eosmat` as
  a typed, provenance-bearing dataset and linked to all three EOS records.
- **Source caveats:** the record explicitly preserves five contradictions:
  BM2 versus BM3 wording; 1425 versus 1450 K; the final article's
  `gamma0 = 0.8, q = 2.1` versus exploratory supplement pairs; and nominal
  300 K versus explicit 298 K reference temperature. In addition, the
  article's displayed exponential-integral equation is not EosFit's Berman
  equation despite the table and prose identifying the fit as Berman.
- **Crystal structure:** the material carries an average disordered
  `R-3m` (#166) hexagonal cell with five occupied asymmetric sites. Their
  multiplicities and split occupancies give B36C9 (`Z = 9` for B4C), and the
  measured `a = 5.601 A`, `c = 12.087 A` cell gives 328.383 A^3, matching the
  EOS reference volume. This is sufficient for PhaseSmith to calculate the
  reflection list and intensities.

### C02 — Cubic boron nitride: add the thermal branch

**Status: accepted and implemented (2026-09-05).**

- **Change:** extend the existing `boron_nitride` (c-BN) material.
- **Source:** Datchi et al. (2007), [Equation of state of cubic boron nitride at
  high pressures and temperatures](https://doi.org/10.1103/PhysRevB.75.214104).
- **Published scope:** 295 K compression to 160 GPa and 500--900 K data to
  80 GPa. The paper gives a Vinet + Debye-Gruneisen parameter set, including
  `theta0 = 1700 K`, `gamma0 = 1.04(2)`, and `q = 4(1.5)`.
- **Resolved mapping:** a distinct 0 K Vinet curve (`V0,0 = 5.9026(4)
  A^3/atom`, `B0,0 = 397(2) GPa`, `B0,0' = 3.62(5)`) plus integrated-
  Gruneisen MGD with `thermal_pressure_reference = absolute_zero`. The existing
  295 K Vinet record is retained separately and is not reused as the cold
  curve.
- **Reproduction:** the published coefficients reproduce Table VI's
  `V(P=0,300 K) = 5.9055 A^3/atom` to the printed last digit and give a 0.584
  GPa pressure RMSE on all 66 Table IV states, matching the reported 0.6 GPa.
  A q-only unweighted pressure-residual refit gives `q = 7.35(1.89)`, less than
  two combined standard deviations from `4(1.5)`; row rounding and source
  weights are not published, so this diagnostic does not replace Table V.
- **Scope and provenance:** Table IV is bundled in full. The record documents
  zinc-blende `F-43m` c-BN, the per-atom-to-eight-atom-cell conversion, H2005
  ruby calibration, modified SrB4O7:Sm2+ high-temperature scale, and the 4He
  EOS route above 100 GPa. Exact observation-level pressure recalculation is
  blocked by absent marker shifts and helium volumes. The printed table states
  no explicit data-reuse license.

### C03 — Epsilon-FeOOH: replace/augment the static record with P-V-T

**Status: accepted and implemented (2026-09-05).**

- **Change:** added a thermal record to the existing low-pressure P21nm
  `e_feooh` material.
- **Source:** Suzuki (2016), [Pressure-volume-temperature equation of state of
  epsilon-FeOOH to 11 GPa and 700 K](https://doi.org/10.2465/jmps.160719c).
- **Published scope:** 33 Table 1 observations from 0.0001 to 11.07 GPa and
  300--700 K; `V0,300 = 66.278(6) A^3`,
  `K0 = 135(3) GPa`, `K0' = 6.1(9)`, `dK/dT = -0.05(2) GPa/K`, and
  `alpha(T) = 2.6(7)e-5 + 1.0(3)e-7 (T - 300) K^-1`.
- **Exact mapping:** `BM3` + `ThermalReferenceStateEOS`, with `Tr = 300 K`,
  linearly temperature-shifted `K0`, and exponential integration of the
  source's reference-temperature linear expansivity. A new
  `linear_reference_temperature` configuration preserves the paper's
  `alpha0 + alpha1*(T-Tr)` convention without translating its published
  `alpha0` into an absolute-temperature intercept.
- **Reproduction:** the stored published curve has 0.10202 GPa pressure RMSE
  and 0.24321 GPa maximum absolute residual on all 33 printed rows. A joint
  Peritheos pressure-volume errors-in-variables refit with `V0,300` fixed
  recovers all five free coefficients within their published one-standard-
  deviation intervals. The precise EosFit weights, unrounded inputs,
  covariance, and scalar fit statistic are not published, so no replacement
  refit record is added. See the
  [dedicated reproduction](epsilon-feooh-suzuki-2016.md).
- **Data and calibration:** all Table 1 lattice/volume observations and printed
  standard deviations are bundled. Pressures use Brown's (1999) NaCl EOS;
  row-wise NaCl volumes are not reported, preventing pressure recalculation.
  The J-STAGE article exposes no supplement or open data license, which is
  recorded explicitly with the transcription.
- **Phase boundary:** this remains epsilon-FeOOH, not starting-material
  goethite, and remains distinct from the separate hydrogen-centered Pnnm
  low-spin material and the approximately 53 GPa spin-transition literature.

### C04 — CaSiO3 perovskite: add a lower-mantle thermal record

**Status: accepted and implemented (2026-09-05).**

- **Change:** extend the existing `ca_perovskite` material.
- **Source:** Noguchi et al. (2013), [High-temperature compression experiments
  of CaSiO3 perovskite to lowermost mantle conditions and its thermal equation
  of state](https://doi.org/10.1007/s00269-012-0549-1).
- **Published scope:** to 127 GPa and 2300 K. The reference isotherm is at
  700 K (`BM2`, `V0 = 46.5(1) A^3`, `K0 = 207(4) GPa`) with an MGD correction
  (`theta0 = 1300(500) K`, `gamma0 = 2.7(3)`, `q = 1.2(8)`).
- **Likely mapping:** `BM2` + `MieGruneisenDebye`, explicitly using
  `Tr = 700 K` rather than silently converting to 300 K.
- **Effort/risk:** low. Validate the cell/formula-unit convention and do not
  merge parameters with the existing Shim or Mao 300 K records.
- **Audit result:** preferred model 1 is represented as one non-default,
  pure-composition cubic record on the Fei et al. (2004) Pt scale. Equations
  (1)--(5), Tables 1--3, the `Z=1` pseudo-cubic volume basis, all coefficient
  errors, and the three rejected broad-peak observations were checked against
  the complete primary article. The Holmes-scale and non-MGD alternatives are
  provenance only and are not separate production records.
- **Reproduction:** an unweighted staged Peritheos refit to the bundled Table 1
  transcription gives `V0 = 46.5005 A^3`, `K0 = 207.369 GPa`,
  `theta0 = 1292.4 K`, `gamma0 = 2.7215`, and `q = 1.2646`, reproducing all
  published model-1 coefficients and their error magnitudes. The complete
  table is bundled with source terms unspecified and scoped CC0 for contributors’
  transcription rights. The reproduction script verifies both audited data
  digests, all four staged fits and historical Pt-scale checks.
  See the [full reproduction](literature-reproductions.md#c04-casio3-perovskite-noguchi-et-al-2013).

### C05 — Mg2SiO4 ringwoodite: add the thermal branch

**Status: blocked after primary-source and numerical audit (2026-09-05); no
production EOS record added.**

- **Change:** extend the existing Mg end-member `ringwoodite` material.
- **Source:** Katsura et al. (2004), [Thermal expansion of Mg2SiO4 ringwoodite
  at high pressures](https://doi.org/10.1029/2004JB003094).
- **Published scope:** 0--24 GPa and 300--2000 K; `V0 = 524.8(1) A^3`, fixed
  `K0 = 182 GPa`, `K0' = 4.6(2)`, `theta0 = 846(26) K`, `gamma0 = 1.93(3)`,
  and `q = 3.5(3)`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`, `Tr = 300 K`.
- **Audit result:** the official Wiley version of record supplies all 127
  Table 2 observations and identifies the Matsui et al. (2000) MgO pressure
  scale. The paper defines `n` as atoms per formula unit but does not print the
  value used or a Debye-temperature volume law. With chemically required
  `n = 7`, the published coefficients miss the table by 1.890 GPa RMSE
  (2.013 GPa for heated rows) under Peritheos' integrated-Gruneisen mapping;
  the alternate variable-exponent law is no better. An independent `n = 7`
  refit moves `gamma0` from 1.93(3) to 1.3996(836), far outside the reported
  error. An unphysical, unreported `n = 5` diagnostic reduces the heated-row
  RMSE to 0.242 GPa, indicating an unresolved hidden normalization rather than
  a defensible parameter choice. The exact blocker and direct table
  transcription are documented in the
  [reproduction](literature-reproductions.md#ringwoodite-katsura-2004).
- **Disposition:** keep the existing Meng static fit unchanged. Add Katsura as
  its own literature record only if an authoritative source resolves the
  numeric atom count, Debye-temperature law, and energy/volume normalization.

### C06 — Bcc vanadium: practical high-temperature marker

- **Implemented (2026-09-19):** main experimental fit in `vanadium_bcc`;
  [source audit](literature-reproductions/crichton-2016-vanadium.md).
  The two fixed-K′ alternatives remain withheld for missing coefficients.
  Plot-only observations support a proxy diagnostic, not the original refit.
- **Change:** add a new elemental material and P-V-T record.
- **Source:** Crichton et al. (2016), [High-temperature equation of state of
  vanadium](https://doi.org/10.1080/08957959.2015.1123256).
- **Published scope:** to 11.5 GPa and 1000 K; BM3 with
  `K0,300 = 150.4(6.2) GPa`, `K0' = 5.5(1.0)`,
  `alpha(T) = 4.8(6)e-5 - 2.4(9)e-8 T`, and
  `dK/dT = -0.0446(7) GPa/K`.
- **Likely mapping:** `BM3` + `ThermalReferenceStateEOS`.
- **Effort/risk:** low-to-medium. Recover the exact `V0` convention from the
  paper, express its alpha polynomial in Peritheos's reference-temperature
  convention, and encode the bcc-to-rhombohedral upper-pressure limit.

### C07 — modern B2 RbCl, RbBr, and RbI thermal pressure-marker family

- **Change:** extend the new Campbell-Heinz (1994) 298 K RbCl-B2 material with
  a modern thermal EOS, and add RbBr and RbI as two new materials.
- **Source:** Farla et al. (2025), [Thermal equations of state of B2-structured
  rubidium halides RbCl, RbBr, and RbI](https://doi.org/10.1063/5.0248905).
- **Published scope:** to 26 GPa and 1800 K. All three use BM3 + MGD. Reported
  `(K0, K0', gamma0, q)` values are RbCl `(19.89(8), 5.00(2), 1.96(4),
  1.05(9))`, RbBr `(16.28(4), 5.28(2), 2.18(14), 1.52(24))`, and RbI
  `(13.69(4), 4.95(1), 2.21(7), 1.42(10))`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye` for each phase.
- **Effort/risk:** medium. Extract `V0`, `theta0`, atom counts, transition
  bounds, and structures from the full paper/supplement; keep pressure-scale
  provenance explicit because CsCl, Mo, and Pt were used in calibration.

## B — useful next

### C08 — Natural siderite Fe0.95Mn0.05CO3

**Current status (2026-10-15): Implemented: published static record and separate joint thermal refit; published thermal alternatives remain withheld. See the current audit.**

- **Source:** Litasov et al. (2013), [P-V-T equation of state of siderite to
  33 GPa and 1673 K](https://doi.org/10.1016/j.pepi.2013.07.011).
- **Published scope:** BM3 with `V0 = 293.4(1) A^3`, `K0 = 120(1) GPa`,
  `K0' = 3.57(9)`, `dK/dT = -0.015(1) GPa/K`, and a reported linear
  ambient-pressure expansivity.
- **Likely mapping:** `BM3` + `ThermalReferenceStateEOS`.
- **Effort/risk:** medium. Name and formula must reflect the measured natural
  composition rather than claim a pure FeCO3 end-member; set the upper range
  below any spin transition or decomposition not represented by the fit.

### C09 — Superhydrous phase B, Mg10Si3O14(OH)4

**Current status (2026-10-15): Implemented: seven records, including five thermal alternatives; no longer an open material gap.**

- **Source:** Litasov et al. (2007), [Thermal equation of state of
  superhydrous phase B to 27 GPa and 1373 K](https://doi.org/10.1016/j.pepi.2007.06.003).
- **Published scope:** BM3 + MGD to 27 GPa and 1373 K. The preferred MGD table
  includes approximately `V0 = 623.53(35) A^3`, `K0 = 135.5(2.3) GPa`,
  `K0' = 5.3(2)`, fixed `theta0 = 860 K`, `gamma0 = 1.33(5)`, and
  `q = 2.03(35)`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`.
- **Effort/risk:** medium. The paper provides fits on two Au pressure scales;
  select one record deliberately or preserve both as separately named scales.

### C10 — Fe2SiO4 ringwoodite

- **Source:** Armentrout and Kavner (2011), [High pressure, high temperature
  equation of state for Fe2SiO4 ringwoodite](https://doi.org/10.1029/2011GL046949).
- **Published scope:** MGD with `V0 = 42.03 cm3/mol`, `K0 = 202(4) GPa`, fixed
  `K0' = 4`, fixed `theta0 = 685 K`, `gamma0 = 1.08(6)`, and fixed `q = 2`.
- **Likely mapping:** `BM2` + `MieGruneisenDebye` as a new Fe end-member
  material, not as another record for Mg-ringwoodite.
- **Effort/risk:** medium. Convert molar to conventional-cell volume exactly
  and preserve that several thermal parameters were adopted/fixed.

### C11 — Dense hydrous phase E

**Integrated 2026-10-16:** `phase_e_fe0072`, the measured Fe-bearing sample.

- **Source:** Crichton and Ross (2000), [Equation of state of phase
  E](https://doi.org/10.1180/002646100549427).
- **Published scope:** single-crystal compression to 6.7 GPa for measured
  `Mg1.96(7)Fe0.072(5)Si1.04(5)H3.7(8)O6`, with BM3
  `K0 = 92.9(7) GPa` and `K0' = 7.3(2)`.
- **Bundled record:** published BM3 `V0 = 105.627(7) Å³`,
  `K0 = 92.9(7) GPa`, `K0' = 7.3(2)`, all freely refined. All eleven
  primary P,a,c,V rows and their reported errors are bundled.
- **Parity:** successfully recovered within reported uncertainty. The weighted
  refit shifts each coefficient by less than 0.123 published error widths;
  the maximum curve difference is about 0.0096 GPa. Published coefficients
  remain unchanged. [Audit and plot](literature-reproductions/phase-e-crichton-ross-2000.md).
- **Structure/calibration:** R-3m, conventional hexagonal cell, O6-normalized
  empirical formula with Z=1; eight indexed phase E reflection positions
  provide the documented diffraction fallback. Target-sample occupied sites
  and H positions are not inferred from different-composition structures.
  Angel (1997) quartz scale is identified; paired quartz readings are absent.
- **Remaining scope:** exact original regression details are optional provenance
  gaps, not failed parity. Wang (2022) thermal data concern a distinct Fe-free
  sample and are not added to this material.

### C12 — Brucite Mg(OH)2

**Current status (2026-10-17): three published isotherms retained; full thermal EOS not recommended.**

- Default: unchanged Fukui (2003) 300 K BM3 through 22.56 GPa; published
  473/673 K slices retain their measured ranges. Hydrogen rearrangement does
  not impose an unsupported 0.5 GPa bulk volume cutoff.
- **Full published thermal EOS not recommended:** parity with the published
  P–V–T data was not achieved. The 873 K normalization and continuous source
  reconstruction remain unresolved. Their unchanged coefficients and all
  Peritheos refits are retained as diagnostic evidence outside executable selection.
- Only the published 300/473/673 K isotherms remain executable, each within its
  source-checked measured pressure span. These do not provide interpolation
  between temperatures. Source observations, errors and primary ambient structure
  are preserved. Further model development is outside this curation task.
  See the [full audit](literature-reproductions/brucite.md).

### C13 — Delta-AlOOH

**Current status (2026-10-15): Implemented: both Vanpeteghem (2002) BM2/BM3 alternatives; thermal coverage remains absent.**

- **Source:** Vanpeteghem et al. (2002), [Equation of state of the hydrous phase
  delta-AlOOH at room temperature up to 22.5 GPa](https://doi.org/10.1029/2001GL014224).
- **Published scope:** two reported BM fits: `K0 = 252(3) GPa` with fixed
  `K0' = 4`, or `K0 = 228(3) GPa` with `K0' = 7(1)`.
- **Likely mapping:** paired BM2 and BM3 records, if the full data audit shows
  both are intended as publishable alternatives.
- **Effort/risk:** medium. Later work reports hydrogen-bond and spin-related
  changes; validity must stop before behavior not captured by a single BM fit.

### C14 — Fe-Si system: B20/B2 FeSi and Fe-9 wt% Si phases

- **Audited (2026-09-19):** seven published records accepted: BM3 and Vinet
  for B20 FeSi, D03 Fe–9Si and hcp Fe–9Si, plus the thermal D03 Fe–16Si
  Vinet reanalysis. All seven supplementary tables and Text S1 were inspected.
- **Source:** Fischer et al. (2014), [Equations of state in the Fe-FeSi system
  at high pressures and temperatures](https://doi.org/10.1002/2013JB010898).
- **Withheld:** both B2 FeSi thermal fits fail against the 114 single-phase
  observations using physical Debye normalization. Doubling thermal-energy
  amplitude recovers the reported residual range, but is not an authoritative
  correction. Two separately labelled, physically normalized Peritheos B2 refits
  are now executable with conditional covariance and sensitivity evidence.
- **Excluded:** three hcp+B2 aggregate fits and eight computational branches.
  Full observations, calibrant inputs, constraints, structure approximations,
  independent diagnostics and all dispositions are in the
  [source audit](literature-reproductions/fischer-2014-fesi.md).

### C15 — Ferropericlase Mg0.75Fe0.25O across spin crossover

**Current status (2026-10-15): Partly implemented: bounded spin branches and a thermal record exist; continuous crossover remains model work.**

- **Source:** Mao et al. (2011), [Thermal equation of state of lower-mantle
  ferropericlase across the spin crossover](https://doi.org/10.1029/2011GL049915).
- **Published scope:** to 140 GPa and 2000 K. Separate high-spin and low-spin
  reference fits are reported; at 300 K, the high-spin state uses
  `V0 = 76.34 A^3`, `K0 = 162(1) GPa`, fixed `K0' = 4`, while the low-spin
  fit uses `V0 = 74.44(60) A^3`, `K0 = 166(7) GPa`, fixed `K0' = 4`.
- **Likely mapping:** two bounded BM2 state records are possible immediately;
  the continuous crossover is not represented by existing Peritheos classes.
- **Effort/risk:** high scientific-review burden. Preserve composition and
  pressure-scale choice, and do not present a piecewise switch as the paper's
  continuous crossover model.

### C16 — Fe- and Al-bearing phase D P-V-T

- **Source:** Litasov et al. (2008), [Thermal equation of state of Al- and
  Fe-bearing phase D](https://doi.org/10.1029/2007JB004937).
- **Published scope:** measured
  `Mg0.99Fe0.12Al0.09Si1.75H2.51O6` to 20.6 GPa and 1273 K, analyzed with both
  high-temperature BM and MGD approaches.
- **Likely mapping:** new composition-specific `BM3` +
  `MieGruneisenDebye` record, rather than attaching thermal parameters to the
  existing antigorite-derived phase-D cards.
- **Effort/risk:** medium. Extract the complete preferred parameter table and
  distinguish scientific-fit provenance from the existing two phase-D
  compositions.

### C17 — Dolomite/ankerite family and dolomite-III

**Current status (2026-10-15): Partly implemented: the Mao (2011) measured-composition cards exist; broader low-pressure dolomite/ankerite coverage remains open.**

- **Sources:** Ross and Reeder (1992), [High-pressure structural study of
  dolomite and ankerite](https://msaweb.org/AmMin/AM77/AM77_412.pdf), and Mao
  et al. (2011), [Dolomite III: a new candidate lower mantle
  carbonate](https://doi.org/10.1029/2011GL049519).
- **Published scope:** the first paper gives low-pressure single-crystal P-V
  data for stoichiometric dolomite and about 70 mol% ankerite; the second gives
  a P-V EOS for Fe-bearing monoclinic dolomite-III observed at 36--83 GPa.
- **Likely mapping:** separate low-pressure dolomite, ferroan-dolomite/ankerite,
  and high-pressure dolomite-III records.
- **Effort/risk:** medium-to-high. Do not splice polymorphs or compositions.
  The old low-pressure paper has no DOI and may require refitting its published
  data; the high-pressure fit may involve a strongly extrapolated zero-pressure
  reference volume.

### C18 — Hcp zinc P-V-T

- **Source:** Errandonea et al. (2018), [High-pressure/high-temperature phase
  diagram of zinc](https://doi.org/10.1088/1361-648X/aacac0).
- **Published scope:** experimental XRD to 16 GPa and 1000 K, with wider
  theoretical coverage; the 300 K BM3 fit gives `V0 = 30.6(1) A^3`,
  `K0 = 63(2) GPa`, and `K0' = 5.6(4)`.
- **Likely mapping:** new hcp-Zn BM3 plus the paper's thermal reference-state
  formulation after equation-level audit.
- **Effort/risk:** medium. Record the strong axial anisotropy and proposed
  isomorphic anomaly near 10 GPa in provenance/validity notes even though the
  scalar volume EOS remains continuous.

### C19 — Molybdenum: independent experimental MGD scale

- **Change:** add an alternate record to existing `molybdenum`.
- **Source:** Huang et al. (2016), [Thermal equation of state of molybdenum
  determined from in situ synchrotron X-ray diffraction](https://doi.org/10.1038/srep19923).
- **Published scope:** 300 K compression to 80 GPa and thermal data to 92 GPa
  and 3470 K. The 300 K BM3 fit gives `V0 = 31.22(8) A^3`,
  `K0 = 273(15) GPa`, and `K0' = 3.6(4)`; the MGD fit fixes the Debye
  temperature near 470 K and fits `q = 0.6`, alongside a separate
  thermodynamic thermal-pressure fit.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`, as an explicitly named
  alternative to the current Sokolova composite scale.
- **Effort/risk:** medium. Resolve which of the paper's two formulations is the
  primary executable record and preserve the fixed-versus-fitted status of the
  Debye temperature.

## C — model work

### C20 — (Al,Fe)-phase H with a spin crossover

- **Source:** Strozewski et al. (2023), [Equation of state and spin crossover
  of (Al,Fe)-phase H](https://doi.org/10.1029/2022JB026291).
- **Published scope:** `Al0.84Fe3+0.07Mg0.02Si0.06OOH` to 125 GPa. The paper
  fits high- and low-spin elastic states simultaneously and places the 300 K
  crossover midpoint at `55.5(9) GPa`, spanning roughly 48--63 GPa.
- **Likely mapping:** a new spin-crossover free-energy EOS, not an ordinary BM3
  record. The paper's MINUTI formulation adds a spin partition-function term
  to the elastic free energy.
- **Effort/risk:** high. This is the clearest candidate for a new general model
  class, but implementation should wait until its thermodynamic variables,
  parameter identifiability, inversion behavior, and serialization contract
  have been designed.

### C21 — KAlSi3O8 liebermannite and K-hollandite II

**Status: audited and deferred (2026-09-04); no executable record created.**

- **Source:** Wang et al. (2026),
  [doi:10.2138/am-2024-9562](https://doi.org/10.2138/am-2024-9562).
- **What worked:** independent BM3 fits reproduce the deposited, rounded
  temperature-indexed pressure-volume grids.
- **Why it remains out:** the grids are derived curves rather than primary fit
  observations; the fitted `V0`, `K0`, and `K0_prime` all vary with
  temperature in a way the current thermal wrapper cannot represent. The
  accessible deposit also does not resolve the exact reference-state method,
  parameter uncertainties/covariance, or authoritative symmetry-reduced
  structures. Encoding individual temperature slices as separate static EOS
  records would misrepresent the paper.
- **Disposition:** retain only as a deferred investigation. Reconsider when a
  source-faithful temperature-continuous model and the missing primary-method
  evidence are available. It is listed explicitly in the
  [paper investigation ledger](paper-investigation-ledger.md).

## Current review order

Use the open additions and existing-record verification queues above. G01–G12
remain stable historical identifiers, not twelve unfinished tasks. Completed
integrations and completed non-reproduction investigations must not be
reintroduced as missing materials or active fitting work. Evidence-dependent
holds remain separate from work that can progress using available inputs.

For each review, the acceptance gate should be: obtain the primary full text
and supplement; identify exact composition, phase, cell and volume convention;
transcribe the published equation rather than only its parameter table;
reproduce at least one published P-V(-T) point; record fixed/fitted parameters
and uncertainty meaning; define the experimentally supported validity domain;
and only then create or modify an `.eosmat` record.

## Accepted direct lead: Ono et al. (2006), CaTi2O4-type MgAl2O4

Accepted 2026-09-08 as `mgal2o4_cati2o4_ono_2006_bm2_2`, a 300 K BM2
with V0=238.9(9) Å³ and K0=219(6) GPa, K0′=4 fixed. All 14 Table 3
observations are bundled. No same-DOI/title row exists in the LitCurate
discovery ledger; no intake identifier was invented. See the
[primary-source audit and reproduction](literature-reproductions/ono-2006-mgal2o4-cati2o4.md)
for uncertainty, fit-selection sensitivity, calibration, and comparison-fit dispositions.


## Direct primary-source additions

- Irifune et al. (2002), DOI `10.1007/s00269-002-0275-1`: accepted one
  constrained CaFe2O4-type MgAl2O4 BM2 record on 2026-09-08. No matching
  LitCurate discovery row exists. See the [sparse-data audit](literature-reproductions/irifune-2002-mgal2o4-cafe2o4.md).
