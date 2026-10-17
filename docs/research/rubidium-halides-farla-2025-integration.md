# Rubidium halide integration handoff

**Scoped commit prepared on 2026-10-16.** The rubidium material records,
datasets, reproduction evidence and documentation are committed separately
from unrelated material and model work. Shared integration was also reviewed
in the local `docs/research/catalog-integration-2026-10-16.md` handoff; that
combined working-tree report and aggregate catalog reconciliation are outside
this commit. Earlier pending-check lists below preserve the scoped handoff
history; scientific evidence gaps remain as stated in the source audit.

Initial date: 2026-10-15; direct EosFit follow-up: 2026-10-16.
Scope G03/C07, Farla et al. (2025),
doi:10.1063/5.0248905. Local additions are complete, subject to central
reconciliation of shared indices and the explicitly qualified source audit.
No changes have been pushed or published, and no author or other chat was
contacted. The separate commit preserves the material-specific scope below.

## Current follow-up decision (2026-10-16)

The remaining differences are documented in the family audit's
[assessment and follow-up decision](../literature-reproductions/rubidium-halides-farla-2025.md#assessment-and-follow-up-decision-2026-10-16).
Preferred Case 1 differs by at most 1.63 published marginal parameter standard
errors; constrained Case 3 by less than 0.085. Five CsCl readings exceed
pressure/volume rounding but differ by at most 0.07090 GPa, comparable to their
0.07 GPa reported errors. The localized RbCl S5 consistency discrepancy is
approximately 0.01% in relative volume. These results do not establish a
substantive scientific error; the parameter comparisons are marginal and the
original covariance remains unavailable.

The user decided to document these limits and **defer author contact for now**.
No author has been contacted. Previously prepared correspondence and the PDF
are unsent drafts, not an active outreach task. Missing original files and
full-precision inputs remain evidence gaps. Published values, source-table
entries and false original-author-fit parity flags remain unchanged.

## Owned files and catalog delta

- `peritheos/data/materials/rbcl.eosmat`: append preferred case-1 thermal
  record and datasets. Existing static record, dataset, structure, all other
  material metadata and default-selection behavior remain unchanged.
- `peritheos/data/materials/rbbr.eosmat` and `rbi.eosmat`: new Z=1 Pm-3m
  B2 cards with ordered Rb/halogen sites and hypothetical reference cells.
- `peritheos/data/datasets/rubidium-halides-farla-2025-*` and the family
  `.LICENSE.md`: six CSV resources, all seven raw numerical source tables,
  source manifest and rights/attribution notice.
- `scripts/recover_rubidium_halides_farla_2025.py` and
  `scripts/reproduce_rubidium_halides_farla_2025.py`.
- `tests/test_rubidium_halides_farla_2025.py` and
  `crates/peritheos/tests/rubidium_halides_farla_2025.rs`.
- `docs/data/rubidium-halides-farla-2025-reproduction.json`,
  `docs/literature-reproductions/rubidium-halides-farla-2025.md`, this note.
- Follow-up calibration/benchmark diagnostic scripts, reports and tests;
  `scripts/run_rubidium_halides_farla_2025_eosfit.py`, its archived-output
  test and the complete `docs/data/rubidium-halides-farla-2025-eosfit7c/`
  evidence directory. A narrowly scoped `.gitattributes` addition protects
  these native files and their hashes from text/line-ending normalization.
- Constrained Case 3 source-input facts, runner, tests and
  `docs/data/rubidium-halides-farla-2025-case3-eosfit7c/` native evidence.
  This includes both RbBr alternatives and explicitly labelled fixed-anchor
  rounding examples. The Case 1 canonical records and source tables retain
  their published values; the new evidence is linked as a separate diagnostic.

This family's isolated delta is **+2 materials, +3 published validated
thermal records, +2 diffraction-ready materials**. RbCl remains present with
its original static default. All three new parameterizations are classified
**hybrid** because their fixed Debye temperatures derive from DFT elastic
constants. New RbBr/RbI cards default to their preferred case-1 records.
Do not apply these deltas blindly to a concurrent catalog count; enumerate
the final reconciled working tree.

Record identifiers:

- `rbcl_b2_farla_2025_bm3_mgd_cscl`
- `rbbr_b2_farla_2025_bm3_mgd_cscl`
- `rbi_b2_farla_2025_bm3_mgd_cscl`

Each links a compound-specific `*_farla_2025_observations` dataset and a
separate `*_farla_2025_benchmarks` dataset. There are 376 retained printed
observation rows, including one blank RbI sample row, and 162 calculated
benchmark ratios. Paired calibrant volumes/errors and literal error tokens
are retained. Source data have CC BY 4.0 attribution; contributor CC0 scope
is restricted to contributor-owned transcription/normalization/metadata.

## Required shared reconciliation

These files were intentionally not edited in this chat:

- `docs/material-eos-candidates.md`: reconcile G03/C07 as implemented with
  source/fit/calibration qualifications; do not describe exact fit parity.
- `docs/catalog-coverage.md`, `docs/data/catalog-coverage.json`, catalog
  counts, README, catalog pages and shared count assertions: regenerate once
  across all concurrent material work.
- `peritheos/data/primary-source-audit.json` and its generated summaries:
  add the three source-validated records with the official accepted-manuscript
  version and public AIP supplement. The publisher PDF was not downloadable;
  do not say the final typeset PDF was inspected.
- Global primary-EOS/refit ledgers, their generator and source-audit index:
  integrate the material-specific reconstruction report. Equation validation
  and coefficient-rounding consistency pass. Conditional independent
  refits are available. Original author fitted-coefficient/covariance and
  workflow parity remain **false**, not unresolved zeros or implicit passes.
- `docs/dataset-requests.md`: if source-workflow recovery is pursued, request
  original EoSFit7 and GSAS-II files, exact row flags, unrounded errors/PVT
  correlations and covariance, raw ruby wavelengths, and saved CsCl runtime
  settings/unrounded inputs behind five pressure discrepancies. Historical
  CsCl coefficients and GUI defaults have now been recovered. No outreach
  was performed; the user explicitly excludes contacting authors.
- Documentation navigation/CHANGELOG: add the family audit when integrating
  the combined additions. This report uses code-formatted repository paths
  for non-doc resources, avoiding broken MkDocs links.

No legacy catalog-compatibility manifest change is required for discovery:
the canonical card loader enumerates the `.eosmat` files, and the public API
loads both new materials successfully. No new equation, schema, native model
or material-family definition is needed.

## Calibration links requiring care

No CsCl, Mo or Pt card and no shared pressure-calibration registry entry was
edited. Do not substitute the existing Campbell-Heinz CsCl static record for
the source's Dewaele-static plus Decker-thermal composite scale.

- DAC ruby lineage is IPPS-Ruby2020. Existing registry identifier
  `ruby_shen_2020` is an appropriate integration link, but original wavelengths
  are absent. The source-reported pressures/errors remain the observations.
- Preferred LVP calibration is a composite CsCl scale: Dewaele (2020) static
  and “Decker et al. 1973” thermal in S2-S4 footnote a; main reference 2 says
  Decker (1971). Historical EosCross v1.3 Vinet-MGD coefficients, a0=4.1221 A
  GUI override, legacy Avogadro constant and BurnMan 1.2.0 dependency are now
  recovered in the family author-inputs JSON. Saved runtime settings remain
  unavailable; 104/109 distinct readings match pressure/volume rounding and
  five exceed that precision. All match reported errors plus rounding.
- Mo comparison pressures: Huang et al. (2016), doi:10.1038/srep19923, as
  stated by footnote c. Do not substitute Litasov/Sokolova alternatives merely
  because the main paper mentions them as comparison calculations.
- Pt comparison pressures: footnote b says Matsui (2007); main reference 30
  is Matsui (2009), doi:10.1063/1.3054331. Existing
  `platinum_matsui_2009_vinet_mgd_electronic` is a possible candidate link
  after checking source version and the footnote discrepancy. It is not
  linked as if that ambiguity were settled. Pt is not the preferred fit gauge.

The cards retain `pressure_calibration.status=partially_resolved` and
`recalculation.status=reference_eos_not_bundled`. Source-reported pressure
columns are executable fit observations; no full pressure-scale replay is
claimed. The independent family script replays historical CsCl defaults;
that implementation is not a registered shared calibrant record. Above about
1100 K the reported Mo/Pt discrepancies and the
possible CsCl BT815 issue remain visible in the audit.

## Scientific qualifications

The 300 K BM3 plus integrated Debye model, reference subtraction, n=2,
fixed V0/theta0, free parameters and two-atom B2 cells are resolved.
The canonical V0 and its error are converted from the molar Table I row;
the inconsistent rounded cell-volume uncertainties are preserved separately.
Independent checks reproduce the fixed Debye temperatures from the
supplemental elastic inputs after rounding, despite apparent formula typos.

All 162 source model ratios agree within an explicit state-specific
coefficient-rounding envelope; many do not agree at the final ratio digit.
Native versus independent quadrature differences are below 1.1e-14 in V/V0.
The initial independent fit uses a documented effective-variance approximation,
diagonal PVT errors and a diagnostic 0.005 GPa pressure-error floor. The
follow-up freezes weights inside each least-squares cycle, updates them
between cycles and removes the floor. This brings marginal parameter
differences below 1.65 reported esds but does not reproduce all coefficients
at publication rounding. A common coefficient vector within rounding bounds
fits all RbBr/RbI calculated ratios. The subsequent single-entry omission
check localizes the remaining RbCl consistency issue to
5 GPa/1500 K; all other 53 entries admit one common rounding-compatible vector.
Constrained Case 3 thermal pairs agree within 0.085 reported marginal esds;
illustrative fixed-anchor rounding tests recover their printed digits without
recovering the authors' actual fixed inputs. The current follow-up decision
above supersedes any earlier implication that these differences require outreach.
Cold-LVP selection sensitivity and unweighted fits are retained. No refit
record was added and published coefficients remain unchanged. The family
audit's recovery and fit parity follow-up gives the source code, five gauge
outliers and limits of the public search for original EoSFit7/GSAS files.

Phase boundaries are not fully quantified by the source. The cards use
marginal observation bounds, with notes that hot data extend only to about
21 GPa, not the 26 GPa room-temperature DAC maximum. The S5-S7 5 GPa /
2500 K molten flags remain benchmarks only. Near-transition RbI has a source
S1 B2-label versus S7 B1-caption ambiguity; no invented reassignment is made.

## Verification log

- Scoped commit verification on 16 October 2026 used an exported snapshot of
  HEAD plus only the selected rubidium changes: **95 family Python tests**
  and **1 dedicated Rust test** passed. The rubidium runners/tests and the
  scoped common-ledger adapter pass lint and formatting; Rust test formatting
  and the staged whitespace check pass. The three adapter results reproduce
  the stored diagnostic coefficients. Existing shared audit/refit entries
  are preserved exactly; only the three Farla records are added. Broader
  aggregate catalog reconciliation and unrelated working-tree changes are
  deliberately excluded. This is not a full-repository test-pass claim.
- Constrained Case 3: **24 native fits**, all converged in two passes.
  With literal published anchors the weighted thermal pairs are
  RbCl (1.96029, 1.04704), RbBr/DAC (2.17459, 1.50152),
  RbBr/K0'=4.98 (2.11447, 1.50637), RbI (2.21407, 1.42931).
  All deviations are below 0.085 reported marginal thermal-parameter esds.
  Illustrative sub-last-digit fitted-DAC anchors make both free coefficients
  round to the publication for all four variants. The deliberately selected
  RbBr K0'=4.98 remains exact; only its V0 changes within printed precision.
  These demonstrate precision consistency without recovering author inputs.
- Case 3 evidence tests: **37 passed**; all family tests combined:
  **95 passed**. Tests validate all four fixed anchors, only gamma0/q free,
  the zero fixed-anchor covariance, native thermal covariance/esds, source
  mapping, hashes, archived macros, and the rounding-example qualifications.
  Cold 300 K row inclusion changes conditional error scaling but not the
  two thermal parameters at native printed precision. The Case 1 runner's
  default macros/data remain byte-identical after adding Case 3 options.
- Direct official EosFit7c 7.60 execution: **18 cases**, all converged in
  both refinement passes. Weighted precise-theta parameter tuples are
  RbCl (19.92227, 4.99511, 1.96401, 1.05929),
  RbBr (16.31590, 5.25587, 2.16658, 1.47412),
  RbI (13.75502, 4.93392, 2.19266, 1.41792).
  Their largest marginal deviations from the publication are 0.40, 1.21
  and 1.63 reported esds. The matching independent IRLS differs by at most
  0.000326 GPa in K0 and 0.000074 in the dimensionless parameters.
- Combined family Python tests after direct-console additions: **58 passed**;
  dedicated Rust integration test: **1 passed**. New runner/tests pass lint
  and formatting. Tests regenerate data and row maps, validate archived
  macros at their original paths, all file hashes, full-MGD/fixed anchors,
  native covariance/errors, two-pass convergence and finite final weights.
  Native source-pressure replay agrees with independent quadrature at its
  printed output precision. The official console executable was not modified.
  Task-local runtime payloads were extracted without installing applications;
  the private display was stopped and temporary authorization removed.
- Actual EosFit warnings are retained. Zero source errors remain zero except
  in the explicitly named pressure-floor sensitivity; all final weighted
  effective variances remain finite. Native refit covariance is retained
  separately and does not become the covariance of the published EOS records.
- Initial material-specific Python suite: **20 passed**, warnings treated as errors.
- Recovery/weighting follow-up plus original material-specific tests:
  **29 passed** on 15 October. This includes independent CsCl/native-model agreement,
  all five pressure discrepancies, iterative weighting and compatible common
  vectors for the rounded RbBr/RbI benchmark tables. Follow-up lint and
  formatting passed. The material-specific Rust test was rerun: **1 passed**.
- The 16 October discrepancy suite plus both earlier family suites:
  **36 passed**. Historical BurnMan 1.2.0 execution agrees with the independent
  CsCl calculation within 6.27e-11 GPa. Half-integer temperature rounding does
  not repair the five pressure discrepancies; the stated 1 K target-control
  tolerance makes one compatible but leaves four. Neither a common a0 nor
  a uniform temperature offset repairs all readings at recovered model settings.
- RbCl's benchmark issue is localized by 54 single-entry omission tests.
  Only omitting 5 GPa/1500 K allows all remaining 53 ratios to fit one vector
  inside the four coefficient half-digit bounds. Its reconstructed ratio rounds
  to 0.9283; S5's original 0.9282 remains unchanged. This supersedes the earlier
  two-outlier localization, not the explicit false original-fit parity status.
  New evidence is in `docs/data/rubidium-halides-farla-2025-discrepancies.json`;
  external historical execution is reproducible through the separate
  `scripts/replay_rubidium_halides_farla_2025_burnman.py` optional-dependency
  command in the family audit. No author contact, shared-ledger regeneration,
  coefficient replacement or source-table correction was performed.
- Material-specific Rust integration test: **1 passed** with all features.
- Targeted Rust Clippy with warnings denied: passed.
- Python lint and formatting: passed repository-wide when run; targeted
  changed scripts/tests also pass after the final edits.
- Mypy (Python 3.9): passed, 56 source files.
- Rust formatting: passed repository-wide after formatting the owned test.
- RbCl preservation check: original record, original dataset and every
  pre-existing material-level field equal their HEAD values exactly.
- Full Rust workspace tests stop at the existing/concurrent Au
  `all_bundled_material_records_load_and_round_trip_through_rust` assertion.
  The failure is in the gold document, not a rubidium card.
- Full workspace Clippy: passed on the final run. An initial run stopped on
  an unreadable literal in the concurrent Diopside test; that was resolved
  outside this scope. No unrelated test was changed here.
- Strict MkDocs build: failed on eight links in concurrent/pre-existing
  Diopside and Au reports. The rubidium report contributes no link warning.
- Shared primary-refit `--check`: failed before comparison with
  `KeyError: 'primary_source_check'` for the pre-existing
  `gold_hirose_2008_fei_2004_speziale_2001_bm3_refit` record. All three
  rubidium records supply that field. No shared ledger was regenerated.
- Adjacent Python regression (`test_eosmat.py`, `test_catalog.py`,
  `test_datasets.py`, `test_materials.py`): **531 passed, 8 failed**. Failures
  are the stale material/record/validation counts, the shared primary-audit
  record set, and the existing Fei Au temperature expectation. At this
  concurrent snapshot assertions expected 227 materials / 650 records /
  640 validated records; enumeration returned 230 / 657 / 647. These are
  combined-checkout counts, not this family's isolated delta.
- Full Python suite with coverage was attempted twice. The initial run was
  interrupted in an unrelated Argon refit after 243 passes and the two
  known `test_catalog.py` failures (454.82 seconds). A retry with numerical
  threads limited and a family-specific coverage file reached the unrelated
  Dewaele six-metal source-refit tests, but was stopped to avoid spending
  further time on expensive unchanged refits after the adjacent suite
  completed. It did not produce a complete coverage result. **No full-suite
  pass or coverage percentage is claimed.**

The caller's baseline already had two `test_catalog.py` expectation failures
(validated-record count and Fei Au temperature bounds). Additional shared
count/index failures can arise from all concurrent additions. Those are
integration tasks, not reasons to edit Au data or other materials in this chat.
