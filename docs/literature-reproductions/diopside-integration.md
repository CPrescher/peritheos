# Diopside integration note

G02 adds the material `diopside` locally. Thompson and Downs (2008) static BM3
is the room-temperature default; Zhao et al. (1998) supplies the published
thermal BM3. Published-model and parameter parity are established within
source precision and reported errors using the graph-corroborated table-error
treatment. The two sources remain independent; all published coefficients,
original observations and source-error provenance are retained. The
[source audit](diopside.md) distinguishes this result from exact recovery of
the original author workflow.

## Owned changes

- `peritheos/data/materials/diopside.eosmat`: one diffraction-ready C2/c,
  Z = 4 CaMgSi2O6 card with two published EOS records.
- `peritheos/data/datasets/diopside-*` and `diopside.LICENSE.md`: both complete
  EOS observation tables (14 static and 32 thermal rows), original numerical
  notation, printed errors, six ambient structure sites and their errors,
  source retrieval receipts, and contributor-only rights dedication.
- `scripts/diopside_reproduce.py` and
  `docs/data/diopside-reproduction.json`: independent equations, deterministic
  transcription/checksum checks, complete-table diagnostic refits and the
  explicitly separate thermal sensitivity calculation.
- `tests/test_diopside.py`, `crates/peritheos/tests/diopside.rs`, and the two
  diopside audit/integration documents.

## Shared reconciliation required

No common backlog, coverage snapshot, source-audit index, scientific ledger,
pressure-calibration library, README or common count assertion was edited.
After reviewing all concurrent material additions together:

1. Add `diopside` and its two records to the material/source audit manifest
   and primary-source validation ledger. The Python catalog discovers the new
   card automatically; no material-loader code change is needed.
2. Register `diopside_thompson_2008_bm3` and
   `diopside_zhao_1998_thermal_bm3` in the numerical validation/refit ledger,
   using both material-specific reports. Mark published-model parity and
   parameter agreement within reported errors **established**, with the
   plot-supported table-error treatment explicitly recorded. Keep exact original
   numerical input/workflow recovery **unrecovered** as a separate field.
   Source-table errors are provenance notes, not an EOS acceptance blocker.
   The audit refit does not replace the published record.
3. If tracking exact original workflow recovery, an optional dataset request
   can cover Zhao's original selected EOS input, V0 treatment, precise row-15
   pressure, weighting and covariance. This is not required for the established
   published-model/parameter parity. Paired NaCl readings would enable
   pressure-scale recalculation. Contact status is **not contacted**; no
   outreach is authorized by this work item.
4. Link the existing `ruby_mao_1978` calibration to the Thompson record in
   any shared consumer/dependency index. The record already identifies it.
   Zhao uses Decker (1971), DOI 10.1063/1.1660714, but no simultaneous NaCl
   readings are printed. No new calibration coefficients or speculative exact
   NaCl record mapping are required. Preserve the missing-observation status.
5. Regenerate coverage, common totals, release notes, navigation and references
   together with the other additions. This item contributes **one material,
   two executable published records, one thermal record**, two linked original
   observation datasets, and one separately registered corrected P–V–T analysis
   dataset. Do not overwrite concurrent materials or Au work.
6. Evaluate a separate curated clinopyroxene family only if warranted by the
   shared family definitions; no family was invented for this single card.

## Verification snapshot from 2026-10-15

The scoped Python source/schema/catalog/physics/data tests and Rust integration
test pass. The targeted combined EOS/catalog/family/thermal check completed with **207
passed and 7 failed**. Its failures are five common aggregate-count/ledger
expectations awaiting reconciliation and the two known Au/catalog expectations;
no material-specific failure remains.

The following broader outcomes apply to the concurrently changing checkout:

- Full Ruff lint and formatting checks passed; mypy passed for 56 source files
  using the Python 3.9 project environment.
- `cargo fmt --all -- --check` and full-workspace Clippy with warnings treated
  as errors passed after formatting only the owned Rust test.
- Full-workspace Rust tests reached an existing Au document equality failure in
  `all_bundled_material_records_load_and_round_trip_through_rust` (38 tests
  passed in that test binary). The dedicated diopside Rust test passes.
- Strict documentation build completed notebook execution but failed on five
  links in the current Hirose/Sakai Au reports. It reported no remaining
  diopside link warnings. Shared navigation is deliberately not edited here.
- The shared primary-refit `--check` run cannot complete: it raises
  `KeyError: primary_source_check` for the current
  `gold_hirose_2008_fei_2004_speziale_2001_bm3_refit` record. Both diopside
  records contain that metadata. Shared Au data and the global script were
  left untouched.
- The full Python coverage run was intentionally interrupted after 489.54 s
  (8 min 9 s) during a SciPy SVD computation, with **243 passed and 2 known
  catalog failures** at that point. It did not finish or establish coverage.
  The completed 214-test targeted run above is the broader integration evidence;
  this work does not claim a clean full-project test run.

No push, publication or author contact is part of this scoped addition.

## Plot check follow-up

The Figure 2 point corresponding to thermal Table 1 row 15 reads approximately
4.58 GPa (graphical bound ±0.05 GPa), corroborating a printed pressure error.
The EOS predicts 4.563 GPa; the table prints 6.95 GPa. The plotted ambient
volume also supports the volume column rather than the inconsistent ambient
lattice dimensions. Original tables, coefficients and defaults are unchanged.

Additional owned files are `scripts/diopside_check_plots.py` and
`docs/data/diopside-plot-check.json`; the existing audit and thermal record now
link this evidence. The annotated PNG in `output/diopside-plot-check.png` is an
inspection artifact containing an attributed source crop, outside the numerical
CC0 dedication. This graphical evidence should accompany shared-ledger
reconciliation without being relabeled as an exact recovered original pressure
or original-author workflow parity.

After the plot follow-up, all **8 material-specific Python tests** pass with
warnings treated as errors. The original reproduction checksum check and
material schema validation also pass. Broader-check limitations above remain
unchanged; no shared ledger or source observation was modified.

## Parity classification

On 2026-10-16 the assessment was updated to record **published-model and
parameter parity within source precision and reported errors**, accounting for
the graph-corroborated source-table error. The five varied thermal coefficients
all pass the comparison; source errors remain documented provenance. The
explicit unrecovered result is the exact original numerical input and fitting
workflow. No source observations, EOS coefficients or catalog defaults changed.

## Corrected dataset for reuse

The user requested a usable corrected dataset on 2026-10-16. The added
`peritheos/data/datasets/diopside-zhao-1998-table1-corrected-pvt.csv` applies
row 15 pressure = 4.58 GPa, with the printed 6.95 GPa retained in an original
pressure column, an explicit correction flag, and a separate ±0.05 GPa
non-statistical graphical reading bound. It keeps all 32 source rows, marks
the repeated ambient entry, and preserves temperatures, printed volumes and
reported errors. Unresolved lattice dimensions remain only in the original
source transcription.

Additional owned files are `scripts/diopside_correct_dataset.py` and the
portable `diopside-zhao-1998-corrected-pvt.README.md` beside the CSV. The
material document registers the corrected resource, checksum, units and
provenance as `diopside_zhao_1998_table1_corrected_pvt` and recommends it for
refitting. It is a supporting analysis dataset, not the original author input
or a new executable refit record. The plot reproduction uses exactly the
rounded exported pressure; original source tables, coefficients and defaults
are unchanged. Shared reconciliation should include this additional derived
dataset alongside the two original observation datasets.

The completed corrected-data checks pass: **9 Python tests and 1 Rust
integration test**, plus schema validation, metadata round trips, source and
corrected-export checksum checks. Shared refit automation should select the
record's `recommended_refit_dataset` for continued analysis, while retaining
`fit_datasets` as the original published-source provenance. The shared
refit/coverage scripts remain untouched in this scoped work.

## Commit readiness on 2026-10-16

The material, original and corrected numerical datasets, reproduction scripts,
source audit and targeted tests are ready for a scoped commit. All nine
material-specific Python tests, the Rust integration test, reproduction and
export checks, scoped linting and metadata/schema checks pass.

Shared source/refit ledgers, aggregate counts, coverage, documentation navigation
and the Unreleased changelog entry still need the combined material-integration
pass described above. These common files remain outside this commit's scope.
This is not a claim that the complete repository is ready to release or that
all shared CI checks pass. The archived full-check snapshot above is historical.
