# Catalog integration: 16 October 2026

This historical first pass reconciles the completed forsterite and diopside additions, the
rubidium-halide family, the existing combined Au refit, and the completed Sakai
Au reconstruction. It describes the local checkout, not a published release.
At that pass, brucite, Fe-ringwoodite and phase E were separately staged and
excluded from these counts. Subsequent phase-E and restricted brucite additions
are reflected in the current [coverage overview](../catalog-coverage.md); see
the [brucite handoff](brucite-integration.md) for its explicit qualifications.

## Reconciled inventory

- 230 material documents and 658 stored EOS records.
- 227 executable materials and 655 executable records.
- 648 primary-source-validated, seven qualified `not_reproduced`, and three deferred records.
- 63 materials with an executable thermal EOS.
- All 658 record identities represented exactly once in both scientific ledgers.
- Refit outcomes: 232 conditional or otherwise scoped parity, 159 similar,
  75 parity not achieved, 176 not directly refittable, 16 source reconstructions,
  and no failed attempts before comparison.

The [generated coverage overview](../catalog-coverage.md) and its JSON snapshot
are the maintained inventory. Dataset inventory wording now explicitly counts
all retained datasets/rows, including calculated checkpoints and diagnostics;
these totals are not counts of independent experimental observations.

## Scientific dispositions

| Contribution | Integrated assessment | Preserved qualification |
|---|---|---|
| Two Katsura (2022) forsterite thermal records | Conditional parameter parity within reported comparison widths | Full-precision author-script vector recommended; printed-table alternative and its errors separate; exact bootstrap and original-2009 combined fit unrecovered |
| Thompson (2008) diopside static EOS | Conditional parameter parity within reported error widths | Source confidence, exact author weights and covariance unspecified |
| Zhao (1998) diopside thermal EOS | Conditional parameter parity with graph-corroborated row-15 treatment | Corrected analysis dataset is separate from the original transcription; 4.58 GPa graphical reading and +/-0.05 GPa non-statistical bound are not recovered original inputs |
| Farla (2025) RbCl/RbBr/RbI thermal EOS | Three conditional numerical-similarity outcomes | Original coefficient-rounding/workflow parity false; independent native EosFit replay, benchmark discrepancies and composite-CsCl calibration qualifications remain explicit |
| Combined Hirose/Fei Au refit on Speziale MgO | Registered Peritheos fit reproduced | Not recovery of Hirose/Fei author coefficients or their historical MgO reducer; predictive/calibrant/interstage uncertainty remains unpropagated |
| Sakai (2025) Au stage | Conditional K0-prime central-value parity at published three-decimal precision, 210 selected states | Other coefficients fixed; uncertainty and full sequential-calibration parity false; Table S9 remains calculated checkpoints |

Only Sakai Au changes an existing common-ledger outcome, from unavailable
full refit to the newly documented conditional one-parameter replay. No existing
published coefficients, thermal conventions, defaults, or source-observation
bytes are changed by this integration.

The missing Au `primary_source_check` is now populated in the record and its
registration routine. The refit report's implementation checksum is refreshed
for that metadata-only code edit; the derived CSV is byte-identical.

## Interchange and maintenance

The Rust JSON reader now enables `float_roundtrip`. The previous reader changed
one small Au coefficient by a floating-point rounding step on repeated JSON
loading. A focused bit-preservation regression and the complete bundled-record
round-trip test protect the fix; source decimals are not rounded as a workaround.

The source-audit manifest, numerical refit and paper ledgers, source-audit index,
coverage JSON/Markdown, Hugoniot inventory, documentation counts and current
backlog are reconciled. Historical release counts remain historical. The API
reference now includes the existing variable-q Debye class. Optional source-file
requests retain their no-outreach status.

Fei Pt remains paused pending an author reply and explicit user resumption.
Yokoo pressure evaluation remains on its evidence-dependent hold. No speculative
CsCl/Pt/Mo calibration link or substitute equation is introduced to clear a
missing-reference status.

## Verification

- All 171 shared catalog/ledger/coverage tests pass; the expanded final regression
  selection passes all 551 tests, including every previously failing module.
- Full record-by-record ledger replay passes `--check`; coverage, source-index
  and paper-ledger freshness checks also pass.
- Full Rust workspace tests, Rust documentation tests, formatting and Clippy pass.
- Strict documentation build passes with notebook execution enabled.
- An isolated wheel install passes the catalog/resource smoke check: 230
  documents, 658 stored records, 227 executable materials, 655 executable records.
- Python lint and mypy pass. Formatting checks for integration-owned files pass;
  separately staged material work is verified in its owning chats.
- The full Python run completed with 2,779 passes and ten failures observed
  while shared files were being reconciled. All ten failures were corrected and
  their modules rerun in a final **551-passing-test** regression selection. The ledger
  check keeps published values, source text, counts and classifications exact,
  while allowing tightly bounded roundoff in calculated optimizer diagnostics.
  Fresh complete replay produced identical rendered Markdown; only numerical
  optimizer details in the Chen argon and Mosenfelder PPv diagnostics differed.
- Dedicated checks retain original published Au equations, coefficients,
  thermal settings and defaults; the combined-refit reproduction check passes.

This integration does not publish a package, push changes, or contact authors.
