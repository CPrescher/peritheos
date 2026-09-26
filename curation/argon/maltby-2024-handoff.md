## Executable unvalidated option (2026-09-26)

User requested a selectable EOS despite unresolved reproduction. Record
`argon_fcc_maltby_2024_published` uses fixed published coefficients, explicit
assumed squared cutoff 64, and geometric fcc shells. It is nondefault and marked
`not_reproduced`; the diagnostic refit is not used. Python and native EOSMAT
execute it, retaining the Table 8 discrepancy and exported qualification.
The source-law characteristic volume is not a zero-pressure volume. Dewaele
remains the default. See the updated reproduction documentation.

# Maltby 2024 integration handoff

Implementation commit: `37ac5b490c2894f7e056087547ef64a21016b923`.
Worktree: `/Users/clemens/.codex/worktrees/9877/peritheos`.
Coordinator: `01a0d259-4995-7a63-a9e4-af02c28c36cc`.

## Outcome and identifiers

- Study identifier: `argon_fcc_maltby_2024`.
- Status: `independent_eos_candidate_not_reproduced`; scientific validation
  `not_reproduced`. This is an independent EOS paper, not measurement-only support.
- No new executable EOSMAT identifiers, schema, dispatcher, material mutation,
  experimental-observation rows, or Studio files. Existing fcc/hcp records stay intact.
- Experimental scalar APIs: Python `peritheos.eos.experimental.Maltby2024Published`
  and Rust `peritheos::experimental::maltby_2024::Maltby2024Published`.
  Both require an explicit occupied fcc `shell_cutoff_squared` (1--256).
- Pressure, unshifted Helmholtz energy, bulk modulus, stable volume inversion,
  and 300 K pressure increments are implemented. Molar volume is J/bar/mol;
  pressure is GPa. Absolute energy/entropy reference offsets are excluded.

The independent SI quadrature agrees with the analytic implementations, but
Table 8's 70 K, 1 MPa volume is 23.97 cm3/mol versus 23.89042753 from the literal
coefficients with cutoff squared 64. The 0.332% discrepancy persists at larger
cutoffs. The paper does not specify its numerical cutoff. No coefficients were
tuned to repair the discrepancy; no validated executable record should be added.

## Source and comparison assets

- `curation/argon/maltby-2024.json`: DOI, coefficient transcription, Table 2
  primary-source counts, theory-only Table 8 values, PDF hashes and Zotero state.
- `docs/literature-reproductions/argon-maltby-2024.md`: source equations, phase,
  units, zero-point/reference choices, uncertainty, pressure scales and limitations.
- `docs/data/argon-maltby-2024-reproduction.json`: independent quadrature,
  cutoff sensitivity, Table 8 mismatch and existing Dewaele comparison metrics.
- `scripts/reproduce_maltby_2024_argon.py`: repeatable independent quadrature and
  diagnostics. The Dewaele comparison retains 130 original rows within 16 GPa,
  not the paper's 38-row primary fit selection. RMS pressure residual 0.48719 GPa;
  mean absolute relative deviation 4.82986%. No full primary-data refit claimed.

## Actual PDFs and Zotero

The following are verified PDF files outside Git, with retrieval provenance in
`/Users/clemens/Documents/Argon-EOS-Papers/maltby-2024/provenance.json`:

1. `/Users/clemens/Documents/Argon-EOS-Papers/maltby-2024/author-manuscript.pdf`
   — NTNU/NVA accepted manuscript, 23 pages, 1,292,435 bytes; exact publication
   DOI `10.1063/5.0237497`; SHA256
   `da12e25c1cfd36917876fc138dcb71b197df3c2ee5c512b56b11c4684214bc4d`.
2. `/Users/clemens/Documents/Argon-EOS-Papers/maltby-2024/supplement.pdf`
   — official AIP Figshare supplement, 215,281 bytes; DOI
   `10.60893/figshare.jpr.27846567.v1`; SHA256
   `08986ed607d464ed16565b6007c0f83ce080f326daa84ddd5260b5c92e086229`.

The coordinator reports the manuscript imported and exact-hash verified in
My Library / Methods / EOS Library (C18): article `CGV58DF9`, attachment
`IDGZPNCM`. The coordinator is attaching the supplement. This task made no
Zotero writes. The main PDF is an accepted manuscript, not the publisher version.

## Verification and integration actions

Passed:

- Python reproduction script, regenerating the committed JSON report.
- `pytest -q tests/test_maltby_2024_argon.py tests/test_argon.py`: 30 passed.
- `cargo test -p peritheos experimental::maltby_2024 --lib`: 2 passed.
- `cargo clippy -p peritheos --lib -- -D warnings`.
- Ruff lint and format checks for all four new Python files; `git diff --check`.

These checks validate numerical implementation and preserve the source mismatch;
they do not establish reproduction of the published EOS. Existing compiled Python
extensions and the pytest interpreter were read from the shared checkout only;
all edits/build output were confined to this worktree or the PDF handoff directory.

Coordinator should cherry-pick the implementation commit and this handoff commit,
finish/record supplement attachment verification, and integrate a pending study
entry into Studio if supported. Display source PDFs and the unresolved status;
do not expose a selectable validated pressure curve or label Table 8 as measured
observations. Shared-checkout landing, Studio integration, and Zotero import
remain coordinator-owned. No push was performed.
