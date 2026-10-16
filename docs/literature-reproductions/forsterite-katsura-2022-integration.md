# Forsterite G01 integration note

2026-10-15: local material-specific work only. No commit, push, publication,
author contact, or message to another chat was performed.

## Completed addition

- Added recommended thermal record
  `forsterite_katsura_2022_author_script_bm3_mgd`, with the full-precision
  geotherm-script coefficients and theta0=760 K.
- Retained `forsterite_katsura_2022_bm3_mgd` as the printed Table 2
  alternative, with theta0=768 K.
- Preserved the entire Finkelstein record and existing static default.
- Linked 51 recovered P–V/V0–T rows with paired MgO observations and errors.
- Archived the original CC BY 4.0 author workbook and relevant MATLAB files.
- Added complete Olivine worksheet CSV and formula/cached-value JSON evidence,
  deterministic recovery/reproduction scripts, and Python/Rust checks.
- Reconstructed a candidate six-column fitting-input CSV from AD/AE/AN/AO/U,
  with sigmaT explicitly missing. The original export's exact pressure
  selection remains unconfirmed; no bootstrap temperature error is invented.
- Recovered all 52 Table 1 rows from the supplied original 2009 PDF,
  including reported lattice, volume and pressure errors. Fifty match the
  later workbook, two are ambient references, and M375006 appears only in
  the workbook. Added the table as a standalone source-comparison dataset.
- Reran the reconstructed-input fit with both the independent optimizer and
  a translation of the author's deterministic optimizer. Both agree within
  1e-5 in all coefficients. The all-51 main fit and 50-row table-matched
  sensitivity fit remain diagnostics; source coefficients are unchanged.
- Classified the recommended script vector's fit reproduction as **`parity`**:
  conditional numerical parameter agreement within printed Table 2 error
  widths (`parity_basis: within_reported_parameter_errors`). Script-specific
  errors remain missing; the exact bootstrap workflow is unverified.
- The 2009 parameterization remains on hold. The deposited executable-source
  vector is recommended. Ordinary refits are diagnostics only; script-specific
  errors and covariance remain missing.

## Global ledger integration completed

On 2026-10-15 the two 2022 records were added to the global primary-source
audit and primary EOS refit ledgers. Both are graded `parity` on the explicit
`within_reported_parameter_errors` basis. The full-precision author-script
record remains the recommended thermal option; no coefficients were replaced
with refit values. The paper investigation ledger marks Katsura (2022)
reproduced with two parity records and documents the comparison-width and
bootstrap qualifications. The original Katsura (2009) article is registered
separately as direct refit unavailable pending the combined-fit inputs and
normalization audit. The corresponding G01 update remains in the shared
backlog working copy, alongside its broader queue rewrite.

The updates preserve every unrelated source-audit and refit entry. The shared
refit coverage check still identifies six unrelated additions absent from
both ledgers: the Au combined refit, two diopside records and three Farla
rubidium-halide records. These require their owners' batch reconciliation.

Reproduce this scoped ledger refresh with:

```console
.venv/bin/python -m scripts.forsterite_reproduce_katsura_2022 --write-report --update-ledgers
```

## Remaining shared updates for batch reconciliation

1. Reconcile the material manifest/index with all parallel additions.
2. Keep the script vector recommended and the printed-table errors attached
   only to their own record. Original author bootstrap reproduction remains
   unresolved; do not describe this as original-2009 fit reproduction.
3. Link the source audit in the source-audit index and docs navigation if
   appropriate. Both global EOS and paper ledgers are now updated.
4. Original Katsura (2009) remains a combined-fit input/normalization hold.
   The full text and original Table 1 have been recovered; do not keep
   those access gaps open.
5. Regenerate the common coverage snapshot and reconcile catalog count
   assertions after all parallel additions. This work adds **two EOS records,
   two thermal records, two primary-source-validated records, zero materials**.
   Existing Au expectation failures must be handled by their owning work.
6. Add a dataset-request entry if tracking the remaining original 2009
   Takeuchi thermal-expansion observations/fit selection, and separate 2022
   `Olivine_PVT2021.xlsx`/temperature errors is useful. Contact status:
   **not contacted**. No author outreach was authorized or attempted.
7. Pressure calibration: the record retains source-reported two-branch MgO
   assignments and does not require an invented shared calibration entry.
   A future executable combined reduction must audit both exact Tange (2009)
   branches and their mixing/error convention first. Do not point this record
   solely at `mgo_b1_tange_2009_vinet` as if it were the complete reduction.

## Verification scope

Material-specific Python and Rust tests cover the source mapping, structure
contents, native pressure evaluation, inverse/derivative/reference identities,
schema/interchange, dataset checksums/transcription, all-state author-code
comparison, and explicit experimental/ordinary-fit diagnostics. Global checks
and their current limitations are recorded at the end of the source audit.
The common refit validator currently stops at the existing combined Au refit's
missing `primary_source_check`; the Au owner must resolve that before a shared
ledger refresh can complete.
