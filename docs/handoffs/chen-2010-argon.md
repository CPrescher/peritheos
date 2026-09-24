# Chen 2010 argon integration handoff

Study ID: `argon_fcc_chen_2010`.
DOI: `10.1103/PhysRevB.81.144110`.
Original-source outcome: **supporting study, `not_reproduced`**. A later,
explicitly requested independent Peritheos refit is described in the addendum
below; it does not turn the source-author parameterization into a validated EOS.
Material: existing `argon_fcc`, conventional cell Z = 4. Wittlinger hcp is
untouched and remains `not_reproduced`.

## Commits to integrate

1. Prerequisite from the Ross worker: `beb665cee539cc8745e5aa38179a6dc726e1a615`
   (standalone datasets, Python/schema/material selection). Its local cherry-pick
   here is `e2b1904`; do not apply both.
2. Rust parity: `e6206c2990e89708aa7b2634f425d1f17e208ec1`.
3. Chen source/data/audit: `b61249e3ec6d08617b159dadbddb145dbc52b8bb`.
4. Alternative pressure-reference audit: `de1f7ba39265a3a200ac9c6d08964be54f19cbc3` (also formats two lines of
   the prerequisite Python changes; no behavior change there).
5. This handoff document is a final separate commit.

Merge `argon_fcc.eosmat` by appending the three Chen datasets, preserving all
other workers' additions. Regenerate the global source-audit index and dataset
inventory counts after integration; this worker's counts describe its isolated
worktree. No shared checkout, Studio checkout, Zotero library, remote branch
or deployment was modified.

## Artifacts and Studio requirements

Packaged metadata: `peritheos/data/studies/argon-fcc-chen-2010.json`. It contains
the complete reference, Zotero-compatible `bibliographic_metadata`, published
constants, range/reference state, calibration, PDF manifest, limitations and
Studio guidance. Reproduction report:
`docs/data/chen-2010-argon-reproduction.json`. Source audit:
`docs/literature-reproductions/chen-2010-argon.md`.

Datasets, all with empty `used_by_eos_records`:

- `argon_fcc_chen_2010_figure5_brillouin`: 80 distinct black-square PDF positions
  with matched graphical P/rho errorbar halfwidths, source density, converted
  four-atom cell volume, PDF coordinates and rendering multiplicity. These are
  **Brillouin-integrated densities**, not direct XRD measurements. Duplicate
  rendering commands do not represent extra measurements.
- `argon_fcc_chen_2010_figure5_curve`: 262 vertices of the published fitted line.
  Kind `published_fit_curve`; plot as a line explicitly labeled published fit.
  Never use these as measured points or a first-principles theory dataset.
- `argon_fcc_chen_2010_reported_values`: 14 printed constants from the own-study
  table row and prose/captions. Row-specific units and uncertainty conventions
  must remain visible; missing uncertainty is not zero.

All resources live in `peritheos/data/datasets/argon-fcc-chen-2010-*.csv` with
SHA256 checksums in the material card. Studio should expose this study in its
supporting-study browser and permit viewing the density points and published
line separately. Do not fabricate a `fit_datasets` link to Dewaele or Ono.
Use 290 K from Table I while showing the room-temperature/300 K prose caveat.
Ruby scale is Mao et al. **1978**, source reference 14, not the paper's
reference 13 (Mao 1986). Raw ruby observations are unavailable.

## PDF and Zotero import

Full paper: *Elasticity, strength, and refractive index of argon at high
pressures*, Bin Chen, A. E. Gleason, J. Y. Yan, K. J. Koski, Simon Clark and
Raymond Jeanloz, Physical Review B 81(14), 144110, 14 April 2010.

- PDF absolute path:
  `/Users/clemens/Documents/Peritheos/papers/chen-2010/chen-2010.pdf`
- Origin URL:
  `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.81.144110/fulltext`
- SHA256:
  `17a0bcd30604bd81f125894a18a45a91255b235fb80643cddd4bcbc6e5bd01f1`
- Verified publisher version, all five pages 144110-1 through 144110-5; title,
  authors, DOI and journal identifiers match. Retrieved 2026-09-24.
- RIS:
  `/Users/clemens/Documents/Peritheos/papers/chen-2010/chen-2010.ris`
- Existing Zotero item/attachment keys: **none found** by title search through
  the Zotero skill helper. Coordinator should check DOI for duplicates at import.
- Destination: **My Library > Methods > EOS Library**, collection **C18**.
- No PDF access blocker remains. The PDF is not committed to git. Import the
  actual local PDF attachment, not just a linked URL or bibliographic record.

## Reproduction status and limitations

The true nonzero-pressure local constraints at 2 GPa are rho = 2.18, KT = 15.1
and Kprime = 5.4. Standard BM3 reconstructed from them gives rho0 = 1.674647
versus printed 1.52 ± 0.05, and 6.880 GPa RMS residual against drawn densities.
The alternative P = 2 + BM3 convention reproduces printed KKsecond (−7.248889
versus −7.3) but still gives rho0 = 1.633737 and 6.251 GPa RMS residual.
The drawn curve itself has rho(2) = 2.18017 but local KT(2) about 11.157 GPa.
Both interpretations are diagnostic and must remain excluded from the
executable catalog and reproduced-fit totals.

The paper does not print an unrounded EOS coefficient set or explicit BM
pressure equation. Original velocity/density tables, interpolation/integration
details, fitting weights, covariance and raw calibrant observations are missing.
Eq. (4) and the published Cp/Cv expression are implemented in the audit script;
no synthetic velocity curve or substitute EOS was invented. No purchase or
author contact was made.

## Validation

- 166 focused Python tests passed (Chen, existing argon, EOSMAT, datasets and
  supporting-dataset behavior), with warnings treated as errors.
- Six documentation consistency checks passed after inventory/index update.
- Four Chen tests passed again after the alternative-reference diagnostic.
- `cargo test -p peritheos`: all Rust unit, integration and doc tests passed,
  including 32 EOSMAT tests and the new empty-link regression.
- Ruff lint and formatting across all 581 Python files passed; `cargo fmt`
  and `git diff --check` passed.
- PDF extraction rerun reproduced 80 marker rows and 262 curve vertices;
  the audit regenerates from bundled CSV without a PDF dependency.
- The optional full Python sweep was attempted with both debug and release
  native builds, then stopped while running unrelated long literature refits.
  It is **not reported as passed**. The completed focused checks above cover
  this study, the shared dataset semantics, all EOSMAT tests and the existing
  argon regressions. Full-suite verification can be performed once on the
  integrated nine-study result.

## User-requested independent refit addendum

After reviewing the source discrepancy, the user explicitly agreed to trying
our own separate digitized-data refit. This adds one executable, nondefault
`record_kind: refit` record: `argon_fcc_chen_2010_bm3_digitized_refit`. Its BM3
model uses the existing `birch_murnaghan_3` adapter. The source study remains
`not_reproduced` and its `eos_record_identifiers` stays empty; the new record is
listed separately in `peritheos_refit_record_identifiers`.

Primary regression: all 80 black-square positions, P/rho errors-in-variables
with graphical coordinate halfwidths as relative scales, three free EOS
parameters and 80 latent densities. No PDF rendering multiplicities, source
curve vertices, printed moduli or extrapolated density constrain the fit.
Use the unchanged source points for plotted measurements; latent coordinates
are diagnostics in the JSON report only.

- V0 = 157.054399126 Å³/four-atom fcc cell, K0 = 4.398863939 GPa,
  K0prime = 4.621762004; zero-pressure extrapolated parameters.
- Original-coordinate pressure RMS = 0.111199 GPa, maximum = 0.209913 GPa.
- Weighting/bin-spacing sensitivity: at most 0.617% in-range volume change.
- Common digitization offsets: at most 0.522% volume change.
- Omit each 1 GPa pressure bin: at most 1.673% volume change and worst
  withheld-bin pressure RMS 0.271 GPa. K0/V0 extrapolations are less stable.
- Qualified range = 1.23188372–26.06236454 GPa at nominal 290 K.
- Parameter errors/covariance are explicitly null. Unknown integration
  covariance and unspecified graphical confidence levels preclude statistical
  uncertainty claims. Sensitivity ranges must not be rendered as confidence
  intervals.

Reproduce: `python -m scripts.refit_chen_2010_argon --plot`. Report and figure:
`docs/data/chen-2010-argon-refit.json` and `.png`. Detailed provenance and bounds
are embedded in `fit_provenance`; the source audit has a new refit section.
Tests: `tests/test_chen_2010_argon_refit.py` plus updated original Chen tests.

Studio should expose the record with its **Peritheos refit** label and fitted
source dataset, keeping the original source discrepancy and published-curve
comparison separate. The 80-point dataset now links to this refit; the other
two Chen datasets remain unlinked. No new equation class is needed.

Per the coordinator's request, global manifests, counts, primary-source/refit
ledgers and the final Studio source pin are coordinated centrally. This adds
one experimental validated *Peritheos refit* record, not one reproduced
source-author fit, and one additional dataset-to-EOS link. The new refit must
not be sent through the old Chen source-parameter failure diagnostic as though
it were a published coefficient record.

Refit validation: 27 focused Python tests passed; the added positive-bulk
modulus and inversion check was rerun separately. Ruff and diff checks passed.
The Rust EOSMAT run passed 31 tests; its sole failure is the centrally owned
hard-coded inventory assertion (617 actual versus 616 expected in this
isolated checkout, `crates/peritheos/tests/eosmat.rs:995`). All records loaded
and round-tripped before that assertion. Update the integrated total rather
than copying this worktree's count. The original completed Rust run above
predates this additional refit record.
