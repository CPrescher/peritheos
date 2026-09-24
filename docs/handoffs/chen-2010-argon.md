# Chen 2010 argon integration handoff

Study ID: `argon_fcc_chen_2010`.
DOI: `10.1103/PhysRevB.81.144110`.
Outcome: **supporting study, `not_reproduced`**, with no executable Chen EOS.
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
