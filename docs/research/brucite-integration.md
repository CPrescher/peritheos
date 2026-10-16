# Brucite G04/C12 catalog integration

2026-10-16. **Published-range lattice-average volume EOSs installed.** This
supersedes the earlier restricted integration. No commit, push, publication or
author contact occurred.

## Curved expansion with adopted published moduli

The user-requested purple-law variant is now installed as five additional,
nondefault `record_kind=refit` records. Four source-temperature isotherms
retain the exact Table 5/6 K0/Kprime values and measured pressure spans. A
separate continuous version uses the same curved expansion, the printed
modulus slopes and explicit Table 5/6 300 K anchors. The continuous version
is not represented as passing through every tabulated modulus.

The selectable continuous identifier is
`brucite_peritheos_fukui_2003_curved_expansion_thermal_refit`, valid at
300–873 K and 5.95–15.48 GPa. Isotherm variants use
`brucite_peritheos_fukui_2003_curved_expansion_{T}k_bm3_refit`.
All three expansion coefficients come from the same unweighted fit to all
11 ambient observations. The full shared covariance and its per-isotherm
propagation are retained; adopted source moduli/slopes and calibration
systematics are outside this conditional prediction covariance.

The original seven records were compared directly with a pre-extension
snapshot and are unchanged, including the published default. Current coverage
is 232 material documents / 671 stored records and 229 executable materials /
668 executable records (659 source-validated, nine qualified, three deferred).
The [curved-expansion report](../data/brucite-curved-expansion.json) distinguishes
ambient fitting from the independent 37-row pressure validation and the
unrecovered original-author workflow. The [full audit](../literature-reproductions/brucite.md#curved-expansion-catalog-refits-with-published-moduli)
provides selection examples and qualifications.

A public uncertainty check exposed a pre-existing finite-difference scale
problem: a unit-sized floor perturbed tiny expansivity-curvature coefficients
nonlocally. The scale now follows the coefficient magnitude or its standard
error, with a fallback only when both are zero. Analytic small-curvature
regressions at zero and nonzero coefficient and correlated-covariance checks
verify the correction.

Extension checks: **632 focused Python tests pass; all 157 workspace Rust
tests pass**. Global Ruff lint/format, mypy (56 library files), Rust formatting
and all-target/all-feature Clippy pass. The isolated built-wheel install
verifies all five new models, source defaults, covariance propagation,
interchange, all 102 dataset rows, CIF inclusion and final-card byte equality.
The normal full Python suite was not rerun for this extension; related
uncertainty/material/catalog/ledger tests are included in the focused run.

Strict documentation build was attempted: seven warnings are unrelated links
in the concurrent Brown-NaCl and Fe-ringwoodite research notes, with no brucite
warning. Those files were preserved. The final whole-catalog numerical replay
(`validate_primary_eos_refits.py --check`) passes for all **671 records**.

## Previous published-range integration

The 0.5 GPa limit was an unsupported audit restriction, not a published bulk
lattice boundary. Fukui's higher-order hydrogen rearrangement and later
Sun (2025) XRD support an empirical layered lattice-average volume EOS through
the source's measured range. H ordering remains annotated; ambient coordinates
are a reference only. No phase-equilibrium or dehydration calculation is implied.

The installed card retains published 300/473/673 K isotherms through
22.56/16.52/16.09 GPa, respectively. The published 873 K curve is preserved
through 15.48 GPa as qualified explicit selection. Its original coefficients
are unchanged. Fig. 6(c)'s plotted markers imply a normalization near 42.70 A^3,
while the printed expansion gives 42.516 A^3; the solid line agrees with the
printed moduli. Exact author intent is not recovered and no source typo is asserted.

Two separately attributed Peritheos unweighted fits make the source tables
usable: an eight-row 873 K BM3 (5.21–15.48 GPa, RMS 0.158 GPa) and a 37-row
continuous thermal fit (5.95–15.48 GPa, 300–873 K, RMS 0.310 GPa).
The latter uses the common pressure interval of all four slices and preserves
linear K0/Kprime temperature laws. Both include full conditional covariance,
fixed V0/alpha assumptions, objective and selection. Published coefficients
remain in separate records. The original source-slope reconstruction also
remains available with `not_reproduced` and explicit-selection qualifications.

The active ambient structure is the directly inspected fluorine-free hydrogenated
neutron refinement in He (2023), Table 4.3, DOI
[10.15083/0002012205](https://doi.org/10.15083/0002012205). The Zigan archive
is retained separately. All 102 Fukui/Xia observations, errors and original
source tokens remain intact. Paired Au readings, original weights, full precision
and author covariance are still unavailable; these do not prevent a clearly
attributed independent table-data regression.

The [full audit](../literature-reproductions/brucite.md) documents source access,
phase evidence, graphical coordinates, numerical results, rights and API usage.
The [coverage report](../catalog-coverage.md) records the current inventory.
The two new refits use the ledger's `registered_peritheos_refit_reproduced`
basis; it does not mean original-author coefficient or workflow parity.

Verification for this extension:

- **559 focused Python tests pass**, covering brucite, catalog/discovery, material
  loading/export, scientific ledgers and documentation consistency.
- **Two native Rust brucite tests pass**, including all 873 K observations and
  continuous-model evaluation/inversion.
- Independent finite-difference uncertainty propagation agrees with the full
  2×2/4×4 covariance matrices. Thermal covariance uses the runtime reference-EOS
  parameter paths and survives interchange.
- Primary-PDF checksum/token checks, deterministic table/report/card builders,
  source index, paper ledger and coverage checks pass. Ruff lint/format,
  Rust formatting and whitespace checks pass.
- An isolated built-wheel installation checks native thermal evaluation,
  uncertainty propagation, all 102 dataset rows, CIF inclusion, qualified
  loading and byte-for-byte equality with the final catalog card.
- Strict documentation build was attempted; only the same two unrelated
  Fe-ringwoodite image links outside the docs tree caused failure. No brucite
  warning was reported; that other task's files were preserved.

The final whole-catalog numerical replay (`validate_primary_eos_refits.py --check`)
passes for all **666** records. All seven brucite ledger outcomes also match
freshly regenerated source/fit evidence.
The preceding narrower integration had a full Python/Rust suite run; those
results are historical and are not presented as a full-suite rerun of this extension.

Logs are `/tmp/peritheos-brucite-broad-tests-final.log`,
`/tmp/peritheos-brucite-broad-ledger-check.log` and
`/tmp/peritheos-brucite-broad-docs.log`.
