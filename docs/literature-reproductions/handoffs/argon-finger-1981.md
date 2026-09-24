# Finger 1981 integration handoff

Implementation commit: **d7870cc1b3f83a5a4b13b58b92beecbe455baa3e**.
Base: 664f97306241bdfde4b0b13cf9bd6b584dd03fe3.
Worktree: `/Users/clemens/.codex/worktrees/0423/peritheos`.
This handoff is a subsequent documentation-only commit; cherry-pick both.
No push, deployment, shared-checkout edit or Zotero write was performed here.

## Records and exact model

- Material: `argon_fcc`; appended record `argon_fcc_finger_1981_murnaghan2_debye`.
- Dataset: `argon_fcc_finger_1981_table1`; 19 measured fcc Ar points, 293 ± 1 K,
  1.28–8.17 GPa. All raw Table I columns retained in
  `peritheos/data/datasets/argon-fcc-finger-1981-table1.csv`.
- New Python/Rust/EOSMAT model: `SecondOrderMurnaghan` / `second_order_murnaghan`.
  Parameters in conventional-cell/GPa units: V0=149.82711895296183, K0=2.3701,
  K0_prime=6.97, K0_double_prime=-0.40 GPa^-1, P0=-0.10289 GPa.
- Thermal model: existing `Dewaele2006` / `dewaele_2006`, exact source law with
  Tr=293, theta0=93.3, gamma0=2.7, gamma_inf=.5, beta=1, n=1; anharmonic_a=0,
  electronic_e=0, inactive exponents=1. Required configuration:
  `thermal_pressure_reference="absolute_zero"`, **`zero_point_pressure="included"`**.
  The new zero-point option defaults to omitted, preserving previous models.
- Python/Rust registries, EOSMAT schema/factories and native fitter support are
  included. Native fitter now also preserves the previously available
  absolute-zero configuration during parameter reconstruction.
- Structure source URL repaired. No edit to Wittlinger hcp or its
  `not_reproduced` status; existing Dewaele/Ono coefficients unchanged.

## Reproduction status

Published equation is independently reproduced, not a BM2 or Vinet surrogate.
`python -m scripts.reproduce_argon_finger1981` regenerates
`docs/data/argon-finger-1981-reproduction.json`; detailed scientific audit is
`docs/literature-reproductions/argon-fcc-finger-1981.md`.

Independent adaptive Debye quadrature versus native pressure: max 3.10e-11 GPa.
Published RMS residual: 0.095355 GPa using printed molar volumes, 0.097220 GPa
using a^3 cell volumes. Three explicit diagnostic objectives yield different
K0'/K0'' estimates. The unweighted refit reproduces both fitted parameters within
the published error bars (0.60 and 0.63 quoted error widths from the source).
This is consistent with an unweighted original fit, but its weighting remains
unconfirmed. Updated status:
`fit_reproduction_status="parameters_reproduced_within_reported_uncertainties"`.
This means parameter-value agreement, not recovery of the unpublished original
weights, covariance or parameter-uncertainty estimates.
Printed Pz0 and zero-point pressure leave +0.001570181 GPa at V0 and 0 K;
coefficients are preserved. Fig.2 non-293 K isotherms are theory, not observations.

## PDF and Zotero

Source manifest: `peritheos/data/datasets/argon-fcc-finger-1981-sources.json`.

- DOI: `10.1063/1.92597`.
- Full citation: L. W. Finger, R. M. Hazen, G. Zou, H. K. Mao, P. M. Bell,
  “Structure and compression of crystalline argon and neon at high pressure
  and room temperature,” Applied Physics Letters **39**(11), 892–894,
  1 December 1981.
- Verified author-hosted publisher PDF: 3 pages, 192142 bytes.
- Origin: https://hazen.carnegiescience.edu/sites/default/files/061-argon-1981.pdf
- Stable PDF: `/Users/clemens/Documents/Peritheos/argon-sources/finger-1981/finger-1981.pdf`.
- SHA256: `a23289268d0a191b894bbeec8920fd4f6ea926847a4a5bc68686aa4a496f8bbb`.
- All pages visually checked for title, authors, journal/year/pages, tables and
  equations. Text extraction alone gives only download watermarks.
- Fresh read-only Zotero lookup: item **2ECG7DTR**, attachment **ZNGG5R8C**.
  Attachment path and identical hash are recorded in the manifest. Initial DOI
  search was empty; item appeared by later `argon` search while integration was
  active. Coordinator must check collection C18 and avoid a duplicate import.
- No access blocker. PDF not committed because copyrighted.

## Studio and shared aggregate work

Rebuild the engine against the integrated Peritheos revision; the generic Rust
EOSMAT loader now executes the exact model. Carry both thermal configuration
fields through any Studio material/draft adapter. Do not map the new mechanism
to BM2 or ordinary first-order Murnaghan.

Map observation columns explicitly to avoid ambiguity between retained source
and derived units:

- pressure: `pressure_gpa` (GPa)
- pressure error: `pressure_uncertainty_kbar / 10` (GPa)
- volume: `volume_a3` (four-atom conventional fcc cell)
- volume error: `volume_uncertainty_a3`
- temperature: `temperature_k` (293 K)
- source molar volume and source a values remain inspectable metadata.

Treat quoted errors as source uncertainties of unspecified confidence level.
The ±1 K is the reported experiment temperature variation, not 19 independently
measured temperature errors. Render measured points separately from EOS curves.
The coexistence point at 1.15 ± .05 GPa has no volume and is documented, not
fabricated as an observation.

Shared aggregate manifest/audit/refit ledgers intentionally await coordinator
regeneration after all nine studies. This addition is +1 EOS / +1 validated
executable record / +1 thermal record, pressure calibration `resolved`,
recalculation `missing_calibrant_observations`. Source audit and fit ledger need
the new ID with honest equation-versus-original-fit distinction.

## Verification

- Focused argon suite: **9 passed**, including new model, independent equation,
  source checksums, inversion, fallback, configuration roundtrip and native fit.
- Broader Python suite: **572 passed, 6 deselected** across argon, fitting,
  materials, datasets, Dewaele iron and EOSMAT.
- Rust EOSMAT/thermal/isothermal/Finger suite: **67 passed, 1 filtered out**.
- `cargo check --workspace`, `cargo clippy --workspace --all-targets -- -D warnings`,
  changed-Python Ruff checks and `git diff --check` passed.
- The six deselected Python tests are aggregate record-count/audit-ledger checks
  in `test_materials.py` and `test_eosmat.py`. They were run and fail solely from
  617/615 actual versus old 616/614 and missing aggregate audit entry.
- Filtered Rust `all_bundled_material_records_load_and_round_trip_through_rust`
  was run: all models loaded/evaluated/round-tripped, then it failed at manifest
  total 617 versus 616. Refresh aggregates and rerun in the integrated tree.
- The normative schema test was updated from 17 to 18 isothermal mechanisms
  and passes. This expectation may need further increments if other studies add
  mechanisms.

Local native build used Python 3.9 from the existing environment (read-only),
`PYO3_PYTHON=/Users/clemens/Programming/peritheos/.venv/bin/python cargo rustc
-p peritheos-python --lib -- -C link-arg=-undefined -C link-arg=dynamic_lookup`,
then copied the resulting dylib only into this isolated worktree as the ignored
`peritheos/_rust.cpython-39-darwin.so`. No binary is committed.
