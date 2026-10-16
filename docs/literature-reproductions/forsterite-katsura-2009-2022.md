# Forsterite: Katsura (2009) source hold and Katsura (2022) thermal EOS

Audit date: 2026-10-15.

## Disposition

The recommended thermal record is
`forsterite_katsura_2022_author_script_bm3_mgd`, using the full-precision
olivine vector embedded in the deposited geotherm script. The printed-table
variant is retained as `forsterite_katsura_2022_bm3_mgd`. The existing
`forsterite_finkelstein_2014_bm3_1` and its default selection are preserved.
Both 2022 records describe Mg-endmember olivine with corrected temperatures
and recalculated MgO pressures. Neither is the original 2009 parameterization
or a new Peritheos refit.

## Recommended author-script parameterization

Recommended evaluation follows the coefficient vector embedded in the
deposited executable source. The printed table is preserved separately for
provenance. The geotherm script's actual coefficient vector is used, rather
than the fitting script's temporary optimizer starting guesses:

| Quantity | Author-script value |
|---|---:|
| V0 (A³/cell) | 289.58 |
| K0 (GPa) | 127.4 |
| K0 prime | 4.21116896988574 |
| theta0 (K) | 760 |
| gamma0 | 0.997258584556838 |
| q | 2.43854091371459 |
| Tr (K), n, Z | 300, 7, 4 |

All values are the first (olivine) elements of the named parameter arrays
in `Katsura Geotherm calculation v2.m`. The independent translation retains
that script's NA=6.02214076e23, R=8.314472 and Debye polynomial. Its pressure
differs from exact Peritheos Debye evaluation by at most **0.0000575 GPa**
over all 51 observations. The author expression gives **12.180116341 GPa**
at V/V0=0.95, T=1700 K.

The script vector's pressure RMSE is **0.2130082 GPa**. Its maximum pressure
difference from the printed-table alternative is **0.0192147 GPa** across
the recovered states. The outlying M375006 observation remains: script
calculated minus observed pressure is **+0.9723329 GPa**, or **6.34** times
the incomplete projected pressure/volume spread. The source choice does not
resolve that observation mismatch.

No script-specific errors or covariance are supplied. They remain missing
on this record; the rounded Table 2 bootstrap uncertainties remain attached
to the separate printed-table record. This avoids assuming that uncertainty
estimates made with theta0=768 K belong unchanged to the 760 K script vector.

## Candidate reconstruction of the fitting input

`forsterite-katsura-2022-candidate-fitting-input.csv` follows the six-column
order read by `xlsread` in the author fitting code:

| Input column | Deposited Olivine column | Recovery |
|---|---|---|
| P | AD | Saved mean pressure |
| sigmaP | AE | Saved branch/bracket sample spread |
| V/V0 | AN | Sample volume ratio |
| sigma(V/V0) | AO | Saved ratio uncertainty |
| T | U | Saved corrected mean temperature |
| sigmaT | absent | Empty, not zero |

All 51 rows survive in their source order. This is an explicitly labeled
candidate reconstruction. The separate original `Olivine_PVT2021.xlsx` is
still absent, so AD/AE versus the individual pressure branches cannot be
confirmed as the exact author's export choice. The full transcription retains
all pressure branches for later comparison.

S and T in the workbook are temperatures corrected at pressure brackets.
Their half-span can be calculated, but the inspected paper and workbook do
not establish that it is the complete sigmaT read by the bootstrap script.
It is therefore not substituted for missing temperature uncertainty. The
candidate permits an ordinary least-squares fit now; it does not support an
unqualified reconstruction of the author's stochastic uncertainty analysis.

The supplied final 2009 PDF resolves the full-text/Table 1 access gap.
All 52 printed rows are now recovered, including lattice, volume and pressure
errors. Fifty match the later workbook; two are ambient references, while
M375006 is present only in the workbook. The original 2009 combined fit still
requires the Takeuchi et al. (1984) thermal-expansion observations and its
source-specific normalization/selection audit. No 2009 executable record or
independent refit record is installed by this follow-up.

## Primary sources and access

- Katsura et al. (2009), *Thermal expansion of forsterite at high pressures
  determined by in situ X-ray diffraction: The adiabatic geotherm in the upper
  mantle*, PEPI 174, 86–92, DOI [10.1016/j.pepi.2008.08.002](https://doi.org/10.1016/j.pepi.2008.08.002).
  Final published PDF supplied and inspected. Table 1 is on journal p.89
  (PDF p.4); all 52 rows are bundled separately with the PDF checksum in
  `forsterite-katsura-2009-source.json`. Sections 3.2–3.3 give adopted
  V0,300=290.79 A³, K0=127.4 GPa, K0 prime=4.2, theta0=760 K, and fitted
  gamma0=1.31(1), q=3.2(1). Section 3.3 additionally uses Takeuchi et al.
  (1984) ambient thermal-expansion data. The original volume normalization
  and pressure/temperature coordinates remain distinct from the 2022 fit.
- Katsura (2022), *A Revised Adiabatic Temperature Profile for the Mantle*,
  JGR Solid Earth 127, e2021JB023562, DOI [10.1029/2021JB023562](https://doi.org/10.1029/2021JB023562).
  The final published version is available from [EPub Bayreuth](https://epub.uni-bayreuth.de/id/eprint/6715/).
  Equation and table layout were visually checked on PDF pages 4–5.
- Author-owned workbook, [Zenodo 5644426](https://doi.org/10.5281/zenodo.5644426),
  `Katsura EMF-corrected PVT.xlsx`, Olivine A1:AO54. All 51 data rows,
  41 columns, 2,091 occupied cells and 663 formulas/cached values are retained.
- Author fitting and geotherm scripts, [Zenodo 5903286](https://doi.org/10.5281/zenodo.5903286).
  The two relevant original files are archived without alteration. The
  pressure and deterministic optimizer functions were translated and run;
  the full MATLAB bootstrap/geotherm programs were not executed, and no
  author's local file path was followed.
- Official Wiley supplements S1 and Table S1 were located, but direct
  downloads returned 403. Their contents are not claimed as inspected.

The final article's Data Availability statement explicitly links the two
Zenodo deposits. Sources, hashes and access limitations are recorded in
`peritheos/data/datasets/forsterite-katsura-2022-source.json`.

## Identity, normalization and printed-table parameters

Section 2.2 represents dry Mg-endmember olivine. This matches the existing
Mg2SiO4, Pbnm (No. 62) forsterite structure. Existing occupied sites give
Mg8Si4O16, hence Z=4 and 28 atoms per conventional cell. The thermal `n=7`
counts atoms per Mg2SiO4 formula unit, not atoms per cell. The EOS volume
289.58 A³ is a conventional-cell volume in the deposited author scripts.
Molar volume is `NA * Vcell * 1e-30 / 4` in m³/mol. The existing ambient
diffraction cell has volume about 290.10 A³ and comes from another sample;
it is not forced to match the thermal EOS reference volume.

| Quantity | Stored value / error | Authority and role |
|---|---|---|
| V0 | 289.58 A³ / missing | Both deposited scripts; fixed; no error invented |
| K0 | 127.4 GPa / missing | Final Table 2; fixed, adopted from Isaak et al. (1989) |
| K0 prime | 4.2 / 0.4 | Final Table 2; bootstrap mean and standard deviation |
| Tr | 300 K / missing | Section 2.2, Equations (7)–(8); fixed |
| theta0 | 768 K / missing | Final Table 2; fixed, adopted from Isaak et al. (1989) |
| gamma0 | 1.00 / 0.02 | Final Table 2; bootstrap result |
| q | 2.4 / 0.5 | Final Table 2; bootstrap result |
| n | 7 / missing | Equation (10) and deposited scripts; exact formula count |

Table errors are bootstrap standard deviations from 1,000 resampled and
perturbed datasets, not known covariance or 95% interval half-widths. No
covariance matrix is available. Adopted parameter errors not supplied in the
represented source are missing, not zero.

## Exact equation mapping and printed-table numerical checks

Equations (7)–(11) define a BM3 **300 K isotherm**, plus

`DeltaP = gamma(V)/Vm * [ED(V,T) - ED(V,300 K)]`,

with SI energies and molar volume (converted to GPa),

`gamma(V) = gamma0 * (V/V0)^q`,

`theta(V) = theta0 * exp[(gamma0 - gamma(V))/q]`, and

`ED = 9*n*R*T/(theta/T)^3 * integral_0^(theta/T) x^3/(exp(x)-1) dx`.

Both energy terms use the same volume-dependent theta. These expressions map
directly to `BM3` plus `MieGruneisenDebye`, with
`debye_temperature_law=integrated_gruneisen` and
`thermal_pressure_reference=reference_temperature`. No new equation class,
absolute-zero cold curve, or variable-exponent theta interpretation is used.

For the printed-table alternative, the independent reference translates the author's deposited `MGD_EOS_RJ`
and `ThermalEnergy` functions, retaining their sixth-order polynomial,
R=8.314472 and NA=6.0221367e23. Peritheos uses exact Debye quadrature and its
current physical constants. Across all 51 states the maximum pressure
difference is **0.0000611 GPa** (test tolerance 0.00007 GPa). At V/V0=0.95,
T=1700 K the independent author expression, configured with the final table
theta0, gives **12.195391775 GPa**. Python and Rust reproduce it within the
same approximation tolerance. This is equation validation against primary
author code, separate from agreement with experimental observations.

Two experimental checkpoints, retained with their source coordinates:

| Workbook row / run | V/V0 | Corrected T (K) | Saved mean P (GPa) | EOS P (GPa) | Projected spread (GPa) |
|---|---:|---:|---:|---:|---:|
| 4 / M371008 | 1.030011491 | 1302.162963 | 0.574524338 | 0.370759833 | 0.074675315 |
| 54 / M375135 | 0.944411233 | 1750.226221 | 13.692764825 | 13.274444221 | 0.214129202 |

Here the diagnostic projected spread is
`sqrt(AE^2 + (KT * AO / AN)^2)`: AE is the deposited pressure spread and AO is
the deposited sample-ratio sigma. These two residuals are below three times
that diagnostic spread. It excludes temperature uncertainty and coefficient
covariance, which are unavailable. It is not a confidence interval for the
complete EOS or a claim that every observation agrees within 3 sigma.

Across all 51 points the published curve has **0.2124994 GPa pressure RMSE**.
The largest normalized residual is **6.38** using the same incomplete spread
projection: run **M375006**, workbook row 14, has calculated minus observed
pressure **+0.9783662 GPa** at 1312.087488 K and saved P=3.735453111 GPa.
This discrepancy is retained; the numerical observations and published
parameters are not adjusted to reduce it. The supplied 2009 Table 1 does
not contain M375006. That supports a separate table-matched sensitivity
selection, but does not establish why the later workbook contains it or
whether the original 2022 fitting export included it.

## Recovered data and pressure calibration

The complete workbook is archived, and a CSV retains all A:AO columns with
their Excel coordinates. A separate readable CSV exposes the following:

| Original coordinates | Meaning |
|---|---|
| A | Run label |
| B | Uncorrected thermocouple temperature |
| D:E | Original Matsui pressure and saved sigma |
| F:H | MgO V/V0 and lower/upper brackets |
| S:T:U | Corrected temperatures at the pressure brackets and their mean |
| X:Y | Corrected-temperature Tange Vinet pressure and saved sigma |
| AB:AC | Corrected-temperature Tange BM3 pressure and saved sigma |
| AD:AE | Mean and sample standard deviation of V,W,Z,AA |
| AN:AO | Forsterite sample V/V0 and saved sigma |

The retained workbook formulas explicitly identify AD as
`AVERAGE(V,W,Z,AA)` and AE as `STDEV(V,W,Z,AA)`. AE includes branch differences
and bracket spread; it must not be treated as a single independent pressure
measurement's standard deviation. The F:H calibrant assignment is additionally
checked against the separate bundled Tange Fit3 Vinet MgO scale: using F at B
reproduces the author's central-pressure AG checkpoints within 0.00008 GPa.
F must not be used as the forsterite volume coordinate.

Section 2.2 specifies the Nishihara et al. (2020) WRe thermocouple EMF
correction and **two Tange et al. (2009) MgO EOSs**. Workbook labels agree.
The Zenodo description instead links the 2012 Tange DOI for its pressure
correction. That metadata conflict is not silently resolved by replacing
saved pressures. The exact combined two-branch reduction is not represented
by the single bundled Vinet calibration, so the record does not claim
observation-level pressure recalculation is ready.

All numerical source values survive unchanged, including negative near-ambient
pressure, repeated runs, source error columns and pressure branches. Formulas
and cached values are retained together, without recalculating the workbook.
All 52 original 2009 lattice/volume/pressure rows and their printed errors
are now bundled separately. Earlier-study combined-fit observations remain
unavailable; none is synthesized from the later EOS. No temperature sigma
is supplied in Table 1.

## Independent refits and author-workflow boundary

Section 2.2 minimizes unweighted squared pressure differences and reports
bootstrap means/standard deviations after resampling and perturbing data.
The deposited fitting script reads **a separate `Olivine_PVT2021.xlsx`** with
six columns (P, sigmaP, V/V0, sigmaV/V0, T, sigmaT). That export, its separate
sigmaT values, and the bootstrap seed/realizations are not deposited.

The audit therefore performs ordinary unweighted pressure refits of all 51
recoverable rows, with P=AD, V/V0=AN, T=U and V0/K0/Tr/n/theta0 fixed:

| theta0 | K0 prime | gamma0 | q | RMSE (GPa) |
|---|---:|---:|---:|---:|
| 768 K (final Table 2) | 4.254577 | 1.000503 | 2.470346 | 0.212395 |
| 760 K (deposited scripts) | 4.249175 | 0.999277 | 2.467664 | 0.212827 |

Both ordinary optima are within the printed bootstrap errors; neither is
the published bootstrap mean. No temperature errors, bootstrap covariance,
or original fitting export is invented. No refit record is installed.

The recommended author-script vector is graded **`parity`** for fit reproduction,
with `parity_basis: within_reported_parameter_errors`. This is conditional
numerical parameter parity: the all-51 reconstructed-input fit recovers all
three free script coefficients within the printed Table 2 error widths.
Those widths accompany the rounded theta0=768 K table parameterization and
remain attached to that record; they are used only as a comparison scale for
the theta0=760 K script vector. Script-specific errors remain missing. Exact
bootstrap reproduction and a formal combined-two-sigma result are unverified.

The explicit rerun now reads the reconstructed six-column CSV itself and
compares an independent optimizer with a translation of the author's own
`NlLSq`, including its temporary olivine guesses, finite-difference
derivatives, damping updates and stopping criterion. Both optimize the same
unweighted squared pressure differences and agree to better than 1e-5 in
every fitted coefficient. See
[the reconstructed-input fit comparison](forsterite-katsura-2022-refit-comparison.md)
for all-51 versus table-matched-50 selections, exact versus polynomial
Debye evaluation, and coefficient differences from the author script.

The table/code theta0 discrepancy is explicit: **the recommended thermal
record follows the scripts at 760 K; the printed-table alternative retains
768 K**. At the same rounded published
free coefficients, the 760 K curve's RMSE is 0.212945 GPa. The alternative
is a diagnostic. The recommended script record uses its full coefficient
vector, rather than changing theta0 alone in the rounded table vector.

## Validity, limitations and reproduction

Saved corrected-data marginal envelopes are approximately -0.047692 to
13.692765 GPa and 304.031 to 1958.093 K. The slightly negative pressure is
near-ambient measurement/reduction scatter, not evidence for tensile phase
stability. The envelopes are not a rectangular guarantee of olivine stability.
The 300 K reference identity is tested with validity checks disabled.

Equation validation, data recovery and ordinary-refit diagnostics are
complete for both 2022 source variants. Original 2009 equation/fit reproduction and
exact 2022 bootstrap parity remain unresolved. No geotherm reproduction is
claimed. In particular, Equation (12) uses a separate empirical ambient
olivine Cp; a generic MGD caloric calculation does not reproduce that whole
geotherm workflow.

```console
.venv/bin/python -m scripts.forsterite_recover_katsura_2022
.venv/bin/python -m scripts.forsterite_reproduce_katsura_2022 --write-report
.venv/bin/python -m pytest -q -W error tests/test_forsterite_katsura_2022.py
cargo test -p peritheos --test forsterite_katsura_2022 --all-features --locked
```

The machine-readable report is
`docs/data/forsterite-katsura-2022-reproduction.json`. Shared ledgers and
catalog totals are intentionally deferred to the material-batch reconciliation;
see `forsterite-katsura-2022-integration.md` for required updates.

## Local check results (shared checkout)

- Forsterite Python checks after the PDF recovery and reconstructed-input rerun: **8 passed**, warnings treated as errors.
- Forsterite Rust integration: **2 passed**; the scoped Clippy check and
  Rust formatting check passed.
- Latest combined forsterite/dataset/pressure-coordinate/supporting-data run
  after the supplied-PDF recovery and reconstructed-input rerun: **81 passed**. All script coefficients
  are checked directly against the archived author parameter arrays.
- Earlier broader `.eosmat`, material and catalog run, before adding the
  author-script alternative: **523 passed, 8 failed**.
  Seven failures are shared total/audit expectations after parallel additions
  (230 versus 227 materials, 657 versus 650 source records, corresponding
  validated/executable counts, and unreconciled audit membership). The eighth
  is the pre-existing Fei Au temperature-coverage assertion. No shared test or
  ledger was changed to remove these failures.
- Dataset, pressure-coordinate, supporting-data, MGD and family suites:
  **138 passed**.
- Repository-wide Ruff lint/format and mypy: **passed**.
- Strict documentation build: **failed on five existing Au-document links**
  in `hirose-2008-gold.md` and `sakai-2025-gold.md`. The forsterite audit links
  introduced no warning; global navigation changes remain deferred.
- Global `validate_primary_eos_refits.py --check`: **failed before comparison**
  with `KeyError: primary_source_check` at line 5237. The current catalog entry
  missing that metadata is the pre-existing
  `gold_hirose_2008_fei_2004_speziale_2001_bm3_refit`. The forsterite record
  includes this field. No Au metadata or shared validator was changed.

This is not a clean complete-project check or a publication-ready batch
integration. The common source/refit ledgers and count assertions still need
the coordinated reconciliation listed in the integration note.
