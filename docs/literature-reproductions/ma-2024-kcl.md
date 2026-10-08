# Ma, Sumita and Murakami (2024): B2-KCl primary pressure scale

The nondefault `kcl_b2_ma_2024_bm3_mgd` record implements the final published
Table 1 BM3 reference isotherm plus Mie-Grüneisen-Debye thermal increment.
The existing Dewaele default and all earlier KCl records are retained.
The executable parameterization reproduces the deposited model grids within
the discrepancy caused by rounded published coefficients. The joint acoustic
and low-pressure diagnostic reproduces the BM3 coefficients within the
published uncertainty widths, but its weighting and covariance must not be
attributed to the authors.

## Final sources and custody

- Ma, N., Sumita, T., and Murakami, M. (2024), *Primary Pressure Scale of KCl
  B2 Phase to the Core-Mantle Boundary*, JGR Solid Earth 129, e2024JB028819,
  [final publication](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JB028819),
  published 29 June 2024. Authority: Table 1, Section 3 and Equations 4-20.
- [Final cited author deposit](https://data.mendeley.com/datasets/7svmv9hvft/1),
  DOI `10.17632/7svmv9hvft.1`, version 1, published 9 April 2024, CC BY 4.0.
  The original `Data Tables.xlsx` is bundled as
  `peritheos/data/datasets/ma_2024_sources/data-tables.xlsx`. Its SHA-256 is
  `e4484174250276f092a8abd0e3b7a4f76d9a7a2cb853b76cf5ca7d54da057b84`,
  matching the public repository file metadata. The original archive preserves
  S1-S4, source text, comparison coefficients, and numerical precision.
- The older deposit `10.17632/mws6hnp49j.1` is **not** used as scientific input.
- The publisher lists `2024JB028819-sup-0001-Supporting Information SI-S01.docx`.
  Its download returned HTTP 403 during this audit. All numerical S1-S4 tables
  were recovered from the final article's cited author deposit; the additional
  supporting text and figures in the publisher DOCX were not independently
  inspected. The DOI and filename remain recorded in the source manifest.

No original source is relabeled CC0. The workbook and its attributed numerical
extractions retain the deposit's CC BY 4.0 terms. Repository code retains the
project license. Acquisition date: 8 October 2026.

## Equation and units

Table 1 gives V0=32.48(9) cm³/mol, KT0=21.33(70) GPa, K0'=4.836(83),
gamma0=1.92(11), and thetaD0=251(22) K. Parentheses are uncertainties in
the last digits; their confidence level is unspecified. The source fixes
q=1 and uses two atoms per KCl formula unit. The source's shear coefficients
are preserved in S2 but are not an additional executable shear model here.

The B2 conventional cubic cell has one formula unit. Therefore

```text
V0[cell Å³] = 32.48 × 10²⁴ / N_A = 53.93430890180654
error[V0 cell Å³] = 0.09 × 10²⁴ / N_A = 0.1494485160456462
N_A = 6.02214076 × 10²³ mol⁻¹ (exact)
```

The public material API uses cell Å³; its thermal implementation converts
internally to J/bar/mol (molar cm³/mol divided by ten). This conversion is
required for the thermal energy in J/mol to yield GPa. The structural lattice
continues to use the existing default record's reference volume.

Equations 13-14 map exactly to the existing `integrated_gruneisen` choice:

```text
gamma(V) = gamma0 × (V/V0)^q
theta(V) = theta0 × exp[(gamma0/q) × (1 - (V/V0)^q)]
P(V,T) = P_BM3(V,300 K) + gamma(V)/V × [E(V,T) - E(V,300 K)]
```

The zero-point term in Equation 20 cancels at fixed volume in this difference.
`Tr=300 K`, `q=1`, and `n=2` are fixed. No constant-exponent or
variable-exponent replacement is made. At 15 cm³/mol the published law gives
about 706 K, consistent with the article's stated Debye temperature.

The EOS is an experimentally constrained thermodynamic synthesis: measured
velocities and volumes determine elastic and Debye inputs, and the inferred
thermal parameters extend the measured room-temperature isotherm within a
quasi-harmonic model. Its determination method is `experimental`; this does
not imply experimental high-temperature coverage.

## Observations, reductions, and model outputs

| Bundled dataset | Rows | Role |
|---|---:|---|
| `kcl_ma_2024_acoustic` | 11 | S1 XRD molar volume, longitudinal/transverse velocities and Raman shift; reduced KT and theta; separately named modeled pressures |
| `kcl_ma_2024_walker` | 8 | S1 original Walker temperature/pressure and author-recalibrated 300 K P-V inputs; separately named modeled pressures |
| `kcl_ma_2024_cold_grid` | 186 | S3 modeled 300 K P-V grid and propagated errors, ending at 144.1826 GPa |
| `kcl_ma_2024_thermal_grid` | 371 | S4 modeled increments at 32.48 cm³/mol, 300-4000 K |

The existing `kcl_walker_2002_table2_pvt` contains the upstream NaCl lattice
readings, their ESDs, KCl cell volumes and their ESDs. The dedicated audit
maps all eight S1 low-pressure rows to those original spectrum identifiers.
The more precise molar volumes in Ma S1 are kept as deposited, rather than
replaced with the original Walker table's rounded molar-volume column.

Only KT-V and author-recalibrated P-V coordinates enter the joint fit. The
diamond Raman pressures monitor and compare with the derived scale; they are
not used to calibrate the source EOS. The Ma primary-scale pressures in S1
and all S3-S4 predictions are excluded from fitting. No pseudo-observations
are generated from the published EOS.

The S1 footnote specifies two standard deviations for diamond pressure and
VL/VT errors. Other row errors have unspecified confidence. The typed columns
preserve them as `uncertainty`, not automatically as one-sigma errors.
The paper propagates fitted covariance into P(V) and KT(V) envelopes, but
does not deposit the numerical covariance matrix. Marginal parameter errors
do not define that covariance. The source pressure-grid errors are preserved
as supplied and do not become independent pressure measurements.

## Joint and staged diagnostic reproduction

Run from the repository root:

```bash
uv run python scripts/reproduce_ma_2024_kcl.py
# Optional lossless re-extraction from the checksummed bundled workbook:
uv run python scripts/reproduce_ma_2024_kcl.py \
  --extract-source peritheos/data/datasets/ma_2024_sources/data-tables.xlsx
```

The standard-library XLSX parser reads cached source values, verifies the
final file checksum, preserves source coordinates, and needs no spreadsheet
software. The result is in `docs/data/ma-2024-kcl-reproduction.json`.

The source specifies simultaneous BM3 fitting of all 11 KT-V and eight
low-pressure P-V points, after acoustic reduction. It does not specify the
complete residual definition, weights, or covariance scaling. Three explicit
objectives therefore quantify the weighting sensitivity:

| Objective | V0 (cm³/mol) | K0 (GPa) | K0' |
|---|---:|---:|---:|
| Published Table 1 | 32.48 ± 0.09 | 21.33 ± 0.70 | 4.836 ± 0.083 |
| Residuals divided by deposited KT/P error widths | 32.479374 | 21.348728 | 4.846212 |
| KT errors also include volume-error propagation; P errors as deposited | 32.472733 | 21.415178 | 4.830811 |
| Unweighted KT/P residuals in GPa | 33.349032 | 16.371147 | 5.303925 |

The first diagnostic has KT RMS residual 7.0510 GPa and low-pressure P RMS
residual 0.02893 GPa. Each coefficient is within its published error width.
The catalog refit ledger nevertheless classifies this as `similar`, because
source weighting and uncertainty confidence are unavailable. A diagnostic
residual-scaled inverse-JᵀJ covariance is saved in the report, with molar V0
basis; it is not attached to the published record as source covariance.
The volume-propagation sensitivity includes the measured acoustic volume
errors but does not invent low-pressure KCl-volume covariance or independent
one-sigma interpretations.

The separate Equation 13 Debye-temperature diagnostic uses all 11 reduced
theta(V) values with q fixed at one. Weighting by the deposited theta-error
widths gives A≈1710.79 K, B≈0.0590200 mol/cm³, gamma0≈1.91697 and
theta0≈251.574 K at the published V0, consistent with Table 1 and the
article's A=1716(117), B=0.0592(35). The unweighted sensitivity is also saved.

### Acoustic Equation 9 discrepancy

As printed in the final HTML article, Equation 9 contains `hbar/(2 kB)` and
Avogadro number without a two-atom number-density factor. Evaluating that
expression gives theta values smaller than S1 by `2 × 2^(1/3)`.
The conventional two-atom expression
`(hbar/kB) × (6 pi² × 2 N_A/V)^(1/3) × v_D` reproduces the S1 theta values
to less than 4 microkelvin. Equations 4 and 7-12 with that convention,
the published rounded gamma(V), and an assumed molar mass 74.5513 g/mol
recover the deposited KT inputs to within 0.021 GPa, far below their
reported error widths. This is a numerical identification of the deposited
convention, not an author-confirmed erratum. The executable thermal EOS uses
the published theta0 directly, so no Equation 9 typo correction changes its
coefficients or the existing reusable thermal model.

### Walker recalibration boundary

Ma's eight deposited reduced pressures span 3.15786-8.07891 GPa at 300 K.
Their ancestry is Walker's paired KCl/NaCl diffraction observations,
originally on the Birch (1986) NaCl B1 scale, recalibrated to the acoustic
NaCl B1 primary scale of Matsui et al. (2012), DOI
[10.2138/am.2012.4136](https://doi.org/10.2138/am.2012.4136).
The [Matsui primary article](https://rruff.info/doclib/am/vol97/AM97_1670.pdf),
Table 2 and Equations 2, 4-11, specifies V0=179.425 Å³ (Z=4),
K0=23.7 GPa, K0'=5.14, K0''=-0.392 GPa⁻¹, gamma0=1.56,
theta0=279 K and q=0.96, with a BM4 isotherm and integrated Debye law.
The source PDF checksum and inspected page locations are in the manifest.

Direct application to the original rounded NaCl lattice readings, evaluating
NaCl pressure at 296.15/297.15 K and adding Walker's KCl
alphaKT=0.00275 GPa/K increment to 300 K, misses Ma's deposited pressures
by 0.0443-0.1471 GPa. Seven Run 1 rows miss by about 0.044-0.064 GPa;
the one Run 2 row misses by about 0.147 GPa. These differences are not
explained by the tiny room-temperature correction alone. The precise
upstream calculation and unrounded calibration inputs are unavailable.
No cause, correction, or exact original-table reproduction is asserted.

Both the author-deposited input and this independently computed diagnostic
remain visible. The joint reproduction uses the deposited input, so it
establishes numerical recovery from Ma's final reduced dataset, not full
reproduction of every upstream calibration step. Matsui's NaCl EOS is used
only in the audit and is not added to the catalog in this change.

## Model-grid checks and coverage

The rounded published coefficients differ from S3 by at most 0.10681 GPa,
below every supplied S3 error width, and from S4 by at most 0.013944 GPa.
Independent reconstruction of the **model grids only** recovers:

```text
S3: V0=32.48441270403208, K0=21.33455152957105, K0'=4.836059266516246
S4 at 32.48 cm³/mol: gamma0=1.9224554719645683, theta0=250.8877146371014
```

These values all round to Table 1 and reconstruct the grids to about
10⁻¹³ GPa. They explain the discrepancies from rounded coefficients, but
are not observation refits, new selectable EOS records, or source covariance.
The catalog stores the published rounded Table 1 values.

Experimental coverage is the room-temperature acoustic scale to about
85 GPa, plus the low-pressure Walker anchors brought to 300 K. Pressure
above that acoustic coverage is extrapolation. The extension to 4000 K and
CMB-like conditions is modeled within the quasi-harmonic approximation,
not direct simultaneous high-P-T measurements in this study. B2 KCl is
unstable at the fictive zero-pressure reference state and transforms to B1
near 2-3 GPa; marginal ranges do not assert a rectangular stability field.

The tests check the original workbook hash and lossless extraction, all row
counts, unit conversion, source-grid benchmarks, the joint and acoustic
reductions, the exposed upstream discrepancy, default preservation,
thermal inversion and arrays, native/Python agreement, and `.eosmat`
round trips. Generic catalog and scientific-ledger checks cover the addition.
