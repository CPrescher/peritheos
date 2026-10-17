# Forsterite: reconstructed-input fit versus the author scripts

2026-10-15. Ordinary least squares on the reconstructed 2022 input gives
coefficients within 1.2% of the full-precision author-script vector. The fit
reproduction grade is **`parity`**, on the explicitly conditional basis
`within_reported_parameter_errors`: the three free coefficients are recovered
within the printed Table 2 error widths. The recommended thermal EOS retains
the author's full-precision script values and theta0=760 K. The
author's own deterministic optimizer, translated literally, converges to
the same result as SciPy. This validates the ordinary-fit result, not the
unavailable bootstrap uncertainty analysis. Published/source-script catalog
parameters are unchanged.

## Input and fixed choices

The fit reads `forsterite-katsura-2022-candidate-fitting-input.csv` directly.
All 51 deposited observations are retained in the main fit. Pressure uses
the author's saved combined AD coordinate; sample V/V0 is AN and corrected
temperature is U. AE and AO errors survive unchanged. SigmaT remains empty.
The separate original `Olivine_PVT2021.xlsx` is unavailable, so its precise
pressure-column selection cannot be certified.

The fixed inputs are V0=289.58 A³/cell, K0=127.4 GPa, Tr=300 K,
theta0=760 K, n=7 atoms/formula unit, Z=4. The three free coefficients are
K0 prime, gamma0 and q. The objective is the sum of squared pressure
residuals, without uncertainty weighting, as in the author's `ChisqPressure`.
Missing sigmaT does not enter this deterministic objective. It prevents an
unqualified reconstruction of the source bootstrap perturbations.

The primary comparison retains the fitting script's Debye polynomial,
R=8.314472 and NA=6.0221367e23. The comparison script vector is the actual
olivine vector in `Geotherm calculation v2.m`, rather than the fitting
program's temporary starting guesses.

## All 51 deposited observations

| Coefficient | Author script | Refit with author equation | Relative difference |
|---|---:|---:|---:|
| K0 prime | 4.211168970 | 4.249221000 | +0.9036% |
| gamma0 | 0.997258585 | 0.999272576 | +0.2020% |
| q | 2.438540914 | 2.467717215 | +1.1965% |

Pressure RMSE is **0.213001256 GPa** at the script coefficients and
**0.212820649 GPa** after refitting. Across these observations the fitted
curve differs from the script curve by at most **0.0203281 GPa**.

The literal author optimizer returns K0 prime=4.249225248,
gamma0=0.999272624 and q=2.467724317. Its largest coefficient difference
from the independent optimizer is below 0.000008, consistent with the
author's relative-step stopping tolerance of 1e-4. Thus optimizer choice
does not explain the remaining difference from the script vector.

The exact Peritheos Debye model gives K0 prime=4.249175494,
gamma0=0.999277223 and q=2.467664523, with RMSE=0.212827064 GPa. These
small equation-evaluation differences are distinct from fitting workflow
or source-data selection differences.

## What the supplied 2009 PDF changes

Table 1 on journal p.89 / PDF p.4 contains **52 observations**:
11 in run M371 and 41 in M375. It includes two ambient reference rows
that are absent from the deposited corrected workbook. Fifty workbook rows
match the printed table by run, uncorrected temperature and original
Matsui pressure within 0.02 GPa. The tolerance accommodates minor
printed/deposited pressure rounding differences, and every match is unique.

The remaining workbook observation is **M375006**, at uncorrected T=1300 K,
original deposited Matsui P=3.897307692 GPa, corrected T=1312.087488 K
and combined corrected P=3.735453111 GPa. It has no matching Table 1 row.
Its pressure residual at the script coefficients is +0.972333 GPa.

The table does not establish whether M375006 was a rejected measurement,
a workbook-only intermediate observation, or part of the original 2022
fitting export. Its absence is evidence for a separate source-table
sensitivity comparison, not permission to delete it silently.

## Fifty rows matched to the published table

This sensitivity fit keeps the same corrected 2022 P, T and V/V0 values,
but selects only observations matched to 2009 Table 1. It does not refit
the original 2009 pressure scale or its combined Takeuchi dataset.

| Coefficient | Author script | Table-matched 50-row refit | Relative difference |
|---|---:|---:|---:|
| K0 prime | 4.211168970 | 4.096436619 | -2.7245% |
| gamma0 | 0.997258585 | 1.008933742 | +1.1707% |
| q | 2.438540914 | 2.468297270 | +1.2203% |

The script curve's RMSE on this selection is **0.165432921 GPa**; the
refitted curve's RMSE is **0.162259167 GPa**. The author's deterministic
optimizer again agrees with the independent optimizer within 0.000004
in every coefficient. Lower scatter comes with a different dataset; it is
not evidence that this selection reproduces the original author export.

## Interpretation and reproducibility

The all-51 ordinary optimum achieves **conditional numerical parameter parity**
under the [catalog's fit-reproduction categories](../catalog.md#fit-reproduction).
Absolute shifts divided by the printed error widths are 0.0951 for K0 prime,
0.1007 for gamma0 and 0.0584 for q, all below one. These widths (0.4, 0.02,
0.5) belong to the rounded Table 2 parameterization at theta0=768 K; they are
used as a comparison scale and are not assigned as script-specific errors.
This grade does not assert a formal combined-two-sigma comparison or reproduce
bootstrap uncertainties. Script-specific errors and covariance remain missing.

Exact equality is not expected from the available evidence because the
publication reports means of fits to 1,000 resampled/perturbed datasets,
whereas this rerun fits the unperturbed reconstructed observations once.
The temperature-error column, exact original export selection and original
bootstrap realizations remain missing. No parameter uncertainty or
covariance is invented, and no refit is installed as a replacement EOS.

The 2009 full text independently specifies theta0=760 K, supporting the
script's adopted value. It also adopts V0,300=290.79 A³ in its own analysis
and combines its measurements with Takeuchi et al. (1984). Those source
choices must not be mixed with the 2022 script's V0=289.58 A³ and corrected
input to claim original-2009 fit reproduction.

```console
.venv/bin/python -m scripts.forsterite_recover_katsura_2009 /path/to/1-s2.0-S0031920108002227-main.pdf
.venv/bin/python -m scripts.forsterite_reproduce_katsura_2022 --write-report
.venv/bin/python -m pytest -q -W error tests/test_forsterite_katsura_2022.py
```

All numerical results and method/selection details are stored in
`docs/data/forsterite-katsura-2022-reproduction.json`, under
`reconstructed_input_refits`. The 52 original table rows and source tokens
are bundled as `peritheos/data/datasets/forsterite-katsura-2009-table1.csv`.

Verification: **81 scoped Python tests and 2 Rust integration tests passed**,
including the original-table transcription, unique workbook join, actual
candidate-input fit, direct author-optimizer comparison, schema, dataset
checksums and source-record preservation. Scoped Ruff lint/format and Clippy
also passed. Shared count/ledger/documentation issues from the earlier audit
were not changed or reclassified by this rerun.
