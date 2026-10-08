# Wang (1996) and Shim (2002) experimental data recovery

Recovery checkpoint: 2026-10-08. This adds published numerical observations and
explicit figure digitizations; it does not change any EOS coefficient, reported
parameter error, phase identity or pressure scale.

## Wang et al. (1996)

The [laboratory-hosted paper](https://millenia.cars.aps.anl.gov/gsecars/LVP/publication/Papers/Wang_et_al95JB03254.pdf)
contains the complete numerical Table 1 on journal page **664**, PDF page 4.
The existing catalog held only its 14 room-temperature Run-13 states. The full
table is now bundled as `ca_perovskite_wang_1996_table1_pvt`:

| Run | Rows | Origin stated in the paper |
| --- | ---: | --- |
| 13 | 34 | Wang et al. (1996) |
| 3 | 15 | Wang and Weidner (1994), Table 1 footnote |
| 4 | 12 | Wang and Weidner (1994), Table 1 footnote |
| 5 | 5 | Wang and Weidner (1994), including two continuation rows |
| Total | 66 | 14 room-temperature and 52 elevated-temperature states |

The transcription preserves printed pressures, temperatures, parenthetical
volume errors, signed differential stresses, run ordering and footnote flags.
The source spans 0.59–12.25 GPa and 301–1594 K. Differential stress is not a
pressure measurement error. No row-specific pressure or temperature uncertainty
is invented, and the full table does not assign a confidence level to the
parenthetical volume errors. Raw parenthetical strings remain available.
The old subset's `volume_sigma_a3` name is retained for compatibility, with its
metadata role qualified as an unspecified-confidence uncertainty rather than
an established Gaussian standard deviation.

Two excluded rows in the old room-temperature subset contained transcription
errors: **1.137 → 1.13 GPa** and **0.595 → 0.59 GPa**. The trailing dagger and
double dagger had been mistaken for digits. Both remain excluded; the 12
fit-included rows are unchanged. The old table location was also corrected from
page 665 to page 664. These corrections and their previous values are recorded
in provenance.

The room-temperature dataset is a subset of the full table, not an independent
experiment to concatenate with it. The inherited 1994 runs are likewise not
new independent 1996 observations. The full table is a standalone source
dataset; no new fit manifest or EOS is assigned by this recovery.

This removes the missing-numerical-table obstacle for examining the preferred
Shim (2000, JGR) combined DAC/LVP fit. Its exact LVP row selection and regression
protocol still need an audit before that combined fit can be claimed reproduced.
The existing DAC-only record remains unchanged numerically.

The subsequent [Wang fit audit](wang-1996-casio3.md) reproduces the rounded
equal-weight Mao reanalysis and checks all direct Table 2 fitting approaches.
Thermal reproduction remains qualified because weights, masks and some
unweighted results cannot be recovered exactly from the printed information.

## Shim et al. (2002)

The [author-hosted PDF](https://duffy.princeton.edu/sites/g/files/toruqf616/files/shim_et_al-2002-grl.pdf)
stores Figure 2 as vector paths. Six unique P–V marker centers and their plotted
error-bar halfwidths are now bundled as
`casio3_perovskite_tetragonal_shim_2002_figure2a_digitized`.

There are **five tetragonal refinements and one cubic refinement at 19.7 GPa**.
The latter is an open symbol, as explained in the caption. A black path beneath
that white symbol is an overdrawn duplicate, not a seventh observation. Legend
symbols, comparison curves and the c/a panel are excluded.

Paragraph 7 prints the six pressures (19.7, 24.2, 25.2, 26.6, 36.1 and 45.8 GPa).
These published values are used as the pressure observations. Independently
digitized horizontal coordinates remain in `plot_pressure_gpa`; they differ
from the printed pressures by up to 0.062 GPa. The coordinate-reading bound
does not account for source plotting offsets or source rounding.

The figure plots V/V0. The conversion `V = 45.58 * (V/V0)` uses the paper's
adopted V0 in Å³ per one-formula-unit pseudocell. The plotted ratio, vector
coordinates, conversion and per-row refinement type remain explicit. A separate
0.1 PDF-point coordinate bound corresponds to approximately 0.024 GPa and
0.0069 Å³. The pressure bound applies to `plot_pressure_gpa`, not the printed
pressures. This is a conservative digitization bound, not a Gaussian sigma.
The caption defines source error bars as 1-sigma, but digitized halfwidths are
still approximate; one volume bar is smaller than the coordinate bound.
Additional stored decimals make extraction reproducible, not experimentally
more precise. These are room-temperature measurements after annealing/quench.

`scripts/digitize_shim_2002_figure2.py` reproduces the CSV from the exact source
PDF, checked by SHA-256, using optional PyMuPDF:

```sh
uv run --with pymupdf python scripts/digitize_shim_2002_figure2.py SOURCE.pdf OUTPUT.csv
```

As a diagnostic, fixing V0=45.58 Å³ and K0′=4 gives:

| Pressure-residual weighting | Fitted K0 (GPa) |
| --- | ---: |
| Unit weights | 262.2289 |
| Digitized pressure error bars | 254.7407 |
| Effective errors in both axes: sqrt(σP² + (dP/dV σV)²) | 255.0081 |

The weighted checks agree with the published **255(5) GPa**, but the original
regression objective, numerical weights and parameter covariance remain
unavailable. The checks use digitized source bars only, not combined Gaussian
digitization errors. No reproduction of the original optimizer or its parameter
errors is asserted, and no published coefficient is replaced. The historical
catalog-wide refit ledger has not been regenerated for this new source recovery.

## Accepted experimental recovery scope

The user accepted the existing qualified Fu (2023) parameterization unchanged,
the **46 printed Mao (1989) rows** despite the prose count mismatch, and the
Wang–Weidner (1994) digitized room-temperature points. No 47th Mao observation
is fabricated; the suggested prose typo is an interpretation, not independently
established evidence. Further original MD/DFT-data recovery is outside this
experimental-data scope.

The recovered datasets supply the published Wang numerical table and the Shim
plotted observations. Unrounded original measurements, original fit protocols
and calibrant observations remain distinct questions. Any CC0 statement for
Peritheos contributors covers transcription, normalization and metadata only;
it does not change third-party source rights.
