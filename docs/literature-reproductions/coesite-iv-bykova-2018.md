# Coesite-IV: conditional refit and independent figure check

Record: `coesite_iv_bykova_2018_am05_static_bm3_refit`. Audit: 23 October 2026.

The stored coefficients are unchanged: **V0 = 438.6702 Å³ per Z = 16 cell,
K0 = 168.52 GPa, K0′ = 3.34**, at 0 K. The record is a **Peritheos refit**
of calculated source states, with both published moduli held fixed. It is not
a fully published parameterization or a reproduction of the authors' regression.

## Source and preserved inputs

Bykova et al. (2018), *Nature Communications* 9, 4789,
[doi:10.1038/s41467-018-07265-z](https://doi.org/10.1038/s41467-018-07265-z).
The [official supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-018-07265-z/MediaObjects/41467_2018_7265_MOESM1_ESM.pdf)
has SHA-256 `9d4bcb36ae6a89b19f5bbcc3fc671f3c1bfc6530a76915382c6132d3d907e459`.

Supplementary Table 11 (PDF page 21) prints K0 and K0′, but **no V0 or
parameter uncertainties**. Supplementary Table 10 on the same page supplies
three calculated states that remain coesite-IV:

| Printed calculated pressure (GPa) | Printed cell volume (Å³) |
| ---: | ---: |
| 39 | 368.5011 |
| 41 | 366.1221 |
| 48 | 357.8274 |

The original CSV, rows, normalization and checksum are preserved. The pressures
are printed as whole GPa; these are not exact P–V states or pressure standard
deviations. Dividing cell volume by Z = 16 gives volume per SiO2 formula unit.
The smaller-volume relaxation at approximately 352.35 Å³ transforms to coesite-V
and remains excluded from this coesite-IV fit.

## What is reproduced

Peritheos fits V0 by minimizing the unweighted sum of squared pressure
residuals, with K0 and K0′ fixed to their printed source values. No digitized
points enter this fit. An independent standard BM3 expression gives
V0 = 438.6702185 Å³ and pressure residuals **+0.108394, −0.013112,
−0.087580 GPa**, where residual means P(model) − P(table). RMS is
**0.0808114 GPa**. The stored rounded coefficients give the same RMS to the
displayed precision.

The numerical assessment therefore checks **reproduction of the stored
Peritheos refit**, with an absolute V0 tolerance of 0.00005 Å³ (half the last
stored decimal). A `parity` result has that scope only. The legacy ledger
field named `published` contains the stored Peritheos target for this refit;
it does not denote a V0 printed by Bykova et al. The previous `similar`
assessment incorrectly implied comparison with a source V0.

The formal V0 standard error, about **0.0843 Å³**, comes from the local
one-parameter Jacobian with covariance scaled by residual variance (two
degrees of freedom). It is a **conditional diagnostic**, not a published
uncertainty or a total uncertainty estimate. It excludes pressure rounding,
uncertainty in the fixed moduli, model inadequacy and DFT systematics. The
unweighted residual sum divided by two is 0.00979573 GPa²; absent a source
error model, it is not a statistical reduced chi-square.

If each printed integer pressure is *assumed* to represent nearest-integer
rounding within ±0.5 GPa, holding printed volumes and moduli fixed, the three
constraints intersect at **V0 = 438.09–439.27 Å³**. This is a
rounding-compatible interval, **not a confidence interval**. Its assumption
does not supply missing source errors.

## Quantitative check against Supplementary Figure 6

The figure (PDF page 7) contains eight yellow-diamond calculated coesite-IV
markers. Their centers were manually read from the two original embedded JPEG
tiles, concatenated in page order without rescaling to a 1430 × 915 pixel
raster. Experimental crosses, coesite-V triangles, the legend and connecting
lines are excluded. Independent affine axis maps use all nine listed pressure
ticks and seven volume ticks. The
[pixel readings and calibration](../data/coesite-iv-bykova-2018-figure6-digitization.json)
preserve the full reading procedure.

A conservative **±4 pixels per coordinate** is allowed: two pixels for center
reading and two for calibration, JPEG antialiasing and tile placement. This
corresponds to approximately **±0.157 GPa** and **±0.0203 Å³ per formula unit**.
These are analyst reading bounds, not published error bars, standard deviations
or confidence intervals. This is an analyst choice, not a measured estimate of
reading repeatability; the coordinate maps themselves come only from axis ticks.
Pressure residual intervals are calculated from the exact BM3 extrema over each
rectangular coordinate box, without statistical quadrature or adding a fitted
parameter error. Zero inside an interval means the curve crosses that box.

Markers 5, 7 and 8 appear to correspond to the three Table 10 states. They are
shown for completeness but excluded from the **five additional-point checks**.
No exact identity between their plotted and printed coordinates is assumed.
The additional markers are independent of the three inputs to the Peritheos
fit; they may have entered the authors' own regression.

| Marker | P(figure), GPa | V/Z, Å³ | Pressure residual, GPa | Residual reading interval, GPa | Additional check |
| ---: | ---: | ---: | ---: | ---: | :---: |
| 1 | 27.01 | 24.086 | +0.013 | [-0.357, +0.383] | Yes |
| 2 | 30.11 | 23.786 | +0.142 | [-0.239, +0.523] | Yes |
| 3 | 34.03 | 23.425 | +0.320 | [-0.074, +0.715] | Yes |
| 4 | 36.11 | 23.247 | +0.354 | [-0.047, +0.757] | Yes |
| 5 | 38.97 | 23.044 | -0.016 | [-0.426, +0.394] | No |
| 6 | 40.19 | 22.952 | -0.083 | [-0.497, +0.331] | Yes |
| 7 | 40.97 | 22.891 | -0.092 | [-0.508, +0.326] | No |
| 8 | 47.95 | 22.352 | +0.134 | [-0.306, +0.576] | No |

All eight markers are compatible with the unchanged curve within these
reading bounds: RMS **0.188 GPa**, maximum absolute residual **0.354 GPa**.
For the five additional markers alone, RMS is **0.226 GPa**, with the same
maximum. This is graphical agreement at the stated resolution. In particular,
the approximately 34 and 36 GPa central readings differ by 0.32 and 0.35 GPa;
the result must not be described as exact marker agreement or validation at
arbitrarily smaller tolerances.
Using a tighter total bound of **±2 pixels** excludes markers **3 and 4**.
This sensitivity is retained in the numerical results: compatibility at the
stated four-pixel resolution does not establish agreement at two-pixel precision.

![Unchanged coesite-IV curve and digitized Figure 6 markers with graphical bounds](assets/coesite-iv-bykova-2018-figure6-check.svg#only-light)
![Unchanged coesite-IV curve and digitized Figure 6 markers with graphical bounds](assets/coesite-iv-bykova-2018-figure6-check-dark.svg#only-dark)

## Interpretation and remaining limits

The 39–48 GPa interval belongs to the numerical fit inputs. Figure readings
support the approximate 27–48 GPa plotted branch at graphical resolution.
Neither this marginal range nor the 0 K curve establishes experimental phase
stability or a finite-temperature EOS. The ambient V0 is an extrapolated model
coefficient, not an independently observed ambient coesite-IV volume.

There is no demonstrated numerical failure here that warrants classifying the
record as `not_reproduced`. Source traceability, repeatability of this
conditional refit and graphical agreement are separate findings. None
establishes recovery of the source author's fitted V0, original E–P–V grid,
objective, weights or covariance, or physical accuracy. Those original outputs
would be needed for author-workflow reproduction. No new usage recommendation
or change to automatic selection follows from this audit.

## Reproduction

These are repository-only maintenance commands, not installed-package tools:

```bash
uv run python -m scripts.reproduce_bykova_2018_coesite_iv --check
uv run --group docs python -m scripts.reproduce_bykova_2018_coesite_iv --check --plot
```

The first command replays the Table 10 fit, the rounding diagnostic, and all
eight graphical coordinate checks; it also verifies that the separate
digitized CSV agrees with the archived pixel readings. `--plot` regenerates the
figure without modifying EOS coefficients or the original Table 10 dataset.
The plot's y-axis is volume and its x-axis is pressure.

[Numerical results](../data/coesite-iv-bykova-2018-audit.json) ·
[Reconstruction report](../reconstruction-reports/coesite_iv.md) ·
[Material and dataset provenance](../material-summaries/coesite_iv.md)
