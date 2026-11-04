# User-requested independent Pd BM3 check with measured V0 fixed

17 October 2026. These fits describe Baty et al.'s experimental 300 K Pd
observations in official supplementary Table S1. They are separate from the
authors' published coefficients and from the separate 0 K DFT EOS. They have
not been added to the Peritheos catalog or made a default.

## Results and interpretation

| Curve | V0 (Å³/conventional fcc cell) | K0 (GPa) | K0' | All-78-row pressure RMS (GPa) |
|---|---:|---:|---:|---:|
| Published experimental BM3 | 58.88 | 190 | 5.3 | 1.273118 |
| Independent BM3, ambient V0 fixed | 58.89565212 | 206.545539 | 4.32161377 | 0.963529 |

**V0 is fixed to the measured ambient volume in the active analysis.** The
updated plots use pressure on the x-axis and atomic volume on the y-axis,
including volume residuals as a function of pressure. The original
pressure-residual objective is unchanged by this plotting convention.

The fixed-V0 fit is a fit-quality diagnostic. Its systematic residuals do
not support adopting a single BM3 EOS for quantitative use. Its fitted K0
and K0' are not independently established elastic constants. The maximum
absolute pressure residual is 2.038892 GPa; the volume RMS is
0.0426465 Å³/atom. Residuals are systematically positive in volume at low
pressures and change sign across the dataset, rather than showing only
random scatter around one smooth BM3 curve. This supports investigating an
apparent compressibility anomaly, without identifying its physical cause or
establishing a phase transition. Source pressure-reduction systematics also
remain unresolved.

The earlier free-V0 fit is preserved in
[the historical numerical report](fit-results-with-free-v0-historical.json)
and [historical residuals](residuals-with-free-v0-historical.csv), but is no
longer recomputed or displayed as a candidate. It predicts 2.876995 GPa at
the measured ambient volume and shifts V0 upward by 1.81%. Allowing that
reference-volume drift does not resolve the scientific problem.

## Explicit fitting choices

- The anchored fit takes V0 from the table's single P=0 observation and fits
  K0 and K0' to all 77 finite-pressure observations with equal pressure
  weights. The ambient row supplies the anchor; no source point is discarded.
- Source dataset-wide uncertainty limits are not assigned as invented
  row-wise Gaussian sigmas. No outlier removal, new EOS family, or thermal
  term is introduced.
- The numerical bounds are broad and do not constrain the solutions. Initial
  values, bounds, software versions, complete conditional covariance,
  parameter correlations, fit statistics, and the data checksum are stored
  in [fit-results.json](fit-results.json).
- Source observations remain unchanged. The checksum is
  `513835bba914f6fa1f02646901e0453edda808c6da7bc425fbb71c8cf1c99829`.

The anchored fit's conditional standard errors are 2.372 GPa for K0 and
0.1263 for K0', with correlation -0.9704. They assume an exact V0 anchor and
independent, equal-variance pressure residuals, with covariance scaled by
SSE/75. They do not include source calibration systematics, measurement
correlations, or model discrepancy. Moving the anchor by the source atomic
volume uncertainty limit (±0.01 Å³/atom, equivalent to ±0.04 Å³/cell) changes
K0 from 204.186 to 208.934 GPa and K0' from 4.25709 to 4.38620. That is a
sensitivity range, not a statistical confidence interval.

## Verification and boundaries

An independent SciPy fit using a literal BM3 equation agrees with the native
Peritheos fixed-V0 curve within 2.7e-8 GPa at every observation. The curve has
positive bulk modulus and passes pressure-volume inversion checks from zero
to 73.7963 GPa. These numerical checks establish implementation consistency,
not scientific accuracy.

Leave-one-out prediction on the same table gives 0.989 GPa for the anchored
fit (77 finite-pressure predictions, retaining the ambient anchor). Run
identities are absent, so this cannot test independent-run or calibration
systematics.

Inspect the diagnostic only at 300 K within the measured compression
range: 0–73.7963 GPa, observed cell volumes 47.24163336–58.89565212 Å³.
There are no measured finite-pressure rows below 5.96952 GPa. There is no
basis here for thermal extrapolation or a precise derivative at ambient
conditions. The unresolved published-fit finding remains in
[the original audit](../../docs/literature-reproductions/palladium-baty-2024-disposition.md).

## What the paper discusses about a possible transition

The complete publisher PDF was checked from the
[University of Valencia repository](https://hdl.handle.net/10550/95926),
bitstream `db263712-69d9-4eff-879e-60eb31ea4343`, SHA-256
`57b16bea12a74cd6bf84b32c26b88b00bc6b14cfb5a21b6e824b905c2d5bdcab`.

On printed page 4 (PDF page 5), the authors claim a close single-BM3 fit to
their 0–80 GPa experimental observations. They attribute the offset between
the 0 K DFT curve and the 300 K measurements to thermal expansion, decreasing
under compression. They do not interpret the room-temperature dataset as
an electronic or isostructural transition. The experimental section on
printed page 2 reports fcc diffraction to the highest measured pressure.
Its introduction mentions hydrogen contamination as a possible cause of
differences between published Pd EOSs, without establishing it for this run.
The electronic-structure discussion on printed page 9 concerns similarities
between the much broader theoretical phase diagrams of Pd and Ag.

A compressibility anomaly is a hypothesis worth investigating. The
unchanged fcc indexing would not by itself exclude an isostructural
electronic transition, but the present residuals do not distinguish that
mechanism from experimental or reduction systematics. A mechanism or
transition pressure has not been established by this analysis.

## Files and reuse

- [Comparison figure](fit-comparison.png); [vector version](fit-comparison.svg).
- [Parameters and fit provenance](fit-results.json).
- [Residuals at every observation](residuals.csv), with row numbers referring
  to data rows in the bundled, pressure-sorted CSV, excluding its header.
- [Volume residuals versus pressure](volume-residuals.csv), using
  V_observed - V_model(P_observed) in Å³/atom.
- [Reproduction script](refit_and_plot.py). It reads the bundled table and
  regenerates only the files in this output folder.

For reproduction of the compression-curve diagnostic with Peritheos, use the
anchored coefficients explicitly (volumes must be conventional-cell Å³):

```python
import json
from pathlib import Path
from peritheos.eos.rt import BM3

folder = Path("output/palladium-baty-2024-refit-2026-10-17")
results = json.loads((folder / "fit-results.json").read_text())
parameters = results["refits"]["ambient_anchored"]["parameters"]
curve = BM3(**parameters)
pressure_gpa = curve.pressure(52.0)
```

The BM3 class itself does not enforce this study's scientific range or usage
qualifications. These must accompany any reused coefficients.

Reproduce from the repository root with:

```console
.venv/bin/python output/palladium-baty-2024-refit-2026-10-17/refit_and_plot.py
```
