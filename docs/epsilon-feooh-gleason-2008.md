# Gleason et al. (2008) epsilon-FeOOH BM2 audit

[Gleason, Jeanloz, and Kunz (2008)](https://doi.org/10.2138/am.2008.2942)
report `V0 = 66.3(5) A^3`, `K0 = 158(5) GPa`, and a fixed
`K0' = 4` for non-quenchable high-pressure epsilon-FeOOH. This audit uses the
[official OSTI manuscript](https://www.osti.gov/servlets/purl/950968) and
[MSA depository item AM-08-056](https://msaweb.org/MSA/AmMin/TOC/2008/ND2008_data/AM-08-056.zip)
rather than a secondary EOS catalog.

The important correction is conceptual as well as numerical: all 49 deposited
epsilon-FeOOH observations were measured between 200 and 400 degC. None is a
room-temperature row. The source nevertheless pooled all of them in its
"simplest analysis" after separate high-temperature and isotherm-by-isotherm
checks found no resolvable softening or volumetric-expansion trend. Peritheos
now reproduces that source-specific reduction instead of discarding 48 rows or
silently pretending that the hot measurements lie on a 300 K isotherm.

## Phase and data selection

The source refined alpha-FeOOH as `Pbnm` and the high-pressure phase as
`P21mn`. The latter is an axis setting of space-group 31; Peritheos stores the
same hydrogen-off-center InOOH-type phase in the axis-permuted `P21nm` setting.
It is distinct from both goethite and the hydrogen-centered, low-spin `Pnnm`
epsilon-FeOOH material.

The EOS selection contains every one of the 49 Table 2 rows:

- pressure: 3.81-21.01 GPa;
- measured temperature: 473.15-673.15 K (200-400 degC);
- conventional-cell volume: 60.40-65.03 A^3;
- phase: only rows refined as high-pressure epsilon-FeOOH;
- exclusions: none.

The deposit does not label compression and decompression row by row, although
Figure 3b distinguishes them graphically. That missing direction label does not
create a hidden selection because the paper explicitly says all epsilon-phase
data were included.

The checked-in CSV is a verbatim, source-order transcription. Its SHA-256 is
`d072bc326b925971e6a79e028016536fa09ff8f3b598864ce194cada0c60d7db`.
The downloaded MSA archive and contained Word document have SHA-256 values
`11cb79a0961ae7664650265d72aadc6f1997caf8c6b2984950c6738aca3aef07`
and `2d4dabfb418a8244316d5b159f3fdd5f401c012ecc350ab25fe66528e14d1e69`,
respectively.

## Pressure and temperature reduction

The experiment mixed sample and gold 7:1 by weight in non-dried 4:1
methanol-ethanol. Ruby, NaCl, and gold standards were used for the separate
ambient-temperature compression experiment. The manuscript states that only
gold, following Shim, Duffy, and Takemura (2002), was used for the
high-temperature runs. Consequently every epsilon-FeOOH row is on that
temperature-aware gold pressure scale.

The depository table contains the already reduced pressures but not the gold
lattice parameters. The pressure calibration is therefore identified but
cannot be recalculated observation by observation.

The epsilon-specific prose does not print a sample-volume temperature
correction. A numerical reconstruction, however, identifies the missing step.
Using the same paper's zero-pressure volumetric expansivity
`alpha0 = 2.3(6)e-5 K^-1` to reduce each measured volume to 300 K,

\[
V_{300}=V(P,T)\exp[-\alpha_0(T-300\ \mathrm{K})],
\]

and then holding the separately extrapolated `V0 = 66.3 A^3` gives
`K0 = 158.099 GPa`, which rounds to the published value. The constant-alpha
exponential is the thermodynamic integral of `alpha = (1/V)(dV/dT)`.

This is strong numerical evidence for the source reduction, but not a printed
epsilon thermal-expansion measurement. Peritheos therefore records it as a
fit-time provenance operation. It does not add `alpha0` to the executable
epsilon EOS or claim a general `P(V,T)` model.

## Staged fit and weighting

The source says that the ambient-pressure room-temperature volume of the
non-quenchable phase was determined with the Jeanloz (1981) `G`-versus-`g`
method. An unweighted linear normalized-stress intercept check on the pooled
raw rows gives `V0 = 66.1628 A^3`, inside the published `66.3(5) A^3`
interval. The rounded published `V0` is then held fixed in the reproduced BM2
stage, while `K0' = 4` is intrinsic to BM2 and `K0` is the only free
coefficient.

The manuscript and deposit do not disclose EOS weights, the fitted residual
variable, covariance, or a scalar fit statistic. Parenthesized pressure,
temperature, lattice, and volume limits in the deposit are explicitly 95%
Rietveld confidence limits; they are not described as EOS weights and are not
silently converted to one-standard-deviation errors. The reproducible fit uses
ordinary unweighted pressure residuals.

## Numerical reproduction

| Calculation | `V0` (A^3) | `K0` (GPa) | Pressure RMSE (GPa) |
|---|---:|---:|---:|
| Publication | 66.3 +/- 0.5 | 158 +/- 5 | - |
| Temperature-reduced published curve | 66.3 fixed | 158 fixed | 0.945732 |
| Temperature-reduced unweighted refit | 66.3 fixed | 158.098698 +/- 1.672307 | 0.945698 |
| Raw hot volumes, unweighted refit | 66.3 fixed | 175.624263 | 0.752364 |
| Raw hot volumes, both coefficients free | 66.139661 | 182.870910 | 0.733335 |

The raw-volume controls matter: they show that direct evaluation of every hot
row against the reference BM2 does not recover the paper. They are retained as
falsification diagnostics, not alternate production parameters.

Varying the adopted `alpha0` across the paper's quoted interval changes the
fixed-`V0` refit from 162.415 GPa at `1.7e-5 K^-1` to 153.951 GPa at
`2.9e-5 K^-1`, consistent with the scale of the published `+/-5 GPa`
uncertainty. Exact coefficient covariance remains unrecoverable.

Run the standalone reproduction with:

```bash
uv run python scripts/reproduce_gleason_2008_epsilon_feooh.py
```

## Use recommendation

Use `e_feooh_gleason_2008_bm2_1` as the paper's 300 K reference BM2 and as a
provenance-preserving early epsilon-FeOOH compression result. Do not use the
stored isothermal object to predict the deposited hot states directly. For an
explicit low-pressure thermal EOS, prefer the separately audited Suzuki (2016)
record, whose five thermal and elastic coefficients were jointly fitted and
printed by its source.
