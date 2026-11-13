# Grossular: datasets driving θ with the 2% Conrad diagnostic

**Subsequent source/covariance test:** the [Du/Bosenick study](grossular-du-bosenick.md)
confirms that Du's temperatures are already NaCl-corrected and propagates
Bosenick's shared measured/calculated anchor. The 32 retained-row controls
do not close θ; full anchor covariance raises it. This page retains the
earlier influence comparisons at fixed diagonal errors.

User-requested investigation, 12 November 2026. The additional **independent
2% Conrad modulus-error allowance is retained as a diagnostic**, added in
quadrature to the formal upper widths. It is not a replacement for source
errors, a catalog default or a recovered author setting. All three selected
Conrad points remain in every case. Other error assumptions stay fixed.

**Du's thermal-expansion series provides the strongest upward pull on θ.
Bosenick's low-temperature data provide most of the local constraint and
oppose that pull.** Their balance depends on the normalization assumptions.
Omitting Du overcorrects θ and moves K′ away from the preceding close match;
it does not reproduce the published joint result.

All 30 declared cases converge from both starts without hitting a bound.
No observation, original error, published coefficient or catalog record was
changed, and no author was contacted.

## Series deletions at fixed diagnostic errors

Published θ = **749.825 K**, K′ = **5.06534**.

| Series omitted | θ, plotted normalization, K | θ, linear 298 K normalization, K | K′, linear 298 K |
| --- | ---: | ---: | ---: |
| None, 219-row baseline | 781.184 | 764.605 | 5.06502 |
| Du V–T, 13 rows | **737.757** | **722.787** | **5.02255** |
| Bosenick V–T, 19 rows | **885.133** | **884.160** | 5.07500 |
| Milani V–T, 35 rows | 787.327 | 764.900 | 5.07654 |
| Rodehorst V–T, 9 rows | 777.266 | 760.325 | 5.06338 |
| Isaak GR-2B V–T, 15 rows | 778.339 | 761.383 | 5.06898 |
| Thiéblot V–T, 7 rows | 783.170 | 766.470 | 5.06583 |
| Kono Kₛ, 36 rows | 781.239 | 764.643 | 5.09558 |
| Gwanmesia Kₛ, 52 rows | 781.983 | 765.337 | 5.07597 |
| Isaak GGD Kₛ, 22 rows | 781.355 | 764.622 | 5.08972 |

The [complete table](../data/garnet-angel-2022/grossular/theta-influence-summary.csv)
also includes Bass, Milani compression and the combined Milani/Isaak
publication deletions. These overlapping publication deletions are not
additive to the series effects. Dropping Kono or Gwanmesia removes its unused
modulus multiplier and remaps the remaining indices without changing inputs.

![Theta and Kprime influence, light theme](../data/garnet-angel-2022/grossular/theta-influence.svg#only-light)
![Theta and Kprime influence, dark theme](../data/garnet-angel-2022/grossular/theta-influence-dark.svg#only-dark)

The top row is the full baseline; every other row omits the named series.
Numbers in parentheses give omitted row counts. Dashed lines are the published
central values, not uncertainty limits. The paired panels show that changing
θ alone is not enough to reproduce the full parameter set.

Du's omission lowers θ by **43.43 K** under plotted normalization and
**41.82 K** under linear 298 K normalization. It crosses below the source
value in both cases, and K′ falls to **5.02568/5.02255**. This identifies
influence rather than a justified exclusion. In contrast, removing Bosenick
raises θ by **103.95/119.56 K** and γ to about **1.196–1.197**. These cold
data restrain the high-θ solution; they are not the main upward driver.

Deleting the elasticity series has small effects on θ at this baseline,
although some strongly affect K′. Milani compression is especially important
for K′: omitting its seven rows gives K′ about **4.86**. The earlier diagnosis
of K′ and this diagnosis of θ therefore concern different information in the
same joint fit.

## A fixed-θ control separates residual tension from information

For each normalization, a control fixes θ at 749.825 K and refits K₀, K′,
γ and both modulus multipliers. For the following comparison, effective
pressure errors are frozen at that control's solution. Curve costs then
use identical rows and errors; they are not calibrated significance tests.

| Quantity | Plotted normalization | Linear 298 K normalization |
| --- | ---: | ---: |
| Actual free-fit θ minus published θ, K | +31.359 | +14.780 |
| Local predicted upward θ shift, K | +29.711 | +14.058 |
| Du local upward contribution, K | **+38.830** | **+38.206** |
| Bosenick local contribution, K | **−11.608** | **−27.003** |
| Rodehorst local contribution, K | +4.208 | +4.254 |
| Bosenick conditional local information | 68.76% | 68.74% |
| Milani V–T conditional local information | 15.73% | 15.75% |
| Du conditional local information | 4.64% | 4.65% |

The θ derivative is projected away from all five scaled nuisance columns
before calculating these signed contributions. They sum over all 219 rows;
the remaining groups partly cancel the contributions shown. The information
fractions sum to one, but are **conditional local information**, not causal
shares or fractions of experimental truth. A series can have limited
information yet a large score if its residual pattern strongly opposes a
particular θ. Individual rows within one series can have opposite signs.

The local prediction differs from the actual refitted shift by **1.65 K**
or **0.72 K**. It is a linearization with frozen errors, not an exact
decomposition of the nonlinear IRLS fit. The baseline θ–γ conditional
correlation is **+0.867/+0.862**, so influences must be interpreted in the
joint parameter space rather than as independent θ measurements.

Fixing θ increases the common-error total squared-residual cost by
**9.595/2.157**. Du contributes **+25.294/+11.762**, while Bosenick contributes
**−14.228/−9.776** and several other groups also improve. The individual
group changes therefore must not be read as independent likelihood ratios.

With θ fixed, K′ remains **5.07477/5.06834**; under linear normalization,
γ becomes **1.14961**, versus 1.15602 published. Thus the source θ can be
imposed with modest change in this conditional objective, but doing so
does not recover the full published coefficient vector or validate the
assumed uncertainties. The absolute residual tension with Du remains.

## Which Du observations matter?

The largest upward row contributions occur at **393.7 K** and **447.1 K**.
Under linear normalization they contribute **+15.65 K** and **+17.03 K** in
the fixed-θ local score. The broader **298 ≤ T < 700 K** Du bin contributes
**+40.29 K**, while its four rows at T ≥ 700 K contribute **−2.09 K**.
This points to the middle-temperature shape/offset rather than simply the
highest-temperature measurements. These are score summaries; no temperature
subsets or individual rows were omitted in additional fits.

| Du temperature, K | Raw cell volume token, Å³ | Observed minus fixed-θ molar volume, cm³/mol, linear normalization |
| ---: | ---: | ---: |
| 393.7 | 1666.0(1) | −0.13060 |
| 447.1 | 1668.3(1) | −0.11719 |

Those temperature, volume and printed-error entries agree with the available
[Du author-preprint Table 2](https://academiccommons.columbia.edu/doi/10.7916/D8V415Z5/download),
manuscript p. 25, checked visually. No transcription correction was found.
The final publisher table remains unverified: current retrieval attempts
were inaccessible, which is not evidence of absence. See the
[source-check receipt](../data/garnet-angel-2022/grossular/theta-influence-source-check.json).

The native-volume departures above use the constrained fit's other fitted
coefficients and the adopted normalization; they are not corrections to
observations. The [native-coordinate table](../data/garnet-angel-2022/grossular/theta-influence-native-comparisons.csv)
retains all V/Kₛ comparisons for the published coefficient vector, diagnostic
baseline and fixed-θ control. Each curve uses its own explicitly recorded
modulus multipliers. These unweighted native residuals differ from the
pressure-coordinate objective.

## Normalization interacts strongly with the cold data

The two baseline normalizations differ by **16.580 K** in θ. When Bosenick
is absent, that difference shrinks to **0.973 K**; when Du is absent, it is
still **14.970 K**. This is a paired interaction from already declared fits,
not a new normalization search. It links much of the normalization sensitivity
to the balance involving Bosenick's cold series. It does not uniquely assign
the change to its anchor: all three prescribed denominator changes and the
remaining parameter correlations are still part of the comparison.

Bosenick's original 295 K anchor and its calculated 298 K alternative remain
distinct. Shared anchor errors, instrument/cycle grouping and cross-point
covariances have not been recovered or propagated. Consequently the present
information fractions and influence ordering describe the declared diagonal
error model, not the unknown original author covariance.

## Retained diagnostic and remaining evidence

The additional independent 2% Conrad allowance is retained in the
[predeclared plan](../data/garnet-angel-2022/grossular/theta-influence-plan.json)
and every retained-row receipt. No original error cell was filled or changed.
This keeps the requested diagnostic reproducible while preserving its
separation from measured uncertainty and published parameterizations.

The strongest remaining input questions concern Du's original fit-file
values, assigned errors, normalization and run grouping, together with the
exact Bosenick reference-anchor treatment and shared covariance. The source
files would distinguish those processing choices from a real tension between
the experiments. These findings do not justify deleting either series.
No new error model, EOS law, specimen correction or replacement coefficient
set was selected, and no further author contact occurred.

## Artifacts and verification

- [Both starts and all fit covariances](../data/garnet-angel-2022/grossular/theta-influence-refits.json).
- [5,980 fitted pressure residuals](../data/garnet-angel-2022/grossular/theta-influence-residuals.csv).
- [438 row-level profile scores](../data/garnet-angel-2022/grossular/theta-influence-profile-rows.csv), [group summaries](../data/garnet-angel-2022/grossular/theta-influence-groups.csv) and [temperature summaries](../data/garnet-angel-2022/grossular/theta-influence-temperature-bins.csv).
- [Paired normalization interactions](../data/garnet-angel-2022/grossular/theta-influence-normalization-interactions.csv).

The 30 cases comprise a baseline, a fixed-θ control, eleven non-Conrad
series deletions and two overlapping whole-publication deletions under each
of two fixed normalizations. Two starts give 60 optimizations. Tests verify
retained values/errors, the unchanged Conrad diagnostic, scale remapping,
baseline parity, constrained stationarity and inverse-information covariance,
case replays, cost and score partitions, native forward/inverse agreement,
normalization contrasts and all eight previous artifact manifests. See the
[verification receipt](../data/garnet-angel-2022/grossular/theta-influence-verification.json).

Repository-only commands:

```bash
python -m scripts.refit_grossular_angel_2022_theta_influence
python -m scripts.plot_grossular_angel_2022_theta_influence
pytest -q tests/test_grossular_theta_influence.py
```
