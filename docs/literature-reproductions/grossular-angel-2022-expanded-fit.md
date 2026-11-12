# Grossular: adding the 52 Gwanmesia moduli to the fit

**Subsequent investigation:** the [remaining-channel and normalization audit](grossular-angel-2022-remaining-channels.md)
adds primary error information and a separately declared 193-row diagnostic.
This page retains the earlier 147-row result and its assumptions.

Diagnostic refit, 12 November 2026. This follows the
[observation recovery](grossular-angel-2022-recovery.md) and answers whether
adding the recovered Gwanmesia data brings the earlier fit closer to Angel
et al. (2022), [*Garnet EoS: a critical review and synthesis*](https://doi.org/10.1007/s00410-022-01918-5).

**The expanded 147-channel fit converges and is close to the published
coefficients, but does not reproduce them.** K₀, K′ and γ move slightly closer
than in the earlier 95-channel fit; θ moves farther away. The recovered
Gwanmesia data were already consistent with the published curve, and their
fitted scale factor almost exactly matches the source value.

## Parameter comparison

The original source coefficients and all earlier receipts remain unchanged.
Values below are central values, displayed to at most six significant digits;
full precision and local covariance are in the
[expanded fit receipt](../data/garnet-angel-2022/grossular/expanded-refits.json).

| Quantity | Published source | Earlier 95 channels | Expanded 147 channels | Expanded minus source | Relative difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| K₀T, GPa | 167.002 | 167.373 | 167.353 | +0.350876 | +0.210% |
| K′₀T | 5.06534 | 4.96395 | 4.98267 | −0.0826686 | −1.632% |
| θ₀, K | 749.825 | 744.384 | 740.183 | −9.64216 | −1.286% |
| γ₀ | 1.15602 | 1.16605 | 1.16467 | +0.00864934 | +0.748% |
| Kono calculated-modulus multiplier | 0.99455 | 0.996923 | 0.996303 | +0.00175334 | +0.176% |
| Gwanmesia calculated-modulus multiplier | 1.00551 | Not fitted | 1.00547 | −0.0000432488 | −0.00430% |

The two fitted volume multipliers are **0.999392 for Milani** and **0.998942
for Du**. They are diagnostic nuisance parameters; the original study instead
pre-normalized the volume datasets. Both starts converge to the same solution
without hitting a bound; their largest θ difference is below 0.00001 K.

For context, Table 1 prints source widths of 0.20 GPa, 0.08, 13 K and 0.006
for the four EOS coefficients. θ's central-value difference is smaller than
its printed width; K₀ and γ differ by more than one printed width, while K′
is about one width away. This is **not a significance test**: source confidence
conventions and the original input covariance remain unresolved.

Our local, unscaled diagonal-error fit widths give K₀ = **167.35 ± 0.26 GPa**,
K′ = **4.98 ± 0.08**, θ₀ = **740 ± 25 K** and γ₀ = **1.165 ± 0.009**.
These are conditional inverse-Jacobian widths, without reduced-objective
rescaling. Their statistical coverage is not established because the adopted
coordinate widths are not all known independent standard deviations.

## Declared protocol

The [plan](../data/garnet-angel-2022/grossular/expanded-refit-plan.json) was saved
before solving. It keeps the same BM3 plus q-compromise Debye formulation,
Vₘ₀ = 125.35 cm³/mol, reference temperature 298 K, 20 atoms per formula unit
and EosFit gas-constant convention as the
[95-channel fit](grossular-angel-2022.md#declared-95-channel-partial-refit).
There is no fitted q or alternative EOS law.

The 95 earlier rows are unchanged: 35 Milani expansion, seven recalibrated
Milani compression, 13 Du expansion, 36 Kono moduli, three Isaak moduli with
explicit errors and one Bass modulus. We add all **52 Gwanmesia Kₛ rows**,
with their original NaCl/Decker pressures and **T = T°C + 273.15**. In
particular, the ten 30 °C rows become **303.15 K**. Gwanmesia XRD volumes
remain excluded, following S1.

For this diagnostic, Gwanmesia pressure and modulus weights use **0.7% of P**
and **0.8% of Kₛ**, respectively, from the primary paper's typical precision
and error estimates. Applying those global percentages as diagonal per-row
weights is an explicit assumption, not recovery of the author's assigned
errors. Confidence levels, calibration covariance and temperature errors are
unknown. Missing temperature-error terms are omitted, not treated as known zero.

The objective remains the earlier iteratively reweighted effective-variance
pressure residual. We fit four EOS coefficients and four nuisance parameters:
two volume multipliers and two calculated-modulus multipliers. Each modulus
multiplier maps **calculated Kₛ onto the observed scale**. Source-scale priors
are not imposed. Existing anchor/calibration covariance omissions remain.
Two predeclared starts use identical rows, weights and bounds. No weights,
models or selections were searched to obtain agreement.

The reduced weighted objective is **2.83756**, with **139 nominal degrees of
freedom**, versus 4.36512 for the earlier 95-channel fit. These statistics
belong to different datasets and parameter counts: the decrease does not
establish a better match to the original study or its unknown χ².

## Residual comparisons on the same 147 rows

The [native-coordinate comparison](../data/garnet-angel-2022/grossular/expanded-native-comparisons.csv)
retains every observed value, original row identity, source and fitted
prediction, nuisance multiplier and residual. The sign is **observed minus
calculated**. Source and fit predictions below use the **same fitted nuisance
multipliers**, isolating changes in the four EOS coefficients. Thus the source
column is not the source's original normalization/scaling workflow. These
forward V/Kₛ residuals are distinct from the minimized pressure objective.

| Dataset | Rows | Coordinate | Source RMS | Expanded-fit RMS |
| --- | ---: | --- | ---: | ---: |
| Milani expansion | 35 | V, cm³/mol | 0.02591 | 0.02890 |
| Du expansion | 13 | V, cm³/mol | 0.04488 | 0.04213 |
| Milani compression | 7 | V, cm³/mol | 0.01318 | 0.01139 |
| Kono | 36 | Kₛ, GPa | 1.66829 | 1.49129 |
| Isaak, explicit-error subset | 3 | Kₛ, GPa | 0.51786 | 0.62871 |
| Bass | 1 | Kₛ, GPa | 0.21966 | 0.15745 |
| Gwanmesia | 52 | Kₛ, GPa | 0.71166 | 0.67753 |

With the **original published Gwanmesia multiplier 1.00551**, the source RMS
is **0.71176 GPa**, as recorded in the recovery audit. The expanded fit's
0.67753 GPa is a modest change. The source's close agreement with the recovered
Gwanmesia series helps explain why adding it causes only small coefficient
shifts rather than closing every difference.

## Limits and verification

This is **147 channels, not the full nominal 220**. The original Thiéblot
table/errors, exact normalization and fit masks, errors for additional Isaak
moduli, Conrad covariance and original fit controls remain incomplete.
Mixed-composition GR-2B expansion is not treated as measured pure-grossular
volume, and digitization bounds are not converted into experimental errors.
No catalog record or source coefficient was changed, and no author was contacted.

Tests replay the expanded optimization, check both starts and covariance rank,
and verify all 52 Gwanmesia Kₛ inversions against independent forward Debye
quadrature for both published and fitted coefficients. The earlier 95-channel
fit still replays with the generalized solver. Historical artifact hashes,
literal source tables and the official EOS remain unchanged. See the
[verification receipt](../data/garnet-angel-2022/grossular/expanded-verification.json).

Repository-only maintenance commands:

```bash
python -m scripts.refit_grossular_angel_2022_expanded
pytest -q tests/test_grossular_expanded_refit.py tests/test_grossular_observation_audit.py tests/test_grossular_followup_recovery.py
```

The result supports **close conditional coefficient agreement**, not exact
parameter parity, original uncertainty reproduction or author-workflow reproduction.
