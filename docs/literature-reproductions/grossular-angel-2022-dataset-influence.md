# Grossular: which datasets drive K′?

User-requested diagnostic, 12 November 2026, following the
[normalization/specimen tests](grossular-angel-2022-normalization-specimen.md).

**The three retained Conrad modulus points are the dominant driver of the
high K′ in our reconstructed fit.** Omitting only those points changes
K′ from **5.23585 to 5.06290**, close to the published **5.06534**.
The same result holds under the independent linear 298 K normalization:
**5.23304 to 5.06004**. This identifies sensitivity to a specific input series;
it does not show that Ross omitted it or that the measurements are wrong.

All 32 declared cases converge from both starts without hitting bounds.
No source observation, assigned error, published coefficient or catalog
record was changed. No author was contacted.

## Deleting series while holding the background fixed

| Series omitted | Rows retained | K′, plotted normalization | K′, linear 298 K normalization |
| --- | ---: | ---: | ---: |
| None | 219 | 5.23585 | 5.23304 |
| Conrad Kₛ | 216 | **5.06290** | **5.06004** |
| Du V–T | 206 | 5.19873 | 5.19587 |
| Kono Kₛ | 183 | 5.25924 | 5.25630 |
| Gwanmesia Kₛ | 167 | 5.26045 | 5.25734 |
| Milani P–V | 212 | 5.25375 | 5.25106 |
| Isaak GGD Kₛ | 197 | 5.40289 | 5.40651 |
| Isaak GR-2B V–T | 204 | 5.24009 | 5.23678 |
| Published K′, for reference | — | 5.06534 | 5.06534 |

The [complete parameter table](../data/garnet-angel-2022/grossular/dataset-influence-summary.csv)
also retains the other series and combined-publication deletions. The latter
overlap the individual series and must not be added to their effects.

![Dataset deletions, light theme](../data/garnet-angel-2022/grossular/dataset-influence.svg#only-light)
![Dataset deletions, dark theme](../data/garnet-angel-2022/grossular/dataset-influence-dark.svg#only-dark)

The top row is the full baseline; every other row omits the named series.
Numbers in parentheses are omitted row counts. These are refitted central
values, not confidence intervals. Similar blue and orange positions show
that the principal influence directions persist under both normalizations.

Du's thermal expansion has a smaller upward influence through the coupled
EOS parameters: its omission lowers K′ by about **0.0371** in either case.
Omitting Kono, Gwanmesia or Milani compression instead raises K′; those series
are not driving the excess in this deletion comparison. Removing Isaak's
GGD moduli raises K′ by **0.1670–0.1735**. These ambient-pressure moduli help
constrain K₀ and the thermal parameters, and thereby restrain K′ in the joint
fit. Their influence is not the same as that of the GR-2B expansion specimen.

## Check the Conrad inputs before interpreting the deletion

The literal C₁₁/C₁₂ tokens and errors were checked against the final
[Conrad et al. (1999)](https://doi.org/10.2138/am-1999-0321), Table 6, p. 379.
There is no transcription mismatch in these values. For cubic grossular,
the derived adiabatic bulk modulus is Kₛ = (C₁₁ + 2C₁₂)/3.

| Pressure, GPa | Derived Kₛ, GPa | Published EOS Kₛ, GPa | Observed minus published, GPa |
| ---: | ---: | ---: | ---: |
| 1.60 | 174.53667 | 176.22819 | −1.69152 |
| 3.50 | 187.09333 | 185.67646 | +1.41688 |
| 6.95 | 205.39333 | 202.56491 | +2.82842 |

The calculated values use the source EOS with its published coefficients,
without a Conrad scale multiplier, at the tabulated pressures and adopted
298 K. The opposite signs at low and high pressure favour a steeper pressure
dependence in this joint fit. At fixed published K′, the 6.95 GPa row supplies
the largest local upward contribution; 3.5 GPa supplies a smaller positive
contribution, while 1.6 GPa opposes it. This is a local diagnostic, not an
additional row-selection search.

The final [Angel et al. (2022)](https://doi.org/10.1007/s00410-022-01918-5),
Fig. 3 caption, p. 8, explicitly excludes the two highest-pressure Conrad
points. Table S1 reports five available and three used, with σ(P) assigned
as 0.1 GPa. Our retained 1.6, 3.5 and 6.95 GPa points follow that mask;
the 7.7 and 10 GPa rows remain excluded in every fit. The ambient Bass row
is not counted again as a Conrad measurement. Thus a mistaken inclusion of
the two published outliers does not explain our K′ discrepancy.

The original C₁₁/C₁₂ covariance and the exact Kₛ errors entered by Ross remain
unrecovered. Our comparison retains the previously declared conservative
covariance-compatible upper widths **0.030, 0.090 and 0.0467 GPa**, plus the
S1-assigned 0.1 GPa pressure width. At the fixed-K′ fit, the pressure term
accounts for **96.8–99.6%** of these three effective pressure-error variances.
Changing C₁₁/C₁₂ correlation within its admissible range can only reduce
the small modulus-error term at that state; it does not provide a larger
uncertainty that would substantially weaken Conrad's weight. This statement
does not exclude additional experimental/reduction errors missing from the
formal table or different author-assigned Kₛ widths.
Unknown inter-point pressure/calibration correlations are also absent from
the present diagonal-error model and could change these conditional influences.

The [native-coordinate receipt](../data/garnet-angel-2022/grossular/dataset-influence-conrad-native.csv)
records predictions for all five original points, with the two excluded rows
explicitly identified, and compares the published EOS, the full fit and the
fit omitting Conrad. The [source-check receipt](../data/garnet-angel-2022/grossular/dataset-influence-source-check.json)
records primary locations and cached-source checksums.

## Constraining K′ isolates the disagreement

For each normalization, a separate control fixes K′ at **5.06534** and
refits K₀, θ, γ and the Kono/Gwanmesia multipliers. At those fitted coefficients,
we freeze the effective errors and compare all curves on exactly the same
219 rows. Under plotted normalization, the sum of squared standardized
residuals for Conrad is **50.38**, versus **28.57** for the free-K′ baseline.
The corresponding difference for Du is **+5.29**; several other series
improve when K′ is fixed. These group costs distinguish coefficient tension
from the much larger total residual cost of Du, which alone is not a measure
of influence on K′.

The local derivative calculation projects the K′ column away from the five
nuisance columns before calculating group scores. Its summed predicted
upward K′ shift is **+0.17045**, close to the actual **+0.17051**. Conrad's
signed contribution is **+0.15366**, including **+0.13237** from the 6.95 GPa
row; Du contributes **+0.02949**, and other groups partly cancel these.
The corresponding linear-298 K results have the same principal ordering.

These are a frozen-error linearization and a conditional information
partition, not a causal decomposition, a statistical significance test or
fractions of experimental truth. A group's local score at constrained K′
need not have the opposite sign of its nonlinear deletion effect. For example,
Isaak GGD has a small positive local score at the constrained fit, but its
deletion raises K′ substantially through parameter correlations. Deletions,
scores, covariance and raw residual sizes answer different questions.
The [group table](../data/garnet-angel-2022/grossular/dataset-influence-groups.csv)
and [row-level profile table](../data/garnet-angel-2022/grossular/dataset-influence-profile-rows.csv)
retain these quantities with their signs and definitions.

## What this resolves and what remains

We have identified the strongest dataset influence on K′ under two fixed
normalization choices. Removing Conrad also moves the fitted modulus
multipliers close to the published values: under plotted normalization,
Kono becomes **0.994824** versus **0.99455**, and Gwanmesia **1.005860** versus
**1.00551**. This is useful diagnostic evidence, not a justified exclusion.

The full source parameter set remains unreproduced. Without Conrad, θ is
still **781.233 K** under plotted normalization and **764.653 K** under the
linear 298 K normalization, versus **749.825 K** published. K₀ also differs
slightly, and the author workflow is not recovered.

Before any later author request, the specific unresolved inputs are now:

- The actual Conrad Kₛ values and errors in the original fit file.
- Any preprocessing, normalization, calibration correction or assigned
  additional error applied to those moduli before regression.
- Confirmation of the pressure values and assigned pressure widths, as well
  as the exact joint fit sequence and point mask.

No new message was drafted or sent. The earlier pyrope correspondence
remains waiting for a reply. We retain the original observations and source
coefficients and do not recommend the deletion fits for prediction.

## Reproduction and verification

The [predeclared plan](../data/garnet-angel-2022/grossular/dataset-influence-plan.json)
contains two normalizations, 12 disjoint series deletions, two overlapping
whole-publication deletions, a baseline and a fixed-K′ control per normalization:
**32 cases and 64 optimizations**. All background observations, estimated errors,
equations and reference conventions are inherited without adjustment from
the preceding study. Deleting Kono or Gwanmesia also removes its unused scale
parameter and remaps the remaining scale indices. The fixed-K′ controls have
five free parameters; the other cases have five or six. Covariances remain
conditional on the diagonal error model and do not establish physical validity.

The [fit receipt](../data/garnet-angel-2022/grossular/dataset-influence-refits.json)
retains both starts, free/fixed parameters, covariances, rank, bounds and
convergence. The [pressure residuals](../data/garnet-angel-2022/grossular/dataset-influence-residuals.csv)
reconcile with every fitted objective. Tests check unchanged retained rows,
scale remapping, baseline parity, constrained-fit stationarity, deletion
replays, covariance rank, local-score partitions, exported costs and
independent forward/inverse agreement on the Conrad moduli. All six previous
artifact manifests remain unchanged. See the
[verification receipt](../data/garnet-angel-2022/grossular/dataset-influence-verification.json).

Repository-only commands:

```bash
python -m scripts.refit_grossular_angel_2022_dataset_influence
python -m scripts.plot_grossular_angel_2022_dataset_influence
pytest -q tests/test_grossular_dataset_influence.py
```
