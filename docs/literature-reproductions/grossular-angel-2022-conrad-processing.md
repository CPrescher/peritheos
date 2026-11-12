# Grossular: Conrad processing and assigned errors

**Retained diagnostic:** the [subsequent θ influence study](grossular-angel-2022-theta-influence.md)
keeps the independent 2% Conrad allowance fixed and all three Conrad points
present. It identifies opposing Du and Bosenick influences without adopting
the estimated errors as source measurements or replacing published coefficients.

User-requested controls, 12 November 2026, following the
[dataset influence study](grossular-angel-2022-dataset-influence.md).

**An additional independent modulus-error allowance can explain the central
K′ discrepancy without removing any Conrad points. A common calibration
error of the same size cannot.** This is conditional sensitivity, not recovery
of Ross's original error assignments. All 26 declared cases retain the same
219 observations, including all three selected Conrad points. All 52
optimizations converge without hitting a bound.

## Results with every selected point retained

Published K′ = **5.06534**. Background errors, specimen treatment and equations
are unchanged; only the declared Conrad treatment changes.

| Conrad treatment | K′, plotted normalization | K′, linear 298 K normalization |
| --- | ---: | ---: |
| Baseline formal-error upper bound | 5.23585 | 5.23304 |
| Formal C₁₁/C₁₂ correlation −1 | 5.23673 | 5.23391 |
| Formal C₁₁/C₁₂ correlation 0 | 5.23629 | 5.23347 |
| Additional independent 0.5% Kₛ width | 5.11815 | 5.11523 |
| Additional independent 1% Kₛ width | 5.08122 | 5.07833 |
| **Additional independent 2% Kₛ width** | **5.06789** | **5.06502** |
| Additional fully shared 2% Kₛ width | 5.18882 | 5.18661 |
| Assign σ(P) = 0.2 GPa | 5.13279 | 5.12993 |
| Assign σ(P) = 0.5 GPa | 5.07682 | 5.07395 |
| Round Kₛ to 0.1 GPa | 5.23632 | 5.23351 |
| Density control from calculated vₚ | 5.24215 | 5.23933 |
| Density control from calculated vₛ | 5.23930 | 5.23649 |
| Freely fit one Conrad scale | 5.18834 | 5.18614 |

The extra modulus widths are added in quadrature to the unchanged formal
widths. The independent 2% case corresponds to approximately **3.49, 3.74
and 4.11 GPa** extra widths for the three points. Those are hypothetical
comparison widths, not recovered measured standard deviations. The pressure
controls explicitly depart from S1's assigned **0.1 GPa**; neither is a
recovered author setting. No fraction or pressure width was adjusted after
seeing the results.

![Conrad controls, light theme](../data/garnet-angel-2022/grossular/conrad-processing.svg#only-light)
![Conrad controls, dark theme](../data/garnet-angel-2022/grossular/conrad-processing-dark.svg#only-dark)

These are fitted central values, not confidence intervals or statistical
significance tests. The [complete parameter table](../data/garnet-angel-2022/grossular/conrad-processing-summary.csv)
also retains K₀, θ, γ and the modulus multipliers.

## What the sources say about preprocessing

[Angel et al. (2022)](https://doi.org/10.1007/s00410-022-01918-5), p. 4,
explicitly state that Kₛ values were **not rescaled before fitting**. They
describe fitting multiplicative scale factors to selected modulus datasets,
with factors multiplying calculated moduli. Their adiabatic–isothermal
conversion occurs inside EosFit, p. 6. Thus an undocumented pre-normalization
of Conrad, or an external conversion to Kₜ followed by treatment as Kₛ, is
not supported by the stated procedure. The extra Conrad scale tested here
is a diagnostic control, not a recovered final author parameter.

[Conrad et al. (1999)](https://doi.org/10.2138/am-1999-0321) describe the
Christoffel inversion and assumed **1% velocity precision**, p. 377. Their
cubic Kₛ formula is (C₁₁ + 2C₁₂)/3, p. 378. Table 6, p. 379, reports formal
one-standard-deviation Cᵢⱼ errors, without their covariance. Changing the
admissible C₁₁/C₁₂ correlation alone moves K′ by less than **0.001** here.
The assigned pressure term dominates these original small formal weights.

The source's 1% velocity precision motivates an error-size control but does
not establish a 2% Kₛ standard deviation. At fixed density, multiplying all
velocities at one pressure by (1+ε) multiplies the inferred stiffness tensor
and Kₛ by (1+ε)², giving a first-order 2ε modulus change. Whether such errors
vary independently between pressures or share a calibration is a separate
assumption. A full inversion error requires the actual velocity covariance,
orientation and density uncertainties, not just this scaling argument.

A later primary experiment,
[Jiang, Speziale and Duffy (2004)](https://doi.org/10.1029/2004JB003081),
paragraph 17, compares its 36 measured directions with Conrad's two and
suggests orientation errors as a possible explanation for their differing
pressure derivatives. This supports examining systematic uncertainty; it
supplies neither a numerical correction to Conrad nor Ross's fitted errors.

The inspected EosFit 7.60 help pages describe input uncertainties and choice
of weighting. They do not recover a particular Conrad assignment. They are
also not the summer-2021 executable used by Angel, so they cannot establish
its exact defaults or the author's chosen input settings.

## Independent scatter and shared calibration behave differently

For the independent cases, the additional Kₛ variances are diagonal. Each
point has more latitude to depart from the EOS, reducing the influence of
Conrad's pressure-dependent residual pattern.

For the shared control, the extra effective pressure covariance is

\[
 C_{ij}=\delta_{ij}\sigma_{P,\mathrm{eff},i}^{2}+u_i u_j,
 \qquad u_i=\left(\frac{\partial P}{\partial K_s}\right)_i 0.02K_{s,i}.
\]

This is a fully correlated 2% modulus calibration allowance, added to the
existing diagonal errors for the three Conrad rows. Each IRLS iteration
freezes this block and whitens its residuals by a Cholesky factor; the other
216 rows keep their existing weights. This generalized covariance control
is not asserted to be the author's fitting procedure.

The shared allowance permits a common shift of the series but does not
release the three pressure-dependent departures independently. Its K′ stays
near **5.189**, similar to the **5.188** obtained with a freely fitted
Conrad multiplier. The latter multiplier is **1.003372** under plotted
normalization. A simple uniform scale adjustment therefore does not resolve
the discrepancy under these controls.

## Density and rounding checks

Table 4's **calculated**, rather than measured, [100] velocities and Table 6's
C₁₁/C₄₄ allow an approximate reconstruction of the density used in the
original velocity calculation:

\[
 \rho_{p}=C_{11}/v_{p,\mathrm{calc}}^2,\qquad
 \rho_{s}=C_{44}/v_{s,\mathrm{calc}}^2.
\]

The controls replace that inferred density with
**3.60 × V₀/Vpublished(P,298)**, using Table 2's reference-density value and
the published EOS volume ratio. Both Kₛ and its formal width are multiplied
by the new/old density ratio. These are calculated central-value controls;
they use the published EOS and cannot independently validate it.

The inferred ratios change Kₛ by only **−0.057% to +0.063%**. Every new density
also lies inside the range allowed by rounding the printed velocities to
0.01 km/s. Hence this does not resolve a density mismatch at the source's
printed precision. The two density controls leave K′ high and slightly
increase it. Rounding Kₛ to 0.1 GPa is likewise insufficient. Unknown
density/calibration covariance is not invented or propagated as recovered
measurement uncertainty.

The [density receipt](../data/garnet-angel-2022/grossular/conrad-processing-density-controls.csv)
retains inferred densities, velocity-rounding ranges and calculated
multipliers. It does not add any new measured channels.

## Interpretation and remaining inputs

The K′ discrepancy can be closed to within about **0.003** while retaining
all three Conrad points, using a predeclared independent 2% modulus-error
allowance. That is evidence for a weighting sensitivity, not proof that Ross
assigned those errors or that the allowance is experimentally correct.

The full source parameter set remains unreproduced. In that case, θ is
**781.184 K** with plotted normalization and **764.605 K** with linear 298 K
normalization, versus **749.825 K** published. The latter fit has
K₀=**166.929 GPa**, γ=**1.15479**, Kono scale=**0.994752** and Gwanmesia
scale=**1.005843**, so a close K′ is not complete coefficient parity.

The most useful missing author input is now the **actual Conrad Kₛ error
column and any additional error assignment**, including whether it represented
independent scatter. The original fit file/log would distinguish that from
an unsupported correction to the central moduli. No further author message
was drafted or sent. All source data, original errors and published coefficients
remain unchanged; none of these controls is a recommended replacement EOS.

## Reproduction and verification

The [plan](../data/garnet-angel-2022/grossular/conrad-processing-plan.json)
was saved before solving: 13 treatments under each of two fixed normalizations,
two starts each. The [fit receipt](../data/garnet-angel-2022/grossular/conrad-processing-refits.json)
retains every case, both starts, covariance, convergence and the three-row
error block. [Input reductions](../data/garnet-angel-2022/grossular/conrad-processing-inputs.csv)
separate literal moduli, rounding controls and calculated density reductions.

The [residual table](../data/garnet-angel-2022/grossular/conrad-processing-residuals.csv)
retains marginal errors and whitened residuals. For shared errors, the sum
of squared whitened residuals is the objective; dividing each residual by
its marginal error would omit the covariance. Tests check this independently
against dense GLS and the rank-one inverse identity, and check stationarity,
replays, retained rows, covariance rank and preservation of all seven previous
artifact manifests. See the [verification receipt](../data/garnet-angel-2022/grossular/conrad-processing-verification.json)
and [source evidence](../data/garnet-angel-2022/grossular/conrad-processing-source-evidence.json).

Repository-only commands:

```bash
python -m scripts.refit_grossular_angel_2022_conrad_processing
python -m scripts.plot_grossular_angel_2022_conrad_processing
pytest -q tests/test_grossular_conrad_processing.py
```
