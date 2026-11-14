# Pyrope: numerical differences and fitting improvements

14 November 2026. Following the [author-input replay](pyrope-author-input-replay.md),
the small Python/EosFit differences are traced to **effective-weight propagation,
modulus inversion and parameter derivatives**, with a smaller convergence/
rounding contribution. Covariance inversion and residual rescaling agree.
No source observations, EOS coefficients or public fitting defaults changed.

## Why Python has a slightly higher χ²

Python reports reduced χ² **2.234252069**. The audited EosFit 7.60 run holds
**2.233663082** internally and prints **2.2337**. The difference is
**0.000588987**, or **0.0264%**. These use different residual/weight calculations;
the ordering does not establish which optimizer found a better minimum of one
shared objective.

We captured all 194 native pre-step pressures, weights and five parameter
derivatives in a debugger-launched copy of the official executable. Numerical
values were read without alteration. Instrumented and ordinary replay results
agree. At the common native-cycle starting state, component swaps give:

| Pressure calculation | Weights | Reduced sum, 189 degrees of freedom |
| --- | --- | ---: |
| Python precise | Python local derivatives | 2.234311942 |
| Python precise | Captured native | 2.232245521 |
| Captured native | Python local derivatives | 2.235743739 |
| Captured native | Captured native | 2.233663874 |

The last value reconstructs the native statistic within **8×10⁻⁷**, with
single-precision arithmetic/accumulation still represented by the captured
value. Replacing weights reduces the sum by about **0.00207**; replacing the
pressure calculation then adds about **0.00142** back. The effects partly
cancel. Python's final statistic also reflects its converged coefficients
rather than this initial-cycle state. The [complete analysis](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-analysis.json)
retains all [row comparisons](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-row-comparison.csv)
and dataset totals.

## Error propagation differs numerically

For x equal to volume, modulus or temperature, Python's audit routine uses
the local first-order contribution

\[
\sigma_{P,x}=|\partial P/\partial x|\sigma_x.
\]

The audited native binary instead uses a finite positive displacement:

\[
\sigma_{P,x}=|P(x+\sigma_x)-P(x)|.
\]

Both add independent contributions in quadrature with σ(P)². The finite step
retains curvature over that particular interval and depends on its direction.
It is not integration over an uncertainty distribution or generally a more
accurate probabilistic propagation. For example, P=x², x=10 and σ(x)=1 give
contributions **20** and **21**. The binary's three-error fallback for a change
below single-precision `tiny` is not needed for the nonzero contributions here.

Changing only this propagation calculation, retaining precise inversion,
moves the converged Python θ from **770.865935 to 770.977834 K**. This changes
the objective definition, rather than improving the original objective.

## Modulus inversion and parameter derivatives

The supporting [CrysFML source snapshot](https://code.ill.fr/scientific-software/crysfml/-/blob/8924b41c13a21f642938e1e5f5961302c3103e41/Src/CFML_EoS.f90)
documents `Get_Pressure_K`: its modulus tolerance is **0.001 times the target
modulus**, with additional pressure-step stopping conditions. Its initial
pressure estimate also affects the approximate result. Our original inversion
solves the same operational KS equation much more tightly. Pressure predictions
differ by about **0.01 GPa** in some modulus rows despite the same printed EOS
coefficients. The native parameter derivatives use adaptive finite steps
starting at about 1% of a parameter; our original replay uses small central
differences. Inversion, spline and single-precision effects enter those steps.

A structural double-precision control closely tracks the native pressure and
weights but is explicitly not a bitwise clone. Four sequential native fits
starting from the author EOS and four starting from Python's optimum give
θ **771.05658–771.06836 K**, with printed χ²=2.2337 throughout. This rerun
variation is smaller than the method-dependent gap. The
[evidence receipt](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-evidence.json)
records executable/source checksums, capture addresses and archived scripts.
The public source snapshot is supporting evidence, not the complete recovered
build source of Ross's executable.

## Explaining the fitted parameters and their errors

One ordered sequence of Gauss–Newton component replacements from the same
native-cycle starting state gives:

| Components used in one step | θ after step, K |
| --- | ---: |
| Python pressure, weights and Jacobian | 770.835786 |
| Replace weights with native values | 771.030699 |
| Also replace calculated pressures | 771.056526 |
| Also replace the parameter Jacobian | 771.058156 |
| Actual native stored coefficient | 771.058105 |

All five native parameter updates reconstruct within **1.23 single-precision
last-bit units**. Further Python IRLS iterations move its first-step θ to
770.865935 K. This replacement order is explicit: interacting components do
not have a unique order-independent attribution.

Holding the residual-scale factor fixed separates covariance effects:

| Jacobian / weights | K₀ error relative to native | K′ error relative to native | Scale error relative to native |
| --- | ---: | ---: | ---: |
| Python / Python | −0.376% | −0.504% | −0.387% |
| Python / native | −0.201% | −0.268% | +0.011% |
| Native / Python | −0.174% | −0.235% | −0.397% |
| Native / native | less than 0.000005% | less than 0.000005% | less than 0.000005% |

With the actual native Jacobian and weights, Python reproduces **all five
native errors within 4.4×10⁻⁸ relative**. The matrix inversion and residual
rescaling therefore agree to the available precision. Weight and derivative
calculations explain the earlier approximately 0.5% error differences. θ's
error is much less sensitive than those of K₀, K′ and the scale. The residual
factor itself differs by only about 0.026% in variance, or 0.013% in widths.

## A tested improvement that preserves the Python replay objective

The repository-only `scripts/refit_pyrope_author_analytic.py` prototype
calculates BM3/Debye derivatives directly and differentiates the KS inversion
implicitly. It retains precise inversion, local derivative weights, all 194
inputs and the same free/fixed parameters. It checks weighted-gradient
stationarity, bounds and the scaled Jacobian condition number.

Independent five-point checks verify every parameter and coordinate derivative.
Coefficients, χ² and errors agree with the existing Python replay; the final
normalized stationarity score is about **3.3×10⁻¹⁰**. Five paired alternating-
order local runs measured medians of **112.4 ms** for the finite-difference
replay and **17.3 ms** for the analytic prototype: about **6.5× faster** for
this complete fit. This is an audit-routine benchmark, not a general public
library or Rust performance claim. [Results and timings](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-analytic-fit.json)
retain both starts and every measurement.

**Scope clarification:** the public Peritheos fitter already adjusts uncertain
coordinates as latent variables and supports observation covariance. It has
an errors-in-variables objective. The pressure-residual IRLS routine used for
this author replay is a separate repository tool. Our earlier conversation's
derivative-propagation description applies to that replay, not the public
API's default objective. No public fitter changed here. An explicit EosFit
compatibility option could help reproduction; these findings do not justify
adopting one-sided propagation or coarser inversion as library defaults.

## Reproducibility and limits

The [native kernel fixture](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-native-kernel.csv)
contains actual pressures, weights and Jacobians. The private archive retains
capture scripts/macros, program outputs, source snapshot and disassembly,
version-locked by hashes and addresses in the evidence receipt. Analysis is
executable from tracked fixtures without a debugger. See the
[verification receipt](../data/garnet-angel-2022/pyrope/numerics-2026-11-14-verification.json).
These numerical findings do not establish raw-coordinate independence, full
preprocessing history, total experimental confidence or a replacement EOS.
Earlier author-input/grossular artifacts remain unchanged. No author message
was sent.

Repository-only commands:

```bash
python -m scripts.audit_pyrope_numerical_differences
python -m scripts.refit_pyrope_author_analytic
pytest -q -W error tests/test_pyrope_numerical_differences.py
```
