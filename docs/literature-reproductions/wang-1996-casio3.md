# Wang et al. (1996): CaSiO3 fit reproduction

Audit: 2026-10-08. The result is **partial reproduction**, with the remaining
fits retained as explicit diagnostics. All published EOS coefficients, errors,
pressure scales and classifications are preserved.

The [source paper](https://millenia.cars.aps.anl.gov/gsecars/LVP/publication/Papers/Wang_et_al95JB03254.pdf)
prints numerical observations in Table 1 (p. 664), the regression equations on
pp. 664 and 666–667, and thermoelastic results in Table 2 (p. 668). The checks
use the bundled 66-row Wang table, its 14-row room-temperature subset and the
46 printed Mao (1989) rows. Checksums are verified before fitting.

Reproduce the row-free audit artifact with:

```sh
.venv/bin/python scripts/reproduce_wang_1996_casio3.py
```

The result is saved to `docs/data/wang-1996-casio3-refit.json`, including fitted
parameters, independent local-Jacobian errors/covariances, source row orders,
fixed parameters, weights, residuals and unresolved source conventions. These
independent error estimates are not the original source covariance.

## Room-temperature compression and Mao reanalyses

| Fit | Source | Recovered result | Assessment |
| --- | --- | --- | --- |
| Preferred new-data BM3, 12 included Run-13 states; K0′=4.8 fixed | V0=45.58(4) Å³, K0=232(8) GPa | Unit-weight pressure residuals: V0=45.5878 Å³, K0=229.8908 GPa | Both parameters agree within source uncertainties; source weighting and exact errors are not established |
| Equal-weight Mao reanalysis, excluding the two sub-1-GPa states | V0=45.71(12), K0=244(12), K0′=4.8(3) | 44 rows: V0=45.7061, K0=243.7279, K0′=4.7695; independent errors below | Rounded coefficients **and errors reproduced** under unit-weight pressure residuals |
| Pressure-weighted Mao reanalysis excluding sub-1-GPa states | V0=45.47(6), K0=268(8), K0′=4.3(2) | Printed-pressure-error weighting: V0=45.3521, K0=279.1467, K0′=4.0748 | **Not reproduced at published precision**; the historical optimizer/residual convention is unresolved |
| Pressure-weighted all-Mao reanalysis | V0=45.3(5), K0=282(8), K0′=4.0(2) | No uniquely defined all-row pressure-error objective | The 0-GPa row has no printed pressure error; no weight is fabricated |

For the equal-weight reanalysis, the independently calculated standard errors
are **0.12178 Å³, 11.53185 GPa and 0.26399**. They round to the printed
0.12 Å³, 12 GPa and 0.3. This is strong numerical agreement, but does not recover
the original unrounded measurements or covariance matrix.

The generic historical catalog-wide refit ledger did not apply these distinct
Mao selection/weight protocols: it used all 46 rows with errors in both axes for
the three reanalysis records. Its old entries are not evidence that the
source-owned regressions were reproduced. This dedicated audit supersedes that
interpretation; it does not regenerate the entire historical ledger.

## Thermal-pressure equation

Equations (3a–b) are implemented directly as

`P(V,T) = BM3(V; V0,K0,K0′=4.8) + [A - B ln(V0/V)] (T-300)`.

Here A=(∂P/∂T)V and B=(∂K/∂T)V; the source's quadratic temperature term is
fixed to zero. Derived quantities are α0=A/K0 and
`(∂K/∂T)P = B - 4.8 A`. This is the source's phenomenological thermal-pressure
model, not a Debye/MGD substitution.

The printed weighting is

`w = 1/[2 |differential stress|]² + 1/[K0 ΔV/V]²`.

It is a **sum of inverse variances**, not the inverse of their sum. Run 3's
first state (full-table source order 35) has printed differential stress
**0.000 GPa**, making the literal full-table weight singular. No hidden stress
floor is introduced; the unrounded stress underlying that printed zero is
unavailable. For Run 13, the literal formula is finite but gives
K0≈207.46 GPa, far from Table 2's 229(4) GPa for weighted fit 2.

The executable alternatives are labelled **proxy weights**: volume-only
`sqrt(w)=V/(232 ΔV)` and quadrature
`sqrt(w)=1/sqrt((2 stress)² + (232 ΔV/V)²)`. The fixed 232 GPa used to convert
volume error to pressure error is explicit. Proxy agreement is not proof of
the original weighting.

| Table 2 thermal fit | Source free coefficients | Explicit diagnostic | Assessment |
| --- | --- | --- | --- |
| Weighted 1 | V0=45.60(5), K0=233(7), A=0.0075(5), B=0.005(22) | All 66 rows, volume-only proxy: V0=45.5855, K0=237.5335, A=0.00768935, B=0.0138604 | All four fitted coefficients within source error bars, under a proxy protocol |
| Weighted 2 | V0=45.61(2), K0=229(4), A=0.0069(1); B=0 fixed | All 34 Run-13 rows, volume-only proxy: V0=45.6040, K0=229.6981, A=0.00698936 | All three fitted coefficients within source error bars, under a proxy protocol |
| Unweighted | V0=45.58(4), K0=233(9), A=0.0071(1); B=0 fixed | All 66 rows: V0=45.5380, K0=238.6139, A=0.00767810 | A differs beyond its printed uncertainty; **not reproduced** |
| Unweighted mask diagnostic | Same published target | Exclude the two sub-2-GPa rows: V0=45.5778, K0=232.3495, A=0.00755724 | V0/K0 closer, but A still differs; **not a successful reproduction** |

The source excludes the low-pressure amorphizing states in its compression
discussion, while the thermal figure captions state that all Table 1 data are
plotted. It does not explicitly supply a thermal regression row manifest.
Both 66- and 64-row diagnostics are retained; neither is silently declared the
original mask.

The prose's three-term and fixed-B fits are also rerun and kept in the JSON.

## Optional independent thermal EOS

The catalog includes **Peritheos refit of Wang (1996), BM3 + linear thermal
pressure (64 points)** under
`ca_perovskite_wang_1996_unweighted_bm3_linear_thermal_refit`. It is selectable
and nondefault; published records remain unchanged.

This joint equal-weight pressure fit uses the 64-row diagnostic above, with
K0′=4.8, Tr=300 K and B=0 fixed. The exact selected CSV retains the original
row/run identities and printed uncertainties. Rows 33 and 34 are excluded
because the source flags their low-pressure amorphization. This is our
explicit selection, not a recovered original thermal regression manifest.

| Fitted parameter | Value | Conditional standard error |
| --- | --- | --- |
| V0 | 45.577827 Å³ per formula unit | 0.043126 Å³ |
| K0 | 232.349472 GPa | 7.413361 GPa |
| A=(∂P/∂T)V | 0.007557237 GPa/K | 0.000150220 GPa/K |

Pressure RMS is **0.356201 GPa**. The full three-parameter covariance is
stored, scaled by RSS/(64−3). Its errors are conditional on the fixed
parameters, independent common-variance pressure residuals and the chosen
model; they omit coordinate errors, shared correlations and calibration
systematics. They are not published Wang errors.

Observed bounds are **2.66–12.25 GPa and 301–1594 K**, an uneven observation
envelope rather than full rectangular coverage. The fitted 300 K,
zero-pressure reference is an extrapolation, not evidence of ambient phase
stability. Runs 3/4/5 are inherited Wang–Weidner measurements; overlapping
source tables must not be combined as independent observations.

Register this record reproducibly with:

```sh
.venv/bin/python -m scripts.register_wang_1996_thermal_refit
```

For Python selection:

```python
from peritheos import get_eos_record

eos = get_eos_record("ca_perovskite_wang_1996_unweighted_bm3_linear_thermal_refit")
pressure_gpa = eos.pressure(45.0, 1000.0)
```
The prose first reports K0=231(7), V0=45.61(4) and α0=3.33(6)×10⁻⁵ K⁻¹;
this differs from the Table 2 weighted-1 column, which has a free B term. Those
are separate fitting targets. Interchanging them would conceal a source
distinction.

## High-temperature Birch–Murnaghan fits

The source's alternative model is implemented with fixed V0=45.58 Å³,
K0=232 GPa and K0′=4.8, and

`V0(T)=45.58 exp[α0 ΔT + b ΔT²/2]`,
`K0(T)=232 + (∂K/∂T)P ΔT`, `ΔT=T-300`.

This integrates α(T)=α0+b(T-300). The quadratic expansivity coefficient c is
zero; the constant-α fits also fix b=0. All 66 states are used in these
pressure-residual diagnostics.

| Table 2 fit | Source | Diagnostic | Assessment |
| --- | --- | --- | --- |
| Weighted | α0=3.55(18)×10⁻⁵; (∂K/∂T)P=−0.036(8) GPa/K | Volume-only proxy: α0=3.66268×10⁻⁵; derivative=−0.0417237 | Both fitted coefficients within source uncertainties; proxy protocol |
| Unweighted 1, constant α | α0=3.06(14)×10⁻⁵; derivative=−0.022(10) | Unit-weight pressure residuals: α0=3.33120×10⁻⁵; derivative=−0.0278507 | α0 lies beyond the source's 1-sigma error; no exact reproduction |
| Unweighted 2, linear α | α0=3.19(48)×10⁻⁵; derivative=−0.035(15); b=0.45(40)×10⁻⁸ K⁻² | Unit-weight pressure residuals: α0=3.47629×10⁻⁵; derivative=−0.0474960; b=0.686512×10⁻⁸ | All three coefficients within source uncertainties; rounded values and errors differ |

The residual coordinate of the source's unweighted nonlinear BM fits is not
specified. Agreement within parameter intervals is a qualified check, not
identification of that optimizer or its error estimates.

## Interpolation and remaining limits

The isochoric and isobaric columns derive from interpolated coordinates shown
in Figures 5–6. The paper does not publish the numerical interpolation output,
algorithm, regression weights or covariance. Their original interpolation fits
cannot be uniquely recovered from Table 1 alone, so no digitized derived
curve is passed off as an independent experimental observation or original
regression input.

All recovered source rows and published fits remain intact. The unresolved
thermal weights, masks and nonlinear objective are separate from missing
unrounded measurements and missing NaCl calibrant lattice parameters. No
pressure-scale correction or preferred replacement EOS follows from this audit.
