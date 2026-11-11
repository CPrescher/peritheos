# Multi-oscillator Gruneisen thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Multi-oscillator Gruneisen thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`MultiOscillatorGruneisenThermalEOS` combines an independently selected
reference isotherm, volume-dependent oscillator temperatures, and quadratic
anharmonic and electronic pressure terms. It needs pressure, bulk modulus, and
$dK/dP$ from the reference component. `EosBase` supplies the last quantity by
central numerical differentiation when the selected isotherm has no specialized
implementation. The published Sokolova et al. catalog records compose the
thermal correction with Holzapfel; another reference isotherm is a valid API
composition, but a new scientific model rather than the published scale.

## Paper versus spreadsheet

> **Reproduction warning:** Peritheos follows the calculation path and
> numerical behavior of the Excel workbook distributed with Sokolova et al.
> (2016), not a literal transcription of the equations printed on page 163 of
> the journal article. The article is the bibliographic and conceptual source;
> the workbook is the computational reference where the two differ.

The differences relevant to this implementation are:

- The printed expression for the pressure derivative of the Holzapfel bulk
  modulus (paper equation 5) is not used directly because it does not reproduce
  the workbook values. Peritheos obtains $K_T'=\partial K_T/\partial P$
  numerically from the implemented Holzapfel pressure and bulk modulus, matching
  the spreadsheet-derived regression values.
- The characteristic temperatures are not evaluated by inserting the printed
  equations 9 and 10 as a closed-form expression. The workbook calculation
  constructs the volume multiplier from the integral
  $\exp[\int_x^1\gamma(u)/u\,du]$. Intrinsic anharmonicity enters the pressure
  through its separate quadratic-temperature term.
- Rather than copying the compact paper equation 12 literally, the workbook
  path evaluates each active oscillator contribution explicitly at $T$ and
  $T_r$, subtracts the reference-temperature contributions, and then adds the
  anharmonic and electronic terms. This ensures
  $\Delta P_{\mathrm{th}}(V,T_r)=0$ for the reference isotherm used by
  Peritheos. The complete workbook expression also permits the optional Bose
  modes represented by `QBo`, `QB1o`, `d`, `d1`, `mb`, and `mb1`.
- The workbook calculation uses kbar for $K_0$ and bar for intermediate
  pressures. Peritheos accepts and returns GPa and applies the conversions
  explicitly.

The Gruneisen expression below is algebraically equivalent to the paper's
equation 11 when `beta=0`; the optional `beta` term retains the generalized
workbook form. Although the article states that its equations 10 and 11 correct
typos in earlier publications, that statement does not make the article and
workbook calculation paths identical.

- `QE1o`, `QE2o` and `mE1`, `mE2` are the two reference Einstein
  temperatures and their multiplicities.
- `QBo`, `QB1o`, `d`, `d1` and `mb`, `mb1` define two optional generalized
  Bose modes. Both multiplicities default to zero.
- `delta`, `t`, and `beta` control the Gruneisen function.
- `a_0`, `m` and `e_0`, `g` control the anharmonic and electronic terms.

With $x=V/V_0$, define the generalized parameter

\[
\widetilde t(x)=t-\beta x^{1/3}
\]

and the model Gruneisen function

\[
\gamma(x)=\delta+
\frac{-3K_T+2P\widetilde t+9K_TK_T'-6\widetilde tK_T}
{6(3K_T-2P\widetilde t)}.
\]

Here $P$, $K_T$, and $K_T'$ come from the selected reference isotherm at
$V=xV_0$. The Sokolova catalog compositions select Holzapfel. Each
characteristic temperature follows

\[
\Theta_i(x)=\Theta_{i0}\exp[I_\gamma(x)],
\qquad
I_\gamma(x)=\int_x^1\frac{\gamma(u)}{u}\,du.
\]

The Einstein-mode energy in kelvin is

\[
\varepsilon_E(\Theta,T)=
\frac{\Theta}{2}+\frac{\Theta}{\exp(\Theta/T)-1},
\]

and the generalized Bose-mode energy with dispersion $d$ is

\[
\varepsilon_B(\Theta,T,d)=
\frac{\Theta(d-1)}{2d}
+\frac{T\Theta d}
{(Td+\Theta)\left[\left(1+\frac{\Theta}{Td}\right)^d-1\right]}.
\]

For a mode with multiplicity $m_i$, its pressure contribution in bar is

\[
P_i(V,T)=\frac{m_iR\,\varepsilon_i[\Theta_i(V),T]\,\gamma(V)}{V}.
\]

The final referenced pressure is

\[
\Delta P_{\mathrm{th}}=10^{-4}\left\{
\sum_i[P_i(V,T)-P_i(V,T_r)]
+C_2(V)(T^2-T_r^2)
\right\},
\]

where

\[
C_2(V)=\frac{3nR}{2\times10^6V}
\left[a_0m\,x^m+e_0g\,x^g\right].
\]

The optional Bose-mode multiplicities default to zero, which recovers the
reduced configuration. This class provides mechanical thermal pressure but not
a complete caloric potential; consequently it does not expose $C_V$, $C_P$,
or $K_S$.

---

**End of Multi-oscillator Gruneisen thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
