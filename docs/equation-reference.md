# Equation reference

This page defines the equations evaluated by Peritheos. Constructor parameter
names are shown in code font; mathematical symbols use the corresponding
subscripts. The definitions here follow the implementation and the numerical
reference cases described under [Validation](validation.md). Original model
publications are collected under [References](references.md).

## Notation and reference state

For an isothermal equation of state (EOS), $V_0$ is the volume at zero model
pressure and $K_0$ is the isothermal bulk modulus there. Peritheos uses

\[
K_T=-V\left(\frac{\partial P}{\partial V}\right)_T,
\qquad
K_0'=\left(\frac{\partial K_T}{\partial P}\right)_{P=0},
\qquad
K_0''=\left(\frac{\partial^2K_T}{\partial P^2}\right)_{P=0}.
\]

Thus $K_0'$ is dimensionless and $K_0''$ has inverse-pressure units. Except
for Holzapfel and thermal models, $V$ and $V_0$ may use any internally
consistent volume unit. Public pressure and modulus values are documented in
GPa. Every model requires $V>0$, with additional domains noted below.

## Isothermal equations

### Birch-Murnaghan family

Define the Eulerian finite strain

\[
f_E=\frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right].
\]

The second-order model (`BM2`) is

\[
P=3K_0 f_E(1+2f_E)^{5/2},
\]

which is the third-order model with the constraint $K_0'=4$. The general
third-order model (`BM3`) is

\[
P=3K_0 f_E(1+2f_E)^{5/2}
\left[1+\frac{3}{2}(K_0'-4)f_E\right].
\]

For the fourth-order model (`BM4`), define

\[
\zeta=\frac{3}{4}(4-K_0'),
\qquad
\xi=\frac{3}{8}\left[
K_0K_0''+(K_0'-4)(K_0'-3)+\frac{35}{9}
\right].
\]

Then

\[
P=3K_0 f_E(1+2f_E)^{5/2}
\left(1-2\zeta f_E+4\xi f_E^2\right).
\]

Higher order does not automatically imply better extrapolation. In particular,
$K_0''$ is often strongly correlated with the lower-order parameters unless
the pressure range constrains curvature well.

### Natural-strain family

The positive compressive natural strain is

\[
f_N=\frac{1}{3}\ln\left(\frac{V_0}{V}\right).
\]

All three natural-strain models use

\[
P=3K_0\frac{V_0}{V}
\left(f_N+A f_N^2+B f_N^3\right).
\]

Their coefficients are

\[
\begin{array}{c|cc}
\text{model} & A & B \\
\hline
\text{NaturalStrain2} & 0 & 0 \\
\text{NaturalStrain3} & \frac{3}{2}(K_0'-2) & 0 \\
\text{NaturalStrain4} & \frac{3}{2}(K_0'-2) &
\frac{3}{2}\left[K_0K_0''+1+(K_0'-2)+(K_0'-2)^2\right]
\end{array}
\]

`NaturalStrain2` therefore has the implied reference derivative $K_0'=2$.

### Murnaghan

The Murnaghan EOS assumes a bulk modulus linear in pressure:

\[
K_T(P)=K_0+K_0'P.
\]

For $K_0'\ne0$, its pressure-volume relation and analytic inverse are

\[
P=\frac{K_0}{K_0'}
\left[\left(\frac{V_0}{V}\right)^{K_0'}-1\right],
\qquad
V=V_0\left(1+\frac{K_0'}{K_0}P\right)^{-1/K_0'}.
\]

Peritheos also supports the continuous $K_0'=0$ limit:

\[
P=K_0\ln\left(\frac{V_0}{V}\right).
\]

### Vinet

With

\[
x=\left(\frac{V}{V_0}\right)^{1/3},
\qquad
\eta=\frac{3}{2}(K_0'-1),
\]

the Vinet pressure is

\[
P=3K_0\frac{1-x}{x^2}\exp[\eta(1-x)].
\]

The analytic bulk modulus used by the implementation is

\[
K_T=K_0x^{-2}\exp[\eta(1-x)]
\left[1+(1+\eta x)(1-x)\right].
\]

### Modified Tait

Let $v=V/V_0$ and define

\[
a=\frac{1+K_0'}{1+K_0'+K_0K_0''},
\]

\[
b=\frac{K_0'}{K_0}-\frac{K_0''}{1+K_0'},
\qquad
c=\frac{1+K_0'+K_0K_0''}
{K_0'^2+K_0'-K_0K_0''}.
\]

The modified Tait pressure and analytic inverse are

\[
P=\frac{1}{b}\left[
\left(\frac{v+a-1}{a}\right)^{-1/c}-1
\right],
\]

\[
V=V_0\left[1-a+a(1+bP)^{-c}\right].
\]

The real-valued implementation requires $(v+a-1)/a>0$ and rejects singular
coefficient combinations. When $K_0''=0$, the equation reduces to the
Murnaghan EOS.

### Holzapfel

The Holzapfel implementation requires $V_0$ and $V$ in
`J bar^-1 mol^-1`, $K_0$ in GPa, $n$ as the number of atoms in
the formula unit, and $Z$ as its atomic-number parameter. For a compound,
$Z$ may be the effective value specified by its source parameterization;
Sokolova et al. equation 3 gives that value for MgO. Define

\[
x=\left(\frac{V}{V_0}\right)^{1/3},
\qquad
P_{\mathrm{FG0}}=1003.6
\left(\frac{Zn}{10V_0}\right)^{5/3},
\]

\[
c_0=-\ln\left(\frac{3K_0}{P_{\mathrm{FG0}}}\right),
\qquad
c_2=\frac{3}{2}(K_0'-3)-c_0.
\]

The pressure is

\[
P=3K_0\exp[c_0(1-x)]
\left(x^{-5}-x^{-4}\right)
\left(1+c_2x-c_2x^2\right).
\]

`bulk_modulus_derivative(V)` evaluates $\partial K_T/\partial P$ numerically
and is used by `MultiOscillatorGruneisenThermalEOS`.

## Thermal equations

Thermal models use $V$ in `J bar^-1 mol^-1`, $T$ in K, energy in
`J mol^-1`, and pressure in GPa. Their total pressure is

\[
P(V,T)=P_{\mathrm{ref}}(V)+\Delta P_{\mathrm{th}}(V,T),
\qquad
\Delta P_{\mathrm{th}}(V,T_r)=0.
\]

Because energy divided by the public molar-volume unit produces bar, the factor
$10^{-4}$ converts thermal pressure to GPa.

### Temperature-dependent reference state

`ThermalReferenceStateEOS` composes an independently selected reference
isotherm that exposes reconstructable `V0` and `K0` parameters. At each
temperature it evaluates that same isotherm after applying

\[
V_0(T)=V_0(T_r)\exp[\alpha_0(T-T_r)],
\qquad
K_0(T)=K_0(T_r)+(T-T_r)\left(\frac{\partial K}{\partial T}\right)_P.
\]

The stored parameters are `Tr`, `alpha0`, and `dK_dT`. Unlike the
energy-based thermal models, this construction inherits the volume unit of its
reference isotherm. The bundled ice VI and ice VII records compose it with BM2,
following equations (1)--(3) and Table II of Bezacier et al. (2014),
[doi:10.1063/1.4894421](https://doi.org/10.1063/1.4894421). The `.eosmat`
interchange type remains `AlphaKT` for Dioptas compatibility, while the stable
model identifier is `thermal_reference_state`.

### Mie-Gruneisen Debye and Einstein

Both Mie-Gruneisen models use

\[
\gamma(V)=\gamma_0\left(\frac{V}{V_0}\right)^q,
\qquad
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

Here `Tr` is $T_r$, `theta0` is $\Theta_0$, `gamma0` is $\gamma_0$,
`q` controls their volume dependence, and `n` scales the number of
vibrational degrees of freedom per formula unit.

`MieGruneisenDebye.debye_temperature_law` selects the characteristic-temperature
relation. If it is omitted, the default is `integrated_gruneisen`:

\[
\Theta(V)=\Theta_0\exp\left\{
-\frac{\gamma_0}{q}
\left[\left(\frac{V}{V_0}\right)^q-1\right]
\right\}, \qquad q\ne0,
\]

with continuous limit

\[
\Theta(V)=\Theta_0\left(\frac{V}{V_0}\right)^{-\gamma_0},
\qquad q=0.
\]

For either oscillator model,

\[
\Delta P_{\mathrm{th}}
=10^{-4}\frac{\gamma(V)}{V}
\left[E(V,T)-E(V,T_r)\right].
\]

The third-order Debye function and Debye energy are

\[
D_3(y)=\frac{3}{y^3}\int_0^y\frac{z^3}{e^z-1}\,dz,
\qquad
E_D(V,T)=3nRTD_3\!\left(\frac{\Theta(V)}{T}\right).
\]

The Einstein energy is

\[
E_E(V,T)=\frac{3nR\Theta(V)}
{\exp[\Theta(V)/T]-1}.
\]

Zero-point energy is omitted from these two public vibrational-energy models;
it cancels from the referenced thermal pressure in any case. Constant-volume
heat capacity is evaluated as $C_V=(\partial E/\partial T)_V$.

### Variable-exponent Debye-temperature law

`MieGruneisenDebye(..., debye_temperature_law="variable_exponent")` retains
the Debye energy and referenced thermal-pressure expression above, but follows
Fei et al. (2007) literally for the characteristic temperature. With
$x=V/V_0$,

\[
\gamma(V)=\gamma_0x^q,
\qquad
\Theta_D(V)=\Theta_0x^{-\gamma(V)}.
\]

This differs from `integrated_gruneisen` when $q\ne0$. The law is explicit
fixed configuration, not a fitted numerical parameter, so reconstruction and
serialization preserve it while uncertainty propagation does not perturb it.
Fei et al. equation 3 defines the referenced thermal pressure; Table 1 supplies
the four catalog parameter sets.

### Tange 2009 MgO thermal model

`Tange2009Debye` retains the Debye energy and referenced Mie-Gruneisen
thermal-pressure expression above, but replaces the power-law Gruneisen model
with Tange et al. (2009), equation 15. With $x=V/V_0$,

\[
\gamma(V)=\gamma_0\left\{1+a\left[x^b-1\right]\right\}.
\]

Integrating $\gamma=-\partial\ln\Theta/\partial\ln V$ gives the paper's
equation 16 in an explicit form:

\[
\Theta(V)=\Theta_0\exp\left[-\gamma_0\left{
(1-a)\ln x+\frac{a}{b}(x^b-1)
\right\}\right], \qquad b\ne0.
\]

The implementation also evaluates the continuous $b=0$ limit. This model is
kept separate from the generic `MieGruneisenDebye` class so that the published
functional form is not approximated by a different $q$ law. The pressure
standard wrapper converts diffraction volumes from conventional-cell
$\mathring{\mathrm A}^3$ to the molar volume required by thermal EOS classes.

### Linear thermal pressure

`LinearThermalPressure` composes any reference isotherm with

\[
P(V,T)=P_{\mathrm{ref}}(V)+\alpha K_T(T-T_r).
\]

Here `alpha_KT` is the fitted product in GPa/K, not a separately constant
thermal expansivity and bulk modulus. Dewaele et al. (2012), equation 2 and
Table V, use this form for B2 KCl and KBr. Walker et al. (2002), equation BE1,
uses the same additive term for KCl; its reported B2 product
`0.0275(9) kbar/K` is represented directly as `0.00275(9) GPa/K` so the
published product uncertainty can be propagated without inventing independent
errors for its correlated factors. Because the correction is independent of
volume, it uses the same volume convention as the composed reference EOS.

### Thermal modified Tait

`ThermalModifiedTait` combines a `ModifiedTait` reference EOS with a fixed
Einstein temperature `theta`. Define

\[
E_E(T)=\frac{3nR\Theta}{\exp(\Theta/T)-1},
\qquad
\phi=\frac{\alpha_0K_0}{C_V(T_r)}.
\]

Then

\[
\Delta P_{\mathrm{th}}(T)=
\phi\left[E_E(T)-E_E(T_r)\right].
\]

This pressure is independent of volume. Its implied Gruneisen parameter is

\[
\gamma(V)=10^4V\phi.
\]

`HollandPowell2011` is an alias for this implementation.

### Multi-oscillator Gruneisen thermal pressure

`MultiOscillatorGruneisenThermalEOS` combines an independently selected
reference isotherm, volume-dependent oscillator temperatures, and quadratic
anharmonic and electronic pressure terms. It needs pressure, bulk modulus, and
$dK/dP$ from the reference component. `EosBase` supplies the last quantity by
central numerical differentiation when the selected isotherm has no specialized
implementation. The published Sokolova et al. catalog records compose the
thermal correction with Holzapfel; another reference isotherm is a valid API
composition, but a new scientific model rather than the published scale.

#### Paper versus spreadsheet

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

## Numerical inversion

`volume(P)` and `volume(P, T)` solve the applicable pressure equation with a
bracketed scalar root for every broadcast point. `temperature(P, V)` solves

\[
\Delta P_{\mathrm{th}}(V,T)=P-P_{\mathrm{ref}}(V)
\]

on the positive-temperature branch nearest $T_r$. The distinct two-volume
DAC inversion and its assumptions are documented under
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md).
