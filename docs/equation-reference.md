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

This reference-isotherm form also applies to the double-Debye Helmholtz models
when their optional `Tr` is set. With `Tr` omitted, those models instead retain
their original absolute simulated free-energy formulation.

Because energy divided by the public molar-volume unit produces bar, the factor
$10^{-4}$ converts thermal pressure to GPa.

### Double-Debye Helmholtz

`DoubleDebyeHelmholtz` supports two explicit reference-state conventions. With
`Tr=None` (the default), it is the complete simulated free-energy EOS and its
`rt_eos` is a Vinet cold curve for a motionless lattice at 0 K. With a numeric
`Tr`, `rt_eos` is instead the complete reference isotherm at `Tr` and the
simulated non-cold free energy is rebased to vanish there. Define

\[
x=(V/V_0)^{1/3},\qquad X=\frac32(K_0'-1)(x-1).
\]

The molar cold energy is

\[
E_{\rm cold}=\phi_0+
\frac{4V_0K_0\,10^4}{(K_0'-1)^2}
\left[1-(1+X)e^{-X}\right],
\]

where $10^4$ converts GPa times J bar$^{-1}$ mol$^{-1}$ to J mol$^{-1}$.
The implementation evaluates the continuous $K_0'\to1$ limit without the
apparent singularity.

Each of $\theta_A$, $\theta_B$, and $\theta_1$ has its own parameter triple:

\[
\theta_i(V)=\theta_{i0}(V/V_p)^{-b_i}
\exp[a_i(V_p-V)],\qquad \gamma_i=a_iV+b_i.
\]

Let $F_D(\theta,T)$ and $U_D(\theta,T)$ be the one-mole-of-atoms Debye
Helmholtz energy and internal energy, including zero point:

\[
F_D=R\left[\frac98\theta+3T\ln(1-e^{-\theta/T})
-TD_3(\theta/T)\right],
\]

\[
U_D=R\left[\frac98\theta+3T D_3(\theta/T)\right].
\]

The volume-dependent weights and ionic free energy are

\[
w_A=\frac{\theta_B-\theta_1}{\theta_B-\theta_A},\qquad
w_B=1-w_A,
\]

\[
F_{\rm ion}=n(w_AF_D(\theta_A,T)+w_BF_D(\theta_B,T)).
\]

Coincident characteristic temperatures are evaluated by their analytic
limiting weights, avoiding a $0/0$ at a shared reference point. The optional
$T^2$ contribution is

\[
\alpha(V)=\alpha_0(V/V_e)^\kappa,\qquad
F_{\rm anh}=-\frac12nR\alpha(V)T^2.
\]

For the absolute (`Tr=None`) convention, the complete energy returned by
`helmholtz_free_energy()` is

\[
F=E_{\rm cold}+F_{\rm ion}+F_{\rm anh}.
\]

Thermodynamic differentiation gives

\[
P_{\rm ion}=\frac{n}{10^4}\left[
\frac{w_A\gamma_AU_A+w_B\gamma_BU_B}{V}
-\frac{dw_A}{dV}(F_A-F_B)\right],
\]

\[
P_{\rm anh}=\frac{nR\kappa\alpha(V)T^2}{2V\,10^4},\qquad
P=P_{\rm Vinet}+P_{\rm ion}+P_{\rm anh}.
\]

For reference-isotherm anchoring at a numeric $T_r$, Peritheos evaluates

\[
F(V,T)=F_{\rm ref}(V)+F_{\rm ion}(V,T)-F_{\rm ion}(V,T_r)
+F_{\rm anh}(V,T)-F_{\rm anh}(V,T_r),
\]

and therefore

\[
P(V,T)=P_{\rm ref}(V)+P_{\rm ion}(V,T)-P_{\rm ion}(V,T_r)
+P_{\rm anh}(V,T)-P_{\rm anh}(V,T_r).
\]

Thus the supplied Vinet curve is recovered exactly at $T=T_r$. This operation
does not replace the simulated thermal model: it preserves its volume-dependent
temperature increment while removing its simulated reference-isotherm bias.
In particular, `Tr=0` is not equivalent to leaving `Tr` unset, because rebasing
at 0 K also subtracts zero-point motion.

The $dw_A/dV$ term is essential when the double-Debye weights vary with
volume. For `Tr=None`, `thermal_pressure()` returns
$P_{\rm ion}+P_{\rm anh}$, including zero-point pressure, while
`thermal_pressure_increment()` returns the pressure above the conventional
300 K comparison isotherm. With numeric `Tr`, both methods return the
reference-relative contribution, which is exactly zero at `Tr`. This lets
`volume_with_dac_confinement()` and `temperature_from_volumes()` use the same
empirical confinement equation as the other thermal models. Use the ordinary
`volume(P,T)` or `temperature(P,V)` whenever total hot pressure is known
independently.

#### Benedict et al. diamond example

The following is the diamond column of Benedict et al. (2014), Table I,
converted from per-atom to Peritheos molar units. The coefficients remain an
example rather than class defaults:

```python
from scipy.constants import Avogadro, electron_volt

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import DoubleDebyeHelmholtz

v_atom = Avogadro * 1e-25  # A^3/atom -> J/bar/mol
e_atom = electron_volt * Avogadro  # eV/atom -> J/mol

diamond = DoubleDebyeHelmholtz(
    rt_eos=Vinet(5.7034 * v_atom, 432.4, 3.793),
    Vp=5.571 * v_atom,
    theta_a0=1887.8,
    a_a=-0.316 / v_atom,
    b_a=0.913,
    theta_b0=1887.8,
    a_b=0.168 / v_atom,
    b_b=0.429,
    theta_1_0=1887.8,
    a_1=0.0846 / v_atom,
    b_1=0.499,
    n=1,
    alpha0=3.79e-5,
    Ve=5.785 * v_atom,
    kappa=0,
    phi0=-9.066 * e_atom,
)

pressure = diamond.pressure(4.6542704116 * v_atom, 3000.0)  # about 150 GPa
volume = diamond.volume(150.0, 3000.0)
temperature = diamond.temperature(150.0, volume)
```

To anchor those thermal coefficients to an independently measured 298 K Vinet
curve, pass that curve as `rt_eos` and add `Tr=298.0`. The bundled
`DIAMOND_BENEDICT_2014_DEWAELE_ANCHORED` and
`DIAMOND_CORREA_2008_DEWAELE_ANCHORED` records do this with Dewaele et al.
(2008). Their unanchored counterparts preserve the original simulation models.

For diamond, the paper identifies the $T^2$ coefficient as an anharmonic ionic
correction rather than an electronic excitation. Because its published
$\kappa=0$, that term changes free energy and heat capacity but not pressure.
The Vinet parameters describe the classical 0 K PBE-DFT lattice; they are not
a 300 K static-compression fit. The paper compares the diamond model directly
with DFT-MD over roughly 3.0--5.6 $\mathrm{\AA^3/atom}$ and 2000--9000 K.
Use outside the diamond stability field, and especially extrapolation of the
quadratic term above its low-temperature purpose, requires an independent
phase-stability and high-temperature plausibility assessment.

There is a factor-of-two inconsistency between the two carbon publications.
Correa et al. (2008), equation 18, use $F_{\rm anh}=-aT^2$ and tabulate
$a=3.8\times10^{-5}\ \mathrm{K^{-1}}$ for diamond. Benedict et al. (2014)
state that they retain that correction, but write
$F_{\rm anh}=-\alpha T^2/2$ and tabulate the essentially unchanged value
$\alpha=3.79\times10^{-5}\ \mathrm{K^{-1}}$. The conventions are equivalent
only when $\alpha=2a$. Peritheos preserves each publication literally rather
than silently choosing one interpretation; this affects caloric properties and
free-energy differences, but not diamond pressure because both coefficients
are volume independent.

#### Logarithmic-moment double-Debye variant

`DoubleDebyeLogMomentHelmholtz` implements the earlier Correa et al. (2008)
diamond branch. It shares the Vinet energy, Debye functions, and optional
reference-isotherm convention above, but conserves
the logarithmic phonon moment $\theta_0$ rather than the arithmetic moment
$\theta_1$:

\[
w_A=\frac{\ln(\theta_B/\theta_0)}{\ln(\theta_B/\theta_A)},\qquad
w_B=1-w_A.
\]

The three characteristic temperatures use independent parameter triples for
$A$, $B$, and $0$. Correa's volume-independent anharmonic term is

\[
F_{\rm anh}=-nR aT^2,
\]

so $P_{\rm anh}=0$ and its heat-capacity contribution is $2nRaT$. The bundled
`DIAMOND_CORREA_2008` record preserves the paper's DFT-GGA value
$V_0=5.785\ \mathrm{\AA^3/atom}$. The authors note that this is about 3% too
large after zero-point and thermal effects; Peritheos does not silently apply
their suggested application-dependent density shift.

The ordinary thermal fitting API can fit pressure-sensitive coefficients while
holding the 0 K Vinet curve fixed. `phi0` is an additive energy zero and cannot
be inferred from pressure-only observations; fix it unless the objective also
contains absolute free-energy data. Likewise, when `kappa=0`, `alpha0` does not
contribute to pressure and must not be treated as pressure-identifiable.

### Temperature-dependent reference state

`ThermalReferenceStateEOS` composes an independently selected reference
isotherm that exposes reconstructable `V0` and `K0` parameters. At each
temperature it evaluates that same isotherm after applying

\[
V_0(T)=V_0(T_r)\exp\left[\int_{T_r}^{T}\alpha(T')\,dT'\right],
\qquad
K_0(T)=K_0(T_r)+(T-T_r)\left(\frac{\partial K}{\partial T}\right)_P.
\]

The default `thermal_expansion_law="constant"` uses
$\alpha(T)=\alpha_0$. The `linear_temperature` law uses

\[
\alpha(T)=\alpha_0+\alpha_1T,
\qquad
V_0(T)=V_0(T_r)\exp\left[
\alpha_0(T-T_r)+\frac{\alpha_1}{2}(T^2-T_r^2)
\right].
\]

The `linear_reference_temperature` law instead defines `alpha0` at the
reference temperature:

\[
\alpha(T)=\alpha_0+\alpha_1(T-T_r),
\qquad
V_0(T)=V_0(T_r)\exp\left[
\alpha_0(T-T_r)+\frac{\alpha_1}{2}(T-T_r)^2
\right].
\]

This is the convention used by the bundled Suzuki (2016) epsilon-FeOOH
thermal EOS. The two linear laws are intentionally distinct because their
stored `alpha0` parameters have different reference points.

Here $\alpha$ is the volumetric expansion coefficient; crystallographic-axis
expansivities are separate quantities. The stored numeric parameters are
`Tr`, `alpha0`, `dK_dT`, and `alpha1` (zero by default). Unlike the energy-based
thermal models, this construction inherits the volume unit of its reference
isotherm.

The separate `reference_volume_law="linear_temperature"` option applies

\[
V_0(T)=V_0(T_r)[1+\alpha_0(T-T_r)].
\]

In this configuration, $\alpha_0$ is the mean expansion coefficient in the
direct reference-volume relation, not a constant instantaneous expansivity.
The bundled Martinez aragonite BM2 record uses this form with linear `K0(T)`,
following equations (2)--(3) and Table 7 of Martinez, Zhang, and Reeder (1996),
[doi:10.2138/am-1996-5-608](https://doi.org/10.2138/am-1996-5-608).

The `reference_volume_law="berman"` option applies

\[
V_0(T)=V_0(T_r)\left[1+\alpha_0(T-T_r)
+\frac{\alpha_1}{2}(T-T_r)^2\right].
\]

This is the truncated Berman (1988) thermal-expansion equation implemented by
EosFit7; `alpha0` is defined at `Tr`. The B4C record from Somayazulu et al.
(2023) uses this form. It must not be replaced by the exponential integral,
which EosFit identifies as the Fei model.

The bundled ice VI and ice VII records use the integrated constant law with BM2,
following equations (1)--(3) and Table II of Bezacier et al. (2014),
[doi:10.1063/1.4894421](https://doi.org/10.1063/1.4894421). The `.eosmat`
interchange type remains `AlphaKT` for Dioptas compatibility, while the stable
model identifier is `thermal_reference_state`. The integrated
linear-expansivity law follows
Martinez, Zhang, and Reeder (1996), equations (2), (4), and (5),
[doi:10.2138/am-1996-5-608](https://doi.org/10.2138/am-1996-5-608).

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

By default, either oscillator model uses the reference-temperature baseline

\[
\Delta P_{\mathrm{th}}
=10^{-4}\frac{\gamma(V)}{V}
\left[E(V,T)-E(V,T_r)\right].
\]

`MieGruneisenDebye` also accepts
`thermal_pressure_reference="absolute_zero"` for publications whose supplied
mechanical curve is explicitly a 0 K cold curve. In that mode,

\[
P_{\mathrm{th}}(V,T)=10^{-4}\frac{\gamma(V)}{V}E_D(V,T),
\qquad P_{\mathrm{th}}(V,0)=0,
\]

and total pressure is $P_c(V)+P_{\mathrm{th}}(V,T)$. `Tr` remains a positive
API reference for temperature inversion and for
`thermal_pressure_increment(V,T)`, which subtracts the absolute pressure at
`Tr`. Datchi et al. (2007), Equations (2)--(4) and Table V, use this convention
for c-BN. The default `reference_temperature` behavior is unchanged.

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
it cancels from referenced thermal pressure and is excluded from the
absolute-zero convention's source definition. Constant-volume
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

### Second-order temperature-compression thermal pressure

`SecondOrderTaylorThermalPressure` composes a cold reference curve exposing
`V0` with an absolute thermal pressure polynomial. Define

\[
\eta=1-\frac{V}{V_0},\qquad \Delta\eta=\eta-\eta_0,
\qquad \Delta T=T-T_r.
\]

Then

\[
P(V,T)=P_c(V)+c_0+c_1\Delta\eta+c_2\Delta T
+\frac{c_3}{2}\Delta\eta^2+\frac{c_4}{2}\Delta T^2
+\frac{c_5}{2}\Delta\eta\Delta T.
\]

The parameter units are GPa for `c0`, `c1`, and `c3`; GPa/K for `c2` and
`c5`; and GPa/K$^2$ for `c4`. Because the thermal term depends only on the
dimensionless ratio $V/V_0$, it inherits the reference EOS volume convention.
Unlike the usual referenced thermal corrections, its absolute thermal pressure
does not vanish at `Tr`. `thermal_pressure_increment(V,T)` remains available
and subtracts the value at `Tr` when a referenced increment is specifically
needed.

Luo et al. (2023), Appendix B, use this form with their zero-temperature
Rydberg--Vinet MgO curve. Treating that curve as a 300 K Vinet isotherm would
discard the nonzero `c0` and compression terms and therefore misrepresent the
published pressure scale.

### Logarithmic-volume linear thermal pressure

`LogVolumeThermalPressure` composes a reference isotherm exposing `V0` with

\[
P(V,T)=P_{\mathrm{ref}}(V)+
\left[\alpha K_{T,r}+
\left(\frac{\partial K_T}{\partial T}\right)_V
\ln\left(\frac{V_0}{V}\right)\right](T-T_r).
\]

The stored parameters are `Tr`, `alpha_KT_ref`, and `dK_dT_V`, in K and
GPa/K. This is Anderson, Isaak, and Yamamoto (1989), Equations (26)--(29),
[doi:10.1063/1.342969](https://doi.org/10.1063/1.342969). Unlike
`ThermalReferenceStateEOS`, it does not shift `V0` or `K0`; it adds a thermal
pressure whose temperature slope changes logarithmically with compression.
The generic class name describes that mechanism rather than the paper or its
use as a pressure standard.

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

### Dorogokupets--Oganov 2007 four-oscillator model

`DorogokupetsOganov2007` is a separate implementation of equations (7)--(14)
of [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/PhysRevB.75.024115).
It must not be substituted with the later Sokolova multi-oscillator model,
whose Gruneisen construction and computational reference differ. The supplied
`rt_eos` is the complete reference isotherm at `Tr` (298.15 K for the bundled
Pt scale).

With $x=V/V_0$, the mode Gruneisen parameter and characteristic temperatures
are

\[
\gamma(x)=\gamma_\infty+(\gamma_0-\gamma_\infty)x^\beta,
\]

\[
\Theta_i(x)=\Theta_{i0}x^{-\gamma_\infty}
\exp\left[\frac{\gamma_0-\gamma_\infty}{\beta}
(1-x^\beta)\right].
\]

This follows directly from integrating
$\gamma=-\partial\ln\Theta/\partial\ln V$. Two generalized-Bose modes use

\[
\varepsilon_B(\Theta,T,d)=\frac{\Theta(d-1)}{2d}
+\frac{T\Theta d}
{(Td+\Theta)[(1+\Theta/(Td))^d-1]},
\]

and two Einstein modes use

\[
\varepsilon_E(\Theta,T)=\Theta\left(\frac12+
\frac{1}{\exp(\Theta/T)-1}\right).
\]

The quasiharmonic pressure contribution is

\[
P_{\mathrm{qh}}=10^{-4}\frac{\gamma(x)}{V}
\sum_i m_iR\varepsilon_i.
\]

For all four modes define $n_i=[\exp(\Theta_i/T)-1]^{-1}$ and
$E_i=\Theta_i(1/2+n_i)$. The remaining molar Helmholtz terms are

\[
F_{\mathrm{anh}}=\sum_i m_iR\frac{a\,10^{-6}x^m}{6}
\left[E_i^2+2\Theta_i^2n_i(n_i+1)\right],
\]

\[
F_{\mathrm{el}}=-\frac32nRe\,10^{-6}x^gT^2,
\qquad
F_{\mathrm{def}}=-\frac32nRT
\exp\left(\frac{S}{x}-\frac{H}{Tx^2}\right).
\]

Peritheos evaluates their analytic volume derivatives and forms the nonreference
pressure $P_*(V,T)=P_{\mathrm{qh}}-10^{-4}\partial
(F_{\mathrm{anh}}+F_{\mathrm{el}}+F_{\mathrm{def}})/\partial V$. The combined
reference-isotherm EOS is

\[
P(V,T)=P_{\mathrm{ref}}(V)+P_*(V,T)-P_*(V,T_r).
\]

The bundled `platinum_dorogokupets_oganov_2007_vinet_4` record uses the Pt
column of Table I. Its published rounded coefficients reproduce all Pt
isochors printed in Table VI to within 0.02 GPa.

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
DAC inversion and forward `volume_with_dac_confinement(P_cold,T,f_dac)` solve
are documented under
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md).
