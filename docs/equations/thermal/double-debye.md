# Double-Debye Helmholtz

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Double-Debye Helmholtz</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`DoubleDebyeHelmholtz` supports two explicit reference-state conventions. With
`Tr=None` (the default), it is the complete simulated free-energy EOS and its
`rt_eos` is a `Vinet`, `BM2`, `BM3`, or `BM4` cold curve for a motionless
lattice at 0 K. With a numeric
`Tr`, `rt_eos` is instead the complete reference isotherm at `Tr` and the
simulated non-cold free energy is rebased to vanish there. For Vinet, define

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

For Birch--Murnaghan references, define the Eulerian strain
$f=[(V_0/V)^{2/3}-1]/2$. Integrating the corresponding pressure gives

\[
E_{\rm ref}=\phi_0+\frac92 V_0K_0\,10^4 f^2
\left[1+c_3 f+c_4 f^2\right],
\]

with $c_3=c_4=0$ for BM2; $c_3=K_0'-4$, $c_4=0$ for BM3; and

\[
c_3=K_0'-4,\qquad
c_4=\frac34\left[K_0K_0''+(K_0'-4)(K_0'-3)+\frac{35}{9}\right]
\]

for BM4. In every case $E_{\rm ref}(V_0)=\phi_0$ and
$P_{\rm ref}=-10^{-4}\,dE_{\rm ref}/dV$. At finite reference temperature,
this energy primitive represents the reference Helmholtz energy difference.
The oscillator terms below are independent of the chosen reference family.

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
P=P_{\rm cold}+P_{\rm ion}+P_{\rm anh}.
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

Thus the supplied reference curve is recovered exactly at $T=T_r$. This operation
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

## Benedict et al. diamond example

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

## Logarithmic-moment double-Debye variant

`DoubleDebyeLogMomentHelmholtz` implements the earlier Correa et al. (2008)
diamond branch. It shares the reference-curve energy, Debye functions, and optional
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

---

**End of Double-Debye Helmholtz.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
