# Mie-Gruneisen Debye and Einstein

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Mie-Gruneisen Debye and Einstein</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

The Mie–Grüneisen–Debye (MGD) model combines a mechanical reference EOS with
a vibrational thermal-pressure contribution. The Debye model describes the
thermal population of lattice vibrations; the Grüneisen parameter describes
how their frequencies change with volume. In the default reference-isotherm
convention, the total pressure is

\[
P(V,T)=P_{\mathrm{ref}}(V,T_r)+\Delta P_{\mathrm{th}}(V,T).
\]

Here $P_{\mathrm{ref}}$ is the curve supplied as `rt_eos`, for example BM3 or
Vinet. Its parameters describe the isotherm at $T_r$, including the thermal
pressure already present there. Subtracting the vibrational contribution at
that same temperature prevents counting it twice.

## Symbols, parameters, and units

| Symbol | Python parameter or method | Meaning and units |
|---|---|---|
| $V$, $V_0$ | `V`, `rt_eos.V0` | Molar volume and reference volume, J bar$^{-1}$ mol$^{-1}$ |
| $T$, $T_r$ | `T`, `Tr` | State and reference temperatures, K |
| $\Theta(V)$, $\Theta_0$ | `characteristic_temperature(V)`, `theta0` | Debye temperature at $V$ and $V_0$, K |
| $\gamma(V)$, $\gamma_0$ | `gruneisen_parameter(V)`, `gamma0` | Grüneisen parameter at $V$ and $V_0$, dimensionless |
| $q$ | `q` | Constant logarithmic volume derivative of $\gamma$, dimensionless |
| $n$ | `n` | Atoms per formula unit defining the molar basis |
| $C_\infty$ | `Cvmax` | High-temperature molar $C_V$ limit; defaults to $3nR$, J mol$^{-1}$ K$^{-1}$ |
| $R$ | — | Molar gas constant, approximately 8.31446 J mol$^{-1}$ K$^{-1}$ |
| $E_D$, $S_D$, $F_D$ | `thermal_energy`, `thermal_entropy`, `thermal_helmholtz_free_energy` | Vibrational energy, entropy, and Helmholtz energy on the same molar basis |

Throughout this section, let $x=V/V_0$ and $y=\Theta(V)/T$. The Debye
temperature is a characteristic phonon energy scale, not the sample
temperature. For fixed $V$, it is independent of $T$ in this quasi-harmonic
model. Volumes, `Tr`, `theta0`, `n`, and an explicit `Cvmax` must be positive;
`gamma0` and `q` must be finite. State temperatures must be strictly positive
in the public API; expressions at 0 K below denote limiting values.

## Grüneisen parameter and Debye temperature

The default volume law is

\[
\gamma(V)=\gamma_0x^q,
\qquad
q=\frac{d\ln\gamma}{d\ln V},
\qquad
\gamma(V)=-\frac{d\ln\Theta}{d\ln V}.
\]

The logarithmic definition of $q$ assumes positive $\gamma$. For the usual
$\gamma_0>0$, positive $q$ makes $\gamma$ decrease on compression ($x<1$),
while positive $\gamma$ makes the Debye temperature increase. Integrating the
last relation from $V_0$ to $V$ gives

\[
\ln\frac{\Theta(V)}{\Theta_0}
=-\int_1^x\gamma_0 u^q\,d\ln u
=-\frac{\gamma_0}{q}(x^q-1),\qquad q\ne0,
\]

and therefore

\[
\boxed{\Theta(V)=\Theta_0
\exp\!\left[-\frac{\gamma_0}{q}(x^q-1)\right].}
\]

The continuous $q=0$ limit has constant $\gamma=\gamma_0$:

\[
\Theta(V)=\Theta_0x^{-\gamma_0}.
\]

This is `debye_temperature_law="integrated_gruneisen"`, the default for
`MieGruneisenDebye` and the law used by `MieGruneisenEinstein`. It ensures
$\Theta(V_0)=\Theta_0$ and $\gamma(V_0)=\gamma_0$. The alternative
[`variable_exponent` law](variable-exponent-debye.md) has a
different volume derivative and must be selected explicitly.

## Debye energy and heat capacity

The third-order Debye function is

\[
D_3(y)=\frac{3}{y^3}\int_0^y\frac{z^3}{e^z-1}\,dz.
\]

Here $z$ is a dimensionless phonon energy divided by $k_BT$, and $y$ is the
Debye cutoff divided by $T$. The factor $1/(e^z-1)$ is the Bose occupation
factor. The Debye density of vibrational states is proportional to frequency
squared, giving the $z^3$ integrand after multiplication by the phonon energy.
The normalization gives $3n$ modes per formula unit when $C_\infty=3nR$.

The molar vibrational energy, excluding zero-point energy, is

\[
E_D(V,T)=C_\infty T D_3(y).
\]

With the default $C_\infty=3nR$, this is equivalently

\[
E_D(V,T)=9nRT\left(\frac{T}{\Theta(V)}\right)^3
\int_0^{\Theta(V)/T}\frac{z^3}{e^z-1}\,dz.
\]

At fixed volume, $\Theta(V)$ is constant. Using

\[
\frac{dD_3}{dy}=\frac{3}{e^y-1}-\frac{3D_3(y)}{y}
\]

gives the constant-volume heat capacity

\[
\boxed{C_V(V,T)=\left(\frac{\partial E_D}{\partial T}\right)_V
=C_\infty\left[4D_3(y)-\frac{3y}{e^y-1}\right].}
\]

`molar_heat_capacity_v()` evaluates the energy derivative numerically; this
analytic expression explains the quantity it returns. An explicit `Cvmax`
replaces $3nR$ throughout energy, entropy, and heat capacity. It is a fitted
normalization when a source uses one, and its molar basis must match $V$.

Useful limiting checks at fixed volume are

\[
\begin{aligned}
T\gg\Theta:&\quad D_3(y)\to1,\qquad
E_D\sim C_\infty T,\qquad C_V\to C_\infty,\\
T\ll\Theta:&\quad D_3(y)\sim\frac{\pi^4}{5y^3},\qquad
E_D\sim\frac{\pi^4 C_\infty}{5}\frac{T^4}{\Theta^3},\qquad
C_V\sim\frac{4\pi^4 C_\infty}{5}\left(\frac{T}{\Theta}\right)^3.
\end{aligned}
\]

The low-temperature $T^3$ heat capacity distinguishes the Debye model from a
single-frequency Einstein oscillator. Numerically, Peritheos uses a small-$y$
series, quadrature at intermediate $y$, and the large-$y$ asymptote for $D_3$.

## Thermal pressure and its reference state

The unreferenced vibrational pressure is

\[
P_{\mathrm{vib}}(V,T)=10^{-4}\frac{\gamma(V)}{V}E_D(V,T).
\]

For the integrated Grüneisen law, this also follows from
$P_{\mathrm{vib}}=-10^{-4}(\partial F_D/\partial V)_T$. The factor $10^{-4}$
converts bar to GPa: $E_D/V$ is in bar when energy is in J mol$^{-1}$ and
volume in J bar$^{-1}$ mol$^{-1}$. If an equation instead uses volume in
cm$^3$ mol$^{-1}$, its numerical prefactor is $10^{-3}$. Convert that volume
before passing it to Peritheos; see [Units and reference states](../../units.md).

For `thermal_pressure_reference="reference_temperature"` (the default),
both energies are evaluated at the **same volume** and the same $\Theta(V)$:

\[
\boxed{\displaystyle
\Delta P_{\mathrm{th}}(V,T)
=10^{-4}\frac{C_\infty\gamma(V)}{V}
\left[
T D_3\!\left(\frac{\Theta(V)}{T}\right)
-T_r D_3\!\left(\frac{\Theta(V)}{T_r}\right)
\right].}
\]

Thus $\Delta P_{\mathrm{th}}(V,T_r)=0$ at every volume and
$P(V,T_r)=P_{\mathrm{ref}}(V,T_r)$. For positive $\gamma$, heating above $T_r$
gives positive thermal pressure and cooling below it gives negative thermal
pressure. If **both** $T$ and $T_r$ greatly exceed $\Theta(V)$,

\[
\Delta P_{\mathrm{th}}\simeq
10^{-4}\frac{C_\infty\gamma(V)}{V}(T-T_r).
\]

A high sample temperature alone does not justify applying this approximation
to a low-temperature reference energy.

`MieGruneisenDebye` also supports two other explicitly chosen baselines:

| `thermal_pressure_reference` | Required meaning of `rt_eos` | Energy subtracted from $E_D(V,T)$ |
|---|---|---|
| `reference_temperature` | Isotherm at $T_r$ | $E_D(V,T_r)$ |
| `reference_isentrope` | Reference isentrope through $(V_0,T_r)$ | $E_D(V,T_S(V))$ |
| `absolute_zero` | 0 K cold curve in the source's zero-point convention | None |

For a source whose mechanical term represents a reference isentrope (for
example a BM3 isentrope), the Debye entropy is constant when $\Theta/T$ is
constant. The reference path therefore obeys

\[
T_S(V)=T_r\frac{\Theta(V)}{\Theta_0},\qquad
\Delta P_{\mathrm{th}}=10^{-4}\frac{\gamma(V)}{V}
\left[E_D(V,T)-E_D(V,T_S(V))\right].
\]

The correction vanishes on $T=T_S(V)$, rather than on the entire $T=T_r$
isotherm. For the integrated law, $d\ln T_S/d\ln V=-\gamma(V)$.

For `absolute_zero`, total pressure is

\[
P(V,T)=P_c(V)+10^{-4}\frac{\gamma(V)}{V}E_D(V,T),
\qquad \lim_{T\to0^+}P_{\mathrm{th}}(V,T)=0.
\]

`Tr` remains a positive API reference for temperature inversion. In **all**
three modes, `thermal_pressure_increment(V,T)` returns
$P(V,T)-P(V,T_r)$ at fixed volume. In the latter two modes this differs from
`thermal_pressure(V,T)`, which is the contribution added to `rt_eos`.
Datchi et al. (2007), Equations (2)--(4) and Table V, use the absolute-zero
convention for c-BN. Select the baseline from the source's definition of its
mechanical curve; changing the option does not convert reference coefficients.

## Thermoelastic derivatives

For any of these temperature-independent reference curves,

\[
\left(\frac{\partial P}{\partial T}\right)_V
=10^{-4}\frac{\gamma(V)C_V(V,T)}{V},
\qquad
\alpha(V,T)=10^{-4}\frac{\gamma(V)C_V(V,T)}{VK_T(V,T)}.
\]

Here $\alpha$ is volumetric expansivity and
$K_T=-V(\partial P/\partial V)_T$ is the bulk modulus of the **total** pressure
surface. It generally differs from the reference-curve modulus. For the
default integrated law with a fixed reference temperature, define
$\Delta E=E_D(V,T)-E_D(V,T_r)$ and $C_{V,r}=C_V(V,T_r)$. Then

\[
K_T(V,T)=K_{\mathrm{ref}}(V)+K_{\mathrm{th}}(V,T),
\]

\[
K_{\mathrm{th}}=(1-q)\Delta P_{\mathrm{th}}
+10^{-4}\frac{\gamma(V)^2}{V}
\left[\Delta E-TC_V(V,T)+T_rC_{V,r}\right].
\]

This follows by differentiating both $\gamma/V$ and the volume-dependent
Debye energies; holding $\Theta$ constant during that derivative would miss
the second term. At $T=T_r$, $K_{\mathrm{th}}=0$. Peritheos evaluates the
mechanical derivatives numerically. Heat-capacity and adiabatic identities,
including $C_P-C_V$ and $K_S=K_TC_P/C_V$, are described under
[Thermoelastic properties](../../thermoelastic-properties.md).

## Entropy, free energy, and normalization

The same Debye oscillator gives

\[
S_D(V,T)=C_\infty\left[\frac{4}{3}D_3(y)
-\ln(1-e^{-y})\right],
\]

\[
F_D(V,T)=E_D-TS_D
=C_\infty T\left[\ln(1-e^{-y})-\frac{1}{3}D_3(y)\right].
\]

These are unreferenced vibrational contributions. `thermal_energy()` and
`thermal_internal_energy()` return $E_D$, not the difference
$E_D(V,T)-E_D(V,T_r)$ used by referenced pressure. Likewise,
`thermal_helmholtz_free_energy()` returns $F_D$, not a total potential for the
mechanical reference curve. For the default convention, the thermal part of
a potential that generates the referenced pressure is $F_D(V,T)-F_D(V,T_r)$.
The [thermoelastic documentation](../../thermoelastic-properties.md#vibrational-thermodynamic-contributions)
explains the returned enthalpy and Gibbs contributions.

Zero-point energy is omitted from these public Debye and Einstein energies.
It cancels from fixed-volume referenced pressure and is excluded from the
absolute-zero mode's thermal contribution. Static formation energies are also
absent; these quantities alone do not define absolute chemical potentials.
The model assumes a quasi-harmonic lattice within a single phase; electronic,
explicitly anharmonic, magnetic, and phase-transition terms require additional
models or source-specific treatment.

Energy and volume must use the same molar basis. `n` is the number of atoms
in that formula unit, not the crystallographic number of formula units per
cell, $Z$. For elemental hcp Fe, use `n=1` with volume per mole of Fe atoms;
$Z=2$ converts the two-atom cell volume to that molar basis. An equivalent
per-mole-cell model uses `n=2` and doubles both $V$ and $V_0$ (and an explicit
`Cvmax`). This preserves $x$, $\gamma(V)$, $\Theta(V)$, and $C_\infty/V$.
Changing only `n` doubles thermal pressure at fixed other coefficients, while
leaving reference-curve pressure unchanged. See
[Atom count and molar-volume basis](../../units.md#atom-count-and-molar-volume-basis).

## Worked MGD calculation

This illustrative parameter set is not a calibrated material or pressure
standard. Its reference molar volume is 10 cm$^3$/mol, entered as
1 J bar$^{-1}$ mol$^{-1}$:

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye

eos = MieGruneisenDebye(
    rt_eos=BM3(V0=1.0, K0=160.0, K0_prime=4.0),
    Tr=300.0,
    theta0=800.0,
    gamma0=1.5,
    q=1.0,
    n=2.0,
)
V, T = 0.9, 1200.0  # J/bar/mol, K
print(f"gamma = {eos.gruneisen_parameter(V):.4f}")
print(f"theta = {eos.characteristic_temperature(V):.2f} K")
print(f"reference pressure = {eos.rt_eos.pressure(V):.4f} GPa")
print(f"thermal pressure = {eos.thermal_pressure(V, T):.4f} GPa")
print(f"total pressure = {eos.pressure(V, T):.4f} GPa")
print(f"C_V = {eos.molar_heat_capacity_v(V, T):.4f} J/mol/K")
```

Expected output:

```text
gamma = 1.3500
theta = 929.47 K
reference pressure = 20.8162 GPa
thermal pressure = 6.0295 GPa
total pressure = 26.8458 GPa
C_V = 48.4218 J/mol/K
```

At this volume, $x=0.9$, $\gamma=1.35$, and
$\Theta=800\exp[1.5(1-0.9)]\simeq929.47$ K. Use this same Debye temperature
in both terms of the referenced pressure. At 300 K, the thermal correction
is exactly zero, although the vibrational energy itself remains positive.

## Einstein alternative

`MieGruneisenEinstein` uses the same integrated volume law and default
reference-temperature pressure subtraction, but replaces the Debye spectrum
by a single frequency:

\[
E_E(V,T)=\frac{3nR\Theta(V)}{e^{\Theta(V)/T}-1},\qquad
C_{V,E}=3nR\frac{y^2e^y}{(e^y-1)^2}.
\]

Its high-temperature heat-capacity limit is also $3nR$, but its
low-temperature heat capacity decays exponentially. `Cvmax` and the alternate
pressure-reference and Debye-temperature options are specific to
`MieGruneisenDebye`.

---

**End of Mie-Gruneisen Debye and Einstein.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
