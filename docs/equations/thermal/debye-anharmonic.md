# Monatomic Debye and anharmonic Helmholtz energy

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Monatomic Debye and anharmonic Helmholtz energy</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`DebyeAnharmonicHelmholtz` implements Xiao et al. (2025),
[doi:10.1007/s10765-024-03469-2](https://doi.org/10.1007/s10765-024-03469-2),
equations (19)–(25). With $x=V/V_0$, $t=T/\theta_0$:

\[
\gamma=\gamma_0x^q,\quad
\theta=\theta_0\exp[-\gamma_0(x^q-1)/q],\quad
A_D=RT\{3\ln[1-e^{-\theta/T}]-D_3(\theta/T)\},
\]

\[
A_a=b_1R\theta_0\frac{t^4}{1+b_2t^2}\exp[b_3(x-1)],\qquad
P_{\rm th}=\frac{\gamma U_D}{10^4V}-\frac{b_3A_a}{10^4V_0}.
\]

The continuous $q=0$ limit is $\theta=\theta_0x^{-\gamma_0}$.
The Debye term excludes zero-point energy and adds its full pressure to a
zero-temperature cold curve. `Tr` only defines pressure increments and an
inversion starting temperature; pressure does not subtract thermal energy at
`Tr`. Positive finite volume and temperature are required, as are positive
`Tr` and `theta0`; `b2` must be nonnegative to avoid a positive-temperature
pole. All other coefficients must be finite. $R=8.31451$ J/mol/K matches the
authors' supplementary workbook, rather than the modern gas constant used
by other Peritheos Debye families.

For the published cold curve, $L=\ln(V_0/V)$ and
$A_0=10^4 V_0(c_1L^2/2+c_2L^3/3+c_3L^4/4)$ when $c_i$ are in GPa and
$V_0$ in J/bar/mol. `NaturalStrain4` reproduces this exactly with
$K_0=c_1$, $K'_0=2+2c_2/c_1$ and
$K''_0=(6c_3/c_1-(K'_0)^2+3K'_0-3)/K_0$.

`helmholtz_free_energy`, `internal_energy`, and `entropy` expose full molar
properties with cold energy zero at $V_0$. `thermal_helmholtz_free_energy` excludes
the cold energy. These functions do not include the source's gas-phase
reference-state offset for phase equilibrium. `molar_heat_capacity_v` includes
the analytic anharmonic derivative; other thermal and caloric methods follow
the common API. Volumes are molar J/bar/mol, pressure is GPa, energy is J/mol,
and heat capacity and entropy are J/mol/K.

---

**End of Monatomic Debye and anharmonic Helmholtz energy.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
