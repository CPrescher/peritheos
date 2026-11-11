# Asymptotic-power-law Debye with T2 excess

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Asymptotic-power-law Debye with T2 excess</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`AsymptoticPowerLawMieGruneisenDebyeExcess` retains the Tange Gruneisen and
Debye-temperature laws and adds the Zhu et al. (2025) excess Helmholtz term

\[
F_{\mathrm{ex}}(V,T)=-\frac{1}{2}\beta_0
\left(\frac{V}{V_0}\right)^mT^2.
\]

Its pressure contribution is

\[
P_{\mathrm{ex}}(V,T)=\frac{\beta_0m}{2V_0}
\left(\frac{V}{V_0}\right)^{m-1}T^2.
\]

Both the Debye and excess pressures are referenced by subtracting their values
at `Tr`. `beta0` is expressed in J mol$^{-1}$ K$^{-2}$ and `m` is
dimensionless. This is the complete thermal pressure form used by the Zhu Au,
Pt, and MgO P-V-T records.

---

**End of Asymptotic-power-law Debye with T2 excess.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
