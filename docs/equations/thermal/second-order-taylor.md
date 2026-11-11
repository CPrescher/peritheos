# Second-order temperature-compression thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Second-order temperature-compression thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

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

---

**End of Second-order temperature-compression thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
