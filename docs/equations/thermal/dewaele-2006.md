# Dewaele 2006 hcp-Fe thermal pressure scale

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Dewaele 2006 hcp-Fe thermal pressure scale</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`Dewaele2006` implements Equations (1)--(2) of
[Dewaele et al. (2006)](https://doi.org/10.1103/PhysRevLett.97.215504).
By default, `rt_eos` is the complete reference isotherm at $T_r$, and every
thermal term is subtracted at that temperature:

\[
P(V,T)=P_{r}(V)+P_{\mathrm{th}}(V,T)-P_{\mathrm{th}}(V,T_r).
\]

With $x=V/V_0$, the model uses

\[
\gamma(x)=\gamma_\infty+(\gamma_0-\gamma_\infty)x^\beta,
\qquad
\Theta(x)=\Theta_0x^{-\gamma_\infty}
\exp\!\left[\frac{\gamma_0-\gamma_\infty}{\beta}(1-x^\beta)\right].
\]

The absolute thermal pressure before reference subtraction is

\[
P_{\mathrm{th}}=
\frac{3nR\gamma}{V}T D_3(\Theta/T)
+\frac{3nR}{2V}m a_0 x^m T^2
+\frac{3nR}{2V}g e_0 x^g T^2.
\]

`anharmonic_a` and `electronic_e` are passed directly in K$^{-1}$, matching
the Dewaele equation. This differs from the micro-K$^{-1}$ parameter convention
used by `DorogokupetsOganov2007`.

With `thermal_pressure_reference="absolute_zero"`, `rt_eos` is instead a
0 K curve and the unreferenced thermal pressure is added without zero-point
pressure. `Tr` remains the positive-temperature anchor for
`thermal_pressure_increment`; it is not subtracted from total pressure.
This represents [Dewaele et al. (2021) fcc argon](../../literature-reproductions/argon-fcc.md)
exactly with gamma0=2.7, gamma_inf=0.5, beta=1 and zero quadratic terms.
The reference-isentrope option is unsupported for this model.

---

**End of Dewaele 2006 hcp-Fe thermal pressure scale.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
