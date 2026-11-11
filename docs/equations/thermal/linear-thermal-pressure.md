# Linear thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Linear thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

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

---

**End of Linear thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
