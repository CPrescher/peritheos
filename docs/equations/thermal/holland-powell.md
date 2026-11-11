# Holland--Powell thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Holland--Powell thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`HollandPowellThermalPressure` combines any reference EOS exposing `K0` with a
fixed Einstein temperature `theta`. Define

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

`ThermalModifiedTait` is the compatibility subclass that requires a
`ModifiedTait` reference EOS. `HollandPowell2011` remains an alias for that
compatibility class. The generic class is needed when a source applies the same
Holland--Powell thermal-pressure term to another reference isotherm, as
Anzellini et al. (2025) do with BM3 for iridium.

---

**End of Holland--Powell thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
