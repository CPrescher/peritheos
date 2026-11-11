# Vinet

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Vinet</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

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

---

**End of Vinet.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
