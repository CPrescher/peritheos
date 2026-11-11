# Murnaghan

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Murnaghan</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

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

---

**End of Murnaghan.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
