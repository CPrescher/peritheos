# Generalized Rydberg-Stacey

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Generalized Rydberg-Stacey</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

The generalized Rydberg-Stacey model (`RydbergStacey`) introduces the
adjustable infinite-pressure bulk-modulus derivative $K'_{\infty}$. Define

\[
x=\left(\frac{V}{V_0}\right)^{1/3},
\qquad
A=\frac{3}{2}K'_0-3K'_{\infty}+\frac{1}{2}.
\]

The pressure is

\[
P=3K_0x^{-3K'_{\infty}}(1-x)\exp[A(1-x)].
\]

The analytic bulk modulus is

\[
K_T=K_0x^{-3K'_{\infty}}\exp[A(1-x)]
\left[x+3K'_{\infty}(1-x)+Ax(1-x)\right].
\]

At $K'_{\infty}=2/3$ this reduces exactly to the Vinet equation. Values such
as $K'_{\infty}=5/3$ instead impose the thermodynamic high-pressure limit used
by Stacey-style pressure scales.

---

**End of Generalized Rydberg-Stacey.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
