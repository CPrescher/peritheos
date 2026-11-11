# Holzapfel

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Holzapfel</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

The Holzapfel implementation requires $V_0$ and $V$ in
`J bar^-1 mol^-1`, $K_0$ in GPa, $n$ as the number of atoms in
the formula unit, and $Z$ as its atomic-number parameter. For a compound,
$Z$ may be the effective value specified by its source parameterization;
Sokolova et al. equation 3 gives that value for MgO. Define

\[
x=\left(\frac{V}{V_0}\right)^{1/3},
\qquad
P_{\mathrm{FG0}}=1003.6
\left(\frac{Zn}{10V_0}\right)^{5/3},
\]

\[
c_0=-\ln\left(\frac{3K_0}{P_{\mathrm{FG0}}}\right),
\qquad
c_2=\frac{3}{2}(K_0'-3)-c_0.
\]

The pressure is

\[
P=3K_0\exp[c_0(1-x)]
\left(x^{-5}-x^{-4}\right)
\left(1+c_2x-c_2x^2\right).
\]

`bulk_modulus_derivative(V)` evaluates $\partial K_T/\partial P$ numerically
and is used by `MultiOscillatorGruneisenThermalEOS`.

---

**End of Holzapfel.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
