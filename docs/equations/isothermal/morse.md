# Morse and Sun-Morse

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Morse and Sun-Morse</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

Sun et al. (2010) define $X=(V/V_0)^{1/n}$. Their ordinary three-dimensional
Morse form (`Morse3`) uses $n=3$, $\beta=K'_0-1$, and

\[
P=\frac{3K_0}{\beta X^2}
\left[\exp(2\beta(1-X))-\exp(\beta(1-X))\right].
\]

The volume-analytic Sun Jiu-Xun--Morse forms (`SunMorse3` and `SunMorse4`)
instead use $n=3$ or $4$, $\alpha=(nK'_0+1)/3$, and

\[
P=\frac{nK_0}{\alpha}
\left[\exp(2\alpha(1-X))-\exp(\alpha(1-X))\right].
\]

Their analytic bulk modulus is

\[
K_T=K_0X\left[2\exp(2\alpha(1-X))-\exp(\alpha(1-X))\right].
\]

These forms are implemented from equations 6--9 and 15 of Sun et al.,
*Zeitschrift für Naturforschung A* **65**, 34--44 (2010),
<https://doi.org/10.1515/zna-2010-1-202>.

---

**End of Morse and Sun-Morse.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
