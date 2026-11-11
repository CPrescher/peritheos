# Dorogokupets--Oganov 2007 four-oscillator model

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Dorogokupets--Oganov 2007 four-oscillator model</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`DorogokupetsOganov2007` is a separate implementation of equations (7)--(14)
of [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/PhysRevB.75.024115).
It must not be substituted with the later Sokolova multi-oscillator model,
whose Gruneisen construction and computational reference differ. The supplied
`rt_eos` is the complete reference isotherm at `Tr` (298.15 K for the bundled
Pt scale).

With $x=V/V_0$, the mode Gruneisen parameter and characteristic temperatures
are

\[
\gamma(x)=\gamma_\infty+(\gamma_0-\gamma_\infty)x^\beta,
\]

\[
\Theta_i(x)=\Theta_{i0}x^{-\gamma_\infty}
\exp\left[\frac{\gamma_0-\gamma_\infty}{\beta}
(1-x^\beta)\right].
\]

This follows directly from integrating
$\gamma=-\partial\ln\Theta/\partial\ln V$. Two generalized-Bose modes use

\[
\varepsilon_B(\Theta,T,d)=\frac{\Theta(d-1)}{2d}
+\frac{T\Theta d}
{(Td+\Theta)[(1+\Theta/(Td))^d-1]},
\]

and two Einstein modes use

\[
\varepsilon_E(\Theta,T)=\Theta\left(\frac12+
\frac{1}{\exp(\Theta/T)-1}\right).
\]

The quasiharmonic pressure contribution is

\[
P_{\mathrm{qh}}=10^{-4}\frac{\gamma(x)}{V}
\sum_i m_iR\varepsilon_i.
\]

For all four modes define $n_i=[\exp(\Theta_i/T)-1]^{-1}$ and
$E_i=\Theta_i(1/2+n_i)$. The remaining molar Helmholtz terms are

\[
F_{\mathrm{anh}}=\sum_i m_iR\frac{a\,10^{-6}x^m}{6}
\left[E_i^2+2\Theta_i^2n_i(n_i+1)\right],
\]

\[
F_{\mathrm{el}}=-\frac32nRe\,10^{-6}x^gT^2,
\qquad
F_{\mathrm{def}}=-\frac32nRT
\exp\left(\frac{S}{x}-\frac{H}{Tx^2}\right).
\]

Peritheos evaluates their analytic volume derivatives and forms the nonreference
pressure $P_*(V,T)=P_{\mathrm{qh}}-10^{-4}\partial
(F_{\mathrm{anh}}+F_{\mathrm{el}}+F_{\mathrm{def}})/\partial V$. The combined
reference-isotherm EOS is

\[
P(V,T)=P_{\mathrm{ref}}(V)+P_*(V,T)-P_*(V,T_r).
\]

The bundled `platinum_dorogokupets_oganov_2007_vinet_4` record uses the Pt
column of Table I. Its published rounded coefficients reproduce all Pt
isochors printed in Table VI to within 0.02 GPa.

---

**End of Dorogokupets--Oganov 2007 four-oscillator model.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
