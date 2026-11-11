# Tange 2009 MgO thermal model

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Tange 2009 MgO thermal model</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`Tange2009Debye` retains the Debye energy and referenced Mie-Gruneisen
thermal-pressure expression above, but replaces the power-law Gruneisen model
with Tange et al. (2009), equation 15. With $x=V/V_0$,

\[
\gamma(V)=\gamma_0\left\{1+a\left[x^b-1\right]\right\}.
\]

Integrating $\gamma=-\partial\ln\Theta/\partial\ln V$ gives the paper's
equation 16 in an explicit form:

\[
\Theta(V)=\Theta_0\exp\left[-\gamma_0\left\{
(1-a)\ln x+\frac{a}{b}(x^b-1)
\right\}\right], \qquad b\ne0.
\]

The implementation also evaluates the continuous $b=0$ limit. This model is
kept separate from the generic `MieGruneisenDebye` class so that the published
functional form is not approximated by a different $q$ law. The pressure
standard wrapper converts diffraction volumes from conventional-cell
$\mathring{\mathrm A}^3$ to the molar volume required by thermal EOS classes.

---

**End of Tange 2009 MgO thermal model.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
