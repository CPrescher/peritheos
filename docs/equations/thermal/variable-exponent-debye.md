# Variable-exponent Debye-temperature law

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Variable-exponent Debye-temperature law</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`MieGruneisenDebye(..., debye_temperature_law="variable_exponent")` retains
the Debye energy and referenced thermal-pressure expression above, but follows
Fei et al. (2007) literally for the characteristic temperature. With
$x=V/V_0$,

\[
\gamma(V)=\gamma_0x^q,
\qquad
\Theta_D(V)=\Theta_0x^{-\gamma(V)}.
\]

Differentiating this expression gives

\[
-\frac{d\ln\Theta_D}{d\ln V}
=\gamma(V)\left[1+q\ln x\right].
\]

Thus it differs from `integrated_gruneisen` when $q\ne0$: the Grüneisen
parameter used in pressure is generally not the logarithmic derivative of
this $\Theta_D$. They coincide at $x=1$ or for $q=0$. The pressure and caloric
functions reproduce the selected expressions, but the identity
$P_{\mathrm{vib}}=-10^{-4}(\partial F_D/\partial V)_T$ and the integrated-law
bulk-modulus equation above must not be assumed for this option. The law is
explicit fixed configuration, not a fitted numerical parameter, so
reconstruction and serialization preserve it while uncertainty propagation
does not perturb it.
Fei et al. equation 3 defines the referenced thermal pressure; Table 1 supplies
the four catalog parameter sets.

---

**End of Variable-exponent Debye-temperature law.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
