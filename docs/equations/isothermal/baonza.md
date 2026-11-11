# Baonza pseudospinodal

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Baonza pseudospinodal</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

The three-parameter Baonza model (`Baonza`) fixes the empirical exponent
$\beta=0.85$. Define

\[
Q=-P_{sp}=\frac{\beta K_0}{K'_0},\qquad
\kappa^*=\frac{Q^\beta}{K_0},
\]

and

\[
A=Q^{1-\beta}-\frac{1-\beta}{\kappa^*}\ln(V/V_0).
\]

The explicit pressure and bulk modulus are

\[
P=A^{1/(1-\beta)}-Q,\qquad
K_T=\frac{A^{\beta/(1-\beta)}}{\kappa^*}.
\]

The model is evaluated only where $A>0$. It exactly recovers $P(V_0)=0$,
$K(V_0)=K_0$, and $(\partial K/\partial P)_0=K'_0$. The implementation follows
the pseudospinodal form of Baonza et al., *Physical Review B* **53**, 5252-5258
(1996), <https://doi.org/10.1103/PhysRevB.53.5252>, as used for the experimental
benchmark fits of Sun et al. (2010).

---

**End of Baonza pseudospinodal.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
