# Modified Tait

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Modified Tait</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

Let $v=V/V_0$ and define

\[
a=\frac{1+K_0'}{1+K_0'+K_0K_0''},
\]

\[
b=\frac{K_0'}{K_0}-\frac{K_0''}{1+K_0'},
\qquad
c=\frac{1+K_0'+K_0K_0''}
{K_0'^2+K_0'-K_0K_0''}.
\]

The modified Tait pressure and analytic inverse are

\[
P=\frac{1}{b}\left[
\left(\frac{v+a-1}{a}\right)^{-1/c}-1
\right],
\]

\[
V=V_0\left[1-a+a(1+bP)^{-c}\right].
\]

The real-valued implementation requires $(v+a-1)/a>0$ and rejects singular
coefficient combinations. When $K_0''=0$, the equation reduces to the
Murnaghan EOS.

---

**End of Modified Tait.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
