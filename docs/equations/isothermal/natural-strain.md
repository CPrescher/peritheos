# Natural-strain family

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Natural-strain family</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

The positive compressive natural strain is

\[
f_N=\frac{1}{3}\ln\left(\frac{V_0}{V}\right).
\]

All three natural-strain models use

\[
P=3K_0\frac{V_0}{V}
\left(f_N+A f_N^2+B f_N^3\right).
\]

Their coefficients are

\[
\begin{array}{c|cc}
\text{model} & A & B \\
\hline
\text{NaturalStrain2} & 0 & 0 \\
\text{NaturalStrain3} & \frac{3}{2}(K_0'-2) & 0 \\
\text{NaturalStrain4} & \frac{3}{2}(K_0'-2) &
\frac{3}{2}\left[K_0K_0''+1+(K_0'-2)+(K_0'-2)^2\right]
\end{array}
\]

`NaturalStrain2` therefore has the implied reference derivative $K_0'=2$.

---

**End of Natural-strain family.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
