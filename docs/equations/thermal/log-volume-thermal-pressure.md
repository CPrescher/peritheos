# Logarithmic-volume linear thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Logarithmic-volume linear thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`LogVolumeThermalPressure` composes a reference isotherm exposing `V0` with

\[
P(V,T)=P_{\mathrm{ref}}(V)+
\left[\alpha K_{T,r}+
\left(\frac{\partial K_T}{\partial T}\right)_V
\ln\left(\frac{V_0}{V}\right)\right](T-T_r).
\]

The stored parameters are `Tr`, `alpha_KT_ref`, and `dK_dT_V`, in K and
GPa/K. This is Anderson, Isaak, and Yamamoto (1989), Equations (26)--(29),
[doi:10.1063/1.342969](https://doi.org/10.1063/1.342969). Unlike
`ThermalReferenceStateEOS`, it does not shift `V0` or `K0`; it adds a thermal
pressure whose temperature slope changes logarithmically with compression.
The generic class name describes that mechanism rather than the paper or its
use as a pressure standard.

---

**End of Logarithmic-volume linear thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
