# Temperature-dependent reference state

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Temperature-dependent reference state</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`ThermalReferenceStateEOS` composes an independently selected reference
isotherm that exposes reconstructable `V0` and `K0` parameters. At each
temperature it evaluates that same isotherm after applying

\[
V_0(T)=V_0(T_r)\exp\left[\int_{T_r}^{T}\alpha(T')\,dT'\right],
\qquad
K_0(T)=K_0(T_r)+(T-T_r)\left(\frac{\partial K}{\partial T}\right)_P.
\]

The default `thermal_expansion_law="constant"` uses
$\alpha(T)=\alpha_0$. The `linear_temperature` law uses

\[
\alpha(T)=\alpha_0+\alpha_1T,
\qquad
V_0(T)=V_0(T_r)\exp\left[
\alpha_0(T-T_r)+\frac{\alpha_1}{2}(T^2-T_r^2)
\right].
\]

The `linear_reference_temperature` law instead defines `alpha0` at the
reference temperature:

\[
\alpha(T)=\alpha_0+\alpha_1(T-T_r),
\qquad
V_0(T)=V_0(T_r)\exp\left[
\alpha_0(T-T_r)+\frac{\alpha_1}{2}(T-T_r)^2
\right].
\]

This is the convention used by the bundled Suzuki (2016) epsilon-FeOOH
thermal EOS. The two linear laws are intentionally distinct because their
stored `alpha0` parameters have different reference points.

Here $\alpha$ is the volumetric expansion coefficient; crystallographic-axis
expansivities are separate quantities. The stored numeric parameters are
`Tr`, `alpha0`, `dK_dT`, and `alpha1` (zero by default). Unlike the energy-based
thermal models, this construction inherits the volume unit of its reference
isotherm.

For references with an independently adjustable pressure derivative, the
optional `dKprime_dT` coefficient (K^-1, default zero) also applies

\[
K'_0(T)=K'_0(T_r)+\mathtt{dKprime\_dT}(T-T_r)
+\mathtt{kprime\_log\_coefficient}(T-T_r)\ln(T/T_r).
\]

The linear term represents the temperature-dependent BM3 formulation described
by Fukui et al. (2003), p. 513,
[doi:10.1007/s00269-003-0353-z](https://doi.org/10.1007/s00269-003-0353-z).
The optional logarithmic term retains the Hirose (2008) law. They add to each
other and both vanish at `Tr`; reference `K0_prime` always means its value at
`Tr`. These shifts support BM3, Baonza, Murnaghan, Morse3, NaturalStrain3,
SunMorse3, SunMorse4 and Vinet. Other reference families reject nonzero shifts.
No `K0_double_prime(T)` law is implied. Both slopes are finite, signed numeric
coefficients. The shifted reference EOS validates its own parameters at each
state, including positive `V0(T)` and `K0(T)`.

For example, the following executes a **conditional** brucite audit surface,
anchoring Fukui's reported temperature slopes to the printed 300 K slice.
The author's separately regressed line intercepts remain unverified; this
example does not establish a validated brucite phase range or reproduce the
author's exact workflow:

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS

eos = ThermalReferenceStateEOS(
    BM3(V0=40.917, K0=41.8, K0_prime=6.6),
    Tr=300.0,
    alpha0=6.69e-5,
    dK_dT=-0.0319,
    dKprime_dT=0.00571,
)
eos.pressure(34.0, 873.0)  # 14.0514545153 GPa
```

The model supports pressure, bulk modulus, state inversion, native fitting and
parameter-covariance propagation for this coefficient. It remains a mechanical
P-V-T model without a specified caloric free energy.

The separate `reference_volume_law="linear_temperature"` option applies

\[
V_0(T)=V_0(T_r)[1+\alpha_0(T-T_r)].
\]

In this configuration, $\alpha_0$ is the mean expansion coefficient in the
direct reference-volume relation, not a constant instantaneous expansivity.
The bundled Martinez aragonite BM2 record uses this form with linear `K0(T)`,
following equations (2)--(3) and Table 7 of Martinez, Zhang, and Reeder (1996),
[doi:10.2138/am-1996-5-608](https://doi.org/10.2138/am-1996-5-608).

The `reference_volume_law="berman"` option applies

\[
V_0(T)=V_0(T_r)\left[1+\alpha_0(T-T_r)
+\frac{\alpha_1}{2}(T-T_r)^2\right].
\]

This is the truncated Berman (1988) thermal-expansion equation implemented by
EosFit7; `alpha0` is defined at `Tr`. The B4C record from Somayazulu et al.
(2023) uses this form. It must not be replaced by the exponential integral,
which EosFit identifies as the Fei model.

The bundled ice VI and ice VII records use the integrated constant law with BM2,
following equations (1)--(3) and Table II of Bezacier et al. (2014),
[doi:10.1063/1.4894421](https://doi.org/10.1063/1.4894421). The `.eosmat`
interchange type remains `AlphaKT` for Dioptas compatibility, while the stable
model identifier is `thermal_reference_state`. The integrated
linear-expansivity law follows
Martinez, Zhang, and Reeder (1996), equations (2), (4), and (5),
[doi:10.2138/am-1996-5-608](https://doi.org/10.2138/am-1996-5-608).

---

**End of Temperature-dependent reference state.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
