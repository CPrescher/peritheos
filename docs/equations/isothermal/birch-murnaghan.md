# Birch-Murnaghan family

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/isothermal/">Isothermal equations</a>
<span aria-hidden="true">›</span>
<span>Birch-Murnaghan family</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

## Finite-pressure BM2 reference

`BM2Reference(V0, K0, P0)` explicitly represents

\[
P(V)=P_0+\frac32 K_0[(V_0/V)^{7/3}-(V_0/V)^{5/3}].
\]

V0 and K0 denote volume and bulk modulus **at P0**, not ambient values unless
P0=0. The reference pressure is required, finite, and uses the same units as
K0. P(V0)=P0, K(V0)=K0 and dK/dP at V0 is 4. Bulk modulus is the usual BM2
derivative, since the offset is independent of V. The molar reference-energy
integral includes `-P0*(V-V0)*1e4` in the library's J/bar/mol volume convention.
Thermal reference-volume/modulus transformations retain P0.

This is a specified reference-pressure law, not an exact way of moving the
reference state of every existing ambient BM2 parameterization. Ordinary
`BM2(V0, K0)` is unchanged. Rust provides
`isothermal::BM2Reference::new(v0, k0, p0)` with scalar/batch evaluation,
inversion, and reference-energy integration.
Python's double-Debye constructors retain their explicit zero-pressure
Vinet/BM2/BM3/BM4 allowlist and reject this new model; no Python caloric
composition support is implied by the isothermal interface.

EOSMAT uses `model: "BM2Reference"` (canonical v3 identifier
`birch_murnaghan_2_reference`), with V0, K0 and P0 parameters. Older readers
without this model must reject the identifier; exporting v3 does not remove
that capability requirement. The source's original parameter names should
also be documented when they differ, such as V32.7 and K32.7 for
[Bykova's joint coesite-IV/V curve](../../literature-reproductions/coesite-iv-v-bykova-2018-experimental.md).

```python
from peritheos.eos.rt import BM2Reference

eos = BM2Reference(V0=375.04, K0=254.0, P0=32.7)
assert eos.pressure(375.04) == 32.7
assert eos.volume(32.7) == 375.04
```

## Standard zero-pressure forms

Define the Eulerian finite strain

\[
f_E=\frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right].
\]

The second-order model (`BM2`) is

\[
P=3K_0 f_E(1+2f_E)^{5/2},
\]

which is the third-order model with the constraint $K_0'=4$. The general
third-order model (`BM3`) is

\[
P=3K_0 f_E(1+2f_E)^{5/2}
\left[1+\frac{3}{2}(K_0'-4)f_E\right].
\]

For the fourth-order model (`BM4`), define

\[
\zeta=\frac{3}{4}(4-K_0'),
\qquad
\xi=\frac{3}{8}\left[
K_0K_0''+(K_0'-4)(K_0'-3)+\frac{35}{9}
\right].
\]

Then

\[
P=3K_0 f_E(1+2f_E)^{5/2}
\left(1-2\zeta f_E+4\xi f_E^2\right).
\]

Higher order does not automatically imply better extrapolation. In particular,
$K_0''$ is often strongly correlated with the lower-order parameters unless
the pressure range constrains curvature well.

---

**End of Birch-Murnaghan family.**

[All isothermal models](../../models/isothermal.md) · [Equation index](../../equation-reference.md)
