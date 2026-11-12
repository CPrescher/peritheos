# Q-compromise Debye thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Q-compromise Debye thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`QCompromiseDebye(rt_eos, Tr, theta0, gamma0, n)` implements the constant-θ,
constant-γ/V approximation of [Angel et al. (2022)](https://doi.org/10.1007/s00410-022-01918-5):

\[
P(V,T)=P_{ref}(V)+\frac{\gamma_0}{10^4 V_0}
\left[E_D(\theta_0,T)-E_D(\theta_0,T_r)\right],
\qquad E_D(\theta,T)=3nRTD_3(\theta/T).
\]

Internal V0 is J/bar/mol of formula units, `n` is atoms per formula unit,
`theta0` and `Tr` are positive finite kelvin, and γ0 must be finite. State
volumes and temperatures must be positive and finite. There is no q parameter.
The correction vanishes at Tr and is independent of V; consequently the
isothermal bulk modulus at a fixed V equals the reference modulus.
This pressure approximation supports native Python/Rust evaluation, fitting,
uncertainty and inversion. It does not define caloric potentials, heat
capacities or a compositional mixing law. See the [garnet source audit](../../literature-reproductions/garnet-angel-2022.md)
and [complete installed-package example](../../examples/garnet-q-compromise.py).

---

**End of Q-compromise Debye thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
