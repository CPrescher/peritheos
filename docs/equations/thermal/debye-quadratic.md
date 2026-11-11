# Debye plus quadratic thermal pressure

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Debye plus quadratic thermal pressure</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`DebyeQuadraticThermalPressure` adds the integrated-law Mie–Grüneisen–Debye
pressure relative to `Tr` and `A (V/V0)^m (T² − Tr²)`. `A` uses GPa/K²,
`m` is dimensionless, and `Tr`, `theta0`, `gamma0`, `q`, `n` follow the Debye
model. Internal volumes are molar J/bar/mol. Positive finite volumes,
temperatures, `Tr`, `theta0`, and `n` are required; other coefficients must
be finite. `A = 0` recovers MGD. The thermal increment vanishes at `Tr`.
It is a pressure surface without caloric observables. The
[Fei (2016) audit](../../literature-reproductions/fei-2016-iron.md) derives the
mapping from the published mass-specific electronic/anharmonic term and
explains the independent electronic Grüneisen convention.

---

**End of Debye plus quadratic thermal pressure.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
