# EOS models

Start with a model category. Each equation family has its own page, including
its formulas, parameter conventions, limits, and source qualifications.

## Isothermal models

[Browse isothermal models](models/isothermal.md) for reference-isotherm equations,
including the Birch–Murnaghan, natural-strain, Vinet, and modified-Tait families.
The inventory lists their constructor parameters and links to each family.

## Thermal models

[Browse thermal models](models/thermal.md) for reference-EOS requirements,
thermal parameters, and caloric capabilities. Follow a model name to its equations
and reference-state conventions.

## Shock Hugoniot models

[`LinearUsUpHugoniot`](hugoniots.md) implements a phase-specific constrained
shock path from `Us = c0 + s * up`. Typed `HugoniotRecord` instances remain in
the material's common `eos_records` collection, marked by
`equation_kind="hugoniot"`, but are deliberately not accepted as reference
isotherms for thermal EOS composition. Hugoniot and equilibrium models share
the `EquationOfState` base interface; the canonical namespaced import is
`peritheos.eos.hugoniot` and `peritheos.hugoniot` remains compatible.

## Selecting a model

- Use BM3 or Vinet for a conventional three-parameter compression fit.
- Use modified Tait or a fourth-order natural-strain EOS when reliable second
  pressure-derivative information exists.
- Use Mie-Gruneisen-Debye when low-temperature heat-capacity behavior matters.
  Its `debye_temperature_law` defaults to `integrated_gruneisen`; select
  `variable_exponent` only for scales that explicitly publish
  $\Theta_D=\Theta_0(V/V_0)^{-\gamma(V)}$.
  `thermal_pressure_reference` defaults to `reference_temperature`; use
  `reference_isentrope` when the wrapped curve is a reference isentrope, and
  `absolute_zero` only when the wrapped curve is explicitly a 0 K cold curve.
- Use Mie-Gruneisen-Einstein for an inexpensive single-frequency approximation.
- Use `Tange2009Debye` for the publication-specific Tange MgO Gruneisen law.
- Use `SoundVelocityDebyeHelmholtz` when a source derives the Debye temperature
  from a cold-curve bulk modulus and published longitudinal/shear
  velocity-density regressions. Those four regression coefficients are part of
  the model definition and must not be silently substituted.
- Use `Dewaele2006` for the source-specific hcp-Fe pressure scale used by
  Dewaele et al. (2006), including its intrinsic-anharmonic and electronic
  pressure terms.
- Use `DorogokupetsOganov2007` only for the paper's complete four-oscillator
  Helmholtz formulation. It is distinct from the later Sokolova workbook
  model even though the authors and physical terms overlap.
- Use linear thermal pressure only when the primary calibration reports a
  constant fitted $\alpha K_T$ product over the intended range.
- Use `SecondOrderTaylorThermalPressure` only for a source that publishes the
  corresponding absolute second-order pressure polynomial, such as Luo et al.'s
  MgO scale; do not reinterpret its cold curve as a reference-temperature isotherm.
- Use thermal modified Tait with Holland-Powell-style datasets.
- Use `MultiOscillatorGruneisenThermalEOS` for the multimode formulation; pass
  `n` explicitly and select the reference isotherm independently. Catalog
  entries pin the exact source-validated combinations.
- Use `DoubleDebyeHelmholtz` when the source supplies volume-dependent
  double-Debye and $T^2$ free-energy terms. Leave `Tr` unset when `rt_eos` is
  the source's motionless-ion 0 K cold curve. Set `Tr` when `rt_eos` is an
  independently measured reference isotherm and only the simulated thermal
  change away from that isotherm should be added.
- Use `DoubleDebyeLogMomentHelmholtz` when the source constrains the two Debye
  weights with the logarithmic phonon moment $\theta_0$, as in Correa et al.
  (2008). It is intentionally distinct from the arithmetic-$\theta_1$
  constraint in `DoubleDebyeHelmholtz`.

Do not choose a higher-order model solely because it has more parameters. Check
parameter correlations and extrapolation behavior with the fitting diagnostics.

`Sokolova2016` remains a compatibility alias. New code and serialized model
identifiers use neutral mechanism-based names; paper names belong to catalog
record identifiers and provenance.
