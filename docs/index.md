# Peritheos

Peritheos is a focused Python library for pressure-volume and
pressure-volume-temperature equations of state for solid materials. It provides
forward evaluation, pressure-volume inversion, thermoelastic properties, and
least-squares parameter fitting.

The documentation is organized around scientific workflows:

- [Getting started](getting-started.md) introduces isothermal and thermal use.
- [EOS models](models.md) lists parameters, assumptions, and useful limits.
- [Equation reference](equation-reference.md) defines the implemented
  isothermal and thermal equations and their coefficients.
- [Thermoelastic properties](thermoelastic-properties.md) defines every derived
  property and its units.
- [Fitting data](fitting.md) covers P-V and P-V-T regression with uncertainty.
- [Calculation uncertainty](uncertainty.md) propagates fitted or published
  parameter errors into EOS predictions.
- [Dioptas and `.eosmat`](dioptas-integration.md) defines the shared material
  ownership and bundled 115-material EOS library.
- [`.eosmat` schema reference](eosmat-schema.md) documents every exchange
  field, model discriminator, default, unit, and compatibility rule.
- [Advanced DAC analysis](dac-thermal-pressure.md) documents the optional
  two-volume thermal-pressure boundary-condition model.
- [Units and reference states](units.md) explains the deliberately strict
  thermal volume convention.
- [Validation](validation.md) describes the numerical and literature checks.
- [API reference](api.md) is a compact import and method reference.
- [API stability](api-stability.md) defines the compatibility contract.

## Scope

Peritheos models homogeneous solid phases on a chosen reference isotherm. It
does not currently calculate phase equilibria, solid solutions, shear moduli,
or fluid fugacities. Thermodynamic free energies exposed by the Mie-Gruneisen
models are vibrational contributions; they are not absolute formation energies.

## Building these docs

```bash
python -m pip install -r docs/requirements.txt
mkdocs serve
```

The site can be built without executing notebook code or downloading data:

```bash
mkdocs build --strict
```

## Citation

Cite Peritheos using the repository's `CITATION.cff` metadata and cite the
original source for each EOS used. See [References](references.md) for the
model-specific literature.
