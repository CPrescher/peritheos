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
- Executable tutorial notebooks provide complete workflows:
  - [Pressure calibration quickstart](notebooks/pressure-calibration-quickstart.ipynb)
    calculates, inverts, validates, and propagates uncertainty for a hot MgO
    pressure marker.
  - [Normalize mixed pressure scales](notebooks/normalize-mixed-pressure-scales.ipynb)
    finds a common target for Au-, Pt-, and ruby-referenced sample EOSs, then
    executes and audits each transformation route.
  - [Comparing room-temperature EOS](notebooks/compare-room-temperature-eos.ipynb)
    contrasts fit diagnostics, parameter correlation, and extrapolation for
    four common equations.
  - [Birch-Murnaghan f-F diagnostic](notebooks/birch-murnaghan-f-f-diagnostic.ipynb)
    compares BM2 and BM3 fits to bundled room-temperature coesite data and
    constructs the conventional normalized-stress plot explicitly.
  - [Comparing gold pressure scales](notebooks/compare-gold-pressure-scales.ipynb)
    evaluates named literature records only on their common validity domains.
  - [Thermal EOS state surfaces](notebooks/thermal-eos-state-surfaces.ipynb)
    uses broadcasting, inversion, and thermoelastic properties on a P-V-T grid.
  - [Fit to prediction uncertainty](notebooks/fit-to-prediction-uncertainty.ipynb)
    carries fitted covariance into correlated linear and Monte Carlo intervals.
  - [Exploring the material library](notebooks/exploring-material-library.ipynb)
    searches the complete bundled catalog and constructs executable records.
  - [EOSMAT round-trip](notebooks/eosmat-roundtrip.ipynb) validates, writes,
    reloads, and safely handles provisional material records.
  - [DAC temperature sensitivity](notebooks/dac-temperature-sensitivity.ipynb)
    studies the two-volume thermal-confinement boundary condition.
  - [Aragonite P-V-T fitting](notebooks/aragonite-eos-fitting.ipynb) reproduces
    a published 64-point staged thermal analysis.
- [Calculation uncertainty](uncertainty.md) propagates fitted or published
  parameter errors into EOS predictions.
- [Dioptas and `.eosmat`](dioptas-integration.md) defines the shared material
  ownership and bundled 139-material/211-record EOS library.
- [`.eosmat` schema reference](eosmat-schema.md) documents every exchange
  field, model discriminator, default, unit, and compatibility rule.
- [Advanced DAC analysis](dac-thermal-pressure.md) documents the optional
  two-volume thermal-pressure boundary-condition model.
- [Units and reference states](units.md) explains the deliberately strict
  thermal volume convention.
- [Validation](validation.md) describes the numerical and literature checks.
- [Primary EOS refit validation](primary-eos-refits.md) reports a fit attempt
  for every bundled material record and lists all coefficient-parity failures.
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

The site executes its bundled fitting notebook against the local CSV dataset;
it does not download scientific data during the build:

```bash
mkdocs build --strict
```

## Citation

Cite Peritheos using the repository's `CITATION.cff` metadata and cite the
original source for each EOS used. See [References](references.md) for the
model-specific literature.
