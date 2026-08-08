# Changelog

All notable changes to Peritheos are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Branch-aware coverage enforcement, Ruff lint and formatting gates, and a
  minimum-supported-dependency CI job.
- Literature-tagged numerical cases for all public EOS families and an
  independently solvable weighted least-squares fitting benchmark.
- GitHub Pages documentation deployment and automatic GitHub Releases with
  source and wheel artifacts.
- Citation metadata, API stability, contribution, security, support, conduct,
  issue, and pull-request policies.
- Versioned Read the Docs builds using the locked uv documentation environment.

### Changed

- Uncertainty propagation tests now cover Monte Carlo state sampling, invalid
  covariance and option handling, one-sided numerical derivatives, wrapper
  behavior, and failed sampling; branch-aware project coverage exceeds 90%.
- Corrected the fourth-order Birch-Murnaghan documentation to match its
  implemented Eulerian-strain signs and exponents.
- Package maturity metadata now identifies the 0.4 development line as beta.

## [0.3.0] - 2026-08-08

### Added

- Murnaghan, modified Tait, and second- through fourth-order natural-strain
  room-temperature equations of state.
- Holland-Powell thermal modified Tait equation of state.
- Bounded P-V and P-V-T errors-in-variables fitting with independently optional
  pressure, volume, and temperature standard errors.
- `FitResult` covariance, correlation, parameter-error, residual, chi-square,
  information-criterion, and adjusted-state diagnostics.
- `EOSUncertainty` for propagating complete covariance matrices, parameter
  errors with correlations, or partial independent parameter errors into EOS
  calculations.
- Linear covariance propagation and reproducible Monte Carlo propagation for
  pressure, volume, bulk modulus, and arbitrary public EOS quantities.
- Optional propagation of pressure, volume, and temperature errors in the
  requested state.
- Explicit combination of separately quantified thermal and reference-EOS
  uncertainty blocks.
- Isothermal compressibility, thermal expansivity, constant-volume and
  constant-pressure heat capacities, thermodynamic Gruneisen parameter, and
  adiabatic bulk modulus APIs.
- Vibrational entropy, internal energy, Helmholtz energy, enthalpy, and Gibbs
  energy for the Mie-Gruneisen-Debye and Mie-Gruneisen-Einstein models.
- Molar-volume and density conversion helpers with explicit unit validation.
- Literature-tagged reference cases and expanded numerical regression tests.
- A multi-page MkDocs documentation site covering models, fitting,
  uncertainty, properties, units, validation, references, and development.

### Changed

- EOS models now expose reconstructable parameter values and safe parameter
  replacement, including dotted reference-EOS parameters such as
  `rt_eos.K0`.
- Thermal EOS calculations consistently document and enforce their molar
  volume, pressure, and temperature units.
- CI now builds the documentation strictly, and source distributions include
  documentation and literature-validation data.
- The legacy fitting argument `sigma` remains supported as an alias for
  `pressure_sigma`.

### Statistical assumptions

- Parameters omitted from an uncertainty specification are treated as exact.
- Individual parameter errors without a correlation matrix are treated as
  mutually independent.
- Combining thermal and reference-EOS results from separate fits requires an
  explicit block-independence assumption.

## [0.2.0] - 2026-08-07

### Added

- Fourth-order Birch-Murnaghan support.
- Complete Sokolova et al. (2016) pressure terms and parameters.
- Automated test and trusted PyPI publication workflows.

### Changed

- EOS pressure and volume validation and inversion were hardened for invalid
  and out-of-domain states.
- Project naming and release metadata were standardized.

[Unreleased]: https://github.com/CPrescher/peritheos/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/CPrescher/peritheos/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/CPrescher/peritheos/releases/tag/v0.2.0
