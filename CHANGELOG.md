# Changelog

All notable changes to Peritheos are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- A Rust workspace containing native EOS, fitting, uncertainty, and private
  PyO3 binding crates, with Rust 1.83 as the library MSRV.
- Shared Python/Rust compatibility fixtures and migration baselines for all
  isothermal and thermal model families.
- Multi-platform native-wheel release jobs for supported CPython versions on
  Linux x86-64/ARM64, macOS Intel/Apple Silicon, and Windows x86-64.

### Changed

- Built-in Python EOS classes now preserve their public API while delegating
  evaluation, inversion, thermoelastic, and caloric calculations to Rust.
- Named bounded robust fitting losses and uncertainty propagation statistics
  now use native numerical kernels. Custom reference EOS classes, callable
  fitting losses, and NumPy-seeded Monte Carlo draws retain documented
  compatibility paths.
- Errors-in-variables fits now use colored latent-coordinate Jacobians and a
  block Schur-complement solve, with stress coverage for large and rank-deficient
  datasets.

## [0.5.0] - 2026-08-30

### Added

- A comprehensive equation reference for every public isothermal and thermal
  EOS, including parameter definitions, model domains, numerical inversion,
  thermoelastic identities, fitting diagnostics, uncertainty propagation, and
  explicit public-unit conversion factors.
- Documentation distinguishing the equations typeset in Sokolova et al. (2016)
  from the accompanying Excel workbook calculation implemented by
  `Sokolova2016`, including differences in the bulk-modulus derivative,
  characteristic-temperature multiplier, oscillator pressure contributions,
  reference-isotherm subtraction, and working pressure units.
- Temperature inversion from pressure and volume for all thermal EOS models,
  including NumPy broadcasting through `temperature()` and
  `calculate_temperature()`.
- Coupled temperature inference from ambient and heated volumes with a
  fractional DAC confinement contribution through `temperature_from_volumes()`.
  The empirical `f_dac * thermal_pressure` pressure increment is solved in its
  algebraically reduced form and requires `0 <= f_dac < 1`.
  Documentation distinguishes this fraction of EOS thermal pressure from a
  fraction of cold pressure and describes its physical limits, calibration,
  identifiability, baseline-drift, and uncertainty constraints.
- Fixed-volume preparation for temperature inversion, avoiding repeated
  Sokolova volume-integral evaluations during root finding.

### Changed

- Refreshed the locked runtime, development, documentation, and release
  dependencies to their latest versions compatible with the supported Python
  and declared package-version ranges.
- The minimum supported SciPy version is now 1.9.3, ensuring a stable binary
  installation and numerical fitting path on Python 3.9 across supported
  platforms.
- Documentation now uses Material for MkDocs with responsive navigation,
  improved search and code presentation, and automatic light and dark themes.
- The model overview now focuses on model selection and reference-EOS
  compatibility, while advanced two-volume DAC analysis has moved out of the
  introductory tutorial into a dedicated guide.

## [0.4.0] - 2026-08-09

### Added

- Branch-aware coverage enforcement, Ruff lint and formatting gates, and a
  minimum-supported-dependency CI job.
- Literature-tagged numerical cases for all public EOS families and an
  independently solvable weighted least-squares fitting benchmark.
- Automatic GitHub Releases with source and wheel artifacts.
- Citation metadata, API stability, contribution, security, support, conduct,
  issue, and pull-request policies.
- Versioned Read the Docs builds using the locked uv documentation environment.
- Joint reference-isotherm and thermal parameter fitting with complete
  cross-covariance through `fit_joint_eos`.
- Per-observation correlated P-V and P-V-T covariance matrices and robust
  least-squares losses for all fitting entry points.
- Human-readable fit summaries and versioned, JSON-safe result export including
  model parameters, covariance, adjusted observations, diagnostics, and solver
  metadata.

### Changed

- Uncertainty propagation tests now cover Monte Carlo state sampling, invalid
  covariance and option handling, one-sided numerical derivatives, wrapper
  behavior, and failed sampling; branch-aware project coverage exceeds 90%.
- Corrected the fourth-order Birch-Murnaghan documentation to match its
  implemented Eulerian-strain signs and exponents.
- Package maturity metadata now identifies the 0.4 development line as beta.
- The minimum supported NumPy version is now 1.21, matching the public typing
  APIs used by Peritheos.
- Distribution package discovery is restricted to `peritheos`, preventing
  generated documentation directories from entering or blocking builds.

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

[Unreleased]: https://github.com/CPrescher/peritheos/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/CPrescher/peritheos/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/CPrescher/peritheos/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/CPrescher/peritheos/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/CPrescher/peritheos/releases/tag/v0.2.0
