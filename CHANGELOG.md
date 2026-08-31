# Changelog

All notable changes to Peritheos are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed

- Marked the published parameter errors for 14 Al, Cu, W, Ni, Ag, diamond,
  alpha/omega-Ti, Si-V/Si-VII/Si-X, Re, corundum, and LiF records as 95%
  confidence half-widths. File-loaded uncertainty propagation now converts
  these intervals to normal-equivalent standard errors instead of treating
  them as one-standard-deviation errors.
- Corrected the imported Dioptas Fei et al. (2007) Au and Ne records to use
  `MieGruneisenDebye` with `debye_temperature_law="variable_exponent"` rather
  than the implicit `integrated_gruneisen` default. The `.eosmat` records
  preserve the original behavior and cite equation 3 in explicit
  migration-correction metadata.
- Corrected the Hazen--Finger (1979) zircon record from an inconsistent BM2
  representation to BM3 with the published assumed `K0' = 6.5`, including the
  reported `V0` uncertainty.
- Corrected the Holmes et al. (1989) platinum record from BM3 to its published
  universal (Vinet) isotherm, restored its model reference volume and
  0--550 GPa static range, and represented Equation (12) with the published
  constant thermal-pressure coefficient.
- Restored the Ross (1997) magnesite fitted `V0` and uncertainty and normalized
  the Haines et al. (2001) Mo2C reference volume from the primary specimen's
  measured ambient lattice parameters.
- Corrected primary-source values and error metadata for CaSiO3, CaO B1/B2,
  rutile GeO2/SnO2, PbS B1, wadsleyite, jadeite, and B2 KCl; corrected the
  migrated B2-KCl thermal component to Walker et al.'s additive BE1 form and
  retained the published uncertainty of its directly fitted `alpha0*K0`
  product. Published errors are retained even where the associated value was
  fixed during a fit, as for Shim et al.'s CaSiO3 `V0`.

### Added

- A reproducible primary-source audit for all 147 migrated `.eosmat` EOS
  records. One hundred sixteen records are now directly validated against original
  publications, official supplements, or stable institutional reports; 31 are explicitly deferred with
  record-level reasons, and no record remains in a generic pending state. The
  bundled machine-readable ledger records source locations, the APS embargo
  affecting Shen--Smith (2026), the B4C order inconsistency, and the restored
  Hanfland graphite `V0` uncertainty.
- A mechanism-oriented `ThermalReferenceStateEOS` implementation for the
  temperature-dependent `V0(T)`/`K0(T)` formulation used by the validated ice
  VI/VII records. The Dioptas `AlphaKT` interchange type maps to the canonical
  `thermal_reference_state` model identifier.
- Primary-audit corrections restore Sokolova `n`/`Z`, silica Debye `n`, and
  ice `Tr` inputs omitted by migration; every validated migrated record is
  constructability-tested.
- A reproducible BurnMan/Pytheos public-API black-box comparison report,
  deliberately separated from primary-source validation and test baselines.
- Executable `Material` conversion through the same canonical `.eosmat`
  format 3 used for Dioptas exchange. Optional symmetry, lattice, space-group,
  atom-site, peak, and unknown extension fields survive a Peritheos round
  trip; cell-to-molar volume conversion is explicit per EOS record. Loading
  uses a fixed model registry, refuses unaudited records by default, and keeps
  snapshot-v2 reading only for compatibility.
- A configurable `debye_temperature_law` on `MieGruneisenDebye`, with
  `integrated_gruneisen` as the backward-compatible default and
  `variable_exponent` for sources that directly publish a volume-dependent
  exponent. Also added the mechanism-named
  `MultiOscillatorGruneisenThermalEOS` class, which accepts any
  isothermal `EosBase` and uses a generic numerical $dK/dP$ fallback where
  needed; the earlier paper-named imports remain compatibility aliases.

- A first-class `Material`/`EOSRecord` catalog API with GPa pressure, conventional-unit-cell
  volumes, scalar/array pressure and volume inversion, explicit material/phase
  and unit metadata, DOI-level parameter provenance, published validity
  envelopes, JSON-safe catalog records, and uncertainty propagation from
  measured volume/temperature and published parameter errors.
- A primary-source-validated catalog: Tange et al. (2009) Fit3-Vinet
  P-V-T MgO B1 and the Dorfman et al. (2012) 300 K Vinet co-compression scales
  for Au, Pt, Mo, NaCl B2, and Ne; Dewaele 2019 LiF and NaCl B1/B2; Dewaele
  2012 KCl and KBr B1/B2; Datchi 2007 c-BN; and Dewaele 2008 diamond, Ag,
  and Ni.
- All eleven Sokolova et al. (2016) thermal pressure markers: MgO, diamond,
  Al, Cu, Ag, Au, Pt, Nb, Ta, Mo, and W, with Table 1 provenance and the
  corrected 2016 equations.
- The Fei et al. (2007) internally consistent Au, Pt, NaCl-B2, and Ne thermal
  scales and a dedicated Debye-temperature convention that preserves the
  paper's equation rather than substituting the generic integrated form.
- The quasi-hydrostatic 300 K hcp Re Vinet scale of Anzellini et al. (2014),
  with Table III lattice-data regressions and its published 95% fit intervals
  retained distinctly from one-standard-deviation errors.
- The Tange et al. volume-dependent Gruneisen Mie-Gruneisen-Debye thermal model,
  with printed Table 5 regression cases and analytic thermodynamic checks.
- A reusable linear thermal-pressure EOS for the Dewaele KCl/KBr equation and
  state-only uncertainty propagation where a source reports no parameter errors.
- A documented catalog inventory, a versioned JSON computation record aligned
  with Dioptas's material-oriented `.eosmat` format, and explicit deferral
  records for Re and other entries where official primary evidence is not yet
  independently available.
- A Peritheos-owned flat `.eosmat` format 3, normative JSON Schema, complete
  116-material/147-record EOS database migrated from Dioptas 0.10.0 with explicit
  validation status and provenance, legacy Dioptas format-2 input, and tested
  Dioptas 0.10.0 read compatibility. A dedicated schema reference documents
  every field, discriminator pairing, default, unit, validation status, and
  consumer compatibility rule.

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
