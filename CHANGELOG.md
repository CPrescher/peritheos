# Changelog

All notable changes to Peritheos are documented here. The project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Added the Luo et al. (2023) B1-MgO pressure scale as a 0 K Vinet curve plus
  the published absolute second-order thermal-pressure polynomial, with all
  five new shock states and the complete 576-value Tables II--III P-V-T grid.

- Added the Campbell and Heinz (1994) RbCl-B2 BM3 material, including all 24
  Table 1 observations and a parity-checked refit. Corrected the CsCl resource,
  which had mistakenly contained an incomplete subset of the RbCl table, and
  added all nine Yagi (1978) CsCl compression ratios needed for the complete
  22-point parity reproduction.

- Added Birch-Murnaghan Eulerian finite-strain/normalized-stress diagnostics,
  including first-order P-V error propagation and fitted-model curves. The
  numerical API is plotting-library independent; an executable notebook shows
  how to produce the conventional F-f plot with Matplotlib.

- Added phase-specific linear `Us`-`up` shock Hugoniot EOSs in Python and Rust,
  including Rankine--Hugoniot pressure-volume inversion, velocity, density,
  energy, tangent-modulus, uncertainty propagation, OLS/WLS/errors-in-variables
  fitting, JSON-serializable fit results, and typed precursor, mass-basis, and
  branch-domain APIs. Hugoniots remain in `Material.eos_records`, with filtered
  `hugoniot_records` and `equilibrium_records` views and category-scoped
  defaults. Loading history and transformed-branch identity are independent;
  record evaluation enforces the declared branch domain, while `V0`, `rho0`,
  formula units, and molar mass are cross-validated. Structured derivation
  metadata covers coefficients obtained from SESAME or published tables.
- Added a bundled executable pressure-calibration library for the Mao (1978),
  Mao--Xu--Bell (1986), Dewaele (2004), Holzapfel (2005), and
  Dorogokupets--Oganov (2007), and IPPS-Ruby2020 ruby R1 scales in Python and
  Rust, together with the Akahama--Kawamura (2006) and Eremets et al. (2023)
  diamond-anvil Raman-edge scales. All 33
  identified ruby-calibrated material records now link to an exact calibration
  identifier. Public APIs convert ruby scales through the R1 wavelength ratio,
  diamond scales through the Raman wavenumber ratio, transform XRD pressure
  scales through virtual same-standard volumes, and recursively normalize a
  material EOS through its recorded calibration graph. An explicit edge
  registry preserves the sources and transformations for simultaneous and
  jointly optimized cross-calibrations. New records add the Fratanduono (2021)
  Au anchor and the Tateno (2019) and Chidester (2021) B2-KCl thermal EOSs; the
  complete 155-row Chidester KCl/Pt-derived-pressure table is bundled. The
  complete Dorogokupets--Oganov (2007) Pt four-oscillator thermal EOS is now
  implemented in Python and Rust and makes the Chidester KCl-to-Pt edge
  executable, including the published KCl-effective-to-Pt-surface temperature
  transformation.
  Observation-level re-reduction from paired measured calibrant volumes remains
  available as a separate operation. A dedicated pressure-scale normalization
  guide documents valid same-standard, ruby, and cross-material paths, the
  qualifying literature, validity and uncertainty rules, and current library
  coverage. A set-level route helper now finds and ranks XRD-standard targets
  reachable from every supplied EOS, with an executable notebook covering
  direct, optical, recursive, and temperature-transforming routes.
- Added the complete 51-row Somayazulu et al. (2023) B4C P-V-T table as a
  provenance-bearing `.eosmat` dataset, linked it to new MGD and Berman thermal
  records, and added a reproducible comparison of published, effective-
  variance, and latent-coordinate fits. The open-source EosFit engine audit
  also adds an explicit truncated-quadratic `berman` reference-volume law in
  Python, Rust, and the interchange schema. An opt-in
  `b4c_somayazulu_2023_berman_refit` record carries the EosFit7c public-data
  refit alongside, without replacing, the published parameterization; generic
  refit lineage and fit-provenance fields make that distinction machine-readable.
- Added forward DAC-confinement prediction through
  `volume_with_dac_confinement(P_cold, T, f_dac=...)`, together with an explicit
  reference-relative `thermal_pressure_increment()` for displaying the full
  thermal increment, retained confinement pressure, and total hot pressure.
  Python material records and the Rust scalar, batch, and `.eosmat` APIs expose
  the same workflow in their public volume conventions.
- Added `DoubleDebyeLogMomentHelmholtz` in Python and Rust together with the
  primary-source-validated Correa et al. (2008) diamond record. The complete
  Helmholtz model implements the Vinet cold curve, logarithmic-moment
  double-Debye weights, zero-point motion, and the published $T^2$
  anharmonic contribution, with pressure, volume, energy, entropy, and heat
  capacity evaluation.
- Added optional reference-isotherm anchoring to both double-Debye Helmholtz
  models. A numeric `Tr` subtracts the simulated ionic and anharmonic
  contribution at that temperature before adding it to `rt_eos`; `Tr=null` in
  `.eosmat` preserves the literal 0 K cold-curve formulation. The catalog now
  includes Correa- and Benedict-thermal diamond variants anchored to the
  experimental Dewaele 298 K Vinet isotherm.
- Added a shared Python/Rust error contract with domain-specific Python
  exceptions, stable machine-readable codes and context, Rust error-kind
  accessors, preserved fitting source chains, and native-extension parity.
- Added interpreted Rust `EosFitResult` parameters, latent states, profiled
  covariance, standard errors, correlation, and fit statistics.
- Added public Rust `.eosmat` validation, serialization, and save APIs with
  extension-preserving canonical and legacy round trips.
- Added PEP 561 typing metadata and a repository-wide `mypy` CI gate.
- Added array-aware pressure and temperature conversions to `peritheos.units`,
  plus conventional-cell to formula-molar-volume conversion in both directions.
- Added separate complete-Python and native Rust coverage gates and
  property-oriented native-versus-compatibility tests.
- Added workflow-oriented Rust API documentation, six executable core and
  fitting examples, docs.rs metadata, and warning-free rustdoc/doctest CI.

### Changed

- Corrected the B2-KCl refit reproductions for Walker (2002), Tateno et al.
  (2019), and Chidester et al. (2021). Walker now follows the source's staged,
  unweighted preferred
  fit. Tateno now uses the final published `gamma0=2.3`, `q=0.8`, integrated-
  Gruneisen Debye law, and the correctly aligned official Supplemental Table S1
  workbook rather than the split accepted-manuscript table. Chidester now uses
  the source's full simultaneous fit scope—123 Dewaele room-temperature rows
  plus 155 new high-temperature rows—and the integrated-Gruneisen Debye law,
  eliminating the false `q=0` boundary and recovering all five coefficients.
- EOS records now evaluate extrapolated states by default. Published ranges are
  treated as calibration/data coverage, with opt-in enforcement through
  `check_validity=True`; `within_calibration_range()` is the preferred coverage
  query and `within_validity()` remains a compatibility alias.
- Documented the factor-of-two normalization conflict between the Correa (2008)
  and Benedict (2014) diamond anharmonic terms while preserving each published
  equation and coefficient literally. For the published diamond parameters the
  coefficient is volume independent, so this discrepancy affects caloric
  quantities but not direct $P(V,T)$ or $V(P,T)$ evaluation.
- Expanded Python fallback, Rust batch, error-category, thermal-domain, and
  `.eosmat` validation tests; raised the branch-aware Python coverage floor to
  90% and the Rust line-coverage floor to 85%.
- Consolidated the public `peritheos-core` and `peritheos-fit` Rust crates into
  one `peritheos` crate. Fitting and uncertainty APIs now live under
  `peritheos::fit`; the Python API is unchanged.
- Deprecated the old `peritheos.utils` pressure and temperature conversion
  imports in favor of the consolidated `peritheos.units` API.
- Reworked the getting-started documentation around validated material records,
  clarified source-build requirements, and grouped contributor-only design
  documents outside the main user navigation.

## [0.6.0] - 2026-09-01

### Added

- Added the generic `DoubleDebyeHelmholtz` full-free-energy EOS with a Vinet
  0 K cold curve, volume-dependent double-Debye modes and weights, zero-point
  motion, a volume-dependent $T^2$ correction, analytic pressure, normal
  pressure/volume/temperature inversion, fitting support, and documented
  Benedict et al. diamond parameters as an example rather than model defaults.
  The public Rust core implements the same Helmholtz, pressure, caloric, and
  inversion equations and loads the model through native `.eosmat` dispatch.
- Added the Benedict et al. (2014) diamond coefficients as an audited material
  record in both the curated catalog and bundled `.eosmat` library, with
  explicit per-atom/conventional-cell conversion and cold-curve caveats.
- Added DAC two-volume temperature inversion for absolute Helmholtz models by
  defining the confined thermal pressure relative to the 300 K isotherm; the
  material-record wrapper accepts conventional-cell volumes directly.
- Added native Rust loading for canonical Peritheos format-3 and legacy
  Dioptas format-2 `.eosmat` files, with executable runtime-dispatched EOS
  records, preserved JSON extensions, and automatic conventional-cell to
  molar-volume conversion for energy-based thermal models.
- Added eight executable notebook tutorials covering pressure calibration,
  material-library exploration, `.eosmat` round trips, room-temperature EOS
  and gold-scale comparisons, thermal state surfaces, DAC temperature
  sensitivity, and fit-to-prediction uncertainty propagation.

### Fixed

- Corrected the scientific provenance of all eleven Sokolova pressure scales:
  their reference inputs and final coefficients originate in Sokolova et al.
  (2013), Tables 1 and 4, while the 2016 paper is the spreadsheet
  implementation/correction source. Each material record now carries structured
  source lineage and fitting-data caveats. The misleading `_sokolova_2016`
  record identifiers and public constants were removed in favor of
  `_sokolova_2013`; `Sokolova2016` remains the corrected calculator formalism.
- Replaced the migrated InN BM3/experimental-volume hybrid with Muñoz and
  Kunc's published theoretical wurtzite Murnaghan fit, including a reference
  volume reconstructed from their Table 1 theoretical lattice constants.
- Made the Campbell--Heinz B2-KCl entry an explicit two-primary-source
  composite: its absolute `V0` and propagated error now follow from the
  published B2/B1 volume ratio and Dewaele et al.'s B1 reference volume.
- Consolidated duplicate `majorite`/`mgsio3-maj` materials, corrected the
  Hanfland lithium equation family from BM3 to its combined-phase Vinet form,
  and removed the unsupported Fei-FeO and Hixson-W standalone BM3 records.
- Resolved the remaining primary-source blockers for CsCl, magnetite, Li,
  majorite, MW60 magnesiowuestite, NiS, phase D, cubic SnO2, and SrO B1/B2.
  Corrected equation orders, cell conventions, phases, fitted-versus-fixed
  flags, ranges, and all printed parameter errors. Phase D now has distinct
  AntA and AntB reference-volume records, and the inherited unsupported
  majorite error and tenfold SnO2 volume-error transcription are removed.
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
- Replaced the migrated Anderson et al. Au `AlphaKT` approximation with the
  exact Equation (29) logarithmic-volume linear thermal pressure, restored its
  density-derived reference volume, Table V domain, and partial published
  `(dKT/dT)V` uncertainty, and promoted the record after primary-source review.
- Removed the Martinez et al. aragonite global HT-BM3 record: Table 7 omits
  its fitted `V0(298 K)`, and the remaining coefficients do not reproduce the
  printed 64-point table under documented pressure- or volume-residual least
  squares. The independently reproducible staged BM2 result is retained and
  extended with its published Equation (2) `K0(T)` slope and Equation (3)
  direct-linear reference-volume law.
- Promoted the Scott et al. (2001) cementite BM3 record after checking the
  complete primary article: the ambient `V0`, weighted-fit coefficients, all
  printed errors, 300 K reference state, and compression interval are now
  traced to pages 1875--1877. The measured `V0` is explicitly fixed in the
  fit metadata, and the unreported covariance/confidence convention remains
  documented rather than inferred.
- Corrected the Clendenen and Drickamer (1966) CoO record from an inherited
  BM3 representation to the published Murnaghan Equation 4, restored the
  Table II ambient cell and Table III 30.8 GPa range, and retained null errors
  because the primary paper reports no parameter uncertainty or covariance.
- Promoted the Noguchi et al. (1999) NiO shock-derived 300 K BM3 isotherm after
  checking the official 1998 primary conference paper for the sample reference
  lattice, its propagated uncertainty, and the Mie--Gruneisen reduction. The
  final journal article supplies the 147.6 GPa range and `K0`/`K0'`; their
  errors remain null because the authors do not report them.

### Added (catalog and native backend)

- An executable documentation notebook using the complete printed Martinez et
  al. (1996) aragonite Table 3 dataset to demonstrate 298 K and staged-isotherm
  BM2 fitting, thermal-trend recovery, scaled joint P-V-T fitting, residual
  visualization, and uncertainty/chi-square interpretation.
- A reproducible primary-source audit for the curated migrated `.eosmat` EOS
  records. All 147 bundled records are now directly validated against original
  publications, official supplements, or stable institutional reports, and no
  record remains pending or deferred. The
  bundled machine-readable ledger records source locations, the independently
  checked Shen--Smith (2026) Vinet fits and errors, the B4C order inconsistency,
  and the restored Hanfland graphite `V0` uncertainty.
- A mechanism-oriented `ThermalReferenceStateEOS` implementation for the
  temperature-dependent `V0(T)`/`K0(T)` formulation used by the validated ice
  VI/VII records. It now supports a generic `thermal_expansion_law`, including
  exact analytical integration of `alpha0+alpha1*T`; the constant law remains
  backward compatible. A separate `reference_volume_law="linear_temperature"`
  represents a directly linear mean-expansion relation without conflating it
  with integrated instantaneous expansivity. The Dioptas `AlphaKT` interchange
  type maps to the canonical `thermal_reference_state` model identifier.
- Primary-source-validated native material records for the Martinez et al.
  (1996) staged aragonite BM2 P-V-T parameterization and the Dewaele et al. (2012) B2-KCl
  P-V-T pressure calibration. KCl uses the paper's Vinet reference isotherm,
  additive thermal-pressure term, fixed fictive `V0`, and explicit
  experimental-versus-computational validity provenance. It is the preferred
  `kcl.eosmat` record. After primary-source corrections and duplicate removal,
  the catalog now contains 147 records, all validated.
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
- All eleven Sokolova thermal pressure markers: MgO, diamond, Al, Cu, Ag, Au,
  Pt, Nb, Ta, Mo, and W, with the original 2013 fit provenance and the
  corrected 2016 workbook equations.
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
- A Peritheos-owned flat `.eosmat` format 3, normative JSON Schema, complete
  115-material/147-record EOS database migrated from Dioptas 0.10.0 with
  explicit validation status and provenance, legacy Dioptas format-2 input,
  and tested Dioptas 0.10.0 read compatibility. A dedicated schema reference
  documents every field, discriminator pairing, default, unit, validation
  status, and consumer compatibility rule.
- A Rust workspace containing native EOS, fitting, uncertainty, and private
  PyO3 binding crates, with Rust 1.83 as the library MSRV.
- Shared Python/Rust compatibility fixtures and migration baselines for all
  isothermal and thermal model families.
- Multi-platform native-wheel release jobs for supported CPython versions on
  Linux x86-64/ARM64, macOS Intel/Apple Silicon, and Windows x86-64.
- Pull-request wheel build and isolated-install smoke tests on Linux, macOS,
  and Windows, complementing the full tagged-release wheel matrix.
- Dependency-free public Rust batch traits, typed joint EOS fitting, and
  model-aware linear and Monte Carlo uncertainty entry points.
- Package-contained scientific fixtures and a two-crate archive verifier that
  tests the required core-before-fit crates.io publication sequence.
- A pinned Rust dependency-source and SPDX-license audit covering the supported
  Linux, macOS, and Windows target graphs.

### Changed

- Built-in Python EOS classes now preserve their public API while delegating
  evaluation, inversion, thermoelastic, and caloric calculations to Rust.
- The material catalog's linear, logarithmic-volume, configurable
  reference-state, variable-exponent Debye, asymptotic-power-law Debye, and
  generic multi-oscillator mechanisms now use the same native evaluation and
  fitting architecture. Thermal fits accept fixed categorical equation choices
  through `configuration`, and native linear uncertainty supports records with
  measurement errors but no published parameter covariance.
- Named bounded robust fitting losses and uncertainty propagation statistics
  now use native numerical kernels. Custom reference EOS classes, callable
  fitting losses, and NumPy-seeded Monte Carlo draws retain documented
  compatibility paths.
- Native fits now return their profiled global-parameter covariance directly;
  Python no longer recomputes it through a separate SciPy/NumPy path.
- Errors-in-variables fits now use colored latent-coordinate Jacobians and a
  block Schur-complement solve, with stress coverage for large and rank-deficient
  datasets.
- Birch-Murnaghan kernels now use an algebraically equivalent cube-root form,
  and large independent EOS arrays use deterministic thresholded parallel
  evaluation while the Python interpreter lock is released.
- Holzapfel bulk-modulus derivatives now execute directly in Rust, obsolete
  Python natural-strain coefficient formulas are removed, and the historical
  coefficient-level Holzapfel helper remains available through a Rust-backed
  compatibility wrapper.
- Native batch calls now have concurrent large-array stress coverage in
  addition to deterministic order, shape, stride, and round-trip checks.
- Native least squares now equilibrates differently scaled Jacobian columns,
  reports failed steps as failures instead of false `xtol` convergence, and
  uses a rank-aware Moore-Penrose covariance calculation.
- Latent-coordinate covariance profiling now preserves observation-local
  blocks, avoiding the previous dense cubic post-fit calculation.
- Native EOS dispatch is restricted to exact built-in Python classes so
  subclass overrides and Debye/Einstein model identity remain authoritative.
- Linear uncertainty kernels now reject non-positive-semidefinite parameter
  covariance and negative state variance instead of clamping invalid inputs.
- Tagged releases now include CPython 3.14 free-threaded wheels alongside the
  standard CPython wheel matrix.

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

[Unreleased]: https://github.com/CPrescher/peritheos/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/CPrescher/peritheos/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/CPrescher/peritheos/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/CPrescher/peritheos/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/CPrescher/peritheos/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/CPrescher/peritheos/releases/tag/v0.2.0
