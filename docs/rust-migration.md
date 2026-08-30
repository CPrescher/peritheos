# Rust migration contract

This document defines the compatibility and validation contract for moving
Peritheos's built-in numerical implementation to Rust. It is an engineering
contract, not a change to the scientific conventions documented in the
[equation reference](equation-reference.md).

The migration is developed on an integration branch. The Python implementation
remains authoritative until every applicable compatibility gate below passes.
Built-in equations must have one implementation after the migration: the Rust
core. Python remains the user-facing compatibility and convenience layer.

## Objectives

- Provide a native Rust API for all built-in equations of state.
- Preserve the documented Python API and NumPy behavior.
- Make the Rust implementation the single source of truth for built-in models.
- Provide native fitting and uncertainty functionality equivalent to the
  current SciPy-backed implementation.
- Improve repeated scalar, inversion, quadrature, and Monte Carlo workloads
  without materially regressing already-vectorized analytical calculations.
- Ship tested Python wheels on every supported platform and Python version.

The migration does not change equations, parameter definitions, working units,
reference states, fitting statistics, covariance conventions, or uncertainty
assumptions. A scientific correction discovered during the migration must be a
separate, explicitly documented change with an independent reference case.

## Crate and binding boundaries

The intended workspace contains three crates:

| Crate | Responsibility | Dependency policy |
|---|---|---|
| `peritheos-core` | EOS models, validation, thermodynamic properties, scalar and batch evaluation, inversion, and required quadrature | No Python dependency; numerical dependencies kept small |
| `peritheos-fit` | Fitting, covariance, and uncertainty kernels | Optional; solver and linear-algebra dependencies remain private implementation details |
| `peritheos-python` | Private PyO3 module used by the existing Python package | May depend on PyO3 and NumPy; not a public Python API |

The public Rust API must not expose SciRS2-specific result, array, error, or
solver types. This preserves the ability to replace individual numerical
components without breaking downstream Rust users.

## Supported model contract

The following built-in isothermal models must be available in both Rust and
Python:

| Model | Constructor parameters | Required quantities |
|---|---|---|
| `BM2` | `V0, K0` | pressure, bulk modulus, volume |
| `BM3` | `V0, K0, K0_prime` | pressure, bulk modulus, volume |
| `BM4` | `V0, K0, K0_prime, K0_double_prime` | pressure, bulk modulus, volume |
| `Murnaghan` | `V0, K0, K0_prime` | pressure, bulk modulus, volume |
| `NaturalStrain2` | `V0, K0` | pressure, bulk modulus, volume |
| `NaturalStrain3` | `V0, K0, K0_prime` | pressure, bulk modulus, volume |
| `NaturalStrain4` | `V0, K0, K0_prime, K0_double_prime` | pressure, bulk modulus, volume |
| `ModifiedTait` | `V0, K0, K0_prime, K0_double_prime` | pressure, bulk modulus, volume |
| `Vinet` | `V0, K0, K0_prime` | pressure, bulk modulus, volume |
| `Holzapfel` | `V0, K0, K0_prime, n, Z` | pressure, bulk modulus, bulk-modulus derivative, volume |

The following thermal models and the documented compatibility alias must be
available:

| Model | Reference EOS constraint | Own parameters |
|---|---|---|
| `MieGruneisenDebye` | any built-in isothermal EOS | `Tr, theta0, gamma0, q, n` |
| `MieGruneisenEinstein` | any built-in isothermal EOS | `Tr, theta0, gamma0, q, n` |
| `ThermalModifiedTait` | `ModifiedTait` | `Tr, theta, alpha0, n` |
| `HollandPowell2011` | alias of `ThermalModifiedTait` | unchanged |
| `Sokolova2016` | `Holzapfel` | all required and optional parameters documented in the API reference |

Every thermal model must preserve total and thermal pressure, P-T to V
inversion, P-V to T inversion, DAC two-volume inversion, isothermal bulk
modulus and compressibility, thermal expansivity, and every caloric quantity
currently supplied by that model.

## Python API contract

The following behavior is public and must remain unchanged:

- Existing import locations and exported class/function names.
- Constructor parameter names, order, keyword support, and defaults.
- Public parameter attributes such as `V0`, `K0`, and `K0_prime`.
- `parameter_values(include_reference=True)` constructor ordering.
- Dotted reference parameters such as `rt_eos.K0`.
- Immutable reconstruction through `with_parameters()`.
- `volume()` and its supported `calculate_volume()` alias.
- `temperature()` and its supported `calculate_temperature()` alias.
- Python `float` output for scalar state input.
- NumPy array output and NumPy broadcasting for array state input.
- Array shape preservation, including zero-dimensional and multidimensional
  inputs accepted by the current API.
- `ValueError`, `TypeError`, or `ArithmeticError` categories for equivalent
  failures. Exact private error and optimizer messages are not stable.
- Support for custom Python `EosBase` subclasses through a compatibility path,
  even when they cannot use every native fast path.

The private extension module will be named `peritheos._rust`. Users must not
need to import it directly.

## Numerical and unit contract

- Public pressure and bulk modulus values use GPa.
- Temperatures use kelvin and must be finite and positive.
- Volumes must be finite and positive.
- Isothermal models other than Holzapfel accept any volume unit consistent with
  `V0`.
- Holzapfel and thermal models use molar volume in `J bar^-1 mol^-1`.
- Thermal energy and heat-capacity quantities retain their documented molar SI
  units.
- The Rust core will freeze each physical constant to a documented value so
  results do not change with the installed SciPy version.
- Reference identities, derivative identities, literature cases, and model
  domain restrictions remain as documented in the validation guide.

Volume inversion must retain the branch nearest the reference volume. It must
search toward compression for pressures above the reference pressure and along
the first expansion branch for lower pressures. Temperature inversion must
return the positive-temperature root nearest `Tr`. Returned roots require a
post-solve residual check.

Compatibility is tolerance based rather than bit-for-bit. Each ported model
must pass the existing tolerances before its Python facade switches to Rust.
Any necessary tolerance change requires an explained numerical analysis.

## Fitting contract

The Rust fitting layer must eventually support the complete signatures of
`fit_rt_eos`, `fit_thermal_eos`, and `fit_joint_eos`, including:

- arbitrary subsets of free and fixed constructor parameters;
- named box bounds;
- independent pressure, volume, and temperature standard deviations;
- the legacy `sigma` alias for pressure uncertainty;
- per-observation correlated covariance matrices;
- latent volume and temperature variables for errors-in-variables fitting;
- `linear`, `soft_l1`, `huber`, `cauchy`, and `arctan` losses;
- `f_scale`, `max_nfev`, and absolute-versus-scaled covariance;
- joint reference/thermal parameters and complete cross-covariance;
- dense small fits and sparse large latent-variable fits; and
- the current validation and failure categories.

`FitResult` must retain parameter ordering, covariance and correlation,
standard errors, adjusted states and corrections, residual arrays, chi-square,
reduced chi-square, degrees of freedom, AIC, BIC, convergence metadata,
human-readable summaries, and schema-version-1 JSON export. Exact optimizer
messages and the last floating-point bit are not compatibility requirements.

## Uncertainty contract

The native implementation must preserve:

- named partial parameter covariance;
- independent parameter-error assumptions;
- supplied correlation matrices;
- symmetric positive-semidefinite covariance validation;
- local linear propagation using parameter and state derivatives;
- optional full output covariance;
- normal confidence limits;
- Monte Carlo parameter and state sampling;
- invalid-sample rejection and the reported rejected fraction;
- deterministic seeded behavior within a backend; and
- explicit block-independence when combining separately fitted models.

Identical pseudorandom samples across NumPy and Rust are not required unless a
shared RNG algorithm is later made part of the public contract. Statistical
results must agree within test-defined sampling tolerances.

## Validation gates

A phase is complete only when its applicable gates pass:

1. Rust unit tests for reference states, limits, invalid domains, and errors.
2. Shared language-independent literature fixtures.
3. Rust/Python deterministic grid comparisons.
4. Analytical-versus-numerical derivative comparisons.
5. Pressure-volume and pressure-volume-temperature round trips.
6. Scalar, array, and broadcasting behavior through Python bindings.
7. Adversarial fitting cases: bounds, outliers, correlated observations,
   ill-conditioning, rank loss, and large latent systems.
8. Linear and Monte Carlo uncertainty comparisons.
9. Linux, macOS, and Windows wheel installation and smoke tests.
10. Supported Python-version tests, including a separately evaluated
    free-threaded build strategy.
11. Benchmarks showing no material regression in common analytical array calls
    and improvements in targeted repeated-scalar workloads.
12. Strict documentation, lint, coverage, Rust formatting, Clippy, and package
    validation.

The legacy built-in Python equations may be removed only after every gate is
satisfied with the Rust backend enabled by default. Until then, migration work
is additive and remains on the integration branch.

## Baseline

The migration begins from Peritheos 0.5.0 at Git commit
`12d033378418cfd6c9ece6050c550fc748ffe02a`. On macOS ARM64 with CPython 3.13,
the clean worktree baseline contains 382 passing tests, passes Ruff lint and
format checks, builds the documentation in strict mode, and produces the
pure-Python wheel and source distribution. Reproducible operation-level timing
is provided by `benchmarks/python_baseline.py`; machine-specific results belong
under `benchmarks/baselines/` and are evidence, not performance promises.

## Integration-branch implementation status

The `codex/rust-core` integration branch now contains all three planned crates.
Built-in isothermal and thermal Python classes route their numerical work to
`peritheos-core`, including inversion, derivatives, caloric properties, and
the full Sokolova expression. Mie-Gruneisen models retain a Python fallback
only when they wrap a user-defined `EosBase` implementation.

The public fitting functions use the `peritheos-fit` bounded solver for all
five named loss functions. Their validation, latent-variable residual model,
correlated-observation whitening, result construction, and JSON schema remain
in the Python compatibility layer. Errors-in-variables fits use a structured
native path: simultaneous finite-difference coloring evaluates all points for
one latent coordinate together, and observation-local normal blocks are
profiled through a Schur complement. Callable Python losses intentionally
retain SciPy's solver. Delta-method matrix propagation and accepted Monte Carlo
sample statistics are native; NumPy retains seeded random draws and Python
retains model reconstruction and invalid-sample rejection.

These boundaries preserve extensibility without duplicating built-in equations
or changing a scientific convention. The branch must remain unmerged until its
multi-platform wheel workflow and all validation gates have completed in CI.

An indicative quick run on macOS ARM64 compared with the committed pre-Rust
baseline improved every recorded workload. Median speedups were approximately
8x for scalar BM3 pressure, 27x for its array path, 1,741x for volume inversion,
226x for Debye pressure, 945x for Sokolova pressure, 1.5-1.7x for fitting, 1.9x
for linear uncertainty, and 24x for Monte Carlo uncertainty. These figures are
machine- and environment-specific evidence, not release performance promises;
CI and release decisions must also consider the deterministic numerical gates.
