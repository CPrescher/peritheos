# Rust migration completion plan

This plan closes the remaining gap between the working native backend and a
compatibility-ready Rust and Python release. Work stays on
`codex/rust-core`; `main` is not changed until every release gate passes.

## Current ownership matrix

| Capability | Current owner | Completion work |
|---|---|---|
| Built-in isothermal, thermal, and caloric equations | Rust | Audit and remove obsolete Python duplicates |
| Scalar inversion and quadrature | Rust for built-ins; Python/SciPy for custom models | Preserve the explicit custom-model fallback |
| NumPy array evaluation | PyO3 with Rust scalar kernels and thresholded Rayon | Add a dependency-free public Rust batch API |
| Named-loss built-in fitting | Rust end to end | Add an ergonomic joint-fitting Rust entry point and broader parity cases |
| Callable-loss and custom-model fitting | SciPy or a Python callback | Retain as documented compatibility behavior |
| Fit covariance | Implemented in Rust but still recomputed by Python/SciPy | Return and consume the native parameter covariance |
| Linear uncertainty matrix propagation | Rust | Add model-aware finite-difference interfaces |
| Monte Carlo evaluation | Python/NumPy with Rust summary statistics | Add a deterministic Rust sampling/evaluation engine without changing NumPy-seeded compatibility |
| Python result objects and serialization | Python | Keep as the stable public presentation layer |
| Crate and wheel packaging | Scaffolded and locally validated | Package crates, validate licenses, and exercise the full wheel matrix |

## Checkpoint sequence

1. **Rust batch API.** Add ordered scalar-equivalent batch methods without an
   array-framework or threading dependency in `peritheos-core`.
2. **Rust fitting and uncertainty ergonomics.** Add explicit joint fitting,
   model-aware finite differences, covariance, and Monte Carlo entry points.
3. **Python integration.** Use Rust covariance and uncertainty kernels for
   exact built-in models; retain custom classes and callable losses as named
   fallbacks.
4. **Single-source cleanup.** Remove unused built-in Python formulas and keep
   custom-model compatibility code isolated and tested.
5. **Adversarial validation.** Cover every loss, correlated observations,
   bounds, rank loss, large latent volume/temperature blocks, deterministic
   ordering, free-threaded use, and concurrent batch calls.
6. **Release engineering.** Validate `cargo package`, API docs, dependency
   licenses, wheel size, source distributions, and supported wheels.
7. **Hosted evidence.** Push only with explicit approval, require every CI and
   wheel job to pass, then review a release candidate before any merge.

## Definition of done

The migration is complete only when the validation gates in the
[Rust migration contract](rust-migration.md#validation-gates) pass, the public
Rust APIs are documented and packageable, intentional Python fallbacks are
listed rather than inferred, and release-wheel tests pass on Linux, macOS, and
Windows. SIMD or SciRS2 adoption is not a completion requirement; either must
independently satisfy the numerical and MSRV gates before use.
