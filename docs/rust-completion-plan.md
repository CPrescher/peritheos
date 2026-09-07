# Rust migration completion record

This historical record describes the final migration from the earlier Python
implementation to the unified native backend. The work is complete on `main`;
the current public API and model inventory are documented in the [Rust API
guide](rust-api.md), [API reference](api.md), and [model list](models.md).

## Final ownership matrix

| Capability | Current owner | Maintained boundary |
|---|---|---|
| Built-in isothermal, thermal, and caloric equations | Rust | Duplicate built-in formulas audited and removed; keep custom-model fallbacks isolated |
| Scalar inversion and quadrature | Rust for built-ins; Python/SciPy for custom models | Preserve the explicit custom-model fallback |
| NumPy array evaluation | PyO3 with Rust scalar kernels and thresholded Rayon | Public dependency-free Rust batch traits complete; broaden stress coverage |
| Named-loss built-in fitting | Rust end to end | Public joint-fitting entry point complete; broaden parity cases |
| Callable-loss and custom-model fitting | SciPy or a Python callback | Retain as documented compatibility behavior |
| Fit covariance | Rust, consumed directly for native fits | Retain the SciPy calculation only for its explicit fallback paths |
| Linear uncertainty matrix propagation | Rust | Public model-aware finite-difference and propagation interfaces complete |
| Monte Carlo evaluation | Public Rust engine; Python/NumPy compatibility path | Preserve NumPy-seeded Python results while validating native sampling |
| Python result objects and serialization | Python | Keep as the stable public presentation layer |
| Crate and wheel packaging | Packageable crates and wheel matrix | Repeat dependency, license, package, and wheel gates for every release |

## Completed checkpoint sequence

1. **Rust batch API.** Added ordered scalar-equivalent batch methods without an
   array-framework or threading dependency in `peritheos`.
2. **Rust fitting and uncertainty ergonomics.** Added explicit joint fitting,
   model-aware finite differences, covariance, and Monte Carlo entry points.
3. **Python integration.** Routed exact built-in models through Rust covariance
   and uncertainty kernels while retaining named custom-model and callable-loss
   fallbacks.
4. **Single-source cleanup.** Removed duplicate built-in Python formulas and
   kept custom-model compatibility code isolated and tested.
5. **Adversarial validation.** Covered losses, correlated observations, bounds,
   rank loss, large latent volume/temperature blocks, deterministic ordering,
   free-threaded use, and concurrent batch calls.
6. **Release engineering.** Validated package metadata, API docs, dependency
   licenses, wheel size, source distributions, and supported wheels.
7. **Hosted evidence.** Required every CI and wheel job to pass before release.

## Definition of done

The migration was considered complete when the validation gates in the
[Rust migration contract](rust-migration.md#validation-gates) passed, the public
Rust APIs were documented and packageable, intentional Python fallbacks were
listed rather than inferred, and release-wheel tests passed on Linux, macOS,
and Windows. Future SIMD or SciRS2 adoption is not part of migration completion
and must independently satisfy the numerical and MSRV gates before use.
