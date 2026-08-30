# Development

## Test suite

```bash
uv sync --all-groups
uv run ruff check .
uv run ruff format --check .
uv run pytest -q -W error --cov --cov-report=term-missing
```

The Python package contains a private PyO3 extension and therefore requires a
Rust toolchain for editable or source installs. The workspace MSRV is Rust
1.83; the release build uses maturin 1.15.

Run the native gates directly with:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features --locked
```

Python coverage omits `peritheos/eos/**` because those modules are facades over
the PyO3 backend and contain retained custom-model compatibility paths. Their
built-in numerical implementation is covered by the Rust unit, fixture,
literature, and integration tests. Python coverage continues to enforce the
90% branch-aware floor for fitting orchestration, uncertainty, units, and the
rest of the Python layer.

If the default uv cache is unavailable in a sandbox:

```bash
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run pytest -q
```

## Documentation

```bash
uv run --group docs mkdocs build --strict
```

## Native architecture and compatibility paths

`peritheos-core` owns all built-in isothermal, thermal, caloric, inversion, and
quadrature calculations. `peritheos-fit` owns bounded robust least squares,
typed EOS residual construction, structured latent-coordinate solving,
covariance profiling, delta-method propagation, and Monte Carlo summary
statistics. `peritheos-python` exposes these through the private
`peritheos._rust` module; application code should continue importing the
documented Python modules.

Named-loss fits of exact built-in Peritheos models cross the PyO3 boundary
once. Model reconstruction, EOS evaluation, residual weighting, correlated
whitening, latent-coordinate assembly, finite differences, and optimization
then remain in Rust, with the Python interpreter lock released for the solve.
Custom Python `EosBase` classes and subclasses retain the callback evaluation
path so overridden scientific behavior is not bypassed. Callable fitting
losses retain the SciPy solver because arbitrary Python loss functions cannot
be represented by the native loss enum. Monte Carlo draws continue to use
NumPy so seeded results and invalid-sample rejection remain compatible;
accepted-sample statistics are native.

Large independent native arrays use Rayon after workload-specific size
thresholds. The binding copies the NumPy inputs before releasing the Python
interpreter lock, preserving NumPy memory safety even if another Python thread
mutates the original array. Parallel results are collected in input order and
errors are converted sequentially, so scheduling does not change output order
or which input error is reported first. Small arrays use the serial path.
`RAYON_NUM_THREADS` may be used by embedding applications that need to limit
the shared Rayon worker pool.

## Building distributions

```bash
uv build
uv run --group release twine check dist/*
```

Release CI builds interpreter-specific wheels for Python 3.9 through 3.14 on
manylinux x86-64 and ARM64, macOS Intel and Apple Silicon, and Windows x86-64.
It builds the source distribution separately and smoke-tests an installed
wheel outside the checkout before publication.

## Adding an isothermal EOS

Subclass `EosBase`, validate `V0` and modulus parameters using the shared
validators, and implement analytic `pressure()` and `bulk_modulus()` methods.
The inherited volume solver supplies array-aware inversion. Add reference-state,
derivative, array, invalid-input, and round-trip tests.

## Adding a thermal EOS

Subclass `ThermalEOS` and implement `thermal_pressure(V, T)`. The base class
provides total pressure, inversion, isothermal bulk modulus, compressibility,
and expansivity. Both volume and temperature inversion are array-aware.
Implement `molar_heat_capacity_v()` only when the model has a defined caloric
potential; the base class then provides `C_P` and `K_S`.

Models with expensive fixed-volume preparation may override the private
`_thermal_pressure_function(V)` hook. It returns a temperature-only callable
used by both inversion methods; Sokolova uses it to evaluate its volume integral
once per solved state.

Document the required molar-volume unit, reference temperature, pressure unit,
parameter domain, source equations, and whether energy methods are absolute or
reference contributions.
