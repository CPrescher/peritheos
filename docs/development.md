# Development

## Test suite

```bash
uv sync --all-groups
uv run ruff check .
uv run ruff format --check .
uv run --python 3.9 mypy
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
cargo llvm-cov --workspace --exclude peritheos-python --all-features --fail-under-lines 80
scripts/verify-rust-packages.sh
```

The package verifier builds the public `peritheos` archive and runs the tests
from its extracted, normalized package. Set `PERITHEOS_ALLOW_DIRTY=1` when
checking an intentional uncommitted change; release verification remains
strict by default.

`deny.toml` restricts normal and platform-specific Rust dependencies to
crates.io and an explicit permissive-license allowlist. The pinned
`cargo-deny` CI action evaluates the complete feature graph on Linux, macOS,
and Windows targets.

Python coverage includes the complete package, including the EOS facades and
the custom-model compatibility paths that do not use the native evaluator. The
branch-aware floor is 89%; this broader denominator replaces the former 90%
figure that excluded the complete `peritheos/eos/**` tree. Deterministic
property-oriented tests compare native and compatibility implementations over
state grids and exercise round-trip and thermodynamic identities.

Rust line coverage covers the unified `peritheos` crate and must remain at
least 80%. `peritheos-python` is excluded from that number because its PyO3
entry points are exercised by the Python suite rather than Rust's test
harness.

The installed package contains `py.typed`; `mypy` checks the inline annotations
on every change. Private native-extension attributes and the intentional
isothermal/thermal override signatures are excluded from internal diagnostics,
while the concrete public classes retain their precise signatures for users.

If the default uv cache is unavailable in a sandbox:

```bash
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run pytest -q
```

## Documentation

```bash
uv run --group docs mkdocs build --strict
```

## Native architecture and compatibility paths

The public `peritheos` crate owns all built-in isothermal, thermal, caloric,
inversion, quadrature, fitting, and uncertainty calculations. Its `fit` module
contains bounded robust least squares, typed EOS residual construction,
structured latent-coordinate solving, covariance profiling, model-aware
delta-method propagation, and deterministic Monte Carlo sampling and summary
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
accepted-sample statistics are native. Native Rust clients can instead use the
public Rust sampler; its deterministic stream intentionally does not claim
bit-for-bit compatibility with NumPy.

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

Maturin is the sole distribution backend. Its include rules explicitly place
the normative schema and complete `peritheos/data/materials` library in both
wheels and source distributions. Documentation sources, including executable
notebooks and their local datasets, are source-distribution inputs; building
the documentation from an sdist therefore has the same content as a checkout.
The inactive setuptools package-data table was removed when the native backend
became authoritative.

Release CI builds interpreter-specific wheels for Python 3.9 through 3.14,
including the separate CPython 3.14 free-threaded ABI, on manylinux x86-64 and
ARM64, macOS Intel and Apple Silicon, and Windows x86-64. It builds the source
distribution separately and smoke-tests every installed wheel outside the
checkout before publication.

Pull-request CI also builds and installs a real Python 3.13 wheel outside the
checkout on Linux, macOS, and Windows. This catches artifact composition,
platform tagging, and dynamic-library loading failures that an editable native
test environment would not expose.

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
