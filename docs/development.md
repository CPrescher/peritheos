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

## Building the two Python distributions

The repository retains complete curation inputs. Editable installs keep those
inputs for maintenance; published wheels and source archives are separate,
explicit exports made by the PEP 517 backend:

```bash
uv build --out-dir dist
uv build packages/peritheos-datasets --out-dir dist
uv run --group release twine check dist/*
```

The first command builds the library, structural CIFs and runtime material
records without experimental observations or research archives. Dataset
declarations retain IDs, columns, sources, checksums and fit links; inline
observations become external resource declarations. Coefficients, model laws,
uncertainties, calibration conventions and scientific statuses are preserved.
Detailed audit inputs stay in the authoritative repository; exported records
link that evidence and record its document fingerprint. No source files are
rewritten by building. Building the core source archive also works independently
of the optional package and the research checkout.

The second command builds a pure Python data wheel and its own standalone
source archive. `packages/peritheos-datasets/datasets.json` is the reviewed
allowlist, with exact declaration fingerprints and source scopes. Changed
declarations or resource hashes fail the export. Adding datasets requires a
source-rights review and an explicit allowlist edit; neighboring artifacts are
not copied. The initial selection has 44 entries. Mixed author/literature,
unresolved and NC/ND sources remain excluded.

`peritheos[datasets]` requests the optional dependency declared under
`[project.optional-dependencies]`; `datasets` is not a second library build. Pin
both project versions and the extra together when preparing a release. Configure
the new `peritheos-datasets` PyPI project with the same trusted-publisher workflow
and environment before publishing the first split release. Local builds do not
publish anything. Use `uv build` for standard builds; direct `maturin build` in
the unexported research checkout bypasses the Python build backend. The release
workflow explicitly stages its runtime tree before invoking maturin for each
platform.

Installed-package verification must run outside the checkout, first with only
the core wheel (`scripts/smoke_test_wheel.py`), then after installing the extra
from local artifacts (`scripts/smoke_test_datasets.py`). The CI packaging job
checks both. The complete repository tests continue to check the preserved
original evidence through an editable installation.

Run the native gates directly with:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features --locked
cargo llvm-cov --workspace --exclude peritheos-python --all-features --fail-under-lines 85
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
branch-aware floor is 90%; this broader denominator includes the complete
`peritheos/eos/**` tree. Deterministic
property-oriented tests compare native and compatibility implementations over
state grids and exercise round-trip and thermodynamic identities.

Rust line coverage covers the unified `peritheos` crate and must remain at
least 85%. `peritheos-python` is excluded from that number because its PyO3
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

## Continuous integration

Routine pushes to `main` and pull requests run four Linux jobs:

- the core compatibility suite on Python 3.13;
- workspace tests, formatting and Clippy on stable Rust;
- Python lint and formatting on Python 3.13;
- typing on Python 3.9, the minimum supported Python version.

The full matrix runs 21 jobs. It runs the complete Python suite on versions
3.9 through 3.14, with the Python 3.13 run enforcing the 90% coverage floor.
It adds Rust 1.83, Rust documentation and crate-package verification, the 85%
Rust coverage floor, dependency license/source policy, strict documentation,
Linux distribution validation, macOS/Windows tests and installed wheels,
free-threaded Python 3.14, and minimum direct dependencies on Linux and macOS.

This full matrix runs weekly on Monday at 05:23 UTC on the default branch.
To run it on a particular branch, open **Actions → CI → Run workflow**, select
the branch, and enable **Run the full compatibility matrix**. Require this full
run to pass on the intended release commit before tagging a release. A manual
run without that option uses the four-job routine matrix. Complete literature
replays and aggregate coverage are checked in full runs, rather than on every
push or pull-request update.

In full runs, Python 3.13 coverage also supplies its regular test result, and
Linux package validation supplies its wheel smoke test, so neither needs a
duplicate job.
Rust dependency builds and the pinned coverage tool are cached, with separate
keys for Python ABIs. Numerical libraries use one thread per process to avoid
CPU oversubscription on small hosted runners. Superseded runs in the same
branch and matrix mode are cancelled; a full compatibility run is kept separate
from routine runs.

### Test selection and profiling

The routine Python job uses the explicit module list in
`.github/pytest-compatibility.ini`. It covers the public API, isothermal and
thermal models, native/Python parity, fitting, uncertainty, units, pressure
calibration, and EOSMAT/catalog/data compatibility. Paper-specific literature
replays and generated-report checks run in the complete weekly and manual
matrix, including the Python 3.13 coverage job. Add tests for new
library interfaces to this module list as part of implementing them.

Run the same core selection locally from the repository root with:

```bash
uv run pytest -c .github/pytest-compatibility.ini --rootdir=. -q -W error --durations=30
```

Normal discovery is limited to `tests/`, so retained snapshots in temporary
directories are not collected. Every CI test job prints its 30 slowest test
phases. To profile the complete suite and retain per-test timing for comparison:

```bash
uv run pytest -q -W error --durations=30 --junitxml=/tmp/peritheos-test-profile.xml
```

Tests can request `reproduce_once` and call `reproduce_once(reproduce)` for an
unchanged, no-argument audit. It computes a fresh result once per callable per
session, then returns a deep copy to each consumer so edits cannot leak between
tests. The cache is held only in test memory; saved reports and production audit
functions are unchanged. Tests that alter input files, patch a solver, or verify
repeatability must call their calculation directly instead.

## Documentation

After editing material documents, keep the readable layout and generated views
current:

```bash
python -m scripts.format_eosmat
python -m scripts.generate_material_summaries
python -m scripts.generate_catalog_coverage
python -m scripts.generate_validation_summary
```

Each command accepts `--check` for read-only verification. See
[Reading and editing EOSMAT](eosmat-authoring.md) for the simple authoring template.

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

The bounded tabulated pressure models use Python interpolation and inversion
with native Debye oscillator kernels. Their Python fitting path therefore uses
callbacks. Rust clients evaluate the same pressure surfaces directly through
the public crate and EOSMAT loader; cross-language checkpoints and bounded
round trips cover both implementations.

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

The full CI matrix also builds and installs a real Python 3.13 wheel outside
the checkout on Linux, macOS, and Windows. This catches artifact composition,
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
expansivity, and forward DAC-confined volume prediction. Both volume and
temperature inversion are array-aware. If `thermal_pressure()` is an absolute
free-energy contribution that does not vanish at the reference temperature,
also override `thermal_pressure_increment(V, T)` with the pressure difference
above the complete reference-temperature isotherm.
Implement `molar_heat_capacity_v()` only when the model has a defined caloric
potential; the base class then provides `C_P` and `K_S`.

Models with expensive fixed-volume preparation may override the private
`_thermal_pressure_function(V)` hook. It returns a temperature-only callable
used by both inversion methods; Sokolova uses it to evaluate its volume integral
once per solved state.

Document the required molar-volume unit, reference temperature, pressure unit,
parameter domain, source equations, and whether energy methods are absolute or
reference contributions.

## Maintaining Rust distribution notices

The [third-party notice audit](research/rust-third-party-notices.md) describes
source verification and toolchain scope. Repository maintenance commands:

```bash
python scripts/generate_rust_notices.py --check
```

Regenerate and review the bundle after Cargo.lock changes, including a
workspace-version bump, before preparing distributions:

```bash
python scripts/generate_rust_notices.py
```

Official publication pins the audited Rust compiler and verifies its notice
bundle. This does not raise the library MSRV. Other-toolchain redistributors
must refresh/review matching standard-library notices as documented in the
audit. License/notice copies are byte-exact; only the Cargo.lock fingerprint
normalizes checkout CRLF to LF.

## Preparing a release

The Python version in `peritheos/__init__.py`, workspace version in
`Cargo.toml`, package entries in `Cargo.lock`, and `CITATION.cff` must agree.
Update the citation release date and move the completed changelog entries to
an explicitly dated version section, leaving `Unreleased` for subsequent work.
Catalog additions also require consistent material manifests, source audits,
refit ledgers, discovery counts, and dataset-resource hashes.

Regenerate the aggregate source audit without rewriting curated records, then
refresh the independent refit and publication ledgers:

```bash
uv run python -m scripts.apply_primary_source_audit --aggregate-only
uv run python -m scripts.validate_primary_eos_refits
uv run python -m scripts.generate_paper_investigation_ledger
uv run python -m scripts.generate_source_audit_index
```

A refit diagnostic must preserve the source's temperature coordinates and
calibration provenance. Material-specific adapters handle Fei Ne's joint PVT
fit and NaCl's staged cold/hot reconstruction; the two NaCl cold figures contain
duplicate observations and must not be pooled. Numerical agreement remains
qualified by the source inputs and uncertainty conventions documented in each
audit.

Run the development, native, documentation and distribution checks above.
Smoke-test the installed wheel outside the checkout, and validate both wheel
and source-distribution metadata with Twine. Commit the completed changes,
push `main`, then manually run **CI** on that branch with **Run the full
compatibility matrix** enabled. Require that full run to pass on the exact
release commit before pushing an annotated `v<version>` tag. Routine CI alone
does not perform the complete pre-release checks. The tag starts
`.github/workflows/publish.yml`, which checks
version/changelog consistency, builds and smoke-tests all supported platform
wheels, publishes them to PyPI, and creates a GitHub release with the changelog
notes and distribution artifacts. Verify both publication destinations before
calling the release complete. A failed run can be restarted with the workflow's
existing-tag `version` input; published version numbers are not reused for
changed source.

## Local research artifacts

Keep downloaded publications, OCR intermediates, working archives, temporary
plots and logs under the ignored `output/` directory or an external research
archive. Ignoring an artifact does not back it up; preserve original evidence
and its checksums separately.

Version curated observations, provenance manifests, reusable replay scripts,
audit notes and small deterministic inputs required by tests. Required checks
must work from a fresh checkout without local publisher PDFs or scratch output.
Published documentation figures remain versioned alongside their reports.
