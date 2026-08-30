# Peritheos benchmarks

`python_baseline.py` records representative end-to-end performance through the
public Python API. The committed `backend: python` file captures Python 0.5.x
before the Rust migration; new runs identify the active native backend and
include conversion and FFI overhead. The script intentionally uses only
runtime dependencies and writes machine-readable JSON:

```bash
uv run python benchmarks/python_baseline.py \
  --output /tmp/peritheos-native-macos-arm64.json
```

Use `--quick` for a smoke run. Full results are machine-specific evidence, not
release performance guarantees. Compare backends on the same machine, Python
environment, workload sizes, and power state. The benchmark covers analytical
array evaluation, scalar calls, volume inversion, costly thermal quadrature,
ordinary and errors-in-variables fitting, and linear and Monte Carlo
uncertainty propagation.

Future Rust Criterion benchmarks should use the same physical parameters and
state ranges. Python-binding benchmarks must continue to measure public Python
calls so data conversion and FFI overhead are included.
