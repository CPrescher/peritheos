# Development

## Test suite

```bash
uv sync --all-groups
uv run ruff check .
uv run ruff format --check .
uv run pytest -q -W error --cov --cov-report=term-missing
```

If the default uv cache is unavailable in a sandbox:

```bash
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run pytest -q
```

## Documentation

```bash
uv run --group docs mkdocs build --strict
```

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
