# Error handling

Peritheos exposes one error contract across its Python implementation, native
extension, and Rust crate. Code should make decisions from an exception type
and, when necessary, its stable error code. Error messages are written for
people and may become more specific without a breaking release.

## Python exception hierarchy

All deliberate library errors derive from `PeritheosError`. They also retain
the built-in categories used before the hierarchy was introduced, so existing
`except ValueError`, `TypeError`, `ArithmeticError`, `RuntimeError`,
`NotImplementedError`, and `KeyError` handlers continue to work.

| Exception | Also behaves as | Meaning |
|---|---|---|
| `ValidationError` | `ValueError` | invalid values, units, arrays, shapes, or state |
| `NumericalError` | `ArithmeticError` | convergence, finiteness, or numerical evaluation failure |
| `ConfigurationError` | `TypeError` | incompatible model or object configuration |
| `UnsupportedOperationError` | `NotImplementedError` | the selected model does not support an operation |
| `EosValidationError` | `EosError`, `ValidationError` | invalid EOS parameter or thermodynamic state |
| `EosNumericalError` | `EosError`, `NumericalError` | EOS evaluation or inversion failed numerically |
| `FitValidationError` | `FitError`, `ValidationError` | invalid observations, bounds, dimensions, or options |
| `FitNumericalError` | `FitError`, `NumericalError` | singular or unstable fit calculation |
| `FitEosValidationError` | `FitValidationError`, `EosValidationError` | an EOS used by a fit rejected a parameter or state |
| `FitEosNumericalError` | `FitNumericalError`, `EosNumericalError` | an EOS used by a fit failed numerically |
| `EosmatError` | `ValidationError` | invalid or unsupported material document |
| `MaterialError` | `ValidationError` | invalid material or EOS record |
| `MaterialLookupError` | `KeyError` | unknown material or record identifier |

The classes are available from `peritheos` and `peritheos.errors`. Every
instance provides:

- `code`: stable machine-readable identifier;
- `operation`: operation family when one is known;
- `field`: public parameter or state-variable name when one is known; and
- `context`: read-only additional values that can help diagnostics.

The Rust extension constructs these same Python classes; there is no separate
native exception hierarchy.

```python
from peritheos import EosValidationError
from peritheos.eos.rt import BM3

try:
    eos = BM3(V0=-1.0, K0=160.0, K0_prime=4.0)
except EosValidationError as error:
    # Prefer the type for broad handling and the code for a specific case.
    if error.code in {"eos.invalid_input", "eos.invalid_parameter"}:
        print(f"invalid {error.field or 'EOS input'}: {error}")
```

External failures retain their native meaning. Filesystem failures from
`.eosmat` load/save operations remain `OSError` subclasses, for example, so a
caller can distinguish permissions or a missing path. Malformed JSON and
document validation are reported as `EosmatError`, with the decoder retained
as `__cause__` where applicable.

## Stable error codes

General pure-Python validation uses the domain default, such as
`eos.invalid_input`, `fit.invalid_input`, or `eosmat.invalid_document`. Native
EOS and Rust APIs provide the more specific codes below.

| Code | Category |
|---|---|
| `eos.invalid_parameter` | model constructor parameter is invalid |
| `eos.invalid_state` | requested volume, pressure, or temperature is invalid |
| `eos.outside_invertible_range` | state is outside the supported inverse branch |
| `eos.bracketing_failed` | a physical root bracket could not be constructed |
| `eos.convergence_failed` | numerical solver did not converge |
| `eos.non_finite_result` | evaluation produced a non-finite value |
| `fit.invalid_input` | fit data, dimensions, bounds, or options are invalid |
| `fit.evaluation_failed` | model factory or callback failed |
| `fit.singular_system` | a required linear system is singular |
| `eosmat.io` | Rust filesystem operation failed |
| `eosmat.json` | JSON decoding or typed deserialization failed |
| `eosmat.invalid_document` | document-level validation failed |
| `eosmat.invalid_record` | an EOS record could not be constructed |

New codes may be added. Existing code should include a broader typed fallback
instead of assuming this table is exhaustive.

## Rust errors and sources

Rust exposes `EosError`, `fit::FitError`, and `EosmatError`. Each has a
non-exhaustive `*ErrorKind` enum plus `kind()` and `code()` accessors. Match
with a wildcard so new categories remain source-compatible.

```rust
use peritheos::{isothermal::BM3, EosErrorKind};

let error = BM3::new(-1.0, 160.0, 4.0).unwrap_err();
match error.kind() {
    EosErrorKind::InvalidParameter => {
        eprintln!("{} in {:?}: {error}", error.code(), error.field());
    }
    _ => eprintln!("{}: {error}", error.code()),
}
```

When a fit factory constructs or evaluates an EOS, use `FitError::from`
instead of converting the EOS error to a string. This preserves the stable EOS
code and makes the original error available through `std::error::Error::source`:

```rust
use peritheos::{fit::FitError, isothermal::BM3};

let factory = |parameters: &[f64]| {
    BM3::new(10.0, parameters[0], parameters[1]).map_err(FitError::from)
};
# let _ = factory(&[160.0, 4.0]);
```

`EosmatError` likewise retains `std::io::Error` and `serde_json::Error` as
sources. `record_identifier()` gives the affected record for
`eosmat.invalid_record`.

## Compatibility rules

The following are public compatibility commitments:

- exported exception classes and their documented built-in base classes;
- documented error-code meanings;
- Rust error variants, kind variants, and source chains; and
- the single-message value in Python `Exception.args`.

Punctuation, capitalization, solver wording, and complete message text are
not stable. Additional context fields, codes, and non-exhaustive Rust variants
may be added. Tests and applications should therefore avoid matching complete
messages.
