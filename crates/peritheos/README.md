# peritheos

`peritheos` is the complete Rust API for equations of state, thermodynamic
properties, inversion, fitting, uncertainty propagation, batch evaluation,
and `.eosmat` loading. It has no Python, NumPy, or array-framework dependency.

The crate is currently consumed from the workspace rather than crates.io:

```toml
[dependencies]
peritheos = { path = "../peritheos/crates/peritheos" }
```

## Choose an entry point

| Need | Start with |
|---|---|
| Known model and parameters | `isothermal` or `thermal` model constructor |
| Pressure, volume, or temperature inversion | `IsothermalEos` or `ThermalEos` |
| Heat capacity and adiabatic properties | `CaloricEos` |
| Several states in input order | extension traits in `batch` |
| Fit P-V or P-V-T observations | high-level routines in `fit` |
| Covariance or prediction uncertainty | propagation routines in `fit` |
| A curated or external material record | `load_eosmat`, `Material`, and `EosRecord` |

## Construct and evaluate a model

```rust
use peritheos::{isothermal::BM3, IsothermalEos};

let eos = BM3::new(10.0, 160.0, 4.0)?;
let pressure = eos.pressure(9.0)?;
let recovered_volume = eos.volume(pressure)?;
assert!((recovered_volume - 9.0).abs() < 1.0e-10);
# Ok::<(), peritheos::EosError>(())
```

Every constructor and calculation returns a checked `Result`. Import the
relevant trait so that its methods participate in Rust method lookup.

## Fit observations

Fitting is a normal module of the same crate. A model factory maps each
ordered parameter slice to a validated EOS:

```rust
use peritheos::fit::{
    fit_isothermal_eos, FitError, IsothermalObservations, SolverOptions,
};
use peritheos::isothermal::BM3;

let pressure = [39.1, 25.3, 15.2, 7.7, 2.2];
let volume = [8.2, 8.6, 9.0, 9.4, 9.8];
let pressure_sigma = [0.2; 5];
let result = fit_isothermal_eos(
    IsothermalObservations {
        pressure: &pressure,
        volume: &volume,
        pressure_sigma: &pressure_sigma,
        volume_sigma: None,
        observation_cholesky: None,
    },
    &[140.0, 4.0],
    &[50.0, 1.0],
    &[300.0, 10.0],
    SolverOptions::default(),
    |parameters| {
        BM3::new(10.0, parameters[0], parameters[1])
            .map_err(|error| FitError::Evaluation(error.to_string()))
    },
)?;
assert_eq!(result.parameters.len(), 2);
assert_eq!(result.covariance.len(), 4);
assert_eq!(result.standard_errors.len(), 2);
# Ok::<(), FitError>(())
```

`fit` also provides joint thermal fitting, robust losses, covariance
estimation, delta-method propagation, and deterministic Monte Carlo sampling.
`EosFitResult` separates model parameters from latent states and directly
contains profiled covariance, standard errors, correlation, and fit
statistics. Its module documentation explains the expected matrix layouts.

## Round-trip `.eosmat` documents

`load_eosmat` validates and constructs every built-in EOS while retaining the
original JSON document, including unknown extension fields. After editing
`material.document`, call `material.validate()`, `material.to_json()`, or
`material.save(path)`. The free functions `validate_eosmat_document`,
`serialize_eosmat`, and `save_eosmat` provide the same workflow for a decoded
`serde_json::Value`.

## Run complete examples

```console
cargo run -p peritheos --example isothermal_workflow
cargo run -p peritheos --example thermal_workflow
cargo run -p peritheos --example material_record
cargo run -p peritheos --example fit_isothermal
cargo run -p peritheos --example fit_thermal
cargo run -p peritheos --example propagate_uncertainty
```

## Units and behavior

- Pressure and bulk modulus: GPa.
- Temperature: K, finite and positive.
- Isothermal volume: any positive unit used consistently with `V0`.
- Thermal energy models: molar volume in `J bar^-1 mol^-1`.
- Molar energy and heat capacity: SI units.
- Batch operations preserve order and stop at the first scalar error.
- Inversion returns the supported branch nearest the reference state.

See the crate-level rustdoc for the guided API tour and the repository's
`docs/rust-api.md` for complete workflows. This API is tested and packaged in
CI but is not yet part of the Python package's public compatibility contract.
