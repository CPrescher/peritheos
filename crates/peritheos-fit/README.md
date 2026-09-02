# peritheos-fit

`peritheos-fit` adds bounded robust nonlinear least squares,
errors-in-variables fitting, covariance estimation, and forward uncertainty
propagation to `peritheos-core` models.

The Rust crates are currently consumed from the workspace rather than crates.io:

```toml
[dependencies]
peritheos-core = { path = "../peritheos/crates/peritheos-core" }
peritheos-fit = { path = "../peritheos/crates/peritheos-fit" }
```

## Choose a workflow

| Data or task | Entry point |
|---|---|
| P-V data and an isothermal model | `fit_isothermal_eos` |
| P-V-T data and a thermal model | `fit_thermal_eos` |
| Reference and thermal parameters together | `fit_joint_eos` |
| Custom residual function | `least_squares` |
| Parameter covariance from a fit | `parameter_covariance` |
| Local linear output uncertainty | `propagate_model_uncertainty` |
| Nonlinear sampled uncertainty | `monte_carlo_model_uncertainty` |

## Fit an EOS

High-level fitting accepts observations, parameter vectors, bounds, solver
options, and a model factory. The order of the factory parameter slice defines
the order everywhere else:

```rust
use peritheos_core::isothermal::BM3;
use peritheos_fit::{
    fit_isothermal_eos, FitError, IsothermalObservations, SolverOptions,
};

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
# Ok::<(), FitError>(())
```

To account for volume error, provide one positive value per observation in
`volume_sigma`; adjusted volumes then become latent solver parameters. Thermal
observations independently support latent volumes and temperatures. Per-point
lower Cholesky factors represent correlated P-V or P-V-T observations.

Run a synthetic recovery, robust loss, and covariance calculation:

```console
cargo run -p peritheos-fit --example fit_isothermal
```

The joint P-V-T example reconstructs both a BM3 reference isotherm and its
thermal-pressure component from each parameter vector:

```console
cargo run -p peritheos-fit --example fit_thermal
```

## Interpret a result

`EosFitResult` contains final predicted pressures plus `SolverResult`:

- `parameters` follows factory order; latent coordinates follow model
  parameters when errors-in-variables fitting is enabled.
- `residuals` are already whitened by the supplied uncertainties.
- `jacobian` is row-major and matches the residual and parameter dimensions.
- `cost`, `optimality`, evaluation counts, status, and message expose solver
  diagnostics.

Use `parameter_covariance` with the complete Jacobian and the count of leading
model parameters. Latent variables are profiled out automatically.

## Propagate parameter uncertainty

`propagate_model_uncertainty` computes an adaptive finite-difference model
Jacobian and applies the delta method. Covariance matrices are flat row-major
slices:

```console
cargo run -p peritheos-fit --example propagate_uncertainty
```

Set `full_covariance` when correlations between output states matter. Pass
independently calculated state-variable variance for every output, or zeros
when propagating parameters only. For strongly nonlinear problems, use seeded
`monte_carlo_model_uncertainty`; it reports rejected samples as well as
confidence intervals.

## Failure handling and robust losses

Model factories should convert `EosError` to `FitError::Evaluation`. Invalid
array dimensions, uncertainties, bounds, or covariance matrices return
`FitError::InvalidInput`. Configure `Linear`, `SoftL1`, `Huber`, `Cauchy`, or
`Arctan` loss through `SolverOptions`.

The crate owns its numerical conventions and does not expose a third-party
solver type. Python callbacks and SciPy remain compatibility fallbacks in the
Python package, not dependencies of this crate.
