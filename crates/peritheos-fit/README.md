# peritheos-fit

`peritheos-fit` provides bounded robust nonlinear least squares, structured
errors-in-variables solving, covariance and uncertainty kernels, and typed
end-to-end EOS fitting for `peritheos-core` models.

The EOS entry points accept physical observations and model-parameter vectors;
they construct weighted residuals and any latent volume or temperature blocks
inside Rust. A typical isothermal fit uses a Rust model factory:

```rust
use peritheos_core::isothermal::BM3;
use peritheos_fit::{
    fit_isothermal_eos, IsothermalObservations, SolverOptions,
};

let pressure = vec![30.0, 15.0, 5.0];
let volume = vec![8.0, 9.0, 9.7];
let sigma = vec![0.1; pressure.len()];
let result = fit_isothermal_eos(
    IsothermalObservations {
        pressure: &pressure,
        volume: &volume,
        pressure_sigma: &sigma,
        volume_sigma: None,
        observation_cholesky: None,
    },
    &[110.0, 4.0],
    &[50.0, 1.0],
    &[250.0, 8.0],
    SolverOptions::default(),
    |parameters| {
        BM3::new(10.0, parameters[0], parameters[1])
            .map_err(|error| peritheos_fit::FitError::Evaluation(error.to_string()))
    },
)?;
# Ok::<(), peritheos_fit::FitError>(())
```

The lower-level `least_squares` entry point remains available for arbitrary
Rust residual functions. Python callbacks and SciPy are compatibility
fallbacks in the Python package, not dependencies of this crate.

`fit_joint_eos` uses the same typed observation interface for a factory whose
ordered parameter slice reconstructs both the reference and thermal model.
After fitting, `parameter_covariance` profiles global parameters from the
complete residual Jacobian, including latent state coordinates.

For forward uncertainty, `propagate_model_uncertainty` combines an adaptive
finite-difference model Jacobian with a complete parameter covariance matrix.
`monte_carlo_model_uncertainty` provides deterministic native sampling,
invalid-model rejection, confidence intervals, and optional output covariance.
Its seeded random stream is stable within the Rust implementation but is not
claimed to reproduce NumPy's random stream.
