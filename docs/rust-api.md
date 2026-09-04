# Rust API guide

Peritheos provides one reusable Rust crate. EOS models and thermodynamic
traits live at the crate root, while fitting and uncertainty live in the
`peritheos::fit` module. This guide focuses on complete user workflows; the
generated rustdoc remains the source for exact types and signatures.

!!! note "Release status"
    The crate is tested as a publishable package but is not yet released on
    crates.io. Depend on it by Git revision or workspace path and pin the
    revision. The Python API's stability contract does not yet cover the Rust
    API.

## Add the crate

From a checkout next to your application:

```toml
[dependencies]
peritheos = { path = "../peritheos/crates/peritheos" }
```

The minimum supported Rust version is 1.83. No feature flag is needed for
fitting or uncertainty work.

## Decide how to represent an EOS

Use a concrete model when its type and parameters are part of your program:

```rust
use peritheos::{isothermal::BM3, IsothermalEos};

let eos = BM3::new(10.0, 160.0, 4.0)?;
let pressure = eos.pressure(9.0)?;
let volume = eos.volume(pressure)?;
# Ok::<(), peritheos::EosError>(())
```

Use an `.eosmat` record when model choice, parameters, reference state,
provenance, and extension metadata should be data:

```rust,no_run
use peritheos::{load_eosmat, load_eosmat_str};

let mut material = load_eosmat("gold.eosmat")?;
let record = material
    .default_equilibrium_record()
    .expect("a default equilibrium EOS record");
let pressure = record.pressure(60.0, 300.0)?;

// Unknown extension fields survive document edits and serialization.
material.document["application_note"] = "calibration copy".into();
material.validate()?;
let serialized = material.to_json()?;
let round_tripped = load_eosmat_str(&serialized)?;
round_tripped.save("gold-copy.eosmat")?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`EosRecord` uses runtime dispatch over the built-in registry. A concrete model
uses static dispatch and gives access to model-specific methods. For decoded
`serde_json::Value` documents, use `validate_eosmat_document`,
`serialize_eosmat`, and `save_eosmat` directly. Serialization validates every
record before writing and retains fields that the Rust model registry does not
interpret.

For a shock path, use `material.default_hugoniot_record()` or iterate
`material.hugoniot_records()`. These return typed views exposing loading path,
branch kind, precursor state, mass basis, branch domain, and
`state_from_particle_velocity()`. The generic `default_record()` remains for
compatibility but always prefers an equilibrium record.

## Work with the common traits

Rust trait methods must be brought into scope:

| Trait | Main operations |
|---|---|
| `IsothermalEos` | pressure, bulk modulus, `dK/dP`, volume inversion |
| `hugoniot::Hugoniot` | constrained-path pressure, volume, velocities, density, and energy |
| `ThermalEos` | thermal/total pressure, P-T-to-V and P-V-to-T inversion, expansivity, compressibility |
| `CaloricEos` | `Cv`, `Cp`, Gruneisen parameter, adiabatic bulk modulus |
| batch extension traits | ordered slice versions of the scalar operations |

Thermal models own their reference isotherm:

```rust
use peritheos::{
    isothermal::BM3,
    thermal::MieGruneisenDebye,
    CaloricEos, ThermalEos,
};

let reference = BM3::new(1.02, 165.0, 5.0)?;
let eos = MieGruneisenDebye::new(reference, 300.0, 170.0, 2.9, 1.0, 1.0)?;
let pressure = eos.pressure(0.95, 1_500.0)?;
let cv = eos.molar_heat_capacity_v(0.95, 1_500.0)?;
# Ok::<(), peritheos::EosError>(())
```

## Fit P-V observations

Fitting uses a closure that reconstructs the model from an ordered parameter
slice. This example fixes `V0` and estimates `K0` and `K0_prime`:

```rust
use peritheos::fit::{fit_isothermal_eos, FitError, IsothermalObservations, SolverOptions};
use peritheos::isothermal::BM3;

let pressure = [39.1, 25.3, 15.2, 7.7, 2.2];
let volume = [8.2, 8.6, 9.0, 9.4, 9.8];
let pressure_sigma = [0.2; 5];
let fit = fit_isothermal_eos(
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
    |p| BM3::new(10.0, p[0], p[1])
        .map_err(FitError::from),
)?;
println!("K0 = {:.3} ± {:.3} GPa", fit.parameters[0], fit.standard_errors[0]);
println!("chi-square = {:.3}", fit.chi_square);
# Ok::<(), FitError>(())
```

Provide `volume_sigma` to use latent adjusted volumes. For P-V-T data,
`ThermalObservations` also supports `temperature_sigma`. Lower-triangular
per-observation Cholesky factors replace independent standard deviations for
correlated measurements. `adjusted_volume` and `adjusted_temperature` expose
the fitted latent states without requiring slices into the solver vector.

## Carry a fit into predictions

`EosFitResult` directly contains the profiled parameter covariance, standard
errors, correlation matrix, and chi-square/AIC/BIC statistics. Its covariance
is unscaled `(J^T J)^+` in flat row-major order; multiply by the reduced
chi-square when estimating variance from residual scatter. Use that matrix
with `propagate_model_uncertainty` for local linear propagation, or
`monte_carlo_model_uncertainty` when invalid sampled states and nonlinear
outputs need to be represented. Both APIs use flat row-major covariance
matrices and return output uncertainty in the evaluator's output order.

## Units and errors

- Pressure and bulk modulus use GPa; temperature uses kelvin.
- Isothermal volume is unit-agnostic when used consistently with `V0`.
- Energy-based thermal equations require molar volume in
  `J bar^-1 mol^-1` and return molar SI energies and heat capacities.
- Model construction and evaluation return `EosError`; fitting returns
  `FitError`. Both provide stable `kind()` and `code()` accessors. Mapping an
  `EosError` with `FitError::from` retains it as the error source.
- Inversion follows the supported branch nearest the model reference state.

See [Error handling](error-handling.md) for the full taxonomy, source chains,
and matching examples.

## Run the examples and build rustdoc

The examples are executable specifications and are compiled in CI:

```console
cargo run -p peritheos --example isothermal_workflow
cargo run -p peritheos --example thermal_workflow
cargo run -p peritheos --example material_record
cargo run -p peritheos --example fit_isothermal
cargo run -p peritheos --example fit_thermal
cargo run -p peritheos --example propagate_uncertainty
cargo doc --workspace --all-features --no-deps --open
```

Start on the `peritheos` landing page in the generated documentation. It
contains runnable quick starts and links to the model and fitting modules.
