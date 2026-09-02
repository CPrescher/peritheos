# Rust API guide

Peritheos provides two reusable Rust crates: `peritheos-core` for equations of
state and `peritheos-fit` for fitting and uncertainty. This guide focuses on
complete user workflows; the generated rustdoc remains the source for exact
types and signatures.

!!! note "Release status"
    The crates are tested as publishable packages but are not yet released on
    crates.io. Depend on them by Git revision or workspace path and pin the
    revision. The Python API's stability contract does not yet cover the Rust
    API.

## Add the crates

From a checkout next to your application:

```toml
[dependencies]
peritheos-core = { path = "../peritheos/crates/peritheos-core" }
peritheos-fit = { path = "../peritheos/crates/peritheos-fit" }
```

The minimum supported Rust version is 1.83. Applications that only evaluate
models need `peritheos-core`; add `peritheos-fit` only for estimation or
uncertainty work.

## Decide how to represent an EOS

Use a concrete model when its type and parameters are part of your program:

```rust
use peritheos_core::{isothermal::BM3, IsothermalEos};

let eos = BM3::new(10.0, 160.0, 4.0)?;
let pressure = eos.pressure(9.0)?;
let volume = eos.volume(pressure)?;
# Ok::<(), peritheos_core::EosError>(())
```

Use an `.eosmat` record when model choice, parameters, reference state,
provenance, and extension metadata should be data:

```rust,no_run
use peritheos_core::load_eosmat;

let material = load_eosmat("gold.eosmat")?;
let record = material.default_record().expect("a default EOS record");
let pressure = record.pressure(60.0, 300.0)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`EosRecord` uses runtime dispatch over the built-in registry. A concrete model
uses static dispatch and gives access to model-specific methods.

## Work with the common traits

Rust trait methods must be brought into scope:

| Trait | Main operations |
|---|---|
| `IsothermalEos` | pressure, bulk modulus, `dK/dP`, volume inversion |
| `ThermalEos` | thermal/total pressure, P-T-to-V and P-V-to-T inversion, expansivity, compressibility |
| `CaloricEos` | `Cv`, `Cp`, Gruneisen parameter, adiabatic bulk modulus |
| batch extension traits | ordered slice versions of the scalar operations |

Thermal models own their reference isotherm:

```rust
use peritheos_core::{
    isothermal::BM3,
    thermal::MieGruneisenDebye,
    CaloricEos, ThermalEos,
};

let reference = BM3::new(1.02, 165.0, 5.0)?;
let eos = MieGruneisenDebye::new(reference, 300.0, 170.0, 2.9, 1.0, 1.0)?;
let pressure = eos.pressure(0.95, 1_500.0)?;
let cv = eos.molar_heat_capacity_v(0.95, 1_500.0)?;
# Ok::<(), peritheos_core::EosError>(())
```

## Fit P-V observations

Fitting uses a closure that reconstructs the model from an ordered parameter
slice. This example fixes `V0` and estimates `K0` and `K0_prime`:

```rust
use peritheos_core::isothermal::BM3;
use peritheos_fit::{fit_isothermal_eos, FitError, IsothermalObservations, SolverOptions};

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
        .map_err(|error| FitError::Evaluation(error.to_string())),
)?;
println!("K0 = {:.3} GPa", fit.solver.parameters[0]);
# Ok::<(), FitError>(())
```

Provide `volume_sigma` to use latent adjusted volumes. For P-V-T data,
`ThermalObservations` also supports `temperature_sigma`. Lower-triangular
per-observation Cholesky factors replace independent standard deviations for
correlated measurements.

## Carry a fit into predictions

The final residual Jacobian can be converted to a parameter covariance with
`parameter_covariance`. Use that matrix with
`propagate_model_uncertainty` for local linear propagation, or
`monte_carlo_model_uncertainty` when invalid sampled states and nonlinear
outputs need to be represented. Both APIs use flat row-major covariance
matrices and return output uncertainty in the evaluator's output order.

## Units and errors

- Pressure and bulk modulus use GPa; temperature uses kelvin.
- Isothermal volume is unit-agnostic when used consistently with `V0`.
- Energy-based thermal equations require molar volume in
  `J bar^-1 mol^-1` and return molar SI energies and heat capacities.
- Model construction and evaluation return `EosError`; fitting returns
  `FitError`. Do not unwrap user-provided parameters or states.
- Inversion follows the supported branch nearest the model reference state.

## Run the examples and build rustdoc

The examples are executable specifications and are compiled in CI:

```console
cargo run -p peritheos-core --example isothermal_workflow
cargo run -p peritheos-core --example thermal_workflow
cargo run -p peritheos-core --example material_record
cargo run -p peritheos-fit --example fit_isothermal
cargo run -p peritheos-fit --example fit_thermal
cargo run -p peritheos-fit --example propagate_uncertainty
cargo doc --workspace --all-features --no-deps --open
```

Start in each crate's landing page in the generated documentation. It contains
runnable quick starts and links to every built-in model family.
