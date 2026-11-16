# Rust API guide

Peritheos provides one reusable Rust crate. EOS models and thermodynamic
traits live at the crate root, while fitting and uncertainty live in the
`peritheos::fit` module. This guide focuses on complete user workflows; the
generated rustdoc remains the source for exact types and signatures.

!!! note "Stability"
    The unified crate is published starting with Peritheos 0.7. The Python
    API's stability contract does not yet cover the Rust API, so applications
    should pin a compatible minor release.

## Add the crate

For the 0.9 release series:

```toml
[dependencies]
peritheos = "0.9"
```

When developing against a local checkout, replace the version with
`{ path = "../peritheos/crates/peritheos" }`.

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

For a second-order Birch–Murnaghan law defined at finite pressure, use
`isothermal::BM2Reference::new(v_ref, k_ref, p_ref)`. It implements
`P(V) = p_ref + 1.5*k_ref*((v_ref/V)^(7/3) - (v_ref/V)^(5/3))`;
its reference volume and modulus belong to `p_ref`. Ordinary `BM2` keeps its
zero-pressure convention. EOSMAT uses the explicit model identifier
`birch_murnaghan_2_reference`, requiring a reader with that capability.

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

## Native EOSMAT v4 files

`load_eosmat`, `load_eosmat_str` and `material_from_value` accept native v4,
including file-local shared metadata, as well as v3/v2. The original native
`Material.document` is retained by `to_json()` and `save()`. Unsupported models,
dangling references, cycles, metadata overrides and invalid reconstruction
assessments fail validation; there is no fallback equation or inferred parity.

```rust,no_run
use peritheos::{load_eosmat, export_eosmat_v3, save_eosmat};

let material = load_eosmat("my-material.eosmat")?;
material.save("native-copy.eosmat")?;
let legacy = export_eosmat_v3(&material.document)?;
save_eosmat("v3-copy.eosmat", &legacy)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

V3 export expands shared metadata while retaining unused supporting definitions.
The [authoring guide](eosmat-authoring.md) includes complete basic, thermal and
database-style v4 files and installed creation/validation/formatting commands.
These teaching files are diagnostic and make no scientific validation claim.

## EOSMAT v5 check headlines and audit sidecars

A v5 record may carry a short `check` block, and bundled materials keep each
record's evidence in an audit sidecar (`../audits/<id>.audit.json`) named by
`check.audit`. The [EOSMAT v5 reference](eosmat-v5.md#record-check-headline-and-audit-sidecars)
defines both. The Rust reader follows the Python reader:

| Input | Sidecars | Result |
| --- | --- | --- |
| `load_eosmat(path)` | Resolved relative to the file and merged | Complete inline document; a missing sidecar is `EosmatError::MissingAudit` (`code() == "eosmat.audit_missing"`) |
| `load_eosmat_without_audit(path)` | Not read | The record file as authored |
| `load_eosmat_str`, `load_eosmat_reader`, `material_from_value` | Not read (no location to resolve against) | The document as given |
| `load_eosmat_str_with_audit(record, sidecar)` | The given sidecar | Complete inline document |

`merge_eosmat_audit(&document, &sidecar)` and
`merge_referenced_eosmat_audits(&document, directory)` return the merged
`serde_json::Value` for callers that decode documents themselves. Merging
regenerates the assessment `report_url`, `report_path` and `report` links that
bundled sidecars omit, and rejects a sidecar for another material, entries for
unknown records and blocks held both inline and in the sidecar.

Every `EosRecord` reports its headline whether or not the sidecar was merged:

```rust,no_run
use peritheos::{load_eosmat, load_eosmat_without_audit, Reproduction};

let material = load_eosmat("peritheos/data/materials/mgo.eosmat")?;
for record in &material.eos_records {
    let tested = record.reproduction() != Reproduction::Untested;
    println!(
        "{}: {} ({}) {:?}",
        record.identifier,
        record.reproduction(),
        record.compared(),
        record.warning()
    );
    // Complete evidence, pressure_calibration and other audit-only blocks.
    let audit = record.audit().expect("load_eosmat merged the sidecar");
    assert!(audit.contains_key("evidence") || !tested);
}

let record_file = load_eosmat_without_audit("peritheos/data/materials/mgo.eosmat")?;
let record = &record_file.eos_records[0];
assert!(record.audit_is_external() && record.audit().is_none());
assert_eq!(record.reproduction(), material.eos_records[0].reproduction());
# Ok::<(), Box<dyn std::error::Error>>(())
```

Each `check` field has an accessor; the [EOS labels](eos-labels.md) define the
values:

| Accessor | `check` field | Returns |
| --- | --- | --- |
| `reproduction()` | `reproduction` | `Reproduction::{Reproduced, Close, Disagrees, Untested}` |
| `compared()` | `compared` | `Compared::{PublishedFit, PeritheosRefit, SourceEquation, Nothing}` |
| `why()` | `why` | The generated reason for `close` or `disagrees`, or `None` |
| `warning()` | `warning` | The stored usage restriction, or `None` |
| `source_check()` | `source_check` | `SourceCheck::{Verified, Pending, Unresolved, Deferred}`, or `None` |
| `source_note()` | `source_note` | The primary-source finding, or `None` |

`as_str()` gives the EOSMAT spelling and `SourceCheck::status_code()` the
`scientific_validation.status` code it renames. A stored `check` value wins;
without one, every field is derived from the record's inline evidence exactly
as the bundled blocks were written, `why` included
(`derive_eosmat_check(&record_document, Some(&material_document))`). A record
without an assessment is `Untested` with `Compared::Nothing`, never a positive
result. `check()` and `audit_reference()` expose the stored block and its
sidecar link.

When `check` and inline evidence are both present, `check.reproduction`,
`check.compared` and `check.source_check` must agree with them or validation
fails. The earlier names `check.finding` and `check.source` are rejected with
an error naming `source_note` and `source_check`. A record file whose
audit is not merged is executable and gives the same pressures as the merged
document. Checks that need audit-only blocks, such as a refit record's
`fit_provenance` or the EOS `pressure_calibration` used by a dataset pressure
reduction, are applied when the complete document is validated. Without its
evidence, such a record reports the `scientific_validation.status` code that
`check.source_check` renames (or `pending_primary_source_check`) in a v3 export.

A v5 material may name its defaults once in `default_records`
(`{"equilibrium": "<record id>", "hugoniot": "<record id>"}`), as bundled
files do; per-record `default`/`default_for` flags remain valid without it.
`EosRecord::is_default` and the default-record methods read either form. A
map naming an unknown record or a record of the other category, an unknown
key, or a record flag that contradicts the map is a validation error.

`export_eosmat_v3` of a merged document carries the full evidence inline and
omits `check`; for a v5 source it writes `default` on every record (and
`default_for` on each designated default) in place of `default_records`; a record whose sidecar was not merged keeps `check` so the
audit link is not lost. `Material::to_json()` and `save()` write the merged
document as one self-contained file with inline evidence. v3 and v4 files are
read unchanged.

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

For analytic derivatives, use `fit_isothermal_eos_with_derivatives` or
`fit_thermal_eos_with_derivatives` (also for joint reference/thermal fits).
The original factory-only functions retain finite differences. A factory can
arbitrarily transform its parameters, so its derivative callback must apply
that same chain rule and return `PressureDerivatives` in factory order.
The helpers `bm3_pressure_derivatives` and
`q_compromise_bm3_pressure_derivatives` support canonical parameter names in
any subset/order. For the BM3 example above, replace the function name and add
this callback after the factory:

```rust
|model, volume| bm3_pressure_derivatives(model, volume, &["K0", "K0_prime"])
```

Import both helpers/functions from `peritheos::fit`. The complete
`fit_isothermal` and `fit_q_compromise` crate examples demonstrate the analytic
APIs, including a joint q-compromise fit with uncertain volume and temperature.
The latter requires volumes in J/bar/mol and pressures in GPa, as for the
thermal model itself. These APIs preserve residual weighting, Cholesky
whitening, robust loss and profiled covariance. Low-level custom residual
models can use `least_squares_with_jacobian` or
`least_squares_structured_with_jacobian`; their callbacks return flat row-major
Jacobians of the weighted residuals, before robust-loss transformation.

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
cargo run -p peritheos --example fit_q_compromise
cargo run -p peritheos --example propagate_uncertainty
cargo doc --workspace --all-features --no-deps --open
```

Start on the `peritheos` landing page in the generated documentation. It
contains runnable quick starts and links to the model and fitting modules.

## Double-Debye reference models

`DoubleDebyeHelmholtz<R = Vinet>` and
`DoubleDebyeLogMomentHelmholtz<R = Vinet>` accept `Vinet`, `BM2`, `BM3`, and
`BM4` references through `isothermal::ReferenceEnergyEos`. The reference
energy is consistent with pressure for both the absolute cold-curve mode and
the reference-isotherm mode. See the
[equation reference](equations/thermal/double-debye.md) for the
energy definitions and reference-temperature conventions.

Existing concrete Vinet models retain their default type parameter. Since
0.9.0, the two corresponding `eosmat::ThermalModel` variants contain
`DoubleDebyeHelmholtz<IsothermalModel>` and
`DoubleDebyeLogMomentHelmholtz<IsothermalModel>`. Code that constructs or
destructures those variants directly must handle `IsothermalModel` references
instead of assuming a concrete `Vinet`.

## Material family metadata

`Material::family_id()` returns `Option<&str>` from the retained `.eosmat`
document. The native reader validates the optional identifier and preserves it
on export, including unknown external families. It does not require a bundled
family registry or change EOS selection. Family registry lookup and grouped
catalog discovery are currently Python APIs.

## Tabulated pressure corrections

`thermal::DebyeTabulatedThermalPressure` composes an integrated Debye model
with a bounded, volume-independent electronic-pressure table. Its pressure
correction subtracts the electronic value at the reference temperature once.
`thermal::AsymptoticDebyeTabulatedPressure` instead adds absolute phonon,
electronic and derived residual pressure to a 0 K cold curve, with separate
cold and ambient reference volumes and explicit numerical bounds. Both models
load from their canonical `.eosmat` records and retain configuration tables
through serialization. The Yokoo reconstructions remain qualified,
explicit-selection Python diagnostics; their Rust executability does not
change their scientific validation status.

The tables use linear interpolation and reject temperature extrapolation.
The absolute model also bounds its volume inversion to the reconstruction
interval and accepts 0 K. These models do not implement `CaloricEos`: a pressure
table does not establish electronic energy, entropy or heat capacity.

As of 0.12.0, the runtime `ThermalModel` and `LoadedEos` enums own table data
and implement `Clone` without `Copy`. Borrow an existing model or use
`.clone()` where an independent owned value is required.
