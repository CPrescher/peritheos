# peritheos-core

`peritheos-core` provides equations of state, thermodynamic properties,
inversion, batch evaluation, and `.eosmat` loading without Python, NumPy, or an
array framework.

The Rust crates are currently consumed from the workspace rather than crates.io:

```toml
[dependencies]
peritheos-core = { path = "../peritheos/crates/peritheos-core" }
```

## Choose an entry point

| Need | Start with |
|---|---|
| Known model and parameters | `isothermal` or `thermal` model constructor |
| Pressure, volume, or temperature inversion | `IsothermalEos` or `ThermalEos` |
| Heat capacity and adiabatic properties | `CaloricEos` |
| Several states in input order | extension traits in `batch` |
| A curated or external material record | `load_eosmat` and `EosRecord` |

## Construct and evaluate a model

```rust
use peritheos_core::{isothermal::BM3, IsothermalEos};

let eos = BM3::new(10.0, 160.0, 4.0)?;
let pressure = eos.pressure(9.0)?;
let recovered_volume = eos.volume(pressure)?;
assert!((recovered_volume - 9.0).abs() < 1.0e-10);
# Ok::<(), peritheos_core::EosError>(())
```

Every constructor and calculation returns a checked `Result`. Importing the
trait is required for its methods to participate in Rust method lookup.

Run the complete example, which also evaluates bulk modulus and a batch:

```console
cargo run -p peritheos-core --example isothermal_workflow
```

## Add temperature and caloric properties

Thermal equations wrap a reference isotherm:

```rust
use peritheos_core::{
    isothermal::BM3,
    thermal::MieGruneisenDebye,
    CaloricEos, ThermalEos,
};

let reference = BM3::new(1.02, 165.0, 5.0)?;
let eos = MieGruneisenDebye::new(reference, 300.0, 170.0, 2.9, 1.0, 1.0)?;
let pressure = eos.pressure(0.95, 1_500.0)?;
let heat_capacity = eos.molar_heat_capacity_v(0.95, 1_500.0)?;
# Ok::<(), peritheos_core::EosError>(())
```

The `thermal_workflow` example demonstrates pressure, P-T-to-V inversion,
P-V-to-T inversion, `Cv`, and `Cp`:

```console
cargo run -p peritheos-core --example thermal_workflow
```

## Load an `.eosmat` material record

Use a material file when the parameter set, reference state, provenance, and
model choice should travel together:

```rust,no_run
use peritheos_core::load_eosmat;

let material = load_eosmat("gold.eosmat")?;
let record = material.default_record().expect("an EOS record");
let pressure = record.pressure(60.0, 300.0)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Canonical Peritheos format-3 and legacy Dioptas format-2 documents are
supported. The loader uses a fixed built-in model registry and preserves
unknown JSON extension fields. A self-contained format-3 example is runnable
with:

```console
cargo run -p peritheos-core --example material_record
```

## Batch behavior

`IsothermalEosBatch`, `ThermalEosBatch`, and `CaloricEosBatch` accept slices,
preserve order, and stop at the first scalar error. They deliberately do not
select a threading or array framework, leaving that decision to applications.

## Units and inversion conventions

- Pressure and bulk modulus: GPa.
- Temperature: K, finite and positive.
- Isothermal volume: any positive unit used consistently with `V0`.
- Thermal energy models: molar volume in `J bar^-1 mol^-1`.
- Molar energy and heat capacity: SI units.
- Volume and temperature inversion return the supported branch nearest the
  reference state.

See the crate-level rustdoc for a guided API tour. The repository's
`docs/rust-api.md` explains the relationship between this crate,
`peritheos-fit`, and the Python package.

This API is tested and packaged in CI but is not yet included in the Python
package's public compatibility contract.
