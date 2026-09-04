# API reference

## Errors

All deliberate Python errors derive from `PeritheosError`; domain classes such
as `EosError` and `FitError` support precise handling, while compatibility with
the previous built-in `ValueError`, `ArithmeticError`, `TypeError`, and related
categories is retained. Instances expose `code`, `operation`, `field`, and
read-only `context` attributes. The native extension raises the same classes.

```python
from peritheos import EosError, EosValidationError, PeritheosError
```

See [Error handling](error-handling.md) for the complete hierarchy, stable
codes, Rust error kinds, and source-chain examples.

## Materials and EOS records

```python
from peritheos import (
    EOSRecord,
    HugoniotBranchDomain,
    HugoniotInitialState,
    HugoniotRecord,
    HugoniotVolumeBasis,
    Material,
    get_eos_record,
    get_material,
    list_eos_records,
    list_materials,
)
from peritheos.materials import DEFERRED_EOS_RECORDS
```

`get_material(identifier)` returns a material phase and
`list_materials(formula=...)` lists or filters the curated pressure-scale
convenience catalog. This compact executable catalog is distinct from the full
115-document shared material library described below. Each
`Material` owns its `eos_records`, supports `get_eos_record(identifier)`, and
provides `to_dict()`/`from_dict()` and `to_eosmat()`/`from_eosmat()` for the
canonical executable format-3 material document. Optional crystallographic
fields are preserved but not interpreted. Loading uses a fixed model registry
and never imports an implementation path. The legacy executable snapshot-v2
reader and `to_snapshot_dict()` remain compatibility-only APIs.

`HugoniotRecord` specializes `EOSRecord` with mandatory loading-path,
branch-kind, initial-state, operational volume-basis, and branch-domain fields.
It stays in the common
`Material.eos_records` collection; `Material.hugoniot_records` and
`Material.equilibrium_records` return convenient filtered views. Hugoniot
records expose `shock_velocity()`, `particle_velocity()`, `density()`, and
`specific_internal_energy_change()` in addition to `pressure()` and `volume()`.

## Shared `.eosmat` material library

```python
from peritheos import (
    eosmat_schema,
    get_material_document,
    list_material_documents,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
```

`list_material_documents()` returns the identifiers of all 116 bundled
materials. `get_material_document(identifier)` returns a defensive copy of one
flat format-3 `.eosmat` document, including optional structure and its raw EOS
records. `load_eosmat()` also accepts native Dioptas 0.10.0 format-2 files;
`save_eosmat()` validates and preserves optional fields. `eosmat_schema()`
returns the bundled normative JSON Schema.
The [`.eosmat` schema reference](eosmat-schema.md) documents the complete
field contract and consumer defaults.

## Pressure-calibration records

```python
from peritheos import (
    find_common_pressure_calibration_routes,
    find_pressure_calibration_path,
    get_cross_calibration_edge,
    get_pressure_calibration,
    list_cross_calibration_edges,
    list_diamond_raman_calibrations,
    list_pressure_calibrations,
    list_ruby_pressure_calibrations,
    list_ruby_xrd_bridges,
    list_xrd_pressure_standards,
    recalculate_diamond_raman_pressure,
    recalculate_eos_pressure_scale,
    recalculate_pressure_calibration_path,
    recalculate_ruby_pressure,
    recalculate_ruby_to_xrd_pressure,
    recalculate_xrd_pressure_scale,
    xrd_standard_pressure,
)
```

Six published ruby R1 calibrations and two diamond-anvil Raman-edge
calibrations are bundled with their equations, coefficients, provenance, and
validity metadata. Ruby supports wavelength, wavelength shift, and normalized
wavelength ratio; diamond Raman supports wavenumber and normalized wavenumber
ratio. Both provide analytical inversion and scale-to-scale conversion. The
explicit cross-calibration edge
registry records simultaneous experiments and jointly optimized families.
XRD-to-XRD conversion inverts the
source EOS to a virtual volume of the common standard and evaluates the target
EOS at that coordinate. Ruby-to-XRD conversion uses a ruby-linked standard EOS
as the first edge. The high-level sample-EOS API recursively discovers a path
from recorded provenance and returns every node and edge used. An edge may
also transform its temperature coordinate; the Chidester KCl-to-Dorogokupets
Pt edge converts the lower KCl effective temperature to the nominal average
Pt-foil temperature and records both edge-end values in `intermediate_states`. See
[Pressure-scale normalization](pressure-scale-normalization.md) for the valid
graph edges, complete examples, cross-calibration literature, and the
difference between a derived transformation and observation-level
re-reduction.

`find_common_pressure_calibration_routes(source_nodes)` discovers XRD-standard
endpoints reachable from every supplied sample EOS and returns the shortest
route for each source. Results are ranked by maximum and total edge count, not
by scientific preference; callers should restrict `target_nodes` to a suitable
standard or internally consistent family and verify the common validity range.

Transferred Dioptas records have completed a primary-source classification,
and native primary-sourced records have been added for aragonite BM2, the
B2-KCl P-V-T pressure calibration, and the Correa and Benedict diamond
Helmholtz models, including derived variants anchored to the experimental
Dewaele 298 K Vinet isotherm, together with a primary-sourced Campbell-Heinz
RbCl-B2 record. All 163 bundled records are
RbCl-B2 record and the Dewaele et al. (2000) MgO BM3-MGD record. All 163 bundled records are
`primary_source_validated`; none remains pending or deferred. `Material.from_eosmat()` constructs
validated records and refuses deferred ones by default; callers can inspect legacy values with
`require_primary_validation=False` and select records with
`record_identifiers=(...)`. This opt-in never changes the stored status.

The Benedict et al. diamond branch is available directly from the convenience
catalog or from the bundled `diamond.eosmat` document. Its public volume is the
eight-atom conventional-cell volume:

```python
from peritheos.materials import DIAMOND_BENEDICT_2014

volume = 8.0 * 4.6542704116  # A^3/conventional cell
pressure = DIAMOND_BENEDICT_2014.pressure(volume, 3000.0)  # about 150 GPa

# Forward phase-line state from a cold pressure, hot temperature, and
# experimentally assumed retained thermal-pressure fraction.
heated_volume = DIAMOND_BENEDICT_2014.volume_with_dac_confinement(
    60.0,
    2000.0,
    f_dac=0.25,
)
thermal_increment = DIAMOND_BENEDICT_2014.thermal_pressure_increment(
    heated_volume, 2000.0
)
confinement_increment = DIAMOND_BENEDICT_2014.dac_thermal_pressure(
    heated_volume, 2000.0, f_dac=0.25
)
total_hot_pressure = 60.0 + confinement_increment

# Both inputs remain public eight-atom conventional-cell volumes.
temperature = DIAMOND_BENEDICT_2014.temperature_from_volumes(
    ambient_volume=39.61987225,
    heated_volume=40.0,
    f_dac=0.25,
)
```

The record's `reference_volume` is the Table I **0 K motionless-ion cold-curve**
volume, not a 300 K zero-total-pressure state. The record also carries the
source's narrower DFT-MD comparison domain and diamond phase-stability caveat.

The earlier Correa logarithmic-moment branch is a separate convenience record:

```python
from peritheos.materials import DIAMOND_CORREA_2008

pressure = DIAMOND_CORREA_2008.pressure(8.0 * 4.43, 5000.0)
# 202.628115 GPa
```

It preserves the published DFT-GGA cold-curve volume and records the authors'
approximately 3% ambient-volume caveat rather than applying an implicit shift.

For example, the independently reproducible staged aragonite BM2 record can be
used at its 298 K reference state or at the represented high temperatures:

```python
from peritheos import Material, get_material_document

document = get_material_document("aragonite")
aragonite = Material.from_eosmat(
    document,
    record_identifiers=["aragonite_martinez_1996_bm2_2"],
).eos_records[0]

pressure = aragonite.pressure(volume=215.0)  # GPa at 298 K
hot_pressure = aragonite.pressure(volume=215.0, temperature=873.0)
result = aragonite.pressure_with_uncertainty(
    volume=215.0,
    volume_sigma=0.05,
)
```

The separate Martinez global thermal BM3 reduction is intentionally absent:
its published fitted reference volume is missing, and its remaining coefficients
do not reproduce the printed dataset under the documented equations.

`get_eos_record(identifier)` and `list_eos_records(formula=...)` are convenient
lookups within that curated pressure-scale set. Use
`list_material_documents()` and `Material.from_eosmat()` for the complete
shared library. Each `EOSRecord` provides:

- `pressure(volume, temperature=None, check_validity=False)` in GPa;
- `volume(pressure, temperature=None, check_validity=False)` in
  angstrom^3/conventional unit cell;
- `pressure_with_uncertainty(...)` and `volume_with_uncertainty(...)`;
- `within_calibration_range(volume, temperature=None)` (preferred) and the
  compatibility alias `within_validity(...)`;
- `reference`, `validity`, `parameter_provenance`, `notes`, and volume metadata.

See [Pressure standards](pressure-standards.md) for the subset and workflows
commonly used for pressure calibration.

## Isothermal equations of state

```python
from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    Vinet,
)
```

Common methods:

- `pressure(V)`
- `bulk_modulus(V)`
- `bulk_modulus_derivative(V)` for $dK/dP$; the base implementation is numerical
- `volume(P)` and `calculate_volume(P)`

Constructor signatures and special requirements are:

| Class | Signature after class name | Special requirement |
|---|---|---|
| `BM2` | `(V0, K0)` | $K_0'=4$ is implied |
| `BM3` | `(V0, K0, K0_prime)` | none |
| `BM4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `Murnaghan` | `(V0, K0, K0_prime)` | supports the `K0_prime=0` limit |
| `NaturalStrain2` | `(V0, K0)` | $K_0'=2$ is implied |
| `NaturalStrain3` | `(V0, K0, K0_prime)` | none |
| `NaturalStrain4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `ModifiedTait` | `(V0, K0, K0_prime, K0_double_prime)` | rejects singular coefficient sets and volumes outside its real domain |
| `Vinet` | `(V0, K0, K0_prime)` | none |
| `Holzapfel` | `(V0, K0, K0_prime, n, Z)` | molar volume in `J bar^-1 mol^-1` |

`Holzapfel` overrides `bulk_modulus_derivative(V, eps=1e-6)` with its existing
specialized numerical form.
See the [equation reference](equation-reference.md#isothermal-equations) for
the mathematical definitions and coefficient domains.

## Thermal equations of state

```python
from peritheos.eos.thermal import (
    DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz,
    HollandPowell2011,
    LinearThermalPressure,
    LogVolumeThermalPressure,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    MultiOscillatorGruneisenThermalEOS,
    SecondOrderTaylorThermalPressure,
    Tange2009Debye,
    ThermalModifiedTait,
    ThermalReferenceStateEOS,
)
```

Thermal constructor signatures are:

| Class | Parameters after `rt_eos` |
|---|---|
| `DoubleDebyeHelmholtz` | `Vp, theta_a0, a_a, b_a, theta_b0, a_b, b_b, theta_1_0, a_1, b_1`, followed by optional `n, alpha0, Ve, kappa, phi0` |
| `DoubleDebyeLogMomentHelmholtz` | `Vp, theta_a0, a_a, b_a, theta_b0, a_b, b_b, theta_0_0, a_0, b_0`, followed by optional `n, anharmonic_a, phi0` |
| `LinearThermalPressure` | `Tr, alpha_KT` |
| `SecondOrderTaylorThermalPressure` | `Tr, eta0, c0, c1, c2, c3, c4, c5` |
| `LogVolumeThermalPressure` | `Tr, alpha_KT_ref, dK_dT_V` |
| `ThermalReferenceStateEOS` | `Tr, alpha0, dK_dT, alpha1=0, thermal_expansion_law="constant", reference_volume_law="integrated_expansivity"`; volume laws also include `linear_temperature` and `berman` |
| `MieGruneisenDebye` | `Tr, theta0, gamma0, q, n, debye_temperature_law="integrated_gruneisen"` |
| `MieGruneisenEinstein` | `Tr, theta0, gamma0, q, n` |
| `ThermalModifiedTait` | `Tr, theta, alpha0, n` |
| `MultiOscillatorGruneisenThermalEOS` | `Tr, QE1o, mE1, QE2o, mE2, delta, t, a_0, m, g, e_0`, followed by optional `beta, QBo, d, mb, QB1o, d1, mb1, n` |
| `Tange2009Debye` | `Tr, theta0, gamma0, a, b, n` |

The Mie-Gruneisen, multi-oscillator, and constant linear thermal-pressure
classes accept any `EosBase` reference. `LogVolumeThermalPressure` and
`SecondOrderTaylorThermalPressure` require `V0`; thermal modified Tait requires `ModifiedTait`; and
`ThermalReferenceStateEOS` requires a reference that reconstructs through
`V0` and `K0`. Energy-based thermal classes require molar volume in
`J bar^-1 mol^-1`; `LinearThermalPressure`,
`SecondOrderTaylorThermalPressure`, and `ThermalReferenceStateEOS` inherit any
volume unit consistent with their reference EOS. The second-order Taylor model
adds absolute thermal pressure to a cold curve; its pressure need not vanish at
`Tr`, although `thermal_pressure_increment()` is zero there by definition.
The double-Debye Helmholtz classes instead require a `Vinet` object representing
the classical 0 K cold curve. Their `thermal_pressure()` is the absolute
non-cold contribution, including zero-point pressure, rather than a difference
from a reference temperature. Each additionally exposes `cold_energy()`,
`zero_point_energy()`, `ion_helmholtz_free_energy()`,
`anharmonic_helmholtz_free_energy()`, `helmholtz_free_energy()`,
`ion_pressure()`, and `anharmonic_pressure()`.
Their ordinary `temperature(P,V)` inversion and DAC
`temperature_from_volumes()` inversion are supported. For the latter,
both classes subtract their pressure on the 300 K isotherm so the
confinement term excludes zero-point and baseline thermal pressure.
`HollandPowell2011` is an alias for
`ThermalModifiedTait`. Exact equations and parameter roles are documented
under [Thermal equations](equation-reference.md#thermal-equations).
`Sokolova2016` is a compatibility alias for
`MultiOscillatorGruneisenThermalEOS`.

`parameter_values()` contains only numeric parameters that fitting and
uncertainty propagation may vary. `configuration_values()` reports fixed
non-numeric constructor choices such as `debye_temperature_law`;
`with_parameters()` preserves those choices when reconstructing an EOS.

Common methods:

- `thermal_pressure(V, T)`
- `thermal_pressure_increment(V, T)`
- `dac_thermal_pressure(V, T, f_dac)`
- `pressure(V, T)`
- `volume(P, T)`
- `volume_with_dac_confinement(P_cold, T, f_dac=...)`
- `temperature(P, V)` and `calculate_temperature(P, V)`
- `temperature_from_volumes(V_ambient, V_heated, f_dac=...)`
- `bulk_modulus(V, T)`
- `isothermal_compressibility(V, T)`
- `thermal_expansivity(V, T)`
- `molar_heat_capacity_v(V, T)` when a caloric model exists
- `molar_heat_capacity_p(V, T)` when a caloric model exists
- `adiabatic_bulk_modulus(V, T)` when a caloric model exists
- `gruneisen_parameter(V, T)` when a caloric model exists

`thermal_pressure_increment()` is the heating pressure above the reference
isotherm. It equals `thermal_pressure()` for reference-relative models and
`pressure(V, T) - pressure(V, Tr)` for either absolute double-Debye Helmholtz
class. `dac_thermal_pressure()` returns only `f_dac` times that increment.
`volume_with_dac_confinement()` solves the forward hot-volume problem from a
cold pressure and temperature; `temperature_from_volumes()` solves its
two-volume inverse. Both follow the empirical confinement model described in
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md);
they require `0 <= f_dac < 1`. In this API, `f_dac` means
`(P_hot - P_ambient) / Delta_P_thermal(V_heated, T)`, where the denominator is
the pressure increment above the reference-temperature isotherm. It is not a
fraction of the cold pressure.

Mie-Gruneisen models additionally expose `gruneisen_parameter()`,
`characteristic_temperature()`, and the vibrational thermodynamic methods
documented under [Thermoelastic properties](thermoelastic-properties.md).

## Shock Hugoniot equations of state

```python
from peritheos import HugoniotBase, HugoniotState, LinearUsUpHugoniot
```

`LinearUsUpHugoniot(V0, rho0, c0, s, P0=0)` implements a phase-specific
`Us = c0 + s * up` relation and the Rankine--Hugoniot jump conditions. See
[Shock Hugoniot equations of state](hugoniots.md) for equations, units,
fitting, EOSMAT records, and phase-branch semantics.

## Fitting

```python
from peritheos.fitting import (
    FitResult,
    HugoniotFitResult,
    fit_joint_eos,
    fit_linear_us_up,
    fit_rt_eos,
    fit_thermal_eos,
)
```

`FitResult` contains the fitted `model`, parameter and uncertainty mappings,
covariance and correlation matrices, raw and weighted residuals, chi-square,
degrees of freedom, AIC, BIC, convergence status, and solver message. It also
reports `adjusted_volume`, `adjusted_temperature`, `volume_corrections`, and
`temperature_corrections` for errors-in-variables fits. See
[Fitting P-V and P-V-T data](fitting.md) for the pressure, volume, and
temperature uncertainty API.

`FitResult.summary()` returns a compact text report. `FitResult.to_dict()` and
`FitResult.to_json()` export a versioned, JSON-safe record containing the model
identity and reconstructable parameters, fit arrays, diagnostics, and solver
configuration. Passing a path to `to_json()` also writes the record to disk.

`peritheos.diagnostics.birch_murnaghan_finite_strain_diagnostic()` constructs
the conventional Birch-Murnaghan normalized-stress diagnostic from P-V
observations. The returned `BirchMurnaghanFiniteStrainDiagnostic` exposes the
transformed arrays and an optional numerical model curve. Plotting is left to
the caller; the fitting guide links to a complete executable example.

`fit_joint_eos()` estimates reference-isotherm and thermal parameters in one
regression. Reference parameters use dotted names such as `rt_eos.V0`; its
covariance includes reference/thermal cross-correlations and is directly
compatible with `FitResult.eos_uncertainty()`.

`fit_linear_us_up()` returns a `HugoniotFitResult`. With no uncertainties it
performs OLS; shock-velocity uncertainties select WLS, and particle-velocity
uncertainties add latent-coordinate errors-in-variables fitting.
`HugoniotFitResult.to_dict()` and `.to_json()` serialize the fitted model,
covariance, adjusted particle velocities, diagnostics, and solver status.

## Units

```python
from peritheos.units import (
    cell_volume_to_molar_volume,
    convert_density,
    convert_molar_volume,
    convert_pressure,
    convert_temperature,
    density_from_molar_volume,
    molar_volume_to_cell_volume,
    molar_volume_from_density,
)
```

Every conversion accepts scalars or NumPy arrays. Cell-volume conversion uses
angstrom cubed per conventional cell and requires the crystallographic number
of formula units per cell; its molar result is per mole of formula units.
Pressure and temperature conversion imports from `peritheos.utils` are
deprecated compatibility aliases.

## Uncertainty propagation

```python
from peritheos import EOSUncertainty, ParameterUncertainty, PredictionUncertainty
```

`EOSUncertainty` wraps a deterministic EOS and provides `pressure()`, `volume()`,
`bulk_modulus()`, and the generic `evaluate()` method. Each returns a
`PredictionUncertainty` containing the nominal value, standard error, confidence
limits, statistical assumptions, and an optional output covariance matrix.

`FitResult.eos_uncertainty()` constructs the wrapper directly from fitted
parameter covariance. See [Uncertainty in EOS calculations](uncertainty.md).
Partial parameter-error sets are supported. Parameters omitted from the error
mapping or covariance ordering are treated as exact, rather than as having an
unknown error that Peritheos will estimate.
