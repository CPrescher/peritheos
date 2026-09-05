# Shock Hugoniot equations of state

A shock Hugoniot is stored as an EOS record because it supplies an executable
pressure-volume relation. It is nevertheless a constrained shock path, not an
equilibrium surface $P(V,T)$. Peritheos makes that distinction explicit with
`equation_kind="hugoniot"` while keeping the record in the material's common
`eos_records` collection.

This representation is phase and branch specific. If a shock crosses a phase
boundary, use a separate EOS record for each independently parameterized stable
branch. Do not silently extend the low-pressure `Us`-`up` relation into the
high-pressure phase.

## Bundled single-phase records

The initial production set deliberately contains only untransformed principal
branches:

- `mgo_b1_duffy_ahrens_1995_hugoniot_5`: published B1 MgO relation over
  14--133 GPa.
- `nickel_oxide_noguchi_1999_linear_hugoniot_2`: Peritheos OLS fit to the eight
  final-state B1 NiO observations over 17.7--147.6 GPa; separately resolved
  elastic-limit states are excluded.

Both records retain their precursor density, phase, branch domain, uncertainty
provenance, and source-page locations. Transformed, liquid, mixed-phase, and
phase-pooled candidates remain outside the production library pending an
explicit phase-transition representation.

## Linear shock-velocity relation

`LinearUsUpHugoniot` implements

\[
U_s=c_0+s u_p,
\]

together with the Rankine-Hugoniot mass and momentum jump conditions

\[
P-P_0=\rho_0 U_su_p,
\qquad
\frac{V}{V_0}=1-\frac{u_p}{U_s}.
\]

Eliminating the velocities gives an executable pressure-volume path. With
$\mu=1-V/V_0$,

\[
P_H(V)=P_0+\rho_0c_0^2
\frac{\mu}{(1-s\mu)^2}.
\]

The implementation also evaluates density and the specific internal-energy
increase implied by the energy jump condition,

\[
\Delta e_H=\frac{1}{2}(P+P_0)
\left(\frac{1}{\rho_0}-\frac{1}{\rho}\right).
\]

Use density in g/cm³, velocities in km/s, and pressure in GPa. These units make
$\rho_0U_su_p$ numerically equal to GPa. `V0` and `V` may use any mutually
consistent volume unit. A material `HugoniotRecord` makes that consistency
enforceable with a formula-unit volume basis: `V0`, `rho0`, the number of
formula units, and the molar mass must describe the same amount of material.
This check is especially important when precursor and represented phases have
different conventional cells.

```python
from peritheos import LinearUsUpHugoniot

hugoniot = LinearUsUpHugoniot(
    V0=60.4,  # angstrom^3/conventional cell
    rho0=21.45,  # g/cm^3
    c0=3.6,  # km/s
    s=1.55,
    P0=0.0,  # GPa
)

pressure = hugoniot.pressure(48.0)
volume = hugoniot.volume(pressure)
state = hugoniot.state(48.0)
```

The full model API is:

- `state_from_particle_velocity(up)` for the complete coupled shock state
- `pressure(V)` and `volume(P)`
- `particle_velocity(V)` and `shock_velocity(V)`
- `density(V)`
- `specific_internal_energy_change(V)` in MJ/kg
- `tangent_modulus(V)`, equal to $-V\,dP_H/dV$
- direct `shock_velocity_from_particle_velocity(up)`,
  `pressure_from_particle_velocity(up)`, and
  `volume_from_particle_velocity(up)` calculations

Only the compressive branch is admitted: $0<V\leq V_0$ and $1-s\mu>0$.

## Fitting `Us`-`up` observations

`fit_linear_us_up` fits only `c0` and `s`. `V0`, `rho0`, and `P0` describe the
initial state and must be supplied independently.

```python
from peritheos.fitting import fit_linear_us_up

fit = fit_linear_us_up(
    particle_velocity=up,  # km/s
    shock_velocity=Us,  # km/s
    V0=60.4,
    rho0=21.45,
    shock_velocity_sigma=sigma_Us,
    particle_velocity_sigma=sigma_up,
    absolute_sigma=True,
)

print(fit.parameters)  # includes fitted c0 and s plus fixed initial state
pressure = fit.model.pressure(48.0)
```

With no uncertainty arguments this is ordinary least squares (OLS): every
point receives equal weight and the solver minimizes the sum of squared
vertical `Us` residuals. Supplying only `shock_velocity_sigma` gives weighted
least squares. Supplying `particle_velocity_sigma` as well gives an
errors-in-variables fit with latent corrected particle velocities. A per-point
2-by-2 `observation_covariance`, ordered as `(Us, up)`, represents correlated
measurement errors.

The result reports the fitted model, covariance and correlation of free
parameters, raw and weighted residuals, adjusted particle velocities, and
standard fit diagnostics. `fit.eos_uncertainty()` propagates the fitted
parameter covariance through any public Hugoniot quantity.

## EOSMAT records and phase transitions

`HugoniotRecord` is a specialized EOS record. It remains in
`material.eos_records`, while `material.hugoniot_records` and
`material.equilibrium_records` provide typed filtered views. Its path,
precursor, operational volume basis, and branch domain are required rather
than inferred.
`record.state_from_particle_velocity(up)` is the preferred evaluation entry
point because it begins from the published `Us`-`up` coordinate and returns a
coupled state whose volume uses the record's public unit.

```python
from peritheos import (
    HugoniotInitialState,
    HugoniotBranchDomain,
    HugoniotRecord,
    HugoniotVolumeBasis,
)
from peritheos.materials import LiteratureReference, ValidityRange

record = HugoniotRecord(
    identifier="example_beta_shock_1",
    name="Example transformed beta-phase Hugoniot",
    material="X",
    phase="beta",  # represented final-state phase
    cell_contents="4 formula units per conventional unit cell",
    eos=hugoniot,
    reference_temperature=298.15,
    reference=LiteratureReference(
        authors="Example et al.",
        year=2026,
        title="Example shock study",
        doi="10.xxxx/example",
        locations=("Table 2",),
    ),
    validity=ValidityRange(
        pressure_gpa=(50.0, 150.0),
        temperature_k=(298.15, 298.15),  # precursor temperature
        volume_ratio=(0.65, 0.85),
    ),
    parameter_provenance={"c0": "Table 2", "s": "Table 2"},
    loading_path="principal",
    branch_kind="transformed",
    initial_state=HugoniotInitialState(
        phase="alpha",
        material_identifier="x_alpha",
        temperature_k=298.15,
        pressure_gpa=hugoniot.P0,
        density_g_cm3=hugoniot.rho0,
    ),
    volume_basis=HugoniotVolumeBasis(
        formula_units=4.0,
        molar_mass_g_mol=195.05412814602,
    ),
    branch_domain=HugoniotBranchDomain(
        particle_velocity_km_s=(1.2, 3.4),
        kind="phase_stability",
        boundary_status="reported_exactly",
    ),
)
```

A canonical record has this shape:

```json
{
  "identifier": "example_alpha_principal_hugoniot_1",
  "label": "Example alpha-phase principal Hugoniot",
  "equation_kind": "hugoniot",
  "loading_path": "principal",
  "branch_kind": "untransformed",
  "initial_state": {
    "phase": "alpha",
    "material_identifier": "example_alpha",
    "temperature_k": 298.15,
    "pressure_gpa": 0.0,
    "density_g_cm3": 8.0
  },
  "volume_basis": {
    "kind": "formula_units",
    "formula_units": 1.0,
    "molar_mass_g_mol": 48.17712608
  },
  "branch_domain": {
    "particle_velocity_km_s": [0.0, 3.0],
    "kind": "phase_stability",
    "boundary_status": "reported_exactly"
  },
  "eos": {
    "type": "LinearUsUpHugoniot",
    "model": "linear_us_up_hugoniot",
    "parameters": {
      "V0": 10.0,
      "rho0": 8.0,
      "c0": 4.0,
      "s": 1.5,
      "P0": 0.0
    }
  }
}
```

`loading_path` records `principal` or `precompressed` loading, independently
of whether `branch_kind` is `untransformed` or `transformed`. Thus a transformed
phase reached from a precompressed state is representable without conflating
the two concepts. Reshocks are deliberately excluded until a model records the
already-shocked material velocity and velocity frame. A
Hugoniot record cannot contain a `thermal` component and cannot serve as the
reference isotherm of a thermal EOS. Temperature is metadata on its initial
state, not an independent input to the path.

When a material has multiple stable shock phases, give each phase-specific
material its own Hugoniot record and coefficients. Multiple parameterizations
of the same represented phase may remain in that material's `eos_records`
array. The material's top-level `phase` identifies the phase represented by
the branch; `initial_state.phase` identifies the precursor state from which
that shock branch begins. Peritheos does not select a phase branch
automatically because a reported transition pressure alone is insufficient to
reconstruct mixed-phase kinetics and path history.

For a principal single shock, phase-specific low- and high-pressure branches
normally share the same precursor and `V0`, `rho0`, and `P0`. The represented
phase changes; the initial state does not. `V0` is the precursor volume
containing `volume_basis.formula_units`, which need not equal one precursor
conventional cell. The basis value must match the represented material's
`formula_units_per_cell`, because public record volumes are represented-phase
conventional-cell volumes. Peritheos verifies this using
`volume_basis.molar_mass_g_mol` and rejects `V0`/`rho0` combinations that
disagree with the stated mass by more than 0.1%.

Every record also declares a `branch_domain` in particle velocity. Record-level
evaluation checks this domain by default, including
`state_from_particle_velocity`, `pressure`, `volume`, and the shock-specific
quantities. Use the raw `record.eos` model only when deliberate mathematical
extrapolation is required. A domain may describe phase stability, experimental
coverage, or a source-recommended range; its kind and boundary status make that
scientific meaning explicit.

## Data provenance

Store or redistribute experimental rows only when their terms permit it.
Published `Us`-`up` coefficients can be entered as parameters with a precise
literature citation and source location; if coefficients are derived from a
table or SESAME calculation, record that transformation and the source table
identifier in a structured `derivation` block and set `record_kind="derived"`.
For SESAME, include the table number/version, sampling domain, transformation or
fit method, software, and an access/licensing note. A fitted relation is not a
substitute for redistribution rights to the underlying table.
