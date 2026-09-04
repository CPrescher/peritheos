# Peritheos

A Python library for thermodynamic equations of state calculations for solid materials.

Built-in EOS evaluation, inversion, bounded robust fitting, and uncertainty
statistics are implemented in Rust and exposed through the unchanged Python API.
NumPy broadcasting and user-defined `EosBase` subclasses remain supported.
The package includes PEP 561 typing metadata for type checkers and IDEs.

Full model, fitting, units, and development documentation is available at
[peritheos.readthedocs.io](https://peritheos.readthedocs.io/).
Rust users can start with the workflow-oriented [Rust API guide](docs/rust-api.md)
and the runnable examples in the unified Rust crate.
Release history is recorded in the [changelog](CHANGELOG.md).

## Features

- Room temperature equations of state (EOS) implementations
  - Birch-Murnaghan
  - Murnaghan
  - Natural strain (orders 2-4)
  - Modified Tait
  - Vinet
  - Holzapfel
- Thermal equations of state (EOS) implementations
  - Mie-Gruneisen-Debye
  - Mie-Gruneisen-Einstein
  - Linear thermal pressure
  - Second-order temperature-compression thermal pressure
  - Holland-Powell thermal modified Tait
  - Multi-oscillator Gruneisen thermal pressure
- P-V and P-V-T parameter fitting with covariance and diagnostics
- Joint reference-isotherm and thermal fitting with cross-covariance
- Correlated observation errors and robust least-squares losses
- Reproducible fit summaries and versioned JSON export
- EOS prediction uncertainty from fitted covariance or published parameter errors
- Thermoelastic derivatives, heat capacities, and vibrational potentials
- Versioned material and EOS-record catalog with explicit literature provenance,
  calibration/data envelopes, extrapolation enabled by default, inversion, and
  measurement/parameter uncertainty
- A Peritheos-owned `.eosmat` schema and 118-material/172-record EOS library with optional
  diffraction structure, stable identifiers, and Dioptas 0.10 storage-read
  compatibility
- Recursive pressure-scale normalization across Au, Pt, KCl, ruby R1, and
  diamond Raman standards with source-documented cross-calibration edges
- Native wheels for supported CPython releases on Linux, macOS, and Windows

## Unit conventions

- Public pressure and bulk-modulus values are in GPa.
- Temperatures are in K.
- Birch-Murnaghan, Murnaghan, modified Tait, and Vinet accept any consistent
  volume unit.
- Holzapfel and energy-based thermal EOS implementations require molar volume in
  J bar^-1 mol^-1, which is equivalent to cm^3/mol divided by 10.
- The linear and second-order Taylor thermal-pressure corrections use the same
  volume convention as their reference isotherm.

## Installation

```bash
pip install peritheos
```

The latest development version can instead be installed directly from GitHub:

```bash
pip install git+https://github.com/CPrescher/peritheos.git
```

Development and other source installations compile the private native
extension and therefore require a Rust toolchain compatible with Rust 1.83 or
newer. Published PyPI wheels do not require Rust.

## Usage

### Materials and EOS records

Materials group one or more literature-specific EOS records. The calculation
API uses GPa, K, and conventional unit-cell volumes in angstrom cubed. Each EOS
record carries its primary reference, parameter provenance, validity range, and
uncertainty assumptions.

```python
from peritheos import get_material

mgo = get_material("mgo_b1")
tange = mgo.get_eos_record("mgo_b1_tange_2009_vinet")
pressure = tange.pressure(volume=60.0, temperature=2000.0)
recovered_volume = tange.volume(pressure, temperature=2000.0)

prediction = tange.pressure_with_uncertainty(
    volume=60.0,
    temperature=2000.0,
    volume_sigma=0.02,
    temperature_sigma=20.0,
)

# Sokolova markers use the same cell-volume API although their composed EOS
# works internally with molar volume.
gold = get_material("au_fcc").get_eos_record("au_fcc_sokolova_2013")
hot_pressure = gold.pressure(volume=55.0, temperature=2000.0)
```

See [Pressure standards](docs/pressure-standards.md) for EOS records commonly
used in that application, and [Dioptas and `.eosmat`](docs/dioptas-integration.md)
for the shared material library. The [`.eosmat` schema reference](docs/eosmat-schema.md)
defines its fields, equation discriminators, defaults, units, and compatibility
rules.

### Room-temperature equations of state

Third-order Birch-Murnaghan equation of state:

```python
from peritheos.eos.rt import BM3

# V0 may use any volume unit for a room-temperature EOS; K0 is in GPa here.
eos = BM3(V0=50, K0=130, K0_prime=4.3)

# Calculate pressure and bulk modulus at a given volume.
pressure = eos.pressure(V=40)
bulk_modulus = eos.bulk_modulus(V=40)

# Invert the EOS to calculate volume at a given pressure.
volume = eos.volume(P=pressure)

print(f"Pressure: {pressure} GPa")
print(f"Bulk modulus: {bulk_modulus} GPa")
print(f"Recovered volume: {volume}")
```

### Thermal equations of state

Mie-Gruneisen-Debye and Mie-Gruneisen-Einstein models can wrap any of the
room-temperature equations of state:

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye

# Thermal models require molar volume in J bar^-1 mol^-1.
rt_eos = BM3(V0=1.0, K0=160.0, K0_prime=4.0)
eos = MieGruneisenDebye(
    rt_eos=rt_eos,
    Tr=300.0,
    theta0=800.0,
    gamma0=1.5,
    q=1.0,
    n=2,
)

pressure = eos.pressure(V=0.9, T=2000.0)
volume = eos.volume(P=pressure, T=2000.0)
temperature = eos.temperature(P=pressure, V=0.9)

# Predict the hot volume and total pressure when a DAC retains 25% of the
# reference-relative thermal pressure above a 40 GPa cold pressure.
cold_pressure = 40.0
hot_temperature = 2000.0
heated_volume = eos.volume_with_dac_confinement(
    cold_pressure,
    hot_temperature,
    f_dac=0.25,
)
thermal_increment = eos.thermal_pressure_increment(heated_volume, hot_temperature)
heated_pressure = cold_pressure + 0.25 * thermal_increment

# Infer temperature from volumes measured before and during DAC heating.
ambient_volume = 0.80000
heated_volume = 0.80001
temperature_with_dac = eos.temperature_from_volumes(
    V_ambient=ambient_volume,
    V_heated=heated_volume,
    f_dac=0.25,
)
ambient_pressure = eos.rt_eos.pressure(ambient_volume)
heated_pressure_from_pair = ambient_pressure + eos.dac_thermal_pressure(
    heated_volume,
    temperature_with_dac,
    0.25,
)
```

Both DAC methods use the empirical `f_dac * thermal_pressure_increment`
confinement term and require `0 <= f_dac < 1`; report and sensitivity-test the
assumed fraction.

Multi-oscillator thermal EOS with a freely chosen reference isotherm

```python
from peritheos.eos.rt.holzapfel import Holzapfel
from peritheos.eos.thermal import MultiOscillatorGruneisenThermalEOS

# Diamond parameters from Sokolova et al. 2016.
# The thermal model requires molar volume in J bar^-1 (= [cm^3/mol] / 10),
# pressure parameters in GPa, and temperatures in K.
V0 = 0.3414
K0 = 441.5
K0_prime = 3.9  # pressure derivative of bulk modulus at reference volume
QE1o = 1561  # first Einstein characteristic temperature
mE1 = 2.436  # first Einstein number
QE2o = 684  # second Einstein characteristic temperature
mE2 = 0.564  # second Einstein number
delta = -0.506  # additive normalizing constant for the Gruneisen parameter
t = 1.085  # generalized Gruneisen parameter
a_0 = 0  # intrinsic anharmonicity parameter
m = 0  # anharmonic analogue of the Grüneisen parameter
e_0 = 0  # free electrons parameter
g = 0  # electronic analogue of the Grüneisen parameter

n = 1  # number of atoms in the formula unit
z = 6  # atomic number of the formula unit
Tr = 298.15  # in K - Reference temperature

# Initialize the Holzapfel EOS
holzapfel = Holzapfel(V0=V0, K0=K0, K0_prime=K0_prime, n=n, Z=z)

# Compose the thermal correction with the source's Holzapfel isotherm.
eos = MultiOscillatorGruneisenThermalEOS(
    rt_eos=holzapfel,
    Tr=Tr,
    QE1o=QE1o,
    mE1=mE1,
    QE2o=QE2o,
    mE2=mE2,
    delta=delta,
    t=t,
    a_0=a_0,
    m=m,
    g=g,
    e_0=e_0,
    n=n,
)

# Calculate the thermal pressure at a given volume and temperature
V = V0 * 0.8
T = 3000  # in K
thermal_pressure = eos.thermal_pressure(V, T)
rt_pressure = holzapfel.pressure(V)
pressure = eos.pressure(V, T)
recovered_volume = eos.volume(pressure, T)
recovered_temperature = eos.temperature(pressure, V)

print(f"Thermal pressure: {thermal_pressure} GPa")
print(f"RT pressure: {rt_pressure} GPa")
print(f"Total pressure: {pressure} GPa")
print(f"Recovered volume: {recovered_volume} J bar^-1")
print(f"Recovered temperature: {recovered_temperature} K")
```

The equation class is named for its mechanism, not a paper. The validated
Sokolova et al. parameterizations remain available through author/year catalog
identifiers such as `diamond_sokolova_2013`; replacing the reference isotherm
creates a new user-composed model and does not inherit that validation claim.

## Citation and support

Use the repository's `CITATION.cff` to cite Peritheos and cite the original
publication for each EOS used. Reproducible bugs and numerical discrepancies
can be reported through [GitHub Issues](https://github.com/CPrescher/peritheos/issues).
See [SUPPORT.md](SUPPORT.md) for the information needed to investigate a result.
