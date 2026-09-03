# Getting started

## Installation

Install a published wheel from PyPI:

```bash
python -m pip install peritheos
```

Wheels do not require a Rust installation. Editable installs and direct source
installs build the private native extension and require Rust 1.83 or newer.

## Start with a literature-backed material record

For published calculations, the record API is the safest starting point. It
fixes the pressure, temperature, and conventional-cell volume units; carries
the primary reference, calibration coverage, and uncertainty metadata. The EOS
remains evaluable outside that coverage; range enforcement is opt-in.

```python
from peritheos import get_material, search_eos_records

mgo = get_material("mgo")
record = mgo.get_eos_record("mgo_sokolova_2013_holzapfel_4")

pressure = record.pressure(volume=60.0, temperature=2000.0)
recovered_volume = record.volume(pressure, temperature=2000.0)
assert record.within_calibration_range(volume=60.0, temperature=2000.0)

prediction = record.pressure_with_uncertainty(
    volume=60.0,
    temperature=2000.0,
    volume_sigma=0.02,
    temperature_sigma=20.0,
)
```

Here volume is in angstrom cubed per conventional unit cell, pressure is in
GPa, and temperature is in K. The returned object came directly from the
bundled `.eosmat` document. Search the same executable catalog without loading
raw documents yourself:

```python
hot_gold_scales = search_eos_records(
    formula="Au",
    thermal=True,
    pressure_gpa=(0.0, 200.0),
    temperature_k=2000.0,
)
```

See [Material catalog](catalog.md) for every discovery filter,
[Pressure standards](pressure-standards.md) for pressure-calibration guidance,
and [Exploring the material library](notebooks/exploring-material-library.ipynb)
for a longer walkthrough.

## Construct a room-temperature model

Use the equation classes directly when fitting a new dataset or reproducing a
parameterization that is not represented by a catalog record.

```python
from peritheos.eos.rt import Vinet

eos = Vinet(V0=10.0, K0=160.0, K0_prime=4.2)
pressure = eos.pressure(8.5)
bulk_modulus = eos.bulk_modulus(8.5)
recovered_volume = eos.volume(pressure)
```

`V0` and `V` may use any consistent volume unit for an isothermal EOS. `K0`,
pressure, and bulk modulus must use the same pressure unit; Peritheos documents
and tests public pressure values as GPa.

All methods accept NumPy arrays:

```python
import numpy as np

volumes = np.linspace(8.0, 10.0, 21)
pressures = eos.pressure(volumes)
```

For the available isothermal formulations and guidance on choosing one, see
[EOS models](models.md#isothermal-models). Their complete mathematical
definitions are in the [isothermal equation reference](equation-reference.md#isothermal-equations),
while constructor signatures and common methods are collected in the
[isothermal API reference](api.md#isothermal-equations-of-state).

## Thermal pressure and properties

Thermal models require molar volume in J bar^-1 mol^-1. This equals molar
cm^3/mol divided by ten.

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye

reference = BM3(V0=1.0, K0=160.0, K0_prime=4.0)
eos = MieGruneisenDebye(
    rt_eos=reference,
    Tr=300.0,
    theta0=800.0,
    gamma0=1.5,
    q=1.0,
    n=2,
)

volume = 0.9
temperature = 1800.0
pressure = eos.pressure(volume, temperature)
recovered_temperature = eos.temperature(pressure, volume)
kt = eos.bulk_modulus(volume, temperature)
alpha = eos.thermal_expansivity(volume, temperature)
cv = eos.molar_heat_capacity_v(volume, temperature)
cp = eos.molar_heat_capacity_p(volume, temperature)
ks = eos.adiabatic_bulk_modulus(volume, temperature)
```

Thermal models combine a reference isothermal EOS with a temperature-dependent
pressure contribution. See [Thermal models](models.md#thermal-models) for the
allowed reference-model combinations, the
[thermal equation reference](equation-reference.md#thermal-equations) for the
implemented formulas, and [Units and reference states](units.md) for the molar
volume convention. Consult [Thermoelastic properties](thermoelastic-properties.md)
before using caloric or free-energy methods.

Use `temperature(P, V)` when the total pressure and volume at the heated state
are both known. If only the reference-temperature volume and heated volume are
measured, the ambient pressure is not the total pressure at high temperature.
That case requires the advanced experimental boundary-condition model described
under [Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md).

The numerical behavior of direct pressure-volume inversion is summarized under
[Numerical inversion](equation-reference.md#numerical-inversion). All thermal
constructors and methods are listed in the
[thermal API reference](api.md#thermal-equations-of-state).

## Broadcasting P and T

Thermal state variables use NumPy broadcasting. A volume column and temperature
row produce a full state grid:

```python
volumes = np.array([0.8, 0.9, 1.0])[:, None]
temperatures = np.array([300.0, 1000.0, 2000.0])[None, :]
pressure_grid = eos.pressure(volumes, temperatures)
```

The same broadcasting rules apply to supported pressure, volume, temperature,
and thermoelastic-property calculations. Method availability by model is listed
in the [thermal API reference](api.md#thermal-equations-of-state).

## Where to go next

| Goal | Detailed documentation |
|---|---|
| Compare EOS families and choose a model | [EOS models](models.md) |
| Inspect the exact implemented formulas | [Equation reference](equation-reference.md) |
| Fit isothermal, thermal, or joint P-V-T data | [Fitting P-V and P-V-T data](fitting.md) |
| Propagate parameter and measurement errors | [Uncertainty in EOS calculations](uncertainty.md) |
| Analyze two-volume laser-heated DAC measurements | [Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md) |
| Calculate moduli, expansivity, heat capacities, or free energies | [Thermoelastic properties](thermoelastic-properties.md) |
| Check units and reference-state conventions | [Units and reference states](units.md) |
| Look up constructors and methods | [API reference](api.md) |
| Review numerical and literature checks | [Validation](validation.md) |
