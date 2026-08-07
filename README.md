# Peritheos

A Python library for thermodynamic equations of state calculations for solid materials.

## Features

- Room temperature equations of state (EOS) implementations
  - Birch-Murnaghan
  - Vinet
  - Holzapfel
- Thermal equations of state (EOS) implementations
  - Sokolova 2016, including its complete thermal-pressure parameter set

## Unit conventions

- Public pressure and bulk-modulus values are in GPa.
- Temperatures are in K.
- Birch-Murnaghan and Vinet accept any consistent volume unit.
- Holzapfel and Sokolova 2016 require molar volume in J bar^-1, which is
  equivalent to cm^3/mol divided by 10.

## Installation

```bash
pip install peritheos
```

The latest development version can instead be installed directly from GitHub:

```bash
pip install git+https://github.com/CPrescher/peritheos.git
```

## Usage

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

Diamond thermal equation of state from sokolova et al. 2016

```python
from peritheos.eos.rt.holzapfel import Holzapfel
from peritheos.eos.thermal.sokolova2016 import Sokolova2016

# Diamond parameters from Sokolova et al. 2016.
# The thermal model requires molar volume in J bar^-1 (= [cm^3/mol] / 10),
# pressure parameters in GPa, and temperatures in K.
V0 = 0.3414
K0 = 441.5
K0_prime = 3.9  # pressure derivative of bulk modulus at reference volume
QE1o = 684  # first Einstein characteristic temperature
mE1 = 0.564  # first Einstein number
QE2o = 1561  # second Einstein characteristic temperature
mE2 = 2.436  # second Einstein number
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

# Initialize the Sokolova 2016 EOS
sokolova = Sokolova2016(
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
)

# Calculate the thermal pressure at a given volume and temperature
V = V0 * 0.8
T = 3000  # in K
thermal_pressure = sokolova.thermal_pressure(V, T)
rt_pressure = holzapfel.pressure(V)
pressure = sokolova.pressure(V, T)
recovered_volume = sokolova.volume(pressure, T)

print(f"Thermal pressure: {thermal_pressure} GPa")
print(f"RT pressure: {rt_pressure} GPa")
print(f"Total pressure: {pressure} GPa")
print(f"Recovered volume: {recovered_volume} J bar^-1")
```
