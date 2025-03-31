# Peritheos

A Python library for thermodynamic equations of state calculations for solid materials.

## Features

- Room temperature equations of state (EOS) implementations
  - Birch-Murnaghan
  - Vinet
  - Holzapfel
- Thermal equations of state (EOS) implementations
  - Sokolova 2016

## Installation

<!-- installation from github -->

```bash
pip install git+https://github.com/yourusername/peritheos.git
```

## Usage

### Room temperature equations of state

Holzapfel equation of state

```python
from peritheos.eos.rt import BM2

# Initialize the BM2 EOS
bm2 = BM2(V0=50, K0=130, K0_prime=4.3)

# Calculate the pressure at a given Volume
pressure = holzapfel.pressure(V=40)

# Calculate the bulk modulus at a given Volume
bulk_modulus = holzapfel.bulk_modulus(V=40)

# Calculate the bulk modulus derivative at a given Volume
bulk_modulus_derivative = holzapfel.bulk_modulus_derivative(V=40)

print(f"Pressure: {pressure} CPa")
print(f"Bulk modulus: {bulk_modulus} GPa")
print(f"Bulk modulus derivative: {bulk_modulus_derivative}")
```

### Thermal equations of state

Diamond thermal equation of state from sokolova et al. 2016

```python
from peritheos.eos.rt.holzapfel import Holzapfel
from peritheos.eos.thermal.sokolova2016 import Sokolova2016

# diamond parameters from sokolova et al. 2016
V0 = 0.3414  # in JBar^-1 (same as [cm^3/mol]/10)
K0 = 4415  # in kbar
K0_prime = 3.9  # pressure derivative of bulk modulus at reference volume
Theta_1 = 684  # Einstein characteristic temperature, Theta_1 in [K]
m1 = 0.564  # The first Einstein number
Theta_2 = 1561  # Einstein characteristic temperature, Theta_02 in [K]
m2 = 2.436  # The second Einstein number
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
sokolova = Sokolova2016(rt_eos=holzapfel, Tr=Tr, QE1o=QE1o, mE1=mE1, QE2o=QE2o, mE2=mE2, delta=delta, t=t, a_0=a_0, m=m, g=g, e_0=e_0)

# Calculate the thermal pressure at a given volume and temperature
V = 0.3414 * 0.8  # in JBar^-1 (same as [cm^3/mol]/10) with 10% compression
T = 3000  # in K
thermal_pressure = sokolova.thermal_pressure(V, T)
rt_pressure = holzapfel.pressure(V)
pressure = sokolova.pressure(V, T)

print(f"Thermal pressure: {thermal_pressure} kbar")
print(f"RT pressure: {rt_pressure} kbar")
print(f"Total pressure: {pressure} kbar")
```
