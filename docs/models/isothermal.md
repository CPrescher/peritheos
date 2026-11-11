# Isothermal models

Select a model name for its equations, domains, and reference-state conventions.

## Models and parameters

| Import | Parameters after `V0`, `K0` | Notes |
|---|---|---|
| [`BM2`](../equations/isothermal/birch-murnaghan.md) | none | Second-order Birch-Murnaghan; implied `K0_prime = 4` |
| [`BM2Reference`](../equations/isothermal/birch-murnaghan.md#finite-pressure-bm2-reference) | required `P0` | BM2 pressure increment above `P0`; `V0` and `K0` refer to that pressure |
| [`BM3`](../equations/isothermal/birch-murnaghan.md) | `K0_prime` | Standard third-order Eulerian finite strain |
| [`BM4`](../equations/isothermal/birch-murnaghan.md) | `K0_prime`, `K0_double_prime` | Fourth-order Eulerian finite strain |
| [`Murnaghan`](../equations/isothermal/murnaghan.md) | `K0_prime` | Assumes bulk modulus linear in pressure |
| [`NaturalStrain2`](../equations/isothermal/natural-strain.md) | none | Hencky strain; implied `K0_prime = 2` |
| [`NaturalStrain3`](../equations/isothermal/natural-strain.md) | `K0_prime` | Third-order Poirier-Tarantola logarithmic EOS |
| [`NaturalStrain4`](../equations/isothermal/natural-strain.md) | `K0_prime`, `K0_double_prime` | Fourth-order natural strain |
| [`ModifiedTait`](../equations/isothermal/modified-tait.md) | `K0_prime`, `K0_double_prime` | Analytically invertible Tait family |
| [`Baonza`](../equations/isothermal/baonza.md) | `K0_prime` | Pseudospinodal form with the conventional fixed exponent beta=0.85 |
| [`Vinet`](../equations/isothermal/vinet.md) | `K0_prime` | Rydberg-Vinet interatomic-potential form |
| [`RydbergStacey`](../equations/isothermal/rydberg-stacey.md) | `K0_prime`, `K_infinity_prime` | Generalized high-compression Rydberg form with an adjustable infinite-pressure derivative |
| [`Morse3`](../equations/isothermal/morse.md) | `K0_prime` | Three-dimensional Morse-potential form |
| [`SunMorse3`](../equations/isothermal/morse.md) | `K0_prime` | Volume-analytic Sun--Morse form with n=3 |
| [`SunMorse4`](../equations/isothermal/morse.md) | `K0_prime` | Volume-analytic Sun--Morse form with n=4 |
| [`Holzapfel`](../equations/isothermal/holzapfel.md) | `K0_prime`, `n`, `Z` | High-compression form with strict molar units |

All isothermal models implement:

```python
pressure(V)
bulk_modulus(V)
volume(P)
```

[Equation index](../equation-reference.md) · [Choosing a model](../models.md#selecting-a-model)
