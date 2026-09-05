# Criniti et al. (2021): MgSiO3 bridgmanite

## Scope and sources

This audit covers all five LitCurate candidates attached to Criniti et al.,
*Single-Crystal Elasticity of MgSiO3 Bridgmanite to Mid-Lower Mantle Pressure*,
*Journal of Geophysical Research: Solid Earth* **126**, e2020JB020967 (2021),
[doi:10.1029/2020JB020967](https://doi.org/10.1029/2020JB020967).
The scientific authorities are the [open-access final article](https://epub.uni-bayreuth.de/5769/1/2020JB020967.pdf)
and the authors' [Figshare measurement workbook](https://doi.org/10.6084/m9.figshare.14265176.v1).
LitCurate was used only for discovery. An exact DOI search found no existing
item in the local Zotero library; this audit did not import it.

## Candidate disposition

| Candidate | Source quantity | Decision |
|---|---|---|
| `litcurate_7c709a358383b534` | Voigt adiabatic bulk-modulus finite-strain fit, `K_S0=257.1(6) GPa`, `K'_S0=3.71(4)` | **Rejected as P(V).** This is a fit of elastic modulus versus strain, not a volumetric pressure equation. |
| `litcurate_1473f5b46c728a08` | Reuss adiabatic bulk-modulus finite-strain fit, `K_SR0=256.7(4) GPa`, `K'_SR0=3.70(3)` | **Rejected as P(V).** It is an elasticity truncation and cannot be substituted for the isothermal pressure integral. |
| `litcurate_9bb38c0a21fcf50b` | Section 3.3 primary absolute scale, `K_TR0=254.5(4) GPa`, `K'_TR0=3.73(2)` | **Accepted** as `bridgmanite_criniti_2021_absolute_bm3`. Equation 3 is an explicit BM3 P(V) relation. |
| `litcurate_19a5bbcd8ab9c0e6` | Independent Jacobsen-ruby P-V BM3, `K_T0=255(10) GPa`, `K'_T0=3.8(2)` | **Hold.** The article does not report this fit's V0, and the public workbook contains no row-wise ruby R1 shifts or explicitly identified Jacobsen-pressure series. Its missing reference state is not inferred. |
| `litcurate_9c32311783c63f9c` | Chantel et al. (2012) ultrasonic comparison value | **Citation trace only.** It must be audited from the underlying primary publication and is not a Criniti source record. |

Production count from this paper is therefore **one**, not four.

## Volume basis and data

The source studies stoichiometric MgSiO3 bridgmanite in the conventional Pbnm
cell (`Z=4`). Table 1 reports 14 room-temperature densities from 0.00010 to
78.8 GPa; the Figshare workbook supplies the underlying sound velocities and
crystal-specific densities. The bundled compact dataset retains each Table 1
pressure, density, and stated one-standard-deviation uncertainty. Cell volumes
are derived with

`Vcell = 4 M(MgSiO3) / (rho NA)`,

using `M=100.3875 g/mol` and the exact Avogadro constant. The ambient
`4.1045(10) g/cm3` gives `V0=162.453274(39579) A3`. That value is the explicit
integration reference for Equation 3; it is not presented as an independently
fitted coefficient.

## Numerical reproduction

The script `scripts/reproduce_criniti_2021_bridgmanite.py` evaluates the
published BM3 (`V0=162.453274 A3`, `K0=254.5 GPa`, `K0'=3.73`) against all 14
rounded Table 1 pressure-density pairs. Pressure RMSE is `0.0658748178 GPa`
and the maximum absolute residual is `0.2434557113 GPa`, at the highest-pressure
row where density is printed to only `0.01 g/cm3`. Recomputing every volume from
density agrees with the stored values within `5e-7 A3`.

The direct curve-data agreement supports the source's executable absolute
P(V) equation. The exact intermediate K_SR-to-K_TR regression is not recreated:
it also depends on supporting-information thermoelastic inputs and the original
global acoustic fit. The published coefficients are retained rather than
replaced by a regression to rounded table values.

## Calibration and limitations

The accepted equation is a self-consistent primary scale: the authors convert
the measured Reuss adiabatic bulk modulus to the isothermal value and integrate
it over measured volume. No secondary pressure marker enters the final P(V)
expression. The separate Jacobsen-ruby comparison remains documented but
non-executable until its exact V0 and point assignment can be recovered. No
parameter covariance is published.
