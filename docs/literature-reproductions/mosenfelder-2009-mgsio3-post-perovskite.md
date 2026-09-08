# Mosenfelder 2009 MgSiO3 post-perovskite

## Audit result

The Mosenfelder et al. (2009) shock-static BM3-MGD record is reproducible. The
previous `not_refittable` classification resulted from treating shock
temperature as a required coordinate and from omitting the 48 static P-V-T
observations of Guignot et al. (2007). The source instead closes each solid
shock state in pressure-density space with Rankine-Hugoniot energy.

The source model is not a 300 K BM3 isotherm. It is a third-order
Birch-Murnaghan reference **isentrope**, with

\[
P(V,T)=P_S(V)+\frac{\gamma(V)}{V}
\{E_D[V,T]-E_D[V,T_S(V)]\},
\quad T_S=T_r\Theta(V)/\Theta_0.
\]

It also fits the high-temperature heat-capacity limit
$C_{V,m}=1.035(50)$ J g$^{-1}$ K$^{-1}$. The executable record stores its exact
molar conversion, `Cvmax = 103.9010625 J mol^-1 K^-1`, using the catalog molar
mass 100.3875 g mol$^{-1}$.

## Inputs and lineage

The two fit inputs are deliberately separate:

- `mgsio3-post-perovskite-guignot-2007-table1-pvt.csv` contains all 48 static
  observations (24 at 300 K and 24 laser-heated to 2535 K).
- `mgsio3-post-perovskite-mosenfelder-2009-table2-shock.csv` contains six
  source-reported PPv shock rows. Flyer velocity, initial density, shock transit
  velocity, and populated pyrometric temperatures are measurements. Particle
  velocity, pressure, and shock density are source-derived by impedance matching
  and Rankine-Hugoniot relations; the phase label is an interpretation.

No derived reference-isentrope state is represented as a measured static datum.
The separate `mgsio3-post-perovskite-mosenfelder-2009-derived-reference-isentrope.csv`
diagnostic contains the reference-isentrope and modeled Hugoniot states computed
from the shock-density inputs. Its dataset kind, column lineage, and notes mark
every modeled quantity as derived, and it is not registered as a fit input.

## Thermal reduction

For a shocked density $\rho_H$, equations (5)-(9) of Mosenfelder et al. (2007)
give

\[
E_H=\tfrac12P_H(1/\rho_0-1/\rho_H),\quad
E_H=E_{tr}+E_S+\frac{P_H-P_S}{\rho_H\gamma},
\]

which is linear in $P_H$ and requires no shock temperature. Transition energies
come from Mosenfelder et al. (2009), Table 3. For the Pv80Mj20 aggregate the
published expression uses majorite mole fraction but the experiment reports
20 vol%; because no conversion used by the authors is published, the audit uses
0.20 and records this as a residual lineage uncertainty.

Table 2 includes three Luo et al. pyrometric temperatures. Section 3 states
explicitly that shock-temperature observations enter the fit for MgSiO3
**liquid**; it does not identify them as constraints for the solid PPv fit.
Accordingly they remain independent diagnostics and missing temperatures are
never imputed.

## Numerical reproduction

The dedicated script minimizes the source-described errors-in-variables
objective, projecting static observations in P, V, and T and shock observations
in P and density. With the published counting convention (54 observations, six
free coefficients), it gives:

| Parameter | Published | Reproduction |
|---|---:|---:|
| $K_{0S}$ (GPa) | 225(2) | 224.558 |
| $K'_{0S}$ | 4.21(7) | 4.21988 |
| $\gamma_0$ | 2.61(67) | 2.51834 |
| $q$ | 2.1(8) | 2.06726 |
| $C_{V,m}$ (J g$^{-1}$ K$^{-1}$) | 1.035(50) | 1.04721 |
| $\Theta_0$ (K) | 990(146) | 998.675 |
| reduced $\chi^2$ | 0.18 | 0.180046 |

The published coefficients themselves reproduce the static pressures with
0.306 GPa RMSE. Their 9.12 GPa pressure RMSE at the rounded reported shock
densities is not the fit residual: the source objective allows P-density
projection within both reported uncertainties.

Run the audit with:

```bash
uv run --frozen python scripts/reproduce_mosenfelder_2009_mgsio3_post_perovskite.py
```

Primary sources: Mosenfelder et al. (2009), DOI
`10.1029/2008JB005900`; Mosenfelder et al. (2007), DOI
`10.1029/2006JB004364`; Guignot et al. (2007), DOI
`10.1016/j.epsl.2007.01.025`.
