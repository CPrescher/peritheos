# Holmes et al. (1989): theoretical platinum pressure scale

## Scope and primary evidence

Holmes, Moriarty, Gathers, and Nellis, “The equation of state of platinum to
660 GPa (6.6 Mbar),” *Journal of Applied Physics* **66**, 2962-2967 (1989),
[doi:10.1063/1.344177](https://doi.org/10.1063/1.344177), combines new shock
experiments with a first-principles theoretical Pt EOS. This audit used the
author-uploaded full article and visually checked Tables II-IV and Figures 3-5.
The article has no supplement or data repository.

The production record is specifically the theoretical 300 K pressure scale.
It is not the linear experimental `Us-up` Hugoniot in Table II, and Table III's
seven shocked states are not 300 K pressure-volume observations.

## Cold and 300 K construction

Equation (7) separates zero-temperature, ion-thermal, and electron-thermal
pressures. The zero-temperature `E0(V)` and `P0(V)` curves were calculated for
fcc Pt over `0.6 < V/V0 < 1.1` with scalar-relativistic Kohn-Sham LDA and the
LMTO atomic-sphere approximation, including the combined correction and
`s`, `p`, `d`, and `f` angular-momentum components. The paper says the only
physical input is atomic number `Z=78`.

That methodological statement is not a deposited numerical calculation. The
paper does not tabulate the LMTO `E0(V)` or `P0(V)` states, k-point or radial
meshes, convergence settings, individual density-of-states values, the 300 K
fit grid, weights, residuals, covariance, or fit statistic. Figure 3 plots only
bcc-fcc and hcp-fcc energy differences; Figures 4 and 5 plot curves. Digitizing
those curves would add graphical pseudo-precision and would not recover the
source regression.

The published equilibrium curve itself is nevertheless exactly reconstructable.
Equation (11) gives

\[
P_{300}(X)=P_T\frac{1-X}{X^2}\exp[\eta(1-X)],\qquad
X=(V/V_0)^{1/3},
\]

and Table IV reports `P_T=798.31 GPa`, `eta=7.2119`, with the experimental
`V0=101.9 bohr^3/atom`. On Peritheos's four-atom conventional fcc-cell basis,
`V0=60.4000884 A^3`. The exactly equivalent Vinet mapping is
`K0=P_T/3=266.103333333 GPa` and
`K0'=1+eta/1.5=5.807933333`. Table IV separately prints the rounded
interpretive values `B_T=266 GPa` and `B_T'=5.81`; using those rounded values
as the executable coefficients does not reproduce the more precise Equation
(11) parameterization exactly.

The bundled eight-row checkpoint grid evaluates Equation (11) directly from
zero pressure through its stated 550 GPa limit. It is source-derived validation
data, not independent fit observations, and is excluded from `fit_datasets`.

## Thermal extension and shock qualification

Equation (12) adds the approximate constant thermal pressure

\[
P(X,T)=P_{300}(X)+\alpha_T B_T(T-300),
\]

using Table IV's `alpha_T=0.261e-4 K^-1` and rounded `B_T=266 GPa`. Thus the
stored slope remains `0.0069426 GPa/K`. The authors limit this approximation to
temperatures below 2000 K and state that its magnitude remains below 12 GPa.
This compact extension is not the full theoretical Hugoniot construction:
Equation (8) uses a Slater-derived volume-dependent ionic Gruneisen parameter,
while Equations (9)-(10) use an LMTO density of states and an approximately
constant electronic Gruneisen parameter near 1.8. Those numerical inputs are
not tabulated, so the theoretical Hugoniot cannot be independently reconstructed.

Table III's seven new shock rows are bundled verbatim, including the source's
two-sigma bounds. Their printed pressures and densities satisfy the
Rankine-Hugoniot momentum and mass identities within rounding. The authors
compare these rows and older McQueen data with their theoretical Hugoniot and
use the agreement to qualify the 300 K curve to at least 10% in pressure. The
older row-level McQueen states and the numerical theoretical Hugoniot are not
printed. Consequently the shock rows are retained as qualification evidence and
are never fitted as a 300 K isotherm.

## Disposition

- The published Equation (11) curve is exactly and independently executable.
- Its underlying LMTO-to-Vinet coefficient fit remains `not_refittable` because
  no numerical theoretical fit grid or fitting protocol is published.
- The seven Table III experiments remain shock-Hugoniot evidence only.
- No Figure 4 or Figure 5 points were digitized, and no shock pressure-density
  pair was relabeled as an equilibrium 300 K observation.
