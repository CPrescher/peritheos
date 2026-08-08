# Validation

Peritheos uses several complementary validation layers.

## Reference-state identities

Every isothermal model is tested for

\[
P(V_0)=0, \qquad K(V_0)=K_0.
\]

Models parameterized by derivatives are also checked numerically for their
supplied (K'_0) and (K''_0). Every thermal model is checked for zero thermal
pressure at its reference temperature.

## Derivative identities

Analytic isothermal bulk moduli are compared with numerical
`-V dP/dV`. Mie-Gruneisen characteristic temperatures are checked against

\[
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

Caloric models are checked for `C_V = dE/dT`, `F = E - TS`,
`C_P - C_V = alpha^2 K_T V T`, and `K_S/K_T = C_P/C_V`.

## Round trips and array behavior

All model families are tested for P-to-V inversion, scalar behavior, NumPy
arrays, and P-V-T broadcasting. Invalid volumes, temperatures, non-finite
parameters, singular parameter sets, and states outside analytic domains are
explicitly rejected.

## Fitting recovery

Synthetic P-V and P-V-T grids with known generating parameters are refitted.
The suite verifies parameter recovery, fixed parameters, bounds, absolute
uncertainty scaling, covariance dimensions, and input errors.

## Literature regressions

The Sokolova diamond model includes published-model numerical regression cases.
Additional equation-level cases and their source DOIs are stored in
`tests/data/literature_reference_cases.json`; keeping the values in a data file
makes changes to scientific baselines visible in review.
The cases cover Birch-Murnaghan orders two through four, Vinet, Holzapfel,
natural strain, modified Tait, Murnaghan, Debye, Einstein, and Holland-Powell
families. Limiting forms and independent derivative identities reduce the risk
of copying the same algebraic mistake into expected values.

The fitting suite also compares a fixed-reference second-order Birch-Murnaghan
fit with its closed-form weighted least-squares solution, including the
absolute covariance. This checks the optimizer and covariance calculation
against an independently solvable statistical result.

Before publishing fitted parameters, users should still compare against the
original fitting program or publication used by their community and report the
model order, reference state, units, weighting, and covariance convention.
