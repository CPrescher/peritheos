# Validation

Peritheos uses several complementary validation layers.

## Reference-state identities

Every isothermal model is tested for

\[
P(V_0)=0, \qquad K(V_0)=K_0.
\]

Models parameterized by derivatives are also checked numerically for their
supplied $K_0'$ and $K_0''$. Every thermal model is checked for zero thermal
pressure at its reference temperature.

## Derivative identities

Analytic isothermal bulk moduli are compared with numerical
`-V dP/dV`. Mie-Gruneisen characteristic temperatures are checked against

\[
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

Caloric models are checked for

\[
C_V=\left(\frac{\partial E}{\partial T}\right)_V,
\qquad F=E-TS,
\qquad \frac{K_S}{K_T}=\frac{C_P}{C_V},
\]

and for the public-unit form

\[
C_P-C_V=10^4\alpha^2K_TVT.
\]

## Round trips and array behavior

All model families are tested for P-to-V and P,V-to-T round trips, scalar
behavior, NumPy arrays, and P-V-T broadcasting. The two-volume DAC inversion is
checked independently against its reduced pressure equation for every thermal
model, including both reduced and complete Sokolova configurations. Its
`f_dac` domain, non-heated volume pairs, non-invertible thermal models, and
fraction sensitivity are also tested. Invalid volumes, temperatures, non-finite
parameters, singular parameter sets, and states outside analytic domains are
explicitly rejected.

## Fitting recovery

Synthetic P-V and P-V-T grids with known generating parameters are refitted.
The suite verifies parameter recovery, fixed parameters, bounds, absolute
uncertainty scaling, covariance dimensions, and input errors.

## Literature regressions

The Sokolova diamond model includes numerical regression cases derived from the
accompanying Excel calculation. These intentionally validate spreadsheet
compatibility rather than a separate literal transcription of the journal
article's printed equations; the distinction is documented under
[Paper versus spreadsheet](equation-reference.md#paper-versus-spreadsheet).
Additional equation-level cases and their source DOIs are stored in
`crates/peritheos-core/tests/data/literature_reference_cases.json`; keeping the
values in the packageable core crate makes changes to scientific baselines
visible in review and lets the published crate run its own tests.
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
