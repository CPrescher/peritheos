# Noguchi et al. (1999) NiO shock-to-isotherm audit

## Outcome

The published 300 K BM3 record is now reproducible as a qualified
shock-to-isotherm fit. The eight bundled Table 1 rows remain Hugoniot states;
none is fitted directly as a room-temperature pressure-volume observation.
The executable reduction first calculates the shock temperature, subtracts
the Mie--Gruneisen Debye thermal pressure at fixed volume, and fits BM3 only to
the resulting 300 K pressures.

The fit to the eight rounded journal rows gives `K0 = 182.4351 GPa` and
`K0' = 4.12465`, compared with the published `191 GPa` and `3.9`. The
residual-scaled refit standard errors are `11.4946 GPa` and `0.35671`, so the
published values differ by only 0.75 and 0.63 refit standard errors. The
pressure RMSE is 2.0843 GPa. This is classified as `similar`, not exact parity,
because the journal does not publish unrounded states, row-wise uncertainties,
regression weights, the numerical integration grid, or the BM3 residual
definition.

## Primary-source chain

The final source is Noguchi et al., *Equation of state of NiO studied by shock
compression*, *Journal of Physics and Chemistry of Solids* **60**, 509--514
(1999),
[doi:10.1016/S0022-3697(98)00296-0](https://doi.org/10.1016/S0022-3697(98)00296-0).
Its Table 1 supplies the eight final shock states through 147.6 GPa, and its
abstract reports the final BM3 coefficients.

The same team's openly readable conference paper, Noguchi et al., *Shock
Compression of NiO to 130 GPa*, *Review of High Pressure Science and
Technology* **7**, 832--834 (1998), supplies the thermodynamic details omitted
from the journal abstract. Pages 833--834 state:

- the continuous relation `Us = 5.36 + 1.19 up`, with velocities in km/s;
- `gamma0 = 1.38`, calculated from `2s = gamma + 1`;
- Debye temperature `theta0 = 390 K`, calculated from sound velocity;
- `gamma/V = constant`; and
- reduction to a 300 K isotherm followed by a Murnaghan--Birch fit.

The conference paper is available from
[J-STAGE](https://www.jstage.jst.go.jp/article/jshpreview1992/7/0/7_0_832/_pdf).
It covers the first seven journal rows, ending at 132.9 GPa, and reports
`K0 = 184 +/- 5 GPa` with `K0' = 4` plus a shock temperature of about 1600 K
at the maximum pressure. The eighth journal experiment extends the range to
147.6 GPa.

## Exact reduction convention

Let `V0` be the initial molar volume calculated from the measured bulk density
`rho0 = 6.781 g/cm3` and the NiO molar mass. The bulk density, not the separately
reported X-ray density, belongs in the Rankine--Hugoniot relations. For the
linear source Hugoniot, Peritheos calculates `PH(V)`, `up(V)`, and `Us(V)`.
The shock energy is

\[
E_H(V)-E_0=\frac{1}{2}\,[P_H(V)+P_0]\,(V_0-V).
\]

The temperature along that path follows directly from the first law:

\[
\frac{dT_H}{dV}=
\frac{dE_H/dV+P_H}{C_V(V,T_H)}-\frac{\gamma(V)T_H}{V}.
\]

The heat capacity is the Debye heat capacity for two atoms per NiO formula
unit. In Peritheos notation `gamma = gamma0*(V/V0)^q`; therefore the paper's
`gamma/V = constant` statement is exactly `q = 1`. Thermodynamic integration
then gives

\[
\theta(V)=\theta_0\exp[\gamma_0(1-V/V_0)].
\]

At every measured final-state volume the 300 K pressure is

\[
P_{300}(V)=P_H^{obs}(V)-\frac{\gamma(V)}{V}
\left[E_D(V,T_H)-E_D(V,300\,\mathrm{K})\right].
\]

The final BM3 regression uses these reduced pressures. This is distinct from
using the raw `PH` column as the fit dependent variable.

## Reconstructed states

| `V/V0` | `PH` (GPa) | `TH` (K) | thermal correction (GPa) | `P300` (GPa) |
|---:|---:|---:|---:|---:|
| 0.930 | 17.7 | 338.324 | 0.219097 | 17.480903 |
| 0.892 | 27.8 | 380.923 | 0.463220 | 27.336780 |
| 0.853 | 42.7 | 460.857 | 0.927298 | 41.772702 |
| 0.766 | 86.6 | 919.197 | 3.676934 | 82.923066 |
| 0.739 | 108.1 | 1203.970 | 5.417770 | 102.682230 |
| 0.722 | 121.2 | 1438.721 | 6.860625 | 114.339375 |
| 0.711 | 132.9 | 1618.660 | 7.969638 | 124.930362 |
| 0.702 | 147.6 | 1784.716 | 8.994801 | 138.605199 |

The decisive method check uses the first seven rows and fixes `K0'=4`, exactly
as the conference paper did. It returns `K0=184.29745 GPa` and
`T_H(132.9 GPa)=1618.66 K`, reproducing both independent published checkpoints.

## Implementation and limits

Run the audit with:

```text
uv run python scripts/reproduce_noguchi_1999_nio.py
```

The implementation uses the existing `LinearUsUpHugoniot`,
`MieGruneisenDebye`, and `BM3` models. No new EOS term is required: the existing
integrated-Gruneisen Debye law with `q=1` is the source's constant-`gamma/V`
case. Tests pin the thermal convention, every reduction direction, the two
conference checkpoints, and the qualified eight-row fit. The source-reported
`191 GPa, 3.9` parameters remain the production record; the rounded-table refit
is validation evidence and does not replace them.
