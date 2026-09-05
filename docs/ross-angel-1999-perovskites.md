# Ross and Angel (1999) CaTiO3 and CaGeO3 perovskites

Peritheos includes two published room-temperature third-order
Birch-Murnaghan records from Ross and Angel, *American Mineralogist* **84**,
277-281 (1999), DOI
[`10.2138/am-1999-0309`](https://doi.org/10.2138/am-1999-0309):

- `calcium_titanate_perovskite_ross_1999_bm3_1`
- `calcium_germanate_perovskite_ross_1999_bm3_1`

Both records use the standard Eulerian finite-strain BM3 pressure equation,

\[
P(V)=\frac{3K_0}{2}\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_0-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\}.
\]

The equation mapping is independently checked against all observations: this
form reproduces the paper's reported maximum pressure residuals and weighted
chi-square values. All three coefficients were fitted; neither fit fixed a
coefficient. The alternative `K0_prime = 4` diagnostics printed by the paper
are not stored as EOS records because they are rejected comparison fits.

## Parameters and scope

| Record | `V0` (A^3) | `K0` (GPa) | `K0_prime` | Table 1 range (GPa) | Rows |
|---|---:|---:|---:|---:|---:|
| CaTiO3 | 223.764(17) | 170.9(14) | 6.6(3) | 0.001-9.700 | 11 |
| CaGeO3 | 206.490(17) | 194.0(21) | 6.1(5) | 0.0001-8.553 | 9 |

Here `V0` is one conventional orthorhombic Pbnm unit cell containing four
formula units (`Z = 4`), not a formula-unit or molar volume. The structures
are GdFeO3-type perovskites in space group Pbnm (No. 62). The material files
use the near-ambient lattice of each EOS specimen and independently sourced
atomic coordinates. The numerical 298 K reference is an operational encoding
of the article's qualitative “room temperature” statement; it is not a
reported thermometer reading.

The reported coefficient errors are retained, but the paper gives neither
their confidence convention nor a parameter covariance matrix. Covariance is
therefore stored as missing. No phase transition was observed over either
fit interval. CaGeO3 transition behavior discussed above the measured range
does not extend the record's validity.

## Observations and pressure calibration

The complete printed Table 1 observations are bundled as:

- `calcium-titanate-perovskite-ross-1999-table1-compression.csv`
- `calcium-germanate-perovskite-ross-1999-table1-compression.csv`

They preserve pressure, all three lattice parameters, conventional-cell
volume, and every printed estimated standard deviation. The near-ambient
pressure values have no printed pressure error and are represented as missing,
not zero.

Pressure was determined from a single-crystal alpha-quartz internal standard
using Angel et al., *Journal of Applied Crystallography* **30**, 461-466
(1997), DOI
[`10.1107/S0021889897000861`](https://doi.org/10.1107/S0021889897000861).
That calibration is the room-temperature BM3 with `V0 = 112.981(2) A^3`,
`K0 = 37.12(9) GPa`, and `K0_prime = 5.99(4)`. It is linked through the
executable Peritheos record `alpha_quartz_angel_1997_bm3_1`.

The calibration equation is fully resolved, but observation-level pressure
recalculation is not possible: Ross and Angel print the reduced pressures and
sample cells, not the row-wise measured quartz volumes. The experiment used a
4:1 methanol-ethanol pressure medium and a Huber four-circle diffractometer.
The publisher article exposes no supplementary file.

## Numerical audit

Direct evaluation of the printed coefficients gives pressure RMSE / maximum
absolute residuals of 0.018946 / 0.029679 GPa for CaTiO3 and
0.019111 / 0.030450 GPa for CaGeO3. These reproduce the paper's rounded
0.03 GPa maximum residual for each material.

An independent effective-variance refit combines the printed pressure errors
with volume errors propagated along the model curve. It recovers:

| Material | Refit `V0` | Refit `K0` | Refit `K0_prime` | Reduced chi-square |
|---|---:|---:|---:|---:|
| CaTiO3 | 223.76464009 | 170.85581942 | 6.57739326 | 1.3808862 |
| CaGeO3 | 206.48942690 | 194.03279818 | 6.09858913 | 1.1028596 |

Every refitted coefficient lies within its printed uncertainty, and the
reduced chi-square values round to the source values 1.4 and 1.1. These are
validation results for the published records, not replacement refit records.
