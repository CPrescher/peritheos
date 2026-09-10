# Redfern et al. (1993): natural magnesite

Primary source: S. A. T. Redfern, B. J. Wood, and C. M. B. Henderson,
“Static compressibility of magnesite to 20 GPa: Implications for MgCO3 in the
lower mantle,” *Geophysical Research Letters* **20**, 2099–2102 (1993),
<https://doi.org/10.1029/93GL02507>.

The original four-page article resolves the earlier data-access limitation.
Table 1 publishes 18 room-temperature static-compression measurements for
natural `Mg0.991Fe0.008Mn0.001CO3` (Harwood specimen 2212): pressure, hexagonal
`a` and `c`, conventional-cell volume, and parenthetical uncertainties. These
rows are transcribed directly in
[`peritheos/data/datasets/magnesite-redfern-1993-table1-pv.csv`](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/magnesite-redfern-1993-table1-pv.csv).
No figure digitization is used.

## Pressure calibration

The Experimental Methods state that powdered NaCl was mixed with the sample as
the internal standard and pressure calibrant. Pressure was calculated from a
least-squares refinement of at least three NaCl peaks using the Decker (1971)
NaCl EOS, <https://doi.org/10.1063/1.1660714>. The reported pressure uncertainty
is ±2%. Table 1 gives the calibrated pressures but not the corresponding NaCl
cell volumes, so the pressure scale is identified while row-wise calibration
recalculation remains unavailable.

## Source fit protocols

The ambient Table 1 value `V0=279.4(2) Å3` defines `V/V0` and the Eulerian
strain. It is therefore treated as a measured, fixed reference value in both
reproductions rather than as a third free EOS coefficient.

The reduced Birch result is described as a least-squares fit with `K0'=4`. The
paper does not state residual direction or weighting for this fit. The
reproduction uses all 18 Table 1 rows and unweighted pressure residuals, with
`V0` and `K0'` fixed.

For the full third-order result, Equation (2) and the Figure 2 discussion state
that only observations with Eulerian strain `fV>0.01` were analyzed. This selects
nine Table 1 rows from 4.6 to 19.7 GPa. The caption calls the calculation a
“weighted regression of the volume data”; the reproduction therefore minimizes
volume residuals using the printed one-sigma volume uncertainties, with `V0`
fixed and `K0`, `K0'` free. The paper does not provide more detailed weighting
or covariance information.

## Direct refit results

| Fit | Published | Refit | Assessment |
|---|---|---|---|
| Reduced Birch, 18 rows | `K0=142(9) GPa`, `K0'=4` | `K0=142.8547 GPa`, pressure RMSE `0.5223 GPa` | parity within the published uncertainty |
| Full BM3, 9 rows | `K0=151(7) GPa`, `K0'=2.5` | `K0=150.7746 GPa`, `K0'=2.7249`, pressure RMSE `0.6133 GPa` | similar; consistent with the published fit, but strict uncertainty parity is unavailable for `K0'` |

The full-fit reproduction gives one-sigma numerical errors of approximately
`7.16 GPa` for `K0` and `1.34` for `K0'`, and a weighted-volume reduced
chi-square of `1.668`. The source itself prints no uncertainty for `K0'`.
The large refit uncertainty shows that `K0'` is weakly constrained, so the
`similar` classification records limited precision rather than a meaningful
disagreement with the published value.

The executable audit is `scripts/reproduce_redfern_1993_magnesite.py`; the
complete machine-readable provenance and checksums are in
[`docs/data/redfern-1993-magnesite-audit.json`](../data/redfern-1993-magnesite-audit.json).

## Secondary-source inconsistency

Yu et al. (2024), Table 1, marks the Redfern `K0'=2.5` value as fixed. The
original paper instead calls this the full Birch–Murnaghan result, derives it
from the slope of Equation (2), and reports `K0'=2.5`. The production record
therefore correctly treats `K0'` as free and only `V0` as fixed.

## LitCurate disposition

| Row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 52 | `litcurate_71122768f103ee77` | ACCEPT | Complete fixed-derivative reduced Birch fit. |
| 53 | `litcurate_bd03302a0d214f4a` | ACCEPT | Complete full third-order Birch–Murnaghan fit. |
| 54 | `litcurate_383afd63c513fd64` | REJECT | The `148/3.2` values come from the separate second-order polynomial in pressure, not the Birch–Murnaghan fit. |
| 55 | `litcurate_b94b146ddc4a305e` | REJECT | The MgO value is adopted comparative input, not a magnesite measurement. |

Result: **2 production records, 2 rejected rows, and 18 authoritative Table 1
observations supporting direct refits**.
