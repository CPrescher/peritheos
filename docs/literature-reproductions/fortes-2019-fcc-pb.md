# Fortes (2019) fcc-Pb pressure-scale audit

## Outcome

The catalog's 300 K BM4 slice is an exact executable slice of the analytical
pressure surface in Fortes (2019), but the reported coefficient fit cannot be
independently reproduced. This is now a source-specific, evidence-backed
`not_refittable` result rather than a generic "parameterization only" label.

The report is open under CC BY 4.0 and has canonical DOI
[10.5286/raltr.2019002](https://doi.org/10.5286/raltr.2019002). It publishes the
complete pressure equation and all nine coefficients. It does **not** publish
the numerical compression dataset used in the final least-squares step. The
decisive high-temperature series was supplied privately by Alexei Kuznetsov;
Kuznetsov et al. (2002) publishes only three fcc transition anchors and a plot
in which the other high-temperature observations have already been reduced to
room temperature.

Two source-published isothermal records are cataloged for fcc Pb:

- `lead_fcc_kuznetsov_2002_bm3_2` is the 296 K BM3 slice of the earlier
  Kuznetsov P-V-T surface (`V0=121.228` Å³ per four-atom cell, `K0=40.5` GPa,
  `K0'=5.74`). It is appropriate for the room-temperature fcc branch: the paper
  observes fcc Pb at ambient temperature, reports transformation onset near
  12 GPa, and prints a 296 K fcc transition anchor at 13.1 GPa.
- `lead_fcc_fortes_2019_bm4_1` remains the default and is the 300 K BM4 slice
  of the revised Fortes surface.

Neither slice is labeled as an independent Peritheos refit. The complete
thermal models remain executable in the audit script, while their polynomial
temperature dependence is outside the current `.eosmat` thermal schema.

The executable audit is
`scripts/audit_fortes_2019_fcc_pb.py`. The three exact public anchors are
bundled as
`peritheos/data/datasets/lead-fcc-kuznetsov-2002-table1-transition-pvt.csv`.
They are retained as a partial source check, not relabeled as the private fit
table.

## Published equation and reference state

Fortes uses fourth-order Eulerian finite strain,

`f = 0.5 [(V/V0)^(-2/3) - 1]`,

`P = 3 K0 f (1+2f)^(5/2) {1 + 3/2 (K0'-4)f + 3/2 [K0 K0'' + (K0'-4)(K0'-3) + 35/9]f^2}`.

With `Delta T = T - 300 K`, the temperature corrections are

`V0(T) = 121.418 + 1.058e-2 Delta T + 3.5e-6 Delta T^2` Å³ per fcc cell,

`K0(T) = 41.73 - 2.544e-5 Delta T - 2.8e-6 Delta T^2` GPa,

`K0'(T) = 5.39 + 1.1e-3 Delta T`, and `K0'' = -0.33 GPa^-1`.

Table 1 labels the polynomial constant terms `V0(0)`, `K0(0)`, and `K'(0)`.
Here zero is the polynomial argument `Delta T=0`, hence 300 K; it is not an
absolute-zero reference state. At 300 K the full surface therefore reduces
exactly to the stored isothermal BM4 record.

## Staged construction

Fortes did not fit all nine coefficients simultaneously.

1. Compiled 10--600 K expansion data were fitted with a Debye internal-energy
   model (Equations 6--7). The exact Debye curve above 100 K was then replaced
   with the printed quadratic `V0(T)`. The report does not tabulate the input
   rows, a complete inclusion mask, or weights.
2. Literature elastic constants supplied `KS=(c11+2c12)/3`. Equation 8 converts
   these to `KT` using the stage-one `V(T)` and expansivity plus a cubic `CP(T)`
   fit. Fortes identifies the selected heat capacities (Meads et al. below room
   temperature and Leadbetter above), but does not print the cubic coefficients
   or input table. The resulting `K0(T)` quadratic is printed completely.
3. Holding `V0(300)`, `a`, `b`, `K0(300)`, `c`, and `d` fixed, the compression
   fit varied `K0'(300)`, `e=dK0'/dT`, and temperature-invariant `K0''` over
   295--788 K.

The material record now preserves that fixed/free split. The parenthetical
errors in Table 1 are retained, but the report gives no covariance matrix and
does not say how the errors were scaled.

## Compression inputs and pressure scales

Section 2.3 names a heterogeneous compilation: Bridgman's mechanical
compression work, shock reductions, Vaidya and Kennedy static compression,
Mao and Bell values shown under both ruby and NaCl calibrations, Mao et al.
(1990), and Kuznetsov et al. (2002). Figure 4 is explicitly only a *selection*
of room-temperature data, not a declaration of the fitted row mask.

The pressure-scale chain is therefore only partially resolved. In particular,
Kuznetsov used NaCl as both medium and calibrant with the Brown (1999) P-V-T
scale ([doi:10.1063/1.371596](https://doi.org/10.1063/1.371596)). The calibration
identity is recoverable, but the row-wise NaCl volumes needed to recalculate
the pressures are not. Fortes does not state that the heterogeneous source
pressures were recalculated to one common scale, nor does it provide the
calibrant observations required to do so. No uniform ruby, NaCl, or other
calibration is inferred in the catalog.

Kuznetsov's own published predecessor model is nevertheless fully executable.
Equations 1--3 and Table 2 give a thermal BM3 referenced to 296 K, with
`V0=30.307 Å³/atom`, `K0=40.5 GPa`, `K0'=5.74`, cubic thermal expansion,
quadratic `K0(T)`, and linear `K0'(T)`. The audit implements all nine printed
fcc coefficients and evaluates them against the three Table 1 transition
anchors (about 0.46 GPa RMSE). This validates the recoverable source model; it
does not turn its fitted curve into the experimental rows shared privately
with Fortes.

The hcp column of the same table is independently executable as well. At
296 K it reduces to BM3 with `V0=29.908 Å³/atom` (`59.816 Å³` per two-atom hcp
cell), `K0=54.2 GPa`, and `K0'=3.61`. Visual inspection of Table 2 confirms the
thermal-coefficient signs that are easily lost by text extraction:
`b1=-5.46e-2 GPa/K`, `b2=+2.34e-5 GPa/K²`, and
`dK0'/dT=+5.01e-3 K^-1`.

The full hcp surface predicts `14.035`, `13.934`, and `13.510 GPa` for the three
rounded Table 1 hcp transition states, an RMSE of `0.754 GPa`. More decisively
for the proposed 296 K slice, it reproduces the seven solid (direct
room-temperature) hcp markers digitized from Figure 3 with `0.266 GPa` RMSE.
Across all 21 solid and source-reduced open markers the RMSE is `0.642 GPa`,
within the limitations of plot digitization, marker overlap, and the source's
unpublished high-temperature reduction details. This verifies the published
hcp slice but does not constitute a new refit. The verified slice is cataloged
as the published record `lead_hcp_kuznetsov_2002_bm3_3`; the two Dewaele (2019)
hcp Vinet pressure-scale reductions remain separate records.

## Fit attempt and why it is not parity

Kuznetsov Table 1 prints three fcc states at the transition boundary:

| T (K) | P (GPa) | atomic V (Å³) | fcc-cell V (Å³) |
|---:|---:|---:|---:|
| 296 | 13.1 | 25.12 | 100.48 |
| 402 | 13.9 | 25.14 | 100.56 |
| 469 | 12.6 | 25.26 | 101.04 |

The audit evaluates the published Fortes coefficients at those states and
runs an unweighted pressure-residual diagnostic with the six preliminary
thermal coefficients fixed. The printed Fortes curve has an RMSE of about
0.81 GPa on these three rounded boundary points. Allowing the three published
compression-stage coefficients to vary lowers the RMSE only to about 0.38 GPa
and produces a severely ill-conditioned Jacobian (condition number above
`10^9`) with coefficients far from Table 1. This is expected: three rounded
phase-boundary checkpoints cannot stand in for the larger private compression
series.

Figure 4 cannot rescue the fit. Every plotted observation there is represented
at room temperature, where `Delta T=0`; analytically `dP/de=0`. Thus the plot
has rank at most two for the three reported free coefficients and cannot
determine `e`. Kuznetsov Figure 3 likewise discards the original temperatures
by reducing open high-temperature markers to room temperature.

An independent reported fit therefore still requires all of the following:

- the numerical Kuznetsov P-V-T table supplied privately to Fortes;
- the exhaustive source/row selection and any phase or quality exclusions;
- the treatment of the heterogeneous source pressure scales;
- the least-squares dependent variable and observation weights; and
- the covariance and error-scaling procedure.

Until those inputs are recovered, the compiled coefficients remain a valid
source-reported pressure surface but not an independently refittable one.
