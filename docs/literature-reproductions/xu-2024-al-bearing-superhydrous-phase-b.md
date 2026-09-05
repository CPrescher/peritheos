# Xu et al. (2024): Al-bearing superhydrous phase B

## Audit outcome

This audit covers all ten LitCurate rows grouped under Chaowen Xu et al.,
*Effect of Al-Incorporation on the Sound Velocities of Superhydrous Phase B at
High Pressure and High Temperature*, **Geophysical Research Letters 51**,
e2023GL107818 (2024),
[doi:10.1029/2023GL107818](https://doi.org/10.1029/2023GL107818).

The source bundle contributes **zero production EOS records**. Four reported
Birch-Murnaghan parameter sets cannot be reproduced from the author's archived
P-V-T rows, three finite-strain rows are coupled acoustic-elasticity models
rather than supported Peritheos P(V) equations, and three rows are comparison
values copied from earlier publications. No `.eosmat` record, dataset CSV,
reproduction script, or test is added for this paper.

LitCurate was used only to identify the candidates. The scientific audit used
the final open publisher article and the corresponding author's public
EarthChem archive,
[doi:10.60520/IEDA/113199](https://doi.org/10.60520/IEDA/113199). The publisher
also registers `2023GL107818-sup-0001-Supporting Information SI-S01.pdf` as
Supporting Information S1. Its Wiley delivery endpoint returned HTTP 403 in the
noninteractive audit environment, so this report does not claim to have read
Text S2 or silently reconstruct its equations.

## Author archive and source integrity

EarthChem record 3199 exposes `3199-1_tabulateddataset.xlsx`. The repository
publishes SHA-1 `a286d63d36735613817056b3297b73b95e717b08`; the retrieved file
matches that checksum and has SHA-256
`b5555ea3f6f2db766c92efc47b08c8786a5caefcd83fe65a7326e1dd804e90b3`.
It is licensed CC-BY-SA-4.0. The archive is a new version of EarthChem
[doi:10.60520/IEDA/113160](https://doi.org/10.60520/IEDA/113160); the two
versions contain the same numerical P-T-V blocks apart from binary floating
point serialization, while the newer version adds the chemical-formula and
ambient-condition sheets.

The workbook contains:

- 33 P-T-V-VP-VS-density rows for sample S3358, spanning 11.4--23.5 GPa and
  300--900 K, including seven 300 K rows;
- 12 rows for sample S3359, spanning 14.82--22.06 GPa and 300--1300 K,
  including three 300 K rows;
- one ambient cell and one analyzed formula for each sample; and
- geotherm, aggregate-velocity, and mineral-proportion tables used for the
  geophysical calculations.

The archive is not internally sufficient to reproduce the published EOS fits.
It provides a single row-wise pressure column without identifying each value as
gold-scale or pressure-scale-free, and it does not include row-wise Au volumes.
Thus the Tsuchiya (2003) gold pressure scale cannot be independently rerun. The
deposited density column also contains nonmonotonic values and exact values
repeated between apparently unrelated rows, preventing it from serving as an
independent reconstruction of the scale-free acoustic pressure route.

Two smaller transcription issues are retained explicitly:

- In `Table 2. ambient condition`, both rows are labeled S3358. The article
  assigns the second row, 621.69(17) A^3, to S3359.
- For the S3359 pressure-scale-free BM fit, Section 3.1 reports
  K'~T0~ = 4.39(2), whereas Table 1 reports 4.38(2). LitCurate stores the table
  value. For the gold-scale fit, the prose reports K~T0~ = 136.6(2) GPa while
  Table 1 rounds it to 137(2) GPa; LitCurate stores the table value.

## Samples, structure, and volume basis

Both specimens are homogeneous, sintered polycrystals synthesized at 20 GPa
and 1373 K for one hour. The paper gives the EDS/mass-balance compositions
Mg~8.93~Si~2.72~Al~0.56~O~18~H~5.57~ for S3358 and
Mg~8.99~Si~2.77~Al~0.54~O~18~H~5.25~ for S3359. Hydrogen is obtained by mass
balance rather than direct site refinement.

The article assigns both samples to orthorhombic Pnn2, space-group number 34,
and reports the following ambient conventional cells:

| sample | a (A) | b (A) | c (A) | V (A^3) |
|---|---:|---:|---:|---:|
| S3358 | 5.0662(9) | 14.0549(23) | 8.7269(12) | 621.41(17) |
| S3359 | 5.0664(9) | 14.0590(24) | 8.7281(28) | 621.69(17) |

The present article does not print Z. The established superhydrous-phase-B
conventional cell contains two O~18~ formula units (Z = 2); this is also the
only physically meaningful basis for a roughly 621 A^3 cell. Because no record
is promoted, the inferred Z is documented rather than used to normalize an
executable EOS.

## Published parameterizations

The article's Table 1 separates three model groups.

### Gold-scale Birch-Murnaghan fits

The paper says that room-temperature volumes and pressures calculated from the
Tsuchiya (2003) Au standard were fit with EosFit to a third-order
Birch-Murnaghan equation. The reported coefficients are:

| sample | stated P-T range | V~0~ (A^3) | K~T0~ (GPa) | K'~T0~ | other thermal terms |
|---|---|---:|---:|---:|---|
| S3358 | 11.4--23.5 GPa, 300--900 K | 621.41 fixed | 137.9(3) | 4.43(4) | dK~T~/dT = -0.018(8) GPa/K; alpha~0~ = 4.3(5) x 10^-5^ K^-1^ |
| S3359 | 14.8--21.2 GPa, 300 K | 621.69 fixed | 137(2) in Table 1; 136.6(2) in prose | 4.43 fixed | dK~T~/dT = -0.018 fixed; alpha~0~ = 4.6(2) x 10^-5^ K^-1^ |

The S3358 row is internally described as both a room-temperature compression
fit and a 300--900 K thermoelastic fit. The S3359 stated interval also conflicts
with the archive: its only three 300 K rows span 14.82--19.89 GPa, not
14.8--21.2 GPa.

### Pressure-scale-free Birch-Murnaghan fits

The authors recalculate pressure from their acoustic finite-strain formalism
and refit the volumes with BM3. These are global thermal fits, not standalone
300 K isotherms:

| sample | stated P-T range | V~0~ (A^3) | K~T0~ (GPa) | K'~T0~ | dK~T~/dT (GPa/K) | alpha~0~ (K^-1^) |
|---|---|---:|---:|---:|---:|---:|
| S3358 | 11.4--23.5 GPa, 300--900 K | 621.41 fixed | 145.07(3) | 4.455(4) | -0.0140(1) | 4.215(7) x 10^-5^ |
| S3359 | 14.8--22.1 GPa, 300--1300 K | 621.69 fixed | 147.7(1) | 4.38(2) table; 4.39(2) prose | -0.0166(2) | 4.25(1) x 10^-5^ |

### Acoustic Eulerian finite-strain fits

The remaining three source-reported candidates fit adiabatic bulk and shear
moduli to Davies-Dziewonski Eulerian finite-strain functions. Table 1 reports
K~S0~, G~0~, their pressure derivatives, and, for the global models,
temperature derivatives and expansivity. These are coupled elasticity models.
They are not interchangeable with a static Birch-Murnaghan P(V) curve merely
because K~S0~ and K'~S~ resemble volumetric EOS coefficients.

| sample/model | K~S0~ (GPa) | K'~S~ | G~0~ (GPa) | G' | dK~S~/dT (GPa/K) | dG/dT (GPa/K) |
|---|---:|---:|---:|---:|---:|---:|
| S3358, 300 K | 148(3) | 4.3(3) | 94(2) | 1.8(8) | -- | -- |
| S3358, 300--900 K | 147(2) | 4.4(1) | 93.9(9) | 1.8(3) | -0.017(1) | -0.019(1) |
| S3359, 300--1300 K | 150(7) | 4.2(6) | 93(4) | 2.0(2) | -0.019(2) | -0.022(2) |

Peritheos has no dedicated coupled acoustic-elasticity record type for this
formalism, and the defining equations remain in inaccessible Supporting Text
S2. Encoding these rows as BM3 would change both the observable and the model.

## Numerical non-reproduction

The standard BM3 pressure equation used for the independent checks was

\[
P(V)=\frac{3K_0}{2}\left[\eta^7-\eta^5\right]
\left\{1+\frac{3}{4}(K'_0-4)[\eta^2-1]\right\},
\qquad \eta=(V_0/V)^{1/3}.
\]

The exact deposited 300 K P-V rows are:

| sample | deposited (P in GPa, V in A^3) pairs |
|---|---|
| S3358 | (11.4, 584.43), (14.5, 574.97), (16.5, 570.06), (18.2, 565.96), (19.7, 561.53), (21.1, 557.89), (22.3, 555.16) |
| S3359 | (14.82, 573.64), (17.55, 566.60), (19.89, 560.84) |

An unweighted pressure-residual refit of the seven S3358 rows, fixing
V~0~ = 621.41 A^3 and varying K~0~ and K'~0~, gives
K~0~ = 170.33642731 GPa and K'~0~ = 2.72989597 with RMSE
0.10205447 GPa. The published gold-scale curve has RMSE 2.10477320 GPa, and
the printed pressure-scale-free reference coefficients have RMSE
1.26700170 GPa. Neither reported parameter pair reproduces the deposit.

For S3359, fixing V~0~ = 621.69 A^3 and the gold-fit K'~0~ = 4.43 while
varying K~0~ gives K~0~ = 153.82068051 GPa and RMSE 0.02385243 GPa. The
published prose value K~0~ = 136.6 GPa gives RMSE 1.96410876 GPa. Allowing
both K~0~ and K'~0~ to vary gives K~0~ = 155.97491844 GPa,
K'~0~ = 4.14012581, and RMSE 0.00471605 GPa. The printed scale-free
coefficients K~0~ = 147.7 GPa and K'~0~ = 4.39 give RMSE 0.73044582 GPa.

The thermal rows do not resolve the discrepancy. Using
K~0~(T) = K~0~ + (dK~0~/dT)(T-300) and either linear or exponential
V~0~(T) with the printed alpha~0~ gives published-curve RMSE values of
1.06272487 and 1.04200217 GPa for S3358, and 0.82029862 and
0.72749756 GPa for S3359. For comparison, a four-parameter exponential-volume
refit to every deposited row gives:

| sample | K~0~ (GPa) | K'~0~ | dK~0~/dT (GPa/K) | alpha~0~ (K^-1^) | RMSE (GPa) |
|---|---:|---:|---:|---:|---:|
| S3358 | 170.257426 | 2.666456 | -0.0120776 | 3.78927 x 10^-5^ | 0.103805 |
| S3359 | 157.899677 | 3.898437 | -0.0172138 | 4.19526 x 10^-5^ | 0.107099 |

These are diagnostics only, not replacement EOS coefficients. The purpose of
the refits is to demonstrate that ordinary rounding, weighting, or the choice
between linear and exponential reference-volume expansion cannot explain the
published/archive mismatch.

## Disposition of every same-DOI LitCurate row

| LitCurate identifier | source row | reported candidate | disposition | reason |
|---|---:|---|---|---|
| `litcurate_cd5dc7c1140d15bf` | 1242 | S3358 Au-scale BM3: `V0=621.41`, `K0=137.9`, `K0'=4.43` | **held / numerically non-reproducible** | The seven deposited 300 K P-V points disagree with the printed curve by 2.105 GPa RMSE. The deposit lacks Au volumes needed to recompute pressure. |
| `litcurate_49be30bb0315f2ed` | 1243 | S3359 Au-scale BM3: `V0=621.69`, `K0=137`, fixed `K0'=4.43` | **held / numerically non-reproducible** | The prose reports 136.6 rather than the table's rounded 137; only three archived 300 K rows exist, and they disagree with the printed curve by 1.964 GPa RMSE. |
| `litcurate_803c07c856498491` | 1244 | S3358 scale-free BM3: `V0=621.41`, `K0=145.07`, `K0'=4.455` | **held / incomplete thermal mapping and non-reproducible** | This is a global 300--900 K model with additional coefficients omitted by the candidate. The complete printed thermal curve still misses the archive by about 1.04 GPa RMSE. |
| `litcurate_6eeacb673fc8fb75` | 1245 | S3359 scale-free BM3: `V0=621.69`, `K0=147.7`, `K0'=4.38` | **held / incomplete thermal mapping and non-reproducible** | This is a global 300--1300 K model. Table 1 and prose disagree on K'~0~, and the complete printed curve misses the archive by at least 0.73 GPa RMSE. |
| `litcurate_403cfa12a19f7103` | 1246 | S3358 room-T finite strain: `KS0=148`, `KS'=4.3` | **held / unsupported coupled elasticity model** | The source fits K~S~ and G simultaneously. It is not a BM P(V) equation and cannot be reduced to one without changing the reported model. |
| `litcurate_2b71322686e77b37` | 1247 | S3358 global finite strain: `KS0=147`, `KS'=4.4` | **held / unsupported coupled elasticity model** | The model includes G~0~, G', both thermal derivatives, and expansivity; the candidate is not a complete executable representation. |
| `litcurate_91a64bb6a2815e53` | 1248 | S3359 global finite strain: `KS0=150`, `KS'=4.2` | **held / unsupported coupled elasticity model** | Same coupled acoustic formalism; the exact equations are delegated to Supporting Text S2 and do not map to a current Peritheos EOS family. |
| `litcurate_e7123de0ecb3ae27` | 1249 | Al-free SuB high-pressure XRD comparison: `K0=132.2` | **rejected as a Xu-source record / citation trace only** | The row lacks V~0~, K'~0~, and an equation identity and summarizes earlier XRD publications. Audit the underlying primary source separately. |
| `litcurate_dcac4085a9de0eeb` | 1250 | SuB ultrasonic/Brillouin comparison: `K0=140.3` | **rejected as a Xu-source record / citation trace only** | It is an incomplete comparison value from earlier elasticity studies, not a Xu et al. EOS fit. |
| `litcurate_b51dd09b4e7b9fc9` | 1251 | SuB first-principles comparison: `K0=154` | **rejected as a Xu-source record / citation trace only** | It is an incomplete lower endpoint of a 154.0--161.8 GPa literature range and has no V~0~, K'~0~, or equation identity. |

## Duplicate and future-work check

No bundled material or primary-source-audit entry uses DOI
`10.1029/2023GL107818`, and no existing material record duplicates either
sample composition. The paper could eventually yield two materials and up to
four BM-family records, but only after a source correction or an unambiguous
row-wise data release reconciles the archived pressures and volumes with the
published coefficient sets. The three acoustic fits additionally require a
dedicated, tested coupled-elasticity framework and the exact Supporting Text S2
equations. Until those conditions are met, this audit remains documentation
only with zero production records.
