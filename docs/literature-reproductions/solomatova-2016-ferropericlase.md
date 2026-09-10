# Solomatova et al. (2016): spin-crossover ferropericlase branches

Primary source: N. V. Solomatova, J. M. Jackson, W. Sturhahn, J. K. Wicks,
J. Zhao, T. S. Toellner, B. Kalkan, and W. M. Steinhardt, “Equation of state and
spin crossover of (Mg,Fe)O at high pressure, with implications for explaining
topographic relief at the core-mantle boundary,” *American Mineralogist* **101**,
1084–1093 (2016), <https://doi.org/10.2138/am-2016-5510>.

## Result

The 16 Table 7 records are paired high-spin (HS) and low-spin (LS) BM3
**reference branches** from eight coupled MINUTI crossover fits. They are not
independent pressure-window fits and neither branch is the physical mixed-spin
curve. This audit concerns the 14 branches fitted to seven earlier studies; the
two Fp48 branches fitted to Solomatova's own observations were already backed
by all 45 Table 1 rows.

Authoritative earlier observations were recovered for Marquardt et al. (29
table rows), Chen et al. (36 table rows), Fei et al. (77 auxiliary-table rows),
and Zhuravlev et al. (71 table rows). Existing Mao et al. digitization provides
42 plotted observations, and the Lin et al. digitization provides 43
normalized-volume observations. Ten of the 14 earlier-study branches now have
executable absolute-volume reconstructions: Marquardt, Chen, Fei, and
Zhuravlev HS/LS are `similar`; Mao HS/LS are retained as
`parity_not_achieved`. The remaining four are `not_refittable`: the two Lin
Mg0.40Fe0.60O branches have no P–V series, while the two Lin Mg0.83Fe0.17O
branches have real normalized observations that cannot independently identify
absolute V0,HS and V0,LS.

## Trace to the earlier observations

| Solomatova branch pair | Earlier primary source | Composition and assignment | Recovered evidence | Audit outcome |
|---|---|---|---|---|
| Mg0.90Fe0.10O HS/LS | Marquardt et al. (2009b), DOI 10.1016/j.epsl.2009.08.017 | Nominal Mg0.90Fe0.10O; Brillouin-constrained HS reference and inferred LS reference | All 29 P–V observations printed in EPSL Table 2, including the source HS/mixed/LS assignment. | Executable; `similar`. The K0S prior remains a constraint, not an observation. |
| Mg0.83Fe0.17O HS/LS | Lin et al. (2005), DOI 10.1038/nature03825 | Mg0.83Fe0.17O; Figure 2 labels HS at ≤55 GPa and LS at ≥74 GPa | All 43 Figure 2 markers recovered as P and V/V0,HS. | `not_refittable` for absolute branches; Table-7-anchored shape diagnostic only. |
| Mg0.75Fe0.25O HS/LS | Mao et al. (2011), DOI 10.1029/2011GL049915 | Mg0.75Fe0.25O; all plotted 300 K states enter the coupled crossover fit | Existing 42-marker Figure 1 vector digitization. | Executable; `parity_not_achieved`. |
| Mg0.65Fe0.35O HS/LS | Chen et al. (2012), DOI 10.1029/2012JB009162 | Catalog/fit label Mg0.65Fe0.35O; measured specimen (Mg0.654(1)Fe0.330(1)Ti0.016(1))O | Complete 36-row Table 1 compression series. | Executable; `similar`. |
| Mg0.61Fe0.39O HS/LS, Fei data | Fei et al. (2007), DOI 10.1029/2007GL030712 | Mg0.61Fe0.39O; source distinguishes compression from decompression | Complete Table S2: 38 compression plus 39 decompression rows. | Executable on the 38 compression rows; `similar`. |
| Mg0.61Fe0.39O HS/LS, Zhuravlev data | Zhuravlev et al. (2010), DOI 10.1007/s00269-009-0347-6 | Measured Mg0.61(2)Fe0.39(3)O; compression collapse at 72.2–76.6 GPa | Complete Table 5: 54 compression plus 17 decompression rows. | Executable on the 54 compression rows; `similar`. |
| Mg0.40Fe0.60O HS/LS | Lin et al. (2005), DOI 10.1038/nature03825 | Mg0.40Fe0.60O; XES crossover constraint | The paper gives an 84–102 GPa spectroscopic interval and approximately 1.6% volume drop near 95 GPa, but no P–V series. | `not_refittable`; transition constraints are not observations. |

Zhuravlev's decompression rows are preserved because they are authoritative
observations, but they are excluded from the reconstruction: the paper reports
metastable LS persistence and hysteresis to about 65 GPa. Mixing them into the
compression fit would erase a documented path dependence.

Fei's Table S2 is handled the same way: the `mw39a` and `mw39b` compression
rows enter the reconstruction, while the 39 `mw39c` decompression rows are
preserved but excluded. The source lattice constants and their one-sigma errors
are retained directly; conventional-cell volumes are calculated as a³. The
column called `error1` belongs to the NaCl lattice constant and is not treated
as pressure uncertainty.

The user-supplied Marquardt supporting file is genuine, but it accompanies the
separate *Science* article “Elastic shear anisotropy of ferropericlase in
Earth's lower mantle” (DOI 10.1126/science.1169365). It contains methods,
figures, and shear/density data—not the P–V observations refitted by
Solomatova. Those 29 observations are printed directly in Table 2 of the
correct EPSL paper, DOI 10.1016/j.epsl.2009.08.017, and are the dataset bundled
here.

Zhuravlev et al. also own a distinct, directly published high-spin EOS. Table 6
and the Figure 9 caption identify the displayed/preferred fit as BM2 with
V0=77.4(2) A^3 per conventional cell and K0=161(3) GPa, using the 43 compression
rows from 15.5 through 70.9 GPa. It is now stored as
`mg061fe039o_zhuravlev_2010_high_spin_bm2_5`. An independent errors-in-variables
refit gives V0=77.2511 A^3 and K0=162.461 GPa, within the combined two-sigma
uncertainties of both published coefficients. The fixed-derivative Vinet and
free-derivative BM3/Vinet rows in Table 6 remain documented as sensitivity fits,
not additional production records.

## Coupled reconstruction

The reproduction follows the public MINUTI description rather than fitting the
HS and LS curves separately. For each trial parameter set it evaluates BM3
elastic energies for both reference states. Fe2+ populations use the published
three-level degeneracies HS:IS:LS = 15:18:1, state indices 0:1:2, and unpaired
fractions 1:1/2:0. The elastic energy difference is normalized by the four
cation sites in the conventional B1 cell and the Fe fraction. Mixed pressure is

`Pmix = PHS + (<alpha>/2) (PLS - PHS)`.

As in Table 7, K0′,LS is fixed at 4. Table 2 supplies priors K0,HS=160±5 GPa,
K0,LS=170±20 GPa, and K0′,HS=4±0.5 for the recovered compositions. The original
MINUTI input files, row weights, residual definition, and covariance were not
published. The reproducible audit therefore minimizes equal-weight approximate
volume residuals plus those priors. It tests whether the observations support a
similar coupled solution; it is not represented as the authors' exact
least-squares calculation.

| Source data | n | Reconstructed HS (V0, K0, K0′) | Reconstructed LS (V0, K0, K0′) | P50 / width20–80 (GPa) | Volume RMSE (Å³) | Ledger |
|---|---:|---|---|---:|---:|---|
| Marquardt Mg0.90Fe0.10O | 29 | 75.343, 161.471, 4.063 | 73.667, 170.250, 4 fixed | 53.821 / 10.663 | 0.1539 | `similar` |
| Mao Mg0.75Fe0.25O | 42 | 76.34 fixed, 160.954, 4.156 | 73.884, 171.981, 4 fixed | 61.515 / 16.091 | 0.0705 | `parity_not_achieved` |
| Chen nominal Mg0.65Fe0.35O | 36 | 77.264, 159.946, 3.993 | 73.841, 170.141, 4 fixed | 64.176 / 14.047 | 0.2059 | `similar` |
| Fei Mg0.61Fe0.39O, compression | 38 | 77.492, 160.256, 4.063 | 73.907, 172.384, 4 fixed | 56.561 / 14.486 | 0.1368 | `similar` |
| Zhuravlev Mg0.61Fe0.39O, compression | 54 | 77.325, 160.404, 4.139 | 73.531, 171.317, 4 fixed | 73.186 / 10.108 | 0.1975 | `similar` |
| Solomatova Fp48 | 45 | 77.29 fixed, 160.696, 4.098 | 73.673, 170.824, 4 fixed | 69.103 / 17.474 | 0.1198 | dedicated coupled diagnostic: `similar`; existing generic branch ledger unchanged |

The Lin Mg0.83Fe0.17O diagnostic fixes V0,HS to the Table 7 value 75.94 Å³.
With that external normalization it fits the 43 ratios with 0.2238 Å³ volume
RMSE and obtains K0,HS=158.755 GPa, K0′,HS=3.934, V0,LS=72.869 Å³,
K0,LS=177.720 GPa, P50=49.859 GPa, and width20–80=10.419 GPa. Those numbers
show that the recovered shape is compatible with a crossover solution; because
the defining absolute-volume anchor comes from the answer being tested, they
are deliberately absent from the successful-refit count.

The implementation and all machine-readable diagnostics are in
`scripts/reproduce_solomatova_2016_ferropericlase.py` and the generated common
ledger `docs/data/primary-eos-refits.json`. The five newly recovered datasets
are checksum-registered in their material documents.

## Observation-integrity rule

No pressure calculated from a Table 7 branch, analytical BM3 checkpoint,
transition midpoint/width, volume-collapse summary, Brillouin prior, or
coefficient repeated by a later paper is treated as an experimental
observation. Such values may constrain or diagnose a reconstruction, but only
source P–V rows or source figure markers contribute observation residuals.

## LitCurate disposition

All 16 same-DOI rows are accepted as the paired Table 7 reference branches:

| Source row | Candidate | Disposition | Branch |
|---:|---|---|---|
| 944 | `litcurate_d84461e22abf5080` | ACCEPT | Fp48 HS |
| 945 | `litcurate_1524ace141e7f7f0` | ACCEPT | Fp48 LS |
| 946 | `litcurate_99bd28480585bdaa` | ACCEPT | Mg0.90Fe0.10O HS |
| 947 | `litcurate_b0f8a936ede6bc5e` | ACCEPT | Mg0.90Fe0.10O LS |
| 948 | `litcurate_ae5f29f994c6d8f3` | ACCEPT | Mg0.83Fe0.17O HS |
| 949 | `litcurate_c7874f561714da3d` | ACCEPT | Mg0.83Fe0.17O LS |
| 950 | `litcurate_eb3e9425325e0c18` | ACCEPT | Mg0.75Fe0.25O HS |
| 951 | `litcurate_7f3cefd535845c56` | ACCEPT | Mg0.75Fe0.25O LS |
| 952 | `litcurate_28e6643f0cf0c399` | ACCEPT | Mg0.65Fe0.35O HS |
| 953 | `litcurate_85a318658060d6c0` | ACCEPT | Mg0.65Fe0.35O LS |
| 954 | `litcurate_67984b8ee2408ca4` | ACCEPT | Mg0.61Fe0.39O HS, Fei data |
| 955 | `litcurate_2fbe372e10a362f1` | ACCEPT | Mg0.61Fe0.39O LS, Fei data |
| 956 | `litcurate_4fe52dbe32e1dfdb` | ACCEPT | Mg0.61Fe0.39O HS, Zhuravlev data |
| 957 | `litcurate_72b31a925a02dab0` | ACCEPT | Mg0.61Fe0.39O LS, Zhuravlev data |
| 958 | `litcurate_1b209d91f037a754` | ACCEPT | Mg0.40Fe0.60O HS |
| 959 | `litcurate_112ff70f49853ccc` | ACCEPT | Mg0.40Fe0.60O LS |

Result: **16 accepted production records**. Acceptance means Table 7 really
publishes that reference branch; it does not imply that the earlier observations
needed to reproduce its coupled fit have been recovered for every pair except
the two Lin compositions described above.
