# Solomatova et al. (2016): spin-crossover ferropericlase branches

Primary source: N. V. Solomatova, J. M. Jackson, W. Sturhahn, J. K. Wicks,
J. Zhao, T. S. Toellner, B. Kalkan, and W. M. Steinhardt, “Equation of state and
spin crossover of (Mg,Fe)O at high pressure, with implications for explaining
topographic relief at the core-mantle boundary,” *American Mineralogist* **101**,
1084–1093 (2016), <https://doi.org/10.2138/am-2016-5510>.

The author-hosted primary PDF was checked against Tables 1–3 and 7 and the
MINUTI fitting discussion. Table 7 reports eight coupled spin-crossover fits,
each as high-spin and low-spin third-order Birch–Murnaghan reference branches;
`K0',LS=4` is fixed. The transition pressure and width belong to the mixer and
are preserved as provenance, but the current Peritheos EOS layer does not
implement MINUTI's population mixer. Consequently, the 16 executable records
are explicitly labelled **reference branches** and must not be mistaken for the
mixed-spin crossover curve.

These are genuine composition/source-data branches, including two independent
fits to Mg0.61Fe0.39O data from Fei et al. (2007) and Zhuravlev et al. (2010).
They are not pressure-window or run pseudo-splits. For the study's exact Fp48
sample `(Mg0.490Fe0.483Ti0.027)O`, all 45 Table 1 P-V observations are bundled.
The published high-spin branch reproduces the 24 measured pre-crossover rows
through 44.4 GPa with 0.52 GPa RMSE; the low-spin branch reproduces the post-crossover
Table 3 grid from 88–140 GPa with 0.30 GPa RMSE. The earlier datasets were not
reprinted numerically, so their source-parameter curves are checked by an
independent BM3 implementation without claiming a new refit.

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

Result: **16 accepted production records**. Each record's validity notes exclude
the mixed-spin interval and direct users to the paired branch provenance.
