# LitCurate tranche C: zero-yield paper audits

This document records every candidate from the thirteen reserved DOI families that produced no executable EOS record. The discovery ledger is not treated as scientific authority. Source-reported rows remain on hold when the primary article does not supply a complete parameterization of an equation implemented by Peritheos; citation rows are rejected from this paper's production scope. No values are completed by borrowing parameters from a different paper or by treating local elastic properties as EOS coefficients.

## 10.1029/JB095iB05p07055 — Isaak et al. (1990), MgO

Primary target: “Calculated elastic and thermal properties of MgO at high pressures and temperatures,” *Journal of Geophysical Research* **95** (1990), DOI [10.1029/JB095iB05p07055](https://doi.org/10.1029/JB095iB05p07055). The publisher record and OpenAlex full-text locations were checked; no open primary article or repository copy was available. More importantly, all fourth-order rows in the discovery extraction omit `K0''`, and its 300 K third-order sensitivity row omits `V0`. They are not executable without inventing coefficients.

| Candidate | Disposition |
|---|---|
| `litcurate_7088e7412608f91a` (8) | hold — source BM4, missing `K0''`. |
| `litcurate_7a436978328704c5` (9) | hold — source BM4, missing `K0''`. |
| `litcurate_0729412c9d11ea76` (10) | hold — source BM4, missing `K0''`. |
| `litcurate_75b8d43422d6da7d` (11) | hold — source BM3 sensitivity row, missing `V0`. |
| `litcurate_a04115f4ab5c2c24` (12) | reject — citation comparison, not source-generated. |
| `litcurate_57e02f0c6ba5c712` (13) | hold — source BM4, missing `K0''`. |
| `litcurate_a1ae1385fbe72b0c` (14) | reject — citation comparison. |
| `litcurate_27883a983f187f54` (15) | hold — source BM4, missing `K0''`. |

Outcome: 0 accepted; 6 incomplete source rows held; 2 citation rows rejected.

## 10.1029/94GL01370 — Matsui, Price, and Patel (1994), MgSiO3 bridgmanite

Primary article: “Comparison between the lattice dynamics and molecular dynamics methods: Calculation results for MgSiO3 perovskite,” *Geophysical Research Letters* **21** (1994), 1659–1662, DOI [10.1029/94GL01370](https://doi.org/10.1029/94GL01370). The open primary PDF was checked directly (SHA-256 `c1d31ad31deffca5b8c472551ca96683f1160ef0b8eeb38adac52eaeae9a1a3f`), especially Table 2 and the method text. Table 2 gives isolated volumes and isothermal moduli for lattice-dynamics and molecular-dynamics calculations at 500 and 1000 K. Although the text says pressure-volume results were fitted with BM3, it reports no `K0'` and no pressure-volume grid. A local `(V,K_T)` pair is not a complete BM3 parameterization.

| Candidate | Disposition |
|---|---|
| `litcurate_c8c2d618f732e173` (92) | hold — 500 K lattice-dynamics local properties; missing `K0'` and fit data. |
| `litcurate_b084b549551e6079` (93) | hold — 500 K molecular-dynamics local properties; missing `K0'` and fit data. |
| `litcurate_18ab289ba61f3b16` (94) | hold — 1000 K lattice-dynamics local properties; missing `K0'` and fit data. |
| `litcurate_43ea05001719e9b0` (95) | hold — 1000 K molecular-dynamics local properties; missing `K0'` and fit data. |

Outcome: 0 accepted; 4 incomplete source rows held.

## 10.2465/jmps.97.13 — Matsui (2002), MgSiO3-Al2O3 bridgmanite

Primary article: “Molecular dynamics simulation of MgSiO3-Al2O3 perovskite,” *Journal of Mineralogical and Petrological Sciences* **97** (2002), DOI [10.2465/jmps.97.13](https://doi.org/10.2465/jmps.97.13). The primary J-STAGE PDF was checked directly (SHA-256 `58a7b8484c467a06c00e6053fdee48db519c4daaf3712ff033801affa418609e`). It states that the results were fitted to a fourth-order Birch-Murnaghan equation but reports only bulk moduli `182(2)`, `179(2)`, and `247.5(1)` GPa for these extracted branches. Reference volumes, `K0'`, `K0''`, and fit observations are absent.

| Candidate | Disposition |
|---|---|
| `litcurate_96d9eb1652572560` (418) | hold — Al-bearing branch, incomplete BM4. |
| `litcurate_1215c131c466928d` (419) | hold — Al-bearing branch, incomplete BM4. |
| `litcurate_562bdfda7b7eecd1` (420) | hold — MgSiO3 branch, incomplete BM4. |

Outcome: 0 accepted; 3 incomplete source rows held.

## 10.3390/min11030322 — Sokolova et al. (2021), Ca-silicates

Open primary article: “Equations of State of Ca-Silicates and Phase Diagram of the CaSiO3 System under Upper Mantle Conditions,” *Minerals* **11** (2021), 322, DOI [10.3390/min11030322](https://doi.org/10.3390/min11030322). The publisher article, its equations, and full parameter table were checked. The paper's EOS is a coupled Helmholtz model: Kunc cold compression (`k=5`), two Einstein oscillators, an Altshuler volume law for the Grüneisen parameter, and an optional anharmonic term. Peritheos does not currently implement that exact combined model. Its existing Holzapfel, multi-oscillator, and thermal wrappers have materially different cold or gamma laws, so translating the table would change the published equation.

| Candidate | Disposition |
|---|---|
| `litcurate_a361f7fdf95fcfaf` (1127), CaSiO3 davemaoite | hold — complete primary parameter table, but exact coupled Kunc/two-Einstein/Altshuler model is unsupported. |

The same primary table was also checked for wollastonite, pseudowollastonite, breyite, beta-larnite, and titanite-CaSi2O5. They are potentially valuable future records but are not separate LitCurate candidates in this ledger, and all require the same model implementation. Outcome: 0 accepted; 1 ledger row plus 5 non-ledger phase parameterizations held for model work.

## 10.1080/08957959.2018.1465056 — Sokolova et al. (2018), MgO-MgSiO3 spreadsheets

Primary citation: T. S. Sokolova, P. I. Dorogokupets, K. D. Litasov, B. S. Danilov, and A. M. Dymshits, “Spreadsheets to calculate P-V-T relations, thermodynamic and thermoelastic properties of silicates in the MgSiO3-MgO system,” *High Pressure Research* **38** (2018), 193–211, DOI [10.1080/08957959.2018.1465056](https://doi.org/10.1080/08957959.2018.1465056). The indexed primary text and equation description were checked. These are not standalone Birch-Murnaghan triplets as the discovery label suggests; they belong to the same Kunc cold compression plus two-Einstein/Altshuler Helmholtz formulation. The complete executable spreadsheet model cannot be mapped exactly to a supported Peritheos model.

| Candidate | Disposition |
|---|---|
| `litcurate_4f0815022700e115` (1005), MgO | hold — exact thermal model unsupported; ledger units are also normalized inconsistently (`K0=1630` corresponds to kbar). |
| `litcurate_3f13a798d5644822` (1006), MgSiO3 bridgmanite | hold — exact thermal model unsupported. |
| `litcurate_6f50fe3425cdda1e` (1007), MgSiO3 post-perovskite | hold — exact thermal model unsupported. |

Outcome: 0 accepted; 3 source parameterizations held for model work.

## 10.1029/97JB03672 — Hama and Suito (1998), Mg-Fe bridgmanite

Primary target: J. Hama and K. Suito, “Equation of state of MgSiO3 perovskite and its thermoelastic properties under lower mantle conditions,” *Journal of Geophysical Research* **103** (1998), 7443–7462, DOI [10.1029/97JB03672](https://doi.org/10.1029/97JB03672). The version of record is closed and no authoritative open full text was located. The source rows describe adopted or alternate theoretical calibrations within a larger thermoelastic construction, not four independently reproducible source datasets. The pure-Mg static coefficients are explicitly traced to Stixrude and Cohen rather than newly fitted here; the Fe-bearing rows lack a primary pressure-volume grid and complete thermal formulation usable in Peritheos.

| Candidate | Disposition |
|---|---|
| `litcurate_7c2e90ce19c25bd5` (178) | reject — citation-reported Vinet comparison. |
| `litcurate_f4f0dd4ab0003536` (179) | reject — citation comparison. |
| `litcurate_c55a1c21c98c355f` (180) | reject — citation comparison. |
| `litcurate_4ed1bea08316b761` (181) | hold — adopted pure-Mg calibration, not a new source fit. |
| `litcurate_c058fa990538a5f7` (182) | hold — alternate calibration/sensitivity row without reproducible primary fit data. |
| `litcurate_1695b56c291d4802` (183) | hold — Fe10 thermoelastic branch lacks reproducible primary fit data/model. |
| `litcurate_83d2833f57fd3634` (184) | hold — alternate Fe10 volume calibration, not an independent fit. |
| `litcurate_2f8531f9536ae515` (185) | reject — cited MgO comparison. |

Outcome: 0 accepted; 4 source rows held; 4 citation rows rejected.

## 10.1103/PhysRevLett.70.3947 — Wentzcovitch et al. (1993), MgSiO3

Primary article: “Ab initio molecular dynamics with variable cell shape: Application to MgSiO3,” *Physical Review Letters* **70** (1993), DOI [10.1103/PhysRevLett.70.3947](https://doi.org/10.1103/PhysRevLett.70.3947). The source reports modulus/derivative pairs for Pbnm and hypothetical cubic structures, but the candidate extraction contains no equilibrium volume and no source numerical pressure-volume table from which it can be recovered without digitizing a mixed comparison plot. A complete EOS is therefore not established.

| Candidate | Disposition |
|---|---|
| `litcurate_05b580c439188d20` (72), Pbnm MgSiO3 | hold — source BM3 missing `V0`. |
| `litcurate_899b5647e58e28a8` (73), hypothetical cubic MgSiO3 | hold — source BM3 missing `V0`; also requires a distinct phase card. |
| `litcurate_99ebd5d905f14ffb` (74) | reject — citation row. |
| `litcurate_ad313770f10a2f4e` (75) | reject — citation row. |
| `litcurate_7604cfa4ca090a22` (76) | reject — citation row. |

Outcome: 0 accepted; 2 incomplete source rows held; 3 citation rows rejected.

## 10.2138/am-2000-2-308 — Matsui, Parker, and Leslie (2000), MgO

Primary target: “The MD simulation of the equation of state of MgO: Application as a pressure calibration standard at high temperature and high pressure,” *American Mineralogist* **85** (2000), 312–316, DOI [10.2138/am-2000-2-308](https://doi.org/10.2138/am-2000-2-308). The publisher abstract/metadata and American Mineralogist issue page were checked, but the primary full text was access-controlled and the public MSA host was unavailable during the audit. The complete-looking BM3 row is therefore held under the primary-authority rule: secondary papers corroborate its use as a pressure scale, but do not substitute for auditing the source equation, reference state, fit basis, and data. The second source row is an elastic-modulus pseudo-split and is incomplete.

| Candidate | Disposition |
|---|---|
| `litcurate_6654fcf45fedcf5f` (283) | hold — complete-looking BM3 triplet, but primary equation/data could not be inspected. |
| `litcurate_3150c4d379bb93b4` (284) | reject pseudo-split — Table 2 elastic-modulus summary, labeled BM4 but missing `V0` and `K0''`. |
| `litcurate_ac66105cbbdf9a84` (285) | reject — Fei citation row. |
| `litcurate_64a2625d3b81147f` (286) | reject — Duffy et al. citation row. |
| `litcurate_e3e1a8a142c3341d` (287) | reject — cited observed elasticity, explicitly not an EOS fit. |

Outcome: 0 accepted; 1 source row held; 1 pseudo-split and 3 citation/non-EOS rows rejected.

## 10.1029/JB095iB13p21671 — Jackson et al. (1990), wuestite

Primary target: “Elasticity, shear-mode softening and high-pressure polymorphism of wüstite (Fe1-xO),” *Journal of Geophysical Research* **95** (1990), DOI [10.1029/JB095iB13p21671](https://doi.org/10.1029/JB095iB13p21671). The three source rows are acoustic/elastic extrapolations for Fe0.943O and omit `V0`; they are not three independently measured compression EOS fits. The remaining five rows quote earlier literature.

| Candidate | Disposition |
|---|---|
| `litcurate_353cbfa64c47a364` (17) | hold — source acoustic BM3-like parameter pair, missing `V0`. |
| `litcurate_4475d71e5f3e8602` (18) | hold — source acoustic sensitivity/extrapolation, missing `V0`. |
| `litcurate_f9fbfe2edd6c1c23` (19) | hold — source acoustic sensitivity/extrapolation, missing `V0`. |
| `litcurate_83aa355a5e33e72b` (20) | reject — citation row. |
| `litcurate_920593beed0790de` (21) | reject — citation row. |
| `litcurate_377006340598e87c` (22) | reject — citation row. |
| `litcurate_496f69bafaadd208` (23) | reject — citation row. |
| `litcurate_07eae2a616424c20` (24) | reject — citation row. |

Outcome: 0 accepted; 3 incomplete source rows held; 5 citation rows rejected.

## 10.1002/2016JB013136 — Pamato et al. (2016), NAL phase

Primary article: “Elasticity of single-crystal NAL phase at high pressure: A potential source of the seismic anisotropy in the lower mantle,” *Journal of Geophysical Research: Solid Earth* **121** (2016), DOI [10.1002/2016JB013136](https://doi.org/10.1002/2016JB013136). Both extracted entries are elastic/acoustic fourth-order finite-strain summaries. Each supplies a reference volume and bulk modulus but omits the pressure derivatives required for BM4; neither can be converted into an executable EOS without inventing coefficients.

| Candidate | Disposition |
|---|---|
| `litcurate_01173f3bae112c77` (889), Fe-free NAL | hold — incomplete BM4 (`K0'` and `K0''` absent). |
| `litcurate_dda72e8732d804bd` (890), Fe-bearing NAL | hold — incomplete BM4 (`K0'` and `K0''` absent). |

Outcome: 0 accepted; 2 incomplete source rows held.

## 10.1029/2023GL107700 — Chen et al. (2024), stishovite sound velocities

Primary article: “Sound Velocities of Stishovite at Simultaneous High Pressure and High Temperature Suggest an Eclogite-Rich Layer Beneath the Hawaii Hotspot,” *Geophysical Research Letters* **51** (2024), DOI [10.1029/2023GL107700](https://doi.org/10.1029/2023GL107700). The paper's new results are sound velocities and acoustic modulus derivatives, not a new pressure-volume fit. The two source rows omit `V0`; the other thirteen rows are literature comparisons. Creating separate EOS records would misclassify elastic parameter sets as compression fits.

| Candidate | Disposition |
|---|---|
| `litcurate_a233355f7e66fafa` (1227) | reject non-EOS — source acoustic parameterization, missing `V0`. |
| `litcurate_bafd34bf85a512d9` (1228) | reject non-EOS — source acoustic sensitivity branch, missing `V0`. |
| `litcurate_dde7f3905ad1e204` (1229) | reject citation row. |
| `litcurate_b211b20d7e6677b8` (1230) | reject citation row. |
| `litcurate_9fae46377e660a6f` (1231) | reject citation row. |
| `litcurate_03878893069c321d` (1232) | reject citation row. |
| `litcurate_b3d51ee59307c54a` (1233) | reject citation row. |
| `litcurate_6305d73204c339e7` (1234) | reject citation row. |
| `litcurate_dd21d3804100769c` (1235) | reject citation row. |
| `litcurate_d8817b7805164ae9` (1236) | reject citation row. |
| `litcurate_503cbaeda56e7b41` (1237) | reject citation row. |
| `litcurate_6c7ee222b246749d` (1238) | reject citation row. |
| `litcurate_944ce973a2fc2874` (1239) | reject citation row. |
| `litcurate_70a5790e03b68ce6` (1240) | reject citation row. |
| `litcurate_ceae2b8bb787b03a` (1241) | reject citation row. |

Outcome: 0 accepted; 2 source non-EOS rows and 13 citation rows rejected.

## 10.4131/jshpreview.21.272 — Ito et al. (2011), Kawai apparatus

Primary article: “High Pressure Generation in a Kawai-Type Apparatus Equipped with Sintered Diamond Anvils,” *Journal of High Pressure Institute of Japan* **49** (2011), DOI [10.4131/jshpreview.21.272](https://doi.org/10.4131/jshpreview.21.272). This apparatus paper quotes two alternative bulk-modulus/derivative pairs for (Mg0.92Fe0.08)SiO3 bridgmanite but does not provide the zero-pressure volume needed by BM3. Treating the two pressure-scale reductions as independent material EOS fits would be a pseudo run split.

| Candidate | Disposition |
|---|---|
| `litcurate_3ec628a4c6956bbb` (756) | reject pseudo-split — source pressure-reduction branch, missing `V0`. |
| `litcurate_5cf1159337ae1fd3` (757) | reject pseudo-split — alternate pressure-reduction branch, missing `V0`. |
| `litcurate_3b7821d637ffbb33` (758) | reject — citation row. |
| `litcurate_5c0608da3c3eecdb` (759) | reject — citation row. |

Outcome: 0 accepted; 2 incomplete pseudo-splits and 2 citation rows rejected.

## 10.1142/S0217979209053047 — Zhang et al. (2009), Mg-Fe silicate phase stability

Primary target: “Pressure-related phase stability of MgSiO3 and (Mg0.75,Fe0.25)SiO3 at lower mantle condition,” *International Journal of Modern Physics B* (2009), DOI [10.1142/S0217979209053047](https://doi.org/10.1142/S0217979209053047). The World Scientific version of record was inaccessible during the audit. The six source rows are labeled fourth-order Birch-Murnaghan but omit `K0''`. The extraction also collapses three distinct calculated structures per composition into “bridgmanite,” so phase identity cannot be repaired from the ledger alone. Five remaining rows are cited comparisons.

| Candidate | Disposition |
|---|---|
| `litcurate_9256ea847204584d` (704) | hold — source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_b6368f71984ac74c` (705) | hold — source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_8a7530a59938b623` (706) | hold — source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_d4d9c902b300d5f1` (707) | hold — Fe25 source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_7de9072c33726687` (708) | hold — Fe25 source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_5df999158b70ca39` (709) | hold — Fe25 source BM4 missing `K0''`; phase assignment unresolved. |
| `litcurate_3ae0fd98c9eff104` (710) | reject — citation row. |
| `litcurate_65c6b9d43f342be4` (711) | reject — citation row. |
| `litcurate_2e8f702f3e06140e` (712) | reject — citation row. |
| `litcurate_a6665f069d6e43c4` (713) | reject — citation row. |
| `litcurate_405cccade894529c` (714) | reject — citation row. |

Outcome: 0 accepted; 6 incomplete/phase-ambiguous source rows held; 5 citation rows rejected.

## 10.2138/am.2008.2717 — Reichmann et al. (2008), MgO elasticity

Primary article: “Single-crystal elastic properties of MgO to 9 GPa,” *American Mineralogist* **93** (2008), DOI [10.2138/am.2008.2717](https://doi.org/10.2138/am.2008.2717); [publisher full-text landing page](https://www.degruyterbrill.com/document/doi/10.2138/am.2008.2717/html). The authors measured pressure-dependent single-crystal elastic moduli and report an adiabatic aggregate bulk modulus, `K_S0 = 161.1(3) GPa`, with pressure derivative `K_S' = 4.2(2)`. This is an acoustic/elastic finite-strain result, not a pressure-volume EOS fit: the paper does not report the reference volume needed to execute BM3. The source row must therefore remain a non-EOS rejection; the other nine rows quote earlier work.

| Candidate | Disposition |
|---|---|
| `litcurate_0d92925ea0fa3e22` (655) | reject non-EOS — source adiabatic elastic-modulus fit, no compression `V0`. |
| `litcurate_a2854facc6d8fc97` (656) | reject — citation row. |
| `litcurate_62ed31ca713127cf` (657) | reject — citation row. |
| `litcurate_dcbd695d937d3c65` (658) | reject — citation row. |
| `litcurate_7f69da78a40a162c` (659) | reject — citation row. |
| `litcurate_485fedd37506ea97` (660) | reject — citation row. |
| `litcurate_d14edc051817fab4` (661) | reject — citation row. |
| `litcurate_ba6f2c963b03a4ce` (662) | reject — citation row. |
| `litcurate_c616011883319a72` (663) | reject — citation row. |
| `litcurate_0f334e36340e51fd` (664) | reject — citation row. |

Outcome: 0 accepted; 1 source non-EOS row and 9 citation rows rejected.

## Batch outcome

These fourteen DOI families contribute **0** production records. Across their 87 LitCurate rows, 41 source rows are held or rejected as incomplete, unsupported-model, pseudo-split, or non-EOS entries, and 46 citation rows are rejected. The separately documented Karki et al. (1997) audit contributes two accepted records from the same reservation set after primary-source resolution of missing ledger fields.
