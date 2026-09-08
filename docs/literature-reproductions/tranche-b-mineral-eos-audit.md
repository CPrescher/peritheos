# Tranche B mineral EOS audit

Audit date: 2026-09-05. This audit applies a source-owned-fit rule: a production record requires an identifiable material, an executable EOS family, a complete coefficient set, and a fit performed in the cited primary paper. Literature-comparison rows are citation traces, not new records. Values transcribed from primary tables are preserved at source precision; no missing coefficients or compositions are inferred.

## Accepted production records

Twenty-four defensible records are retained from ten papers. The analytical curves are reproduced by `scripts/reproduce_tranche_b_mineral_eos.py` and tested by `tests/test_tranche_b_mineral_eos.py`.

| Primary source | Accepted records | Evidence |
|---|---:|---|
| Sherman (1993), doi:10.1029/93JB02175 | 1 | Abstract/source fit discussion: complete basis-set-B static BM3. |
| Wang and Weidner (1994), doi:10.1029/94GL00976 | 1 | Primary abstract and compression discussion: one room-temperature BM2 fit above 2 GPa. |
| Akber-Knutson, Bukowinski, and Matas (2002), doi:10.1029/2001GL013523 | 4 | Primary Tables 1–2: one stishovite and three CaSiO3 source-owned BM3 fits. |
| Tsuchiya and Mookherjee (2015), doi:10.1038/srep15534 | 1 | Primary Results, Figure 2, Table 1, and Methods: ordered model-1 phase-H BM3. |
| Sagatova et al. (2021), doi:10.1134/S0016702921080073 | 7 | Official full text, Methods and Table 2: seven source-owned room-temperature Vinet parameterizations. |
| Liu et al. (2011), doi:10.1088/1674-0068/24/06/703-710 | 1 | Primary Sections II–III and Table I: source-owned static PBE-GGA BM3. |
| Ismailova et al. (2016), doi:10.1126/sciadv.1600427 | 1 | Primary article and supplement: one measured Fe-deficient bridgmanite BM2. |
| Ricolleau et al. (2009), doi:10.1029/2008GL036759 | 4 | Official full text, Section 3 and Table 1: three independent thermal BM2 fits and one room-temperature low-spin branch. |
| Liu et al. (2007), doi:10.1088/0953-8984/19/24/246103 | 1 | Primary calculation section and Table 1: source-owned cubic CaSiO3 LDA BM3. |
| Bykova et al. (2018), doi:10.1038/s41467-018-07265-z | 3 | Official open article and supplement: one published coesite-I/II/III BM3 plus two explicitly derived AM05 BM3 reconstructions from source Tables 10–11. |

The phase-H dataset contains all 12 rows of primary Table 1. Its lattice-derived volumes reproduce the published BM3 with a pressure RMS residual of 0.5462582 GPa and maximum absolute residual of 1.1329591 GPa; the residual is consistent with three-decimal lattice-parameter rounding. The Sagatova Ca-perovskite values are printed per formula unit and multiplied by four for the conventional I4/mcm material-card basis; all other Sagatova Table 2 volumes already use the respective conventional-cell basis. Ricolleau's three thermal records use the stated integrated linear-in-temperature expansivity and linear bulk-modulus temperature derivative, with a 300 K BM2 reference curve. For Bykova's two theoretical high-pressure phases, Supplementary Table 11 prints K0 and K0′ but omits V0. The derived records hold those coefficients fixed and reconstruct only V0 from exact Supplementary Table 10 calculated P–V states. Coesite-IV uses three unambiguously phase-retaining states and reproduces them to 0.080812 GPa RMS; coesite-V is exactly anchored by its tabulated 57 GPa state and checked against the plotted curve. They are marked `derived`, not misrepresented as fully printed parameter triples.

## Row-by-row disposition

### Sherman (1993), 10.1029/93JB02175

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_7d3b7069671bd006 | Source | Accept | Complete basis-set-B BM3 → `ca_perovskite_sherman_1993_basis_b_static_bm3`. |
| litcurate_b27524b22bd0085c | Source | Hold | The row says “finite strain” but the exact analytical form is unresolved; it cannot be mapped safely from the coefficient triple alone. |
| litcurate_63c16f922d10a453 | Source | Reject | Small-basis calculation omits K0′ and is incomplete for BM3. |
| litcurate_55fe40018a8643d7 | Citation | Reject | Mao et al. comparison value, incomplete and not fit in this paper. |

### Wang and Weidner (1994), 10.1029/94GL00976

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_0de2ffa4e2cfb6f6 | Source | Accept | Independent room-temperature BM2 → `ca_perovskite_wang_weidner_1994_bm2`. BM2 encodes K0′=4 structurally, so K0′ is documented but not passed as a model parameter. |
| litcurate_d683536d7eee5d30 | Citation | Reject | Comparison modulus; V0 and EOS family absent. |
| litcurate_d3f946dd383b2e3c | Citation | Reject | Comparison row, not fit in the source. |
| litcurate_c277deef8ed8e045 | Citation | Reject | Comparison modulus; V0 and EOS family absent. |
| litcurate_894404ba8787ed75 | Citation | Reject | Mg-Fe bridgmanite comparison, wrong material and not source-owned. |
| litcurate_0f857c12911d9183 | Source | Reject duplicate | Molar-volume restatement/rounding of the accepted CaSiO3 fit, not an independent fit. |

### Karki et al. (2001), 10.1029/2001GL012910

The primary Table 1 reports thermodynamic derivatives obtained from quasiharmonic free energies, but does not identify an executable volume-EOS analytical form. The complete-looking V0/K0/K0′ triples therefore remain on hold rather than being guessed as BM3 or Vinet.

| LitCurate candidate | Origin | Disposition | Reason |
|---|---|---|---|
| litcurate_02bcac74964d7a33 | Source | Hold | 300 K derivative triple; EOS family not specified. |
| litcurate_cfe89bed76d0e0f3 | Source | Hold | Static derivative triple; EOS family not specified. |
| litcurate_68023a6beb314d4f | Citation | Reject | External comparison row. |
| litcurate_481eff2d0b46c7c3 | Citation | Reject | External Fe-bearing comparison row. |
| litcurate_e621755283ff16bd | Citation | Reject | External comparison and incomplete V0. |
| litcurate_2ad76a15fa5c8c46 | Citation | Reject | External Al-bearing comparison and incomplete coefficients. |

### Akber-Knutson et al. (2002), 10.1029/2001GL013523

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_591dba9831c14ec5 | Source | Accept | Thermally corrected VIBC stishovite BM3 → `stishovite_akber_knutson_2002_vibc_300k_bm3`. |
| litcurate_672e78b5a3a6d678 | Citation | Reject | Stishovite comparison and incomplete K0′. |
| litcurate_1353e38c544a4b77 | Citation | Reject | Stishovite comparison, not source-owned. |
| litcurate_b378bfb34ea70e3a | Source | Accept | VIBC Pnma CaSiO3 BM3 → `ca_perovskite_akber_knutson_2002_vibc_pnma_300k_bm3`. |
| litcurate_86be9b10a1cde79f | Source | Accept | VIBC Pm3m CaSiO3 BM3 → `ca_perovskite_akber_knutson_2002_vibc_pm3m_300k_bm3`. |
| litcurate_07e49dfd60cc66b9 | Source | Accept | VIB Pm3m CaSiO3 BM3 → `ca_perovskite_akber_knutson_2002_vib_pm3m_300k_bm3`. |
| litcurate_d9d3f0397cd488b4 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_3e3ccee45a3b0658 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_7177e71c5a4b194e | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_0226a4a31e93f474 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_e2af4883e4dfbdeb | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_2f1b5e23cf321c5b | Citation | Reject | Literature comparison, not source-owned; anomalous basis is not repaired. |
| litcurate_9f0a381357727abb | Citation | Reject | Literature comparison, not source-owned; anomalous basis is not repaired. |
| litcurate_c3a768d08977ad22 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_412c69bc20cee269 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_7a48df5a6e272006 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_051fd43130f4a8c9 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_ab0303b377a2c179 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_5d5b88e06aa9cbe7 | Citation | Reject | Literature comparison, not source-owned. |

### Faust and Knittle (1994), 10.1029/94GL01592

The source used a natural F-rich chondrodite (Berkeley specimen 13056). The reusable chemical formula and exact site occupancies of that specific specimen cannot be reconstructed confidently from the accessible primary evidence. The second row reuses the same mechanical coefficients in a compositional-density calculation and is not another EOS fit.

| LitCurate candidate | Origin | Disposition | Reason |
|---|---|---|---|
| litcurate_adb4b30081847c84 | Source | Hold | Real BM3 fit, but exact specimen composition/material identity remains unresolved. |
| litcurate_faae4e18134071da | Source | Reject pseudo-split | Assumed-composition density model reusing the same K0 and K0′, not an independent compression fit. |

### Tsuchiya and Mookherjee (2015), 10.1038/srep15534

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_1b09c1018d6788fa | Source | Accept | Freely fitted K0′ identifies BM3; ordered P2/m model-1 → `phase_h_tsuchiya_mookherjee_2015_gga_static_bm3`. |

### Sagatova et al. (2021), 10.1134/S0016702921080073

The English translation DOI is used as the canonical identifier. The official Russian full text also displays 10.31857/S0016752521080070. Table 2 contains seven complete source-owned Vinet rows. LitCurate exposes only the two Ca-perovskite rows; the five additional phase rows are accepted as a same-paper primary-table expansion: `wollastonite_sagatova_2021_gga_300k_vinet`, `pseudowollastonite_sagatova_2021_gga_300k_vinet`, `breyite_sagatova_2021_gga_300k_vinet`, `casio2o5_titanite_sagatova_2021_gga_300k_vinet`, and `larnite_sagatova_2021_gga_300k_vinet`.

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_e35512b57defe823 | Source | Accept | LDA tetragonal CaSiO3 Vinet → `ca_perovskite_tetragonal_sagatova_2021_lda_300k_vinet`. |
| litcurate_00a482d0d9c2f37d | Source | Accept | GGA tetragonal CaSiO3 Vinet → `ca_perovskite_tetragonal_sagatova_2021_gga_300k_vinet`. |
| litcurate_2626bbd8e5f3794a | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_2cce7733b9056366 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_0e811e23b63c35c2 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_501e43b82eaac5ba | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_a7b2acb62404c667 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_9f43e7114494d3cb | Citation | Reject | Literature comparison, not source-owned; anomalous basis is not repaired. |
| litcurate_c1fc812a85823b68 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_d9bef591c80666a8 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_30312d58adb66eb3 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_7dae2aa15793f3d9 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_bd59ae56d23dcd7c | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_1f72ca5438d8da85 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_1b0f0ed4f1cab43d | Citation | Reject | Literature comparison, not source-owned. |

### Schoelmerich et al. (2020), 10.1038/s41598-020-66340-y

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_399b5f6c219e8346 | Source | Reject | The reported BM3 coefficients cannot be independently reproduced: corrected 300 K states are graphical and the complete correction and weighting protocol is absent. Primary Table 1 shock states remain as an audit fixture. |
| litcurate_6afe201aeae497fc | Citation | Reject | Prior shock-compression comparison, not fit in this source and lacks V0. |
| litcurate_71b55cdf996ff53b | Citation | Reject | Prior DAC comparison, not source-owned. |
| litcurate_8fae504a33ad0c78 | Citation | Reject | Literature comparison, not source-owned. |
| litcurate_fdbd04aa602db711 | Citation | Reject | Literature comparison, not source-owned. |

### Liu et al. (2011), 10.1088/1674-0068/24/06/703-710

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_93f5bf006271340f | Source | Accept | Static PBE-GGA orthorhombic MgSiO3 BM3 → `bridgmanite_liu_2011_gga_static_bm3`. |
| litcurate_cdb244bac79e77a9 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_e16b865a5f3f0cee | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_8cee3e29abed2229 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_adf74af6252e7375 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_22995f1a5b7229a9 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_a5650795ab13d13d | Citation | Reject | Table-I literature comparison, not source-owned. |

### Ismailova et al. (2016), 10.1126/sciadv.1600427

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_607d82ad74297299 | Source | Accept | Independent BM2 for measured `(Fe2+0.64(2)Fe3+0.24(2))Si1.00(3)O3` with about 12% A-site vacancies → `fe088sio3_bridgmanite_ismailova_2016_300k_bm2`. |

### Ricolleau et al. (2009), 10.1029/2008GL036759

The records describe coexisting phases separated from natural KLB-1 pyrolite. Their compositions evolve slightly with pressure and temperature, so the material cards retain sample-specific symbolic formulas instead of inventing fixed endmember stoichiometries.

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_8bcba2433560341a | Source | Accept | KLB-1 Mg-perovskite thermal BM2 → `klb1_mg_perovskite_ricolleau_2009_bm2_alphakt`. |
| litcurate_3cb77d9522d73b1a | Source | Accept | KLB-1 Ca-perovskite thermal BM2 → `klb1_ca_perovskite_ricolleau_2009_bm2_alphakt`. |
| litcurate_2a71d6f77a5e3c68 | Source | Accept | High-spin KLB-1 ferropericlase thermal BM2 → `klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt`. |
| litcurate_3d74852af1005f2f | Source | Accept | Distinct 300 K low-spin ferropericlase BM2 branch → `klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2`. |

### Liu et al. (2007), 10.1088/0953-8984/19/24/246103

| LitCurate candidate | Origin | Disposition | Reason / production mapping |
|---|---|---|---|
| litcurate_4314c5f756c351cf | Source | Accept | Static LDA cubic CaSiO3 BM3 → `ca_perovskite_liu_2007_lda_static_bm3`. |
| litcurate_98b4ea6477e4ecb0 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_4dabf08d83f12857 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_f440dcb27777d70b | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_0f0c3363fca0f155 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_202ccf5f7c369b2a | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_605467c894143180 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_a9d66e249068b2a3 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_32e7ac9bd60655f3 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_b3e298bc8a7859e7 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_5aa2d6e3abc959c5 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_8fd9b28b9b180f0a | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_746dc0bcfbf16267 | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_66272925694c4b0c | Citation | Reject | Table-I literature comparison, not source-owned. |
| litcurate_10c181958cdb6ec5 | Citation | Reject | Table-I literature comparison, not source-owned. |

### Bykova et al. (2018), 10.1038/s41467-018-07265-z

LitCurate contains no candidate rows for this DOI. It was investigated as an explicitly documented same-topic primary-source expansion. Three independent and executable curves are retained: the published combined experimental coesite-I/II/III BM3, and separately calculated coesite-IV and coesite-V AM05 BM3 curves. The latter are transparent source-table derivations because the paper prints K0 and K0′ but not V0. The experimental combined coesite-IV/V BM2 is not added: it is parameterized at a 32.7 GPa reference state (`V32.7`, `K32.7`), which the current zero-pressure BM2 implementation cannot represent exactly without an unsupported pressure offset or a new refit.

### Shieh, Duffy, and Li (2002), 10.1103/PhysRevLett.89.255507

The complete APS primary PDF was checked. The source fits the hydrostatic-orientation points below 50 GPa with a Birch–Murnaghan equation, fixes K0=304 GPa, and reports fitted K0′=4.6(1), but never gives the V0 used in that fit. An executable EOS cannot be reconstructed without importing an ambient volume from another source, so the single genuine fit remains on hold.

| LitCurate candidate | Origin | Disposition | Reason |
|---|---|---|---|
| litcurate_c47b28c0e27876d7 | Source | Hold | Real source fit, but V0 is absent; the candidate also missed the primary text's K0′=4.6(1). |
| litcurate_f0f71337af5d108b | Citation | Reject | Prior He-medium compression comparison, not source-owned. |
| litcurate_5f4c7da0e39bed86 | Citation | Reject | Prior ultrasonic comparison, not source-owned. |
| litcurate_d0516ae794933ab0 | Citation | Reject | Prior shock comparison, not source-owned. |
| litcurate_c5e67abcbf8d9c1b | Citation | Reject | Prior quenched laser-heating comparison, not source-owned. |

### Sinogeikin and Bass (1999), 10.1103/PhysRevB.59.R14141

The complete APS primary PDF was checked. Table I reports a source-owned fourth-order finite-strain bulk-modulus parameterization (`K0=163.2(10) GPa`, `K0′=4.0(1)`, `K0″=-0.04(2) GPa^-1`), but does not publish V0. The measured density is obtained iteratively inside the elasticity inversion. Without the reference volume, no executable BM4 record is defensible.

| LitCurate candidate | Origin | Disposition | Reason |
|---|---|---|---|
| litcurate_757f9b8c49305ab4 | Source | Hold | Genuine fourth-order finite-strain result, but V0 is not published; LitCurate also omits K0″. |
| litcurate_d5dc28e03b8b138b | Citation | Reject | Jackson and Niesler comparison values, not source-owned. |
| litcurate_a3b6c75bf2d11f1b | Citation | Reject | Yoneda comparison values, not source-owned. |

## Pressure calibration and limitations

The theoretical records use `not_applicable` pressure calibration with an ab-initio method. Wang and Weidner identify NaCl as their pressure standard, but the exact calibration equation and complete paired calibrant observations were not resolved, so recalculation is explicitly marked unavailable. Ismailova's exact pressure calibration was likewise left explicitly unresolved. Ricolleau used the Fei et al. (2007) Au scale; recalculation remains unavailable because the paired Table-S1 calibrant observations were not bundled. No covariance matrix is invented where the papers do not publish one. The Sagatova source plots relative-volume curves but does not publish the calculated P–V grids, so numerical checking is limited to analytical anchors and monotonic compressed-volume checkpoints.

## Zotero-ready primary-source metadata

- Sherman, David M. (1993). “Equation of state, elastic properties, and stability of CaSiO3 perovskite: First principles (periodic Hartree-Fock) results.” *Journal of Geophysical Research: Solid Earth* 98(B11), 19795–19805. doi:10.1029/93JB02175.
- Wang, Yanbin; Weidner, Donald J. (1994). “Thermoelasticity of CaSiO3 perovskite and implications for the lower mantle.” *Geophysical Research Letters* 21(10), 895–898. doi:10.1029/94GL00976.
- Akber-Knutson, Sofia; Bukowinski, Mark S. T.; Matas, Jan (2002). “On the structure and compressibility of CaSiO3 perovskite.” *Geophysical Research Letters* 29(3), 4-1–4-4. doi:10.1029/2001GL013523.
- Tsuchiya, Jun; Mookherjee, Mainak (2015). “Crystal structure, equation of state, and elasticity of phase H (MgSiO4H2) at Earth's lower mantle pressures.” *Scientific Reports* 5, 15534. doi:10.1038/srep15534.
- Sagatova, Dinara N.; Shatskiy, Anton F.; Sagatov, Nursultan E.; Litasov, Konstantin D. (2021). “Phase Relations in CaSiO3 System up to 100 GPa and 2500 K.” *Geochemistry International* 59(8), 791–800. doi:10.1134/S0016702921080073.
- Schoelmerich, M. O.; Tschentscher, T.; Bhat, S.; Bolme, C. A.; Cunningham, E.; Farla, R.; Galtier, E.; Gleason, A. E.; Harmand, M.; Inubushi, Y.; Katagiri, K.; Miyanishi, K.; Nagler, B.; Ozaki, N.; Preston, T. R.; Redmer, R.; Smith, R. F.; Tobase, T.; Togashi, T.; Tracy, S. J.; Umeda, Y.; Wollenweber, L.; Yabuuchi, T.; Zastrau, U.; Appel, K. (2020). “Evidence of shock-compressed stishovite above 300 GPa.” *Scientific Reports* 10, 10197. doi:10.1038/s41598-020-66340-y.
- Liu, Zi-jiang; Sun, Xiao-wei; Zhang, Cai-rong; Hu, Jian-bo; Song, Ting; Qi, Jian-hong (2011). “Elastic Tensor and Thermodynamic Property of Magnesium Silicate Perovskite from First-principles Calculations.” *Chinese Journal of Chemical Physics* 24(6), 703–710. doi:10.1088/1674-0068/24/06/703-710.
- Ismailova, Leyla; Bykova, Elena; Bykov, Maxim; Cerantola, Valerio; McCammon, Catherine; Boffa Ballaran, Tiziana; Bobrov, Andrei; Sinmyo, Ryosuke; Dubrovinskaia, Natalia; Glazyrin, Konstantin; Liermann, Hanns-Peter; Kupenko, Ilya; Hanfland, Michael; Prescher, Clemens; Prakapenka, Vitali; Svitlyk, Volodymyr; Dubrovinsky, Leonid (2016). “Stability of Fe,Al-bearing bridgmanite in the lower mantle and synthesis of pure Fe-bridgmanite.” *Science Advances* 2(7), e1600427. doi:10.1126/sciadv.1600427.
- Ricolleau, Angèle; Fei, Yingwei; Cottrell, Elizabeth; Watson, Heather; Deng, Liwei; Zhang, Li; Fiquet, Guillaume; Auzende, Anne-Line; Roskosz, Mathieu; Morard, Guillaume; Prakapenka, Vitali (2009). “Density profile of pyrolite under the lower mantle conditions.” *Geophysical Research Letters* 36(6), L06302. doi:10.1029/2008GL036759.
- Liu, Z. J.; Sun, X. W.; Chen, Q. F.; Cai, L. C.; Wu, H. Y.; Ge, S. H. (2007). “First-principles study of the elastic and thermodynamic properties of CaSiO3 perovskite.” *Journal of Physics: Condensed Matter* 19(24), 246103. doi:10.1088/0953-8984/19/24/246103.
- Bykova, E.; Bykov, M.; Černok, A.; Tidholm, J.; Simak, S. I.; Hellman, O.; Belov, M. P.; Abrikosov, I. A.; Liermann, H.-P.; Hanfland, M.; Prakapenka, V. B.; Prescher, C.; Dubrovinskaia, N.; Dubrovinsky, L. (2018). “Metastable silica high pressure polymorphs as structural proxies of deep Earth silicate melts.” *Nature Communications* 9, 4789. doi:10.1038/s41467-018-07265-z.

Investigated zero-yield primary sources (not import targets): Shieh, Sean R.; Duffy, Thomas S.; Li, Baosheng (2002), “Strength and Elasticity of SiO2 across the Stishovite–CaCl2-type Structural Phase Boundary,” *Physical Review Letters* 89, 255507, doi:10.1103/PhysRevLett.89.255507; Sinogeikin, S. V.; Bass, J. D. (1999), “Single-crystal elasticity of MgO at high pressure,” *Physical Review B* 59, R14141–R14144, doi:10.1103/PhysRevB.59.R14141.

The two held papers are also fully documented above, but are not presented as import targets for accepted production records.
