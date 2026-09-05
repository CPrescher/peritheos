# Paper investigation ledger

Generated from the primary-source audit, the record-level refit ledger, and
the explicit nonproduction investigation register. This page answers a
different question from the EOS catalog: it records what happened to every
primary paper that was actually investigated, including papers that did not
produce an executable record.

## Status definitions

- **Reproduced:** every executable record from the paper reached `parity` or
  `similar` in the documented independent check.
- **Partly reproduced:** at least one record was reproduced, while another
  could not be refitted directly from available row-level evidence.
- **Coefficient parity not achieved:** the refit ran, but at least one
  published coefficient was outside both the uncertainty and numerical
  similarity criteria. These are source-fit discrepancies, not software-run
  failures; the record may remain for faithful published-curve provenance.
- **Direct refit unavailable:** the equation and parameters were audited, but
  independent coefficient recovery was impossible because primary rows, an
  executable calibration, or the original reduction were unavailable or
  circular.
- **Withheld/deferred:** investigation did not pass the executable-record
  acceptance gate, so no production EOS was added.

## Summary

The register covers **150 primary papers**: **148** support the 233 audited catalog records and **2** were investigated without adding a production record.
No numerical refit attempt failed before producing a comparison. The adverse
outcomes are instead explicit coefficient discrepancies, unavailable direct
refits, or acceptance-gate holds.

| Paper-level outcome | Papers |
|---|---:|
| Reproduced | 110 |
| Partly reproduced | 3 |
| Mixed: reproduced and discrepant records | 3 |
| Coefficient parity not achieved | 9 |
| Direct refit unavailable | 23 |
| Withheld: could not reproduce | 1 |
| Deferred: incomplete source/model mapping | 1 |

## Withheld or deferred papers

### [Katsura et al. (2004), Mg2SiO4 ringwoodite](https://doi.org/10.1029/2004JB003094)

**Outcome:** Withheld: could not reproduce (2026-09-05).

All 127 official Table 2 observations were transcribed, but the published BM3-MGD coefficients miss them by 1.890 GPa RMSE when the chemically required n=7 atoms per formula unit is used. A refit moves gamma0 far outside its reported uncertainty, while an unphysical and unpublished n=5 normalization fits much better. The hidden normalization, Debye-temperature law, and energy/volume convention remain unresolved, so no production EOS was added.

Evidence: [literature-reproductions.md#ringwoodite-katsura-2004](literature-reproductions.md#ringwoodite-katsura-2004).

### [Wang et al. (2026), KAlSi3O8 liebermannite and K-hollandite II](https://doi.org/10.2138/am-2024-9562)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-04).

The deposited temperature-indexed BM3 curves can be reproduced, but they are derived grids rather than primary fit observations. Their V0(T), K0(T), and K0-prime(T) behavior is not faithfully represented by the current thermal wrapper, and the primary methods, uncertainty/covariance information, exact reference-state convention, and authoritative symmetry-reduced structures were unavailable. The candidate remains deliberately non-executable.

Evidence: [material-eos-candidates.md](material-eos-candidates.md).

## Papers with coefficient discrepancies

These **12 papers** account for all 18 records
classified as `parity_not_achieved`. Papers with other successful records
are marked as mixed in the complete register.

| Paper | Affected record | Published-to-refit discrepancy |
|---|---|---|
| [Anzellini et al. (2019)](https://doi.org/10.1038/s41598-019-51931-1) | `silicon_vii_anzellini_2019_vinet_1` | K0 96.9 -> 4.845 |
| [Baty et al. (2024)](https://doi.org/10.1063/5.0179469) | `palladium_baty_2024_bm3_1` | K0 190 -> 152.057 |
| [Clendenen and Drickamer (1966)](https://doi.org/10.1063/1.1726610) | `coo_clendenen_1966_murnaghan_1` | K0_prime 3.9 -> 5.1481 |
| [Dobrosavljevic et al. (2019)](https://doi.org/10.3390/min9120762) | `mgfe94o_b1_dobrosavljevic_2019_bm3_1` | K0_prime 3.79 -> 2.71444 |
|  | `mgfe94o_rhombohedral_dobrosavljevic_2019_bm3_1` | K0 217 -> 168.771 |
| [Finkelstein et al. (2017)](https://doi.org/10.2138/am-2017-5966) | `mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1` | K0 148 -> 337.791; K0_prime 4.09 -> 9.6577 |
|  | `mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2` | K0 163 -> 360.366; K0_prime 4.02 -> 9.41478 |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | `goethite_gleason_2008_bm3_1` | rt_eos.K0 140.3 -> 183.338; rt_eos.K0_prime 4.6 -> 0 |
| [Jacobsen et al. (2002)](https://doi.org/10.1029/2001jb000490) | `ferropericlase_fe27_jacobsen_2002_bm3_1` | V0 76.336 -> 113.015; K0 158.4 -> 37.1717 |
|  | `magnesiowustite_fe56_jacobsen_2002_bm3_1` | V0 77.453 -> 113.002; K0 155.8 -> 37.1669 |
|  | `magnesiowustite_fe75_jacobsen_2002_bm3_1` | V0 78.082 -> 113.014; K0 151.3 -> 37.188 |
| [Jacobsen et al. (2005)](https://doi.org/10.1107/s0909049505022326) | `fe093o_b1_jacobsen_2005_bm3_1` | V0 79.41 -> 59.7395 |
|  | `mg073fe027o_jacobsen_2005_bm3_1` | V0 77.3 -> 57.963; K0_prime 4 -> 2.38146 |
| [Katsura et al. (2009)](https://doi.org/10.1029/2009gl038107) | `wadsleyite_katsura_2009_bm3_1` | gamma0 1.64 -> 1.12567 |
| [Ono et al. (2000)](https://doi.org/10.1007/s002690000108) | `sno2_cubic_27gpa_ono_2000_bm3_1` | K0 252 -> 379.59 |
|  | `sno2_pa_3_at_48gpa_ono_2000_bm3_1` | K0 252 -> 379.59 |
| [Somayazulu et al. (2023)](https://doi.org/10.1098/rsta.2022.0331) | `b4c_somayazulu_2023_bm3_1` | q 2.1 -> 1.04991 |
| [Thompson et al. (2017)](https://doi.org/10.1002/2017jb014168) | `e_feooh_hc_low_spin_thompson_2017_bm3_1` | K0 223 -> 185.921 |

The full data selection, model mapping, residuals, and bounded explanation
for every row above are in the
[primary EOS refit validation](primary-eos-refits.md#detailed-non-parity-investigations).

## Papers with unavailable direct refits

These **28 papers** contain 49 records for which a
source-faithful coefficient refit could not be performed. A paper can also
have other records that were reproduced.

| Paper | Affected records | Why direct refitting was unavailable |
|---|---|---|
| [Anderson et al. (1989)](https://doi.org/10.1063/1.342969) | `gold_anderson_1989_bm3_1` | The paper derives a thermodynamic gold parameterization from literature properties and publishes coefficient tables, not a new row-level experimental compression dataset. |
| [Anzellini et al. (2025)](https://doi.org/10.1038/s43246-025-00963-4) | `iridium_anzellini_2025_bm3_1` | The bundled rows are all heated states. The stored coefficients are the 300 K reference part of a combined thermal fit, but the record does not represent the source's thermal correction needed to refit those rows. |
| [B1 Fe0.94O, Fischer et al. (2011)](https://doi.org/10.1016/j.epsl.2011.02.025) | `feo_b8_2_fischer_2011_bm3_1`, `feo_fischer_2011_bm3_2` | Only 1 observation(s) lie at the reference temperature for 2 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |
| [Baty et al. (2024)](https://doi.org/10.1063/5.0179469) | `palladium_baty_2024_bm3_dft_2` | The calculated P(V) grid is not published as independent row-level fit input. Table S3 contains pressures generated from the already fitted EOS at selected volumes. |
| [Benedict et al. (2014)](https://doi.org/10.1103/physrevb.89.224109) | `diamond_benedict_2014_double_debye_4` | This is a theoretical multiphase carbon EOS. It publishes the fitted diamond model coefficients, but not the underlying electronic-structure grid as row-level data. |
| [Caracas et al. (2005)](https://doi.org/10.1029/2004gl022144) | `ca_perovskite_caracas_2005_bm3_3` | The article publishes fitted EOS coefficients and relative energies but not the first-principles E(V) observations, fit weights, residuals, or covariance. No supporting-information or official data attachment is listed on the publisher article page or the UCL deposit, so an independent coefficient refit is not possible. |
| [Correa et al. (2008)](https://doi.org/10.1103/physrevb.78.024101) | `diamond_correa_2008_double_debye_log_moment_5` | This is a theoretical multiphase carbon EOS. It publishes fitted cold-curve and thermal-model coefficients, but not the underlying DFT energy-volume grid as row-level data. |
| [Datchi et al. (2007)](https://doi.org/10.1103/physrevb.75.214104) | `diamond_datchi_2007_vinet_1` | For the diamond record, Datchi et al. reanalyze the previously published Occelli et al. diamond compression data on the H05 pressure scale and report the resulting EOS parameters in Table II; they do not republish the row-level diamond observations. The paper's new c-BN observations are separately bundled with the c-BN record. |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.77.094106) | `diamond_benedict_2014_dewaele_anchored`, `diamond_correa_2008_dewaele_anchored` | The linked diffraction rows constrain only the Dewaele reference isotherm; the Benedict thermal term is a separately published theoretical model. The linked diffraction rows constrain only the Dewaele reference isotherm; the Correa thermal term is a separately published theoretical model. |
| [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/physrevb.75.024115) | `platinum_dorogokupets_oganov_2007_vinet_4` | This is a semiempirical multi-material pressure-scale construction from published shock, ultrasonic, X-ray, and thermochemical literature. It publishes EOS coefficients and calculated calibration values, but no new row-level experimental platinum observations. |
| [Duffy and Ahrens (1995)](https://doi.org/10.1029/94jb02065) | `mgo_b1_duffy_ahrens_1995_hugoniot_5` | The published phase-specific coefficients are transcribed directly. The article's observation table is not redistributed because no open table-data license was identified. |
| [Fortes (2019)](https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf) | `lead_fcc_fortes_2019_bm4_1` | Fortes (2019) derives an fcc-Pb pressure scale from published literature data and tabulates model coefficients and comparisons, but no new row-level experimental P-V-T observations. |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | `e_feooh_gleason_2008_bm2_1` | Only 1 observation(s) lie at the reference temperature for 2 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |
| [Guigue et al. (2020)](https://doi.org/10.1063/1.5138697) | `palladium_guigue_2020_vinet_1` | The underlying pure-Pd observations are plotted but not tabulated in the accessible primary article; no numerical refit is claimed. |
| [Holmes et al. (1989)](https://doi.org/10.1063/1.344177) | `platinum_holmes_1989_vinet_1` | The bundled rows are shock-Hugoniot qualification experiments; the stored equilibrium Vinet curve is a theoretical 300 K isotherm and cannot be refitted directly to those rows. |
| [Kawai and Tsuchiya (2014)](https://doi.org/10.1002/2013jb010905) | `ca_perovskite_kawai_2014_vinet_mgd_3` | The publisher page exposes no supporting-information or data file, and the article plots but does not tabulate the underlying FPMD P-V-T stress averages. A direct refit is therefore impossible. All 60 printed Table 1 fitted-isotherm benchmark states are bundled separately for numerical reproduction; they are model values, not primary observations. |
| [Li et al. (2006)](https://doi.org/10.1029/2005jb004251) | `mgo_li_2006_bm3_absolute_acoustic` | The Table 1 pressures are outputs of the stored acoustic-derived BM3, not independent pressure-volume observations. The source-derived isothermal coefficients are instead validated by the bundled velocity-density data and the dedicated acoustic finite-strain reproduction. |
| [Luo et al. (2023)](https://doi.org/10.1103/physrevb.107.134116) | `mgo_b1_luo_2023_vinet_thermal_5` | The five bundled Table I rows are only the new shock subset of a global quasi-Debye fit. The complete earlier-study observations, numerical sound-velocity-density fits, objective weights, and covariance are not published; Tables II-III are derived EOS output and cannot serve as independent refit observations. |
| [Mosenfelder et al. (2009)](https://doi.org/10.1029/2008jb005900) | `mgsio3_post_perovskite_mosenfelder_2009_bm3_1` | The bundled rows are shock states and the source's thermal reduction cannot be reconstructed as a direct P-V-T least-squares fit because most rows do not report temperature. |
| [Muñoz and Kunc (1993)](https://doi.org/10.1088/0953-8984/5/33/010) | `indium_nitride_munoz_1993_murnaghan_1` | This is a first-principles study. The calculated E(V) points are plotted but not tabulated; Table 1 contains only the fitted theoretical parameters. |
| [Noguchi et al. (1999)](https://doi.org/10.1016/s0022-3697(98)00296-0) | `nickel_oxide_noguchi_1999_bm3_1` | The bundled rows are Hugoniot states; the stored 300 K isotherm is the source's Mie-Gruneisen reduction, not a direct fit to Hugoniot P-V pairs. |
| [Noguchi et al. (2013)](https://doi.org/10.1007/s00269-012-0549-1) | `ca_perovskite_noguchi_2013_bm2_mgd_1` | Table 1 contains 54 P-V-T rows, paired Fei/Holmes pressures, Pt lattice parameters, and three explicit fit exclusions. A complete local transcription was used for the audit and independent refit, but is not redistributed because the subscription article states no reusable data license. |
| [Shen and Smith (2026)](https://doi.org/10.1103/fxgq-96sg) | `fe_shen_2026_vinet_1`, `gold_shen_2026_vinet_3`, `iron_shen_2026_vinet_2`, `mgo_shen_2026_vinet_3`, `molybdenum_shen_2026_vinet_1`, `nacl_b1_shen_2026_vinet_1`, `nacl_b2_shen_2026_vinet_2`, `platinum_shen_2026_vinet_2`, `tantalum_shen_2026_vinet_2`, `tungsten_shen_2026_vinet_3` | The workbook contains simultaneous volumes but no pressures, and the record declares its Cu anchor as reference_model_not_supported. |
| [Sokolova et al. (2013)](https://doi.org/10.1016/j.rgg.2013.01.005) | `aluminum_sokolova_2013_holzapfel_2`, `copper_sokolova_2013_holzapfel_2`, `diamond_sokolova_2013_holzapfel_3`, `gold_sokolova_2013_holzapfel_4`, `mgo_sokolova_2013_holzapfel_4`, `molybdenum_sokolova_2013_holzapfel_2`, `niobium_sokolova_2013_holzapfel_2`, `platinum_sokolova_2013_holzapfel_3`, `silver_sokolova_2013_holzapfel_2`, `tantalum_sokolova_2013_holzapfel_3`, `tungsten_sokolova_2013_holzapfel_4` | This is an internally consistent multi-marker optimization. It publishes input constants and optimized EOS coefficients, but no new row-level experimental P-V-T observations; the calibration comparisons are graphical. |
| [Sun et al. (2016)](https://doi.org/10.1002/2016jb013062) | `ca_perovskite_sun_2016_bm3_3` | The published thermal-EOS coefficients are transcribed directly. The article's P-V-T table is not redistributed because no open table-data license was identified. |
| [Sun et al. (2022)](https://doi.org/10.2138/am-2021-7913) | `ca_perovskite_tetragonal_sun_2022_bm3_1` | The published fixed-derivative BM3 coefficients are transcribed directly. The article's P-V table is not redistributed because no open table-data license was identified. |
| [Tange et al. (2009)](https://doi.org/10.1029/2008jb005813) | `mgo_b1_tange_2009_vinet` | This is a unified least-squares analysis of previously published pressure-scale-free thermal, elastic, and shock datasets. It reports optimized MgO EOS parameters and residuals, but no new row-level experimental observations. |
| [Zhao et al. (1997)](https://doi.org/10.1029/96gl03769) | `naalsi2o6_zhao_1997_bm3_1` | Only 1 observation(s) lie at the reference temperature for 1 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |

## Complete investigated-paper register

This is the exhaustive paper-level index. `Bundled` means numerical primary
rows are stored; `plot only` means observations were digitized;
`parameterization only` means only the published equation/coefficients can
be checked. Record-level links, fit metrics, and evidence locations are in
the primary-source and refit ledgers.

| Paper | Final disposition | Catalog records | Record-level results | Primary-data form |
|---|---|---:|---|---|
| [Anderson et al. (1989)](https://doi.org/10.1063/1.342969) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Andrault et al. (2003)](https://doi.org/10.2138/am-2003-2-307) | Reproduced | 1 | 1 parity | 1 bundled |
| [Angel and Jackson (2002)](https://doi.org/10.2138/am-2002-0419) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Angel et al. (1997)](https://doi.org/10.1107/s0021889897000861) | Reproduced | 1 | 1 parity | 1 bundled |
| [Anzellini et al. (2014)](https://doi.org/10.1063/1.4863300) | Reproduced | 1 | 1 parity | 1 bundled |
| [Anzellini et al. (2019)](https://doi.org/10.1038/s41598-019-51037-8) | Reproduced | 1 | 1 similar | 1 bundled |
| [Anzellini et al. (2019)](https://doi.org/10.1038/s41598-019-51931-1) | Mixed: reproduced and discrepant records | 4 | 1 parity; 2 similar; 1 parity not achieved | 4 bundled |
| [Anzellini et al. (2022)](https://doi.org/10.1038/s41598-022-10523-2) | Reproduced | 1 | 1 parity | 1 bundled |
| [Anzellini et al. (2025)](https://doi.org/10.1038/s43246-025-00963-4) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [B1 Fe0.94O, Fischer et al. (2011)](https://doi.org/10.1016/j.epsl.2011.02.025) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 bundled |
| [Baty et al. (2024)](https://doi.org/10.1063/5.0179469) | Coefficient parity not achieved | 2 | 1 parity not achieved; 1 direct refit unavailable | 1 bundled; 1 theoretical parameterization only |
| [Bejina et al. (2021)](https://doi.org/10.5194/ejm-33-519-2021) | Reproduced | 1 | 1 parity | 1 bundled |
| [Benedict et al. (2014)](https://doi.org/10.1103/physrevb.89.224109) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Besson et al. (1994)](https://doi.org/10.1103/physrevb.49.12540) | Reproduced | 1 | 1 similar | 1 bundled |
| [Bezacier et al. (2014)](https://doi.org/10.1063/1.4894421) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Boffa Ballaran et al. (2007)](https://doi.org/10.2138/am.2007.2715) | Reproduced | 2 | 2 parity | 2 bundled |
| [Campbell and Heinz (1991)](https://doi.org/10.1016/0022-3697(91)90181-x) | Reproduced | 1 | 1 parity | 1 bundled |
| [Campbell and Heinz (1993)](https://doi.org/10.1016/0022-3697(93)90106-2) | Reproduced | 1 | 1 parity | 1 bundled |
| [Campbell and Heinz (1994)](https://doi.org/10.1029/94jb00127) | Reproduced | 2 | 2 parity | 2 bundled |
| [Caracas et al. (2005)](https://doi.org/10.1029/2004gl022144) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Chidester et al. (2021)](https://doi.org/10.1103/physrevb.104.094107) | Reproduced | 1 | 1 parity | 1 bundled |
| [Clendenen and Drickamer (1966)](https://doi.org/10.1063/1.1726610) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Correa et al. (2008)](https://doi.org/10.1103/physrevb.78.024101) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Crichton et al. (2002)](https://doi.org/10.2138/am-2002-2-316) | Reproduced | 1 | 1 parity | 1 bundled |
| [Criniti et al. (2023)](https://doi.org/10.2138/am-2022-8559) | Reproduced | 2 | 2 similar | 2 bundled |
| [Cynn and Yoo (1999)](https://doi.org/10.1103/physrevb.59.8526) | Reproduced | 1 | 1 parity | 1 bundled |
| [Daniel et al. (2004)](https://doi.org/10.1029/2004gl020213) | Reproduced | 2 | 2 parity | 2 bundled |
| [data from JCPDS and Levien and Prewitt, 1981](https://msaweb.org/AmMin/AM66/AM66_324.pdf) | Reproduced | 1 | 1 similar | 1 bundled |
| [Datchi et al. (2007)](https://doi.org/10.1103/physrevb.75.214104) | Partly reproduced | 3 | 1 parity; 1 similar; 1 direct refit unavailable | 2 bundled; 1 parameterization only |
| [Dewaele (2019)](https://doi.org/10.3390/min9110684) | Reproduced | 1 | 1 parity | 1 bundled |
| [Dewaele and Torrent (2013)](https://doi.org/10.1103/physrevb.88.064107) | Reproduced | 1 | 1 parity | 1 bundled |
| [Dewaele et al. (2000)](https://doi.org/10.1029/1999jb900364) | Reproduced | 1 | 1 parity | 1 bundled |
| [Dewaele et al. (2004)](https://doi.org/10.1103/physrevb.70.094112) | Reproduced | 4 | 4 parity | 4 bundled |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.77.094106) | Partly reproduced | 3 | 1 parity; 2 direct refit unavailable | 3 bundled |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.78.104102) | Reproduced | 2 | 2 similar | 2 bundled |
| [Dewaele et al. (2012)](https://doi.org/10.1103/physrevb.85.214105) | Reproduced | 4 | 4 similar | 4 bundled |
| [Dewaele et al. (2015)](https://doi.org/10.1103/physrevb.91.134108) | Reproduced | 2 | 2 parity | 2 bundled |
| [Dobrosavljevic et al. (2019)](https://doi.org/10.3390/min9120762) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/physrevb.75.024115) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Dubrovinsky et al. (2002)](https://doi.org/10.1080/08957950212807) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Duffy and Ahrens (1995)](https://doi.org/10.1029/94jb02065) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Duffy et al. (1995)](https://doi.org/10.1103/physrevlett.74.1371) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Fedotenko et al. (2020)](https://doi.org/10.1016/j.jallcom.2020.156179) | Reproduced | 2 | 2 parity | 2 bundled |
| [Fei et al. (2000)](https://doi.org/10.2138/am-2000-11-1229) | Reproduced | 1 | 1 similar | 1 plot only/digitized |
| [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) | Reproduced | 3 | 2 parity; 1 similar | 3 plot only/digitized |
| [Finkelstein et al. (2014)](https://doi.org/10.2138/am.2014.4526) | Reproduced | 1 | 1 parity | 1 bundled |
| [Finkelstein et al. (2017)](https://doi.org/10.2138/am-2017-5966) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Fortes (2019)](https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Frank et al. (2004)](https://doi.org/10.1016/j.gca.2003.12.007) | Reproduced | 1 | 1 parity | 1 bundled |
| [Fratanduono et al. (2021)](https://doi.org/10.1126/science.abh0364) | Reproduced | 1 | 1 parity | 1 bundled |
| [Frost et al. (2023)](https://doi.org/10.1063/5.0161038) | Reproduced | 2 | 2 similar | 2 bundled |
| [Fu et al. (2024)](https://doi.org/10.2138/am-2023-8969) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Fujihisa and Takemura (1996)](https://doi.org/10.1103/physrevb.54.5) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Gerward et al. (2005)](https://doi.org/10.1016/j.jallcom.2005.04.008) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | Coefficient parity not achieved | 2 | 1 parity not achieved; 1 direct refit unavailable | 2 bundled |
| [Guigue et al. (2020)](https://doi.org/10.1063/1.5138697) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 plot only/digitized |
| [Haines et al. (2001)](https://doi.org/10.1088/0953-8984/13/11/303) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Hanfland et al. (1989)](https://doi.org/10.1103/physrevb.39.12598) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Hanfland et al. (1999)](https://doi.org/10.1016/s0038-1098(99)00322-1) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hanna et al. (2011)](https://doi.org/10.1063/1.3644969) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Hazen and Finger (1979)](https://msaweb.org/AmMin/AM64/AM64_196.pdf) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hazen and Finger (1981)](https://doi.org/10.1016/0022-3697(81)90074-3) | Reproduced | 2 | 2 similar | 2 bundled |
| [Heinz and Jeanloz (1984)](https://doi.org/10.1103/physrevb.30.6045) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hemley et al. (1989)](https://doi.org/10.1103/physrevb.39.11820) | Reproduced | 1 | 1 parity | 1 bundled |
| [Holmes et al. (1989)](https://doi.org/10.1063/1.344177) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Jacobsen et al. (2002)](https://doi.org/10.1029/2001jb000490) | Coefficient parity not achieved | 3 | 3 parity not achieved | 3 bundled |
| [Jacobsen et al. (2005)](https://doi.org/10.1107/s0909049505022326) | Mixed: reproduced and discrepant records | 3 | 1 parity; 2 parity not achieved | 3 bundled |
| [Jacobsen et al. (2008)](https://doi.org/10.2138/am.2008.2988) | Reproduced | 2 | 2 parity | 2 bundled |
| [Katsura et al. (2004), Mg2SiO4 ringwoodite](https://doi.org/10.1029/2004JB003094) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Katsura et al. (2009)](https://doi.org/10.1029/2009gl038107) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Kawai and Tsuchiya (2014)](https://doi.org/10.1002/2013jb010905) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Knittle and Jeanloz (1987)](https://doi.org/10.1126/science.235.4789.668) | Reproduced | 1 | 1 similar | 1 bundled |
| [Knorr et al. (2003)](https://doi.org/10.1140/epjb/e2003-00034-6) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Koemets et al. (2023)](https://doi.org/10.3389/fchem.2023.1258389) | Reproduced | 1 | 1 parity | 1 bundled |
| [Kubo et al. (2000)](https://doi.org/10.2183/pjab.76.103) | Reproduced | 2 | 2 parity | 2 bundled |
| [Kubo et al. (2006)](https://doi.org/10.1029/2006gl025686) | Reproduced | 2 | 2 parity | 2 bundled |
| [Le Godec et al. (2014)](https://doi.org/10.3103/s1063457614010092) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Li et al. (2006)](https://doi.org/10.1029/2005jb004251) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Liu and Bassett (1973)](https://doi.org/10.1029/jb078i035p08470) | Reproduced | 1 | 1 parity | 1 bundled |
| [Luo et al. (2023)](https://doi.org/10.1103/physrevb.107.134116) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Lv et al. (2020)](https://doi.org/10.2138/am-2020-7279) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Magad-Weiss et al. (2021)](https://doi.org/10.1103/physrevb.103.014101) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Mao et al. (1974)](https://doi.org/10.1029/jb079i008p01165) | Reproduced | 1 | 1 parity | 1 bundled |
| [Mao et al. (1989)](https://doi.org/10.1029/jb094ib12p17889) | Reproduced | 1 | 1 parity | 1 bundled |
| [Mao et al. (1990)](https://doi.org/10.1029/jb095ib13p21737) | Reproduced | 1 | 1 similar | 1 bundled |
| [Mao et al. (2015)](https://doi.org/10.1002/2015gl064400) | Reproduced | 1 | 1 similar | 1 bundled |
| [Martinez et al. (1996)](https://doi.org/10.2138/am-1996-5-608) | Reproduced | 1 | 1 parity | 1 bundled |
| [Matsui et al. (2012)](https://doi.org/10.2138/am.2012.3937) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [McHardy et al. (2026)](https://doi.org/10.1103/zp3m-kjpc) | Reproduced | 3 | 3 parity | 3 bundled |
| [Meng et al. (1994)](https://doi.org/10.1007/bf00203299) | Reproduced | 1 | 1 parity | 1 bundled |
| [Milani et al. (2015)](https://doi.org/10.1016/j.lithos.2015.03.017) | Reproduced | 2 | 2 parity | 2 bundled |
| [Miozzi et al. (2018)](https://doi.org/10.1029/2018je005582) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Mookherjee et al. (2019)](https://doi.org/10.2138/am-2019-6694) | Reproduced | 1 | 1 parity | 1 bundled |
| [Mosenfelder et al. (2009)](https://doi.org/10.1029/2008jb005900) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Muhammad et al. (2024)](https://doi.org/10.1039/d4nr00093e) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Muñoz and Kunc (1993)](https://doi.org/10.1088/0953-8984/5/33/010) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Noguchi et al. (1999)](https://doi.org/10.1016/s0022-3697(98)00296-0) | Partly reproduced | 2 | 1 parity; 1 direct refit unavailable | 2 bundled |
| [Noguchi et al. (2013)](https://doi.org/10.1007/s00269-012-0549-1) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Ono et al. (2000)](https://doi.org/10.1007/s002690000108) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Ono et al. (2006)](https://doi.org/10.2138/am.2006.2347) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ono et al. (2006)](https://doi.org/10.2138/am.2006.2118) | Reproduced | 3 | 3 parity | 3 bundled |
| [Pepin et al. (2014)](https://doi.org/10.1103/physrevlett.113.265504) | Reproduced | 2 | 2 parity | 2 plot only/digitized |
| [Prescher et al. (2015)](https://doi.org/10.1038/ngeo2370) | Reproduced | 1 | 1 similar | 1 bundled |
| [Qin et al. (2023)](https://doi.org/10.2138/am-2022-8432) | Reproduced | 2 | 2 parity | 2 bundled |
| [Redfern and Angel (1999)](https://doi.org/10.1007/s004100050471) | Reproduced | 1 | 1 parity | 1 bundled |
| [Reynard et al. (1996)](https://doi.org/10.2138/am-1996-1-206) | Reproduced | 2 | 2 parity | 2 bundled |
| [Richet et al. (1988)](https://doi.org/10.1029/jb093ib12p15279) | Reproduced | 2 | 2 parity | 2 bundled |
| [Richet et al. (1989)](https://doi.org/10.1029/jb094ib03p03037) | Reproduced | 1 | 1 parity | 1 bundled |
| [Rodrigo-Ramon et al. (2024)](https://doi.org/10.1038/s41598-024-78006-0) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ross (1997)](https://doi.org/10.2138/am-1997-7-805) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ross and Angel (1999)](https://doi.org/10.2138/am-1999-0309) | Reproduced | 2 | 2 parity | 2 bundled |
| [Sakai et al. (2016)](https://doi.org/10.1038/srep22652) | Reproduced | 1 | 1 parity | 1 bundled |
| [Sato and Jeanloz (1981)](https://doi.org/10.1029/jb086ib12p11773) | Reproduced | 1 | 1 parity | 1 bundled |
| [Saxena et al. (1999)](https://doi.org/10.2138/am-1999-0303) | Reproduced | 1 | 1 similar | 1 bundled |
| [Schouwink et al. (2011)](https://doi.org/10.2138/am.2011.3775) | Reproduced | 1 | 1 parity | 1 bundled |
| [Schulze et al. (2018)](https://doi.org/10.2138/am-2018-6562) | Reproduced | 1 | 1 parity | 1 bundled |
| [Scott et al. (2001)](https://doi.org/10.1029/2000gl012606) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Shen and Smith (2026)](https://doi.org/10.1103/fxgq-96sg) | Direct refit unavailable | 10 | 10 direct refit unavailable | 10 bundled |
| [Shi et al. (2022)](https://doi.org/10.1029/2021jb023805) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Shieh et al. (2000)](https://doi.org/10.1016/s0012-821x(00)00033-9) | Reproduced | 2 | 2 parity | 2 bundled |
| [Shim et al. (2000)](https://doi.org/10.1029/2000jb900183) | Reproduced | 1 | 1 similar | 1 bundled |
| [Shim et al. (2000)](https://doi.org/10.1016/s0031-9201(00)00154-0) | Reproduced | 1 | 1 similar | 1 bundled |
| [Siersch et al. (2021)](https://doi.org/10.1016/j.pepi.2021.106786) | Reproduced | 1 | 1 parity | 1 bundled |
| [Sokolova et al. (2013)](https://doi.org/10.1016/j.rgg.2013.01.005) | Direct refit unavailable | 11 | 11 direct refit unavailable | 11 parameterization only |
| [Somayazulu et al. (2023)](https://doi.org/10.1098/rsta.2022.0331) | Mixed: reproduced and discrepant records | 3 | 1 parity; 1 similar; 1 parity not achieved | 3 bundled |
| [Speziale et al. (2001)](https://doi.org/10.1029/2000jb900318) | Reproduced | 1 | 1 parity | 1 bundled |
| [Stinton et al. (2014)](https://doi.org/10.1103/physrevb.90.134105) | Reproduced | 2 | 2 parity | 2 plot only/digitized |
| [Sun et al. (2016)](https://doi.org/10.1002/2016jb013062) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Sun et al. (2019)](https://doi.org/10.1029/2019jb017853) | Reproduced | 2 | 2 similar | 2 bundled |
| [Sun et al. (2022)](https://doi.org/10.2138/am-2021-7913) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Suzuki (2016)](https://doi.org/10.2465/jmps.160719c) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura (2004)](https://doi.org/10.1103/physrevb.70.012101) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura and Dewaele (2008)](https://doi.org/10.1103/physrevb.78.104119) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura and Singh (2006)](https://doi.org/10.1103/physrevb.73.224119) | Reproduced | 1 | 1 parity | 1 bundled |
| [Tange et al. (2009)](https://doi.org/10.1029/2008jb005813) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Tange et al. (2012)](https://doi.org/10.1029/2011jb008988) | Reproduced | 2 | 2 similar | 2 bundled |
| [Tateno et al. (2019)](https://doi.org/10.2138/am-2019-6779) | Reproduced | 1 | 1 parity | 1 bundled |
| [Thompson et al. (2017)](https://doi.org/10.1002/2017jb014168) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Vanpeteghem et al. (2006)](https://doi.org/10.1029/2005gl024955) | Reproduced | 1 | 1 parity | 1 bundled |
| [Vocadlo (1999)](https://doi.org/10.2138/am-1999-1017) | Reproduced | 1 | 1 similar | 1 bundled |
| [Walker et al. (2002)](https://doi.org/10.2138/am-2002-0701) | Reproduced | 1 | 1 similar | 1 bundled |
| [Wang et al. (2012)](https://doi.org/10.1029/2011jb009100) | Reproduced | 1 | 1 parity | 1 bundled |
| [Wang et al. (2026), KAlSi3O8 liebermannite and K-hollandite II](https://doi.org/10.2138/am-2024-9562) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Wittlinger et al. (1997)](https://doi.org/10.1107/s0108768197005739) | Reproduced | 1 | 1 similar | 1 plot only/digitized |
| [Wolf et al. (2015)](https://doi.org/10.1002/2015jb012108) | Reproduced | 2 | 2 parity | 2 bundled |
| [Xu et al. (2020)](https://doi.org/10.1029/2020gl088877) | Reproduced | 1 | 1 parity | 1 bundled |
| [Yagi et al. (1992)](https://doi.org/10.1016/0031-9201(92)90063-2) | Reproduced | 1 | 1 similar | 1 bundled |
| [Yu et al. (2024)](https://doi.org/10.1029/2023jb028026) | Reproduced | 1 | 1 similar | 1 bundled |
| [Zhao et al. (1997)](https://doi.org/10.1029/96gl03769) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Zhu et al. (2020)](https://doi.org/10.1029/2020jb019964) | Reproduced | 3 | 3 parity | 3 bundled |

## Maintenance

This page is generated by `scripts/generate_paper_investigation_ledger.py`.
Update record-level evidence first, add nonproduction investigations to
`docs/data/nonproduction-paper-investigations.json`, regenerate this page,
and run the generator with `--check` in validation.
