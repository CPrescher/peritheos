# References

- Correa, A. A., Benedict, L. X., Young, D. A., Schwegler, E. & Bonev, S. A.
  (2008). First-principles multiphase
  equation of state of carbon under extreme conditions. *Physical Review B*,
  78, 024101. [doi:10.1103/PhysRevB.78.024101](https://doi.org/10.1103/PhysRevB.78.024101).
  `DoubleDebyeLogMomentHelmholtz` implements equations 2--7 and 13--18 for the
  diamond branch, including logarithmic-moment weights and the constant
  anharmonic coefficient from Table I.

- Benedict, L. X., Driver, K. P., Hamel, S., Militzer, B., Qi, T., Correa,
  A. A., Saul, A. & Schwegler, E. (2014). A multiphase equation of state for
  carbon addressing high pressures and temperatures. *Physical Review B*,
  89, 224109. [doi:10.1103/PhysRevB.89.224109](https://doi.org/10.1103/PhysRevB.89.224109).
  `DoubleDebyeHelmholtz` implements equations 3--7 for the double-Debye solid,
  including the Vinet cold energy, zero point, and the Table I diamond
  anharmonic coefficient.

- Bezacier, L., Journaux, B., Perrillat, J.-P., Cardon, H., Hanfland, M. &
  Daniel, I. (2014). Equations of state of ice VI and ice VII at high pressure
  and high temperature. *Journal of Chemical Physics*, 141, 104505.
  doi:10.1063/1.4894421. Peritheos uses equations 1--3, the H2O PVT rows in
  Table II, and representative measured states from Table I.
- Anzellini, S., Dewaele, A., Occelli, F., Loubeyre, P. & Mezouar, M.
  (2014). Equation of state of rhenium and application for ultra high pressure
  calibration. *Journal of Applied Physics*, 115, 043511.
  doi:10.1063/1.4863300. Peritheos uses equation 6, the fitted values and 95%
  confidence intervals following it, Table IV's 0.64-144 GPa domain, and
  lattice-data regressions from Table III.
- Anderson, O. L., Isaak, D. G. & Yamamoto, S. (1989). Anharmonicity and the
  equation of state for gold. *Journal of Applied Physics*, 65, 1534--1543.
  [doi:10.1063/1.342969](https://doi.org/10.1063/1.342969). Peritheos uses
  Equation 29 and the top rows of Table V, including the logarithmic-volume
  correction to the thermal-pressure slope. The ambient reference volume is
  converted from the density and atomic mass stated on page 1535.
- Birch, F. (1947). Finite elastic strain of cubic crystals. *Physical Review*,
  71, 809-824.
- Clendenen, R. L. & Drickamer, H. G. (1966). Lattice parameters of nine
  oxides and sulfides as a function of pressure. *Journal of Chemical Physics*,
  44, 4223--4228. [doi:10.1063/1.1726610](https://doi.org/10.1063/1.1726610).
  The CoO record uses Equation 4 and Tables II, III, and VI. The published
  Murnaghan constants are a whole-range empirical fit; no parameter errors or
  covariance are printed.
- Dewaele, A., Fiquet, G. & Gillet, P. (1998). Temperature and pressure
  distribution in the laser-heated diamond-anvil cell. *Review of Scientific
  Instruments*, 69, 2421-2426. doi:10.1063/1.1148970.
- Datchi, F., Dewaele, A., Le Godec, Y. & Loubeyre, P. (2007). Equation of
  state of cubic boron nitride at high pressures and temperatures. *Physical
  Review B*, 75, 214104. doi:10.1103/PhysRevB.75.214104. Peritheos uses the
  295 K Vinet fit from Table I and the volume convention in Figure 1.
- Dewaele, A., Datchi, F., Loubeyre, P. & Mezouar, M. (2008). High
  pressure-high temperature equations of state of neon and diamond. *Physical
  Review B*, 77, 094106. doi:10.1103/PhysRevB.77.094106. The diamond entry
  follows equations 2, 3, and 6 and Table III.
- Dewaele, A., Torrent, M., Loubeyre, P. & Mezouar, M. (2008). Compression
  curves of transition metals in the Mbar range: Experiments and projector
  augmented-wave calculations. *Physical Review B*, 78, 104102.
  doi:10.1103/PhysRevB.78.104102. Peritheos uses the unconstrained Ag and Ni
  fits in Table IV and regression rows from Table II.
- Dewaele, A., Belonoshko, A. B., Garbarino, G., Occelli, F., Bouvier, P.,
  Hanfland, M. & Mezouar, M. (2012). High-pressure-high-temperature equation
  of state of KCl and KBr. *Physical Review B*, 85, 214105.
  doi:10.1103/PhysRevB.85.214105. Peritheos uses equation 2 and Tables III and
  V, with separate B1 and B2 records.
- Campbell, A. J. & Heinz, D. L. (1991). Compression of KCl in the B2
  structure to 56 GPa. *Journal of Physics and Chemistry of Solids*, 52,
  495--499. doi:10.1016/0022-3697(91)90181-X. The material-library composite
  uses the primary abstract's B2/B1 `V0` ratio, `K0`, fixed `K0'`, and range;
  its absolute `V0` is explicitly combined with Dewaele et al.'s B1 value.
- Campbell, A. J. & Heinz, D. L. (1993). Equation of state and high pressure
  phase transition of NiS in the NiAs structure. *Journal of Physics and
  Chemistry of Solids*, 54, 5--7.
  [doi:10.1016/0022-3697(93)90106-2](https://doi.org/10.1016/0022-3697(93)90106-2).
  Peritheos uses the ambient lattice constants, Table 1 compression data, and
  the BM3 coefficients printed on page 6, including all reported errors.
- Campbell, A. J. & Heinz, D. L. (1994). High-pressure acoustic wave
  velocities and equations of state of the alkali chlorides. *Journal of
  Geophysical Research*, 99, 11765--11774.
  [doi:10.1029/94JB00127](https://doi.org/10.1029/94JB00127). The CsCl and
  RbCl-B2 records use the separate Table 1 compression blocks and finite-strain
  fits on pages 11767--11768. CsCl uses the accepted ambient lattice parameter
  and Yagi's corrected Table 1 ratios; non-quenchable RbCl-B2 uses the paper's
  hypothetical zero-pressure density.
- Yagi, T. (1978). Experimental determination of thermal expansivity of several
  alkali halides at high pressures. *Journal of Physics and Chemistry of
  Solids*, 39, 563--571.
  [doi:10.1016/0022-3697(78)90037-9](https://doi.org/10.1016/0022-3697(78)90037-9).
  The CsCl validation uses the nine 25 degC `V/V0` values in Table 1 and applies
  Campbell and Heinz's documented 4.118-to-4.123 angstrom reference correction.
- Dewaele, A. (2019). Equations of State of Simple Solids (Including Pb, NaCl
  and LiF) Compressed in Helium or Neon in the Mbar Range. *Minerals*, 9, 684.
  doi:10.3390/min9110684. Peritheos uses equation 1, the unified
  Dorogokupets-calibrated column of Table 1, and Tables 3 and 4.
- Dorogokupets, P. I. (2010). P-V-T equations of state of MgO and
  thermodynamics. *Physics and Chemistry of Minerals*, 37, 677-684.
  doi:10.1007/s00269-010-0367-2.
- Dorogokupets, P. I., Sokolova, T. S., Danilov, B. S. & Litasov, K. D.
  (2012). Near-absolute equations of state of diamond, Ag, Al, Au, Cu, Mo, Nb,
  Pt, Ta, and W for quasi-hydrostatic conditions. *Geodynamics &
  Tectonophysics*, 3, 129--166.
  [doi:10.5800/GT-2012-3-2-0067](https://doi.org/10.5800/GT-2012-3-2-0067).
  This paper documents the simultaneous-optimization methodology that precedes
  the final Sokolova et al. diamond pressure scale.
- Dorfman, S. M., Prakapenka, V. B., Meng, Y. & Duffy, T. S. (2012).
  Intercomparison of pressure standards (Au, Pt, Mo, MgO, NaCl and Ne) to
  2.5 Mbar. *Journal of Geophysical Research: Solid Earth*, 117, B08210.
  doi:10.1029/2012JB009292. Peritheos uses equation 2, Tables 1 and 2, and
  sections 3.1, 3.2, and 4 for the 300 K Vinet catalog entries.
- Fei, Y., Ricolleau, A., Frank, M., Mibe, K., Shen, G. & Prakapenka, V.
  (2007). Toward an internally consistent pressure scale. *Proceedings of the
  National Academy of Sciences*, 104, 9182-9186.
  doi:10.1073/pnas.0609013104. Peritheos uses equation 3, its stated
  Debye-temperature convention, Table 1, and the data envelopes in Figures
  1-5 for the Au, Pt, NaCl-B2, and Ne thermal scales.
- Muñoz, A. & Kunc, K. (1993). Structure and static properties of indium
  nitride at low and moderate pressures. *Journal of Physics: Condensed
  Matter*, 5, 6015--6022. doi:10.1088/0953-8984/5/33/010. Peritheos uses the
  Murnaghan energy-volume formulation in section 3 and the theoretical
  wurtzite lattice and EOS coefficients in Table 1.
- Shen, G. & Smith, J. S. (2026). Simultaneous x-ray diffraction measurements
  of nine pressure calibrants to 140 GPa. *Physical Review B*, 113, 144113.
  [doi:10.1103/fxgq-96sg](https://doi.org/10.1103/fxgq-96sg). Peritheos uses
  the reduced-300 K Vinet form in Equation 4, the fixed reference volumes and
  fitted parameters in Table II, and the phase/range qualifications in Table I
  and Section III.E. The printed parameter errors are retained without an
  inferred confidence level or covariance.
- Frank, M. R., Fei, Y. & Hu, J. (2004). Constraining the equation of state of
  fluid H2O to 80 GPa using the melting curve, bulk modulus, and thermal
  expansivity of ice VII. *Geochimica et Cosmochimica Acta*, 68, 2781--2790.
  doi:10.1016/j.gca.2003.12.007. The catalog uses Equation 2 and the
  simultaneous three-parameter 300 K fit in section 3.1; Table 1 supplies an
  independent measured-state regression.
- Gerward, L., Olsen, J. S., Petit, L., Vaitheeswaran, G., Kanchana, V. &
  Svane, A. (2005). Bulk modulus of CeO2 and PrO2—an experimental and
  theoretical study. *Journal of Alloys and Compounds*, 400, 56--61.
  doi:10.1016/j.jallcom.2005.04.008. Peritheos uses Equation 1 and the
  experimental fluorite-phase BM3 fits in Tables 1 and 2.
- Fortes, A. D. (2019). *A revised equation of state for in situ pressure
  determination using fcc-Pb (0 < P < 13 GPa, T > 100 K)*. STFC Rutherford
  Appleton Laboratory Technical Report RAL-TR-2019-002. Peritheos uses
  equations 1--5 and the 300 K values in Table 1. The stable institutional
  primary copy is [RAL-TR-2019-002](https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf).
- Hazen, R. M. & Finger, L. W. (1979). Crystal structure and compressibility
  of zircon at high pressure. *American Mineralogist*, 64, 196--201.
  Peritheos uses Table 1 and the Birch--Murnaghan result on page 198:
  `V0 = 260.79(4) angstrom^3`, `K0 = 227(2) GPa`, with `K0' = 6.5` assumed.
  [Primary article](https://msaweb.org/AmMin/AM64/AM64_196.pdf).
- Lin, J.-F., Degtyareva, O., Prewitt, C. T., Dera, P., Sata, N., Gregoryanz,
  E., Mao, H.-K. & Hemley, R. J. (2004). Crystal structure of a
  high-pressure/high-temperature phase of alumina by in situ X-ray
  diffraction. *Nature Materials*, 3, 389--393.
  [doi:10.1038/nmat1121](https://doi.org/10.1038/nmat1121). The refined
  113 GPa, 300 K `Pbcn` cell, coordinates, and four-formula-unit contents are
  the structural source for the Rh2O3(II)-type alumina material.
- Holland, T. J. B. & Powell, R. (2011). An improved and extended internally
  consistent thermodynamic dataset for phases of petrological interest,
  involving a new equation of state for solids. *Journal of Metamorphic
  Geology*, 29, 333-383. doi:10.1111/j.1525-1314.2010.00923.x.
- Heinz, D. L. (1990). Thermal pressure in the laser-heated diamond anvil cell.
  *Geophysical Research Letters*, 17, 1161-1164.
  doi:10.1029/GL017i008p01161.
- Haines, J., Léger, J. M., Chateau, C. & Lowther, J. E. (2001). Experimental
  and theoretical investigation of Mo2C at high pressure. *Journal of Physics:
  Condensed Matter*, 13, 2447--2454.
  doi:10.1088/0953-8984/13/11/303. The Mo2C record uses the BM3 fit in
  section 4.1 and normalizes its reference volume from the specimen lattice
  parameters in section 2.
- Holmes, N. C., Moriarty, J. A., Gathers, G. R. & Nellis, W. J. (1989). The
  equation of state of platinum to 660 GPa (6.6 Mbar). *Journal of Applied
  Physics*, 66, 2962--2967. doi:10.1063/1.344177. Peritheos uses the
  theoretical 300 K universal isotherm in Equation 11, Table IV, and the
  approximate finite-temperature extension in Equation 12; the shock-Hugoniot
  data are retained as validation context rather than mislabeled as a static
  range.
- Holzapfel, W. B. (2001). Equations of state for solids under strong
  compression. *Zeitschrift fuer Kristallographie*, 216, 473-488.
  doi:10.1524/zkri.216.9.473.20346.
- Jackson, I. & Rigden, S. M. (1996). Analysis of P-V-T data: constraints on
  the thermoelastic properties of high-pressure minerals. *Physics of the
  Earth and Planetary Interiors*, 96, 85-112.
  doi:10.1016/0031-9201(96)03143-3.
- Levien, L. & Prewitt, C. T. (1981). High-pressure crystal structure and
  compressibility of coesite. *American Mineralogist*, 66, 324--333.
  Peritheos uses Table 7 and the unweighted Birch--Murnaghan fit on page 328;
  the conventional-cell `V0` follows the ambient lattice parameters.
  [Primary article](https://msaweb.org/AmMin/AM66/AM66_324.pdf).
- Liu, L.-G. & Bassett, W. A. (1973). Changes of the crystal structure and the
  lattice parameter of SrO at high pressure. *Journal of Geophysical Research*,
  78, 8470--8473.
  [doi:10.1029/JB078i035p08470](https://doi.org/10.1029/JB078i035p08470).
  Peritheos uses the printed Birch equation, reported ambient cell, coefficient
  errors, and all B1 data through 34.05 GPa.
- Mao, H.-K., Takahashi, T., Bassett, W. A., Kinsland, G. L. & Merrill, L.
  (1974). Isothermal compression of magnetite to 320 kbar and pressure-induced
  phase transformation. *Journal of Geophysical Research*, 79, 1165--1170.
  [doi:10.1029/JB079i008p01165](https://doi.org/10.1029/JB079i008p01165).
  Peritheos uses the ambient lattice result and third-order finite-strain fit;
  the assumed derivative and the paper's composite error convention remain
  explicit.
- Martinez, I., Zhang, J. & Reeder, R. J. (1996). In situ X-ray diffraction of
  aragonite and dolomite at high pressure and high temperature: Evidence for
  dolomite breakdown to aragonite and magnesite. *American Mineralogist*, 81,
  611--624. doi:10.2138/am-1996-5-608. Peritheos uses Equations 1--3 and Tables
  6--7 for the executable staged aragonite BM2 P-V-T parameterization. The
  separate global HT-BM3 reduction is excluded because its fitted reference
  volume is omitted and the remaining coefficients are not reproducible from
  the printed dataset. [Primary article](https://rruff.info/doclib/am/vol81/AM81_611.pdf).
- Hanfland, M., Loa, I., Syassen, K., Schwarz, U. & Takemura, K. (1999).
  Equation of state of lithium to 21 GPa. *Solid State Communications*, 112,
  123--127.
  [doi:10.1016/S0038-1098(99)00322-1](https://doi.org/10.1016/S0038-1098(99)00322-1).
  The catalog retains the paper's single empirical Vinet fit across bcc and
  fcc data and labels that cross-phase scope explicitly.
- Murnaghan, F. D. (1944). The compressibility of media under extreme
  pressures. *Proceedings of the National Academy of Sciences*, 30, 244-247.
- Noguchi, Y., Uchino, M., Hikosaka, H., Kusaba, K., Fukuoka, K., Mashimo, T.
  & Syono, Y. (1998). Shock compression of NiO to 130 GPa. *The Review of
  High Pressure Science and Technology*, 7, 832--834.
  doi:10.4131/jshpreview.7.832. Peritheos uses the ambient pseudo-cubic lattice
  parameter on page 832 and the shock-to-300 K Mie--Gruneisen reduction and
  Murnaghan--Birch fit on pages 833--834.
- Noguchi, Y., Uchino, M., Hikosaka, H., Atou, T., Kusaba, K., Fukuoka, K.,
  Mashimo, T. & Syono, Y. (1999). Equation of state of NiO studied by shock
  compression. *Journal of Physics and Chemistry of Solids*, 60, 509--514.
  doi:10.1016/S0022-3697(98)00296-0. The final article extends the isotherm to
  147.6 GPa and reports `K0 = 191 GPa` and `K0' = 3.9` without coefficient
  errors.
- Occelli, F., Loubeyre, P. & LeToullec, R. (2003). Properties of diamond
  under hydrostatic pressures up to 140 GPa. *Nature Materials*, 2, 151--154.
  [doi:10.1038/nmat831](https://doi.org/10.1038/nmat831). The 300 K
  compression data are a principal experimental constraint behind the later
  Sokolova diamond pressure scale.
- Poirier, J.-P. & Tarantola, A. (1998). A logarithmic equation of state.
  *Physics of the Earth and Planetary Interiors*, 109, 1-8.
  doi:10.1016/S0031-9201(98)00112-5.
- Ross, N. L. (1997). The equation of state and high-pressure behavior of
  magnesite. *American Mineralogist*, 82, 682--688.
  doi:10.2138/am-1997-7-805. Peritheos uses the unconstrained BM3 fit on pages
  684--685 and a measured P--V state from Table 1 as a regression case.
- Ono, S., Ito, E., Katsura, T., Yoneda, A., Walter, M. J., Urakawa, S.,
  Utsumi, W. & Funakoshi, K. (2000). Thermoelastic properties of the
  high-pressure phase of SnO2 determined by in situ X-ray observations up to
  30 GPa and 1400 K. *Physics and Chemistry of Minerals*, 27, 618--622.
  [doi:10.1007/S002690000108](https://doi.org/10.1007/S002690000108).
  The two cubic diffraction cards share the paper's 300 K BM3 fit. Its
  separate 25 GPa expansivity is documented, not promoted to a full thermal
  EOS.
- Richet, P., Mao, H.-K. & Bell, P. M. (1989). Bulk moduli of
  magnesiowuestites from static compression measurements. *Journal of
  Geophysical Research*, 94, 3037--3045.
  [doi:10.1029/JB094iB03p03037](https://doi.org/10.1029/JB094iB03p03037).
  Peritheos uses the MW60 ambient lattice parameter and second-order Birch fit.
- Sato, Y. & Jeanloz, R. (1981). Phase transition in SrO. *Journal of
  Geophysical Research*, 86, 11773--11778.
  [doi:10.1029/JB086iB12p11773](https://doi.org/10.1029/JB086iB12p11773).
  The B2 reference volume and error are converted from the reported
  extrapolated density; the derivative is fixed as in the source.
- Shi, W., Wei, W., Sun, N., Mao, Z. & Prakapenka, V. B. (2022). Thermal
  equations of state of corundum and Rh2O3 (II)-type Al2O3 up to 153 GPa and
  3400 K. *Journal of Geophysical Research: Solid Earth*, 127,
  e2021JB023805. [doi:10.1029/2021JB023805](https://doi.org/10.1029/2021JB023805).
  Peritheos uses Equations 1--6, the selected Rh2O3(II) model in Tables 1 and
  3, and all 75 P--T--V rows in official Supporting Information Table S2.
- Shieh, S. R., Mao, H.-K., Hemley, R. J. & Ming, L. C. (2000). In situ X-ray
  diffraction studies of dense hydrous magnesium silicates at mantle
  conditions. *Earth and Planetary Science Letters*, 177, 69--80.
  [doi:10.1016/S0012-821X(00)00033-9](https://doi.org/10.1016/S0012-821X(00)00033-9).
  Peritheos represents the joint BM2 modulus with separate AntA and AntB
  measured ambient volumes from Table 2.
- Scott, H. P., Williams, Q. & Knittle, E. (2001). Stability and equation of
  state of Fe3C to 73 GPa: Implications for carbon in the Earth's core.
  *Geophysical Research Letters*, 28, 1875--1878.
  doi:10.1029/2000GL012606. Peritheos uses the ambient reference volume on
  page 1875 and the weighted third-order Birch--Murnaghan fit on page 1876;
  Figure 3 documents the fitted 300 K compression data.
- Yagi, T., Uchiyama, Y., Akaogi, M. & Ito, E. (1992). Isothermal compression
  curve of MgSiO3 tetragonal garnet. *Physics of the Earth and Planetary
  Interiors*, 74, 1--7.
  [doi:10.1016/0031-9201(92)90063-2](https://doi.org/10.1016/0031-9201(92)90063-2).
  Peritheos preserves the printed fixed reference volume and derivative and
  does not manufacture parameter errors absent from the paper.
- Reeber, R. R. & Wang, K. (1996). Thermal expansion, molar volume and specific
  heat of diamond from 0 to 3000 K. *Journal of Electronic Materials*, 25,
  63--67. [doi:10.1007/BF02666175](https://doi.org/10.1007/BF02666175).
- Sokolova, T. S., Dorogokupets, P. I. & Litasov, K. D. (2013).
  Self-consistent pressure scales based on the equations of state for ruby,
  diamond, MgO, B2-NaCl, as well as Au, Pt, and other metals to 4 Mbar and
  3000 K. *Russian Geology and Geophysics*, 54, 181--199.
  [doi:10.1016/j.rgg.2013.01.005](https://doi.org/10.1016/j.rgg.2013.01.005).
  Tables 1 and 4 are the shared scientific source of Peritheos' eleven
  Sokolova marker coefficient sets.
- Sokolova, T. S., Dorogokupets, P. I., Dymshits, A. M., Danilov, B. S. &
  Litasov, K. D. (2016). Microsoft
  Excel spreadsheets for calculation of
  P-V-T relations and thermodynamic properties from equations of state of MgO,
  diamond and nine metals as pressure markers in high-pressure and
  high-temperature experiments. *Computers & Geosciences*, 94, 162-169.
  doi:10.1016/j.cageo.2016.06.002.
  The article and workbooks are the executable implementation source for the
  eleven catalog records; it is not represented as a new experimental fit.
  The exact MgO implementation also uses the 2016 correction to its earlier
  anharmonic coefficient. Peritheos follows the accompanying workbook
  calculation path when it differs from the typeset equations; see
  [Paper versus spreadsheet](equation-reference.md#paper-versus-spreadsheet).
- Tange, Y., Nishihara, Y. & Tsuchiya, T. (2009). Unified analyses for P-V-T
  equation of state of MgO: A solution for pressure-scale problems in high
  P-T experiments. *Journal of Geophysical Research: Solid Earth*, 114,
  B03208. doi:10.1029/2008JB005813. Peritheos uses equations 2, 4, 5, 15,
  and 16 and Tables 1, 2, 4, and 5 for the Fit3-Vinet pressure standard.
- Tange, Y., Nishihara, Y. & Tsuchiya, T. (2010). Correction to “Unified
  analyses for P-V-T equation of state of MgO: A solution for pressure-scale
  problems in high P-T experiments.” *Journal of Geophysical Research: Solid
  Earth*, 115, B12203. doi:10.1029/2010JB007959. The correction changes Figure
  11 only and does not alter the EOS equations, parameters, or regression table
  used by Peritheos.
- Victor, A. C. (1962). Heat capacity of diamond at high temperatures.
  *Journal of Chemical Physics*, 36, 1903--1911.
  [doi:10.1063/1.1701288](https://doi.org/10.1063/1.1701288).
- Vinet, P. and coauthors (1986, 1987). Universal/Rydberg equation of state for
  compressed solids.
- Yen, C. E., Williams, Q. & Kunz, M. (2020).
  Thermal pressure in the laser-heated diamond anvil cell: A quantitative study
  and implications for the density versus mineralogy correlation of the mantle.
  *Journal of Geophysical Research: Solid Earth*, 125, e2020JB020006.
  doi:10.1029/2020JB020006.
- Zouboulis, E. S., Grimsditch, M., Ramdas, A. K. & Rodriguez, S. (1998).
  Temperature dependence of the elastic moduli of diamond: A
  Brillouin-scattering study. *Physical Review B*, 57, 2889--2896.
  [doi:10.1103/PhysRevB.57.2889](https://doi.org/10.1103/PhysRevB.57.2889).

The DAC papers use several different normalizations for reported pressure
increases. A percentage of constant-volume thermodynamic thermal pressure is
not interchangeable with a percentage of cold pressure. Peritheos defines
`f_dac` using the former normalization, as described in
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md).

The source docstrings contain the model-specific reference nearest each
implementation. Bibliographic precision should be checked against the original
publication when citing a model in scientific work.
