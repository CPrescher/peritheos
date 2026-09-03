# Candidate material equations of state

Literature search date: 2026-09-02

This is a review backlog, not a set of validated Peritheos records. A candidate
is included when a primary paper reports a pressure-volume or
pressure-volume-temperature relation for a homogeneous solid phase, the phase
is absent from the bundled library or materially extends an existing record,
and the published formulation is plausibly representable in Peritheos.

Every accepted candidate must also be diffraction-ready for PhaseSmith and
Dioptas. Its material file must contain the phase-specific crystal system,
conventional unit cell, formula-unit count, space group, and occupied atomic
sites (or a documented fallback peak list when a defensible atomic model is
not published). Crystal structure and EOS provenance are audited separately;
an EOS is not accepted into the library with only a phase name or nominal
formula.

The current baseline is 115 material files and 159 EOS records. The model
inventory already includes BM2/BM3/BM4, Murnaghan, natural-strain, modified
Tait, Vinet, Holzapfel, referenced Mie-Gruneisen-Debye/Einstein, linear and
log-volume thermal pressure, a temperature-dependent reference-state model,
thermal modified Tait, a multi-oscillator model, and a double-Debye Helmholtz
model. The list below therefore prioritizes material parameterizations over new
generic equation families.

## Ranking key

- **A — review first:** strong scientific or calibration value, explicit
  parameters, and likely implementable with existing classes.
- **B — useful next:** worthwhile coverage, but composition, crystallography,
  reference-state conversion, or source reconciliation needs more work.
- **C — model work:** scientifically attractive, but faithful implementation
  likely needs a new phase-transition or spin-crossover formulation.

"Low effort" means no new public EOS class is expected. It does not waive the
normal primary-source, units, uncertainty, phase, and numerical reproduction
audit.

## A — review first

### C01 — Boron carbide B4C: add the thermal branch

**Status: accepted and implemented (2026-09-03).**

- **Change:** extend the existing `b4c` material rather than add a duplicate;
  add both the MGD and Berman thermal parameterizations.
- **Source:** Somayazulu et al. (2023), [P-V-T equation of state of boron
  carbide](https://doi.org/10.1098/rsta.2022.0331).
- **Published scope:** 8--50 GPa and 300--2500 K. The implemented record uses
  the quench-data BM3 (`V0 = 328.4 A^3` fixed, `K0 = 221(2) GPa`,
  `K0' = 3.3(1)`) plus the final article's MGD pair (`gamma0 = 0.8` fixed,
  fitted `q = 2.1`). `Tr = 298 K` and `theta0 = 1425 K` follow the explicit
  official supplement.
- **Mapping:** `BM3` + `MieGruneisenDebye` with the explicit
  `integrated_gruneisen` Debye-temperature law. The paper's generalized Tange
  relation reduces exactly to this model for its fixed `a = 1`.
- **Reproduction:** the official supplement includes 51 complete P-V-T rows
  with uncertainties. At their printed V and T, the final parameters give a
  1.39 GPa pressure RMSE and chi-squared per heated observation of 0.961. A
  Peritheos errors-in-variables refit of the 41 heated rows with `gamma0 = 0.8`
  fixed gives `q = 1.050(519)` and reduced chi-squared 0.767; it does not
  independently recover the article's `q = 2.1`. See the
  [full reproduction record](literature-reproductions.md#c01-boron-carbide-somayazulu-et-al-2023).
  The 289.1 A^3, 2023 K row is reproduced as 40.45 GPa versus 40.4(2.3) GPa
  observed.
- **Berman alternative:** the second record maps the paper's EosFit7c fit to
  `BM3 + ThermalReferenceStateEOS` with EosFit's truncated-quadratic Berman
  volume law. The published `alpha0 = 1.94(16)e-5 K^-1` and
  `dK0/dT = -0.008(3) GPa/K` are retained. Matching EosFit's effective-variance
  objective on the 41 printed heated rows gives `alpha0 = 1.811e-5 K^-1` and
  `dK0/dT = -0.01312 GPa/K`; Peritheos's latent-coordinate fit gives
  `1.781e-5 K^-1` and `-0.01232 GPa/K`. The remaining discrepancy may reflect
  a different internal data revision, unavailable unrounded inputs, or an
  undocumented fit selection. The direct EosFit7c result is available as the
  opt-in `b4c_somayazulu_2023_berman_refit` record; the published record is
  unchanged.
- **Primary data:** all 51 supplementary rows are embedded in `b4c.eosmat` as
  a typed, provenance-bearing dataset and linked to all three EOS records.
- **Source caveats:** the record explicitly preserves five contradictions:
  BM2 versus BM3 wording; 1425 versus 1450 K; the final article's
  `gamma0 = 0.8, q = 2.1` versus exploratory supplement pairs; and nominal
  300 K versus explicit 298 K reference temperature. In addition, the
  article's displayed exponential-integral equation is not EosFit's Berman
  equation despite the table and prose identifying the fit as Berman.
- **Crystal structure:** the material carries an average disordered
  `R-3m` (#166) hexagonal cell with five occupied asymmetric sites. Their
  multiplicities and split occupancies give B36C9 (`Z = 9` for B4C), and the
  measured `a = 5.601 A`, `c = 12.087 A` cell gives 328.383 A^3, matching the
  EOS reference volume. This is sufficient for PhaseSmith to calculate the
  reflection list and intensities.

### C02 — Cubic boron nitride: add the thermal branch

- **Change:** extend the existing `boron_nitride` (c-BN) material.
- **Source:** Datchi et al. (2007), [Equation of state of cubic boron nitride at
  high pressures and temperatures](https://doi.org/10.1103/PhysRevB.75.214104).
- **Published scope:** 295 K compression to 160 GPa and 500--900 K data to
  80 GPa. The paper gives a Vinet + Debye-Gruneisen parameter set, including
  `theta0 = 1700 K`, `gamma0 = 1.04(2)`, and `q = 4(1.5)`.
- **Likely mapping:** the existing 295 K `Vinet` record plus
  `MieGruneisenDebye` referenced to 295 K.
- **Effort/risk:** low-to-medium. The paper tabulates its Vinet component as a
  0 K/static curve, whereas Peritheos's generic MGD wrapper is referenced to a
  measured isotherm. Reproduce the paper numerically before deciding whether
  the existing 295 K curve can be composed directly or whether a cold-curve
  adapter is needed.

### C03 — Epsilon-FeOOH: replace/augment the static record with P-V-T

- **Change:** add a thermal record to the existing `e_feooh` material.
- **Source:** Suzuki (2016), [Pressure-volume-temperature equation of state of
  epsilon-FeOOH to 11 GPa and 700 K](https://doi.org/10.2465/jmps.160719c).
- **Published scope:** to 11 GPa and 700 K; `V0,300 = 66.278(6) A^3`,
  `K0 = 135(3) GPa`, `K0' = 6.1(9)`, `dK/dT = -0.05(2) GPa/K`, and
  `alpha(T) = 2.6(7)e-5 + 1.0(3)e-7 (T - 300) K^-1`.
- **Likely mapping:** `BM3` + `ThermalReferenceStateEOS` with integrated linear
  expansivity and `Tr = 300 K`.
- **Effort/risk:** low. Keep this composition/phase distinct from goethite and
  from high-pressure spin-transition descriptions of epsilon-FeOOH.

### C04 — CaSiO3 perovskite: add a lower-mantle thermal record

- **Change:** extend the existing `ca_perovskite` material.
- **Source:** Noguchi et al. (2013), [High-temperature compression experiments
  of CaSiO3 perovskite to lowermost mantle conditions and its thermal equation
  of state](https://doi.org/10.1007/s00269-012-0549-1).
- **Published scope:** to 127 GPa and 2300 K. The reference isotherm is at
  700 K (`BM2`, `V0 = 46.5(1) A^3`, `K0 = 207(4) GPa`) with an MGD correction
  (`theta0 = 1300(500) K`, `gamma0 = 2.7(3)`, `q = 1.2(8)`).
- **Likely mapping:** `BM2` + `MieGruneisenDebye`, explicitly using
  `Tr = 700 K` rather than silently converting to 300 K.
- **Effort/risk:** low. Validate the cell/formula-unit convention and do not
  merge parameters with the existing Shim or Mao 300 K records.

### C05 — Mg2SiO4 ringwoodite: add the thermal branch

- **Change:** extend the existing Mg end-member `ringwoodite` material.
- **Source:** Katsura et al. (2004), [Thermal expansion of Mg2SiO4 ringwoodite
  at high pressures](https://doi.org/10.1029/2004JB003094).
- **Published scope:** 0--24 GPa and 300--2000 K; `V0 = 524.8(1) A^3`, fixed
  `K0 = 182 GPa`, `K0' = 4.6(2)`, `theta0 = 846(26) K`, `gamma0 = 1.93(3)`,
  and `q = 3.5(3)`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`, `Tr = 300 K`.
- **Effort/risk:** low. This should be a new literature record, not a thermal
  component attached to the existing Meng et al. static fit.

### C06 — Bcc vanadium: practical high-temperature marker

- **Change:** add a new elemental material and P-V-T record.
- **Source:** Crichton et al. (2016), [High-temperature equation of state of
  vanadium](https://doi.org/10.1080/08957959.2015.1123256).
- **Published scope:** to 11.5 GPa and 1000 K; BM3 with
  `K0,300 = 150.4(6.2) GPa`, `K0' = 5.5(1.0)`,
  `alpha(T) = 4.8(6)e-5 - 2.4(9)e-8 T`, and
  `dK/dT = -0.0446(7) GPa/K`.
- **Likely mapping:** `BM3` + `ThermalReferenceStateEOS`.
- **Effort/risk:** low-to-medium. Recover the exact `V0` convention from the
  paper, express its alpha polynomial in Peritheos's reference-temperature
  convention, and encode the bcc-to-rhombohedral upper-pressure limit.

### C07 — B2 RbCl, RbBr, and RbI: a pressure-marker family

- **Change:** add three new alkali-halide materials as one coordinated review.
- **Source:** Farla et al. (2025), [Thermal equations of state of B2-structured
  rubidium halides RbCl, RbBr, and RbI](https://doi.org/10.1063/5.0248905).
- **Published scope:** to 26 GPa and 1800 K. All three use BM3 + MGD. Reported
  `(K0, K0', gamma0, q)` values are RbCl `(19.89(8), 5.00(2), 1.96(4),
  1.05(9))`, RbBr `(16.28(4), 5.28(2), 2.18(14), 1.52(24))`, and RbI
  `(13.69(4), 4.95(1), 2.21(7), 1.42(10))`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye` for each phase.
- **Effort/risk:** medium. Extract `V0`, `theta0`, atom counts, transition
  bounds, and structures from the full paper/supplement; keep pressure-scale
  provenance explicit because CsCl, Mo, and Pt were used in calibration.

## B — useful next

### C08 — Natural siderite Fe0.95Mn0.05CO3

- **Source:** Litasov et al. (2013), [P-V-T equation of state of siderite to
  33 GPa and 1673 K](https://doi.org/10.1016/j.pepi.2013.07.011).
- **Published scope:** BM3 with `V0 = 293.4(1) A^3`, `K0 = 120(1) GPa`,
  `K0' = 3.57(9)`, `dK/dT = -0.015(1) GPa/K`, and a reported linear
  ambient-pressure expansivity.
- **Likely mapping:** `BM3` + `ThermalReferenceStateEOS`.
- **Effort/risk:** medium. Name and formula must reflect the measured natural
  composition rather than claim a pure FeCO3 end-member; set the upper range
  below any spin transition or decomposition not represented by the fit.

### C09 — Superhydrous phase B, Mg10Si3O14(OH)4

- **Source:** Litasov et al. (2007), [Thermal equation of state of
  superhydrous phase B to 27 GPa and 1373 K](https://doi.org/10.1016/j.pepi.2007.06.003).
- **Published scope:** BM3 + MGD to 27 GPa and 1373 K. The preferred MGD table
  includes approximately `V0 = 623.53(35) A^3`, `K0 = 135.5(2.3) GPa`,
  `K0' = 5.3(2)`, fixed `theta0 = 860 K`, `gamma0 = 1.33(5)`, and
  `q = 2.03(35)`.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`.
- **Effort/risk:** medium. The paper provides fits on two Au pressure scales;
  select one record deliberately or preserve both as separately named scales.

### C10 — Fe2SiO4 ringwoodite

- **Source:** Armentrout and Kavner (2011), [High pressure, high temperature
  equation of state for Fe2SiO4 ringwoodite](https://doi.org/10.1029/2011GL046949).
- **Published scope:** MGD with `V0 = 42.03 cm3/mol`, `K0 = 202(4) GPa`, fixed
  `K0' = 4`, fixed `theta0 = 685 K`, `gamma0 = 1.08(6)`, and fixed `q = 2`.
- **Likely mapping:** `BM2` + `MieGruneisenDebye` as a new Fe end-member
  material, not as another record for Mg-ringwoodite.
- **Effort/risk:** medium. Convert molar to conventional-cell volume exactly
  and preserve that several thermal parameters were adopted/fixed.

### C11 — Dense hydrous phase E

- **Source:** Crichton and Ross (2000), [Equation of state of phase
  E](https://doi.org/10.1180/002646100549427).
- **Published scope:** single-crystal compression to 6.7 GPa for measured
  `Mg1.96(7)Fe0.072(5)Si1.04(5)H3.7(8)O6`, with BM3
  `K0 = 92.9(7) GPa` and `K0' = 7.3(2)`.
- **Likely mapping:** a new BM3 material record.
- **Effort/risk:** low-to-medium. Extract `V0` and atomic structure from the
  paper or a primary crystallographic source and retain the measured formula
  with its analytical uncertainties.

### C12 — Brucite Mg(OH)2

- **Sources:** Fei and Mao (1993), [Static compression of Mg(OH)2 to 78 GPa at
  high temperature](https://doi.org/10.1029/93JB00701), and Fukui et al.
  (2003), [Thermal expansion of Mg(OH)2 brucite under high
  pressure](https://doi.org/10.1007/s00269-003-0353-z).
- **Published scope:** the 1993 study reports 300 K `K0 = 54.3(1.5) GPa`,
  `K0' = 4.7(2)` and a separate 600 K fit; the 2003 work supplies a fuller
  P-V-T treatment to 16 GPa and 873 K.
- **Likely mapping:** start with a defensible 300 K BM3 record; review whether
  the later P-V-T formulation maps to `ThermalReferenceStateEOS`.
- **Effort/risk:** medium. Do not mix static and shock-derived parameters, and
  encode stability/dehydration limitations rather than extrapolating blindly.

### C13 — Delta-AlOOH

- **Source:** Vanpeteghem et al. (2002), [Equation of state of the hydrous phase
  delta-AlOOH at room temperature up to 22.5 GPa](https://doi.org/10.1029/2001GL014224).
- **Published scope:** two reported BM fits: `K0 = 252(3) GPa` with fixed
  `K0' = 4`, or `K0 = 228(3) GPa` with `K0' = 7(1)`.
- **Likely mapping:** paired BM2 and BM3 records, if the full data audit shows
  both are intended as publishable alternatives.
- **Effort/risk:** medium. Later work reports hydrogen-bond and spin-related
  changes; validity must stop before behavior not captured by a single BM fit.

### C14 — Fe-Si system: B20/B2 FeSi and Fe-9 wt% Si phases

- **Source:** Fischer et al. (2014), [Equations of state in the Fe-FeSi system
  at high pressures and temperatures](https://doi.org/10.1002/2013JB010898).
- **Published scope:** room-temperature and thermal BM3 + MGD fits for
  stoichiometric B20 and B2 FeSi and multiple Fe-9Si structures, with data to
  about 145 GPa and 3400 K for FeSi.
- **Likely mapping:** several phase- and composition-specific material records,
  each using `BM3` + `MieGruneisenDebye` where the paper supplies a complete
  thermal fit.
- **Effort/risk:** medium-to-high. This must be split by exact composition,
  structure, and phase field. Volumes are reported per mole of atoms in parts
  of the paper, so conversion errors are a major audit target.

### C15 — Ferropericlase Mg0.75Fe0.25O across spin crossover

- **Source:** Mao et al. (2011), [Thermal equation of state of lower-mantle
  ferropericlase across the spin crossover](https://doi.org/10.1029/2011GL049915).
- **Published scope:** to 140 GPa and 2000 K. Separate high-spin and low-spin
  reference fits are reported; at 300 K, the high-spin state uses
  `V0 = 76.34 A^3`, `K0 = 162(1) GPa`, fixed `K0' = 4`, while the low-spin
  fit uses `V0 = 74.44(60) A^3`, `K0 = 166(7) GPa`, fixed `K0' = 4`.
- **Likely mapping:** two bounded BM2 state records are possible immediately;
  the continuous crossover is not represented by existing Peritheos classes.
- **Effort/risk:** high scientific-review burden. Preserve composition and
  pressure-scale choice, and do not present a piecewise switch as the paper's
  continuous crossover model.

### C16 — Fe- and Al-bearing phase D P-V-T

- **Source:** Litasov et al. (2008), [Thermal equation of state of Al- and
  Fe-bearing phase D](https://doi.org/10.1029/2007JB004937).
- **Published scope:** measured
  `Mg0.99Fe0.12Al0.09Si1.75H2.51O6` to 20.6 GPa and 1273 K, analyzed with both
  high-temperature BM and MGD approaches.
- **Likely mapping:** new composition-specific `BM3` +
  `MieGruneisenDebye` record, rather than attaching thermal parameters to the
  existing antigorite-derived phase-D cards.
- **Effort/risk:** medium. Extract the complete preferred parameter table and
  distinguish scientific-fit provenance from the existing two phase-D
  compositions.

### C17 — Dolomite/ankerite family and dolomite-III

- **Sources:** Ross and Reeder (1992), [High-pressure structural study of
  dolomite and ankerite](https://msaweb.org/AmMin/AM77/AM77_412.pdf), and Mao
  et al. (2011), [Dolomite III: a new candidate lower mantle
  carbonate](https://doi.org/10.1029/2011GL049519).
- **Published scope:** the first paper gives low-pressure single-crystal P-V
  data for stoichiometric dolomite and about 70 mol% ankerite; the second gives
  a P-V EOS for Fe-bearing monoclinic dolomite-III observed at 36--83 GPa.
- **Likely mapping:** separate low-pressure dolomite, ferroan-dolomite/ankerite,
  and high-pressure dolomite-III records.
- **Effort/risk:** medium-to-high. Do not splice polymorphs or compositions.
  The old low-pressure paper has no DOI and may require refitting its published
  data; the high-pressure fit may involve a strongly extrapolated zero-pressure
  reference volume.

### C18 — Hcp zinc P-V-T

- **Source:** Errandonea et al. (2018), [High-pressure/high-temperature phase
  diagram of zinc](https://doi.org/10.1088/1361-648X/aacac0).
- **Published scope:** experimental XRD to 16 GPa and 1000 K, with wider
  theoretical coverage; the 300 K BM3 fit gives `V0 = 30.6(1) A^3`,
  `K0 = 63(2) GPa`, and `K0' = 5.6(4)`.
- **Likely mapping:** new hcp-Zn BM3 plus the paper's thermal reference-state
  formulation after equation-level audit.
- **Effort/risk:** medium. Record the strong axial anisotropy and proposed
  isomorphic anomaly near 10 GPa in provenance/validity notes even though the
  scalar volume EOS remains continuous.

### C19 — Molybdenum: independent experimental MGD scale

- **Change:** add an alternate record to existing `molybdenum`.
- **Source:** Huang et al. (2016), [Thermal equation of state of molybdenum
  determined from in situ synchrotron X-ray diffraction](https://doi.org/10.1038/srep19923).
- **Published scope:** 300 K compression to 80 GPa and thermal data to 92 GPa
  and 3470 K. The 300 K BM3 fit gives `V0 = 31.22(8) A^3`,
  `K0 = 273(15) GPa`, and `K0' = 3.6(4)`; the MGD fit fixes the Debye
  temperature near 470 K and fits `q = 0.6`, alongside a separate
  thermodynamic thermal-pressure fit.
- **Likely mapping:** `BM3` + `MieGruneisenDebye`, as an explicitly named
  alternative to the current Sokolova composite scale.
- **Effort/risk:** medium. Resolve which of the paper's two formulations is the
  primary executable record and preserve the fixed-versus-fitted status of the
  Debye temperature.

## C — model work

### C20 — (Al,Fe)-phase H with a spin crossover

- **Source:** Strozewski et al. (2023), [Equation of state and spin crossover
  of (Al,Fe)-phase H](https://doi.org/10.1029/2022JB026291).
- **Published scope:** `Al0.84Fe3+0.07Mg0.02Si0.06OOH` to 125 GPa. The paper
  fits high- and low-spin elastic states simultaneously and places the 300 K
  crossover midpoint at `55.5(9) GPa`, spanning roughly 48--63 GPa.
- **Likely mapping:** a new spin-crossover free-energy EOS, not an ordinary BM3
  record. The paper's MINUTI formulation adds a spin partition-function term
  to the elastic free energy.
- **Effort/risk:** high. This is the clearest candidate for a new general model
  class, but implementation should wait until its thermodynamic variables,
  parameter identifiability, inversion behavior, and serialization contract
  have been designed.

## Proposed one-by-one review order

1. C01 B4C thermal extension
2. C03 epsilon-FeOOH thermal extension
3. C05 Mg-ringwoodite thermal extension
4. C04 CaSiO3-perovskite thermal extension
5. C06 bcc vanadium
6. C07 B2 rubidium halides
7. C02 c-BN thermal extension
8. C08 siderite
9. C09 superhydrous phase B
10. C10 Fe-ringwoodite
11. C11 phase E
12. C12 brucite
13. C16 Fe-Al phase D
14. C14 Fe-Si phase suite
15. C18 zinc
16. C19 molybdenum alternate scale
17. C13 delta-AlOOH
18. C17 dolomite/ankerite polymorph suite
19. C15 ferropericlase spin states
20. C20 spin-crossover phase H

For each review, the acceptance gate should be: obtain the primary full text
and supplement; identify exact composition, phase, cell and volume convention;
transcribe the published equation rather than only its parameter table;
reproduce at least one published P-V(-T) point; record fixed/fitted parameters
and uncertainty meaning; define the experimentally supported validity domain;
and only then create or modify an `.eosmat` record.
