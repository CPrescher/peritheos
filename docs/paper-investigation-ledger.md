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
  similarity criteria, or a coupled source-level objective had a demonstrably
  different optimum. These are source-fit discrepancies, not software-run
  failures; the record may remain for faithful published-curve provenance.
- **Direct refit unavailable:** the equation and parameters were audited, but
  independent coefficient recovery was impossible because primary rows, an
  executable calibration, or the original reduction were unavailable or
  circular.
- **Withheld/deferred:** investigation did not pass the executable-record
  acceptance gate, so no production EOS was added.

## Summary

The register covers **318 primary papers**: **220** support the 814 audited catalog records and **98** were investigated without adding a production record.
No numerical refit attempt failed before producing a comparison. The adverse
outcomes are instead explicit coefficient discrepancies, unavailable direct
refits, or acceptance-gate holds.

| Paper-level outcome | Papers |
|---|---:|
| Reproduced | 126 |
| Partly reproduced | 8 |
| Mixed: reproduced and discrepant records | 4 |
| Coefficient parity not achieved | 12 |
| Direct refit unavailable | 74 |
| Withheld: could not reproduce | 6 |
| Deferred: incomplete source/model mapping | 87 |

## Withheld or deferred papers

### [A modified Anderson–Grüneisen model for the pressure dependence of thermal expansivity (2019)](https://doi.org/10.1139/cjp-2019-0326)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Ahrens and Gaffney (1975), Mg-Fe silicate shock interpretation](https://doi.org/10.1111/j.1365-246X.1975.tb06463.x)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The legacy compressibility and modulus values lack reference volumes and mostly omit pressure derivatives; the pure-Mg endpoint is extrapolated or cited rather than an independent complete fit.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Akaogi et al. (1996), natural-basalt majoritic garnet](https://doi.org/10.1029/96GL03027)

**Outcome:** Direct refit unavailable (2026-09-06).

The two candidate fits require the complex sample's reference volume and exact fit convention, but the final primary article was not recoverable; two further rows are citations.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [ANALYSIS OF THERMAL EXPANSIVITY OF SOLIDS UNDER HIGH PRESSURES (2012)](https://doi.org/10.1142/s0217984912501461)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Anderson and Zou (1990), MgO thermodynamic functions](https://doi.org/10.1063/1.555873)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The temperature-dependent thermodynamic construction is not a standalone complete BM record; the source row omits V0, K0 prime, equation order, and the full thermal parameterization.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Braithwaite (2002), high-pressure solid/melt representations](https://doi.org/10.1063/1.1483512)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The primary Vinet/thermal table omits V0. Its liquid-MgO row is therefore non-executable, while the solid row contains an obvious table-column extraction error.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Caracas and Cohen (2005), MgSiO3-FeSiO3-Al2O3 pv/ppv chemistry](https://doi.org/10.1029/2005GL023164)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

All 16 LitCurate rows are held because the primary parameter table was unavailable and the accessible abstract does not establish the EOS equation, volume basis, or magnetic/electronic branch identities. No coefficients were promoted from discovery-only evidence.

Evidence: [literature-reproductions/caracas-cohen-2005-chemistry.md](literature-reproductions/caracas-cohen-2005-chemistry.md).

### [Chen et al. (2024), stishovite velocities](https://doi.org/10.1029/2023GL107700)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The new results are acoustic velocities and modulus derivatives rather than a pressure-volume fit, and the extracted rows omit V0.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Comparative Analysis of Grüneisen Parameters for Selected Geophysical Minerals Using Advanced Equations of State (2024)](https://doi.org/10.69626/sea.2024.0152)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All three MgO rows compare alternative laws using the same incomplete adopted input triplet and omit V0; the paper supplies no new independently fitted EOS.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Compositional constraints on the equation of state and thermal properties of the lower mantle (2001)](https://doi.org/10.1046/j.1365-246x.2001.00437.x)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 17 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Constraints on lower mantle composition and temperature from density and bulk sound velocity profiles (1990)](https://doi.org/10.1029/gl017i008p01153)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 6 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Dewaele et al. (2005), high-pressure metrology abstract](https://doi.org/10.1107/S0108767305096972)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

Both LitCurate rows are publication-level misattributions: this DOI is a pressure-metrology abstract and does not report either claimed silicate EOS.

Evidence: [literature-reproductions/iucr-2005-s0108767305096972.md](literature-reproductions/iucr-2005-s0108767305096972.md).

### [Effect of Pressure on the Composition of the Lower Mantle End Member Fe x O (1993)](https://doi.org/10.1126/science.259.5091.66)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 5 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Elastic properties of Fe-bearing Akimotoite at mantle conditions: Implications for composition and temperature in lower mantle transition zone (2022)](https://doi.org/10.1016/j.fmre.2021.12.013)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 3 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Equation of State, Phase Stability of (Mg0.92, Fe0.08)SiO3 Perovskite from Shock Wave Study and Its Geophysical Implications (2004)](https://doi.org/10.1063/1.1780510)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The source-labelled proceedings row is a publication duplicate of the data-backed Gong et al. shock EOS already accepted under the primary GRL DOI 10.1029/2003GL019132; it is not an independent fit.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Exploring the High-Pressure Equation of State in Earth’s Mantle with a Focus on the MgSiO3−MgO System (2025)](https://doi.org/10.15407/mfint.47.06.0601)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [EXTREME COMPRESSION BEHAVIOUR OF SOLIDS BASED ON THE ROY-ROY INVERTED EQUATION OF STATE (2008)](https://doi.org/10.1142/s0217979208038910)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Faust and Knittle (1994), natural chondrodite](https://doi.org/10.1029/94GL01592)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The real BM3 fit belongs to a natural F-rich sample whose exact reusable composition and site occupancies could not be established.

Evidence: [literature-reproductions/tranche-b-mineral-eos-audit.md](literature-reproductions/tranche-b-mineral-eos-audit.md).

### [Ferre et al. (2009), dislocations in CaSiO3 perovskite](https://doi.org/10.2138/am.2009.3003)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The dislocation and Peierls-Nabarro study's source and comparison pairs omit V0 and do not define an executable pressure-volume fit.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Finding the isentropic density of perovskite: Implications for iron concentration in the lower mantle (1997)](https://doi.org/10.1029/96gl03951)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 5 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Fuchizaki (2019), updated MgO melting curve](https://doi.org/10.7566/JPSJ.88.065003)

**Outcome:** Direct refit unavailable (2026-09-06).

The source-labelled row has no V0 and an anomalous fixed derivative; the access-controlled primary equation/table could not be inspected to repair the suspected extraction error.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Funamori et al. (1996), MgSiO3 perovskite thermoelasticity](https://doi.org/10.1029/95JB03732)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The two apparent LitCurate records are run-specific ambient-volume normalizations of a single combined thermal EOS, not two independently fitted equations. Both adopt K0=261 GPa and K0-prime=4 from Mao et al. (1991), while the source's executable thermal parameterization is not represented by either split candidate. No production EOS was added.

Evidence: [literature-reproductions/funamori-1996-mgsio3-perovskite.md](literature-reproductions/funamori-1996-mgsio3-perovskite.md).

### [Fundamental thermodynamic relations and silicate melting with implications for the constitution of D″ (1990)](https://doi.org/10.1029/jb095ib12p19311)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Haavik et al. (2000), defect-clustered wuestite](https://doi.org/10.1039/B006026G)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The accessible primary evidence does not supply a complete executable coefficient set and reference volume for any source-owned branch.

Evidence: [literature-reproductions/haavik-2000-wustite.md](literature-reproductions/haavik-2000-wustite.md).

### [Hama and Suito (1998), Mg-Fe bridgmanite](https://doi.org/10.1029/97JB03672)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The closed primary and incomplete thermal/model lineage do not support the apparent independent source branches.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Hama and Suito (2000), Vinet-Debye mineral parameterizations](https://doi.org/10.2138/am-2000-2-310)

**Outcome:** Direct refit unavailable (2026-09-06).

The final article was not recoverable, so exact same-DOI coefficients, compositions, and reference states for the three apparent Vinet inputs could not be independently established.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [High pressure and high temperature in situ X‐ray observation of MgSiO3 Perovskite under lower mantle conditions (1993)](https://doi.org/10.1029/92gl02960)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Inbar and Cohen (1995), MgO thermal properties](https://doi.org/10.1029/95GL01086)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source reports thermodynamic derivatives at four temperatures but no analytical pressure-volume EOS family or numerical P-V grid; BM3 would be an unsupported assignment.

Evidence: [literature-reproductions/inbar-cohen-1995-mgo.md](literature-reproductions/inbar-cohen-1995-mgo.md).

### [Indoor seismology by probing the Earth's interior by using sound velocity measurements at high pressures and temperatures (2007)](https://doi.org/10.1073/pnas.0608609104)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 5 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Isaak et al. (1990), calculated MgO properties](https://doi.org/10.1029/JB095iB05p07055)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

All fourth-order discovery rows omit K0-double-prime and the third-order sensitivity omits V0; the inaccessible primary could not repair them.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Ito et al. (2011), Kawai apparatus](https://doi.org/10.4131/jshpreview.21.272)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The apparatus paper quotes two alternative modulus/derivative pairs but no V0 and does not fit an independent material EOS.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Jackson et al. (1990), wuestite elasticity](https://doi.org/10.1029/JB095iB13p21671)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source rows are acoustic extrapolations that omit V0, not three independent pressure-volume fits.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Jackson et al. (2006), ferropericlase elasticity](https://doi.org/10.1029/2005JB004052)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source-owned quantities are acoustic moduli and derivatives; the apparent EOS rows do not define independent pressure-volume fits.

Evidence: [literature-reproductions/jackson-2006-ferropericlase-elasticity.md](literature-reproductions/jackson-2006-ferropericlase-elasticity.md).

### [Jacobs and Oonk (2000), GGK MgO EOS](https://doi.org/10.1039/A910247G)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source uses an unsupported GGK pressure law; its ambient derivative triple cannot be re-labelled as BM3.

Evidence: [literature-reproductions/jacobs-oonk-2000-mgo.md](literature-reproductions/jacobs-oonk-2000-mgo.md).

### [Karki and Stixrude (1999), lower-mantle elastic moduli](https://doi.org/10.1029/1999JB900069)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The five source rows are finite-strain fits of elastic moduli M(P), not pressure-volume equations of state, and they contain no V0. Encoding them as volumetric BM records would be scientifically false, so all five were rejected from production.

Evidence: [literature-reproductions/oganov-wang-chizmeshya-karki-eos-audit.md#karki-and-stixrude-1999-0-accepted-5-rejected](literature-reproductions/oganov-wang-chizmeshya-karki-eos-audit.md#karki-and-stixrude-1999-0-accepted-5-rejected).

### [Karki et al. (2000), theoretical MgSiO3 akimotoite](https://doi.org/10.2138/am-2000-2-309)

**Outcome:** Direct refit unavailable (2026-09-06).

The complete-looking static BM3 row could not be checked because the final primary article was unavailable; the remaining three rows are earlier experimental citations.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Karki et al. (2001), MgSiO3 thermodynamic derivatives](https://doi.org/10.1029/2001GL012910)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source prints derivative triples from quasiharmonic free energies but does not specify an executable analytical volume-EOS family.

Evidence: [literature-reproductions/tranche-b-mineral-eos-audit.md](literature-reproductions/tranche-b-mineral-eos-audit.md).

### [Katsura et al. (2004), Mg2SiO4 ringwoodite](https://doi.org/10.1029/2004JB003094)

**Outcome:** Withheld: could not reproduce (2026-09-05).

All 127 official Table 2 observations were transcribed, but the published BM3-MGD coefficients miss them by 1.890 GPa RMSE when the chemically required n=7 atoms per formula unit is used. A refit moves gamma0 far outside its reported uncertainty, while an unphysical and unpublished n=5 normalization fits much better. The hidden normalization, Debye-temperature law, and energy/volume convention remain unresolved, so no production EOS was added.

Evidence: [literature-reproductions.md#ringwoodite-katsura-2004](literature-reproductions.md#ringwoodite-katsura-2004).

### [Kawai and Tsuchiya (2015), CaSiO3 thermoelasticity](https://doi.org/10.1002/2015GL063446)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The candidate is a citation-reported earlier thermoelastic parameterization already represented from primary DOI 10.1002/2013JB010905; this extracting paper adds no independent EOS.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Kholiya et al. (2014), MgO EOS comparison](https://doi.org/10.1155/2014/289353)

**Outcome:** Withheld: could not reproduce (2026-09-06).

All eight rows quote earlier publications; the source compares normalized-volume predictions and reports no absolute V0 or new fit.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Knittle and Jeanloz (1991), MgSiO3 perovskite transition and thermal expansion](https://doi.org/10.1126/science.251.4992.410)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The paper reports a phase transition and thermal expansion; its assumed modulus and derivative do not form a complete source-owned EOS with a stated equation.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Kumari and Dass (1990), EOS applied to 50 solids II](https://doi.org/10.1088/0953-8984/2/39/003)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The MgO source row omits the absolute V0 and does not resolve the equation identity; its second row is citation-reported.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Lakshtanov et al. (2007), Al-H stishovite elasticity](https://doi.org/10.2138/am.2007.2294)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source's composition-specific elastic results do not include a complete source-defined volume EOS model.

Evidence: [literature-reproductions/lakshtanov-2007-alh-stishovite-elasticity.md](literature-reproductions/lakshtanov-2007-alh-stishovite-elasticity.md).

### [Li-and-Zhang-data elastic refit (2010)](https://doi.org/10.1016/j.pnsc.2009.09.002)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The source row re-fits earlier elastic data without a new pressure-volume dataset, exact Fe content, or reference volume; the other two rows are citations.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Liu (2008), MgO bulk modulus method](https://doi.org/10.1515/zna-2008-1-209)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The method paper adopts its ambient modulus from a cited result and does not provide V0, a derivative, or a complete source-owned executable EOS.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Liu (2011), two-parameter MgO EOS analysis](https://doi.org/10.1139/p11-040)

**Outcome:** Withheld: could not reproduce (2026-09-06).

All five temperature-dependent modulus pairs are citation-reported adopted inputs and omit V0; the paper generates no independent fit.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Lower-mantle MgSiO3 thermoelastic parameter review (1996)](https://doi.org/10.1098/rsta.1996.0053)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The apparent bridgmanite parameters cannot be separated into a complete pure-endmember source fit from accessible primary evidence; the ferropericlase row is explicitly cited.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Matsui (1993), molecular dynamics of high-pressure silicates](https://doi.org/10.5940/jcrsj.35.190)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The source molecular-dynamics row reports V0 and K0 but neither an analytical EOS family nor K0-prime; the experimental row is a cited comparison.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Matsui (2002), MgSiO3-Al2O3 simulations](https://doi.org/10.2465/jmps.97.13)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source calls the fits BM4 but reports only bulk moduli; V0, K0-prime, K0-double-prime, and fit observations are absent.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Matsui et al. (1994), MgSiO3 simulations](https://doi.org/10.1029/94GL01370)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The primary reports isolated volumes and moduli but no K0-prime or numerical pressure-volume grid for its stated BM3 reductions.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Matsui et al. (2000), MgO pressure standard](https://doi.org/10.2138/am-2000-2-308)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The complete-looking row remains held because the primary full text was unavailable and equation, reference state, and fit basis could not be audited.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [McCarthy and Harrison (1994), MgO bulk properties](https://doi.org/10.1103/PhysRevB.49.8574)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The reported values are zero-pressure elastic bulk moduli from local distortion fits, not pressure-volume EOS parameterizations.

Evidence: [literature-reproductions/mccarthy-harrison-1994-mgo.md](literature-reproductions/mccarthy-harrison-1994-mgo.md).

### [MgSiO3 elasticity temperature-derivative study (2004)](https://doi.org/10.1029/2003GL018762)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All five candidate rows are comparison or citation values and none supplies a source-owned reference volume and executable pressure-volume equation.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [MgSiO3 interaction-potential study (2007)](https://doi.org/10.1088/1674-0068/20/05/547-551)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The source reports only an equilibrium volume for a new interaction potential, not a complete EOS parameterization.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [MgSiO3 molecular-dynamics compression study (2006)](https://doi.org/10.1360/CJCP2006.19(4).311.4)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The article reports an ambient molar volume and plots pressure-volume behavior but gives no complete analytical EOS and coefficient set; its derivative-only row is cited.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Miyajima et al. (2025), electron diffraction of a dense hydrous magnesium silicate](https://doi.org/10.1029/2025GL115280)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The paper reports a single-state electron-diffraction structure and ambient volume plus literature comparisons, not a fitted EOS.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Mookherjee et al. (2019), phase Egg accepted-manuscript alias](https://doi.org/10.2138/am-2018-6694)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

This noncanonical accepted-manuscript alias directs citation to final DOI 10.2138/am-2019-6694, whose two records are already bundled; no duplicate was created.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Morishima et al. (1994), fixed-pressure CaSiO3 thermal expansion](https://doi.org/10.1029/94GL00844)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The experiment measures thermal expansion at fixed pressure; the candidate compression coefficients are assumed inputs and V0 is absent.

Evidence: [literature-reproductions/ca-other-source-exhaustion-2026-09-06.md](literature-reproductions/ca-other-source-exhaustion-2026-09-06.md).

### [Myhill (2022), anisotropic high-P-T EOS](https://doi.org/10.1093/gji/ggac180)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source's coupled anisotropic SLB2005 tensor/thermal formulation is not representable by a scalar Peritheos volume EOS without changing its meaning.

Evidence: [literature-reproductions/myhill-2022-anisotropic-eos.md](literature-reproductions/myhill-2022-anisotropic-eos.md).

### [Pamato et al. (2016), NAL elasticity](https://doi.org/10.1002/2016JB013136)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

Both rows are elastic finite-strain summaries missing the pressure derivatives required for the claimed BM4 mapping.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Performance of an ab initio equation of state for magnesium oxide (2004)](https://doi.org/10.1088/0953-8984/16/30/006)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Pressure Dependence of Interatomic Separation and Thermal Expansivity for Alkali Halides and Periclase (MgO) (2009)](https://doi.org/10.12693/aphyspola.115.709)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Pressure-Induced Magnetization in FeO: Evidence from Elasticity and Mössbauer Spectroscopy (2004)](https://doi.org/10.1103/physrevlett.93.215502)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Pressure-volume equation of state of the high-pressureB2phase of NaCl (2002)](https://doi.org/10.1103/physrevb.65.104114)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Refinement of enthalpy measurement of MgSiO3 perovskite and negative pressure‐temperature slopes for Perovskite‐forming reactions (1993)](https://doi.org/10.1029/93gl01265)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Reichmann et al. (2008), MgO elasticity](https://doi.org/10.2138/am.2008.2717)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The sole source row is an adiabatic acoustic bulk-modulus fit without compression V0, not a volumetric EOS.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Research progresses on the structure, phase transition and physical properties of MgSiO&amp;lt;sub&amp;gt;3&amp;lt;/sub&amp;gt; under high temperature and high pressure (2025)](https://doi.org/10.3724/j.issn.1007-2802.20240151)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 22 rows are review-level summaries of earlier MgSiO3 and Fe/Al-bearing primary studies, so their source lineage remains with those cited papers.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Satta et al. (2025), Fe-bearing delta-AlOOH elasticity](https://doi.org/10.1007/s00269-025-01319-7)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source reports elastic-modulus finite-strain fits, not standalone volumetric P(V) equations with a complete reference state.

Evidence: [literature-reproductions/satta-2025-delta-alooh-elasticity.md](literature-reproductions/satta-2025-delta-alooh-elasticity.md).

### [Sherman et al. (1993), stishovite and modified-fluorite SiO2](https://doi.org/10.1029/93JB00783)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The accessible primary abstract reports a static HF-LCAO stishovite BM3 parameterization, but the calculated energy-volume grid is confined to an inaccessible figure and the abstract conflicts with body values extracted by LitCurate. The modified-fluorite candidate also has an unresolved conventional-cell basis. Both source fits are therefore documented but withheld, and ten comparison values are retained only as citation traces.

Evidence: [literature-reproductions/sherman-1993-stishovite.md](literature-reproductions/sherman-1993-stishovite.md).

### [Shieh et al. (2002), MgSiO3 post-perovskite](https://doi.org/10.1103/PhysRevLett.89.255507)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source gives fixed K0 and fitted K0-prime but omits V0, preventing an executable EOS without importing a value from elsewhere.

Evidence: [literature-reproductions/tranche-b-mineral-eos-audit.md](literature-reproductions/tranche-b-mineral-eos-audit.md).

### [Shukla et al. (2016), Fe3+- and Al-bearing bridgmanite](https://doi.org/10.1002/2016GL069332)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

Three composition-specific calculated rows are held because the primary article prints V0, K0, and K0-prime but never identifies the analytical cold-compression equation. Three further rows are comparisons owned by cited primary papers.

Evidence: [literature-reproductions/shukla-2016-ferric-al-bridgmanite.md](literature-reproductions/shukla-2016-ferric-al-bridgmanite.md).

### [Singh and Singh (2021), alkaline-earth oxide EOS formulation](https://doi.org/10.12693/aphyspola.140.131)

**Outcome:** Withheld: could not reproduce (2026-09-06).

The MgO and CaO modulus pairs are explicitly cited inputs; plots use only V/V0, and no new absolute-volume fit is reported.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Sinogeikin and Bass (1999), MgO elasticity](https://doi.org/10.1103/PhysRevB.59.R14141)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The complete fourth-order elastic finite-strain coefficients omit V0; the reference volume cannot be reconstructed without invention.

Evidence: [literature-reproductions/tranche-b-mineral-eos-audit.md](literature-reproductions/tranche-b-mineral-eos-audit.md).

### [Sinogeikin et al. (2004), MgSiO3 elasticity](https://doi.org/10.1029/2004GL019559)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The paper measures elastic properties and does not publish the complete independent volumetric EOS implied by the discovery rows.

Evidence: [literature-reproductions/sinogeikin-2004-mgsio3-elasticity.md](literature-reproductions/sinogeikin-2004-mgsio3-elasticity.md).

### [Sokolova et al. (2018), MgO-MgSiO3 spreadsheets](https://doi.org/10.1080/08957959.2018.1465056)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The spreadsheet parameterizations use the same unsupported coupled Kunc/two-Einstein/Altshuler model and cannot be reduced to BM triplets.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Sokolova et al. (2021), Ca-silicate EOS](https://doi.org/10.3390/min11030322)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The coupled Kunc/two-Einstein/Altshuler Helmholtz model is not currently supported as one faithful Peritheos EOS.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Sound-velocity temperature coefficients of MgSiO3 perovskite (2000)](https://doi.org/10.1088/0256-307X/17/3/022)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The sole 263 GPa value is citation-reported and has neither a reference volume nor pressure derivative, so it is not an executable source EOS.

Evidence: [literature-reproductions/mgsio3-source-exhaustion.md](literature-reproductions/mgsio3-source-exhaustion.md).

### [Spin crossover and Mott—Hubbard transition under high pressure and high temperature in the low mantle of the Earth (2015)](https://doi.org/10.1088/1742-6596/653/1/012095)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Stixrude et al. (1992), thermoelasticity and mantle stratification](https://doi.org/10.1126/science.257.5073.1099)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The Fe-bearing assemblage rows do not establish exact Birch-Murnaghan order and reference state; one is additionally missing both modulus coefficients.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Taniguchi et al. (1995), Ca-silicate calculation models](https://doi.org/10.2465/minerj.17.290)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

A calculated V0-K0-K0-prime triple is printed without identifying the analytical EOS or reference temperature, so no model is inferred.

Evidence: [literature-reproductions/taniguchi-1995-casilicates.md](literature-reproductions/taniguchi-1995-casilicates.md).

### [The effect of temperature on the product of bulk modulus and volume thermal expansion coefficient, and its application to the thermal expansion of MgO and other minerals (2004)](https://doi.org/10.1002/pssb.200302047)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [The sound velocity of wüstite at high pressures: implications for low-velocity anomalies at the base of the lower mantle (2020)](https://doi.org/10.1186/s40645-020-00333-3)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [The texture of the post-perovskite phase controls the characteristics of the D” seismic discontinuity (2025)](https://doi.org/10.1038/s43247-025-02383-1)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 2 LitCurate candidates are explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Thermal pressure in MgO and MgSiO 3 perovskite at lower mantle conditions (2000)](https://doi.org/10.2138/am-2000-1013)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Thermo-chemical and thermo-physical properties of the high-pressure phase anhydrous B (Mg14Si5O24): An ab-initio all-electron investigation (2010)](https://doi.org/10.2138/am.2010.3368)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Thermodynamic estimation the compressibility of ferropericlase under high pressure (2016)](https://doi.org/10.1063/1.4967779)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 31 rows are a compilation of earlier MgO, FeO, and ferropericlase moduli across composition rather than source-owned complete EOS fits.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Thermoelasticity of perovskite: An emerging consensus (1994)](https://doi.org/10.1029/94eo01093)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

All 1 LitCurate candidate is explicitly citation-reported or adopted inputs rather than a source-owned EOS fit under this DOI.

Evidence: [literature-reproductions/litcurate-source-exhaustion-citation-audit.md](literature-reproductions/litcurate-source-exhaustion-citation-audit.md).

### [Tsuchiya and Kawamura (2001), B1 oxide elasticity](https://doi.org/10.1063/1.1371498)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The source MgO BM3 row omits V0 and primary observations; the other five entries are citation comparisons.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Tsuchiya et al. (2004), MgSiO3 phase transition](https://doi.org/10.1016/j.epsl.2004.05.017)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

Table 1 merges the uncertainty ranges of BM3 and Vinet fits without publishing either model's individual coefficients or the calculated P-V rows.

Evidence: [literature-reproductions/tsuchiya-2004-mgsio3-phase-transition.md](literature-reproductions/tsuchiya-2004-mgsio3-phase-transition.md).

### [Tsuchiya et al. (2004), post-perovskite elasticity](https://doi.org/10.1029/2004GL020278)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The paper reports elastic tensors and quotes EOS summaries from a separate study; it does not own a complete executable P(V) fit.

Evidence: [literature-reproductions/tsuchiya-2004-post-perovskite-elasticity.md](literature-reproductions/tsuchiya-2004-post-perovskite-elasticity.md).

### [Vijay (2024), generalized Rydberg-Vinet and Stacey thermoelasticity](https://doi.org/10.32908/hthp.v53.1503)

**Outcome:** Withheld: could not reproduce (2026-09-06).

Both MgO rows reuse the same citation-reported modulus pair under comparison equations and omit V0 and reference temperature.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Wang et al. (2026), KAlSi3O8 liebermannite and K-hollandite II](https://doi.org/10.2138/am-2024-9562)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-04).

The deposited temperature-indexed BM3 curves can be reproduced, but they are derived grids rather than primary fit observations. Their V0(T), K0(T), and K0-prime(T) behavior is not faithfully represented by the current thermal wrapper, and the primary methods, uncertainty/covariance information, exact reference-state convention, and authoritative symmetry-reduced structures were unavailable. The candidate remains deliberately non-executable.

Evidence: [material-eos-candidates.md](material-eos-candidates.md).

### [Wentzcovitch et al. (1993), MgSiO3 molecular dynamics](https://doi.org/10.1103/PhysRevLett.70.3947)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The source reports modulus/derivative pairs but no equilibrium volume or numerical pressure-volume grid for an executable EOS.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

### [Wu et al. (2016), Fe-Al phase D elasticity](https://doi.org/10.1002/2016JB013209)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

The extracted rows are acoustic finite-strain elasticity results rather than complete pressure-volume EOS parameterizations.

Evidence: [literature-reproductions/wu-2016-feal-phase-d.md](literature-reproductions/wu-2016-feal-phase-d.md).

### [Xu et al. (2024), Al-bearing superhydrous phase B](https://doi.org/10.1029/2023GL107818)

**Outcome:** Withheld: could not reproduce (2026-09-05).

The official EarthChem archive contains all 45 P-T-V rows, but none of the four published Birch-Murnaghan parameterizations reproduces those observations: direct refits shift K0 by roughly 18-32 GPa and K0-prime substantially. Three additional rows are coupled acoustic-elasticity fits rather than standalone P(V) equations, and three are citation-only comparisons. No production EOS was added.

Evidence: [literature-reproductions/xu-2024-al-bearing-superhydrous-phase-b.md](literature-reproductions/xu-2024-al-bearing-superhydrous-phase-b.md).

### [Zha, Mao, and Hemley (2000), MgO elasticity pressure scale](https://doi.org/10.1073/pnas.240466697)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-06).

The primary source reports adiabatic and corrected isothermal modulus pairs but no absolute V0 in the EOS parameter table; the discovery row also mislabels the adiabatic fit order.

Evidence: [literature-reproductions/mgo-oxide-litcurate-exhaustion.md](literature-reproductions/mgo-oxide-litcurate-exhaustion.md).

### [Zhang and Weidner (1999), Al-enriched silicate perovskite](https://doi.org/10.1126/science.284.5415.782)

**Outcome:** Direct refit unavailable (2026-09-05).

The three apparent records describe separate runs on the same 5 mol% Al2O3 sample, correctly normalized as Mg0.95Al0.10Si0.95O3, but the closed primary article exposes neither a numerical P-V-T table nor a recoverable pressure-volume figure. The extracted coefficients have no source uncertainties or independent observations with which to verify the run fits, so all three remain on hold and five comparison rows are citation traces.

Evidence: [literature-reproductions/zhang-weidner-1999-aluminous-perovskite.md](literature-reproductions/zhang-weidner-1999-aluminous-perovskite.md).

### [Zhang et al. (2009), Mg-Fe silicate phase stability](https://doi.org/10.1142/S0217979209053047)

**Outcome:** Deferred: incomplete source/model mapping (2026-09-05).

All purported BM4 branches omit K0-double-prime and the inaccessible primary could not resolve collapsed phase identities.

Evidence: [literature-reproductions/tranche-c-zero-yield-audits.md](literature-reproductions/tranche-c-zero-yield-audits.md).

## Papers with coefficient discrepancies

These **16 papers** account for all 31 records
classified as `parity_not_achieved`. Papers with other successful records
are marked as mixed in the complete register.

| Paper | Affected record | Published-to-refit discrepancy |
|---|---|---|
| [Anzellini et al. (2019)](https://doi.org/10.1038/s41598-019-51931-1) | `silicon_vii_anzellini_2019_vinet_1` | K0 96.9 -> 4.845 |
| [Baty et al. (2024)](https://doi.org/10.1063/5.0179469) | `palladium_baty_2024_bm3_1` | K0 190 -> 152.057 |
| [Clendenen and Drickamer (1966)](https://doi.org/10.1063/1.1726610) | `coo_clendenen_1966_murnaghan_1` | K0_prime 3.9 -> 5.1481 |
| [Dewaele (2019)](https://doi.org/10.3390/min9110684) | `iron_dewaele_2019_dor_vinet` | V0 22.354 -> 23.5216; K0 168.4 -> 101.058; K0_prime 5.33 -> 6.63486 |
|  | `iron_dewaele_2019_mao_vinet` | K0 164.5 -> 101.058; K0_prime 4.96 -> 6.63486 |
|  | `lead_hcp_dewaele_2019_dor_vinet` | K0_prime 4.77 -> 3.30377 |
|  | `lead_hcp_dewaele_2019_mao_vinet` | K0_prime 4.4 -> 3.30377 |
| [Dobrosavljevic et al. (2019)](https://doi.org/10.3390/min9120762) | `mgfe94o_b1_dobrosavljevic_2019_bm3_1` | K0_prime 3.79 -> 2.71444 |
|  | `mgfe94o_rhombohedral_dobrosavljevic_2019_bm3_1` | K0 217 -> 168.771 |
| [Dorfman et al. (2012)](https://doi.org/10.1029/2012jb009292) | `gold_dorfman_2012_tange_mgo_k0_fixed_vinet` | K0_prime 5.88 -> 5.86027 |
|  | `gold_dorfman_2012_tange_mgo_k0_free_vinet` | K0 167 -> 178.713; K0_prime 5.84 -> 5.46833 |
|  | `molybdenum_dorfman_2012_tange_mgo_k0_fixed_vinet` | K0_prime 4.19 -> 4.1381 |
|  | `molybdenum_dorfman_2012_tange_mgo_k0_free_vinet` | K0 271 -> 255.46; K0_prime 3.89 -> 4.2951 |
|  | `platinum_dorfman_2012_tange_mgo_k0_fixed_vinet` | K0_prime 5.43 -> 5.36541 |
|  | `platinum_dorfman_2012_tange_mgo_k0_free_vinet` | K0 280 -> 294.309; K0_prime 5.29 -> 4.89884 |
| [Finkelstein et al. (2017)](https://doi.org/10.2138/am-2017-5966) | `mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1` | K0 148 -> 337.791; K0_prime 4.09 -> 9.6577 |
|  | `mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2` | K0 163 -> 360.366; K0_prime 4.02 -> 9.41478 |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | `goethite_gleason_2008_bm3_1` | rt_eos.K0 140.3 -> 183.338; rt_eos.K0_prime 4.6 -> 0 |
| [Gong et al. (2004)](https://doi.org/10.1029/2003gl019132) | `mg092fe008sio3_bridgmanite_gong_2004_adiabatic_bm3` | V0 163.137 -> 231.512 |
| [Jacobsen et al. (2002)](https://doi.org/10.1029/2001jb000490) | `ferropericlase_fe27_jacobsen_2002_bm3_1` | V0 76.336 -> 113.015; K0 158.4 -> 37.1717 |
|  | `magnesiowustite_fe56_jacobsen_2002_bm3_1` | V0 77.453 -> 113.002; K0 155.8 -> 37.1669 |
|  | `magnesiowustite_fe75_jacobsen_2002_bm3_1` | V0 78.082 -> 113.014; K0 151.3 -> 37.188 |
| [Jacobsen et al. (2005)](https://doi.org/10.1107/s0909049505022326) | `fe093o_b1_jacobsen_2005_bm3_1` | V0 79.41 -> 59.7395 |
|  | `mg073fe027o_jacobsen_2005_bm3_1` | V0 77.3 -> 57.963; K0_prime 4 -> 2.38146 |
| [Katsura et al. (2009)](https://doi.org/10.1029/2009gl038107) | `wadsleyite_katsura_2009_bm3_1` | gamma0 1.64 -> 1.12567 |
| [Ono et al. (2000)](https://doi.org/10.1007/s002690000108) | `sno2_cubic_27gpa_ono_2000_bm3_1` | K0 252 -> 379.59 |
|  | `sno2_pa_3_at_48gpa_ono_2000_bm3_1` | K0 252 -> 379.59 |
| [Solomatova et al. (2016)](https://doi.org/10.2138/am-2016-5510) | `mg0490fe0483ti0027o_solomatova_2016_hs_bm3_reference_1` | K0 160 -> 202.391; K0_prime 4.12 -> 2.31737 |
|  | `mg0490fe0483ti0027o_solomatova_2016_ls_bm3_reference_2` | V0 73.64 -> 78.5805; K0 173 -> 135.945 |
| [Somayazulu et al. (2023)](https://doi.org/10.1098/rsta.2022.0331) | `b4c_somayazulu_2023_bm3_1` | q 2.1 -> 1.04991 |
| [Thompson et al. (2017)](https://doi.org/10.1002/2017jb014168) | `e_feooh_hc_low_spin_thompson_2017_bm3_1` | K0 223 -> 185.921 |

The full data selection, model mapping, residuals, and bounded explanation
for every row above are in the
[primary EOS refit validation](primary-eos-refits.md#detailed-non-parity-investigations).

## Papers with unavailable direct refits

These **81 papers** contain 564 records for which a
source-faithful coefficient refit could not be performed. A paper can also
have other records that were reproduced.

| Paper | Affected records | Why direct refitting was unavailable |
|---|---|---|
| [Akber-Knutson et al. (2002)](https://doi.org/10.1029/2001gl013523) | `ca_perovskite_akber_knutson_2002_vib_pm3m_300k_bm3`, `ca_perovskite_akber_knutson_2002_vibc_pm3m_300k_bm3`, `ca_perovskite_akber_knutson_2002_vibc_pnma_300k_bm3`, `stishovite_akber_knutson_2002_vibc_300k_bm3` | The paper does not tabulate the calculated energy-volume grid. The paper tabulates coefficients but not the calculated energy-volume grid. |
| [Akber-Knutson et al. (2005)](https://doi.org/10.1029/2005gl023192) | `al2o3_perovskite_akber_knutson_2005_gga_bm3`, `al2o3_post_perovskite_akber_knutson_2005_gga_bm3`, `bridgmanite_akber_knutson_2005_gga_bm3`, `mg09375al0125si09375o3_bridgmanite_akber_knutson_2005_gga_bm3`, `mg09375al0125si09375o3_post_perovskite_akber_knutson_2005_gga_bm3`, `mgal00625h00625si09375o3_bridgmanite_akber_knutson_2005_gga_bm3`, `mgal00625h00625si09375o3_post_perovskite_akber_knutson_2005_gga_bm3`, `mgal0125si0875o29375_bridgmanite_akber_knutson_2005_gga_bm3`, `mgal0125si0875o29375_post_perovskite_akber_knutson_2005_gga_bm3`, `mgsio3_post_perovskite_akber_knutson_2005_gga_bm3` | The source states the E(V) volume envelope but does not tabulate the individual calculated E(V) observations, fit weights, residuals, or covariance. The E(V) rows, weights, residuals, and covariance are not tabulated. |
| [Akins et al. (2004)](https://doi.org/10.1029/2004gl020237) | `mgsio3_liquid_akins_2004_adiabatic_bm3` | The publisher deposits the new shock states and the derivation used to calculate theoretical Hugoniots; this record preserves the printed candidate parameters without refitting. |
| [Anderson et al. (1989)](https://doi.org/10.1063/1.342969) | `gold_anderson_1989_bm3_1` | Tables I-IV contain heterogeneous literature properties, separately regressed coefficients, and derived thermodynamic diagnostics, while Table V is output from Equation (29). The source performs staged smoothing, one-dimensional regressions, numerical integrations, and qualitative K0' trials; it does not define a global observation matrix, objective, weights, integration protocol, or covariance that could be reproduced as a direct EOS coefficient refit. |
| [Anzellini et al. (2025)](https://doi.org/10.1038/s43246-025-00963-4) | `iridium_anzellini_2025_bm3_1` | The bundled rows are all heated states. The stored coefficients are the 300 K reference part of a combined thermal fit, but the record does not represent the source's thermal correction needed to refit those rows. |
| [B1 Fe0.94O, Fischer et al. (2011)](https://doi.org/10.1016/j.epsl.2011.02.025) | `feo_b8_2_fischer_2011_bm3_1`, `feo_fischer_2011_bm3_2` | Only 1 observation(s) lie at the reference temperature for 2 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |
| [Baty et al. (2024)](https://doi.org/10.1063/5.0179469) | `palladium_baty_2024_bm3_dft_2` | The calculated P(V) grid is not published as independent row-level fit input. Table S3 contains pressures generated from the already fitted EOS at selected volumes. |
| [Benedict et al. (2014)](https://doi.org/10.1103/physrevb.89.224109) | `diamond_benedict_2014_double_debye_4` | This is a theoretical multiphase carbon EOS. It publishes the fitted diamond model coefficients, but not the underlying electronic-structure grid as row-level data. |
| [Bykova et al. (2018)](https://doi.org/10.1038/s41467-018-07265-z) | `coesite_i_iii_bykova_2018_300k_bm3`, `coesite_v_bykova_2018_am05_static_bm3_refit` | The combined observations are plotted in Figure 1 but not supplied as a numerical table. Table 10 contains only one coesite-V pressure-volume anchor, which is sufficient to verify the published-parameter reconstruction but not to independently refit its three BM3 coefficients. |
| [Caracas et al. (2005)](https://doi.org/10.1029/2004gl022144) | `ca_perovskite_caracas_2005_bm3_3`, `ca_perovskite_caracas_2005_bm4_4`, `ca_perovskite_tetragonal_caracas_2005_bm3_1`, `ca_perovskite_tetragonal_caracas_2005_bm4_2`, `casio3_perovskite_i4mmm_caracas_2005_bm3_1`, `casio3_perovskite_i4mmm_caracas_2005_bm4_2`, `casio3_perovskite_im3_caracas_2005_bm3_1`, `casio3_perovskite_im3_caracas_2005_bm4_2`, `casio3_perovskite_imma_caracas_2005_bm3_1`, `casio3_perovskite_imma_caracas_2005_bm4_2`, `casio3_perovskite_p42nmc_caracas_2005_bm3_1`, `casio3_perovskite_p42nmc_caracas_2005_bm4_2`, `casio3_perovskite_p4mbm_caracas_2005_bm3_1`, `casio3_perovskite_p4mbm_caracas_2005_bm4_2`, `casio3_perovskite_pnma_caracas_2005_bm3_1`, `casio3_perovskite_pnma_caracas_2005_bm4_2`, `casio3_perovskite_r3c_caracas_2005_bm3_1`, `casio3_perovskite_r3c_caracas_2005_bm4_2` | The article publishes fitted EOS coefficients and relative energies but not the first-principles E(V) observations, fit weights, residuals, or covariance. No supporting-information or official data attachment is listed on the publisher article page or the UCL deposit, so an independent coefficient refit is not possible. Complete coefficients and density checkpoints are published; the underlying E(V) grid is not. The source publishes complete coefficients and independent density checkpoints but not the underlying E(V) grid or fit covariance. |
| [Chantel et al. (2012)](https://doi.org/10.1029/2012gl053075) | `bridgmanite_chantel_2012_bm3_mgd` | The bundled density and acoustic-velocity observations validate the published thermoelastic pressure surface in the dedicated Chantel reproduction. The stored K0 and K0-prime come from the source's combined acoustic fit, so these rows are not independent observations for a generic pressure-volume coefficient refit. |
| [Chen et al. (2018)](https://doi.org/10.2138/am-2018-6087) | `ca_perovskite_tetragonal_chen_2018_bm2` | The article states that all EOS observations are in Table 1, but no authoritative open table artifact suitable for lossless redistribution was retrieved in this audit; no digitized pseudo-table was created. |
| [Chizmeshya et al. (1996)](https://doi.org/10.1029/96gl02624) | `ca_perovskite_chizmeshya_1996_lapw7_static_bm3`, `ca_perovskite_chizmeshya_1996_lapw8_static_bm3`, `ca_perovskite_chizmeshya_1996_lapw9_300k_bm3`, `ca_perovskite_chizmeshya_1996_lapw9_300k_kp4_bm3`, `ca_perovskite_chizmeshya_1996_lapw9_static_bm3` | Coefficients are tabulated; energy-volume points and weights are unavailable. No E-V grid is published. The coefficients are complete but the corrected E-V grid is not published. Thermally corrected E-V points are not published. |
| [Cohen and Lin (2014)](https://doi.org/10.1103/physrevb.90.140102) | `fesio3_bridgmanite_cohen_lin_2014_vinet_1`, `fesio3_post_perovskite_cohen_lin_2014_vinet_1`, `fesio3_post_perovskite_ii_cohen_lin_2014_vinet_1` | The eight energy-volume observations are plot-only, but Table III independently prints V100=34.27 A3/FeSiO3, K100=597 GPa, and K100'=3.34; the stored curve reproduces all three after rounding. The energy-volume observations are plot-only; Table III independently prints V100=33.98 A3/FeSiO3, K100=579 GPa, and K100'=3.47. The energy-volume grid is plot-only; Table III independently prints V100=34.49 A3/FeSiO3, K100=580 GPa, and K100'=3.44. |
| [Correa et al. (2008)](https://doi.org/10.1103/physrevb.78.024101) | `diamond_correa_2008_double_debye_log_moment_5` | This is a theoretical multiphase carbon EOS. It publishes fitted cold-curve and thermal-model coefficients, but not the underlying DFT energy-volume grid as row-level data. |
| [Datchi et al. (2007)](https://doi.org/10.1103/physrevb.75.214104) | `diamond_datchi_2007_vinet_1` | For the diamond record, Datchi et al. reanalyze the previously published Occelli et al. diamond compression data on the H05 pressure scale and report the resulting EOS parameters in Table II; they do not republish the row-level diamond observations. The paper's new c-BN observations are separately bundled with the c-BN record. |
| [Delta archive (2023)](https://doi.org/10.24435/materialscloud:5e-mv) | `aluminum_delta_archive_experimental_reference_bm3`, `barium_bcc_delta_archive_experimental_reference_bm3`, `bismuth_a7_delta_archive_experimental_reference_bm3`, `cadmium_hcp_delta_archive_experimental_reference_bm3`, `calcium_fcc_delta_archive_experimental_reference_bm3`, `cesium_bcc_delta_archive_experimental_reference_bm3`, `chromium_delta_archive_experimental_reference_bm3`, `cobalt_hcp_delta_archive_experimental_reference_bm3`, `copper_delta_archive_experimental_reference_bm3`, `fe_delta_archive_experimental_reference_bm3`, `gold_delta_archive_experimental_reference_bm3`, `hafnium_hcp_delta_archive_experimental_reference_bm3`, `indium_bct_delta_archive_experimental_reference_bm3`, `iridium_delta_archive_experimental_reference_bm3`, `lead_fcc_delta_archive_experimental_reference_bm3`, `lithium_9r_delta_archive_experimental_reference_bm3`, `magnesium_hcp_delta_archive_experimental_reference_bm3`, `manganese_alpha_delta_archive_experimental_reference_bm3`, `molybdenum_delta_archive_experimental_reference_bm3`, `nickel_delta_archive_experimental_reference_bm3`, `niobium_delta_archive_experimental_reference_bm3`, `osmium_delta_archive_experimental_reference_bm3`, `palladium_delta_archive_experimental_reference_bm3`, `platinum_delta_archive_experimental_reference_bm3`, `potassium_bcc_delta_archive_experimental_reference_bm3`, `rhenium_delta_archive_experimental_reference_bm3`, `rhodium_delta_archive_experimental_reference_bm3`, `rubidium_bcc_delta_archive_experimental_reference_bm3`, `ruthenium_delta_archive_experimental_reference_bm3`, `scandium_hcp_delta_archive_experimental_reference_bm3`, `silver_delta_archive_experimental_reference_bm3`, `sodium_9r_delta_archive_experimental_reference_bm3`, `strontium_fcc_delta_archive_experimental_reference_bm3`, `tantalum_delta_archive_experimental_reference_bm3`, `thallium_hcp_delta_archive_experimental_reference_bm3`, `tin_alpha_delta_archive_experimental_reference_bm3`, `titanium_alpha_delta_archive_experimental_reference_bm3`, `tungsten_delta_archive_experimental_reference_bm3`, `vanadium_bcc_delta_archive_experimental_reference_bm3`, `yttrium_hcp_delta_archive_experimental_reference_bm3`, `zinc_hcp_delta_archive_experimental_reference_bm3`, `zirconium_alpha_delta_archive_experimental_reference_bm3` | The archive supplies only a constructed coefficient table. Exact pre-correction inputs, a unified P-V dataset, pressure calibration, uncertainties, weights, and selection rules are not available, so an independent fit of the composite triplet is impossible. Property-level upstream recovery is recorded separately and includes row-level data only where the cited source actually prints it. |
| [Deng et al. (2006)](https://doi.org/10.1088/0256-307x/23/8/101) | `bridgmanite_deng_2006_lda_bm3` | Calculations were performed every 10 GPa through 120 GPa, but their numerical volume grid is plotted rather than tabulated. |
| [Dewaele (2019)](https://doi.org/10.3390/min9110684) | `aluminum_dewaele_2019_dor_vinet`, `aluminum_dewaele_2019_mao_vinet`, `beryllium_hcp_dewaele_2019_dor_vinet`, `beryllium_hcp_dewaele_2019_mao_vinet`, `cobalt_hcp_dewaele_2019_dor_vinet`, `cobalt_hcp_dewaele_2019_mao_vinet`, `copper_dewaele_2019_dor_vinet`, `copper_dewaele_2019_mao_vinet`, `gold_dewaele_2019_dor_vinet`, `gold_dewaele_2019_mao_vinet`, `molybdenum_dewaele_2019_dor_vinet`, `molybdenum_dewaele_2019_mao_vinet`, `nickel_dewaele_2019_dor_vinet`, `nickel_dewaele_2019_mao_vinet`, `platinum_dewaele_2019_dor_vinet`, `platinum_dewaele_2019_mao_vinet`, `rhenium_dewaele_2019_dor_vinet`, `rhenium_dewaele_2019_mao_vinet`, `silver_dewaele_2019_dor_vinet`, `silver_dewaele_2019_mao_vinet`, `tantalum_dewaele_2019_dor_vinet`, `tantalum_dewaele_2019_mao_vinet`, `tungsten_dewaele_2019_dor_vinet`, `tungsten_dewaele_2019_mao_vinet`, `zinc_hcp_dewaele_2019_dor_vinet`, `zinc_hcp_dewaele_2019_mao_vinet` | Published coefficients are executable; this import does not claim a new refit of row-level observations. |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.77.094106) | `diamond_benedict_2014_dewaele_anchored`, `diamond_correa_2008_dewaele_anchored` | The linked diffraction rows constrain only the Dewaele reference isotherm; the Benedict thermal term is a separately published theoretical model. The linked diffraction rows constrain only the Dewaele reference isotherm; the Correa thermal term is a separately published theoretical model. |
| [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/physrevb.75.024115) | `platinum_dorogokupets_oganov_2007_vinet_4` | The 36 bundled Dewaele static rows are a recoverable stage-2 subset, and the seven Holmes shots are an independent Equation-15 diagnostic. The published result is nevertheless a two-stage weighted optimization of six metals: the complete observation inventory and weights, selected Shock Wave Database rows, Collard-McLellan Figure 1 K_S(T) coordinates, and platinum-specific thermochemical row mapping are not published. The dedicated audit therefore validates exact calculated outputs and partial source subsets without claiming a numerical reconstruction of the global objective. |
| [Dorogokupets et al. (2015)](https://doi.org/10.1016/j.rgg.2015.01.011) | `akimotoite_dorogokupets_2015_298k_rydberg_stacey`, `bridgmanite_dorogokupets_2015_298k_rydberg_stacey`, `mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey` | Complete coefficients are printed, but no consolidated observation-level fit dataset and weights are deposited. The complete optimized coefficients are tabulated, while the heterogeneous literature fit observations and weights are not deposited as one machine-readable dataset. Complete optimized coefficients are tabulated; the underlying heterogeneous P-V-T observations and fitting weights are not deposited together. |
| [Driver et al. (2010)](https://doi.org/10.1073/pnas.0912130107) | `alpha_quartz_driver_2010_qmc_300k_vinet`, `seifertite_driver_2010_qmc_300k_vinet`, `sio2_stv_andr_driver_2010_qmc_300k_vinet` | The source describes approximately six volumes spanning +/-10% and plots the statistical envelope, but does not tabulate the individual QMC energies or pressures. |
| [Fortes (2019)](https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf) | `lead_fcc_fortes_2019_bm4_1` | Fortes (2019) derives an fcc-Pb pressure scale from published literature data and tabulates model coefficients and comparisons, but no new row-level experimental P-V-T observations. |
| [Fu et al. (2023)](https://doi.org/10.2138/am-2022-8435) | `ca_perovskite_fu_2023_bm3_mgd_refit` | The current paper does not reprint all literature P-V-T rows or the fitting weights; executable analytical checkpoints verify the published composite coefficients without claiming exact refit parity. |
| [Funamori et al. (1998)](https://doi.org/10.1029/98jb01575) | `mgal2o4_cafe2o4_funamori_1998_bm2_1`, `mgal2o4_cati2o4_funamori_1998_bm2_1` | The primary article reports only the ambient and compressed endpoint for this polymorph. Those two states reproduce the published fixed-V0, fixed-K0-prime curve in the dedicated Funamori reproduction, but do not provide enough degrees of freedom for the generic refit campaign. |
| [Ghosh and Karki (2016)](https://doi.org/10.1038/srep37269) | `mgo_liquid_ghosh_karki_2016_3000k_bm3_1` | The source plots but does not tabulate the pure-liquid P-V simulation states; no graphical pseudo-precision was introduced. |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | `e_feooh_gleason_2008_bm2_1` | Only 1 observation(s) lie at the reference temperature for 2 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |
| [Guigue et al. (2020)](https://doi.org/10.1063/1.5138697) | `palladium_guigue_2020_vinet_1` | The underlying pure-Pd observations are plotted but not tabulated in the accessible primary article; no numerical refit is claimed. |
| [Hama and Suito (1996)](https://doi.org/10.1088/0953-8984/8/1/008) | `mgo_hama_suito_1996_qsm_static_vinet` | The theoretical comparison curves are plotted but no pressure-volume calculation table is published; validation uses the exact printed equation and coefficients without inventing pseudo-observations. |
| [Hamahata et al. (2000)](https://doi.org/10.2465/jmps.95.236) | `bridgmanite_hamahata_2000_md_300k_bm3` | The complete fit coefficients are printed; the simulated 300 K volumes are plotted but not tabulated, so no pseudo-observations were digitized. |
| [Holland et al. (2013)](https://doi.org/10.1093/petrology/egt035) | `al2o3_perovskite_holland_2013_apv_modified_tait`, `bridgmanite_holland_2013_mpv_modified_tait`, `ca_perovskite_holland_2013_cpv_modified_tait`, `feo_holland_2013_fper_modified_tait`, `fesio3_bridgmanite_holland_2013_fpv_modified_tait`, `mgo_holland_2013_per_modified_tait`, `sio2_stv_andr_holland_2013_stv_modified_tait` | The official tc-ds62 apv block provides all four coefficients. The official tc-ds62 mpv block provides V0, K0, K0', and K0'' without model-family conversion. The official tc-ds62 cpv block provides all four coefficients. The official tc-ds62 fper block provides all four coefficients. The official tc-ds62 fpv block provides all four coefficients. The official tc-ds62 per block provides all four coefficients. The official tc-ds62 stv block provides all four coefficients. |
| [Holmes et al. (1989)](https://doi.org/10.1063/1.344177) | `platinum_holmes_1989_vinet_1` | The bundled rows are shock-Hugoniot qualification experiments; the stored equilibrium Vinet curve is a theoretical 300 K isotherm and cannot be refitted directly to those rows. |
| [Ismailova et al. (2016)](https://doi.org/10.1126/sciadv.1600427) | `fe088sio3_bridgmanite_ismailova_2016_300k_bm2` | The source supplies crystallographic and compressibility data in its supplement; the production record preserves the directly reported fitted parameterization. |
| [Karki and Crain (1998)](https://doi.org/10.1029/98gl51952) | `ca_perovskite_karki_crain_1998_static_bm3` | The source reports the complete EOS coefficients and plotted calculated curves but no reusable energy-volume table or fit covariance. |
| [Karki and Wentzcovitch (2002)](https://doi.org/10.1029/2001jb000702) | `akimotoite_karki_2002_1000k_bm4_3`, `akimotoite_karki_2002_2000k_bm4_4`, `akimotoite_karki_2002_300k_bm4_2`, `akimotoite_karki_2002_static_bm4_1` | The fitted coefficients and plotted curve are published, but the numerical free-energy-volume grid is unavailable. The fitted coefficients and curve are published, but the numerical free-energy-volume grid is unavailable. Complete BM4 coefficients are tabulated, but the underlying free-energy-volume grid is shown only graphically and is not deposited. |
| [Karki et al. (1997)](https://doi.org/10.2138/am-1997-1-207) | `mgo_karki_1997_lda_static_bm3`, `mgo_karki_1997_lda_static_bm4` | The calculated states are plotted but not tabulated; no pseudo-observations were constructed. The calculated states are plotted in Figure 1 but not tabulated; validation therefore checks the published equation, coefficients, derivative identities, and inverse curve rather than inventing row-level data. |
| [Kawai and Tsuchiya (2014)](https://doi.org/10.1002/2013jb010905) | `ca_perovskite_kawai_2014_vinet_mgd_3` | The publisher page exposes no supporting-information or data file, and the article plots but does not tabulate the underlying FPMD P-V-T stress averages. A direct refit is therefore impossible. All 60 printed Table 1 fitted-isotherm benchmark states are bundled separately for numerical reproduction; they are model values, not primary observations. |
| [Kiefer et al. (2002)](https://doi.org/10.1029/2002gl014683) | `bridgmanite_kiefer_2002_gga_bm3`, `mg075fe025sio3_bridgmanite_kiefer_2002_gga_bm3` | Complete coefficients are printed, but the underlying energy-volume grid is not tabulated. |
| [Lee and Wan (2008)](https://doi.org/10.1103/physrevb.78.224103) | `mgo_lee_wan_2008_gga_static_bm3`, `mgo_lee_wan_2008_lda_static_bm3` | The exact BM3 coefficients and plotted curves are published, but the underlying pressure-volume calculation table and fitting weights are not. |
| [Lejaeghere et al. (2016)](https://doi.org/10.1126/science.aad3000) | `aluminum_lejaeghere_2016_fleur_pbe_bm3`, `aluminum_lejaeghere_2016_wien2k_pbe_bm3`, `barium_bcc_lejaeghere_2016_wien2k_pbe_bm3`, `beryllium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `bismuth_a7_lejaeghere_2016_wien2k_pbe_bm3`, `cadmium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `calcium_fcc_lejaeghere_2016_wien2k_pbe_bm3`, `cesium_bcc_lejaeghere_2016_wien2k_pbe_bm3`, `chromium_lejaeghere_2016_wien2k_pbe_bm3`, `cobalt_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `copper_lejaeghere_2016_wien2k_pbe_bm3`, `fe_lejaeghere_2016_fleur_pbe_bm3`, `fe_lejaeghere_2016_wien2k_pbe_bm3`, `gallium_alpha_lejaeghere_2016_wien2k_pbe_bm3`, `gold_lejaeghere_2016_wien2k_pbe_bm3`, `hafnium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `indium_bct_lejaeghere_2016_wien2k_pbe_bm3`, `iridium_lejaeghere_2016_wien2k_pbe_bm3`, `lead_fcc_lejaeghere_2016_wien2k_pbe_bm3`, `lithium_9r_lejaeghere_2016_wien2k_pbe_bm3`, `lutetium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `magnesium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `manganese_fcc_afm_lejaeghere_2016_wien2k_pbe_bm3`, `mercury_bct_dft_lejaeghere_2016_wien2k_pbe_bm3`, `molybdenum_lejaeghere_2016_wien2k_pbe_bm3`, `nickel_lejaeghere_2016_wien2k_pbe_bm3`, `niobium_lejaeghere_2016_wien2k_pbe_bm3`, `osmium_lejaeghere_2016_wien2k_pbe_bm3`, `palladium_lejaeghere_2016_wien2k_pbe_bm3`, `platinum_lejaeghere_2016_wien2k_pbe_bm3`, `polonium_simple_cubic_lejaeghere_2016_wien2k_pbe_bm3`, `potassium_bcc_lejaeghere_2016_wien2k_pbe_bm3`, `rhenium_lejaeghere_2016_wien2k_pbe_bm3`, `rhodium_lejaeghere_2016_wien2k_pbe_bm3`, `rubidium_bcc_lejaeghere_2016_wien2k_pbe_bm3`, `ruthenium_lejaeghere_2016_wien2k_pbe_bm3`, `scandium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `silver_lejaeghere_2016_wien2k_pbe_bm3`, `sodium_9r_lejaeghere_2016_wien2k_pbe_bm3`, `strontium_fcc_lejaeghere_2016_wien2k_pbe_bm3`, `tantalum_lejaeghere_2016_wien2k_pbe_bm3`, `technetium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `thallium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `tin_alpha_lejaeghere_2016_wien2k_pbe_bm3`, `titanium_alpha_lejaeghere_2016_wien2k_pbe_bm3`, `tungsten_lejaeghere_2016_wien2k_pbe_bm3`, `vanadium_bcc_lejaeghere_2016_wien2k_pbe_bm3`, `yttrium_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `zinc_hcp_lejaeghere_2016_wien2k_pbe_bm3`, `zirconium_alpha_lejaeghere_2016_wien2k_pbe_bm3` | The official archive preserves the complete EOS coefficients and frozen structure, but not the seven row-level E(V) values for this code; direct coefficient refitting is unavailable. |
| [Leonov et al. (2017)](https://doi.org/10.1103/physrevb.96.075136) | `feo_leonov_2017_hs_bm3_1`, `mg0125fe0875o_leonov_2017_hs_bm3_1`, `mg025fe075o_leonov_2017_hs_bm3_1`, `mg0375fe0625o_leonov_2017_hs_bm3_1`, `mg05fe05o_leonov_2017_hs_bm3_1`, `mg0625fe0375o_leonov_2017_hs_bm3_1`, `mg075fe025o_leonov_2017_hs_bm3_1`, `mg0875fe0125o_leonov_2017_hs_bm3_1` | Independent equation checkpoints verify transcription without claiming a refit. Coefficients and plotted energy-volume curves are published, but the numerical grid is not deposited; independent equation checkpoints verify transcription without claiming a refit. Independent checkpoints only; no grid deposited. |
| [Li and Zeng (2009)](https://doi.org/10.1142/s0129183109014242) | `bridgmanite_li_zeng_2009_gga_bm3`, `bridgmanite_li_zeng_2009_gga_natural_strain3`, `bridgmanite_li_zeng_2009_gga_vinet`, `mgsio3_post_perovskite_li_zeng_2009_gga_bm3`, `mgsio3_post_perovskite_li_zeng_2009_gga_natural_strain3`, `mgsio3_post_perovskite_li_zeng_2009_gga_vinet` | Figure 1 plots E(V), but the individual energy grid is not tabulated. |
| [Li et al. (2006)](https://doi.org/10.1029/2005jb004251) | `mgo_li_2006_bm3_absolute_acoustic` | The Table 1 pressures are outputs of the stored acoustic-derived BM3, not independent pressure-volume observations. The source-derived isothermal coefficients are instead validated by the bundled velocity-density data and the dedicated acoustic finite-strain reproduction. |
| [Liu et al. (2007)](https://doi.org/10.1088/0953-8984/19/24/246103) | `ca_perovskite_liu_2007_lda_static_bm3` | The source plots its calculated EOS but does not tabulate the energy-volume grid. |
| [Liu et al. (2010)](https://doi.org/10.1142/s0217984910022391) | `mgsio3_post_perovskite_liu_2010_lda_static_bm3` | Figure 1 plots the EOS, but the underlying calculated E-V points and regression covariance are not tabulated; validation uses the exact printed equation and coefficients. |
| [Liu et al. (2011)](https://doi.org/10.1088/1674-0068/24/06/703-710) | `bridgmanite_liu_2011_gga_static_bm3` | The calculated 0-150 GPa volume series is plotted in Figure 1 but not tabulated; Table I provides the complete fitted coefficients. |
| [Luo et al. (2023)](https://doi.org/10.1103/physrevb.107.134116) | `mgo_b1_luo_2023_vinet_thermal_5` | The five bundled Table I rows are only the new shock subset of a global quasi-Debye fit. The complete earlier-study observations, numerical sound-velocity-density fits, objective weights, and covariance are not published; Tables II-III are derived EOS output and cannot serve as independent refit observations. |
| [Mao et al. (1991)](https://doi.org/10.1029/91jb00176) | `bridgmanite_mao_1991_bm2_1`, `mg080fe020sio3_bridgmanite_mao_1991_bm2_1`, `mg09fe01sio3_bridgmanite_mao_1991_bm2_2` | No accessible numerical pressure-volume table was available; no graphical points were fabricated. |
| [Marcondes et al. (2020)](https://doi.org/10.1103/physrevb.102.104112) | `mg09375fe00625o_marcondes_2020_11nn_hs_bm3_1`, `mg09375fe00625o_marcondes_2020_11nn_ls_bm3_2`, `mg09375fe00625o_marcondes_2020_11nn_ms_bm3_3`, `mg09375fe00625o_marcondes_2020_2nn_hs_bm3_4`, `mg09375fe00625o_marcondes_2020_2nn_ls_bm3_5`, `mg09375fe00625o_marcondes_2020_2nn_ms_bm3_6`, `mg096875fe003125o_marcondes_2020_hs_bm3_1`, `mg096875fe003125o_marcondes_2020_ls_bm3_2` | The source publishes complete coefficients but no numerical energy-volume grid; independent BM3 checkpoints verify every stored curve. |
| [Metsue and Tsuchiya (2012)](https://doi.org/10.1111/j.1365-246x.2012.05511.x) | `bridgmanite_metsue_2012_static_bm3_1`, `mg09375fe00625sio3_bridgmanite_metsue_2012_hs_model1_bm3`, `mg09375fe00625sio3_bridgmanite_metsue_2012_hs_model2_bm3`, `mg09375fe00625sio3_bridgmanite_metsue_2012_hs_model3_bm3`, `mg09375fe00625sio3_bridgmanite_metsue_2012_ls_model1_bm3`, `mg09375fe00625sio3_bridgmanite_metsue_2012_ls_model2_bm3`, `mg09375fe00625sio3_bridgmanite_metsue_2012_ls_model3_bm3` | The coefficients and six calculation pressures are stated; row-wise calculated P-V values are not tabulated. The six calculation pressures and full coefficients are stated; row-wise calculated P-V values are not tabulated. Coefficients and calculation pressures are given; row-wise P-V values are not tabulated. |
| [Mookherjee et al. (2015)](https://doi.org/10.2138/am-2015-5312) | `mgsioh6_365a_phase_mookherjee_2015_gga_bm4_model_crystal_3` | The full four-coefficient BM4 is tabulated and the computed P-V markers are plotted, but the numerical energy-volume grid is not deposited; no false-precision digitization is bundled. |
| [Mosenfelder et al. (2009)](https://doi.org/10.1029/2008jb005900) | `mgsio3_post_perovskite_mosenfelder_2009_bm3_1` | The bundled rows are shock states and the source's thermal reduction cannot be reconstructed as a direct P-V-T least-squares fit because most rows do not report temperature. |
| [Muñoz and Kunc (1993)](https://doi.org/10.1088/0953-8984/5/33/010) | `indium_nitride_munoz_1993_murnaghan_1` | This is a first-principles study. The calculated E(V) points are plotted but not tabulated; Table 1 contains only the fitted theoretical parameters. |
| [Noguchi et al. (1999)](https://doi.org/10.1016/s0022-3697(98)00296-0) | `nickel_oxide_noguchi_1999_bm3_1` | The bundled rows are Hugoniot states; the stored 300 K isotherm is the source's Mie-Gruneisen reduction, not a direct fit to Hugoniot P-V pairs. |
| [Noguchi et al. (2013)](https://doi.org/10.1007/s00269-012-0549-1) | `ca_perovskite_noguchi_2013_bm2_mgd_1` | Table 1 contains 54 P-V-T rows, paired Fei/Holmes pressures, Pt lattice parameters, and three explicit fit exclusions. A complete local transcription was used for the audit and independent refit, but is not redistributed because the subscription article states no reusable data license. |
| [Oganov and Dorogokupets (2003)](https://doi.org/10.1103/physrevb.67.224110) | `mgo_oganov_2003_ecp_large_core_static_bm3`, `mgo_oganov_2003_ecp_small_core_static_bm3`, `mgo_oganov_2003_paw_large_core_static_bm3`, `mgo_oganov_2003_pressure_corrected_0k_bm3`, `mgo_oganov_2003_pressure_corrected_1000k_bm3`, `mgo_oganov_2003_pressure_corrected_2000k_bm3`, `mgo_oganov_2003_pressure_corrected_298k_bm3`, `mgo_oganov_2003_pressure_corrected_3000k_bm3`, `mgo_oganov_2003_pressure_corrected_4000k_bm3` | The complete coefficients are tabulated but the fitted energy-volume grid is not published. No numerical energy-volume grid is published. The coefficients and plotted curve are published; no row-level thermal grid is deposited. No row-level thermal grid is published. Coefficients and curve are published; underlying thermal grid is not. |
| [Ono (2013)](https://doi.org/10.3390/e15104300) | `ca_perovskite_ono_2013_bm3_log_thermal` | The article states that 27 high-temperature AIMD states were fitted but does not tabulate their P-V-T values, fit weights, residual statistic, or covariance. |
| [Redfern et al. (1993)](https://doi.org/10.1029/93gl02507) | `mg0991fe0008mn0001co3_redfern_1993_bm2_1`, `mg0991fe0008mn0001co3_redfern_1993_bm3_2` | The accessible primary record did not expose a numerical pressure-volume table; no figure points were invented. |
| [Ricolleau et al. (2009)](https://doi.org/10.1029/2008gl036759) | `klb1_ca_perovskite_ricolleau_2009_bm2_alphakt`, `klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt`, `klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2`, `klb1_mg_perovskite_ricolleau_2009_bm2_alphakt` | The source reports the P-V-T observations in Table S1; this record preserves the Table 1 fit. The source reports the P-V observations in Table S1; this record preserves the Table 1 fit. The source reports 17 room-temperature and 136 high-temperature patterns in Table S1; this record preserves the Table 1 fit. |
| [Sagatova et al. (2021)](https://doi.org/10.1134/s0016702921080073) | `breyite_sagatova_2021_gga_300k_vinet`, `ca_perovskite_tetragonal_sagatova_2021_gga_300k_vinet`, `ca_perovskite_tetragonal_sagatova_2021_lda_300k_vinet`, `casio2o5_titanite_sagatova_2021_gga_300k_vinet`, `larnite_sagatova_2021_gga_300k_vinet`, `pseudowollastonite_sagatova_2021_gga_300k_vinet`, `wollastonite_sagatova_2021_gga_300k_vinet` | Relative-volume curves are plotted in Figure 6; calculated P-V points are not tabulated. The 300 and 2000 K P-V curves are plotted in Figure 7, but calculated P-V points are not tabulated. |
| [Sakai et al. (2025)](https://doi.org/10.1038/s43246-025-00792-5) | `copper_sakai_2025_rydberg_stacey_1`, `gold_sakai_2025_rydberg_stacey_1`, `iron_sakai_2025_rydberg_stacey_1`, `mgo_sakai_2025_rydberg_stacey_1`, `molybdenum_sakai_2025_rydberg_stacey_1`, `nacl_b2_sakai_2025_rydberg_stacey_1`, `platinum_sakai_2025_rydberg_stacey_1`, `rhenium_sakai_2025_rydberg_stacey_1`, `tungsten_sakai_2025_rydberg_stacey_1` | The bundled grid is calculated from the reported coefficients. Raw simultaneous-volume observations exist in Tables S1-S8, but the multi-stage fit also incorporates earlier studies and does not publish one flat pressure-volume regression table or covariance. |
| [Schoelmerich et al. (2020)](https://doi.org/10.1038/s41598-020-66340-y) | `stishovite_schoelmerich_2020_shock_300k_bm3` | Primary Table 1 reports the shock observables and derived states; the production record preserves the published fitted parameterization, while a separate transcription was not needed for executable verification. |
| [Shen and Smith (2026)](https://doi.org/10.1103/fxgq-96sg) | `fe_shen_2026_vinet_1`, `gold_shen_2026_vinet_3`, `iron_shen_2026_vinet_2`, `mgo_shen_2026_vinet_3`, `molybdenum_shen_2026_vinet_1`, `nacl_b1_shen_2026_vinet_1`, `nacl_b2_shen_2026_vinet_2`, `platinum_shen_2026_vinet_2`, `tantalum_shen_2026_vinet_2`, `tungsten_shen_2026_vinet_3` | The workbook contains simultaneous volumes but no pressures, and the record declares its Cu anchor as reference_model_not_supported. |
| [Sherman (1993)](https://doi.org/10.1029/93jb02175) | `ca_perovskite_sherman_1993_basis_b_static_bm3` | The complete calculated E-V grid and weighting are not tabulated. |
| [Shim et al. (2002)](https://doi.org/10.1029/2002gl016148) | `casio3_perovskite_tetragonal_shim_2002_bm3_1` | The six pressures and plotted volume/c-axis-ratio observations appear in Figure 2, but numerical cell volumes are not tabulated; no pseudo-precision was introduced by digitizing the small four-page rendering. |
| [Sokolova et al. (2013)](https://doi.org/10.1016/j.rgg.2013.01.005) | `aluminum_sokolova_2013_holzapfel_2`, `copper_sokolova_2013_holzapfel_2`, `diamond_sokolova_2013_holzapfel_3`, `gold_sokolova_2013_holzapfel_4`, `mgo_sokolova_2013_holzapfel_4`, `molybdenum_sokolova_2013_holzapfel_2`, `niobium_sokolova_2013_holzapfel_2`, `platinum_sokolova_2013_holzapfel_3`, `silver_sokolova_2013_holzapfel_2`, `tantalum_sokolova_2013_holzapfel_3`, `tungsten_sokolova_2013_holzapfel_4` | This is an internally consistent multi-marker optimization. It publishes input constants and optimized EOS coefficients, but no new row-level experimental P-V-T observations; the calibration comparisons are graphical. |
| [Solomatova et al. (2016)](https://doi.org/10.2138/am-2016-5510) | `mg061fe039o_solomatova_2016_fei_hs_bm3_reference_1`, `mg061fe039o_solomatova_2016_fei_ls_bm3_reference_2`, `mg061fe039o_solomatova_2016_zhuravlev_hs_bm3_reference_3`, `mg061fe039o_solomatova_2016_zhuravlev_ls_bm3_reference_4`, `mg065fe035o_solomatova_2016_hs_bm3_reference_1`, `mg065fe035o_solomatova_2016_ls_bm3_reference_2`, `mg075fe025o_solomatova_2016_hs_bm3_reference_1`, `mg075fe025o_solomatova_2016_ls_bm3_reference_2`, `mg083fe017o_solomatova_2016_hs_bm3_reference_1`, `mg083fe017o_solomatova_2016_ls_bm3_reference_2`, `mg090fe010o_solomatova_2016_hs_bm3_reference_1`, `mg090fe010o_solomatova_2016_ls_bm3_reference_2`, `mgfe60o_solomatova_2016_hs_bm3_reference_1`, `mgfe60o_solomatova_2016_ls_bm3_reference_2` | The paper publishes the complete refitted branch coefficients and transition parameters but does not reprint the underlying earlier P-V rows; independent BM3 checkpoints verify coefficient transcription. Independent BM3 checkpoints verify transcription; underlying rows are not reprinted. |
| [Sun et al. (2010)](https://doi.org/10.1515/zna-2010-1-202) | `ag_sun_2010_low_bn`, `ag_sun_2010_low_mrs3`, `ag_sun_2010_low_sms3`, `ag_sun_2010_low_sms4`, `ag_sun_2010_low_vn`, `ag_sun_2010_mrs3_1`, `ag_sun_2010_sms3_1`, `ag_sun_2010_sms4_1`, `al2o3_sun_2010_sms3_1`, `al2o3_sun_2010_sms4_1`, `al_sun_2010_low_mrs3`, `al_sun_2010_low_sms3`, `al_sun_2010_low_sms4`, `al_sun_2010_low_vn`, `al_sun_2010_sms3_1`, `al_sun_2010_sms4_1`, `au_sun_2010_low_bn`, `au_sun_2010_low_mrs3`, `au_sun_2010_low_sms3`, `au_sun_2010_low_sms4`, `au_sun_2010_low_vn`, `au_sun_2010_mrs3_1`, `au_sun_2010_sms3_1`, `au_sun_2010_sms4_1`, `be_sun_2010_low_mrs3`, `be_sun_2010_low_sms3`, `be_sun_2010_low_sms4`, `be_sun_2010_low_vn`, `be_sun_2010_sms3_1`, `be_sun_2010_sms4_1`, `ca_sun_2010_low_mrs3`, `ca_sun_2010_low_sms3`, `ca_sun_2010_low_sms4`, `ca_sun_2010_low_vn`, `ca_sun_2010_sms3_1`, `ca_sun_2010_sms4_1`, `cd_sun_2010_low_mrs3`, `cd_sun_2010_low_sms3`, `cd_sun_2010_low_sms4`, `cd_sun_2010_low_vn`, `cd_sun_2010_sms3_1`, `cd_sun_2010_sms4_1`, `co_sun_2010_low_mrs3`, `co_sun_2010_low_sms3`, `co_sun_2010_low_sms4`, `co_sun_2010_low_vn`, `co_sun_2010_sms3_1`, `co_sun_2010_sms4_1`, `cr_sun_2010_low_mrs3`, `cr_sun_2010_low_sms3`, `cr_sun_2010_low_sms4`, `cr_sun_2010_low_vn`, `cr_sun_2010_sms3_1`, `cr_sun_2010_sms4_1`, `csbr_sun_2010_sms3_1`, `csbr_sun_2010_sms4_1`, `cu_sun_2010_low_bn`, `cu_sun_2010_low_mrs3`, `cu_sun_2010_low_sms3`, `cu_sun_2010_low_sms4`, `cu_sun_2010_low_vn`, `cu_sun_2010_mrs3_1`, `cu_sun_2010_sms3_1`, `cu_sun_2010_sms4_1`, `in_sun_2010_low_mrs3`, `in_sun_2010_low_sms3`, `in_sun_2010_low_sms4`, `in_sun_2010_low_vn`, `in_sun_2010_sms3_1`, `in_sun_2010_sms4_1`, `k_sun_2010_sms3_1`, `k_sun_2010_sms4_1`, `kf_sun_2010_sms3_1`, `kf_sun_2010_sms4_1`, `ki_sun_2010_sms3_1`, `ki_sun_2010_sms4_1`, `li_sun_2010_sms3_1`, `li_sun_2010_sms4_1`, `libr_sun_2010_sms3_1`, `libr_sun_2010_sms4_1`, `licl_sun_2010_sms3_1`, `licl_sun_2010_sms4_1`, `lif_sun_2010_sms3_1`, `lif_sun_2010_sms4_1`, `lii_sun_2010_sms3_1`, `lii_sun_2010_sms4_1`, `mg_sun_2010_low_mrs3`, `mg_sun_2010_low_sms3`, `mg_sun_2010_low_sms4`, `mg_sun_2010_low_vn`, `mg_sun_2010_sms3_1`, `mg_sun_2010_sms4_1`, `mgo_sun_2010_sms3_1`, `mgo_sun_2010_sms4_1`, `mo_sun_2010_low_bn`, `mo_sun_2010_low_mrs3`, `mo_sun_2010_low_sms3`, `mo_sun_2010_low_sms4`, `mo_sun_2010_low_vn`, `mo_sun_2010_mrs3_1`, `mo_sun_2010_sms3_1`, `mo_sun_2010_sms4_1`, `na_sun_2010_sms3_1`, `na_sun_2010_sms4_1`, `nabr_sun_2010_sms3_1`, `nabr_sun_2010_sms4_1`, `nacl_sun_2010_sms3_1`, `nacl_sun_2010_sms4_1`, `naf_sun_2010_sms3_1`, `naf_sun_2010_sms4_1`, `nai_sun_2010_sms3_1`, `nai_sun_2010_sms4_1`, `nb_sun_2010_low_mrs3`, `nb_sun_2010_low_sms3`, `nb_sun_2010_low_sms4`, `nb_sun_2010_low_vn`, `nb_sun_2010_sms3_1`, `nb_sun_2010_sms4_1`, `nd_sun_2010_sms3_1`, `nd_sun_2010_sms4_1`, `ni_sun_2010_low_mrs3`, `ni_sun_2010_low_sms3`, `ni_sun_2010_low_sms4`, `ni_sun_2010_low_vn`, `ni_sun_2010_sms3_1`, `ni_sun_2010_sms4_1`, `pb_sun_2010_low_mrs3`, `pb_sun_2010_low_sms3`, `pb_sun_2010_low_sms4`, `pb_sun_2010_low_vn`, `pb_sun_2010_sms3_1`, `pb_sun_2010_sms4_1`, `pd_sun_2010_low_mrs3`, `pd_sun_2010_low_sms3`, `pd_sun_2010_low_sms4`, `pd_sun_2010_low_vn`, `pd_sun_2010_mrs3_1`, `pd_sun_2010_sms3_1`, `pd_sun_2010_sms4_1`, `pt_sun_2010_low_bn`, `pt_sun_2010_low_mrs3`, `pt_sun_2010_low_sms3`, `pt_sun_2010_low_sms4`, `pt_sun_2010_low_vn`, `pt_sun_2010_mrs3_1`, `pt_sun_2010_sms3_1`, `pt_sun_2010_sms4_1`, `rb_sun_2010_sms3_1`, `rb_sun_2010_sms4_1`, `rbbr_sun_2010_sms3_1`, `rbbr_sun_2010_sms4_1`, `rbcl_sun_2010_sms3_1`, `rbcl_sun_2010_sms4_1`, `rbf_sun_2010_sms3_1`, `rbf_sun_2010_sms4_1`, `rbi_sun_2010_sms3_1`, `rbi_sun_2010_sms4_1`, `sn_sun_2010_low_mrs3`, `sn_sun_2010_low_sms3`, `sn_sun_2010_low_sms4`, `sn_sun_2010_low_vn`, `sn_sun_2010_sms3_1`, `sn_sun_2010_sms4_1`, `solid_h2_sun_2010_mrs3_1`, `solid_h2_sun_2010_sms3_1`, `solid_h2_sun_2010_sms4_1`, `ta_sun_2010_low_bn`, `ta_sun_2010_low_mrs3`, `ta_sun_2010_low_sms3`, `ta_sun_2010_low_sms4`, `ta_sun_2010_low_vn`, `ta_sun_2010_mrs3_1`, `ta_sun_2010_sms3_1`, `ta_sun_2010_sms4_1`, `th_sun_2010_low_mrs3`, `th_sun_2010_low_sms3`, `th_sun_2010_low_sms4`, `th_sun_2010_low_vn`, `th_sun_2010_sms3_1`, `th_sun_2010_sms4_1`, `ti_sun_2010_low_mrs3`, `ti_sun_2010_low_sms3`, `ti_sun_2010_low_sms4`, `ti_sun_2010_low_vn`, `ti_sun_2010_mrs3_1`, `ti_sun_2010_sms3_1`, `ti_sun_2010_sms4_1`, `tl_sun_2010_low_mrs3`, `tl_sun_2010_low_sms3`, `tl_sun_2010_low_sms4`, `tl_sun_2010_low_vn`, `tl_sun_2010_sms3_1`, `tl_sun_2010_sms4_1`, `v_sun_2010_low_mrs3`, `v_sun_2010_low_sms3`, `v_sun_2010_low_sms4`, `v_sun_2010_low_vn`, `v_sun_2010_sms3_1`, `v_sun_2010_sms4_1`, `w_sun_2010_low_bn`, `w_sun_2010_low_mrs3`, `w_sun_2010_low_sms3`, `w_sun_2010_low_sms4`, `w_sun_2010_low_vn`, `w_sun_2010_mrs3_1`, `w_sun_2010_sms3_1`, `w_sun_2010_sms4_1`, `zn_sun_2010_low_mrs3`, `zn_sun_2010_low_sms3`, `zn_sun_2010_low_sms4`, `zn_sun_2010_low_vn`, `zn_sun_2010_mrs3_1`, `zn_sun_2010_sms3_1`, `zn_sun_2010_sms4_1`, `zr_sun_2010_low_mrs3`, `zr_sun_2010_low_sms3`, `zr_sun_2010_low_sms4`, `zr_sun_2010_low_vn`, `zr_sun_2010_sms3_1`, `zr_sun_2010_sms4_1` | Published coefficients are executable; this import does not claim a new refit of row-level observations. Source coefficients are complete and executable, but the fitted row-level compression observations and regression weights are not republished. The published coefficients are executable, but an independent refit was not attempted: no reusable license was identified for the all-rights-reserved handbook table, its OCR is not reliable enough to transcribe silently, and the original experiment-level rows and numerical implementation are not deposited. This record preserves Sun et al.'s published coefficients and is not itself refitted. A separate dataset-backed Peritheos record, cu_sun_2010_low_vn_refit, uses the manually verified handbook Cu cells; experiment-level rows and the source authors' numerical implementation remain unavailable. Source-author coefficients are complete and executable, but the exact fitted isotherm rows and numerical implementation are not republished. This catalog action is a transcription, not a Peritheos refit. |
| [Sun et al. (2016)](https://doi.org/10.1002/2016jb013062) | `ca_perovskite_sun_2016_bm3_3` | The published thermal-EOS coefficients are transcribed directly. The article's P-V-T table is not redistributed because no open table-data license was identified. |
| [Sun et al. (2019)](https://doi.org/10.1029/2018gl081421) | `fesio3_liquid_sun_2019_2500k_bm4_1` | The bundled Table 1 grid mixes liquid and nonliquid simulations. Figure 1, rather than the numerical table, identifies the liquid states used by the source fit, so the table is a checkpoint resource and not an asserted exact regression input. |
| [Sun et al. (2022)](https://doi.org/10.2138/am-2021-7913) | `ca_perovskite_tetragonal_sun_2022_bm3_1` | The published fixed-derivative BM3 coefficients are transcribed directly. The article's P-V table is not redistributed because no open table-data license was identified. |
| [Tange et al. (2009)](https://doi.org/10.1029/2008jb005813) | `mgo_b1_tange_2009_vinet` | This is a unified least-squares analysis of previously published pressure-scale-free thermal, elastic, and shock datasets. It reports optimized MgO EOS parameters and residuals, but no new row-level experimental observations. |
| [Wu et al. (2013)](https://doi.org/10.7498/aps.62.049101) | `bridgmanite_wu_2013_gga_bm3`, `mg075fe025sio3_bridgmanite_wu_2013_gga_bm3` | Pressure-dependent calculations are plotted, but the numerical P-V grid is not tabulated. |
| [Xiao et al. (2013)](https://doi.org/10.2138/am.2013.4470) | `srsio3_6h_xiao_2013_gga_bm2_1`, `srsio3_cubic_xiao_2013_gga_bm2_2` | Calculated states are plotted but not tabulated; no graphical pseudo-precision was introduced. |
| [Yang et al. (2015)](https://doi.org/10.1038/srep17188) | `mg092fe008o_yang_2015_hs_bm3_reference` | No machine-readable P-V table is published. Analytical BM3 checkpoints and inverse-volume round trips independently verify executable transcription of the source coefficients. |
| [Zhang and Bukowinski (1991)](https://doi.org/10.1103/physrevb.44.2495) | `mgo_b1_zhang_bukowinski_1991_mpib_bm3`, `mgo_b2_zhang_bukowinski_1991_mpib_bm3`, `stishovite_zhang_bukowinski_1991_mpib_bm3` | Calculated states are plotted but not tabulated. |
| [Zhang and Wentzcovitch (2022)](https://doi.org/10.1103/physrevb.106.054103) | `bridgmanite_zhang_wentzcovitch_2022_phq_lda_300k_bm3`, `bridgmanite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3`, `mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_lda_300k_bm3`, `mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3` | Five F(V) states per temperature are described and the curves are plotted, but numerical free-energy rows, weights, residuals, and covariance are not deposited. Five F(V) states per temperature are described and plotted but not numerically deposited. |
| [Zhao et al. (1997)](https://doi.org/10.1029/96gl03769) | `naalsi2o6_zhao_1997_bm3_1` | Only 1 observation(s) lie at the reference temperature for 1 free isothermal coefficients; the other rows require a thermal relation that this record does not represent. |
| [Zhu et al. (2025)](https://doi.org/10.22541/essoar.176236186.65259830/v1) | `gold_zhu_2025_vinet_300k`, `mgo_zhu_2025_vinet_300k`, `platinum_zhu_2025_vinet_300k` | The record cites an external repository or source dataset that is not bundled as a row-level material dataset. |

## Complete investigated-paper register

This is the exhaustive paper-level index. `Bundled` means numerical primary
rows are stored; `bundled indirect` means source coefficient or derived-output
tables are stored but do not form a common direct-fit observation matrix;
`plot only` means observations were digitized;
`parameterization only` means only the published equation/coefficients can
be checked. Record-level links, fit metrics, and evidence locations are in
the primary-source and refit ledgers.

| Paper | Final disposition | Catalog records | Record-level results | Primary-data form |
|---|---|---:|---|---|
| [A modified Anderson–Grüneisen model for the pressure dependence of thermal expansivity (2019)](https://doi.org/10.1139/cjp-2019-0326) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Ahrens and Gaffney (1975), Mg-Fe silicate shock interpretation](https://doi.org/10.1111/j.1365-246X.1975.tb06463.x) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [AIP Handbook Cu low-pressure Vinet m=8 refit](https://doi.org/10.1515/zna-2010-1-202) | Reproduced | 1 | 1 similar |  |
| [Akaogi et al. (1996), natural-basalt majoritic garnet](https://doi.org/10.1029/96GL03027) | Direct refit unavailable | 0 | no production record | investigation evidence only |
| [Akber-Knutson et al. (2002)](https://doi.org/10.1029/2001gl013523) | Direct refit unavailable | 4 | 4 direct refit unavailable | 4 parameterization only |
| [Akber-Knutson et al. (2005)](https://doi.org/10.1029/2005gl023192) | Direct refit unavailable | 10 | 10 direct refit unavailable | 10 parameterization only |
| [Akins et al. (2004)](https://doi.org/10.1029/2004gl020237) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [ANALYSIS OF THERMAL EXPANSIVITY OF SOLIDS UNDER HIGH PRESSURES (2012)](https://doi.org/10.1142/s0217984912501461) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Anderson and Zou (1990), MgO thermodynamic functions](https://doi.org/10.1063/1.555873) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Anderson et al. (1989)](https://doi.org/10.1063/1.342969) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled indirect |
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
| [Belmonte (2017)](https://doi.org/10.3390/min7100183) | Reproduced | 1 | 1 similar | 1 bundled |
| [Benedict et al. (2014)](https://doi.org/10.1103/physrevb.89.224109) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Besson et al. (1994)](https://doi.org/10.1103/physrevb.49.12540) | Reproduced | 1 | 1 similar | 1 bundled |
| [Bezacier et al. (2014)](https://doi.org/10.1063/1.4894421) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Boffa Ballaran et al. (2007)](https://doi.org/10.2138/am.2007.2715) | Reproduced | 2 | 2 parity | 2 bundled |
| [Braithwaite (2002), high-pressure solid/melt representations](https://doi.org/10.1063/1.1483512) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Bykova et al. (2018)](https://doi.org/10.1038/s41467-018-07265-z) | Partly reproduced | 3 | 1 similar; 2 direct refit unavailable | 2 bundled; 1 plot only/digitized |
| [Campbell and Heinz (1991)](https://doi.org/10.1016/0022-3697(91)90181-x) | Reproduced | 1 | 1 parity | 1 bundled |
| [Campbell and Heinz (1993)](https://doi.org/10.1016/0022-3697(93)90106-2) | Reproduced | 1 | 1 parity | 1 bundled |
| [Campbell and Heinz (1994)](https://doi.org/10.1029/94jb00127) | Reproduced | 2 | 2 parity | 2 bundled |
| [Caracas and Cohen (2005), MgSiO3-FeSiO3-Al2O3 pv/ppv chemistry](https://doi.org/10.1029/2005GL023164) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Caracas et al. (2005)](https://doi.org/10.1029/2004gl022144) | Direct refit unavailable | 18 | 18 direct refit unavailable | 18 theoretical parameterization only |
| [Chantel et al. (2012)](https://doi.org/10.1029/2012gl053075) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Chen et al. (2018)](https://doi.org/10.2138/am-2018-6087) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Chen et al. (2024), stishovite velocities](https://doi.org/10.1029/2023GL107700) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Chidester et al. (2021)](https://doi.org/10.1103/physrevb.104.094107) | Reproduced | 1 | 1 parity | 1 bundled |
| [Chizmeshya et al. (1996)](https://doi.org/10.1029/96gl02624) | Direct refit unavailable | 5 | 5 direct refit unavailable | 5 theoretical parameterization only |
| [Clendenen and Drickamer (1966)](https://doi.org/10.1063/1.1726610) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Cohen and Lin (2014)](https://doi.org/10.1103/physrevb.90.140102) | Direct refit unavailable | 3 | 3 direct refit unavailable |  |
| [Comparative Analysis of Grüneisen Parameters for Selected Geophysical Minerals Using Advanced Equations of State (2024)](https://doi.org/10.69626/sea.2024.0152) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Compositional constraints on the equation of state and thermal properties of the lower mantle (2001)](https://doi.org/10.1046/j.1365-246x.2001.00437.x) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Constraints on lower mantle composition and temperature from density and bulk sound velocity profiles (1990)](https://doi.org/10.1029/gl017i008p01153) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Correa et al. (2008)](https://doi.org/10.1103/physrevb.78.024101) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Crichton et al. (2002)](https://doi.org/10.2138/am-2002-2-316) | Reproduced | 1 | 1 parity | 1 bundled |
| [Criniti et al. (2021)](https://doi.org/10.1029/2020jb020967) | Reproduced | 1 | 1 parity | 1 bundled |
| [Criniti et al. (2023)](https://doi.org/10.2138/am-2022-8559) | Reproduced | 2 | 2 similar | 2 bundled |
| [Cynn and Yoo (1999)](https://doi.org/10.1103/physrevb.59.8526) | Reproduced | 1 | 1 parity | 1 bundled |
| [Daniel et al. (2004)](https://doi.org/10.1029/2004gl020213) | Reproduced | 2 | 2 parity | 2 bundled |
| [data from JCPDS and Levien and Prewitt, 1981](https://msaweb.org/AmMin/AM66/AM66_324.pdf) | Reproduced | 1 | 1 similar | 1 bundled |
| [Datchi et al. (2007)](https://doi.org/10.1103/physrevb.75.214104) | Partly reproduced | 3 | 1 parity; 1 similar; 1 direct refit unavailable | 2 bundled; 1 parameterization only |
| [Delta archive (2023)](https://doi.org/10.24435/materialscloud:5e-mv) | Direct refit unavailable | 42 | 42 direct refit unavailable | 42 parameterization only |
| [Deng et al. (2006)](https://doi.org/10.1088/0256-307x/23/8/101) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Dewaele (2019)](https://doi.org/10.3390/min9110684) | Mixed: reproduced and discrepant records | 31 | 1 parity; 4 parity not achieved; 26 direct refit unavailable | 3 bundled; 26 parameterization only |
| [Dewaele and Torrent (2013)](https://doi.org/10.1103/physrevb.88.064107) | Reproduced | 1 | 1 parity | 1 bundled |
| [Dewaele et al. (2000)](https://doi.org/10.1029/1999jb900364) | Reproduced | 5 | 5 parity | 1 bundled |
| [Dewaele et al. (2004)](https://doi.org/10.1103/physrevb.70.094112) | Reproduced | 12 | 12 parity | 12 bundled |
| [Dewaele et al. (2005), high-pressure metrology abstract](https://doi.org/10.1107/S0108767305096972) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.77.094106) | Partly reproduced | 3 | 1 parity; 2 direct refit unavailable | 3 bundled |
| [Dewaele et al. (2008)](https://doi.org/10.1103/physrevb.78.104102) | Reproduced | 2 | 2 similar | 2 bundled |
| [Dewaele et al. (2012)](https://doi.org/10.1103/physrevb.85.214105) | Reproduced | 4 | 4 similar | 4 bundled |
| [Dewaele et al. (2015)](https://doi.org/10.1103/physrevb.91.134108) | Reproduced | 2 | 2 parity | 2 bundled |
| [Dobrosavljevic et al. (2019)](https://doi.org/10.3390/min9120762) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Dorfman et al. (2012)](https://doi.org/10.1029/2012jb009292) | Coefficient parity not achieved | 6 | 6 parity not achieved | 6 bundled |
| [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/physrevb.75.024115) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 partial primary rows bundled |
| [Dorogokupets et al. (2015)](https://doi.org/10.1016/j.rgg.2015.01.011) | Direct refit unavailable | 3 | 3 direct refit unavailable | 3 parameterization only |
| [Driver et al. (2010)](https://doi.org/10.1073/pnas.0912130107) | Direct refit unavailable | 3 | 3 direct refit unavailable |  |
| [Dubrovinsky et al. (2002)](https://doi.org/10.1080/08957950212807) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Duffy and Ahrens (1995)](https://doi.org/10.1029/94jb02065) | Reproduced | 1 | 1 parity | 1 bundled |
| [Duffy et al. (1995)](https://doi.org/10.1103/physrevlett.74.1371) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Effect of Pressure on the Composition of the Lower Mantle End Member Fe x O (1993)](https://doi.org/10.1126/science.259.5091.66) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Elastic properties of Fe-bearing Akimotoite at mantle conditions: Implications for composition and temperature in lower mantle transition zone (2022)](https://doi.org/10.1016/j.fmre.2021.12.013) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Equation of State, Phase Stability of (Mg0.92, Fe0.08)SiO3 Perovskite from Shock Wave Study and Its Geophysical Implications (2004)](https://doi.org/10.1063/1.1780510) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Exploring the High-Pressure Equation of State in Earth’s Mantle with a Focus on the MgSiO3−MgO System (2025)](https://doi.org/10.15407/mfint.47.06.0601) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [EXTREME COMPRESSION BEHAVIOUR OF SOLIDS BASED ON THE ROY-ROY INVERTED EQUATION OF STATE (2008)](https://doi.org/10.1142/s0217979208038910) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Faust and Knittle (1994), natural chondrodite](https://doi.org/10.1029/94GL01592) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Fedotenko et al. (2020)](https://doi.org/10.1016/j.jallcom.2020.156179) | Reproduced | 2 | 2 parity | 2 bundled |
| [Fei et al. (2000)](https://doi.org/10.2138/am-2000-11-1229) | Reproduced | 1 | 1 similar | 1 plot only/digitized |
| [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) | Reproduced | 3 | 2 parity; 1 similar | 3 plot only/digitized |
| [Ferre et al. (2009), dislocations in CaSiO3 perovskite](https://doi.org/10.2138/am.2009.3003) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Finding the isentropic density of perovskite: Implications for iron concentration in the lower mantle (1997)](https://doi.org/10.1029/96gl03951) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Finkelstein et al. (2014)](https://doi.org/10.2138/am.2014.4526) | Reproduced | 1 | 1 parity | 1 bundled |
| [Finkelstein et al. (2017)](https://doi.org/10.2138/am-2017-5966) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Fiquet et al. (2000)](https://doi.org/10.1029/1999gl008397) | Reproduced | 1 | 1 parity | 1 bundled |
| [Fortes (2019)](https://epubs.stfc.ac.uk/manifestation/40740885/RAL-TR-2019-002.pdf) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Frank et al. (2004)](https://doi.org/10.1016/j.gca.2003.12.007) | Reproduced | 1 | 1 parity | 1 bundled |
| [Fratanduono et al. (2021)](https://doi.org/10.1126/science.abh0364) | Reproduced | 1 | 1 parity | 1 bundled |
| [Frost et al. (2023)](https://doi.org/10.1063/5.0161038) | Reproduced | 2 | 2 similar | 2 bundled |
| [Fu et al. (2023)](https://doi.org/10.2138/am-2022-8435) | Partly reproduced | 3 | 2 similar; 1 direct refit unavailable | 2 bundled |
| [Fu et al. (2024)](https://doi.org/10.2138/am-2023-8969) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Fuchizaki (2019), updated MgO melting curve](https://doi.org/10.7566/JPSJ.88.065003) | Direct refit unavailable | 0 | no production record | investigation evidence only |
| [Fujihisa and Takemura (1996)](https://doi.org/10.1103/physrevb.54.5) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Funamori et al. (1996), MgSiO3 perovskite thermoelasticity](https://doi.org/10.1029/95JB03732) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Funamori et al. (1998)](https://doi.org/10.1029/98jb01575) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 bundled |
| [Fundamental thermodynamic relations and silicate melting with implications for the constitution of D″ (1990)](https://doi.org/10.1029/jb095ib12p19311) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Gerward et al. (2005)](https://doi.org/10.1016/j.jallcom.2005.04.008) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Ghosh and Karki (2016)](https://doi.org/10.1038/srep37269) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) | Coefficient parity not achieved | 2 | 1 parity not achieved; 1 direct refit unavailable | 2 bundled |
| [Gong et al. (2004)](https://doi.org/10.1029/2003gl019132) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Guigue et al. (2020)](https://doi.org/10.1063/1.5138697) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 plot only/digitized |
| [Haavik et al. (2000), defect-clustered wuestite](https://doi.org/10.1039/B006026G) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Haines et al. (2001)](https://doi.org/10.1088/0953-8984/13/11/303) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Hama and Suito (1996)](https://doi.org/10.1088/0953-8984/8/1/008) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Hama and Suito (1998), Mg-Fe bridgmanite](https://doi.org/10.1029/97JB03672) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Hama and Suito (2000), Vinet-Debye mineral parameterizations](https://doi.org/10.2138/am-2000-2-310) | Direct refit unavailable | 0 | no production record | investigation evidence only |
| [Hamahata et al. (2000)](https://doi.org/10.2465/jmps.95.236) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Hanfland et al. (1989)](https://doi.org/10.1103/physrevb.39.12598) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Hanfland et al. (1999)](https://doi.org/10.1016/s0038-1098(99)00322-1) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hanna et al. (2011)](https://doi.org/10.1063/1.3644969) | Reproduced | 2 | 1 parity; 1 similar | 2 plot only/digitized |
| [Hazen and Finger (1979)](https://msaweb.org/AmMin/AM64/AM64_196.pdf) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hazen and Finger (1981)](https://doi.org/10.1016/0022-3697(81)90074-3) | Reproduced | 2 | 2 similar | 2 bundled |
| [Heinz and Jeanloz (1984)](https://doi.org/10.1103/physrevb.30.6045) | Reproduced | 1 | 1 parity | 1 bundled |
| [Hemley et al. (1989)](https://doi.org/10.1103/physrevb.39.11820) | Reproduced | 1 | 1 parity | 1 bundled |
| [High pressure and high temperature in situ X‐ray observation of MgSiO3 Perovskite under lower mantle conditions (1993)](https://doi.org/10.1029/92gl02960) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Hirose et al. (2005)](https://doi.org/10.2138/am.2005.1702) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Holland et al. (2013)](https://doi.org/10.1093/petrology/egt035) | Direct refit unavailable | 7 | 7 direct refit unavailable |  |
| [Holmes et al. (1989)](https://doi.org/10.1063/1.344177) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Inbar and Cohen (1995), MgO thermal properties](https://doi.org/10.1029/95GL01086) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Indoor seismology by probing the Earth's interior by using sound velocity measurements at high pressures and temperatures (2007)](https://doi.org/10.1073/pnas.0608609104) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Irifune et al. (2002)](https://doi.org/10.1007/s00269-002-0275-1) | Reproduced | 1 | 1 parity | 1 bundled |
| [Isaak et al. (1990), calculated MgO properties](https://doi.org/10.1029/JB095iB05p07055) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Ismailova et al. (2016)](https://doi.org/10.1126/sciadv.1600427) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Ito et al. (2011), Kawai apparatus](https://doi.org/10.4131/jshpreview.21.272) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Jackson et al. (1990), wuestite elasticity](https://doi.org/10.1029/JB095iB13p21671) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Jackson et al. (2006), ferropericlase elasticity](https://doi.org/10.1029/2005JB004052) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Jacobs and Oonk (2000), GGK MgO EOS](https://doi.org/10.1039/A910247G) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Jacobsen et al. (2002)](https://doi.org/10.1029/2001jb000490) | Coefficient parity not achieved | 3 | 3 parity not achieved | 3 bundled |
| [Jacobsen et al. (2005)](https://doi.org/10.1107/s0909049505022326) | Mixed: reproduced and discrepant records | 3 | 1 parity; 2 parity not achieved | 3 bundled |
| [Jacobsen et al. (2008)](https://doi.org/10.2138/am.2008.2988) | Reproduced | 2 | 2 parity | 2 bundled |
| [Karki and Crain (1998)](https://doi.org/10.1029/98gl51952) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Karki and Stixrude (1999), lower-mantle elastic moduli](https://doi.org/10.1029/1999JB900069) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Karki and Wentzcovitch (2002)](https://doi.org/10.1029/2001jb000702) | Direct refit unavailable | 4 | 4 direct refit unavailable |  |
| [Karki et al. (1997)](https://doi.org/10.2138/am-1997-1-207) | Direct refit unavailable | 2 | 2 direct refit unavailable |  |
| [Karki et al. (2000), theoretical MgSiO3 akimotoite](https://doi.org/10.2138/am-2000-2-309) | Direct refit unavailable | 0 | no production record | investigation evidence only |
| [Karki et al. (2001), MgSiO3 thermodynamic derivatives](https://doi.org/10.1029/2001GL012910) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Katsura et al. (2004), Mg2SiO4 ringwoodite](https://doi.org/10.1029/2004JB003094) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Katsura et al. (2009)](https://doi.org/10.1029/2008gl035658) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Katsura et al. (2009)](https://doi.org/10.1029/2009gl038107) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Kawai and Tsuchiya (2012)](https://doi.org/10.2138/am.2012.3915) | Reproduced | 2 | 2 similar | 2 plot only/digitized |
| [Kawai and Tsuchiya (2014)](https://doi.org/10.1002/2013jb010905) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Kawai and Tsuchiya (2015), CaSiO3 thermoelasticity](https://doi.org/10.1002/2015GL063446) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Kholiya et al. (2014), MgO EOS comparison](https://doi.org/10.1155/2014/289353) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Kiefer et al. (2002)](https://doi.org/10.1029/2002gl014683) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 theoretical parameterization only |
| [Knittle and Jeanloz (1987)](https://doi.org/10.1126/science.235.4789.668) | Reproduced | 1 | 1 similar | 1 bundled |
| [Knittle and Jeanloz (1991), MgSiO3 perovskite transition and thermal expansion](https://doi.org/10.1126/science.251.4992.410) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Knorr et al. (2003)](https://doi.org/10.1140/epjb/e2003-00034-6) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Koemets et al. (2023)](https://doi.org/10.3389/fchem.2023.1258389) | Reproduced | 1 | 1 parity | 1 bundled |
| [Kubo et al. (2000)](https://doi.org/10.2183/pjab.76.103) | Reproduced | 2 | 2 parity | 2 bundled |
| [Kubo et al. (2006)](https://doi.org/10.1029/2006gl025686) | Reproduced | 2 | 2 parity | 2 bundled |
| [Kumari and Dass (1990), EOS applied to 50 solids II](https://doi.org/10.1088/0953-8984/2/39/003) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Lakshtanov et al. (2007), Al-H stishovite elasticity](https://doi.org/10.2138/am.2007.2294) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Le Godec et al. (2014)](https://doi.org/10.3103/s1063457614010092) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Lee and Wan (2008)](https://doi.org/10.1103/physrevb.78.224103) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 theoretical parameterization only |
| [Lejaeghere et al. (2016)](https://doi.org/10.1126/science.aad3000) | Direct refit unavailable | 50 | 50 direct refit unavailable | 50 parameterization only |
| [Leonov et al. (2017)](https://doi.org/10.1103/physrevb.96.075136) | Direct refit unavailable | 8 | 8 direct refit unavailable |  |
| [Li and Zeng (2009)](https://doi.org/10.1142/s0129183109014242) | Direct refit unavailable | 6 | 6 direct refit unavailable |  |
| [Li et al. (2006)](https://doi.org/10.1029/2005jb004251) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Li-and-Zhang-data elastic refit (2010)](https://doi.org/10.1016/j.pnsc.2009.09.002) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Liu (2008), MgO bulk modulus method](https://doi.org/10.1515/zna-2008-1-209) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Liu (2011), two-parameter MgO EOS analysis](https://doi.org/10.1139/p11-040) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Liu and Bassett (1973)](https://doi.org/10.1029/jb078i035p08470) | Reproduced | 1 | 1 parity | 1 bundled |
| [Liu et al. (2007)](https://doi.org/10.1088/0953-8984/19/24/246103) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Liu et al. (2010)](https://doi.org/10.1142/s0217984910022391) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Liu et al. (2011)](https://doi.org/10.1088/1674-0068/24/06/703-710) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Lower-mantle MgSiO3 thermoelastic parameter review (1996)](https://doi.org/10.1098/rsta.1996.0053) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Luo et al. (2023)](https://doi.org/10.1103/physrevb.107.134116) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Lv et al. (2020)](https://doi.org/10.2138/am-2020-7279) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Magad-Weiss et al. (2021)](https://doi.org/10.1103/physrevb.103.014101) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Mao et al. (1974)](https://doi.org/10.1029/jb079i008p01165) | Reproduced | 1 | 1 parity | 1 bundled |
| [Mao et al. (1989)](https://doi.org/10.1029/jb094ib12p17889) | Reproduced | 1 | 1 parity | 1 bundled |
| [Mao et al. (1990)](https://doi.org/10.1029/jb095ib13p21737) | Reproduced | 1 | 1 similar | 1 bundled |
| [Mao et al. (1991)](https://doi.org/10.1029/91jb00176) | Direct refit unavailable | 3 | 3 direct refit unavailable | 3 parameterization only |
| [Mao et al. (2011)](https://doi.org/10.1029/2011gl049519) | Reproduced | 3 | 2 parity; 1 similar | 3 plot only/digitized |
| [Mao et al. (2011)](https://doi.org/10.1029/2011gl049915) | Reproduced | 2 | 2 parity |  |
| [Mao et al. (2015)](https://doi.org/10.1002/2015gl064400) | Reproduced | 1 | 1 similar | 1 bundled |
| [Marcondes et al. (2020)](https://doi.org/10.1103/physrevb.102.104112) | Direct refit unavailable | 8 | 8 direct refit unavailable |  |
| [Martin et al. (2007)](https://doi.org/10.2138/am.2007.2473) | Reproduced | 1 | 1 parity | 1 bundled |
| [Martinez et al. (1996)](https://doi.org/10.2138/am-1996-5-608) | Reproduced | 1 | 1 parity | 1 bundled |
| [Matsui (1993), molecular dynamics of high-pressure silicates](https://doi.org/10.5940/jcrsj.35.190) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Matsui (2002), MgSiO3-Al2O3 simulations](https://doi.org/10.2465/jmps.97.13) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Matsui et al. (1994), MgSiO3 simulations](https://doi.org/10.1029/94GL01370) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Matsui et al. (2000), MgO pressure standard](https://doi.org/10.2138/am-2000-2-308) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Matsui et al. (2012)](https://doi.org/10.2138/am.2012.3937) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [McCarthy and Harrison (1994), MgO bulk properties](https://doi.org/10.1103/PhysRevB.49.8574) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [McHardy et al. (2026)](https://doi.org/10.1103/zp3m-kjpc) | Reproduced | 3 | 3 parity | 3 bundled |
| [Meng et al. (1994)](https://doi.org/10.1007/bf00203299) | Reproduced | 1 | 1 parity | 1 bundled |
| [Metsue and Tsuchiya (2012)](https://doi.org/10.1111/j.1365-246x.2012.05511.x) | Direct refit unavailable | 7 | 7 direct refit unavailable |  |
| [MgSiO3 elasticity temperature-derivative study (2004)](https://doi.org/10.1029/2003GL018762) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [MgSiO3 interaction-potential study (2007)](https://doi.org/10.1088/1674-0068/20/05/547-551) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [MgSiO3 molecular-dynamics compression study (2006)](https://doi.org/10.1360/CJCP2006.19(4).311.4) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Milani et al. (2015)](https://doi.org/10.1016/j.lithos.2015.03.017) | Reproduced | 2 | 2 parity | 2 bundled |
| [Miozzi et al. (2018)](https://doi.org/10.1029/2018je005582) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Miyajima et al. (2025), electron diffraction of a dense hydrous magnesium silicate](https://doi.org/10.1029/2025GL115280) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Mookherjee et al. (2015)](https://doi.org/10.2138/am-2015-5312) | Partly reproduced | 3 | 1 parity; 1 similar; 1 direct refit unavailable |  |
| [Mookherjee et al. (2019)](https://doi.org/10.2138/am-2019-6694) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Mookherjee et al. (2019), phase Egg accepted-manuscript alias](https://doi.org/10.2138/am-2018-6694) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Morishima et al. (1994), fixed-pressure CaSiO3 thermal expansion](https://doi.org/10.1029/94GL00844) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Mosenfelder et al. (2009)](https://doi.org/10.1029/2008jb005900) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Muhammad et al. (2024)](https://doi.org/10.1039/d4nr00093e) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Muñoz and Kunc (1993)](https://doi.org/10.1088/0953-8984/5/33/010) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 theoretical parameterization only |
| [Myhill (2022), anisotropic high-P-T EOS](https://doi.org/10.1093/gji/ggac180) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Noguchi et al. (1999)](https://doi.org/10.1016/s0022-3697(98)00296-0) | Partly reproduced | 2 | 1 parity; 1 direct refit unavailable | 2 bundled |
| [Noguchi et al. (2013)](https://doi.org/10.1007/s00269-012-0549-1) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Oganov and Dorogokupets (2003)](https://doi.org/10.1103/physrevb.67.224110) | Partly reproduced | 11 | 2 similar; 9 direct refit unavailable | 9 theoretical parameterization only |
| [Ono (2013)](https://doi.org/10.3390/e15104300) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Ono et al. (2000)](https://doi.org/10.1007/s002690000108) | Coefficient parity not achieved | 2 | 2 parity not achieved | 2 bundled |
| [Ono et al. (2006)](https://doi.org/10.2138/am.2006.2347) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ono et al. (2006)](https://doi.org/10.1007/s00269-006-0068-z) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ono et al. (2006)](https://doi.org/10.2138/am.2006.2118) | Reproduced | 3 | 3 parity | 3 bundled |
| [Pamato et al. (2016), NAL elasticity](https://doi.org/10.1002/2016JB013136) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Pepin et al. (2014)](https://doi.org/10.1103/physrevlett.113.265504) | Reproduced | 2 | 2 parity | 2 plot only/digitized |
| [Performance of an ab initio equation of state for magnesium oxide (2004)](https://doi.org/10.1088/0953-8984/16/30/006) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Prescher et al. (2015)](https://doi.org/10.1038/ngeo2370) | Reproduced | 1 | 1 similar | 1 bundled |
| [Pressure Dependence of Interatomic Separation and Thermal Expansivity for Alkali Halides and Periclase (MgO) (2009)](https://doi.org/10.12693/aphyspola.115.709) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Pressure-Induced Magnetization in FeO: Evidence from Elasticity and Mössbauer Spectroscopy (2004)](https://doi.org/10.1103/physrevlett.93.215502) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Pressure-volume equation of state of the high-pressureB2phase of NaCl (2002)](https://doi.org/10.1103/physrevb.65.104114) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Qin et al. (2023)](https://doi.org/10.2138/am-2022-8432) | Reproduced | 2 | 2 parity | 2 bundled |
| [Redfern and Angel (1999)](https://doi.org/10.1007/s004100050471) | Reproduced | 1 | 1 parity | 1 bundled |
| [Redfern et al. (1993)](https://doi.org/10.1029/93gl02507) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 parameterization only |
| [Refinement of enthalpy measurement of MgSiO3 perovskite and negative pressure‐temperature slopes for Perovskite‐forming reactions (1993)](https://doi.org/10.1029/93gl01265) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Reichmann et al. (2008), MgO elasticity](https://doi.org/10.2138/am.2008.2717) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Research progresses on the structure, phase transition and physical properties of MgSiO&amp;lt;sub&amp;gt;3&amp;lt;/sub&amp;gt; under high temperature and high pressure (2025)](https://doi.org/10.3724/j.issn.1007-2802.20240151) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Reynard et al. (1996)](https://doi.org/10.2138/am-1996-1-206) | Reproduced | 2 | 2 parity | 2 bundled |
| [Richet et al. (1988)](https://doi.org/10.1029/jb093ib12p15279) | Reproduced | 2 | 2 parity | 2 bundled |
| [Richet et al. (1989)](https://doi.org/10.1029/jb094ib03p03037) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ricolleau et al. (2009)](https://doi.org/10.1029/2008gl036759) | Direct refit unavailable | 4 | 4 direct refit unavailable |  |
| [Rodrigo-Ramon et al. (2024)](https://doi.org/10.1038/s41598-024-78006-0) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ross (1997)](https://doi.org/10.2138/am-1997-7-805) | Reproduced | 1 | 1 parity | 1 bundled |
| [Ross and Angel (1999)](https://doi.org/10.2138/am-1999-0309) | Reproduced | 2 | 2 parity | 2 bundled |
| [Sagatova et al. (2021)](https://doi.org/10.1134/s0016702921080073) | Direct refit unavailable | 7 | 7 direct refit unavailable |  |
| [Sakai et al. (2016)](https://doi.org/10.1038/srep22652) | Reproduced | 1 | 1 parity | 1 bundled |
| [Sakai et al. (2025)](https://doi.org/10.1038/s43246-025-00792-5) | Direct refit unavailable | 9 | 9 direct refit unavailable |  |
| [Sato and Jeanloz (1981)](https://doi.org/10.1029/jb086ib12p11773) | Reproduced | 1 | 1 parity | 1 bundled |
| [Satta et al. (2025), Fe-bearing delta-AlOOH elasticity](https://doi.org/10.1007/s00269-025-01319-7) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Saxena et al. (1999)](https://doi.org/10.2138/am-1999-0303) | Reproduced | 1 | 1 similar | 1 bundled |
| [Schoelmerich et al. (2020)](https://doi.org/10.1038/s41598-020-66340-y) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Schouwink et al. (2011)](https://doi.org/10.2138/am.2011.3775) | Reproduced | 1 | 1 parity | 1 bundled |
| [Schulze et al. (2018)](https://doi.org/10.2138/am-2018-6562) | Reproduced | 1 | 1 parity | 1 bundled |
| [Scott et al. (2001)](https://doi.org/10.1029/2000gl012606) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Shen and Smith (2026)](https://doi.org/10.1103/fxgq-96sg) | Direct refit unavailable | 10 | 10 direct refit unavailable | 10 bundled |
| [Sherman (1993)](https://doi.org/10.1029/93jb02175) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Sherman et al. (1993), stishovite and modified-fluorite SiO2](https://doi.org/10.1029/93JB00783) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Shi et al. (2022)](https://doi.org/10.1029/2021jb023805) | Reproduced | 2 | 1 parity; 1 similar | 2 bundled |
| [Shieh et al. (2000)](https://doi.org/10.1016/s0012-821x(00)00033-9) | Reproduced | 2 | 2 parity | 2 bundled |
| [Shieh et al. (2002), MgSiO3 post-perovskite](https://doi.org/10.1103/PhysRevLett.89.255507) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Shieh et al. (2006)](https://doi.org/10.1073/pnas.0506811103) | Reproduced | 3 | 3 parity | 3 plot only/digitized |
| [Shim et al. (2000)](https://doi.org/10.1029/2000jb900183) | Reproduced | 1 | 1 similar | 1 bundled |
| [Shim et al. (2000)](https://doi.org/10.1016/s0031-9201(00)00154-0) | Reproduced | 1 | 1 similar | 1 bundled |
| [Shim et al. (2002)](https://doi.org/10.1029/2002gl016148) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Shukla et al. (2016), Fe3+- and Al-bearing bridgmanite](https://doi.org/10.1002/2016GL069332) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Siersch et al. (2021)](https://doi.org/10.1016/j.pepi.2021.106786) | Reproduced | 1 | 1 parity | 1 bundled |
| [Singh and Singh (2021), alkaline-earth oxide EOS formulation](https://doi.org/10.12693/aphyspola.140.131) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Sinogeikin and Bass (1999), MgO elasticity](https://doi.org/10.1103/PhysRevB.59.R14141) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Sinogeikin et al. (2004), MgSiO3 elasticity](https://doi.org/10.1029/2004GL019559) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Sokolova et al. (2013)](https://doi.org/10.1016/j.rgg.2013.01.005) | Direct refit unavailable | 11 | 11 direct refit unavailable | 11 parameterization only |
| [Sokolova et al. (2018), MgO-MgSiO3 spreadsheets](https://doi.org/10.1080/08957959.2018.1465056) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Sokolova et al. (2021), Ca-silicate EOS](https://doi.org/10.3390/min11030322) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Solomatova et al. (2016)](https://doi.org/10.2138/am-2016-5510) | Coefficient parity not achieved | 16 | 2 parity not achieved; 14 direct refit unavailable |  |
| [Somayazulu et al. (2023)](https://doi.org/10.1098/rsta.2022.0331) | Mixed: reproduced and discrepant records | 3 | 1 parity; 1 similar; 1 parity not achieved | 3 bundled |
| [Sound-velocity temperature coefficients of MgSiO3 perovskite (2000)](https://doi.org/10.1088/0256-307X/17/3/022) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Speziale et al. (2001)](https://doi.org/10.1029/2000jb900318) | Reproduced | 1 | 1 parity | 1 bundled |
| [Speziale et al. (2007)](https://doi.org/10.1029/2006jb004730) | Reproduced | 1 | 1 similar | 1 bundled |
| [Spin crossover and Mott—Hubbard transition under high pressure and high temperature in the low mantle of the Earth (2015)](https://doi.org/10.1088/1742-6596/653/1/012095) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Stinton et al. (2014)](https://doi.org/10.1103/physrevb.90.134105) | Reproduced | 2 | 2 parity | 2 plot only/digitized |
| [Stixrude et al. (1992), thermoelasticity and mantle stratification](https://doi.org/10.1126/science.257.5073.1099) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Sun et al. (2010)](https://doi.org/10.1515/zna-2010-1-202) | Direct refit unavailable | 220 | 220 direct refit unavailable | 220 parameterization only |
| [Sun et al. (2016)](https://doi.org/10.1002/2016jb013062) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Sun et al. (2019)](https://doi.org/10.1029/2018gl081421) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Sun et al. (2019)](https://doi.org/10.1029/2019jb017853) | Reproduced | 2 | 2 similar | 2 bundled |
| [Sun et al. (2022)](https://doi.org/10.2138/am-2021-7913) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Suzuki (2016)](https://doi.org/10.2465/jmps.160719c) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura (2004)](https://doi.org/10.1103/physrevb.70.012101) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura and Dewaele (2008)](https://doi.org/10.1103/physrevb.78.104119) | Reproduced | 1 | 1 parity | 1 bundled |
| [Takemura and Singh (2006)](https://doi.org/10.1103/physrevb.73.224119) | Reproduced | 1 | 1 parity | 1 bundled |
| [Tange et al. (2009)](https://doi.org/10.1029/2008jb005813) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 parameterization only |
| [Tange et al. (2012)](https://doi.org/10.1029/2011jb008988) | Reproduced | 2 | 2 similar | 2 bundled |
| [Taniguchi et al. (1995), Ca-silicate calculation models](https://doi.org/10.2465/minerj.17.290) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Tateno et al. (2019)](https://doi.org/10.2138/am-2019-6779) | Reproduced | 1 | 1 parity | 1 bundled |
| [The effect of temperature on the product of bulk modulus and volume thermal expansion coefficient, and its application to the thermal expansion of MgO and other minerals (2004)](https://doi.org/10.1002/pssb.200302047) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [The sound velocity of wüstite at high pressures: implications for low-velocity anomalies at the base of the lower mantle (2020)](https://doi.org/10.1186/s40645-020-00333-3) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [The texture of the post-perovskite phase controls the characteristics of the D” seismic discontinuity (2025)](https://doi.org/10.1038/s43247-025-02383-1) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Thermal pressure in MgO and MgSiO 3 perovskite at lower mantle conditions (2000)](https://doi.org/10.2138/am-2000-1013) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Thermo-chemical and thermo-physical properties of the high-pressure phase anhydrous B (Mg14Si5O24): An ab-initio all-electron investigation (2010)](https://doi.org/10.2138/am.2010.3368) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Thermodynamic estimation the compressibility of ferropericlase under high pressure (2016)](https://doi.org/10.1063/1.4967779) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Thermoelasticity of perovskite: An emerging consensus (1994)](https://doi.org/10.1029/94eo01093) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Thompson et al. (2017)](https://doi.org/10.1002/2017jb014168) | Coefficient parity not achieved | 1 | 1 parity not achieved | 1 bundled |
| [Tsuchiya and Kawamura (2001), B1 oxide elasticity](https://doi.org/10.1063/1.1371498) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Tsuchiya and Mookherjee (2015)](https://doi.org/10.1038/srep15534) | Reproduced | 1 | 1 parity | 1 bundled |
| [Tsuchiya et al. (2004), MgSiO3 phase transition](https://doi.org/10.1016/j.epsl.2004.05.017) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Tsuchiya et al. (2004), post-perovskite elasticity](https://doi.org/10.1029/2004GL020278) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Vanpeteghem et al. (2002)](https://doi.org/10.1029/2001gl014224) | Reproduced | 2 | 2 parity | 2 bundled |
| [Vanpeteghem et al. (2006)](https://doi.org/10.1029/2005gl024955) | Reproduced | 1 | 1 parity | 1 bundled |
| [Vijay (2024), generalized Rydberg-Vinet and Stacey thermoelasticity](https://doi.org/10.32908/hthp.v53.1503) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Vocadlo (1999)](https://doi.org/10.2138/am-1999-1017) | Reproduced | 1 | 1 similar | 1 bundled |
| [Walker et al. (2002)](https://doi.org/10.2138/am-2002-0701) | Reproduced | 1 | 1 similar | 1 bundled |
| [Wang and Weidner (1994)](https://doi.org/10.1029/94gl00976) | Reproduced | 1 | 1 parity | 1 plot only/digitized |
| [Wang et al. (1996)](https://doi.org/10.1029/95jb03254) | Reproduced | 4 | 4 parity |  |
| [Wang et al. (2012)](https://doi.org/10.1029/2011jb009100) | Reproduced | 1 | 1 parity | 1 bundled |
| [Wang et al. (2026), KAlSi3O8 liebermannite and K-hollandite II](https://doi.org/10.2138/am-2024-9562) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Wentzcovitch et al. (1993), MgSiO3 molecular dynamics](https://doi.org/10.1103/PhysRevLett.70.3947) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Wittlinger et al. (1997)](https://doi.org/10.1107/s0108768197005739) | Reproduced | 1 | 1 similar | 1 plot only/digitized |
| [Wolf et al. (2015)](https://doi.org/10.1002/2015jb012108) | Reproduced | 2 | 2 parity | 2 bundled |
| [Wu et al. (2013)](https://doi.org/10.7498/aps.62.049101) | Direct refit unavailable | 2 | 2 direct refit unavailable | 2 theoretical parameterization only |
| [Wu et al. (2016), Fe-Al phase D elasticity](https://doi.org/10.1002/2016JB013209) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Xiao et al. (2013)](https://doi.org/10.2138/am.2013.4470) | Partly reproduced | 3 | 1 similar; 2 direct refit unavailable | 1 bundled |
| [Xu et al. (2020)](https://doi.org/10.1029/2020gl088877) | Reproduced | 1 | 1 parity | 1 bundled |
| [Xu et al. (2024), Al-bearing superhydrous phase B](https://doi.org/10.1029/2023GL107818) | Withheld: could not reproduce | 0 | no production record | investigation evidence only |
| [Yagi et al. (1992)](https://doi.org/10.1016/0031-9201(92)90063-2) | Reproduced | 1 | 1 similar | 1 bundled |
| [Yang et al. (2015)](https://doi.org/10.1038/srep17188) | Direct refit unavailable | 1 | 1 direct refit unavailable |  |
| [Yu et al. (2024)](https://doi.org/10.1029/2023jb028026) | Reproduced | 1 | 1 similar | 1 bundled |
| [Zha, Mao, and Hemley (2000), MgO elasticity pressure scale](https://doi.org/10.1073/pnas.240466697) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Zhang and Bukowinski (1991)](https://doi.org/10.1103/physrevb.44.2495) | Direct refit unavailable | 3 | 3 direct refit unavailable |  |
| [Zhang and Weidner (1999), Al-enriched silicate perovskite](https://doi.org/10.1126/science.284.5415.782) | Direct refit unavailable | 0 | no production record | investigation evidence only |
| [Zhang and Wentzcovitch (2022)](https://doi.org/10.1103/physrevb.106.054103) | Direct refit unavailable | 4 | 4 direct refit unavailable |  |
| [Zhang et al. (2009), Mg-Fe silicate phase stability](https://doi.org/10.1142/S0217979209053047) | Deferred: incomplete source/model mapping | 0 | no production record | investigation evidence only |
| [Zhang et al. (2025)](https://doi.org/10.3390/cryst15030221) | Reproduced | 3 | 3 parity | 3 bundled; 3 final-input parity, upstream reduction partial |
| [Zhao et al. (1997)](https://doi.org/10.1029/96gl03769) | Direct refit unavailable | 1 | 1 direct refit unavailable | 1 bundled |
| [Zhu et al. (2020)](https://doi.org/10.1029/2020jb019964) | Reproduced | 3 | 3 parity | 3 bundled |
| [Zhu et al. (2025)](https://doi.org/10.22541/essoar.176236186.65259830/v1) | Direct refit unavailable | 3 | 3 direct refit unavailable |  |

## Maintenance

This page is generated by `scripts/generate_paper_investigation_ledger.py`.
Update record-level evidence first, add nonproduction investigations to
`docs/data/nonproduction-paper-investigations.json`, regenerate this page,
and run the generator with `--check` in validation.
