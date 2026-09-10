# MgSiO3-family LitCurate source-exhaustion audit

Audit date: 2026-09-06

This audit covers every LitCurate row attached to the fourteen assigned publication DOIs. LitCurate is discovery-only: a row was promoted only when the primary publication established composition and phase, an executable equation family, a complete reference state, and all required coefficients. Values merely quoted from another paper are not duplicated. Values that cannot be checked in accessible primary full text remain nonproduction even when LitCurate looks numerically complete.

Five records pass that standard:

- `bridgmanite_dorogokupets_2015_298k_rydberg_stacey`
- `akimotoite_dorogokupets_2015_298k_rydberg_stacey`
- `mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey`
- `bridgmanite_hamahata_2000_md_300k_bm3`
- `mgsio3_liquid_akins_2004_adiabatic_bm3`

## Source-level conclusions

| DOI | Primary-source finding | Production disposition |
| --- | --- | --- |
| `10.1126/science.251.4992.410` | The Science article reports a phase transition and thermal expansion; the candidate table's assumed modulus and derivative are not a complete, source-owned EOS fit with a stated equation. | No record. Register as `not_an_eos_fit`. |
| `10.1029/96GL03027` | The source paper concerns a compositionally complex natural-basalt majoritic garnet. The two fitted alternatives require the reference volume and exact fit convention. The primary full article was not openly recoverable in this audit; a later secondary quotation is insufficient authority. The stishovite and CaSiO3 values are citations. | No record. Register as `primary_source_inaccessible`; revisit only with the final article and a dedicated complex-garnet material identity. |
| `10.1098/rsta.1996.0053` | The apparent MgSiO3-perovskite parameters are not independently separable into a complete pure-endmember source fit from the accessible source evidence, and one row is explicitly ferropericlase cited from elsewhere. | No record. Register as `insufficient_source_owned_parameterization`. |
| `10.1088/0256-307X/17/3/022` | Sound-velocity temperature-coefficient study; the lone 263 GPa number is citation-reported and has neither reference volume nor pressure derivative. | No record. Register as `citation_only`. |
| `10.2138/am-2000-2-309` | LitCurate prints a numerically complete static BM3 for MgSiO3 ilmenite/akimotoite, but the final primary article was not recoverable from the publisher or repository endpoints during this audit. The remaining three rows are earlier experimental citations. | Hold row 288; no record until the equation, volume basis, computational method, and fit interval are checked in primary full text. Register as `primary_source_inaccessible`. |
| `10.2138/am-2000-2-310` | LitCurate prints source-owned Vinet inputs for bridgmanite, CaSiO3 perovskite, and magnesiowustite, but the final 2000 article was not recoverable. A related 1998 proceedings paper confirms the general Vinet-plus-Debye method but cannot establish the exact same-DOI coefficients or compositions. | Hold row 292; rows 293-294 are outside this task's file ownership. Register this DOI as `primary_source_inaccessible`. |
| `10.2465/jmps.95.236` | The open final PDF explicitly gives a third-order Birch-Murnaghan density equation and, for the 300 K MgSiO3-perovskite MD compression, rho0=4.030 g cm-3, KT0=271 GPa, and KT0'=4.4. | Accept row 314 after transparent density-to-Z=4 volume conversion. Row 315 is a Funamori et al. citation. |
| `10.1029/2003GL018762` | This elasticity-temperature-derivative paper's five LitCurate entries are comparison/citation values. None has a source-owned reference volume or executable equation. | No record. Register as `citation_only`. |
| `10.1029/2004GL020237` | Paragraph 4 prints a source-owned candidate high-pressure MgSiO3-melt EOS: rho0=3.68 g cm-3, K0S=125 GPa, K0S'=4.0. The author derivation shows that theoretical Hugoniots are Mie-Gruneisen energy-balance offsets from a BM3 isentrope, not direct BM3 fits to shock P-rho pairs. A bounded reconstruction from the three later re-reduced enstatite melt states recovers rho0=3.6891 g cm-3 and K0S=125.463 GPa when K0S'=4 and the printed thermal terms are fixed. | Retain row 458 in production as a qualified published candidate reference isentrope. The reconstruction is `similar`, not strict parity or an independent EOS determination. |
| `10.1360/CJCP2006.19(4).311.4` | The open article reports the ambient MD molar volume and plots P-V behavior, but gives no complete fitted EOS equation and coefficients. The derivative-only row is a citation. | No record. Register as `incomplete_parameterization`. |
| `10.1088/1674-0068/20/05/547-551` | The article reports an equilibrium volume for a new interaction potential, not a complete EOS parameterization. | No record. Register as `incomplete_parameterization`. |
| `10.1111/j.1365-246X.1975.tb06463.x` | The 1975 shock-interpretation article is incorrectly year-labelled 2007 in LitCurate. Its tabulated compressibility/modulus values lack the reference volumes and, except one mixed-composition estimate, pressure derivatives needed for executable Murnaghan records. The pure-Mg endpoint is extrapolated/cited rather than an independent complete fit. | No record. Register as `incomplete_legacy_parameterization`. |
| `10.1016/j.pnsc.2009.09.002` | The source row is the authors' re-fit of Li and Zhang elastic data, not a new P-V dataset; exact Fe content and reference volume are absent. The pure-MgSiO3 and ferropericlase rows are explicit citations. | No record. Register as `secondary_refit_without_identity_or_v0`. |
| `10.1016/j.rgg.2015.01.011` | Equation (2) and Table 1 define a complete generalized Rydberg reference isotherm. It maps exactly to Peritheos Rydberg-Stacey with `K_infinity_prime=k/3`; the source fixes `k=5`. Its full Helmholtz thermal model uses two Einstein modes not represented by these records. | Accept rows 867-869 as 298.15 K reference isotherms. Row 870 is an MgO citation. |

## Row-level disposition

| Source row | LitCurate identifier | Material / reported content | Disposition |
| ---: | --- | --- | --- |
| 44 | `litcurate_1bb0547da02a92e1` | MgSiO3 bridgmanite; V0=24.725 cm3/mol, K0=246.5 GPa, K'=4, all marked assumed and equation unknown | Reject: no complete source-owned EOS. |
| 45 | `litcurate_b7ac10e23ee2799a` | MgSiO3 bridgmanite; no EOS coefficients | Reject: transition/thermal-expansion result, not an EOS fit. |
| 150 | `litcurate_941ab5276debb9a2` | Complex basaltic majorite; BM3 K0=226.2 GPa, K'=4 | Hold: primary V0/context not independently recovered. |
| 151 | `litcurate_4ee8a39e356e6c1b` | Same majorite; alternate BM3 K0=180 GPa, K'=7 | Hold: primary V0/context not independently recovered. |
| 152 | `litcurate_e3d86575aafbd3c5` | Stishovite K0=342 GPa, K'=4 | Reject: citation and incomplete. |
| 153 | `litcurate_820577b09e8377c9` | CaSiO3 perovskite K0=281 GPa, K'=4 | Reject: citation, incomplete, and outside assigned materials. |
| 156 | `litcurate_679d39dcc8734634` | (Mg,Fe)SiO3 bridgmanite BM3, V0=24.46 cm3/mol, K0=261 GPa, K'=3.75 | Reject: composition is not specified and independence from the cited comparison fit is not defensible. |
| 157 | `litcurate_de8b4d32fd97998f` | MgSiO3 bridgmanite alternate K0=261 GPa, K'=4, no V0 | Reject: incomplete alternate/citation fit. |
| 158 | `litcurate_4e1dd9a14f024f75` | (Mg0.6Fe0.4)O BM3 | Reject: citation and outside assigned materials. |
| 281 | `litcurate_97bae74a2c52109d` | (Mg,Fe)SiO3 bridgmanite K0=263 GPa only | Reject: citation and incomplete. |
| 288 | `litcurate_d4b33f075a456c53` | MgSiO3 akimotoite static BM3, V0=252.75 A3 (Z=6), K0=224 GPa, K'=4.18 | Hold: plausible and numerically complete, but primary full text was unavailable. |
| 289 | `litcurate_b5db8183c3658c78` | MgSiO3 akimotoite experimental V0=262.54 A3, K0=212 GPa | Reject: citation and incomplete. |
| 290 | `litcurate_3c972fc0cc8aa1d6` | MgSiO3 akimotoite K0=212 GPa, K'=5.6 | Reject: cited Reynard ice-scale fit, missing V0. |
| 291 | `litcurate_a15a1c169a510297` | MgSiO3 akimotoite K0=212 GPa, K'=7.5 | Reject: cited Reynard ruby-scale fit, missing V0. |
| 292 | `litcurate_ab6531269ad0df87` | nominal (Mg,Fe)SiO3 perovskite Vinet V0=40.185 A3/f.u., K0=266 GPa, K'=4.2 | Hold: final primary full text and exact composition/context unavailable. |
| 293 | `litcurate_8a49fd0cc8b3621c` | CaSiO3 perovskite Vinet | Not handled: source hold plus outside assigned materials. |
| 294 | `litcurate_a21656f2d1f9011c` | magnesiowustite Vinet | Not handled: source hold plus outside assigned files. |
| 314 | `litcurate_d6b5cf68e4b9f829` | MgSiO3 bridgmanite 300 K MD BM3, K0=271 GPa, K'=4.4 | Accept; primary rho0 supplies V0. |
| 315 | `litcurate_0a8330d27f2e169d` | MgSiO3 bridgmanite BM3 K0=261 GPa, K'=4 | Reject: cited Funamori et al. value. |
| 442 | `litcurate_69a3d8609023f4de` | MgSiO3 bridgmanite K0=249 GPa, K'=4 | Reject: citation and no V0/equation. |
| 443 | `litcurate_4d235e278b4f8253` | MgSiO3 bridgmanite K0=253 GPa, K'=3.9 | Reject: citation and no V0/equation. |
| 444 | `litcurate_07be5b78c3a71fe7` | (Mg,Fe)SiO3 bridgmanite K0=261 GPa, K'=4 | Reject: citation, composition unspecified, no V0. |
| 445 | `litcurate_ace0da5a84beb5b0` | MgSiO3 bridgmanite K0=261 GPa, K'=4 | Reject: citation and no V0/equation. |
| 446 | `litcurate_fcb3515d6ac5f0ab` | MgSiO3 bridgmanite K0=261 GPa, K'=4 | Reject: citation and no V0/equation. |
| 458 | `litcurate_042ccd0ac8e39a2b` | MgSiO3 high-pressure melt, K0S=125 GPa, K0S'=4 | Accept after recovering V0=45.29797155879917 A3/f.u. from rho0=3.68 g cm-3 and confirming the source's BM3-isentrope construction. |
| 581 | `litcurate_1ab29fe7e3d0d463` | MgSiO3 perovskite V0=24.414 cm3/mol only | Reject: no equation, K0, or K'. |
| 582 | `litcurate_93fb6ba2631fe798` | MgSiO3 perovskite K'=4.92 only | Reject: citation and incomplete. |
| 634 | `litcurate_d3b79918949ac92d` | MgSiO3 perovskite V0=25.071 cm3/mol only | Reject: equilibrium volume, not a complete EOS. |
| 635 | `litcurate_c2faa35e1545f401` | MgSiO3 akimotoite legacy modulus/compressibility | Reject: no V0 or K'. |
| 636 | `litcurate_5e358f21fa9797b4` | FeSiO3 akimotoite legacy modulus/compressibility | Reject: no V0 or K'. |
| 637 | `litcurate_0af2ce351344d13e` | En90 bronzitite, K'=2 plus legacy modulus/compressibility | Reject: no V0 and phase/composition endpoint is not a distinct complete EOS. |
| 638 | `litcurate_e4944226cefe3b1c` | (Mg0.32Fe0.68)SiO3 legacy modulus/compressibility | Reject: no V0 or K'. |
| 639 | `litcurate_cf5bd7af5ebcabb3` | (Mg0.44Fe0.56)SiO3 legacy modulus/compressibility | Reject: no V0 or K'. |
| 640 | `litcurate_5f8f168f71634c5d` | (Mg0.72Fe0.28)SiO3 legacy modulus/compressibility | Reject: no V0 or K'. |
| 641 | `litcurate_4e5a26bc205113c9` | second (Mg0.72Fe0.28)SiO3 legacy value | Reject: no V0 or K'. |
| 642 | `litcurate_d6dcebea49ff0757` | (Mg0.15Fe0.85)SiO3 legacy modulus/compressibility | Reject: no V0 or K'. |
| 643 | `litcurate_5421eb06202b1f14` | extrapolated/cited MgSiO3 endpoint | Reject: citation and incomplete. |
| 676 | `litcurate_8e787df413112a0d` | authors' BM3 re-fit of (Mg,Fe)SiO3 elastic data, K0S=253 GPa, K'=4.3 | Reject: no new P-V data, exact x and V0 absent. |
| 677 | `litcurate_d5af969d188c7efe` | MgSiO3 K0S=253 GPa, K'=4.4 | Reject: cited Li and Zhang value, no V0. |
| 678 | `litcurate_af4797b895f5d868` | ferropericlase BM3 | Reject: citation, incomplete, outside assigned files. |
| 867 | `litcurate_7828334f1795f436` | MgSiO3 perovskite V0=24.45 cm3/mol, K0=252.0 GPa, K'=4.38 | Accept as corrected generalized Rydberg-Stacey, not LitCurate's BM label. |
| 868 | `litcurate_ef9e7ca080ec8309` | MgSiO3 post-perovskite V0=24.2 cm3/mol, K0=253.7 GPa, K'=4.03 | Accept as corrected generalized Rydberg-Stacey. |
| 869 | `litcurate_8acd9ecf5854a501` | MgSiO3 akimotoite V0=26.35 cm3/mol, K0=215.3 GPa, K'=4.91 | Accept as corrected generalized Rydberg-Stacey. |
| 870 | `litcurate_64cfd8b62909fcb2` | MgO periclase comparison EOS | Reject: citation and outside assigned files. |

## Reproduction details

Dorogokupets et al. give

`P0(V) = 3 K0 X^(-k) (1-X) exp[(3 K0'/2-k+1/2)(1-X)]`, with `X=(V/V0)^(1/3)`.

Peritheos `RydbergStacey` uses the same expression with `k=3 K_infinity_prime`, so the source's fixed `k=5` is exactly `K_infinity_prime=5/3`. Molar volumes are converted to conventional-cell volumes using exact `N_A` and Z=4, 6, and 4 respectively.

For Hamahata et al., `V0 = Z M/(rho0 N_A)` with Z=4, M=100.387 g/mol, and rho0=4.030 g/cm3 gives 165.45561819988183 A3. For Akins et al., the same calculation with Z=1 and rho0=3.68 g/cm3 gives 45.29797155879917 A3 per formula unit.

### Akins et al. candidate-melt reconstruction

The 2004 publisher page still lists separate supporting files for the new shock
states, model parameters, and derivation. Its automated download endpoint
returned a challenge during this audit, so no byte identity with the final EDS
text files is claimed. The audit instead triangulates the source content through
the official Caltech copy of Akins's thesis, which contains the underlying 14-row
state table and equations 2.6-2.18, and Mosenfelder et al. (2009), which
re-tabulates the Akins shots and explicitly documents the reduction.

That distinction matters. Initial density, flyer velocity, and shock velocity
are experimental observables. Particle velocity, peak pressure, shock density,
and internal energy are calculated by impedance matching. Phase labels are
interpretations. Finally, the candidate pressure at a given shock density is a
calculated Hugoniot obtained from

`Delta_EH = E_transition + Delta_ES + Delta_EV`

and `Delta_EV = V_H (P_H-P_S)/gamma`, with the Rankine-Hugoniot energy relation
and `gamma=gamma0(V/V0)^q`. Solving the source equations gives

`P_H = [P_S - gamma(E_transition+E_S)/V_H] / [1 - gamma(V_initial-V_H)/(2V_H)]`,

where `P_S` and `E_S` are the BM3 reference-isentrope pressure and energy. Thus
regressing the three liquid P-rho states directly with a BM3 would be circular
and physically wrong.

`scripts/reproduce_akins_2004_mgsio3_liquid.py` implements the complete pressure
construction. Using the three Akins enstatite states classified as melt in the
official 2009 re-reduction (shots 318, 322, and 319), an unweighted bounded
pressure-residual fit refines only `rho0` and `K0S`. It fixes `K0S'=4`,
`gamma0=2.4`, `q=1`, and `E_transition=2.4 MJ/kg` to the final 2004 values.
The printed `Cv=0.92(3nR)` term is retained as provenance but does not enter
this pressure-only energy-balance reconstruction; it is needed to calculate
temperature along the candidate Hugoniot. The pressure-fit result is:

| Quantity | Published candidate | Bounded reconstruction |
|---|---:|---:|
| `rho0` (g cm-3) | 3.68 | 3.68912 |
| `V0` (A3/formula unit) | 45.29797 | 45.18598 |
| `K0S` (GPa) | 125 | 125.46303 |
| pressure RMSE (GPa) | 6.84602 | 6.79454 |

This is strong curve-level consistency but weak parameter identification. Three
states and two free coefficients leave one residual degree of freedom; the
`rho0`-`K0S` correlation is 0.99964 and the diagnostic one-sigma errors are much
larger than the tiny difference between the point estimates. A pressure-error-
weighted sensitivity fit moves to about 3.39 g cm-3 and 97 GPa, confirming that
the unpublished source objective matters. `K0S'` and the thermal terms cannot be
independently fitted from this subset. The published values therefore remain
unchanged, carry no invented uncertainties, and stay in production only as an
explicitly labeled candidate reference isentrope—not as an ambient liquid
isotherm, a direct shock-data BM3, or a complete executable Hugoniot.

`scripts/reproduce_mgsio3_source_exhaustion.py` evaluates the serialized models against independent direct equations at several compression ratios. `scripts/reproduce_akins_2004_mgsio3_liquid.py` performs the separate source-equation Hugoniot validation described above.

## Zotero-ready metadata

No Zotero write was performed in this task. The five production records require three bibliography items:

- Dorogokupets, P. I.; Dymshits, A. M.; Sokolova, T. S.; Danilov, B. S.; Litasov, K. D. (2015). *The equations of state of forsterite, wadsleyite, ringwoodite, akimotoite, MgSiO3-perovskite, and postperovskite and phase diagram for the Mg2SiO4 system at pressures of up to 130 GPa*. Russian Geology and Geophysics 56, 172-189. DOI `10.1016/j.rgg.2015.01.011`.
- Hamahata, K.; Ohtani, E.; Kawamura, K. (2000). *Molecular Dynamics (MD) simulation of the elasticity of MgSiO3 perovskite and the temperature anomaly in the lower mantle*. Journal of Mineralogical and Petrological Sciences 95, 236-244. DOI `10.2465/jmps.95.236`.
- Akins, J. A.; Luo, S.-N.; Asimow, P. D.; Ahrens, T. J. (2004). *Shock-induced melting of MgSiO3 perovskite and implications for melts in Earth's lowermost mantle*. Geophysical Research Letters 31, L14612. DOI `10.1029/2004GL020237`.
