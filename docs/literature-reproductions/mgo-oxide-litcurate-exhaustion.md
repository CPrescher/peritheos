# MgO and related-oxide LitCurate source-exhaustion audit

Audit date: 2026-09-06.

This audit covers all 43 LitCurate candidates whose extracting-paper DOI is
one of the twelve DOI families below. It applies the source-ownership and
minimum-parameter rules in `docs/adding-materials-and-eos.md`: a comparison or
adopted-input value remains attributed to its generating primary source, and
an EOS is not made executable by borrowing a missing reference volume,
pressure derivative, equation identity, or reference state from another
paper.

Result: **0 new production records; 11 incomplete source rows held; 32
citation, adopted-input, malformed, or non-EOS rows rejected.** No material
card, numerical dataset, or global ledger was changed by this audit.

## Source-level conclusions and candidate inventory

### 10.1063/1.555873 — Anderson and Zou (1990)

Primary article: “Thermodynamic Functions and Properties of MgO at High
Compression and High Temperature,” *Journal of Physical and Chemical
Reference Data* **19**, 69–83 (1990), DOI
[10.1063/1.555873](https://doi.org/10.1063/1.555873). The NIST record and
indexed article description were checked. This is a temperature-dependent
thermodynamic representation over 300–2000 K and up to 150 GPa; Birch–Murnaghan
terms are part of that construction rather than a standalone ambient BM
triplet. The extraction supplies neither `V0` nor `K0'` for its source row.

| Candidate | Disposition |
|---|---|
| `litcurate_a8448a761108a6f3` (row 25) | hold — source thermodynamic/Birch–Murnaghan construction, but equation order, `V0`, `K0'`, and the complete thermal parameterization are absent. |
| `litcurate_55c765b09fe97d93` (row 26) | reject — citation-reported modulus (`163.9 GPa`) with no complete source-owned EOS. |

Outcome: 0 accepted; 1 source row held; 1 citation row rejected.

### 10.1088/0953-8984/2/39/003 — Kumari and Dass (1990)

Primary article: “An equation of state applied to 50 solids. II,” *Journal of
Physics: Condensed Matter* **2**, 7891–7895 (1990), DOI
[10.1088/0953-8984/2/39/003](https://doi.org/10.1088/0953-8984/2/39/003).
The article metadata and its connection to the separately published equation
derivation in part I were checked. The MgO fit is reported as
`K0 = 1487.08 kbar` (`148.708 GPa`) and `K0' = 5.519`, but the LitCurate row
does not establish the equation or an absolute `V0`. The comparison row is
another paper's result.

| Candidate | Disposition |
|---|---|
| `litcurate_9044c5d727e832bd` (row 27) | hold — source fit, but equation identity and absolute `V0` are unresolved. |
| `litcurate_7bf02f430efa4b1b` (row 28) | reject — citation comparison (`1691.8 kbar`, `K0'=3.95`), not a source-owned fit. |

Outcome: 0 accepted; 1 source row held; 1 citation row rejected.

### 10.1126/science.257.5073.1099 — Stixrude et al. (1992)

Primary article: “Thermoelasticity of Silicate Perovskite and
Magnesiowüstite and Stratification of the Earth's Mantle,” *Science* **257**,
1099–1101 (1992), DOI
[10.1126/science.257.5073.1099](https://doi.org/10.1126/science.257.5073.1099).
The version-of-record metadata and indexed abstract were checked. These rows
belong to mantle-composition calculations for Fe-bearing bridgmanite and
magnesiowüstite, not pure MgO. The extraction does not identify the exact
Birch–Murnaghan order or reference temperature/pressure; one row additionally
lacks both modulus coefficients. Assigning BM3 merely because three numbers
are present would be an unsupported equation mapping.

| Candidate | Disposition |
|---|---|
| `litcurate_c7f61c55b52e2798` (row 46), `(Mg0.9Fe0.1)SiO3` | hold — `V0`, `K0`, and `K0'` are present, but BM order and reference state are not established and the composition needs its own defensible phase record. |
| `litcurate_20802bef7aad971e` (row 47), `(Mg,Fe)SiO3` | hold — assemblage/model input with only `V0=24.46 cm3/mol`; coefficients and composition are incomplete. |
| `litcurate_dd9afac1717aca07` (row 48), `(Mg0.6Fe0.4)O` | hold — BM order/reference state are unresolved for this composition-specific model input. |

Outcome: 0 accepted; 3 incomplete/identity-ambiguous source rows held.

### 10.1073/pnas.240466697 — Zha, Mao, and Hemley (2000)

Open primary article: “Elasticity of MgO and a primary pressure scale to
55 GPa,” *Proceedings of the National Academy of Sciences* **97**,
13494–13499 (2000), DOI
[10.1073/pnas.240466697](https://doi.org/10.1073/pnas.240466697). The primary
HTML, equations, Table 1, Table 2, and method text were checked. Table 2 gives
the source elastic result `K0S=162.5(7) GPa`, `K0S'=3.99(3)` and the corrected
isothermal values `K0T=160.2(7) GPa`, `K0T'=4.03(3)`. The text says the
adiabatic bulk modulus was fitted with a **third-order** finite-strain
equation; the LitCurate “fourth-order” label for row 275 is therefore wrong.
The paper does not provide an absolute reference volume in that parameter
table, so neither modulus pair alone defines an executable pressure-volume
EOS. Rows 278–279 are explicitly elastic measurements rather than EOS fits.

| Candidate | Disposition |
|---|---|
| `litcurate_7b046223b3101448` (row 275) | hold — source adiabatic elastic fit, `V0` absent and equation order mis-extracted. |
| `litcurate_c1738b470211be79` (row 276) | hold — source isothermal conversion, but `V0` absent. |
| `litcurate_3f9996ab942c1ba5` (row 277) | reject — cited comparison and incomplete fourth-order label. |
| `litcurate_727e2d6c5fd818c5` (row 278) | reject non-EOS — cited elastic result. |
| `litcurate_897da44bf5e7e06a` (row 279) | reject non-EOS — cited elastic result. |
| `litcurate_b09be3bd1df58998` (row 280) | reject — cited theoretical comparison, incomplete for the claimed fourth-order equation. |

Outcome: 0 accepted; 2 source rows held; 4 citation/non-EOS rows rejected.

### 10.1063/1.1371498 — Tsuchiya and Kawamura (2001)

Primary target: “Systematics of elasticity: Ab initio study in B1-type
alkaline earth oxides,” *The Journal of Chemical Physics* **114**,
10086–10093 (2001), DOI
[10.1063/1.1371498](https://doi.org/10.1063/1.1371498). The version-of-record
metadata and indexed method/abstract were checked. The source MgO calculation
is labelled BM3 and reports `K0=160 GPa`, `K0'=4.01`, but the extraction and
accessible primary material do not supply its `V0` or a numerical P–V table.
The other five rows are literature comparisons using several methods and no
resolved common equation.

| Candidate | Disposition |
|---|---|
| `litcurate_55c3bb2842521015` (row 357) | hold — source BM3 calculation missing `V0` and primary observations. |
| `litcurate_4b1de32da216cbb2` (row 358) | reject — citation comparison. |
| `litcurate_8abfa93e0da535bd` (row 359) | reject — citation comparison, also missing `K0'`. |
| `litcurate_ec2e580699ce5ff3` (row 360) | reject — citation comparison. |
| `litcurate_743717339d758ee9` (row 361) | reject — citation comparison. |
| `litcurate_ba7dd4a0b054456c` (row 362) | reject — experimental citation comparison. |

Outcome: 0 accepted; 1 source row held; 5 citation rows rejected.

### 10.1515/zna-2008-1-209 — Liu (2008)

Primary article: “Bulk Modulus and Equation of State Under the Effect of High
Temperature and High Pressure for MgO,” *Zeitschrift für Naturforschung A*
**63**, 53–56 (2008), DOI
[10.1515/zna-2008-1-209](https://doi.org/10.1515/zna-2008-1-209). The indexed
primary text, equations, and table were checked. This is a Harrison/Hildebrand
method paper; the `180.1 GPa` ambient modulus is adopted from the cited
pressure-induced-brillouin result, not newly fitted here. The extraction has
no `V0`, derivative, or complete source-owned parameterization.

| Candidate | Disposition |
|---|---|
| `litcurate_34bf9c16c3c77453` (row 653) | hold — method output is incomplete as an executable EOS and uses an adopted ambient modulus. |
| `litcurate_4f1bfcf7d0cbdbb4` (row 654) | reject — citation-reported modulus only. |

Outcome: 0 accepted; 1 incomplete/adopted-input source row held; 1 citation
row rejected.

### 10.1139/p11-040 — Liu (2011)

Primary article: “Analysis of a new two-parameter equation of state for MgO,”
*Canadian Journal of Physics* **89**, 709–712 (2011), DOI
[10.1139/p11-040](https://doi.org/10.1139/p11-040). The primary method paper
uses temperature-specific modulus pairs from earlier literature to exercise a
normalized-volume two-parameter equation. All five LitCurate entries are
correctly marked citation-reported, omit `V0`, and are not independent fits by
this paper.

| Candidate | Disposition |
|---|---|
| `litcurate_00497aab579b8f6d` (row 751), 300 K | reject — citation/adopted input (`K0=180 GPa`, `K0'=4.15`), `V0` absent. |
| `litcurate_86314accee4bc183` (row 752), 500 K | reject — citation/adopted input, `V0` absent. |
| `litcurate_b483f92ebb7b5815` (row 753), 1000 K | reject — citation/adopted input, `V0` absent. |
| `litcurate_23b7037d8b70e62f` (row 754), 1500 K | reject — citation/adopted input, `V0` absent. |
| `litcurate_08212e46ee5707a1` (row 755), 2000 K | reject — citation/adopted input, `V0` absent. |

Outcome: 0 accepted; 5 citation/adopted-input rows rejected.

### 10.1155/2014/289353 — Kholiya, Chandra, and Verma (2014)

Open primary article: “Analysis of Equation of States for the Suitability at
High Pressure: MgO as an Example,” *The Scientific World Journal* **2014**,
1–5, DOI [10.1155/2014/289353](https://doi.org/10.1155/2014/289353). The
primary PDF was checked directly (SHA-256
`ebe7fb5fbf2b34a27944ae8432b5ccbdb01a2ce1ee01c44f7a33e5600e6448b8`),
including equations (2)–(6), Table 1, and the normalized compression plots.
The paper compares its own algebraic form with Shanker, Tait, Vinet, and BM3,
but obtains all eight modulus pairs from cited publications. It adopts Li et
al.'s `K0=161.3 GPa`, `K0'=4.24` for its comparisons. Every plotted abscissa is
`V/V0`; no absolute MgO reference volume is reported. Digitizing those curves
would only reproduce predictions generated from cited inputs and would still
not determine the missing absolute `V0`.

| Candidate | Disposition |
|---|---|
| `litcurate_2bd1963c975f50d4` (row 849) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_d892dd55bca8897a` (row 850) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_4489bbcda3858896` (row 851) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_1ae823e505dcb265` (row 852) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_89fa059945c900b0` (row 853) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_ceeb093971b691fa` (row 854) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_29ebf9a63b275875` (row 855) | reject — citation row, equation unresolved, `V0` absent. |
| `litcurate_c1e3dd00ab220a85` (row 856) | reject — Li et al. citation adopted for a BM3 comparison; `V0` absent. |

Outcome: 0 accepted; 8 citation/adopted-input rows rejected.

### 10.7566/JPSJ.88.065003 — Fuchizaki (2019)

Primary target: “Predicting the Melting Curve of MgO: An Essential Update,”
*Journal of the Physical Society of Japan* **88**, 065003 (2019), DOI
[10.7566/JPSJ.88.065003](https://doi.org/10.7566/JPSJ.88.065003). The
publisher metadata, abstract, related project report, and reference inventory
were checked. The two-page version of record is access-controlled and no
authoritative open manuscript was located. The source-labelled extraction
(`K0=141 GPa`, fixed `K0'=20`) is anomalous relative to both comparison rows
and lacks `V0`; without primary equation/table access it cannot be determined
whether `20` is an EOS derivative, a range/point count, or a column-shift
artifact. It is not safe to map this row to Peritheos's generalized
Rydberg–Stacey model.

| Candidate | Disposition |
|---|---|
| `litcurate_416aac9efea84337` (row 1040) | hold — source label, but primary equation/parameter mapping is inaccessible, `V0` is absent, and the extracted fixed derivative is suspect. |
| `litcurate_bd84babf3daf78ee` (row 1041) | reject — Kono et al. citation comparison, equation and `V0` unresolved. |
| `litcurate_5aee4c20a9e5a489` (row 1042) | reject — Okamoto/Fuchizaki 2017 citation comparison, equation and `V0` unresolved. |

Outcome: 0 accepted; 1 source row held; 2 citation rows rejected.

### 10.12693/aphyspola.140.131 — Singh and Singh (2021)

Open primary article: “New Formulation of Equation of State and Study of
Elastic Properties of Alkaline Earth Oxides under High Pressure,” *Acta
Physica Polonica A* **140**, 131–137 (2021), DOI
[10.12693/aphyspola.140.131](https://doi.org/10.12693/aphyspola.140.131).
The primary PDF was checked directly (SHA-256
`50fa9e5ba56b061ae5d13637e94b587a442339113d57beeaa7ddc39c0a4b5b33`),
especially equations (30), (36), (39), (42), Table I, and Figure 2. The paper
derives an `n`-powered Eulerian-strain family; `n=2` recovers BM3. Table I
explicitly says its MgO (`161.6 GPa`, `4.13`) and CaO (`110.6 GPa`, `4.05`)
inputs are recorded in the cited references. Results are plotted only against
`V/V0`, so the paper supplies no absolute reference volume and no new fit.

| Candidate | Disposition |
|---|---|
| `litcurate_056ba335b5e3a8b0` (row 1116), MgO | reject — citation/adopted input, `V0` absent. |
| `litcurate_ed293122ae8307df` (row 1117), CaO | reject — citation/adopted input, `V0` absent. |

Outcome: 0 accepted; 2 citation/adopted-input rows rejected.

### 10.32908/hthp.v53.1503 — Vijay (2024)

Primary target: “Thermoelastic properties of materials based on the
Generalized Rydberg-Vinet and the Stacey reciprocal K-primed equations of
state,” *High Temperatures–High Pressures* **53**, 219–232 (2024), DOI
[10.32908/hthp.v53.1503](https://doi.org/10.32908/hthp.v53.1503). The
publisher issue record and indexed article metadata were checked. Both MgO
rows are marked citation-reported and repeat the same adopted input pair
(`K0=162 GPa`, `K0'=4.15`) under two comparison equations. Neither has `V0`
or a reference temperature; they are not two independent primary fits.

| Candidate | Disposition |
|---|---|
| `litcurate_90732f574a1451fc` (row 1263) | reject — citation input used in a generalized Rydberg–Vinet calculation, `V0`/temperature absent. |
| `litcurate_5678d4f58b736bf0` (row 1264) | reject — same citation input reused in a Stacey calculation, not an independent EOS. |

Outcome: 0 accepted; 2 citation/adopted-input rows rejected.

### 10.1063/1.1483512 — Braithwaite (2002)

Primary article: “Thermodynamic Representations for Solid/Melt Systems at
High Pressure and Temperature,” *AIP Conference Proceedings* **620**,
185–190 (2002), DOI
[10.1063/1.1483512](https://doi.org/10.1063/1.1483512). The indexed complete
primary text, equations (4)–(5), Table 2 header, and surrounding method text
were checked. Equation (4) is a Vinet reference-isotherm relation; equation
(5) adds thermal pressure relative to `TR`. Table 2's published columns are
phase, `TR`, `B0`, `B0'`, and `alpha`; there is **no `V0` column**. Consequently
the source liquid-MgO pair at `TR=3105 K` cannot define the Vinet curve in
absolute volume. The available text extraction also interleaves the solid and
liquid MgO table lines, explaining the impossible `3147 GPa` solid-row value.
The plotted curves do not expose a tabulated, phase- and basis-resolved
zero-pressure liquid volume suitable for transparent digitization.

| Candidate | Disposition |
|---|---|
| `litcurate_d9ae64f866b3a0bd` (row 410), solid MgO | reject malformed citation row — `3147` is a table/OCR column shift, not a 3147 GPa modulus; `V0` is also absent. |
| `litcurate_f9d91dc933b8d3fe` (row 411), liquid MgO | hold — source Vinet/thermal parameters include `K0=37.5 GPa`, `K0'=3.9`, `TR=3105 K`, but omit the required liquid `V0`; the complete thermal mapping is also not represented by a standalone Vinet record. |

Outcome: 0 accepted; 1 source row held; 1 malformed citation row rejected.

## Why digitization did not create a record

Digitization is allowed when a primary plot contains identifiable source
observations. It does not cure a missing reference state when the plot uses
only normalized `V/V0`, nor does it turn a model-prediction curve generated
from cited coefficients into new source observations. The Kholiya and Singh
papers show only normalized predictions; the PNAS paper's modulus table omits
the reference volume; the Braithwaite table and accessible plot text do not
provide a phase- and basis-resolved liquid `V0`; and the JPSJ primary figure is
not openly inspectable. No defensible digitization target remained.

## Proposed nonproduction-paper investigation entries

The following entries are proposed for
`docs/data/nonproduction-paper-investigations.json`. They are recorded here so
this scoped audit does not modify that shared file.

```json
[
  {"citation":"Anderson and Zou (1990), MgO thermodynamic functions","doi":"10.1063/1.555873","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The temperature-dependent thermodynamic construction is not a standalone complete BM record; the source row omits V0, K0 prime, equation order, and the full thermal parameterization."},
  {"citation":"Kumari and Dass (1990), EOS applied to 50 solids II","doi":"10.1088/0953-8984/2/39/003","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The MgO source row omits the absolute V0 and does not resolve the equation identity; its second row is citation-reported."},
  {"citation":"Stixrude et al. (1992), thermoelasticity and mantle stratification","doi":"10.1126/science.257.5073.1099","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The Fe-bearing assemblage rows do not establish exact Birch-Murnaghan order and reference state; one is additionally missing both modulus coefficients."},
  {"citation":"Zha, Mao, and Hemley (2000), MgO elasticity pressure scale","doi":"10.1073/pnas.240466697","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The primary source reports adiabatic and corrected isothermal modulus pairs but no absolute V0 in the EOS parameter table; the discovery row also mislabels the adiabatic fit order."},
  {"citation":"Tsuchiya and Kawamura (2001), B1 oxide elasticity","doi":"10.1063/1.1371498","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The source MgO BM3 row omits V0 and primary observations; the other five entries are citation comparisons."},
  {"citation":"Liu (2008), MgO bulk modulus method","doi":"10.1515/zna-2008-1-209","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The method paper adopts its ambient modulus from a cited result and does not provide V0, a derivative, or a complete source-owned executable EOS."},
  {"citation":"Liu (2011), two-parameter MgO EOS analysis","doi":"10.1139/p11-040","outcome":"withheld_unreproduced","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"All five temperature-dependent modulus pairs are citation-reported adopted inputs and omit V0; the paper generates no independent fit."},
  {"citation":"Kholiya et al. (2014), MgO EOS comparison","doi":"10.1155/2014/289353","outcome":"withheld_unreproduced","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"All eight rows quote earlier publications; the source compares normalized-volume predictions and reports no absolute V0 or new fit."},
  {"citation":"Fuchizaki (2019), updated MgO melting curve","doi":"10.7566/JPSJ.88.065003","outcome":"direct_refit_unavailable","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The source-labelled row has no V0 and an anomalous fixed derivative; the access-controlled primary equation/table could not be inspected to repair the suspected extraction error."},
  {"citation":"Singh and Singh (2021), alkaline-earth oxide EOS formulation","doi":"10.12693/aphyspola.140.131","outcome":"withheld_unreproduced","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The MgO and CaO modulus pairs are explicitly cited inputs; plots use only V/V0, and no new absolute-volume fit is reported."},
  {"citation":"Vijay (2024), generalized Rydberg-Vinet and Stacey thermoelasticity","doi":"10.32908/hthp.v53.1503","outcome":"withheld_unreproduced","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"Both MgO rows reuse the same citation-reported modulus pair under comparison equations and omit V0 and reference temperature."},
  {"citation":"Braithwaite (2002), high-pressure solid/melt representations","doi":"10.1063/1.1483512","outcome":"deferred_incomplete_model","investigation_date":"2026-09-06","evidence":"literature-reproductions/mgo-oxide-litcurate-exhaustion.md","reason":"The primary Vinet/thermal table omits V0. Its liquid-MgO row is therefore non-executable, while the solid row contains an obvious table-column extraction error."}
]
```

## Zotero-ready metadata

Metadata below was resolved from the DOI registration records. It is ready for
DOI-based Zotero import; this audit did not modify Zotero.

| DOI | Authors | Year | Journal / volume / pages |
|---|---|---:|---|
| `10.1063/1.555873` | Orson L. Anderson; Keshan Zou | 1990 | *Journal of Physical and Chemical Reference Data* 19(1), 69–83 |
| `10.1088/0953-8984/2/39/003` | M. Kumari; N. Dass | 1990 | *Journal of Physics: Condensed Matter* 2(39), 7891–7895 |
| `10.1126/science.257.5073.1099` | Lars Stixrude; R. J. Hemley; Y. Fei; H. K. Mao | 1992 | *Science* 257(5073), 1099–1101 |
| `10.1073/pnas.240466697` | Chang-Sheng Zha; Ho-kwang Mao; Russell J. Hemley | 2000 | *PNAS* 97(25), 13494–13499 |
| `10.1063/1.1371498` | T. Tsuchiya; K. Kawamura | 2001 | *The Journal of Chemical Physics* 114(22), 10086–10093 |
| `10.1515/zna-2008-1-209` | Quan Liu | 2008 | *Zeitschrift für Naturforschung A* 63(1–2), 53–56 |
| `10.1139/p11-040` | Quan Liu | 2011 | *Canadian Journal of Physics* 89(6), 709–712 |
| `10.1155/2014/289353` | Kuldeep Kholiya; Jeewan Chandra; Swati Verma | 2014 | *The Scientific World Journal* 2014, 1–5 |
| `10.7566/JPSJ.88.065003` | Kazuhiro Fuchizaki | 2019 | *Journal of the Physical Society of Japan* 88(6), 065003 |
| `10.12693/aphyspola.140.131` | S. P. Singh; D. Singh | 2021 | *Acta Physica Polonica A* 140(2), 131–137 |
| `10.32908/hthp.v53.1503` | A. Vijay | 2024 | *High Temperatures–High Pressures* 53(3), 219–232 |
| `10.1063/1.1483512` | M. Braithwaite | 2002 | *AIP Conference Proceedings* 620, 185–190 |

## Batch accounting

| DOI | Rows | Accepted | Held | Rejected |
|---|---:|---:|---:|---:|
| `10.1063/1.555873` | 2 | 0 | 1 | 1 |
| `10.1088/0953-8984/2/39/003` | 2 | 0 | 1 | 1 |
| `10.1126/science.257.5073.1099` | 3 | 0 | 3 | 0 |
| `10.1073/pnas.240466697` | 6 | 0 | 2 | 4 |
| `10.1063/1.1371498` | 6 | 0 | 1 | 5 |
| `10.1515/zna-2008-1-209` | 2 | 0 | 1 | 1 |
| `10.1139/p11-040` | 5 | 0 | 0 | 5 |
| `10.1155/2014/289353` | 8 | 0 | 0 | 8 |
| `10.7566/jpsj.88.065003` | 3 | 0 | 1 | 2 |
| `10.12693/aphyspola.140.131` | 2 | 0 | 0 | 2 |
| `10.32908/hthp.v53.1503` | 2 | 0 | 0 | 2 |
| `10.1063/1.1483512` | 2 | 0 | 1 | 1 |
| **Total** | **43** | **0** | **11** | **32** |
