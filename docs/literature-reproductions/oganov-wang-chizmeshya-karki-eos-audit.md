# MgO and CaSiO3 EOS audit: four candidate papers

## Outcome

Twenty source-owned BM3 parameterizations are accepted: eleven MgO fits from
Oganov and Dorogokupets (2003), four CaSiO3-perovskite fits or explicit
reanalyses from Wang et al. (1996), and five CaSiO3-perovskite fits from
Chizmeshya et al. (1996). None of the five Karki and Stixrude (1999) candidate
rows is a volume equation of state. This document accounts for every one of
the 29 assigned source rows; comparison values belonging to cited papers are
also identified and excluded.

The acceptance rule is source ownership plus a complete, explicitly identified
volume-EOS parameter set. Alternative fits are retained only where the source
actually performed and interpreted them. Repeated values embedded in a larger
thermal model, a unit-converted restatement, elastic-modulus fits, and citation
traces are not promoted to records.

## Oganov and Dorogokupets (2003): 11 accepted

Artem R. Oganov and Peter I. Dorogokupets, “All-electron and pseudopotential
study of MgO: Equation of state, anharmonicity, and stability,” *Physical
Review B* **67**, 224110 (2003),
doi:[10.1103/PhysRevB.67.224110](https://doi.org/10.1103/PhysRevB.67.224110).
The audit used the
[APS primary article](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.67.224110/fulltext).

Equations (3)-(4) and the text immediately below Table I identify the
third-order Birch-Murnaghan form. Table I reports four independent static GGA
fits that vary the electronic-structure treatment. Table V reports seven
independent fits to the pressure-corrected static/thermal EOS. The stated
-7.736 GPa pressure shift is part of the latter model family, so those rows are
not duplicates of Table I.

| Assigned source row | Source location | V0 (A3, Z=4) | K0 (GPa) | K0' | Disposition |
|---|---|---:|---:|---:|---|
| ECP, large-core Mg, static | Table I | 77.629 | 151.707 | 4.212 | ACCEPT: source-computed static BM3. |
| ECP, small-core Mg, static | Table I | 76.595 | 150.839 | 4.052 | ACCEPT: source-computed static BM3. |
| PAW, large-core Mg, static | Table I | 76.049 | 154.183 | 4.141 | ACCEPT: source-computed static BM3. |
| PAW, small-core Mg, static | Table I | 76.947 | 150.597 | 4.103 | ACCEPT: source-computed static BM3. |
| Pressure-corrected static | Table V | 73.425 | 181.240 | 3.997 | ACCEPT: explicitly refitted corrected static curve. |
| Pressure-corrected 0 K | Table V | 74.439 | 173.480 | 4.014 | ACCEPT: explicitly tabulated isotherm fit. |
| Pressure-corrected 298 K | Table V | 74.670 | 170.530 | 4.036 | ACCEPT: explicitly tabulated isotherm fit. |
| Pressure-corrected 1000 K | Table V | 76.549 | 152.595 | 4.130 | ACCEPT: explicitly tabulated isotherm fit. |
| Pressure-corrected 2000 K | Table V | 79.915 | 127.719 | 4.244 | ACCEPT: explicitly tabulated isotherm fit. |
| Pressure-corrected 3000 K | Table V | 83.772 | 106.110 | 4.331 | ACCEPT: explicitly tabulated isotherm fit. |
| Pressure-corrected 4000 K | Table V | 88.006 | 88.473 | 4.385 | ACCEPT: explicitly tabulated isotherm fit. |

Table I's additional `Experiment` row (74.712, 159.94, 4.112) is rejected as a
citation trace: its footnote says it was obtained from Refs. 1 and 22-26 rather
than from observations generated in this work.

Table III supplies eight small-core PAW pressure-volume checkpoints. They are
transcribed in `mgo-oganov-2003-table3-pv.csv` and used to test both the
uncorrected Table-I curve and the pressure-shifted Table-V static curve. They
are not claimed to be the original energy-volume fit grid.

## Wang, Weidner, and Guyot (1996): 4 accepted, 4 rejected

Yanbin Wang, Donald J. Weidner, and Francois Guyot, “Thermal equation of state
of CaSiO3 perovskite,” *Journal of Geophysical Research: Solid Earth* **101**,
661-672 (1996),
doi:[10.1029/95JB03254](https://doi.org/10.1029/95JB03254).
The full primary article was checked, including equation (1), Tables 1, 2 and
4, and the room-temperature EOS discussion on pages 664-665.

The preferred room-temperature BM3 is fitted to the twelve run-13 observations
at 2.66-10.07 GPa with K0'=4.8 adopted from the authors' reanalysis of Mao et
al. (1989). The two lower-pressure run-13 observations are published in the
dataset but excluded from the fit because the source reports amorphization and
large differential stress. The three Mao-data rows below are genuine
source-performed sensitivity reanalyses, not attributions to Mao's own fit.

| Assigned source row | V0 (A3) | K0 (GPa) | K0' | Disposition |
|---|---:|---:|---:|---|
| Preferred run-13 room-temperature BM3 | 45.58 | 232 | 4.8 | ACCEPT: source's preferred low-pressure fit; K0' fixed. |
| Mao data, equal weights, sub-1-GPa values excluded | 45.71 | 244 | 4.8 | ACCEPT: explicit source reanalysis and selected compatibility curve. |
| Mao data, pressure-uncertainty weighted, all values | 45.3 | 282 | 4.0 | ACCEPT: explicit source sensitivity fit; contaminated low-P warning retained. |
| Mao data, pressure-uncertainty weighted, sub-1-GPa values excluded | 45.47 | 268 | 4.3 | ACCEPT: explicit source sensitivity fit. |
| Thermal-pressure fit, weighted 1 | 45.60 | 233 | 4.8 | REJECT / INCOMPLETE THERMAL MODEL: these three values are fixed/static terms inside a multicoefficient thermal-pressure fit in Table 2, not another independently reported room-temperature EOS. |
| Thermal-pressure fit, weighted 2 | 45.61 | 229 | 4.8 | REJECT / INCOMPLETE THERMAL MODEL: alternative weighting of the composite thermal model; extracting only V0/K0/K0' would discard its thermal coefficients. |
| Thermal-pressure fit, unweighted | 45.58 | 233 | 4.8 | REJECT / INCOMPLETE THERMAL MODEL: alternative weighting of the same composite thermal model. |
| Lower-mantle parameter summary | 27.45 cm3/mol | 232 | 4.8 | REJECT / DUPLICATE: Table 4 restates the preferred fit, converting 45.58 A3 per formula unit to molar volume. |

Four further comparison values discussed by Wang et al. are citation traces
and remain under their primary sources: Mao et al. (1989)'s own constrained
fit (45.37 A3, 281 GPa, K0'=4), Tarrida and Richet (1989), Tamai and Yagi
(1989), and the prior Wang and Weidner (1994) dataset.

The reproduction script refits the twelve published run-13 points with K0'=4.8
and separately reruns the three Mao-data objectives using the existing complete
Mao Table-I transcription. The preferred run-13 coefficients are recovered
closely. Exact recovery of the pressure-weighted legacy fits is not claimed:
the article does not specify enough objective-function details to reproduce its
rounding and treatment of the reported P/V uncertainties uniquely.

## Chizmeshya, Wolf, and McMillan (1996): 5 accepted

Alex Chizmeshya, George H. Wolf, and Paul F. McMillan, “First-principles
calculation of the equation-of-state, stability, and polar optic modes of
CaSiO3 perovskite,” *Geophysical Research Letters* **23**, 2725-2728 (1996),
doi:[10.1029/96GL02624](https://doi.org/10.1029/96GL02624).

The source reports three successively converged static LAPW BM3 fits and two
300 K thermally corrected fits. The constrained thermal row is retained as an
explicit sensitivity fit rather than being silently merged with the preferred
unconstrained row.

| Assigned source row | Source location | V0 (A3) | K0 (GPa) | K0' | Disposition |
|---|---|---:|---:|---:|---|
| LAPW(7), static | Table 1 | 45.04 | 241.8 | 4.15 | ACCEPT: independent seven-point convergence fit. |
| LAPW(8), static | Table 1 | 45.02 | 241.0 | 4.16 | ACCEPT: independent eight-point convergence fit. |
| LAPW(9), static | Table 1 | 45.06 | 238.2 | 4.18 | ACCEPT: independent final static fit. |
| LAPW(9), 300 K, constrained | Table 2 | 45.55 | 237.9 | 4.00 | ACCEPT: explicitly reported K0'=4 sensitivity fit. |
| LAPW(9), 300 K, preferred | Abstract and Table 2 | 45.62 | 227 | 4.29 | ACCEPT: source's preferred thermally corrected fit. |

Five additional PHF/experimental comparison rows in the candidate extraction
are rejected as citation traces because the calculations or observations were
not generated by this paper. The preferred 45.62/227/4.29 set is also repeated
by Akber-Knutson et al. (2002), Table 2, doi:10.1029/2001GL013523.

The publisher links a formal one-page correction, “Correction to
First-principles calculation ...,” *GRL* **25**(5), 711 (1998),
doi:[10.1029/97GL02812](https://doi.org/10.1029/97GL02812). Its full correction
text was not independently retrievable during this audit. Accordingly every
record carries the correction DOI and the audit does not infer undocumented
changes; the later published comparison retaining the preferred set is the
available cross-check. No underlying E-V grid or fit weights are published, so
the deterministic check is limited to exact transcription and independent
curve evaluation.

## Karki and Stixrude (1999): 0 accepted, 5 rejected

Bijaya B. Karki and Lars Stixrude, “Seismic velocities of major silicate and
oxide phases of the lower mantle,” *Journal of Geophysical Research: Solid
Earth* **104**, 13025-13033 (1999),
doi:[10.1029/1999JB900069](https://doi.org/10.1029/1999JB900069).
The audit used an archived scan of the full primary article.

The five candidate rows arise from Table 1, whose heading is “Calculated
Athermal Elastic (M) Moduli and Their Pressure Derivatives (M') at Zero
Pressure.” Here `K` is the aggregate elastic bulk modulus and `K'` its elastic
pressure derivative. The third-order Eulerian finite-strain equations mentioned
in the table note are fits of each elastic modulus `M(P)`, not fits of volume
versus pressure. No V0 is supplied. The underlying elastic calculations are
also attributed to the authors' earlier papers.

| Assigned source row | Table-I K (GPa) | K' | Disposition |
|---|---:|---:|---|
| MgSiO3 perovskite, calculated | 260 | 4.02 | REJECT / NOT A VOLUME EOS: elastic-modulus row; V0 absent. |
| CaSiO3 perovskite, calculated | 236 | 4.42 | REJECT / NOT A VOLUME EOS: elastic-modulus row; V0 absent. |
| MgO, calculated | 157 | 4.27 | REJECT / NOT A VOLUME EOS: elastic-modulus row; V0 absent. |
| CaO, calculated | 115 | 4.48 | REJECT / NOT A VOLUME EOS: elastic-modulus row; V0 absent. |
| Stishovite, calculated | 310 | 4.24 | REJECT / NOT A VOLUME EOS: elastic-modulus row; V0 absent. |

The five experimental rows immediately below them are likewise elastic
comparisons sourced from earlier experimental papers. They are rejected both
as non-volume-EOS rows and as citation traces.

## Complete LitCurate candidate inventory

The tables above account for the 29 source-owned rows assigned to this batch.
For machine-auditable completeness, every extracted candidate under the four
DOIs is listed exactly once below, including all citation traces.

| Candidate | Disposition |
|---|---|
| `litcurate_495276da66c58057` | ACCEPT: Oganov Table-I ECP large-core. |
| `litcurate_1f5654f2f00f9f23` | ACCEPT: Oganov Table-I ECP small-core. |
| `litcurate_d49496a812e42045` | ACCEPT: Oganov Table-I PAW large-core. |
| `litcurate_f0dc0ccdf1cf08a5` | ACCEPT: Oganov Table-I PAW small-core. |
| `litcurate_9e84630d76a0440b` | REJECT / CITATION TRACE: aggregated experiment from Refs. 1 and 22-26. |
| `litcurate_8fe272d95d738aaf` | ACCEPT: Oganov Table-V corrected static. |
| `litcurate_bbc8127f4a8a8553` | ACCEPT: Oganov Table-V 0 K. |
| `litcurate_eeb995324a9e42ec` | ACCEPT: Oganov Table-V 298 K. |
| `litcurate_2a37886b36f784e4` | ACCEPT: Oganov Table-V 1000 K. |
| `litcurate_27542a1e72e4db93` | ACCEPT: Oganov Table-V 2000 K. |
| `litcurate_75a5273b7f434d75` | ACCEPT: Oganov Table-V 3000 K. |
| `litcurate_9b4cc923f365760b` | ACCEPT: Oganov Table-V 4000 K. |
| `litcurate_c8dd7ce993040ef3` | ACCEPT: Wang preferred run-13 fit. |
| `litcurate_429d4b0a6abc8ee1` | ACCEPT: Wang equal-weight Mao-data reanalysis. |
| `litcurate_228b6c1c84dc3934` | ACCEPT: Wang pressure-weighted all-data reanalysis. |
| `litcurate_ec8391e08cf13b41` | ACCEPT: Wang pressure-weighted reanalysis excluding sub-1-GPa points. |
| `litcurate_b96e6df515e03b31` | REJECT / CITATION TRACE: Mao et al. original fit. |
| `litcurate_f624d81245b377a1` | REJECT / CITATION TRACE: Tarrida and Richet result. |
| `litcurate_1b934f78ddb1e995` | REJECT / CITATION TRACE: Tamai and Yagi result. |
| `litcurate_a5214b57abb7e030` | REJECT / INCOMPLETE THERMAL MODEL: weighted-1 Table-II coefficients. |
| `litcurate_d3236851bdaabbb5` | REJECT / INCOMPLETE THERMAL MODEL: weighted-2 Table-II coefficients. |
| `litcurate_c860b4736d0a9c29` | REJECT / INCOMPLETE THERMAL MODEL: unweighted Table-II coefficients. |
| `litcurate_4bb90558a9c2e3fa` | REJECT / DUPLICATE: CaSiO3 Table-IV unit-converted preferred fit. |
| `litcurate_8eea64a6ac473deb` | REJECT / CITATION TRACE: MgSiO3 Table-IV data from Wang et al. (1994). |
| `litcurate_7f7acbac7cdbd1d8` | ACCEPT: Chizmeshya LAPW(7). |
| `litcurate_6a7d6917dfc7efdf` | ACCEPT: Chizmeshya LAPW(8). |
| `litcurate_c69b68c6ef4ad4e2` | ACCEPT: Chizmeshya LAPW(9). |
| `litcurate_c875851dc0cf2eba` | REJECT / CITATION TRACE: Sherman (1993) PHF static. |
| `litcurate_ebd3f572b1dd4d16` | ACCEPT: Chizmeshya constrained 300 K LAPW(9). |
| `litcurate_bbc8817045d7f711` | ACCEPT: Chizmeshya preferred 300 K LAPW(9). |
| `litcurate_c1856048ca8b675c` | REJECT / CITATION TRACE: Sherman (1993) constrained PHF. |
| `litcurate_92eff845f89656c1` | REJECT / CITATION TRACE: Sherman (1993) unconstrained PHF. |
| `litcurate_7055b8fa3b4f06f9` | REJECT / CITATION TRACE: Wang et al. experiment. |
| `litcurate_c2e7cca741b63d3b` | REJECT / CITATION TRACE: Mao et al. experiment. |
| `litcurate_098a976f2dafe671` | REJECT / NOT A VOLUME EOS: Karki MgSiO3 calculated elasticity. |
| `litcurate_db2a25c90398e1bf` | REJECT / CITATION TRACE / NOT A VOLUME EOS: MgSiO3 experiment. |
| `litcurate_19300ac5d1cb95e9` | REJECT / NOT A VOLUME EOS: Karki CaSiO3 calculated elasticity. |
| `litcurate_d69256693dfe5100` | REJECT / CITATION TRACE / NOT A VOLUME EOS: CaSiO3 experiment. |
| `litcurate_2ebe5e3ef9c993a5` | REJECT / NOT A VOLUME EOS: Karki MgO calculated elasticity. |
| `litcurate_8796ee52ed73dc76` | REJECT / CITATION TRACE / NOT A VOLUME EOS: MgO experiment. |
| `litcurate_258ef4fdbeb963b7` | REJECT / NOT A VOLUME EOS: Karki CaO calculated elasticity. |
| `litcurate_abc783c960b5b17d` | REJECT / CITATION TRACE / NOT A VOLUME EOS: CaO experiment. |
| `litcurate_9d5d64d73a5fafb8` | REJECT / NOT A VOLUME EOS: Karki stishovite calculated elasticity. |
| `litcurate_1dfe4b0e87d2ec69` | REJECT / CITATION TRACE / NOT A VOLUME EOS: stishovite experiment. |

## Reproduction artifacts

- `scripts/reproduce_oganov_wang_chizmeshya_eos.py` evaluates every accepted
  curve and performs all fits possible from tabulated primary data.
- `peritheos/data/datasets/mgo-oganov-2003-table3-pv.csv` preserves all eight
  Oganov Table-III checkpoints.
- `peritheos/data/datasets/ca-perovskite-wang-1996-table1-room-temperature.csv`
  preserves all fourteen Wang run-13 observations, including the two excluded
  low-pressure states and the source's inclusion decision.

No records are fabricated for missing raw grids, incomplete thermal
parameterizations, duplicate summaries, citation traces, or elastic moduli.
