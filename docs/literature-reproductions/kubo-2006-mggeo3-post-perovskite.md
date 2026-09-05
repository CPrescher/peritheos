# Kubo et al. (2006): MgGeO3 post-perovskite

## Primary source and scope

Atsushi Kubo, Boris Kiefer, Guoyin Shen, Vitali B. Prakapenka, Robert J.
Cava, and Thomas S. Duffy, *Stability and equation of state of the
post-perovskite phase in MgGeO3 to 2 Mbar*, **Geophysical Research Letters
33**, L12S12 (2006),
[doi:10.1029/2006GL025686](https://doi.org/10.1029/2006GL025686).
The audit used the [author-hosted final article](https://duffy.princeton.edu/sites/g/files/toruqf616/files/media/kubo-2006-georeslet.pdf)
and the publisher's plain-text
[Supporting Information Table S1](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2006GL025686&file=grl21266-sup-0002-ts01.txt).

The phase is pure MgGeO3 post-perovskite in the CaIrO3-type `Cmcm`
structure. The conventional C-centered orthorhombic cell contains four
formula units (`Z=4`). The material card's diffraction structure is the
source's 87.8 GPa Rietveld refinement; it is not an ambient structure.

The room-temperature EOS observations came from three laser-heated
diamond-anvil-cell runs. Pt was both absorber and pressure standard. Pressures
were reduced with the Holmes et al. (1989) Pt EOS using only the Pt 111
reflection to reduce nonhydrostatic bias. Table S1 contains 25 exact rows with
`a`, `b`, and `c` and their 1-sigma errors, Pt 111 and 200 spacings, reduced
pressure, differential stress, and exposure timing.

## Fit selection and numerical reproduction

The source used the standard third-order Birch-Murnaghan equation. Because no
reliable observations were available near ambient pressure, it also invoked
the semi-empirical `K*V=constant` systematics for a common crystal structure
using MgSiO3 post-perovskite references. The selected coefficients are

| coefficient | selected value | status |
|---|---:|---|
| `V0` | 179.2(7) A^3/cell | fitted/extrapolated |
| `K0` | 207(5) GPa | fitted with the source constraint |
| `K0'` | 4.4 | fixed |

The paper excludes its three rows below 47 GPa because the weakened and
broadened diffraction peaks yielded volumes above the smooth high-pressure
trend. Those 45.0, 36.3, and 7.4 GPa rows remain in the bundled dataset with
`fit_included=0`. The selected fit therefore uses 22 rows over
47.3--196.3 GPa at 300 K.

Recomputing every cell volume as `V=a*b*c` and evaluating the printed curve
gives a pressure RMSE of 1.394506953 GPa and a maximum absolute residual of
3.937891209 GPa. An independent unweighted pressure-residual refit with
`K0'=4.4` fixed returns `V0=179.20054626 A^3` and `K0=206.71122613 GPa`.
Both coefficients reproduce the published values well inside their reported
uncertainties. The reported Pt differential stress is not a pressure standard
deviation and was not used as one.

Paragraph 17 also reports the same-data sensitivity result
`V0=175.9(6) A^3`, `K0=245(5) GPa`, and fixed `K0'=4`. Although the prose
retains the third-order name, `K0'=4` removes the third-order correction, so
Peritheos stores this curve canonically as BM2. It has a pressure RMSE of
1.356128559 GPa; an independent BM2 refit returns `V0=175.92371146 A^3` and
`K0=244.61773418 GPa`. Its slightly smaller unweighted RMSE does not make it
the source-selected result: the paper presents it specifically to demonstrate
the large `V0`-`K0`-`K0'` trade-off. It is executable but non-default.

Recalculating all 25 pressures from the rounded Pt 111 spacings with the
bundled `platinum_holmes_1989_vinet_1` record gives a 0.253508148 GPa RMSE and
a 0.445212011 GPa maximum difference from the printed pressures. This confirms
the calibration lineage to the precision permitted by the four-decimal-place
Pt spacings.

No digitization was needed. The machine-readable publisher table is complete.

## Disposition of every same-DOI LitCurate row

The source DOI generated eight LitCurate candidates. Only two are executable
experimental curves from this paper. The other six are either the paper's own
unverifiable DFT summary or comparison values copied from other primary
publications.

| LitCurate identifier | reported row | disposition | reason |
|---|---|---|---|
| `litcurate_257ec8acbb2a2877` | MgGeO3 experiment: `V0=179.2`, `K0=207`, fixed `K0'=4.4` | **accepted, preferred/default** as `mggeo3_post_perovskite_kubo_2006_bm3_1` | Source-selected constrained fit; all 22 selected observations reproduce it. |
| `litcurate_973395417d1d42d7` | MgGeO3 experiment: `V0=175.9`, `K0=245`, fixed `K0'=4` | **accepted, nonpreferred sensitivity** as `mggeo3_post_perovskite_kubo_2006_bm2_sensitivity_2` | Reproduces the same 22 rows, but the source uses it only to demonstrate parameter trade-off. |
| `litcurate_24cd91aaacc9abd7` | This study's static LDA: `V0=178.02`, `K0=201.9`, `K0'=4.34` at 0 K | **held / not refittable** | Table 2 reports the coefficients, but neither the article nor its supplement publishes the DFT energy-volume or pressure-volume grid. It may become a separate theoretical record, never an experimental alternative, if parameter-only theory records are intentionally added. |
| `litcurate_49869576e372eb74` | Hirose et al. (2005) MgGeO3 experiment: `V0=183.1(8)`, `K0=192(5)`, `K0'=4` | **rejected as a Kubo-source record; primary-source audit required** | Kubo Table 2 explicitly identifies this as Hirose et al.'s earlier experiment. It belongs to DOI `10.2138/am.2005.1702`; importing it from the Kubo comparison table would break source lineage. |
| `litcurate_21f4bb278eb404db` | Tsuchiya et al. (2004) MgSiO3 DFT: `V0=163.81(5)`, `K0=222(1)`, `K0'=4.2(1)` | **rejected as a Kubo-source record; citation trace only** | Different composition and underlying primary DFT publication. |
| `litcurate_0a9cb64914d6f525` | Oganov and Ono (2004) MgSiO3 DFT: `V0=162.86`, `K0=231.9`, `K0'=4.43` | **rejected as a Kubo-source record; citation trace only** | Different composition and underlying primary DFT publication. |
| `litcurate_89524810ffc6d316` | Shieh et al. (2006) (Mg,Fe)SiO3 experiment: `V0=164.9(6)`, `K0=219(5)`, `K0'=4` | **rejected as a Kubo-source record; citation trace only** | Different composition and sample; its primary source is DOI `10.1073/pnas.0506811103`. |
| `litcurate_7000ad6e95c15aea` | Ono et al. (2006) MgSiO3 experiment: `V0=162.86`, `K0=225`, `K0'=4` | **rejected as incomplete and duplicate** | Kubo Table 2 actually reports `K0=225--249`, not one universal 225 GPa value. The primary paper is already represented by the three pressure-scale-specific `mgsio3_post_perovskite_ono_2006_*` records. |

The preliminary 2005 IUCr conference-abstract candidate
`litcurate_84fcb601deba874d` (`V0=179.7(9) A^3`, `K0=203(6) GPa`,
`K0'=4.4`) is not the same DOI, but it is the same research program and was
superseded by this full paper's constrained result. It should not be added as a
second production EOS without an explicit historical reason.

## Other reported source fits that are not material EOS records

The paper's own 0 K LDA fit is retained in the disposition table above. The
underlying calculated states are unavailable, so it cannot be independently
refitted without digitizing the sparse blue points in Figure 3; such a
digitization would be a low-information curve check rather than recovery of the
authors' calculation grid.

Paragraph 19 additionally fits each lattice axis with a modified axial
Birch-Murnaghan form: `a0=2.803(2) A`, `Ka0=230(4) GPa`, `Ka0'=4.6`;
`b0=9.292(27) A`, `Kb0=161(7) GPa`, `Kb0'=4.0`; and `c0=6.882(5) A`,
`Kc0=247(3) GPa`, `Kc0'=4.7`. These are directional incompressibility models,
not volumetric material EOS records, and Peritheos does not promote them into
the equilibrium EOS catalog. The article appears to print cubic-angstrom units
for two zero-pressure axis lengths; dimensional consistency requires angstroms.

## Reproduction command

```bash
python scripts/reproduce_kubo_2006_mggeo3_post_perovskite.py
```
