# Shieh et al. (2006): natural-pyroxene pv and ppv

## Scope and primary source

This audit covers all nine LitCurate candidates attached to Shieh et al.,
*Equation of state of the postperovskite phase synthesized from a natural
(Mg,Fe)SiO3 orthopyroxene*, *Proceedings of the National Academy of Sciences*
**103**, 3039--3043 (2006),
[doi:10.1073/pnas.0506811103](https://doi.org/10.1073/pnas.0506811103).
The open-access primary article is the scientific authority. LitCurate was used
only for discovery and candidate enumeration.

The electron-microprobe composition of the starting orthopyroxene is
`(Mg1.80Fe0.18Al0.01Ca0.01)Si2O6`. Peritheos normalizes this exactly per SiO3 as
`Mg0.90Fe0.09Al0.005Ca0.005SiO3`; it does not round the minor Al and Ca away.
Laser heating produced coexisting Pbnm perovskite (pv, bridgmanite) and Cmcm
post-perovskite (ppv). The room-temperature volume observations were collected
on decompression without additional heating.

## Candidate disposition

| Candidate | Quantity represented | Decision |
|---|---|---|
| `litcurate_40925dc9160b7bdd` | This-study ppv second-order Eulerian finite-strain fit: `V0=164.9(6) A^3`, `K0=219(5) GPa`, implicit fixed `K0'=4` | **Accepted / preferred** as `mg090fe009al0005ca0005sio3_post_perovskite_shieh_2006_bm2_1`. |
| `litcurate_a84d5bbee1a84283` | This-study ppv sensitivity fit: `V0=166.2(7) A^3`, `K0=198(5) GPa`, fixed `K0'=4.4` | **Accepted / nonpreferred sensitivity** as `mg090fe009al0005ca0005sio3_post_perovskite_shieh_2006_bm3_sensitivity_2`. LitCurate's BM2 suggestion is corrected to BM3 because `K0'` is not 4. |
| `litcurate_c2542f0b2490af33` | This-study coexisting pv fit: measured `V0=163.3(1) A^3`, `K0=255(10) GPa`, fixed `K0'=3.7` | **Accepted** as `mg090fe009al0005ca0005sio3_bridgmanite_shieh_2006_bm3_1`. |
| `litcurate_b75920864fe802c4` | Oganov and Ono (2004) MgSiO3 ppv LDA values in comparison Table 1 | **Rejected as a Shieh-source record / citation trace only.** Audit the underlying primary theoretical paper rather than importing a secondary table row. |
| `litcurate_9a4875b670cf3c0e` | Oganov and Ono (2004) MgSiO3 ppv GGA values in comparison Table 1 | **Rejected as a Shieh-source record / citation trace only.** Same underlying-paper requirement. |
| `litcurate_a439638bf2273241` | Tsuchiya et al. (2005) MgSiO3 ppv LDA values in comparison Table 1 | **Rejected as a Shieh-source record / citation trace only.** Different composition and primary DFT publication. |
| `litcurate_d32fda0dff76cc71` | Fiquet et al. (2000) MgSiO3 pv experimental values in comparison Table 1 | **Rejected as a Shieh-source record / citation trace only.** The underlying experimental paper controls its data, pressure scale, and fit definition. |
| `litcurate_e7fdbac2dda84b1e` | Andrault et al. (2001) Mg0.95Fe0.05SiO3 pv experimental values in comparison Table 1 | **Rejected as a Shieh-source record / citation trace only.** Different composition and underlying publication. |
| `litcurate_766ca79517d1bb73` | Sinogeikin et al. (2004) MgSiO3 pv `K0=253(3) GPa` comparison value | **Rejected as incomplete / citation trace only.** The row has neither V0 nor K0' and is not a Shieh pressure-volume fit. |

The three accepted records are separate because they represent two crystal
structures and, for ppv, an explicitly nonpreferred parameter-tradeoff curve.
The six comparison rows are not counted as EOS additions from this paper.

## Structure and volume basis

Both phases use conventional orthorhombic cells with `Z=4`, so all published
volumes are stored directly in A^3 per conventional cell. Shieh et al. identify
Pbnm and Cmcm but do not publish composition-specific refined fractional
coordinates or trace-element site occupancies. The two material cards therefore
use explicitly labeled pure-MgSiO3 positional-topology proxies:

- Pbnm positions from the Yagi et al. (1978) MgSiO3 topology already bundled by
  Peritheos, with the lattice scaled to the measured pv `V0=163.3 A^3`;
- Cmcm positions from the Murakami et al. (2004) MgSiO3 topology already bundled
  by Peritheos, with the lattice scaled to the preferred ppv `V0=164.9 A^3`.

The proxy labels are deliberate: they support phase-specific diffraction-card
geometry without pretending that Fe, Al, and Ca partition or scattering factors
were refined. In particular, the ppv lattice is not described as a recovered
ambient structure; the paper states that a reliable low-pressure ppv volume
could not always be obtained.

## Figure 2 digitization

No numerical P-V table accompanies the article. The open-access PDF was
rendered at 400 dpi and Figure 2 was calibrated against its linear 20 GPa and
10 A^3 tick spacing. Filled circles were assigned to ppv and open circles to the
coexisting pv phase; theoretical/literature symbols and the ppv arrow labeled
“Theory” were excluded. The datasets contain:

- 13 separable ppv markers from 12.07 to 105.78 GPa; and
- 11 separable high-pressure pv markers from 82.91 to 105.78 GPa, plus the
  independently printed recovered volume `163.3(1) A^3` at ambient pressure.

The raster-derived points carry conservative uniform center-picking bounds of
0.50 GPa and 0.15 A^3. Those bounds are digitization uncertainty, not the
authors' experimental standard deviations. Some experimental error bars are
visible in the plot but cannot be transcribed consistently at the available
resolution, so the reproduction does not invent them. The dataset provenance
stores hashes of the source PDF and rendered figure crop together with the
pixel-axis calibration.

## Numerical reproduction

The deterministic script
`scripts/reproduce_shieh_2006_natural_pyroxene.py` evaluates the standard
Eulerian Birch--Murnaghan equation and refits unweighted pressure residuals.
The digitized figure does not recover the source's exact regression weights,
so parameter parity within printed uncertainty is the acceptance criterion.

| Record | Published-curve pressure RMSE (GPa) | Refit constraints | Independent refit | Assessment |
|---|---:|---|---|---|
| preferred ppv BM2 | 1.295150 | implicit `K0'=4`; fit V0, K0 | `V0=164.787751 A^3`, `K0=219.925996 GPa` | Both coefficients lie within the reported one-sigma intervals. |
| ppv BM3 sensitivity | 1.368471 | fixed `K0'=4.4`; fit V0, K0 | `V0=165.646948 A^3`, `K0=201.693422 GPa` | Both coefficients lie within the reported one-sigma intervals. |
| pv BM3 | 1.224195 | fixed measured `V0=163.3 A^3` and `K0'=3.7`; fit K0 | `K0=253.838363 GPa` | K0 lies within the reported one-sigma interval. |

The modest residuals include experimental scatter, figure center-picking, and
unknown source weights. The independent refits validate the printed curves but
do not replace them.

## Pressure calibration and validity

The three experimental configurations used internal Pt, Au, or NaCl standards,
citing Holmes et al. (1989), Shim et al. (2002), and Sata et al. (2002),
respectively. Argon or NaCl served as pressure medium/insulation depending on
the run. The article does not map each Figure 2 marker to one calibrant or print
row-wise calibrant volumes, so the reference scales are resolved but exact
observation-by-observation pressure recalculation is not possible.

The ppv records are bounded to the source-reported 12--106 GPa, 300 K
decompression range. Their persistence to low pressure is metastability, not a
thermodynamic stability claim. The pv fit is anchored by recovered V0 and only
high-pressure observations from 83--106 GPa; its intervening low-pressure curve
is correspondingly an interpolation across a sparsely observed interval.
