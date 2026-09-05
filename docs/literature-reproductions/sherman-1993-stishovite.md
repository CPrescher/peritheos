# Sherman (1993): stishovite and hypothetical post-stishovite SiO2

## Primary source and audit boundary

This audit covers all 12 LitCurate rows grouped under David M. Sherman,
*Equation of state and high-pressure phase transitions of stishovite (SiO2):
Ab initio (periodic Hartree-Fock) results*, **Journal of Geophysical Research:
Solid Earth 98**, 11865--11873 (1993),
[doi:10.1029/93JB00783](https://doi.org/10.1029/93JB00783).

The accessible publisher page exposes the primary abstract and reference list,
but not the article body or Figure 4. OpenAlex and Semantic Scholar report the
article as closed and identify no repository full text, and a title/DOI search
of the local Zotero library found no copy. No supplement is registered. The
raw LitCurate workbook from its
[Zenodo deposit](https://doi.org/10.5281/zenodo.22118629) was inspected to
understand the extraction, but it remains discovery evidence rather than
primary scientific authority.

Consequently this paper contributes no production EOS record in the present
batch. The source-reported fits are documented as holds, and values copied by
Sherman from older publications remain citation traces.

## Source-reported stishovite fit

The primary abstract identifies the material as SiO2 stishovite and the method
as a static periodic Hartree-Fock calculation. It states that the calculated
total energies were fit with a third-order Birch--Murnaghan equation and
reports

| coefficient | primary abstract value |
|---|---:|
| `V0` | 46.1 A^3 |
| `K0` | 328 GPa |
| `K0'` | 4.0 |

The approximately 46 A^3 volume is the conventional rutile-type stishovite
cell, `P42/mnm`, `Z=2`, matching the existing `sio2_stv_andr` material. This is
a zero-Kelvin theoretical fit, so an experimental pressure calibration does
not apply. The abstract further says that compression to 80% of the
zero-pressure volume corresponds to a pressure near 110 GPa.

For the standard BM3 pressure form,

\[
P(V)=\frac{3K_0}{2}
\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_0-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\},
\]

the abstract coefficients give 114.466551 GPa at `V/V0=0.8`. This is
consistent with the qualitative phrase "near 110 GPa" and confirms that the
reported coefficients map plausibly to the standard BM3 convention. It is a
curve spot check, not a coefficient refit.

LitCurate's evidence field says that the article body and Figure 4 instead
report `V0=46.2(1) A^3`, `K0=328(5) GPa`, and `K0'=4.0(5)`. That creates a
46.1-versus-46.2 A^3 abstract/body discrepancy. Because neither the body nor
the energy-volume points could be recovered from a primary copy, this audit
does not select one value silently, promote the secondary extraction's errors,
or synthesize an observation grid. Candidate
`litcurate_6e1a8d4040649782` is therefore **held / not refittable**. It may be
revisited if the article body or an author-provided calculation grid becomes
available.

## Source-reported modified-fluorite fit

The primary abstract says that Sherman also calculated hypothetical
modified-fluorite and alpha-PbO2 structures. It reports only qualitative
energetic and density conclusions for them, not a complete EOS
parameterization.

LitCurate candidate `litcurate_55f140f68d70d27a` reports a source fit with
`V0=42.9 A^3`, `K0=387 GPa`, and `K0'=4.4`, plus an internal energy 131 kJ/mol
above stishovite. Its evidence describes the volume as a cell containing two
SiO2 formula units, but the candidate leaves the phase and structure
unspecified. The likely modified-fluorite or pyrite-type `Pa-3` conventional
cell contains four formula units, which would require an explicit conversion
to 85.8 A^3 for Peritheos's conventional-cell volume convention. The
inaccessible article body is needed to establish whether 42.9 A^3 is a
two-formula-unit normalization, a computational primitive-cell convention, or
something else, and to recover the structural coordinates and calculated
states.

This candidate is **held / unresolved cell basis**. The alpha-PbO2 discussion
does not itself form another candidate: the accessible primary abstract gives
only an energy difference within 10 kJ/mol and a density about 2% greater than
stishovite, not a complete executable EOS.

## Disposition of every same-DOI LitCurate row

| LitCurate identifier | reported row | disposition | reason |
|---|---|---|---|
| `litcurate_6e1a8d4040649782` | This work, rutile stishovite HF-LCAO BM3: `V0=46.2`, `K0=328`, `K0'=4.0` | **held / not refittable** | Figure 4's calculated energy-volume grid is unavailable, and the accessible primary abstract instead prints `V0=46.1 A^3`. |
| `litcurate_55f140f68d70d27a` | This work, unspecified hypothetical SiO2 BM3: `V0=42.9`, `K0=387`, `K0'=4.4` | **held / unresolved cell basis** | The primary abstract does not print these coefficients; the likely modified-fluorite structure and the reported two-formula-unit volume do not establish a conventional-cell record without the article body. |
| `litcurate_2e1a2869e78573c6` | Bassett and Barnett (1970): `K0=300`, `K0'=4` | **citation trace only** | Earlier experiment; `V0` and the equation definition are absent from Sherman's comparison row. Audit the underlying paper, DOI `10.1016/0031-9201(70)90044-0`. |
| `litcurate_d62422e996bba1c5` | Liu et al. (1974): `K0=344(27)`, `K0'=2--7` | **citation trace only** | Earlier nonhydrostatic experiment; LitCurate collapses the derivative range to its lower endpoint. Audit DOI `10.1029/JB079i008p01160`. |
| `litcurate_058788c236ed797e` | Sato (1977): `K0=298`, `K0'=0.7` | **citation trace only** | Earlier hydrostatic-compression study; audit DOI `10.1016/0012-821X(77)90015-2`. |
| `litcurate_2987efb3bdee53f4` | Weidner et al. (1982): `K0=306`, no `K0'` | **citation trace only** | Brillouin elasticity result, not a complete Sherman EOS. Audit DOI `10.1029/JB087iB06p04740`. |
| `litcurate_7426a266150e18e2` | Sugiyama et al. (1987): `V0=46.591`, `K0=313`, `K0'=6.0` | **citation trace only** | Earlier diffraction fit; audit DOI `10.2465/minerj.13.455`. |
| `litcurate_511d9d47e2806b38` | Tsuchida and Yagi (1989): `K0=376`, `K0'=4.0` | **citation trace only** | Earlier nonhydrostatic experiment; `V0` is absent from the comparison row. Audit DOI `10.1038/340217a0`. |
| `litcurate_a69e1118096651c2` | Ross et al. (1990): `V0=46.615`, `K0=313`, `K0'=2.8` | **citation trace only** | Earlier experimental publication; source lineage must remain with Ross et al., not Sherman. |
| `litcurate_3d0214da7860bef3` | Park et al. (1988): `V0=47.5`, `K0=288`, `K0'=3.14` | **citation trace only** | Earlier LAPW calculation; audit DOI `10.1038/336670a0`. |
| `litcurate_9067bba915292c59` | Cohen (1991): `V0=46.16`, `K0=324`, `K0'=4.04` | **citation trace only** | Earlier LAPW calculation; audit the underlying *American Mineralogist* article. |
| `litcurate_47b4ad839b80c80e` | Keskar et al. (1991): `V0=45.64`, `K0=292`, `K0'=5.86` | **citation trace only** | Earlier plane-wave pseudopotential calculation; audit DOI `10.1103/PhysRevB.44.4081`. |

## Duplicate and future-work check

No bundled record uses DOI `10.1029/93JB00783`, and no SiO2 record duplicates
either source-reported coefficient set. The existing stishovite material
contains independently sourced experimental and thermal EOS records, so a
future Sherman static-theory record would be a distinct source
parameterization rather than a duplicate. It should be reconsidered only when
the primary energy-volume curve or an equivalent independently specified
source benchmark is recoverable.
