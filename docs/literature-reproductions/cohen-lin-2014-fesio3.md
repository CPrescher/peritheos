# Cohen and Lin (2014): static FeSiO3 Pv, PPv, and PPv-II

## Outcome

All three source-reported LitCurate candidates are accepted as distinct static
Vinet records. They describe three different FeSiO3 structures, not repeated
fits or experimental-run splits. The records are stored in phase-specific
material cards:

- `fesio3_bridgmanite_cohen_lin_2014_vinet_1`
- `fesio3_post_perovskite_cohen_lin_2014_vinet_1`
- `fesio3_post_perovskite_ii_cohen_lin_2014_vinet_1`

The disposition after source-exhaustion review is **retain all three published
parameterizations, but retain `not_refittable` for all three**. No authoritative
row-level energy-volume data were found, and a graphical extraction cannot
reconstruct the source's total energies, energy zero, weights, or E(V) fitting
protocol.

## Primary source

R. E. Cohen and Y. Lin, “Prediction of a potential high-pressure structure of
FeSiO3,” *Physical Review B* **90**, 140102(R) (2014),
doi:[10.1103/PhysRevB.90.140102](https://doi.org/10.1103/PhysRevB.90.140102).
The audit used the [authoritative UCL repository
copy](https://discovery.ucl.ac.uk/id/eprint/1461221/1/PhysRevB.90.140102.pdf),
SHA-256 `c4bd0265b83f2859c3bb2158f91dd20be064771a03cabbcb6bb164ad4ea5a711`.

The calculations are static PAW-PBE GGA+U (`U=6.0 eV`) Quantum ESPRESSO
calculations in 20-atom cells, using a 4x4x4 k-point mesh and 80 Ry cutoff. The
paper states that all three phases use high-spin antiferromagnetic states.

## Source-exhaustion audit

The following independent primary-source surfaces were inspected on 2026-09-08:

- the live APS version-of-record page, which exposes the article PDF and no
  supplemental-material or data attachment;
- Crossref metadata, which lists only APS version-of-record links and no related
  object, and OpenAlex metadata, which lists only APS and arXiv locations and no
  dataset;
- the UCL Discovery record, which contains exactly one file, the 708 kB article
  PDF cited above;
- the Carnegie Science publication record, which is bibliographic only;
- a DataCite related-identifier query, which returns only the arXiv DOI
  `10.48550/arXiv.1407.0361`, and an exact-DOI Zenodo query, which returns zero
  records;
- both arXiv submissions for [`1407.0361`](https://arxiv.org/abs/1407.0361).
  Version 1 is an author-generated PDF (SHA-256
  `368f7774c0ec004e1c99f9100a927c1559998eb2849047c7856277470430378c`).
  The version 2 source archive (SHA-256
  `460ed6935a6cadcd97bdf4fddfb5b74fd7593cd9fea631dc93d020c0611ec4e0`)
  contains one TeX file, one BibTeX output file, and eight EPS figures—no table,
  calculation input/output, spreadsheet, or other numerical supplement.

The APS, UCL, arXiv v1, and arXiv v2 PDFs each contain zero embedded file
attachments.

The `Figures_03.eps` asset is not a vector graph with recoverable plotted
coordinates. Its header identifies `bmeps 1.2.7`, followed by a single
ASCII85/run-length encoded 1606 x 978 RGB raster (`colorimage`); its SHA-256 is
`d1780680166a0234ffbfe75d8cb04eb0e7708efedcc921459e12c7c1172a0372`.
Exact source rows therefore remain unavailable after checking the publisher,
the two institutional records, both author submissions, DOI relations, and
research-data indexes.

## EOS basis and reproduction

Table III explicitly labels the equation Vinet and prints `V0`, `K0`, and
`K0'` for Pv, PPv, and PPv-II. Its volumes are per FeSiO3. The material cards
use conventional 20-atom cells (`Z=4`), so only volume is multiplied by four;
bulk modulus and its derivatives are unchanged.

The underlying eight static energies correspond to -10, 0, 25, 50, 75, 100,
125, and 150 GPa and appear only in Figure 3. Digitizing that small energy plot
would add false precision. More usefully, Table III provides three independent
100 GPa checks for every fit: `V100`, `K100`, and `K100'`.

| Phase | V0 (A3/formula) | K0 (GPa) | K0' | Source V100 | Reproduced V100 | Source K100 | Reproduced K100 | Source K100' | Reproduced K100' |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pv | 44.31 | 225 | 4.42 | 34.27 | 34.28075 | 597 | 596.8896 | 3.34 | 3.33886 |
| PPv | 44.90 | 189 | 4.73 | 33.98 | 33.94965 | 579 | 577.8706 | 3.47 | 3.46808 |
| PPv-II | 45.45 | 195 | 4.67 | 34.49 | 34.50403 | 580 | 580.5804 | 3.44 | 3.44258 |

Across the nine final-table checks, the maximum differences are 0.031
A3/FeSiO3 in volume, 1.13 GPa (0.20%) in modulus, and 0.0026 in `K'`. These
small nonzero offsets reflect the coarse two-decimal or integer coefficients
printed in final Table III; they should not be described as exact rounding
parity. `scripts/reproduce_cohen_lin_2014_fesio3.py` evaluates the Vinet
equation independently and performs these checks deterministically.

That rounding interpretation has a direct primary-source cross-check. The
superseded arXiv v1 Table 2 prints `K0=224.8`, `189.7`, and `194.7 GPa`, and
`K100=596.5`, `579.0`, and `580.1 GPa` for Pv, PPv, and PPv-II. Replaying those
extra-decimal coefficients gives, respectively:

| Phase | Replayed V100 | v1 V100 | Replayed K100 | v1 K100 | Replayed K100' | v1 K100' |
|---|---:|---:|---:|---:|---:|---:|
| Pv | 34.27564 | 34.27 | 596.6558 | 596.5 | 3.33845 | 3.34 |
| PPv | 33.97136 | 33.98 | 578.7261 | 579.0 | 3.46988 | 3.47 |
| PPv-II | 34.49488 | 34.49 | 580.2168 | 580.1 | 3.44184 | 3.44 |

Version 1 is used only to explain final-publication rounding; it does not
replace the final Table III coefficients in the production records.

## Approximate plot diagnostic

For completeness, marker centers were read from the volume axis of the arXiv
v2 Figure 3 raster. This is explicitly an **approximate plot diagnostic, not a
primary observation dataset or an independent refit**. Only marker x positions
were used; the plotted energies were not promoted to numerical data. At the
nominal 100 GPa point the approximate marker volumes are 34.2267, 33.9255, and
34.4562 A3/FeSiO3 for Pv, PPv, and PPv-II, within 0.055 A3/FeSiO3 of the three
Table III values. Across all eight pressure labels, the largest marker-to-
published-curve volume difference is 0.217 A3/FeSiO3. This is adequate as a
visual transcription check and inadequate as source-faithful E(V) evidence.

## Structure and scope

Table I supplies complete 100 GPa Cmcm PPv and Cmmm PPv-II lattices and
fractional coordinates, which are stored without modification. It does not
tabulate the Pv structure, so no Pv lattice, space-group setting, or coordinates
are invented. Table I site multiplicities establish `Z=4` for both tabulated
cells. The printed lattice products give 33.9343 and 34.4492 A3/FeSiO3,
respectively, within 0.046 A3/FeSiO3 of the independent Table III 100 GPa
volumes after the three-decimal lattice rounding.

The mathematical fit spans the eight static calculation pressures, but this is
not a phase-stability interval. In particular, the authors find that Cmmm
PPv-II distorts to C2/m at low pressure. The Cmmm structure card is therefore
anchored to the explicitly tabulated 100 GPa state.

## Candidate dispositions

| LitCurate identifier | Disposition | Reason |
|---|---|---|
| `litcurate_9b4780a09d66d440` | ACCEPT | Source Table III Pv Vinet; independent V100/K100/K100' agreement within printed-coefficient precision. |
| `litcurate_ceaea2f0eec4299a` | ACCEPT | Distinct source Table III PPv Vinet; Cmcm structure tabulated. |
| `litcurate_7596b5eb1cb53baf` | ACCEPT | Distinct source Table III PPv-II Vinet; Cmmm structure tabulated. |

There are no citation-reported comparison rows under this DOI in the current
LitCurate ledger.

No removal is recommended. Each record is a directly published, equation-
identified parameterization with multiple independent table and structure
checks. Conversely, no status upgrade is warranted: a source-faithful refit
requires authoritative numerical E(V) rows (or author-supplied calculation
outputs), not the raster diagnostic above.
