# Delta-project elemental-metal EOS

This addition contributes 50 static PBE third-order Birch-Murnaghan (BM3)
parameterizations for metals: the 48 metallic elements in the archived WIEN2k
reference table, plus independent FLEUR parameterizations for Al and Fe.  It
therefore contains 50 EOS records for 48 distinct metals.  The two duplicates
are retained as genuinely separate code implementations, not averaged or
silently selected.

## Primary sources and archived evidence

The third-order Birch-Murnaghan equation and original seven-volume benchmark
protocol are documented by Lejaeghere, Van Speybroeck, Van Oost, and
Cottenier, *Critical Reviews in Solid State and Materials Sciences* **39**,
1-24 (2014),
[doi:10.1080/10408436.2013.772503](https://doi.org/10.1080/10408436.2013.772503).
Section IV.A describes the seven-point elemental-crystal EOS workflow.  An
[author-deposited manuscript](https://backoffice.biblio.ugent.be/download/4188194/4188200)
was used for the page-level method audit. This 2014 paper used an earlier
WIEN2k 11.1 coefficient set; its DOI is therefore retained only as
methodological lineage, not attributed to the archived WIEN2k 13.1 records.

Both the WIEN2k 13.1 and FLEUR 0.26 results belong to the multi-code
reproducibility campaign of Lejaeghere et al., *Science* **351**, aad3000
(2016),
[doi:10.1126/science.aad3000](https://doi.org/10.1126/science.aad3000).
The original Delta-project website and supplementary files are preserved in
the CC BY 4.0 [Materials Cloud archive](https://archive.materialscloud.org/record/2023.133),
[dataset doi:10.24435/materialscloud:5e-mv](https://doi.org/10.24435/materialscloud:5e-mv).

The checked files are `README.txt`, `WIEN2k.txt`,
`history/history/FLEUR.txt`, `eosfit.py`, and the element CIFs in
`CIFs.tar.gz`. Their SHA-256 digests are preserved in
`peritheos/data/datasets/deltacodesdft-metal-eos-source.json`. The exact
selected rows and cell-volume conversions are in
`peritheos/data/datasets/deltacodesdft-metal-eos-parameters.csv`.

## Equation and units

The official `eosfit.py` fits a third-order Birch-Murnaghan energy-volume
curve and reports (V_0) in Å³/atom, (B_0) in GPa, and (B'_0) as a
dimensionless derivative.  Peritheos evaluates the corresponding pressure:

\[
P(V)=\frac{3B_0}{2}(\eta^7-\eta^5)
\left[1+\frac{3}{4}(B'_0-4)(\eta^2-1)\right],\qquad
\eta=(V_0/V)^{1/3}.
\]

Each archived atomic (V_0) is multiplied by the target material card's
formula-unit count.  This changes only the public volume basis: the pressure
curve is invariant.  The CSV records both the source atomic value and the
executable cell value, making every conversion auditable.

The workflow uses seven fixed geometries made by scaling the archived CIF to
94%, 96%, 98%, 100%, 102%, 104%, and 106% of its volume.  The cell shape and
internal coordinates remain frozen.  The stored validity interval is that
exact source-CIF envelope expressed as (V/V_0), rather than an assumed
phase-stability range.  Searchable pressure bounds are derived by evaluating
the stored BM3 curve at those two volume endpoints and are labeled theoretical,
not experimental.  All records are static 0 K PBE models and omit spin-orbit
coupling. The WIEN2k rows use version 13.1, all-electron APW+lo; the FLEUR
rows use version 0.26, all-electron LAPW (+lo). Both are scalar-relativistic.
The source requires antiferromagnetic spin polarization for Cr and Mn, and
ferromagnetic spin polarization for Fe, Co, and Ni; these states are explicit
in each affected EOS record.

## Selection and phase handling

The WIEN2k selection is:

| Elements | Count |
|---|---:|
| Li, Be, Na, Mg, Al, K, Ca, Sc, Ti, V, Cr, Mn | 12 |
| Fe, Co, Ni, Cu, Zn, Ga, Rb, Sr, Y, Zr, Nb, Mo | 12 |
| Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Cs, Ba, Lu, Hf | 12 |
| Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po | 12 |
| FLEUR Al and Fe cross-code parameterizations | 2 |
| **Total EOS records** | **50** |

Twenty-three exact composition/phase matches extend existing material cards.
Twenty-five structures absent from the catalog receive new cards.  Those new
cards retain the official archive's explicit P1 representation, complete
atomic positions, and source cell, so diffraction and volume conversion do
not depend on inferred symmetry operations.

Three assignments require special care:

- Li and Na use the source's three-atom rhombohedral 9R cells and therefore do
  not extend the ambient bcc/legacy material cards.
- Mn is the source's deliberately simplified antiferromagnetic fcc model, not
  complex alpha-Mn.  It is a distinct `manganese_fcc_afm` material.
- The WIEN2k Mn fit reports (B'_0=-0.21).  The value is transcribed exactly,
  flagged prominently, and restricted to the seven-volume fitting envelope;
  it is unsuitable for broad extrapolation.

## Numerical reproduction

`scripts/reproduce_deltacodesdft_metal_eos.py`
reads the independent CSV transcription, verifies all 50 stored coefficient
triplets and cell conversions, evaluates a separately written BM3 pressure
expression at (V/V_0=0.94), 1, and 1.06, and compares it with the Peritheos
implementation.  The current result is:

```text
records: 50
elements: 48
wien2k_records: 48
fleur_records: 2
largest_parameter_error: 0.0
largest_p0_gpa: 0.0
largest_independent_pressure_error_gpa: 3.055333763768431e-13
```

The archive does not preserve the seven row-level E(V) values for these code
tables, so a direct coefficient refit is unavailable.  No synthetic energy
rows, uncertainties, covariance, or experimental pressure range are claimed.
