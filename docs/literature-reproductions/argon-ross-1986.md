# Ross et al. (1986): published tables recovered; Monte Carlo reproduction pending

M. Ross, H. K. Mao, P. M. Bell and J. A. Xu,
*The equation of state of dense argon: A comparison of shock and static studies*,
J. Chem. Phys. **85**(2), 1028-1033 (1986),
[doi:10.1063/1.451346](https://doi.org/10.1063/1.451346).

**Primary full text verified.** The earlier access blocker is resolved by the
user-supplied journal PDF. The `not_reproduced` status means that an independent
finite-temperature Monte Carlo EOS has not yet been implemented and validated;
it does not mean the published work has been shown to be irreproducible.
No executable Ross EOS or substitute Birch-Murnaghan/Vinet fit is registered.

## Verified primary source and Zotero

The original PDF is `/Users/clemens/Downloads/1028_1_online.pdf`; the stable copy
is `/Users/clemens/Documents/Peritheos-sources/argon-ross-1986/ross-1986-jcp.pdf`.
It has seven PDF pages: publisher cover plus journal pages 1028-1033, 869034 bytes,
SHA256 `8f5f50ab785adb5f03c3d5769d1c3b91a0ec7c3e5eb3a1374faa273112aaf9c3`.
Title, authors, DOI, tables and page numbers were checked against the PDF.
The file is outside git. Local paths and checksums are in the
[handoff manifest](../data/argon-ross-1986-handoff.json).

Zotero: **My Library > Methods > EOS Library**, collection `JT8V6LUL`,
parent `QWLP4V5H`, PDF `8FAK8WKL`. The stored attachment's SHA256 matches the
source. Older metadata-only item `TJWGXSZR` remains a duplicate; its historical
access-blocker note is obsolete. No merge or deletion was performed.

## Journal datasets

| Dataset | Source and interpretation | Temperature |
| --- | --- | --- |
| `argon_ross_1986_table1` | All 42 Table I static measurements, page 1029; pressure, parenthetical error, lattice edge and molar volume | 298 K |
| `argon_ross_1986_table3` | All 16 Table III published calculated isotherm states, page 1030; pressure and Gruneisen gamma; alpha=13.2 | Isotherm 293 K; gamma 298 K per prose |
| `argon_ross_1986_table2` | All six Table II calculated liquid Hugoniot states, including initial state; alpha=13.2 | 87-11963 K |

The three CSV assets are in `peritheos/data/datasets/argon-ross-1986-table*.csv`.
Each has a checksum, DOI and source-PDF checksum. Theory datasets are explicitly
marked theoretical and are not experimental observations. All have empty
`used_by_eos_records`; none were used to fit Dewaele or Ono.

Pressure conversion is 1 kbar = 0.1 GPa. Source volumes are cm3/mol of atoms.
`volume_a3 = Vmolar * 4e24 / NA` uses four atoms per conventional fcc cell.
For the liquid only this is an equivalent four-atom normalization, not a crystal
cell. Printed molar volumes determine the converted volumes; they are not
recomputed from rounded lattice edges. Extra converted decimals preserve
arithmetic and do not imply experimental precision.

Table I covers 1.6-80.6 GPa. Parenthetical pressure errors are preserved as
**generic uncertainties**, because the paper does not specify a confidence
level. No volume uncertainty is invented. Repeated and anomalous rows remain:

- Row 17 prints 247(13) kbar, a=4.047 A, V=9.98 cm3/mol, a nonmonotonic pressure.
- Row 32 prints a=3.685 A at 568(17) kbar and V=8.83 cm3/mol. The lattice edge
  is inconsistent with the printed volume. Both printed values are preserved.

These flags identify conspicuous source issues, not an exhaustive error model.
Table I temperature is 298 K; Table III's heading and the discussion on page
1032 specify 293 K for the calculated isotherm. Page 1030 explicitly gives 298 K
for gamma. Separate columns preserve these source statements without silently
harmonizing them. Table II initial pressure is 0 as printed, although the prose
also reports a model initial pressure of -28 bar for alpha=13.2 including quantum
corrections. The table and prose values are not silently substituted.

The ruby pressure scale, Eq. (1), is
`P(Mbar)=19.04/7.665*((1+delta_lambda/lambda0)^7.665-1)`.
Original ruby wavelengths are not provided, so reported pressures are retained.

## Model availability and reproduction boundary

Journal Eq. (2) prints the complete effective exp-6 pair potential:

```text
phi(r)/kB = 122/(alpha-6) *
            [6*exp(alpha*(1-r/3.85)) - alpha*(3.85/r)^6].
```

Distances are A; epsilon/kB=122 K, r*=3.85 A. The preferred alpha is 13.2;
13.0 is a comparison, not a statistical parameter uncertainty. The audit script
implements this pair energy on the discussed r>=2 A interval and verifies its
minimum. The short-distance exp-6 catastrophe is outside that audit domain.

Pressure and energy were calculated by Monte Carlo; Eq. (3) defines the shock
Hugoniot relation and the temperature iteration is described on page 1030.
Evaluating a pair energy alone does not provide the finite-temperature pressure
EOS. A full independent reproduction still requires implementing the statistical
mechanics calculation and resolving numerical conventions, including the cited
methodology and quantum corrections. This work does not claim those steps were
completed or that the study is irreproducible.

Table III can already be plotted directly. Interpolation between its states
would represent the published tabulation, not an independently regenerated Monte
Carlo calculation. A conventional fit to Table I would be a new fit. Neither is
added here. Table III's high pressures are theoretical extrapolation beyond the
80.6 GPa measured range; the paper proposes use as a standard to about 3-4 Mbar.

## Distinct 1985 precursor retained

The [UCRL-93030 conference precursor](https://www.osti.gov/biblio/5471606)
remains a separate source, not the journal article. Its PDF checksum is
`fbc357ed1cd980fd495f26b22932adb3540d3d36d9c10349b32899fc42800302`.
`argon_ross_1985_precursor_figure1_subset` retains the previous 19 digitized
293 K static triangles, about 1.81-57.66 GPa, with their original pixel
coordinates and calibration. The subset is incomplete and has no journal DOI.
Its conservative coordinate-reading bounds are not experimental uncertainties.
Do not relabel these points 298 K or replace their version provenance.

## Validation and display

Run `python -m scripts.reproduce_argon_ross` to regenerate the
[audit report](../data/argon-ross-1986-reproduction.json). It checks table counts,
unit conversions, precursor digitization arithmetic and pair potential values;
fit residuals remain null. Tests also check dataset loading, schema, temperature
separation, printed anomalies and generic uncertainty semantics.

Studio should show primary source verified and published tables recovered,
with independent Monte Carlo reproduction pending. It should distinguish
measurements, calculated solid isotherm, calculated liquid Hugoniot and the older
precursor subset. No Ross executable EOS curve is supplied by this change.
