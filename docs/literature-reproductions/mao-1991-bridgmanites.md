# Mao et al. (1991): (Fe,Mg)SiO3 perovskites

Primary source: H. K. Mao, R. J. Hemley, Y. Fei, J. F. Shu, L. C. Chen,
A. P. Jephcoat, Y. Wu, and W. A. Bassett, “Effect of pressure, temperature,
and composition on lattice parameters and density of (Fe,Mg)SiO3-perovskites
to 30 GPa,” *Journal of Geophysical Research* **96**, 8069–8079 (1991),
<https://doi.org/10.1029/91JB00176>.

The final article reports indistinguishable room-temperature compression for
three compositions. Table 1 gives distinct conventional-cell reference volumes
of `162.49(7)`, `162.79(8)`, and `163.53(10) Å3` for MgSiO3,
Mg0.9Fe0.1SiO3, and Mg0.8Fe0.2SiO3. Crucially, the fit paragraph says that only
the Mg0.9Fe0.1SiO3 data were used for the least-squares Birch-Murnaghan fit,
which yielded the shared `K0=261(4) GPa` with `K0'=4`. These are therefore one
elastic fit applied with three independently measured ambient volumes, not a
joint regression of all three composition tables.

## Final article data and direct reproduction

Tables 2 and 3 contain 16 room-temperature rows: nine Fe10 fit observations,
three MgSiO3 comparison observations, and four Fe20 comparison observations.
The complete transcription, including the printed gold lattice parameter,
perovskite volume and lattice parameters, all parenthetical uncertainties, run
numbers, and explicit fit-selection flag, is in
`peritheos/data/datasets/mao-1991-tables2-3-bridgmanite-pv.csv`.

With the Table 1 Fe10 value `V0=162.79 Å3` and `K0'=4` fixed, an unweighted
least-squares fit of pressure residuals to all nine 300 K Fe10 rows gives:

| Quantity | Published | Reproduction |
|---|---:|---:|
| `K0` | `261 ± 4 GPa` | `260.800 ± 1.734 GPa` |
| pressure RMSE | not reported | `0.285 GPa` |

The source says “least squares” but does not state its weights or
dependent-variable convention. The unweighted pressure objective is recorded
as the reproducible convention because it recovers the printed coefficient
from rounded rows; volume-weighted alternatives are unnecessary assumptions.

The paper also resolves the pressure calibration. Ruby R1 fluorescence and
neon compression were cross-checks. Neon became unsuitable for accurate
calibration because it produced single-crystal or spotty diffraction, and rows
where the sample bridged the anvils and the gauges disagreed by several GPa
were rejected. Only pressures derived from the Anderson, Isaak, and Yamamoto
(1989) gold EOS were placed in the tables and used for the perovskite EOS. That
calibration is already executable as `gold_anderson_1989_bm3_1`, and all 16
rows preserve the required Au lattice parameter and temperature.

Direct evaluation of the bundled Anderson curve gives pressures systematically
`0.125-0.261 GPa` below the printed Mao values (`0.193 GPa` RMSE). This is a
rounding-level reference-state issue, not a missing implementation: the bundled
Anderson `V0=67.79 Å3` is converted from the source's rounded ambient density.
Holding Anderson's `K0=166.65 GPa` and `K0'=5.4823` fixed while fitting only the
effective reference volume gives `V0=67.8454 Å3` and reduces the pressure RMSE
to `0.017 GPa`. The main gold pressure path is therefore executable; the
overall calibration remains `partially_resolved` because the article does not
print the row-wise ruby and neon cross-check readings or the unrounded Anderson
reference state used in its reduction.

## Recovered authoritative precursor data

The [Geophysical Laboratory's 1988–1989 annual report](https://archive.org/details/annualreportofd198889carn)
(also discoverable through Carnegie Science's [official year-book index](https://carnegiescience.edu/about/history/publications/carnegie-year-books))
contains an eight-page report by seven of the final article's eight authors
under essentially the same title. Its Table 8 supplies 12 numerical
observations at 298 K: three MgSiO3,
five Mg0.9Fe0.1SiO3, and four Mg0.8Fe0.2SiO3 rows. The table prints pressure,
all three orthorhombic axes, conventional-cell volume, three normalized axes,
and V/V0. All 12 rows are transcribed in
`peritheos/data/datasets/mao-1989-carnegie-table8-bridgmanite-pv.csv`; the
companion source JSON records the scan hash, page location, search trail, and
the distinction between this report and the final paper.

The precursor describes the experiment and pressure reduction unusually well:

- the three samples were synthesized from synthetic pyroxenes by laser heating
  at 40 GPa;
- neon loaded at 200 MPa supplied the hydrostatic medium;
- 16.1 keV synchrotron powder diffraction was recorded on film at CHESS;
- gold (Jamieson et al., 1982; Ming et al., 1983) and Mao et al.'s 1986 ruby
  scale were the primary gauges; their pressures were averaged, with a stated
  0.2 GPa standard deviation;
- the Hemley et al. (1989) neon EOS supplied a third pressure check.

The report says the complete three composition sets were combined. It states
no exclusions, volume uncertainties, or regression weights. An unweighted
least-squares fit of pressure residuals to the printed V/V0 values recovers all
three reported precursor results from rounding-level data:

| Precursor fit | Report | Reproduction |
|---|---:|---:|
| Murnaghan, K0 and K0' free | `275 ± 8 GPa`, `3.7 ± 0.8` | `275.61 ± 8.79 GPa`, `3.669 ± 0.828` |
| Murnaghan, K0'=4 | `272.5 ± 2.4 GPa` | `272.26 ± 2.40 GPa` |
| Birch, K0'=4 | `273.4 ± 2.4 GPa` | `273.21 ± 2.39 GPa` |

This also resolves the model identity: the source's second-order Birch form is
the `K0'=4` special case of the repository's executable BM3 implementation.

## Why the precursor and final moduli differ

The final article makes the earlier `273.4(2.4)` versus final `261(4) GPa`
discrepancy explainable. It adds the four 300 K Fe10 runs D1, D3, D12, and D20,
fits only the nine Fe10 rows rather than combining all three compositions, and
uses Anderson et al. gold pressures instead of the precursor's average of the
older gold and ruby scales. Several row pressures change, and P10A5 changes
from `153.51` to `153.29 Å3`. The precursor remains a useful reduction
diagnostic, but it is no longer mistaken for the final fit input.

The three production records are now `bundled` / `parity`. Their `fit_datasets`
point to the final table; the earlier Carnegie Table 8 resource remains listed
only as a diagnostic precursor dataset.

The executable audit is `scripts/reproduce_mao_1991_bridgmanites.py`. It performs
the final nine-row parity refit, checks the three production curves at equal
compression, and retains all three precursor fits for comparison.

## LitCurate disposition

| Row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 29 | `litcurate_3a319d747245f4a1` | ACCEPT | Source-reported MgSiO3 fit with composition-specific V0. |
| 30 | `litcurate_d89301cb5e4082f8` | ACCEPT | Source-reported Fe10 fit with composition-specific V0. |
| 31 | `litcurate_bd7bfa49a3abafa2` | ACCEPT | Source-reported Fe20 fit with composition-specific V0. |
| 32 | `litcurate_072d0c07b96dcd17` | REJECT | Earlier MgSiO3 comparison; incomplete citation trace. |
| 33 | `litcurate_d945bf809c25ef83` | REJECT | Earlier MgSiO3 comparison; incomplete citation trace. |
| 34 | `litcurate_0d34446b774288f5` | REJECT | Earlier MgSiO3 comparison; incomplete citation trace. |
| 35 | `litcurate_ff4d5fbec277f38e` | REJECT | Earlier MgSiO3 comparison; incomplete citation trace. |
| 36 | `litcurate_9060bce3627f7052` | REJECT | Mg0.88Fe0.12SiO3 comparison belongs to its underlying cited paper. |

Result: **3 production records, 5 rejected citation traces**.
