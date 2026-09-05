# Mao et al. (2011) ferropericlase spin-branch audit

## Sources and correction

- Primary paper: Zhu Mao, Jung-Fu Lin, Jin Liu, and Vitali B. Prakapenka,
  “Thermal equation of state of lower-mantle ferropericlase across the spin
  crossover,” *Geophysical Research Letters* 38, L23308 (2011),
  [doi:10.1029/2011GL049915](https://doi.org/10.1029/2011GL049915).
- Corrected author-hosted final PDF:
  [MaoFp-EoSGRL2011.pdf](https://www.jsg.utexas.edu/lin/files/MaoFp-EoSGRL2011.pdf),
  SHA-256 `54d40b7391f64761543038396473fd062f43e57172034d6a3b0d5b12e2627981`
  (audited 2026-09-05).
- Published correction: Mao et al. (2012), *Geophysical Research Letters* 39,
  L02399, [doi:10.1029/2011GL050814](https://doi.org/10.1029/2011GL050814).
- Pressure scale: Fei et al. (2007), “Toward an internally consistent pressure
  scale,” [doi:10.1073/pnas.0609013104](https://doi.org/10.1073/pnas.0609013104).

The corrected PDF explicitly says that errors throughout the article have
been corrected there and that the originally published text remains in the
HTML. This audit therefore takes the corrected PDF as authority. In
particular, Results paragraph 7 says that the fit used the **300 K** P-V data,
and reports the low-spin extrapolated volume as `74.4(6) A^3`, not the
`74.44 A^3` value still exposed by some uncorrected metadata.

## Material, experiment, and pressure scale

The sample is polycrystalline `(Mg0.75Fe0.25)O`, mechanically ground with
5 wt.% Au for approximately 10 hours and pressed into 10-15 micrometre disks.
The B1 rocksalt structure is cubic `Fm-3m` with four formula units per
conventional cell. The mixture was loaded in rhenium-gasketed diamond-anvil
cells between dried KCl foils below 60 GPa and dried NaCl above 60 GPa.

Diffraction was measured at GSECARS, Advanced Photon Source. The experiments
used double-sided laser heating and collected patterns at 300, 1200, 1500,
1800, and 2000 K. Pressures came from in-situ Au using the thermal scale of
Fei et al. (2007). The published figure and accessible article do not contain
the row-wise Au lattice parameters, so the reduced pressure coordinate cannot
be independently recalculated.

The exact composition already exists in `mg075fe025o.eosmat` for the
independent Matsui et al. (2012) high-spin thermal EOS. The Mao fits are not
duplicates: they come from a different DAC experiment and add explicitly
separated 300 K high-spin and low-spin reference branches. The top-level
material phase label is now spin-neutral; spin identity and validity remain
record-specific.

## Corrected EOS parameterizations

Results paragraph 7 calls both fits third-order isothermal Birch-Murnaghan
EOS with fixed `K0'=4`. Peritheos preserves that source family label as BM3;
the fixed derivative makes each curve algebraically identical to the BM2
special case.

| Branch | V0 (A^3/conventional cell) | K0 (GPa) | K0' | Executable validity |
| --- | ---: | ---: | ---: | --- |
| High spin | 76.34(1), fixed | 162(1) | 4 fixed | 300 K, 0-50 GPa |
| Low spin | 74.4(6), fitted extrapolation | 166(7) | 4 fixed | 300 K, 76.1-129.6 GPa |

The source places the 300 K spin crossover between 50 and 75 GPa. Neither
reference curve is a model of that continuous mixed-spin interval. The
low-spin `V0` is an extrapolation of the high-pressure branch and must not be
presented as an observed ambient low-spin cell.

## Primary data and digitization

The publisher page lists official Text S1 and Table S1 attachments, including
`grl28755-sup-0006-ts01.txt`, but each tested official Wiley download route
returned HTTP 403 during this audit. The main figure contains clear vector-like
open-circle markers, and digitization is permitted by the project protocol.

Corrected Figure 1 was rendered at 500 dpi. Linear pixel-to-axis calibration
used every major tick from 0 to 140 GPa and 50 to 80 A^3. All 42 black 300 K
marker centers were extracted and rounded to 0.1 GPa and 0.01 A^3. The
resulting dataset preserves:

- 15 high-spin observations from 0.0 to 48.4 GPa;
- nine mixed-spin observations from 50.8 to 70.8 GPa, excluded from both
  branch fits; and
- 18 low-spin observations from 76.1 to 129.6 GPa.

The dataset's `0.10 GPa` and `0.03 A^3` uncertainty columns describe one-pixel
digitization resolution, not experimental standard deviations. No row-wise
experimental uncertainty is invented.

## Numerical reproduction

Run:

```bash
.venv/bin/python scripts/reproduce_mao_2011_ferropericlase.py
```

An unweighted pressure-space refit of the 15 high-spin points, with
`V0=76.34 A^3` and `K0'=4` fixed as in the source, gives
`K0=163.0461 GPa`. This is 1.05 GPa above the reported `162(1) GPa`, comparable
to the printed uncertainty and finite Figure 1 marker/curve width. Restricting
the diagnostic to points at or below 45 GPa, farther from the crossover onset,
gives `K0=162.70 GPa`. The rounded published curve has a 0.2492 GPa pressure
RMSE against the selected digitized points.

The corresponding 18-point low-spin refit gives `V0=74.4715 A^3` and
`K0=165.6619 GPa`, both well inside the printed `0.6 A^3` and `7 GPa`
uncertainties. The rounded published curve has a 0.4943 GPa pressure RMSE.

These results establish equation convention, conventional-cell volume basis,
branch assignment, and coefficient parity. They do not reconstruct the
authors' unreported objective, weights, or covariance.

## LitCurate disposition: both same-DOI candidates

| LitCurate identifier | Candidate | Disposition | Mapping |
| --- | --- | --- | --- |
| `litcurate_0017737ebf338199` | `(Mg0.75Fe0.25)O` high-spin reference | **Accepted** | `mg075fe025o_mao_2011_high_spin_bm3_1`; corrected `V0=76.34 A^3` fixed, `K0=162(1) GPa`, fixed `K0'=4`; 15 digitized 300 K points. |
| `litcurate_33fe12afaa4758a7` | `(Mg0.75Fe0.25)O` low-spin reference | **Accepted with corrected coefficient** | `mg075fe025o_mao_2011_low_spin_bm3_1`; corrected `V0=74.4(6) A^3`, `K0=166(7) GPa`, fixed `K0'=4`; 18 digitized 300 K points. |

Zotero's local API and connector were running during the audit. An exact DOI
search for `10.1029/2011GL049915` returned no existing item; this scoped audit
did not import a new record.
