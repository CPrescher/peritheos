# Boffa Ballaran et al. (2007) CaIrO3 perovskite and post-perovskite audit

## Scope and source

This audit covers every LitCurate candidate attached to Boffa Ballaran,
Trønnes, and Frost, “Equations of state of CaIrO3 perovskite and
post-perovskite phases,” *American Mineralogist* **92**, 1760–1763 (2007),
[doi:10.2138/am.2007.2715](https://doi.org/10.2138/am.2007.2715).

Both candidates were already present in the Peritheos catalog before this
batch.  Their source audit used an author-uploaded copy of the final published
pages.  Table 1 supplies eleven pressure–lattice–volume observations for each
polymorph, including one decompression observation that the text explicitly
excludes from fitting.  The bundled CSVs preserve all eleven rows and mark the
ten fitted compression rows.  Table 2 and the abstract report the BM3
coefficients and standard deviations.  The experiment used a 4:1
methanol–ethanol pressure medium and the Mao et al. (1986) ruby calibration.

## Numerical verification and low pressure derivatives

The source explicitly defines the standard Eulerian finite strain and
normalized stress used by third-order Birch–Murnaghan EOS fits.  It treats
`V0`, `K0`, and `K0_prime` as fitted parameters; neither low derivative is a
transcription of a fixed coefficient.  The authors also publish normalized
stress slopes and confidence ellipses and caution against extrapolating the
short, 0–7.79 GPa fits beyond roughly 15–20 GPa.

| Phase | Published `(V0 A3, K0 GPa, K0')` | Independent errors-in-variables refit | Published-curve pressure RMSE |
|---|---|---|---:|
| Pbnm perovskite | `(229.463, 198, 1.2)` | `(229.463435, 198.370553, 1.233101)` | 0.027725 GPa |
| Cmcm post-perovskite | `(226.38, 181, 2.3)` | `(226.377363, 181.466972, 2.276317)` | 0.023982 GPa |

Every independently fitted coefficient lies within the corresponding printed
standard deviation.  The BM3-implied normalized-stress slopes, −831.6 GPa and
−461.55 GPa, also agree with the source values of −825(244) GPa and −433(214)
GPa.  The low `K0_prime` values are therefore retained as short-range fitted
results, with the source's extrapolation warning, rather than normalized to a
more typical value.

## Disposition of every same-DOI candidate

| LitCurate identifier | Candidate | Disposition |
|---|---|---|
| `litcurate_88696b48a0b7225b` | CaIrO3 Pbnm/Pnma perovskite BM3 | **Already accepted** as `cairo3_perovskite_boffa_ballaran_2007_bm3_1`; phase, conventional-cell basis, complete Table 1 data, free parameters, and numerical parity verified. |
| `litcurate_e17fb8989f90c501` | CaIrO3 Cmcm post-perovskite BM3 | **Already accepted** as `cairo3_post_perovskite_boffa_ballaran_2007_bm3_1`; phase, conventional-cell basis, complete Table 1 data, free parameters, and numerical parity verified. |

Final catalog count for DOI `10.2138/am.2007.2715`: **2 production records**.
Net new production count in this batch: **0 records**.  No same-DOI candidate
is rejected or held.

## Reproduction assets

- `peritheos/data/materials/cairo3_perovskite.eosmat`
- `peritheos/data/materials/cairo3_post_perovskite.eosmat`
- `peritheos/data/datasets/cairo3-perovskite-boffa-ballaran-2007-table1-compression.csv`
- `peritheos/data/datasets/cairo3-post-perovskite-boffa-ballaran-2007-table1-compression.csv`
- `tests/test_cairo3_boffa_ballaran_2007.py`

Run `python -m pytest -q tests/test_cairo3_boffa_ballaran_2007.py` to repeat
the source-table checksum, BM3 curve, normalized-stress, independent-refit,
phase-identity, and pressure-calibration checks.
