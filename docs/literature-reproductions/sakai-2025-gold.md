# Sakai et al. (2025): Au source observations and conditional calibration

Audit date: 2026-10-12. The published Au Rydberg–Stacey coefficients remain
`V0=67.716 Å³/cell`, `K0=167 GPa`, `K0′=5.780(4)`, and fixed `K∞′=5/3`.
**Conditional numerical parity is achieved for Au's independently refitted
K0′ central value at the paper's published precision.** V0, K0 and K∞′ are
fixed at the paper's values and are not independently recovered. Uncertainty
parity and full sequential-calibration parity remain unverified: the original
weights, covariance and exact preceding regressions are not established.

## Source coverage

The actual [Supplementary Data 1 workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43246-025-00792-5/MediaObjects/43246_2025_792_MOESM5_ESM.xlsx)
was inspected. Its SHA-256 is
`c559b29fea9c793831e4bc6deef4d62712b5e78821fc420969f709cbd1060ff9`.
Tables S1–S7 contain observed conventional cell volumes, not only ratios.
Table S8 is a 35-row **reference and calibrant index**, not a numerical
compilation of the preceding studies. Table S9 contains calculated final
curve nodes and remains separate and unchanged.

| Source | Rows retained | Role for Au |
|---|---:|---|
| S1, Cu–Re, runs 1a/1b | 139 | Cu → Re ancestor |
| S2, Pt–Re, run 2 | 33 | Re → Pt ancestor |
| S3, Cu–W, run 3 | 54 | Cu → W ancestor |
| S4, Au–W, run 4 | 19 | Direct W → Au input |
| S7, Pt–MgO, runs 7c/7b/7a in source order | 81 | Pt → MgO ancestor |
| Dewaele et al. (2018), supplementary Tables 4–5 | 84 | Direct Re → Au input |
| Dewaele et al. (2004), Table I, existing resource | 37 | Ruby → Au input |
| Takemura and Dewaele (2008), Table III, existing resource | 86 | Ruby → Au input |
| Hirose et al. (2008), Table 1, existing resource | 21 | MgO → Au input; nine rows at 300 K |

The five direct Au source families retain **247 rows**. There are **235 RT
rows**, of which **210** enter the documented diagnostic after excluding 25
Dewaele2018 states below 200 GPa. The 326-row network file includes the 19
S4 rows exposed in the separate Au–W file: those resources overlap. New
Sakai-specific dataset aliases reuse the complete existing resources without
changing their other record links. None of these links creates new measurements.

The [source manifest](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/sakai-2025-au-source.json)
records source/product hashes, all upstream numerical-resource pointers,
normalizations, rights, and scope. The
[cell archive](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/sakai-2025-au-workbook-cells.json)
preserves every cell on S1–S8 and S14, including literal OOXML tokens,
coordinates, strings, formulas and cached values. S5 (Mo–W), S6 (Fe–W), and
downstream S8 branches are archival context, not Au fit inputs. The normalized
network includes all numeric rows on the five Au ancestor/direct sheets.

No errors are imputed. All 19 Au errors in S4 are blank; the corresponding W
errors are retained. The source says `err.` without a confidence or standard
deviation convention. Re lattice a/c, their errors and c/a are retained where
provided. The S7 stress estimates and calculated c/a values remain identified
as calculations rather than extra independent observations. Units and paired
phases are explicit: volumes in Å³/cell, lattice parameters in Å, stresses and
pressures in GPa; fcc Au/Cu/Pt, hcp Re, bcc W, and B1 MgO.

## Earlier Au–Re measurements

[Dewaele et al. (2018)](https://doi.org/10.1038/s41467-018-05294-2)
provides numerical Au–Re lattice parameters in supplementary Tables 4–5,
pages 9–10. Both pages were visually checked against extracted values.
The [normalized table](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/gold-dewaele-2018-tables-4-5-au-re.csv)
retains all 36/31/17 rows of runs 1/2/3, including duplicate membrane settings,
Au111 and Au200, Re a/c, original Re-scale pressures, paired Al lattice
values in run 3, and every printed numeric token. No row-wise errors are
printed. The source constrained Re c/a during lattice reduction; these are
reduced diffraction observations, not independent measurements of both Re axes.

Au volumes are `a111³` and Re volumes are `sqrt(3)*a²*c/2` for conventional
four-atom fcc and two-atom hcp cells. The source's displayed atomic-volume
formulas have a factor of four in the denominator; the present archive uses
conventional cells to match Sakai Table 2. Au200 is retained but is not used
to infer a hydrostatic correction. Au111 is the stress-insensitive volume
coordinate specified in the experimental discussions.

The [Sakai supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43246-025-00792-5/MediaObjects/43246_2025_792_MOESM2_ESM.pdf)
Fig. S6 caption explicitly excludes Dewaele2018 data below 200 GPa from the
Au fit. These rows remain in the source dataset. The replay marks selection
using recalculated Sakai Re pressure; all 25 excluded rows are also below
200 GPa on the original reported Re coordinate, so that choice introduces
no threshold ambiguity in this recovered table.

## Numerical reconstruction and its limits

Run `python -m scripts.reproduce_sakai_2025_gold --output
docs/data/sakai-2025-gold-reconstruction.json` from the checkout. The
[report](../data/sakai-2025-gold-reconstruction.json) archives every reconstructed
RT pressure and selection flag, including the excluded rows. It is distinct
from the unchanged Table S9 calculated grid.

1. Reduce S4 using published W, Hirose's nine 300 K rows using published MgO,
   and Dewaele2018 using published Re. These calibrants are Sakai's rounded
   Table 2 coefficients, so this step is conditional on their earlier fits.
2. Infer effective ruby shifts from Dewaele2004's rounded Mao1986 pressures
   and Takemura2008's rounded Dorogokupets–Oganov pressures. Apply
   [Shen et al. (2020), Equation 3](https://www.osti.gov/servlets/purl/1861244),
   `P=1870*s*(1+5.63*s)` GPa. The inverse coordinates are the source-reported
   Mao equation `P=1904/7.665*((1+s)^7.665-1)` and Dorogokupets–Oganov
   Equation 16 `P=1884*s*(1+5.5*s)`. Preserve both original ruby columns.
   This recovers pressure coordinates from rounded reductions, not raw
   fluorescence spectra, wavelength references, or experimental errors.
3. Keep V0, K0, and K∞′ fixed and fit K0′ to the resulting 210 states.

With **explicitly assumed equal pressure residuals**, K0′ is **5.780096584**
and pressure RMS is **3.141600 GPa**. The successful independent K0′ refit
rounds to the paper's **5.780** at three decimal places, achieving conditional
numerical parity. The central-value difference is **0.000096584**. This
classification is conditional on the published upstream calibrants, fixed
Au parameters and assumed pressure residuals; it does not require recovery
of original weights or covariance. Equal volume residuals give **5.787771542**
with volume RMS **0.184403 Å³**. This sensitivity prevents identifying the
authors' objective solely from coefficient agreement. The reported **±0.004**
has not been reproduced, so uncertainty parity is explicitly false. Neither
result replaces the production coefficients or supplies author covariance.
No uncertainty is computed by treating unspecified source errors as
statistical sigma. In machine-readable metadata, `author_fit_parity=false`
refers to the exact original procedure and uncertainty; the separate
`conditional_numerical_parity_achieved=true` records central-coefficient parity.

The independent Equation 1 pressure implementation agrees with the native
model on all recovered W/MgO/Re reduction states. The published Au curve
against S4 alone has RMS 5.721623 GPa and maximum 17.644817 GPa; the larger
residual at some observed states is not a Table S9 implementation error.
Near-ambient W states yield slightly negative calculated pressures and are
retained without clipping.

A second diagnostic actually carries out the measured-only Cu → W → Au
sequence: fit W using all 54 S3 pairs, then fit Au using all 19 S4 pairs.
Equal pressure residuals give W K0′=4.417912338 and Au K0′=5.715517182.
Equal volume residuals give 4.423582056 and 5.709307773. This verifies that
the paired-volume data support sequential calibration, while demonstrating
the effect of omitting Table S8's earlier-study constraints.

S14 E/F also permits a conditional refit of the Cu scale conversion, using
35 **calculated external-scale input** rows at ≤1200 GPa. Equal pressure
residuals give K0′=5.202434878, K∞′=2.003716330, RMS 0.054102 GPa and maximum
0.138359 GPa. The 1209.215 GPa node remains archived but lies outside that
selection. This calculation is a diagnostic of the Cu conversion, not an
independent Cu measurement fit or proof of the original conversion objective.

The full calibration still requires the exact original reductions and
objectives of the W, Re, Pt and MgO stages, including Table S8's prior
constraints and the W → Re auxiliary input. Existing primary resources are
linked in the manifest. The Sakai2018 resource is figure-digitized rather
than recovered original Pt–Re pairs. Original fit weights, cross-material
error correlations, coefficient covariance, unrounded ruby reductions and
confidence conventions remain unrecovered. Absence of one flat P–V table
does not make this source unreconstructable: **all five Au input families
now support conditional numerical parity for the Au K0′ central value**.
Exact full-network calibration, original fitting-procedure parity and
uncertainty parity remain unestablished.

## Rights and integration scope

CC0 applies only to Peritheos contributors' transcription, normalization code
and metadata. It grants no rights in third-party source material. Sakai's
article carries CC BY-NC-ND 4.0 and Dewaele2018 carries CC BY 4.0; original
rights, attribution and exclusions remain applicable. Existing earlier-source
rights declarations are preserved. See the
[rights note](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/sakai-2025-au.LICENSE.md).

Only Sakai Au data coverage/source links and eight new dataset entries are
patched in `gold.eosmat`. Other records, existing datasets and coefficients
are preserved. Shared audit indices, fit-reproduction ledgers and generated
catalog manifests are not regenerated; their Au coverage summaries need
reconciliation after the concurrent Au follow-ups are integrated.
