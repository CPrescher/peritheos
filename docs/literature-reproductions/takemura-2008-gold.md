# Takemura and Dewaele (2008) Au recalculation metadata

Audit date: 2026-10-12. Scope: `gold_takemura_2008_vinet_6` recalculation
readiness, pressure coordinate, and a fresh conditional refit. Published
coefficients, uncertainties, fixed parameters, and source observations are unchanged.

The bundled `gold_takemura_2008_table3_compression` dataset contains 86 observed
rows across seven runs (13, 21, 6, 10, 9, 9, and 18 rows). Its packaged loader
verifies SHA-256 `04e2015f65f5203d89966ebf4b08eb3000d9a3fb7bfc7e3fd032e253e72b3858`.
It retains the two reported ruby pressure reductions and the measured
`am(111)` lattice parameter. Conventional fcc cell volume is `am(111)^3`;
atomic volume is `am(111)^3/4`. These are observed compression data, not
calculated EOS output or digitized curve coordinates.

Recalculation status is now `ready` for row-wise pressure-scale conversion and
diagnostic refitting. The linked executable `ruby_dorogokupets_oganov_2007`
implements `P = 1884 s (1 + 5.5 s)` GPa, with `s = lambda/lambda0 - 1`, from
[Dorogokupets and Oganov (2007), Eq. (16)](https://uspex-team.org/static/file/PressureScales-PRB-2007.pdf).
All 86 reported source pressures can be inverted to effective R1 wavelength
ratios and converted to another executable ruby scale. Focused tests check the
inversion against an independent analytic root, conversion against the published
Mao equation, and the return conversion against the reported source pressures.

This readiness does not recover raw ruby spectra, row-wise R1 wavelengths or
zero-pressure wavelength references, unrounded pressure/lattice reductions,
pointwise pressure/lattice uncertainties, or original fit weights/covariance.
The source table does not provide pointwise uncertainties. Deterministic
conversion of rounded reported pressures cannot reconstruct those missing
inputs or establish exact author-fit reproduction. The stored constrained fit
remains `V0 = 67.848` cell Å³, `K0 = 167` GPa, and `K0_prime = 5.9(1)`.

[Takemura and Dewaele (2008)](https://doi.org/10.1103/PhysRevB.78.104119)
summarize their experimental limit as 123 GPa in the abstract. Table III's
highest-pressure row reports `P = 122.5` GPa on the Zha coordinate and
`P-prime = 130.8` GPa on the Dorogokupets-Oganov coordinate, with
`am(111) = 3.6699` Å. This Dorogokupets record's experimental range is therefore
`[0, 130.8]` GPa. The paper's rounded 123 GPa limit remains valid as a description
of the Zha-coordinate experiment; changing coordinates does not add observations
or extend the experiment.

Only factual source-table transcriptions are bundled. This metadata correction
does not grant rights to the third-party article, spectra, or table images.

## Fresh refit and parity assessment

`python -m scripts.reproduce_takemura_2008_gold` writes the dedicated report
`docs/data/takemura-2008-gold-reproduction.json`; `--check` verifies it.
The Vinet equation and compression-branch volume inversion are implemented
independently of the library evaluator. All 86 observations are used with
the published `V0 = 67.848` cell Å³ and `K0 = 167` GPa fixed. No rows are
excluded and no source-pressure conversion is needed for either printed column.

| Ruby coordinate | Published K0-prime | Equal-pressure residual fit | Equal-volume residual fit |
|---|---:|---:|---:|
| Dorogokupets-Oganov P-prime | 5.9(1) | 5.92871 | 5.92928 |
| Zha P | 5.5(1) | 5.51340 | 5.53013 |

**Conditional numerical parameter parity is achieved:** all four diagnostic
fits lie inside the respective printed ±0.1 parameter widths. On the stored
Dorogokupets coordinate, the published curve's pressure RMSE is 0.799948 GPa;
the equal-pressure refit's RMSE is 0.777873 GPa. The latter reproduces the
earlier coefficient estimate, now with an explicit scope and correct statistical
qualification. This comparison uses the source error width without assigning
it a standard-deviation or confidence-level interpretation.

The previous common audit set `absolute_sigma=True` despite the absence of
pointwise errors and reported a refit error and reduced chi-square from the
implicit unit pressure scale. Those values cannot establish statistical parity.
The corrected entry leaves refit error, combined-two-sigma comparison, and
reduced chi-square unavailable. Equal-pressure and equal-volume residuals are
explicit diagnostic assumptions; exact author weights, objective, unrounded
reductions, uncertainty estimation, and covariance remain unreproduced.

The shared refit JSON, its Markdown row and qualification, and the paper-ledger
Takemura row were reconciled with targeted edits. Their generators now preserve
this conditional scope. No complete global regeneration was performed; other
chats' changes and unrelated records remain untouched. Derived catalog exports
still require any later coordinated project-wide refresh.
