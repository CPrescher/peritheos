# Diopside source audit and equation reproduction

Diopside is added as nominal CaMgSi2O6, C2/c clinopyroxene, with conventional
Z = 4 cell volumes. The default is the independent room-temperature BM3 of
Thompson and Downs (2008). Zhao et al. (1998) supplies the published thermal
BM3 covering the reported 300–1280 K observations. Published coefficients
and original observations are preserved.

**Published-model and parameter parity are established.** The thermal curves
agree at graphical resolution, and all five varied coefficients of the
independent refit lie within their reported errors when the Figure 2 pressure
reading is used for the inconsistent Table 1 row 15. The static record also
passes equation, observation and parameter-agreement checks. Source-table
errors are recorded as provenance and do not block acceptance of the published
EOS. Exact recovery of the authors' original numerical inputs and regression
workflow is a separate, unrecovered result.

## Primary authority and source selection

- Zhao, Y., Von Dreele, R. B., Zhang, J. Z., and Weidner, D. J. (1998).
  *Thermoelastic Equation of State of Monoclinic Pyroxene: CaMgSi2O6 Diopside*.
  The Review of High Pressure Science and Technology 7, 25–27.
  [Publication and primary full text](https://doi.org/10.4131/jshpreview.7.25).
  All three pages were inspected, including Equation (1), Table 1 and Figure 2
  on page 26. The J-STAGE article has no separate numerical supplement.
- Thompson, R. M., and Downs, R. T. (2008). *The crystal structure of diopside
  at pressure to 10 GPa*. American Mineralogist 93, 177–186.
  [Publication](https://doi.org/10.2138/am.2008.2684),
  [primary article archived by RRUFF](https://rruff.net/doclib/am/vol93/AM93_177.pdf).
  Tables 1–3 and the experimental/results text on pages 179–180 define the
  adopted static EOS and observed ambient structure. The article identifies
  structural CIF deposit AM-08-006; its relevant ambient structure is fully
  printed in Table 2. The deposit is not a numerical EOS regression supplement.

The review also inspected Cameron et al. (1973), American Mineralogist 58,
594–618, and both Raudsepp, Hawthorne and Turnock (1990) studies, American
Mineralogist 75, 1274–1281 and Canadian Mineralogist 28, 93–109. The synthetic
D0/Ni0 model is unsuitable as a literal fully occupied endmember structure:
its refined M2 population is Ca0.89Mg0.11. It was not silently idealized.
Thompson's observed nominal endmember refinement provides a simpler,
independent structural reference. No acoustic or elastic-tensor modulus is
used as an isothermal EOS coefficient.

The source receipt in
`peritheos/data/datasets/diopside-source-audit.json`
records retrieval URLs, byte counts and SHA-256 hashes of the inspected PDFs.
The publisher files are not redistributed with the numerical data. No authors
were contacted.

## Material and structural identity

The material identifier is `diopside`; the formula is CaMgSi2O6. The conventional
cell is unique-b C2/c, No. 15, with four formula units and 40 atoms. The stored
observed room-temperature cell is Thompson Table 1 run 0:

| a in Å | b in Å | c in Å | β in degrees | Printed volume in Å³ |
|---|---|---|---|---|
| 9.7397(3) | 8.9174(3) | 5.2503(2) | 105.866(3) | 438.63(3) |

The metric gives 438.63123 Å³. The atom positions are from the **observed**
1-atm column of Table 2, not the adjacent ideal geometric model. Mg1 and Ca2
occupy 4e (0,y,¼); Si and three independent oxygens occupy 8f. Full nominal
populations give Ca4Mg4Si8O24. Position errors and isotropic B values are retained
in `diopside-thompson-2008-table2-ambient-structure.csv`.

The natural specimen is clear, colorless DeKalb diopside, Smithsonian R18685.
The source uses a nominal full-site CaMgSi2O6 model without a new analytical
purity determination or refined occupancy errors. This provenance does not
claim that natural and synthetic specimens have identical trace composition.
The structure is a diffraction reference, not a measured high-temperature
atomic-coordinate or axial-strain model.

The structural volume 438.63(3), fitted static V0 = 438.66(2), and Zhao's
printed ambient volume 438.67(6) Å³ are distinct source quantities. They have
not been forced to agree by rescaling the cell or combining coefficients.

## Static EOS and source reproduction

`diopside_thompson_2008_bm3` is the default published static record. Page 179
reports V0 = 438.66(2) Å³, KT0 = 118(1) GPa and K′ = 3.8(3), fitted to the
weighted Table 1 data using Levenberg–Marquardt. Table 3 instead prints an
error of 0.2 on K′; the adopted record retains the prose error of 0.3 and
records the discrepancy. No covariance matrix or exact weighting objective
is printed. The source's separate fixed-K′ = 4 alternative has KT0 = 117.1(2)
GPa and V0 = 438.66(2) Å³; it is not combined with the free-K′ coefficients.

The source is room-temperature compression, without an exact numerical
experimental temperature. The operational reference is nominally 300 K and
has no thermal dependence. This chosen label is explicit; no measured
300 K temperatures have been added to Table 1. The static validity check
accepts that reference temperature only.

All 14 observations, their run identifiers, intensity-collection flags,
pressures, lattice parameters, β angles, cell volumes and printed errors are
retained. Their pressure range is 0.0001–10.16 GPa. The pressure medium was
15:4:1 methanol:ethanol:water, with Mao et al. (1978) ruby calibration
([DOI](https://doi.org/10.1063/1.325277), existing `ruby_mao_1978` calibration).
The source estimates 0.04 GPa pressure error. Ruby wavelengths are not printed,
so observation-level pressure recalculation is unavailable. The same crystal
was reloaded after gasket failure near 5 GPa.

Independent finite-strain evaluation gives 10.12944 GPa for the last Table 1
volume, 407.16 Å³, compared with 10.16(4) GPa. Across all 14 rows the published
curve has pressure RMS residual 0.037397 GPa and maximum absolute residual
0.068227 GPa. Cell metrics agree with the printed volumes to within 0.008 Å³.

A diagnostic fit with effective variance σP² + (KσV/V)² gives approximately
V0 = 438.6554 Å³, KT0 = 117.8642 GPa and K′ = 3.7963. Each is within its
published error. Ambient pressure is treated as an imposed reference in this
diagnostic; its missing error is not reclassified as a measured zero error.
This demonstrates parameter agreement under a stated independent objective.
It does not recover the unreported original weighting implementation,
covariance scaling, regression code or software version.

## Thermal equation and exact model mapping

`diopside_zhao_1998_thermal_bm3` preserves the published thermal parameterization.
Equation (1) is

\[
f=\tfrac12[(V_T/V)^{2/3}-1],\qquad
P=3K_T f(1+2f)^{5/2}[1+\tfrac32(K'-4)f].
\]

The source defines K_T = K_T0 + (∂K/∂T)(T−300), α(T) = a+bT and
V_T = V0 exp(∫α dT). Anchoring the integral at 300 K gives

\[
V_T=438.67\exp[a(T-300)+\tfrac12b(T^2-300^2)]\;\mathrm{Å^3}.
\]

This is exactly `ThermalReferenceStateEOS` with BM3,
`thermal_expansion_law="linear_temperature"` and
`reference_volume_law="integrated_expansivity"`. The intercept a is defined
at T = 0, not at 300 K. In particular α(300) = 2.884×10⁻⁵ K⁻¹. A direct
linear-volume law, a quadratic truncated volume, and an expansivity intercept
at 300 K would represent different equations.

| Parameter | Adopted published value | Authority |
|---|---|---|
| V0 | 438.67 Å³ ambient anchor | Table 1 ambient rows and Equation (1)'s ambient-volume definition |
| KT0 | 109.1 GPa; reported error 4 GPa | Abstract and Figure 2; prose 109(4) |
| K′ | 4.84; reported error 0.6 | Abstract and Figure 2; prose 4.8(6) |
| ∂K/∂T | −0.0205(4) GPa K⁻¹ | Page 26: −2.05(4)×10⁻² |
| a | 2.32(5)×10⁻⁵ K⁻¹ | Page 26 |
| b | 1.88(7)×10⁻⁸ K⁻² | Page 26 |

V0 is the source's printed ambient cell-volume anchor. The paper does not
separately print a fitted V0 or state whether it was varied. Its observation
error 0.06 Å³ is retained in the table, and is **not** presented as a fitted
parameter error. The stored V0 parameter error is missing. Fixed V0 in the
diagnostic refits is an explicitly stated reconstruction choice, not a
claim about the original workflow. Error confidence and covariance are
unreported. The source neglects the c/T² expansivity term and higher temperature
derivatives; K′ is constant. Page 27's separate pressure-correction example
uses K′ = 4.85 and ∂K/∂T = −0.0206; these rounded alternatives are not mixed
into the Figure 2 parameterization.

The 300–1300 K mean of α(T) is 3.824×10⁻⁵ K⁻¹, agreeing with the independent
page 27 result 3.82(8)×10⁻⁵ K⁻¹. Off-reference Table 1 checks give:

| Observed pressure in GPa | Temperature in K | Volume in Å³ | Model pressure in GPa |
|---|---|---|---|
| 8.18 | 1272 | 421.43 | 8.27721 |
| 6.88 | 1280 | 426.65 | 6.78220 |
| 2.17 | 1272 | 444.61 | 2.21465 |

The first two differences are smaller than the source's quoted aggregate
0.08–0.12 GPa pressure precision scale. This is an equation/curve consistency
check, not proof that every table row generated the reported fit. The actual
experimental envelope is 0–8.18 GPa and 300–1280 K; the 1300 K mean-expansion
comparison is not an extension of that envelope or a phase-stability guarantee.

The thermal model defines mechanical P–V–T behavior. It does not supply a
caloric free-energy model. Neither source provides an axial thermal strain
law for diffraction.

## Thermal observations and recorded source inconsistencies

All 32 printed Table 1 rows are preserved in order, including the duplicate
ambient rows 1 and 32. No pressure, volume, lattice dimension, error or
calibrant observation is synthesized. The source reports layer-packed NaCl
on Decker (1971), DOI 10.1063/1.1660714. NaCl lattice readings are absent, so
pressure recalculation is blocked. The 0.08–0.12 GPa aggregate pressure
precision is not assigned as individual row errors. Temperature errors are
not printed. The stated thermocouple is W–26%Rh/W, retained literally despite
the unusual composition. The source describes heating after cold compression
and fitting only the subsequent hydrostatic paths; exact selected row indices,
regression weights, treatment of repeated ambient data and original fit code
are unavailable.

The following discrepancies remain visible rather than corrected:

- Source row 15 is printed at 6.95 GPa, 1072 K, V = 431.66(6) Å³. The
  published EOS predicts 4.56302 GPa there, a −2.38698 GPa residual. Neither
  a substitute pressure nor an author exclusion is asserted.
- Both ambient rows print a = 9.757, b = 8.939, c = 5.276 Å and β = 105.43°.
  These imply 443.57529 Å³, whereas their volume column prints 438.67(6) Å³.
  The stored structure therefore comes from the independent Thompson source.
- Rows 6, 11 and 16 have metric-minus-printed-volume discrepancies of about
  0.352, 0.211 and 0.665 Å³. Row 4 prints a = 9.584(7), which is retained
  with error 0.007 Å. None of these entries is normalized to an inferred value.
- Figure 1's plot heading says 6.88 GPa/1280 K; the page 26 discussion of that
  pattern instead says 6.95 GPa/1280 K. This cannot resolve row 15 at 1072 K.
- Figure 4 labels the expansion-pressure slope −1.54(2)×10⁻⁶ K⁻¹ GPa⁻¹,
  whereas the prose prints −1.52(2)×10⁻⁶. This secondary regression is not a
  coefficient of the accepted mechanical EOS and is not fitted here.

The published curve's RMS residual is 0.426861 GPa across all 32 rows and
0.065521 GPa in a diagnostic calculation omitting only row 15. The complete
unweighted pressure refit, with the ambient volume fixed, yields roughly
KT0 = 118.93 GPa, K′ = 0.396, a = 2.3783×10⁻⁵, b = 1.2913×10⁻⁸ and
∂K/∂T = −0.00523 GPa K⁻¹. It fails published coefficient agreement.
The sensitivity fit excluding row 15 yields KT0 = 107.48 GPa, K′ = 5.235,
a = 2.2928×10⁻⁵, b = 1.9187×10⁻⁸ and ∂K/∂T = −0.020478 GPa K⁻¹,
within the individual reported errors. Both fits, their selections and
objectives are retained in the machine-readable report. The sensitivity result
suggests a source-table problem, but does not establish a correction or
recover the authors' original workflow.

The original numerical input, exact row-15 pressure, V0 treatment, row
selection, weights, covariance and regression code would be needed to recover
the authors' exact fitting workflow. Their absence does not invalidate the
published-model and parameter parity demonstrated below. Paired NaCl readings
would additionally enable pressure-scale reproduction. Contact status remains
**not contacted**.

## Plot audit of the thermal discrepancy

The high-resolution Figure 2 check independently supports a **printed pressure
error** in Table 1 row 15. Axis calibration from the printed tick marks places
the matching hollow-circle marker at approximately **4.58 GPa and 431.63 Å³**.
Conservative graphical reading bounds are ±0.05 GPa and ±0.2 Å³; these are
reading bounds, not experimental standard errors. This agrees with the row's
printed volume 431.66 Å³ and the published EOS pressure **4.56302 GPa** at
1072 K. The table position (6.95 GPa, 431.66 Å³) has no corresponding marker.
The plot therefore changes the assessment from a suspected table problem to
one with direct graphical corroboration. It does not recover an exact original
pressure to the precision of a numerical source table.

The plotted ambient marker is about 438.8 Å³ within the same reading bound,
corroborating the printed 438.67 Å³ volume rather than the 443.575 Å³ implied
by the inconsistent ambient lattice dimensions. It cannot establish corrected
a, b, c or β values, which remain untouched.

All six Figure 2 isotherms have the expected ordering. A point-free cross
section near 11.292 GPa gives graphical-minus-model differences from −0.249
to −0.060 Å³. This is agreement at the figure's graphical scale, rather than
proof of identity with the original plotting program. The plotted temperature
labels are 300, 473, 673, 873, 1073 and 1273 K; the caption's wording about
200 K steps above 300 K is inexact.

Figure 3 displays raw V–T points without identifying their pressures, so it
cannot independently recover the questionable pressure. Figure 4's
zero-pressure mean expansivity is consistent with the adopted thermal
coefficients; its slope-label/prose discrepancy remains as described above.

As a separate sensitivity calculation, replacing only the analysis pressure
of row 15 by the graphical reading rounded to 4.58 GPa gives a published-curve pressure
RMS residual of **0.06456 GPa** over all 32 rows. An independent unweighted
fit then gives K0 ≈ 107.73 GPa, K′ ≈ 5.122, a ≈ 2.2973×10⁻⁵ K⁻¹,
b ≈ 1.9025×10⁻⁸ K⁻² and dK/dT ≈ −0.020206 GPa K⁻¹, each within the
individual reported parameter error. This establishes parameter parity within
the reported errors with the documented table-error treatment. The original observations, published
coefficients and room-temperature default remain unchanged. The catalog
records this positive parity assessment separately from the unrecovered exact
original workflow.

The analysis pixel coordinates, affine calibration, reading bounds and
numerical checks are retained in
[the plot-check report](../data/diopside-plot-check.json), generated by
`python -m scripts.diopside_check_plots`. Its PDF hash identifies the inspected
source version. The annotated comparison PNG under `output/` is an attributed
inspection artifact; the source crop is not included in the numerical CC0
dedication or redistributed as a bundled dataset.

## Corrected dataset for continued analysis

Use `peritheos/data/datasets/diopside-zhao-1998-table1-corrected-pvt.csv`
for subsequent P–V–T work. The material document registers it as
`diopside_zhao_1998_table1_corrected_pvt` and identifies it as the recommended
refitting dataset. It retains all 32 source rows and applies only the
plot-supported row-15 pressure correction: **6.95 → 4.58 GPa**. The exported
value is rounded to 0.01 GPa, and the reproduction calculation uses precisely
that value rather than the additional digits from axis calibration.

The original pressure appears in `source_pressure_gpa`; correction flags,
pressure origin and a row-specific note identify the change. A separate
`pressure_graphical_reading_bound_gpa` records the ±0.05 GPa graphical bound.
The experimental pressure-error field remains blank, because that uncertainty
was not reported. Temperatures, volumes and their source errors are unchanged.
Row 32 is marked as a repeated ambient entry of row 1, without automatically
excluding it.

This is a corrected P–V–T view. Unresolved lattice dimensions remain in the
original transcription; they are not copied as supposedly corrected values
into this view. Quality notes mark affected source rows. The portable guide
`peritheos/data/datasets/diopside-zhao-1998-corrected-pvt.README.md` describes
usage, limits and provenance. The original source CSV and numerical-notation
transcription remain unchanged, and no published EOS coefficient is replaced.

Regenerate or check the corrected export with
`python -m scripts.diopside_correct_dataset` or its `--check` option.

## Executable checks and integration

Run `python -m scripts.diopside_reproduce --check` for source-notation/CSV
agreement, resource hashes and published-curve results. The report is
[`diopside-reproduction.json`](../data/diopside-reproduction.json).
`tests/test_diopside.py` checks the catalog, schema, observation transcription,
source checkpoints, exact thermal conventions, cell contents, inversion,
bulk-modulus derivatives, native execution and metadata retention.
`crates/peritheos/tests/diopside.rs` checks Rust loading, evaluation, inversion
and document round trips. Validity checks use interior forward states because
an exactly bounded inverse can produce a pressure a few floating-point units
outside the closed range; the boundary round trips are still tested numerically.

Shared ledgers and aggregate counts must be reconciled with the other material
additions; see the [integration note](diopside-integration.md). Source licensing
is unspecified. CC0 covers contributor transcription, normalization and metadata
rights only, as stated in
`peritheos/data/datasets/diopside.LICENSE.md`.
