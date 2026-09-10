# Iron: Brown (2000), Dubrovinsky (2000), Yamazaki (2012), Sakai (2014), and Dewaele (2006)

The five primary papers and available official supplements were read from the
user's Zotero library on 2026-09-10. Attachment and item keys, filenames and
SHA-256 hashes are recorded in
`peritheos/data/datasets/iron-source-papers-zotero.json`. Only factual numeric
transcriptions are bundled; publisher PDFs are not redistributed.

## Catalog additions

| Primary paper | New records | Scope |
|---|---:|---|
| [Brown, Fritz and Hixson (2000)](https://doi.org/10.1063/1.1319320) | 1 | Published linear Us–up Hugoniot, with an explicit solid-branch operational restriction |
| [Dubrovinsky et al. (2000)](https://doi.org/10.1103/PhysRevLett.84.1720) | 1 | Table II thermal BM3 with reciprocal quadratic compressibility |
| [Yamazaki et al. (2012)](https://doi.org/10.1029/2012GL053540) | 2 | Table 1 BM3 and Vinet, each with complete vibrational, anharmonic and electronic pressure |
| [Sakai et al. (2014)](https://doi.org/10.1016/j.pepi.2013.12.010) | 38 | New NaCl-B2 calibration; six pressure scales and three P1 reference-volume choices: 16 isothermal fits each for Fe and Fe0.9Ni0.1, plus five alloy thermal fits |
| [Dewaele et al. (2006)](https://doi.org/10.1103/PhysRevLett.97.215504) | 1 | Distinct Table I BM3 alternative; the existing Vinet thermal record is preserved |

These 43 published records extend `iron` and introduce `fe09ni01_hcp`; the new P5 calibration extends `nacl_b2`.
The latter is a substitutional hcp alloy with Fe/Ni occupancies 0.9/0.1 on
Wyckoff 2c, space group P63/mmc, and two atoms per conventional cell. Its
reference lattice is the measured F10N06_012 cell at 27.1 GPa and 300 K;
it is distinct from the extrapolated zero-pressure EOS volume. Existing
catalog defaults are preserved. No local regression coefficients replace a
published record.

## Primary datasets and pressure scales

Three new datasets preserve 348 observations:

- Yamazaki Table S1: all 207 rows, including Au volumes, Fe lattice constants,
  temperatures, pressures and quoted errors. Pressures retain **Tsuchiya
  (2003) Au**, unlike the later Fei-scale re-reduction in Zhang's compilation.
- Sakai Table 2: all 104 rows, comprising 27 pure-Fe and 77 alloy measurements
  (45 room-temperature and 32 heated alloy points). Separate P1/P4/P6/P7
  columns, missing cells, NaCl/MgO marker volumes and uncertainties are retained.
  P1, P2, P3, P4, P5 and P6 correspond respectively to NaCl calibrated against
  Matsui Pt, Fei Pt, Dorogokupets–Oganov Pt, Tange MgO, Yokoo Pt and Holmes Pt.
  P7 is the direct Tange MgO thermal scale. Low-pressure NaCl-B1 rows use
  **Brown (1999)**, as specified by Table 2's footnote.
- Brown Tables I–II: all 37 Us–up–P–density rows and their one-sigma errors.
  No shock temperatures are supplied. The source samples are approximately
  99% low-carbon steel, with up to 0.7% Mn; the source initial density
  **7.850 g/cm³** is preserved rather than substituted with pure-Fe density.

The existing Dewaele EPAPS dataset already contains all 63 rows, including
bcc/hcp phase and original-run selection flags. The BM3 diagnostic uses the
37 helium hcp rows marked as belonging to the original 2006 fit. Dewaele's
BM3 atomic V0 = 11.234(12) Å³ becomes cell V0 = **22.468(24)** Å³. The source
errors are 95% confidence intervals; a larger uncertainty copied in a later
comparison table is not substituted.

## Exact equation mappings and independent checks

Yamazaki's power-law Grüneisen expression is exactly the existing
`Dewaele2006` model with `gamma_inf=0` and `beta=q`. The Python and native
constructors now permit this nonnegative boundary; negative values remain
invalid. No approximate small positive replacement is used. Reference-subtracted
Debye, anharmonic and electronic terms all use the one-atom molar basis,
converted from the two-atom conventional cell.

Dubrovinsky's equation is implemented as
`K0(T)=1/(b1+b2*T+b3*T²)` and
`V0(T)=V0(300)*exp(alpha*(T-300))`, using the existing `AlphaKT` reciprocal
compressibility law. K0(300) is derived from the printed b coefficients.
The single tabulated expansivity is treated as constant; no higher-order
expansivity coefficients are invented. The independent calculation gives
K(211 GPa,300 K)=1111.44 GPa against the printed 1110 GPa and mean expansivity
9.205e-6/K at 202 GPa,5200 K against 9.13e-6/K. Original coefficient covariance
is unavailable, so no rigorous derived K0 error is claimed.

Sakai's five Table 7 thermal fits use the P4 BM3 reference and P7 hot pressures
(Table 8). Type 2 uses zero gamma_inf; Type 1 uses constant gamma, with gamma_inf
tied to gamma0 and theta proportional to V^(-gamma0), following the limiting
relation stated in Section 3.5. The otherwise untabulated q in the Table 6 Type 1
heading is not treated as an extra fitted parameter. The predicted Table 8
volumes at 329 GPa,5000 K agree within 0.014 Å³ of the printed values. These
states are explicit extrapolation checkpoints, outside thermal observations.
Yamazaki's two 330 GPa,6000 K density checkpoints agree within 0.025 g/cm³ of
the source values; rounded coefficients need not reproduce every final digit.

All 43 equations are independently evaluated using NumPy and Gauss–Legendre
Debye quadrature, compared with the executable native records, and checked
through `.eosmat` round trips. The alloy structure and all tabulated unit-cell
geometries are also checked.

## Refits and remaining limits

`scripts/reproduce_iron_source_papers.py` regenerates the full numeric report in
`docs/data/iron-source-papers-reproduction.json`. Its results feed the normal
primary-refit ledger. Every fit specifies its selected rows, pressure column,
free parameters and objective.

Brown's linear and quadratic Us–up coefficients reproduce at their printed
precision. The linear RMS with denominator N is 60.8 m/s; with N−2 it is
62.5 m/s, near the paper's 62 m/s. The quadratic RMS is 38.84 m/s against
39 m/s. The existing catalog supports a linear Hugoniot, so the quadratic
coefficients and their reproduction remain explicit source metadata and audit
results rather than being disguised as another linear model. The executable
linear record is restricted to the below-200-GPa solid subset used by Fei
(2016); this operational bound is **not** a measured phase boundary. The
published linear regression still uses all 37 points to 442.1 GPa. Neither
that full pressure range nor the quadratic model is presented as an hcp
isothermal or general thermal EOS.

Dewaele's conditional BM3 refit agrees within reported errors. For Sakai,
independent unweighted pressure fits conditional on the published g-G V0 are
possible for P1/P4/P6. P5 pressures are reconstructed from the paper’s new
NaCl-B2 BM3 calibration (V0=38.34(4.69) Å³, K0=45.18(48) GPa, K0′=4.22(2))
and the measured NaCl volumes, preserving the separate Brown B1 pressures.
The unusually large printed V0 uncertainty is retained. The original paired
2011 NaCl/Pt observations are absent, so this calibration itself is not refitted.
P2/P3 pressures are not tabulated; those published
records remain executable but lack a complete independent scale re-reduction.
The preceding g-G regression, weighting and covariance are not reproduced.
Sakai's five thermal fits yield approximately 7.1 GPa RMS against the tabulated
hot points using the documented P4/P7 combination, versus about 1.1 GPa reported
in Table 7. Their Table 8 checkpoints nevertheless agree. This unresolved
source-reduction/fit discrepancy is retained. The thermal Type 4-1 diagnostic
reaches its declared gamma0 upper bound of 8; it is not an unconstrained
parameter determination. Bounds are included in the machine-readable report.
The discrepancy is not corrected by changing
pressures or replacing source coefficients.

Yamazaki's 207-row unweighted refits yield slightly lower pressure RMS than the
published curves (published RMS 0.661 and 0.690 GPa), but do not recover all
parameters within their quoted errors. The original weights and parameter
covariance are absent. Dubrovinsky's 109 new plus 79 earlier P-V-T observations
are not tabulated in the attached four-page paper; its dhcp phase discussion
explicitly says that a dhcp EOS has not been determined. No dhcp EOS or invented
188-row dataset is added.

These additions improve the available primary inputs for Fei (2016), but they
do not establish parity for Fei's combined regression. Remaining gaps include
Dubrovinsky's actual observations, Fei's exact pressure intercalibration and row
selection, Brown shock temperatures used in 2016, fitting weights, and the
unpublished matched Pt measurements for Fei Table S2.
