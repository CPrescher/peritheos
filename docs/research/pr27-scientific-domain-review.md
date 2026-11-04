# PR 27: scientific transcription and provenance review

Reviewed 4 November 2026 against revision `55dadf85` and base
`e2d4b163`. This review checks the actual PR difference; the finite-reference
model and pressure-conversion implementation already shipped in earlier
releases. The remaining scientific additions are six cementite EOS records,
recovered Coesite and Fe7C3 inputs, and qualified research evidence.

**No coefficient or observation-token discrepancy was found in the comparisons
below.** Current-summary errors were found and corrected. Source observations,
published errors, fixed/adopted parameters, original files, historical fits and
use restrictions are preserved. This result does not establish physical
validation or full reproduction of the original author workflows.

## Sources and comparisons

The [machine-readable review](../data/pr27-scientific-domain-review.json)
records the checked source-file hashes, comparison counts and limitations.
The independent review calculations and page renders remain in the local
research archive rather than adding publisher files to the repository.
No workbook was recalculated or saved.

| Case | Source examined | Result and retained scope |
| --- | --- | --- |
| Li (2002), DOI [10.1007/s00269-001-0224-4](https://doi.org/10.1007/s00269-001-0224-4) | Preserved final four-page article, visually inspected Table 1 and equation; original token CSV | All 14 rows and 155 central/error tokens agree. The reference-card row is not a measurement; the ambient measurement fixes V0; 12 nonambient rows and ten paired neon readings remain distinct. K0=174(6), K0-prime=4.8(8), V0=155.28(11) are retained. Ruby B=7.715 and averaged ruby/neon pressures are distinguished from later recalibrations. |
| Ono–Mibe (2010), DOI [10.1016/j.pepi.2010.03.008](https://doi.org/10.1016/j.pepi.2010.03.008) | Final article Table 1, Table 3 and Figure 4; preserved token CSV | All 14 P/a/b/c/V rows agree. The source is experimental BM3, 167/6.7, below 35 GPa; actual included maximum is 34.3 GPa. Measured V0=155.4(2) is explicitly adopted, without inventing a source fixed-V0 constraint. The 0 K calculations and later Vinet recalibration remain separate. |
| Sata (2010), DOI [10.1029/2009JB006975](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009JB006975) | Final publisher HTML Table 1, Table 6 and paragraph 15; independently preserved 2009 manuscript Table 1 | All nine Fe3C/MgO rows and 72 numeric/error tokens agree. Published 9.341 A3/atom becomes 149.456 A3 for 16 atoms, without replacing it by rounded 149.46. Pr=0.0(1.6), Kr=290(13), K-prime=3.76(18), the three overlapping Li rows and the unresolved MgO-reducer precision remain explicit. The manuscript is corroboration, not the publication of record. |
| Takahashi (2019), DOI [10.3390/min9120744](https://doi.org/10.3390/min9120744) | Preserved final article Tables 1–2 and Equations 1–6; official supplement XLSX | All retained raw cell values/formats and all 60 normalized rows agree. The cell P5 value 121.3 has format 0.00, so Q5 `(8)` means 0.08 A3. Both pressure-scale coefficient/error sets, fixed theta0=314 K, Z=4 and four atoms per formula are retained. The integrated Gruneisen Debye law matches Equation 6. The 2600 K supplementary row and SiO2-marker rows are explicitly included; this does not recover the authors' exact selection or weighting. |
| Hu (2019), DOI [10.1029/2019GL084545](https://doi.org/10.1029/2019GL084545) | Official [PANGAEA workbook](https://doi.pangaea.de/10.1594/PANGAEA.905474), all Sheet1 A1:I8 source cells | All 70 retained header/shot cells agree, including individual initial densities and calculated temperatures. Published c0=4.627(73) km/s and slope=1.614(28) are retained over all six shots, 85.6–245 GPa. The mean porous bulk precursor density and bulk mass-equivalent V0 are not diffraction-cell measurements. The empirical principal Hugoniot is not a 300 K isotherm or proof of solid phase identity through 245 GPa. |
| Coesite, Bykova (2018), DOI [10.1038/s41467-018-07265-z](https://doi.org/10.1038/s41467-018-07265-z) | Original author DAT and original XLSX, including cached values, shared formulas, number formats and source addresses | All 22 input rows, 100 separate calculated-curve rows and retained workbook cells agree. There are 21 measured entries and one derived ambient EOS anchor. Literal DAT replay gives approximately 34.204033/102.729498/3.024481. The documented precision reconstruction retains the demonstrated one-decimal current-study V/Z while using workbook P/sigma_P and literature precision; it gives 34.203783/102.750605/3.022928. This reproduces the communicated saved curve at printed precision. Using full precision for *all* workbook volumes is a different control and does not reproduce that curve. Neither agreement proves original solver/covariance identity. |
| Fe7C3, Prescher (2015), DOI [10.1038/ngeo2370](https://doi.org/10.1038/ngeo2370) | Preserved original 28-row author text, saved fit export, normalized CSV, and preserved cell-addressed workbook extraction | All 112 original input values agree with the normalized CSV and with extracted `gesamt_sync new!E4:H31` cached values. An independent NumPy/SciPy effective-variance calculation reproduces all six printed export coefficients/errors: 748.09/1.05, 168.15/4.11, 6.06/0.13. Zero pressure widths remain zero; confidence semantics and the exact author solver remain unknown. Table S4 is retained separately under its existing provenance-only restriction. |
| Prescher (2012), DOI [10.1103/PhysRevB.85.140402](https://doi.org/10.1103/PhysRevB.85.140402) | Final user-supplied supplement, all 19 numbered Table S1 pressures/volumes | The 19 states agree, including p16=133.25 A3. There is no measured ambient row in this table. The current primary comparison is the 19-row fixed-adopted-V0 replay; its K0 error is not reproduced. The 22 graphical markers and split-fit diagnostics remain separate historical evidence. |
| Finger (1979), American Mineralogist 64, 1002–1009 | Preserved scanned article and pages 1004–1005 | Table 2 Fe2SiO4 volumes 558.7(2), 553.5(7), 547.7(8) A3 at 0.001/20(1)/40(1) kbar agree. The text fixes K-prime=4 and gives 1.96(6) Mbar, or 196(6) GPa. The curve remains an explicitly derived static reconstruction with a measured adopted ambient volume; no author-fitted intercept or thermal validation is inferred. |

The new nonmagnetic cementite structure was separately checked against
Prescher Table S1 p19 (47.4 GPa, 293 K): all three lattice values and nine
fractional coordinates follow the declared cyclic axis mapping, with no origin
shift. Ideal full occupancies give Fe12C4 and Z=4; they are not described as
refined occupancies. The finite-pressure structural volume is not an EOS
ambient-volume measurement. The structural fields of all 233 pre-existing
material documents are unchanged.

## Preservation exceptions and independent protection

The supplemental merge bindings permit exact, reviewed metadata transitions;
they are not evidence of scientific correctness by themselves. Independently
of those exceptions, this review compared all 679 existing records with the
PR base: model, parameters/errors/fixed flags, thermal definition, reference
temperature and parameter-error confidence are unchanged. Existing dataset
rows, resource declarations and column definitions are unchanged as well.
Source-resource byte checks remain in the frozen release-source test.

The new `scientific_source_snapshot` maintenance check freezes all current
scientific coefficients, units, crystallographic identity/coordinates,
reference-state and thermal conventions,
source arrays/resource bytes, observation views and uncertainty semantics.
Its reviewed fingerprint is separate from the merge metadata exceptions.
Changing a numerical source value cannot be excused by merely accepting a new
document metadata fingerprint.

## Corrections made

- Correct the current cementite summary and changelog: the primary Prescher
  comparison uses 19 numbered states and does **not** reproduce the K0 error.
  The earlier 22-marker result remains preserved and explicitly historical.
- Correct the summary's obsolete claim of three installed Prescher records:
  only the full-range record is installed; the splits remain audit evidence.
- Qualify the Fe7C3 synchronization explanation: different cell representations
  are documented, but the original cause is unresolved.
- Update stale API-stability discovery totals to 234 documents, 231 executable
  materials, 630 source choices and 682 executable records.
- Name the actual nested native-v4 thermal fields in the merge-preservation
  review. No material coefficients, arrays or scientific dispositions change.

## Remaining limits

The original wider Fe7C3 workbook/instrument archive is absent at its historical
path in this checkout. Its original input and fit-export copies are available;
the preserved workbook extraction identifies original workbook SHA256
`75b93cac2a5a912af33ddd7e9f10a4195b942b95f18d4044fa75cfdfec0c4b8b`.
This review verifies those preserved copies and extraction, not fresh access
to that missing XLSX, raw detector frames or complete instrument logs.
The Coesite objective is attributed to author communication as supplied by
the user; no independently authenticated correspondence file was available.
These boundaries prevent a claim of complete original-workflow reproduction.

All original observations, transcriptions, coefficients, historical diagnostics
and existing recommendation holds remain available. The review supplies no
replacement EOS, resolves no previously unresolved physical phase boundary,
and does not promote a conditional numerical comparison to experimental
validation. Release readiness additionally requires the full compatibility
workflow on the exact intended release commit.
