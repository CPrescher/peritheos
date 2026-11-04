# Prescher (2015) Fe7C3: recovered original fit input

The published BM3 coefficients achieve numerical parameter parity on the
recovered **28-row author fit input**, supplied by Clemens Prescher on
25 October 2026. The former Peritheos audit used different, CIF-derived
observations from Supplementary Table S4.

**Table S4 should not be used for scientific calculations, EOS fitting or
pressure calibration.** It is retained verbatim for provenance only because
its cell values differ from the revised author inputs and the archived
processing/export files contain inconsistent cell representations. The author
input is a separate dataset; it does not silently replace the published table.

Source publication: Prescher et al. (2015), *High Poisson's ratio of Earth's
inner core explained by carbon alloying*, Nature Geoscience 8, 220–223,
[doi:10.1038/ngeo2370](https://doi.org/10.1038/ngeo2370).

## Evidence and data identity

The checksummed resources bundled with the material are:

- `iron-carbide-fe7c3-prescher-2015-author-input.txt`: original author input,
  including all printed P, P_err, V and V_err tokens and original row order.
- `iron-carbide-fe7c3-prescher-2015-author-fit-export.txt`: original saved
  coefficient/error export, kept separately from observations.
- `iron-carbide-fe7c3-prescher-2015-author-fit.csv`: the same 28 observations
  with named columns and run identifiers.
- `iron-carbide-fe7c3-prescher-2015-source.json`: source hashes, original
  tokens, workbook identity, processing-log excerpts, internal-CIF checks,
  and the historical Table S4 reconstruction.
- `iron-carbide-fe7c3-prescher-2015-table-s4-compression.csv`: the unchanged
  published 27-row structural table, explicitly prohibited for scientific use.

All resource names above are relative to the installed package's
`peritheos/data/datasets` directory. Numerical reproduction requires no Dropbox
access or workbook-reading dependency.

The author archive's `Fe7C3 EOS V-P data.xlsx` contains distinct `series d` and
`series d new` sheets. Nine rows from `series a new` and 19 from `series d new`
match every coordinate and width in the author input. The same 28 rows occur in
`gesamt_sync new!E4:H31`, explicitly excluding series b. An additional rotated
series-a observation at 3.92 GPa explains the 28-versus-27 row count.

Separately, `create_table.py` reads structural CIF cell quantities into
`data_table.txt`. Every lattice length and volume in those 27 rows matches the
published Table S4 transcription. Both datasets use conventional-cell volumes
in Å³ with Z=8; an axis permutation or normalization change cannot explain
their differences.

## Where the two data paths diverge

March 14, 2013 CrysAlis logs explicitly contain both representations for the
same runs:

| Run | Author pressure (GPa) | Revised cell used by EOS (Å³) | Reduction/export cell (Å³) | Table S4 volume (Å³) |
| --- | ---: | ---: | ---: | ---: |
| d_p2 | 46.03 | 634.6(2) | 635.392994 | 635.4(2) |
| d_p9 | 101.215 | 571.8(3) | 564.726811 | 564.7(3) |
| d_p9_90 | 101.0 | 571.3(2) | 567.195500 | 567.2(3) |
| d_p11_90 | 111.675 | 562.4(3) | 563.484938 | 563.5(4) |

For p9, `REDLOGThu-Mar-14-15-44-03-2013.txt` reports the orthorhombically
constrained cell 4.0637(6), 11.013(7), 12.778(5) Å, V=571.8(3) Å³ at lines
1979–1982 and 2340–2343, matching the revised workbook. The same session
reports the other cell and V=564.726811 Å³ in its reduction output. The latter
cell is retained in the exported CrysAlis CIF and subsequent SHELXL CIF used
by the structural table. The log records a precise peak-table update using
three-dimensional profile analysis before displaying the revised cell.

The p9 CrysAlis CIF is internally inconsistent: its stored orientation matrix
and wavelength imply approximately 571.575 Å³, while its cell-volume field
prints 564.7(2) Å³. The p9_90 CIF similarly implies approximately 571.052 Å³
but lists 567.19(18) Å³. These calculations use the CrysAlis wavelength-scaled
UB relation V=λ³/|det(UB)|; finite printed wavelength precision is retained.
They diagnose inconsistent representations, not replacement measurements.

Thus the split already exists in the processing/export state. Interpreting
the structural export as retaining stale cell information is supported by the
refinement sequence, but the precise software behavior or missing update step
is unrecovered. The table already has the differing values in inspected
February 2014 drafts and retains them through the final supplementary files.

Other differences, including printed widths, remain intact. For d_p18_90 the
author input has 153.43 GPa whereas Table S4 prints 151.6 GPa. The acquisition
logbook contains distinct 151.4/151.7 and 153.0/153.86 GPa reading pairs, with
crossed-out run labels; the exact naming/transcription step is unresolved.

## Numerical parameter reconstruction

Fit a three-parameter BM3 using all 28 author rows, with no fixed coefficient:

\[
\sum_i \frac{[P_{\mathrm{BM3}}(V_i)-P_i]^2}
{s_{P,i}^2+[K(V_i)s_{V,i}/V_i]^2}.
\]

The denominator changes with the fitted parameters. Printed zero pressure
widths are preserved; positive volume widths keep it finite. Workbook
pressure widths often represent half the before/after reading difference,
so they are not established independent statistical standard deviations.
Covariance is reconstructed as inverse(JᵀJ) multiplied by the weighted
residual sum of squares divided by 25 degrees of freedom. The original
confidence semantics and exact solver implementation remain unknown.

| Parameter | Published | Saved author export | Reconstructed |
| --- | ---: | ---: | ---: |
| V0 (Å³/cell) | 748(1) | 748.09 ± 1.05 | 748.087675 ± 1.047441 |
| K0 (GPa) | 168(4) | 168.15 ± 4.11 | 168.154116 ± 4.110102 |
| K0′ | 6.1(1) | 6.06 ± 0.13 | 6.063294 ± 0.132209 |

All six archived coefficient/error numbers are reproduced at their printed
two-decimal precision. The published coefficients agree within their reported
error widths. This is numerical parameter parity on recovered author inputs,
not proof of exact original-workflow identity or independent physical validation.
Published coefficients and their source errors remain unchanged in the catalog.

The published rounded triplet has pressure RMS residual 0.776259 GPa on the
author inputs, compared with 3.587693 GPa on Table S4. Both are evaluated at
recorded volumes, with residual P(model)−P(input). The historical Table S4
refit (V0=745.114, K0=188.419, K0′=5.22807) is preserved as an obsolete
different-data diagnostic. Its earlier 2.04963 GPa optimizer residual used
adjusted volumes; its RMS at recorded Table S4 volumes is about 2.283 GPa.

## Reproduction and limits

The repository-only `scripts.reproduce_prescher_2015_fe7c3` module reads the
bundled author CSV and saved fit export. The primary-fit audit and generated
reconstruction plots use that input. Table S4 is excluded from `fit_datasets`
and retained as archival context with an explicit do-not-use qualification.

The source investigation inventoried 10,375 archive files, extracted 182
document/workbook/notebook files and 65 PDFs, inspected the relevant scanned
SXRD logbook pages, and searched instrument logs, CIFs, summaries, backups and
manuscript revisions. Proprietary project objects were not fully decoded and
raw detector frames were not reintegrated. No note explaining the exact
synchronization failure was found. Independent selection of physically
preferable cells would require further raw-data analysis.
