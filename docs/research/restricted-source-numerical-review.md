# Sakai, Ismailova and iridium numerical redistribution review

Reviewed 5 November 2026 against `294a651f`. This review concerns what can
reasonably be redistributed in the optional **`peritheos-datasets`** package
(`peritheos[datasets]`) in general. The core `peritheos` package separately
supplies calculations, structures and EOS parameters without observation values;
it is the intended starting point for consumers such as Dioptas.

**Admit 21 existing numerical entries and retain one hold.** The optional
selection becomes **448 of 453 material dataset entries**, backed by 388 CSV
resources. No source grant has changed. This is a reasoned maintainer assessment
of the rights applicable to these particular numerical extractions, not a claim
that the original publications are freely licensed or a legal guarantee.

| Source / content | Entries | Decision and basis |
| --- | ---: | --- |
| Sakai Tables S1–S7 and overlapping Au calibration subset | 8 | Admit primary numerical measurements and functional run labels |
| Sakai Table S9, nine materials | 9 | Admit calculated numerical checkpoints, not observations |
| Sakai Table S14, calculated Cu-scale columns | 1 | Admit mathematical scale-conversion results |
| Sakai Table S8 literature/calibrant lineage | 1 | Hold recent assembled literature index |
| Ismailova Table S2 numerical excerpt | 1 | Admit four selected pressure/lattice/volume states |
| Anzellini iridium Tables S1–S3 | 1 | Admit 122 primary experimental numerical rows |
| Anzellini iridium Figure 4 | 1 | Admit 107 quantitative marker-center readouts, not artwork |

The [bound machine-readable receipt](../data/restricted-source-numerical-review.json)
records all 22 exact declarations/resources, source-copy hashes, decisions and
limitations. It is embedded as `restricted_source_review` in the installed
`redistribution-review.json`. Original CSVs, declarations, source-rights flags,
coefficients and scientific classifications remain byte-for-byte unchanged.
Original source PDFs/XLSX, rendered review pages and a pre-change snapshot are
preserved privately; they are excluded from both distributions.

## Rights basis and limits

[Sakai (2025)](https://doi.org/10.1038/s43246-025-00792-5) and
[Anzellini (2025)](https://doi.org/10.1038/s43246-025-00963-4) carry
**CC BY-NC-ND 4.0** source notices.
[Ismailova (2016)](https://doi.org/10.1126/sciadv.1600427) carries
**CC BY-NC 4.0**. These declarations remain explicit, including all 22
`restricted_source` inventory flags. None is changed to CC0, CC BY or MIT.

The assessment distinguishes numerical facts/mathematical results from protected
expression and database rights. The [CC legal code](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.en)
limits its grant to applicable copyright/similar rights; it does not create
rights in otherwise unprotected content. Conversely, its Section 4 covers
qualifying database extraction. The [CC database guidance](https://creativecommons.org/faq/#which-components-of-databases-are-protected-by-copyright)
and [German database law](https://www.gesetze-im-internet.de/urhg/__87a.html)
require considering collection rights separately from individual values.

These admissions do **not** use NC permission as a basis for commercial
redistribution or ND permission as a basis for sharing protected adaptations.
They rest on the specific numerical content reproduced and the inspected
source methods. There is no article prose, copied caption, figure artwork,
source-sheet layout or author software in the admitted resources. Ordinary
scientific fields, source tokens and run/phase labels remain for reproducibility.
Contributor CC0 dedications apply only where actually made and within rights held.

The primary experimental tables describe creation of measurements and their
in-measurement processing. No separate substantial collection investment
independent of that creation was established for the admitted resources. This
is an assessment of the evidence, **not proof that qualifying verification or
presentation investment cannot exist**. Whole numerical-table reuse is
acknowledged; neither a small row count nor experimental origin automatically
exempts a collection. Acquisition contracts were not exhaustively audited.
A valid third-party restriction is not overridden by this review. A separate,
file-specific CC BY 4.0 or CC0 grant covering relevant held rights would provide
stronger evidence for general reuse. Distributors can identify and exclude
assessed resources through the installed ledger; their policies may differ.

## Sakai: new experiments, calculated grids and an assembled index

The fresh official workbook matches the earlier checksum:
`c559b29fea9c793831e4bc6deef4d62712b5e78821fc420969f709cbd1060ff9`.
All nine preserved raw worksheet objects match the fresh workbook exactly.
Regenerating the 326-row Au calibration input reproduces its CSV bytes exactly.
That file contains runs S1/S2/S3/S4/S7 from the present experiments; it is **not**
a collection of older literature measurements. The seven individual measured
tables contain 493 paired-volume states in total. The aggregate and individual
entries overlap and must not be treated as independent observations.

The article's Methods, p.7, describes the diffraction/volume generation and
stress analysis. Numerical volumes/errors, literal OOXML tokens, cell locations,
run/apparatus identities and missing errors remain intact. There are no
invented pressure observations or fit weights. Calibration reconstruction,
uncertainty parity and full author-workflow limitations remain unchanged.

Table S9 contains formula-generated checkpoints from Table 2 parameters.
All 234 volume/ratio cells and the 225 noninitial pressure cells match the
fresh workbook at serialization precision. An existing normalization is now
explicit: each CSV starts at **0 GPa**, while workbook A6 contains
**−0.0002242214081267206 GPa**. Ambient volume ratios also retain the source's
rounding deviations from one. The original workbook and CSV values are both
preserved; exact source-byte equality is not claimed for this normalization.
These are calculated checkpoints, not independent fit observations or a newly
regenerated replacement dataset. Existing equation-reproduction tolerances stay.

Table S14 contributes only the volume and two calculated pressure-scale columns
(E/F/G), not the historical-observation columns A/B. It remains a calculated
scale-conversion input; no source workbook formulas or layout are distributed.

**Table S8 remains excluded.** Its 35 rows assemble pre-existing literature
sources, calibration relationships, selection and notes. This presents a
plausible separate obtaining/verification/presentation and selection issue that
the inspected evidence does not resolve. To include it, obtain an appropriate
file-specific grant from the relevant rightsholders, or reconstruct factual
links independently from primary sources with a new selection/presentation and
an explicit collection-rights review. Reformatting the current index alone
would not resolve that issue. This is not a claim that each bibliographic fact
is protected.

## Ismailova: four selected numerical checkpoints

The fresh institutional supplementary PDF matches the preserved SHA-256
`5ab1c1c60e7d5486299cf5f0653d05b6b8bed6957a382e839ba292e2456ecd34`
and MD5. All pressure, a/b/c/V and reported temperature tokens were compared
with rendered Table S2, PDF p.18. The CSV extracts functional numerical fields,
not the complete crystallographic refinement table, layout or prose.

The missing first temperature remains missing; reported high-temperature
conditions are not substituted for the 300 K EOS regression temperature.
All printed parenthesized errors remain beside their normalized values.
These four checkpoints are not the complete compression/decompression series.
The missing row-wise neon observations still prevent recalibration; the EOS
remains `not_refittable`. Admission makes no new scientific validation claim.
The earlier blanket commercial-use wording in the source audit is superseded
by this scoped assessment; the original source notice is retained.

## Iridium: experimental tables and quantitative digitization

The fresh supplement matches SHA-256
`f661efccbb391e92cb6adf3070c4a21a01f5fba0a2306075175948f344980f73`.
Tables S1–S3 contain the authors' eight laser-heating runs. The 122-row CSV
retains the source MgO/KCl availability, both pyrometry measurements, printed
Delta-T values and missing uncertainties. The mean temperature is identified
arithmetic, not an additional measurement. Table structure/calibrant roles and
source tokens were inspected; this is not represented as a second independent
manual transcription of every cell.

Rerunning the maintained Figure 4 vector-center extractor on the exact source
PDF reproduces the 107-row CSV **byte-for-byte**. The CSV contains only numerical
centers, mapped pressure/volume coordinates and functional source labels.
It contains no circle/vector paths, image pixels, fitted curves, legend or caption.
The 0.75 GPa and 0.08 Å³ bounds remain graphical curation bounds, not experimental
errors or source regression weights.

The 91 black 300 K readouts refer to Monteseguro (2019); 16 dark-red readouts
refer to present-study 830 K observations. The fresh
[Monteseguro source PDF](https://doi.org/10.1038/s41598-019-45401-x) matches
its prior hash and prints CC BY 4.0. That is ancestry evidence, **not permission
for the 2025 replot**. The present extraction is admitted on the numerical-coordinate
assessment, not by asserting that the entire 2025 figure became CC BY.
Scientific reconstruction limitations and published EOS coefficients remain.

## Remaining exclusions

Five material entries remain held: Sakai Table S8 (one CSV) and four views of
the original Campbell integration-results collection (one shared CSV). The two
Barker supporting-study entries remain outside the material-dataset interface;
that is not a restrictive-license finding. No additional exclusion for Debian
is presumed here: core `peritheos` can be packaged without this optional dataset
archive at all. The split remains unreleased and published 0.17.4 is unchanged.
