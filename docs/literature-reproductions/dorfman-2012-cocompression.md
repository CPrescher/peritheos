# Dorfman et al. (2012) co-compression audit

This audit covers the six published Au, Mo, and Pt Vinet records derived from
Dorfman et al. (2012), DOI
[`10.1029/2012JB009292`](https://doi.org/10.1029/2012JB009292). The record
coefficients remain the printed Table 2 values. They are not replaced by the
diagnostic refits below.

## Primary evidence and custody

- The article was checked from the publisher page and the author-hosted PDF.
  The inspected PDF has SHA-256
  `1c9fe0f884fdb61f9d52d60308276d63dd434f27e59bcea2c74471566d0ec5e0`.
- The numerical observations are in the publisher's auxiliary PDF
  `jgrb17272-sup-0002-txts01.pdf`, Tables S1-S6. The inspected file has SHA-256
  `d7f7214cae4dabb3aae8a8585df4af7ada7990cb36b83a42ed7ed97e9eef7f2a`.
- A lossless long-form transcription of every nonblank numerical table cell is
  bundled as
  [`dorfman-2012-tables-s1-s6-cocompression.csv`](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/dorfman-2012-tables-s1-s6-cocompression.csv).
  Exact printed tokens are retained beside normalized numeric values, and the
  CSV is checksummed in the material metadata.
- Peritheos dedicates the rights its contributors hold in that factual CSV
  transcription, normalization, column naming, and arrangement under CC0 1.0.
  The scoped
  [license notice](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/dorfman-2012-tables-s1-s6-cocompression.LICENSE.md)
  expressly excludes the article, publisher PDF, and third-party rights. The
  publisher PDF itself is not redistributed.
- The publisher correction, DOI
  [`10.1029/2012JB009800`](https://doi.org/10.1029/2012JB009800), changes only
  the Figure S3 caption and does not change the EOS equation, table, or rows.

The reproducible parser requires the exact publisher auxiliary PDF, verifies
its checksum, extracts the printed table cells with `pdftotext -tsv`, and can
regenerate both the bundled CSV and the fit summary. There is no plot
digitization or hand-entered observation table. The derived fit result is
stored in
[`docs/data/dorfman-2012-cocompression-refit.json`](../data/dorfman-2012-cocompression-refit.json).

## Data scope and normalization

Tables S1-S6 yield 165 aligned observation rows, 368 conventional-cell volume
values, and 241 unordered material pairs. The runs represented are ANN11,
MMH6, MMH12, MNN11, MNN13, PMH6, PMN3, PMNN7, and PNH7. Three MgO entries are
marked as determined from the (200) reflection alone; they are retained because
the source does not instruct readers to exclude them.

Article Table 1 also lists run AN012 (Au + NaCl, no pressure medium, 1-164 GPa),
but no AN012 rows occur in the auxiliary PDF. Thus the official numerical
deposit is not the complete run inventory described by the article. This does
not establish that AN012 was included in either Table 2 fit, or that its absence
causes the coefficient discrepancy.

All tabulated volumes are in angstrom cubed per conventional cell. The cell
normalizations are Au 4 atoms, MgO 4 formula units, Mo 2 atoms, B2 NaCl 1 formula
unit, Ne 4 atoms, and Pt 4 atoms. The fit is the 300 K reference isotherm.

## Published coupled fit

Equation (2) is the standard Vinet pressure:

```text
x = (V/V0)^(1/3)
P = 3 K0 (1-x) exp[3/2 (K0'-1)(1-x)] / x^2
```

Equation (3) minimizes all co-compression pressure differences together:

```text
chi^2 = sum_i (P_a,i - P_b,i)^2 / mean(P_a,i, P_b,i)
```

In this replay, every unordered pair present in a row contributes once; the
authors' exact pair inventory remains unrecovered. Printed pressure
columns are outputs of the Table 2 EOSs and are used only as a transcription
check, never as independent fit observations. Recomputing all 368 printed
pressures from the published fixed-metal coefficients gives a maximum absolute
difference of 0.1824 GPa and an RMS difference of 0.0439 GPa, consistent with
the table's printed precision.

All `V0` values are fixed. MgO is fixed to the 300 K Vinet part of Tange et al.
(2009) Fit 3: `V0=74.698 A^3`, `K0=160.6 GPa`, `K0'=4.37`; the source is DOI
[`10.1029/2008JB005813`](https://doi.org/10.1029/2008JB005813), Table 4. In the
fixed-metal solution, Au, Mo, and Pt `K0` are fixed while their `K0'` values and
both NaCl and Ne coefficients vary. In the free-metal solution, `K0` and `K0'`
vary for all five non-MgO standards.

The article reports formal coefficient errors but warns that they are
artificially small and recommends approximately 2% realistic uncertainties for
`K0` and `K0'`. No covariance matrix or fit code was recovered from the inspected
primary sources. The
parenthetical volume precision in Tables S1-S6 is preserved by the parser but
is not introduced as a new weight because Equation (3) specifies the weight.

## Independent result

The literal simultaneous refit converges in both source-defined constraint
modes, but it does not reproduce Table 2. In both modes the independently
optimized Equation (3) objective is lower than the objective at the published
coefficients.

| Mode | Equation (3) objective, published | objective, refit | Au `K0`, `K0'` | Mo `K0`, `K0'` | Pt `K0`, `K0'` |
|---|---:|---:|---|---|---|
| fixed metal `K0` | 19.01399 | 14.37466 | 167 fixed, 5.86027 | 261 fixed, 4.13810 | 277 fixed, 5.36541 |
| free metal `K0` | 15.45416 | 12.19784 | 178.713, 5.46833 | 255.460, 4.29510 | 294.309, 4.89884 |

For completeness, the simultaneously fitted auxiliary standards are:

| Mode | NaCl `K0`, `K0'` | Ne `K0`, `K0'` |
|---|---|---|
| fixed metal `K0` | 27.2342, 5.42500 | 0.986588, 8.57976 |
| free metal `K0` | 29.0987, 5.24252 | 1.11008, 8.36343 |

This is a non-parity finding, not a substitute parameterization. Plausible
unpublished inputs include the missing AN012 rows, unrounded volumes, a row
mask, or implementation details not stated in the paper; the available evidence
does not establish which explanation is correct. The exact remaining scientific
step is to obtain the authors' complete numerical input and fitting code (or a
source statement resolving the row selection and weighting), then rerun the
checksummed coupled fit.

## Unweighted first approximation

At the user's request, an unweighted baseline was checked on 12 October 2026.
Its objective is `sum_i (P_a,i - P_b,i)^2`, with one equally weighted residual
per available unordered material pair. It uses the same 165 aligned rows, 241
pairs, fixed volumes, fixed MgO reference, and two metal-constraint modes as the
archived replay. Only the mean-pressure denominator is removed. The
[diagnostic report](../data/dorfman-2012-au-unweighted-baseline-2026-10-12.json)
records source hashes, optimizer settings, versions, and results.

| Mode | Published Au `K0`, `K0'` | Unweighted Au `K0`, `K0'` | Unweighted objective at published solution (GPa squared) | Optimized objective (GPa squared) |
|---|---|---|---:|---:|
| fixed metal `K0` | 167 fixed, 5.88 | 167 fixed, 5.85617 | 1991.63898 | 1301.47840 |
| free metal `K0` | 167, 5.84 | 173.75091, 5.63894 | 1519.37413 | 1181.28222 |

Both fits converge, but neither reproduces the published Au coefficients.
These objective values cannot be compared directly with Equation (3)'s values:
the unweighted objective has units of GPa squared, while the printed objective
divides by mean pressure and has units of GPa. Its weight is `1/mean(P_a,P_b)`,
not an inverse experimental variance or `1/mean(P_a,P_b)^2`.

An unweighted fit is the baseline diagnostic; the printed objective remains a
separate reproduction target. Neither calculation resolves the incomplete run
inventory, the original fit selection, or the cause of the discrepancy. The
numerical source CSV and published material coefficients are unchanged.

## Au input follow-up, 12 October 2026

The user also asked whether Au `K0=167 GPa` could have remained fixed while Mo
and Pt `K0` were freed. Table 2 marks Au's first value as fixed (`f`), but the
second as `167 +/- 3 GPa` without `f`, under the free-metal heading. A
[targeted constraint check](../data/dorfman-2012-au-k0-constraint-check-2026-10-12.json)
on the same rows gives Au `K0'=5.86414` with Equation (3) and `5.85896` without
weighting when only Au `K0` stays fixed among the metals. Thus that constraint
does not recover the printed free-mode `5.84` to displayed precision. The
original constraint implementation still requires author confirmation.

No additional numerical observations, fit code, unrounded inputs, or authoritative
row-selection specification were recovered in this follow-up. Both Au records
therefore retain the existing non-parity finding and published coefficients.
The [source-search report](../data/dorfman-2012-au-source-input-search-2026-10-12.json)
records the routes, queries, attachment inventory, access failures, selected
repository metadata, response hashes, and local verification. Its source
metadata does not extend the CSV's CC0 dedication to third-party material.

The [publisher article](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JB009292)
lists nine attachments: a readme, the existing numerical PDF, three EPS figures,
and tab-delimited versions of article Tables 1-4. The latter are described as
article tables, not another measurement deposit. The numerical PDF remains the
only attachment described as unit-cell observations. Attachment contents could
not be downloaded anew: the web reader returned internal errors and direct
publisher requests returned HTTP 403. The existing CSV and archived result were
checked locally; the publisher PDF checksum was not reverified in this follow-up.

The [publisher correction](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JB009800)
was rechecked and concerns only the Figure S3 caption. The publisher attachment
listing still carries the old description; this is not evidence of a correction
to the fitted rows or coefficients. Crossref's records for the article and
correction expose no related-data or code link; their empty relation fields
cannot rule out other deposits or corrections.

The [Duffy publication archive](https://duffy.princeton.edu/publications) links
the article PDF, including the correction as its final page. The
[Dorfman publication list](https://sukidorfman.rocks/publications/),
[Princeton publication record](https://collaborate.princeton.edu/en/publications/intercomparison-of-pressure-standards-au-pt-mo-mgo-nacl-and-ne-to/),
[Carnegie publication list](https://carnegiescience.edu/epl/epl-publications/geophysical-lab-contributions-2012-2013),
and [CDAC publication list](https://cdac.phys.uic.edu/publications/cdac-publications-2012/)
provided no additional data/code link for this article. The exact
[EPFL repository record](https://infoscience.epfl.ch/handle/20.500.14299/89629)
was recovered through its public API: it is marked `metadata-only`, and the
bundle listing contains only `LICENSE`, with no article/data/code bundle.

Exact DOI, run-name, title, and correction searches, including searches scoped
to Princeton DataSpace, EPFL, Zenodo, Figshare, GitHub, and EarthRef, yielded no
usable additional inputs. These were indexed-web searches, not exhaustive
repository inventories. A possible further source is Dorfman's 2012 Princeton
dissertation, *Effects of iron enrichment on the chemistry and physical
properties of deep lower mantle silicates*, listed by the
[Princeton Mineral and Gem Collection](https://minerals.princeton.edu/about-collection/publications).
Its full text was not obtained: direct DataSpace search returned HTTP 401, and
an exact-title EPFL API query returned zero records. No claim is made about the
dissertation's contents or private, restricted, or unindexed author files.

The bundled Au observations all belong to ANN11: 22 Au volume cells across 26
aligned rows. Au is blank at row-in-run 10, 16, 23, and 26; NaCl B2 is blank at
row 2. The reproducer uses the pairs actually present and does not fill these
blanks. Six existing Dorfman tests passed, including replay of both archived
constraint modes from the CSV. No alternative objective or row mask was tried
during the source search; the later user-requested unweighted baseline is
documented separately above.

### Exact questions for the authors

The subsequent AN012 exclusion check found a consistent pattern: Table 1's nine
He/Ne-medium runs are exactly the nine runs present in the CSV; AN012 is both
the only no-medium run and the only absent run. Article paragraphs 9 and 31
frame the study and pressure-scale result around noble-gas/quasi-hydrostatic
compression. Figure 1 and the Au part of Figure 11 identify ANN11. These are
consistent with intentional exclusion of AN012, but do not establish the
original fit mask or the reason for its exclusion. Paragraph 13 describes a
simultaneous fit of the material data sets without a per-run inclusion list;
Figure 6's caption likewise supplies no per-run selection. No explicit AN012
exclusion statement was found in the inspected article or correction. AN012 is
therefore an absent experimental run, not a demonstrated missing fit input,
and must not be assumed to explain non-parity. The first author question below
remains necessary to confirm this distinction.

1. For each Table 2 constraint mode, can you supply the complete, unrounded
   simultaneous-fit input with run/point identifiers, phase assignments, and
   inclusion flags? Was AN012 included, excluded, or used only for another
   analysis? If included, please supply its paired Au/NaCl measurements and
   specify the treatment of NaCl B1 versus B2.
2. Which rows and material pairs entered each fit? Please specify the handling
   of ANN11's blank cells, the three single-(200)-reflection MgO values, any
   pressure or stress exclusions, and whether literature co-compression data
   shown in the figures also entered the objective.
3. Can you supply the fitting code or the exact residual definition: all three
   unordered pairs for a three-material point or another pairing; the mean
   pressure in Equation (3) evaluated from current model parameters or held
   fixed; and any additional uncertainty weights, correlations, or penalties?
4. What exact fixed volumes and MgO parameters, free/fixed parameter flags,
   bounds or priors, optimizer settings, and stopping rule were used? Please
   include unrounded best-fit coefficients and the final objective for each
   mode so that the 167 GPa / 5.84 free-Au result can be checked against the
   archived literal-equation diagnostic (178.713 GPa / 5.46833).
5. Can you supply the covariance matrix and the rule used for formal errors
   (including any residual or degrees-of-freedom scaling), separately from the
   recommended approximately 2% realistic coefficient uncertainties?

These questions have not been sent. No material records, observations, source
coefficients, reproducer, archived coupled-fit result, or shared generated
ledgers/manifests were changed by the source-search follow-up. That stage wrote
only this audit and its Dorfman-specific source-search report; a later documentation build
can publish them without refitting or reconciling generated numerical data.

## Reproduction

With a locally obtained copy of the publisher auxiliary file:

```console
python scripts/reproduce_dorfman_2012_cocompression.py \
  --source-pdf /path/to/jgrb17272-sup-0002-txts01.pdf \
  --output docs/data/dorfman-2012-cocompression-refit.json \
  --dataset-output peritheos/data/datasets/dorfman-2012-tables-s1-s6-cocompression.csv
```

The script refuses a file whose checksum differs. It emits aggregate counts,
source hashes, selection and weighting rules, transcription diagnostics, both
coupled fits, coefficient-by-coefficient comparisons, and the lossless
long-form CSV. The source PDF remains external.
