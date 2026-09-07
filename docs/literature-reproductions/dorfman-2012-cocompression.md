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
  [`dorfman-2012-tables-s1-s6-cocompression.csv`](../../peritheos/data/datasets/dorfman-2012-tables-s1-s6-cocompression.csv).
  Exact printed tokens are retained beside normalized numeric values, and the
  CSV is checksummed in the material metadata.
- Peritheos dedicates the rights its contributors hold in that factual CSV
  transcription, normalization, column naming, and arrangement under CC0 1.0.
  The scoped
  [license notice](../../peritheos/data/datasets/dorfman-2012-tables-s1-s6-cocompression.LICENSE.md)
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
deposit is not the complete run inventory described by the article.

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

Every unordered pair present in a row contributes once. Printed pressure
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
`K0` and `K0'`. No covariance matrix or fit code is published. The
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
