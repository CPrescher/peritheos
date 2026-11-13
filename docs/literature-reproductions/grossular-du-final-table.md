# Grossular: Du final Table 2 versus the preprint

Source comparison, 13 November 2026. **The final typeset grossular Table 2
matches the preprint in every row and printed uncertainty.** A change between
these versions cannot explain the remaining θ discrepancy in our grossular
reconstruction. The additional independent 2% Conrad error allowance remains
a diagnostic, with no change to any source error or published coefficient.

## Recovered source and comparison scope

Du, Clark and Walker (2015), *Thermo-compression of pyrope-grossular garnet
solid solutions: Non-linear compositional dependence*, American Mineralogist
100, 215–222, [DOI](https://doi.org/10.2138/am-2015-4752).

The eight-page final journal typesetting was recovered through the public
[author-uploaded ResearchGate copy](https://www.researchgate.net/publication/271207857_Thermo-compression_of_pyrope-grossular_garnet_solid_solutions_Non-linear_compositional_dependence).
The page identifies Wei Du as the uploader on 21 October 2016. The PDF has
the journal imprint, DOI and final pagination; Table 2 is on **p. 217 (PDF
page 3)**. This is an author-hosted copy of the typeset article, not an
authenticated publisher download or a demonstrated byte-identical publisher
file. The full article remains in the private source archive.

We visually compared the complete relevant pages with the
[Columbia preprint](https://academiccommons.columbia.edu/doi/10.7916/D8V415Z5/download):
Table 2's grossular block and footnote on manuscript/PDF p. 25, and the
temperature calibration on p. 7. The scope is **the 13 grossular thermal
rows and their definitions**, not identity of the entire paper or its other
composition/pressure tables.

| Compared item | Result |
| --- | --- |
| 13 temperatures, 296.1–859.6 K | Identical |
| 13 lattice parameters and parenthesized errors | Identical |
| 13 unit-cell volumes and parenthesized errors | Identical |
| 13 calculated α entries and parenthesized errors | Identical; not independent fit channels |
| NaCl temperature calibration, units and 0.2 °C error statement | Same definitions |

In particular, **393.7 K → 1666.0(1) Å³** and
**447.1 K → 1668.3(1) Å³** are present unchanged in the final table. They
are not preprint transcription mistakes. The separate
[final transcription](../data/garnet-angel-2022/grossular/du-2015-final-table-2.csv)
preserves all printed grossular entries; it does not overwrite the
[previous preprint transcription](../data/garnet-angel-2022/grossular/du-2015-preprint-table-2.csv).

## Calibration and uncertainty meaning

Final p. 217, Eq. 1, retains

\[
T_2=-3.6\times10^{-5}T_1^2+0.982T_1-0.565,
\]

where both variables are in **°C**, T₁ is the thermocouple temperature,
and T₂ is derived from the NaCl internal standard. The methods use T₂
throughout; the table already reports corrected temperatures in kelvin.
Applying the mapping again is not the stated source procedure.

Both table footnotes give a NaCl temperature error bar of 0.2 °C
(the same increment as 0.2 K), and describe parenthesized errors as estimated
standard deviations from unit-cell least-squares refinement. Neither version
establishes that all temperature errors are independent Gaussian deviations,
or supplies calibration covariance. The α column is calculated from the
thermal expansion analysis and is not an additional measured data series.
Both footnotes refer to Eq. 1 for α, although the derivative definition is
Eq. 2 in the body; this shared cross-reference issue does not alter the table
or temperature mapping.

## Consequence for the existing fits

All nine fields used from each preprint row, including parsed widths and
original lattice/volume tokens, are exactly equal in the new transcription.
The [machine-readable comparison](../data/garnet-angel-2022/grossular/du-final-table-comparison.json)
reports zero differing rows and zero numeric differences. Substituting the
final source into the [32 Du/Bosenick diagnostic cases](grossular-du-bosenick.md)
leaves every prepared observation, error, denominator, start and fit control
unchanged. All 219 fit rows and the three Conrad points remain included.
There is consequently no new fit to run from this source-version check.

This closes the final-table availability question at the stated scope.
It does not recover raw NaCl readings, calibration/run covariance, Angel's
exact normalization denominators or original EosFit errors and weights.
The existing conditional refits and unresolved θ finding remain unchanged;
this is not author-workflow reproduction or validation of a replacement EOS.
No author message was sent. The earlier pyrope correspondence remains
waiting for Ross's reply.

## Evidence and repeatable checks

The [source receipt](../data/garnet-angel-2022/grossular/du-final-table-source-evidence.json)
records private PDF paths, SHA-256 checksums, source versions, page locations,
retrieval limitations and visual checks. The previous fit receipts and all
ten preceding artifact manifests remain unchanged, including their historical
statement that the final table had not yet been recovered. See the new
[verification receipt](../data/garnet-angel-2022/grossular/du-final-table-verification.json).

Repository-only commands:

```bash
python -m scripts.compare_grossular_du_final_table
pytest -q -W error tests/test_grossular_du_final_table.py
```
