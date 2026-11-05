# Campbell (2009): independently evaluated buffer polynomials

`campbell-2009-independent-polynomial-grid.csv` is newly calculated by
`scripts/generate_campbell_2009_polynomial_tables.py` from the mathematical
formulas and numerical coefficients in the official supplement, Tables S4-S5,
rows 26-32. The source parameters and reproducible grid are preserved in
`campbell-2009-polynomial-source.json` and in each dataset declaration.

Attribute Campbell, Danielson, Righter, Seagle, Wang and Prakapenka (2009),
*Earth and Planetary Science Letters* 286, 556-564,
<https://doi.org/10.1016/j.epsl.2009.07.022>. The workbook was retrieved from
<https://ars.els-cdn.com/content/image/1-s2.0-S0012821X09004373-mmc2.xls>;
SHA-256 `b857b9e40f76389dd78da9b757c504b17607d69e4cdc55087fff6e2cb8566f30`.
The full workbook remains private source evidence and is not added to packages.

Peritheos contributors dedicate their held rights in the new calculated grid
and its arrangement under CC0 1.0. The independent generator is MIT. This is a
maintainer assessment of mathematical parameter evaluation; it does not grant
rights in the original workbook, source layout or tabulated integration results.
Source-data terms remain unspecified. No source prose, artwork or author code
is included. Downstream acceptance remains a separate assessment.

These are **polynomial approximations**, not measurements, exact recovered
author outputs, new EOS fits or a thermochemical reintegration. The original
3434 integration results remain unchanged and excluded from the optional
package. Using the printed coefficients gives IW RMS/max differences
0.01757084/0.05136113 log units and NNO 0.01304041/0.05389993. The IW maximum
exceeds the workbook's stated 0.046; that discrepancy is unresolved, not repaired
by a new fit or attributed to coefficient rounding without evidence.
All prior scientific qualifications, including deferred Fe/FeO EOS records,
remain in force. The declared calculation grid is 0-100 GPa and 1000-2600 K;
it does not establish phase stability or validated prediction throughout it.
