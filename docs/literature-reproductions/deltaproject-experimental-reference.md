# Delta-project experimental-reference parameterizations

The catalog contains 42 third-order Birch--Murnaghan comparison curves selected
from the CC BY 4.0
[Delta-project archive](https://archive.materialscloud.org/record/2023.133),
[doi:10.24435/materialscloud:5e-mv](https://doi.org/10.24435/materialscloud:5e-mv).
The exact selected coefficients are bundled in
`deltaproject-experimental-reference-parameters.csv` and the source hierarchy,
checksums, and recovery status in
`deltaproject-experimental-reference-source.json`.

## Scientific classification

These records are derived reference parameterizations, not raw pressure-volume
observations or 42 fits to unified experimental datasets. Each combines a
crystallographic volume, bulk modulus, and pressure derivative from
heterogeneous sources. Lejaeghere et al.,
[doi:10.1080/10408436.2013.772503](https://doi.org/10.1080/10408436.2013.772503),
approximately corrected the first two properties to a static-lattice 0 K
comparison state and left the derivative uncorrected.

The archived triplet constructs a BM3 curve over the historical Delta
0.94--1.06 volume-ratio integration window. That interval is not experimental
coverage, a phase-stability range, or an extrapolation claim. A unified
experimental refit is impossible because the archive supplies no common
row-level observations, pressure calibration, uncertainties, covariance,
weights, or selection rules.

## Original-data recovery

The property-level citations were followed upstream wherever a lawful source
copy was available:

- Knittle's complete element table and references recover 30 of the 32 cited
  default B0-prime values directly or by the stated arithmetic mean. Fe=4.6
  and Au=6.4 remain editorially unexplained. The exact reconstruction retains
  Knittle's reference numbers and phase notes.
- Anderson et al. (1990) prints 2.41(5) for Sr, whereas Knittle prints 2.47.
  Delta's 2.485 exactly averages Knittle's 2.5 and 2.47, so the immediate-source
  arithmetic is preserved while the upstream conflict is flagged.
- Fujihisa and Takemura's alpha-Mn paper reports a Vinet fit of B0=158+/-3 GPa
  and B0-prime=4.6+/-0.1, but publishes its observations only as a plot. Delta
  combines that B0 with a separate Knittle B0-prime=6.6. The plot was not
  digitized.
- Price, Rowe, and Nicklow's alpha-Sn elastic constants reproduce B0=42.5 GPa
  through `(C11 + 2*C12)/3`. Delta combines that alpha-Sn property with
  B0-prime=4.0 averaged from three beta-Sn rows, so the resulting curve is
  explicitly flagged as phase-mixed.
- Takemura's complete 28-row Os Table I dataset is already bundled. An
  unweighted fixed-V0 BM3 fit to its 27 increasing-pressure run-A/run-B rows
  gives B0=395.566 GPa and B0-prime=4.49887, reproducing the published central
  result from the rounded table. Delta uses only that derivative and combines
  it with a separately compiled Kittel B0.
- For Tl, archived `exp.txt` gives B0-prime=3.0, while the cited 2014
  supplementary table prints 5.8. The undocumented archive replacement
  matches Knittle's Tl row and is retained with an explicit conflict warning.

Only source-printed tables or author/publisher-hosted numerical observations
are treated as recovered data. Plots are not digitized, and inaccessible tables
are not inferred.

## Verification

`scripts/reproduce_deltaproject_experimental_reference.py` verifies the 42-row
transcription, per-atom-to-cell normalization, independent BM3 pressures, all
30 recoverable Knittle reductions, the alpha-Sn elastic derivation, and the
fixed-V0 Os refit. Supplying the extracted archive `history/history/exp.txt`
also verifies its SHA-256 checksum and every selected source coefficient.
