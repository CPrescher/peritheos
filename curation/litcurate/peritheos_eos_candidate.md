Extract every distinct equation-of-state parameterization reported in this
paper. The output is a Peritheos review candidate, not an executable record.

Use only statements, equations, tables, captions, and supplements in the
provided document. Do not infer a conventional equation family from familiar
parameter names. Do not convert units, repair values, choose between conflicting
source locations, or fill missing parameters from another composition or paper.

Create a separate entry for every scientifically distinct combination of
phase, composition, sample, equation, reference state, pressure scale, and fit.
Classify a value as `source_reported` only when the paper presents it as a new
determination. Values copied from earlier literature are `citation_reported`.

For the equation:

- record the printed name, order, equation number, and complete equation text;
- distinguish a zero-temperature/static cold curve from a finite-temperature
  reference isotherm;
- include every coefficient required by the printed equation, including
  higher pressure derivatives and thermal or electronic parameters;
- use `existing_family_candidate` only for an explicit conventional BM2, BM3,
  BM4, Vinet, Murnaghan, modified-Tait, natural-strain, Holzapfel, or other
  family already named in the paper;
- use `unresolved` when the paper does not identify the formulation.

For each parameter, copy the reported value and unit verbatim. Record whether
it was fitted, fixed, adopted, measured, or derived, together with uncertainty,
uncertainty convention, and an exact table/equation/page locator.

Also extract the material's phase and composition, structure and space group,
formula units represented by the reported volume, reference pressure and
temperature, experimental/computational range, pressure calibrants, primary
data locations, and data license when stated.

Evidence excerpts must be short and paired with precise locators. Put every
contradiction between prose, equations, tables, and supplements in `conflicts`.
Put every fact still needed for an executable Peritheos record in `blockers`.
Use null or `unknown` when the paper does not state a value. Return
`{"eos_entries": []}` when no numerical EOS parameterization is reported.
