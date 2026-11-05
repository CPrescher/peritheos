# Independently calculated reference patterns

Reviewed 5 November 2026 against `42abd741`. All ten inherited beamline
reflection arrays are replaced by independent kinematic X-ray calculations.
The original **452 rows** and complete original EOSMAT records remain privately
preserved, with checksums. The current arrays contain **1,154 calculated lines**.
These are calculated reference sticks, not measured peak strengths or a
reproduction of the inherited card selection.

| Material | Calculated lines | Structural input and limits |
| --- | ---: | --- |
| fe3s | 123 | CC0 COD 9007338 schreibersite; Fe/Ni -> Fe, P -> S. Approximate Fe3P-type topology proxy; no Fe3S refinement. Root atom sites now use the same new proxy. |
| feh2 | 36 | Pépin (2014) numerical I4/mmm model: refined Fe z=0.853, DFT H 4c/4d. Nominal occupancy; retained reference cell differs from the 67 GPa source cell. No separate source-publication grant established. |
| iceviii | 41 | CC0 COD 1566661, D2O ice VIII. X-ray intensities use the H electron form factor for D; these are not neutron intensities. |
| majorite | 410 | CC0 COD 9010029, including Mg/Si disorder and the source's explicit symmetry operations/origin. |
| mg7si2o8 | 241 | CC0 COD 9006822 Phase A, including H sites; retain the existing Horiuchi reference cell rather than the slightly smaller Kagi CIF cell. |
| mgfe60o | 10 | Independently authored ideal rocksalt model: statistical Mg0.4Fe0.6 / O. Reference cell attributed to Richet (1989). |
| perovskite_orthorhombic | 181 | CC0 COD 9005757 pure MgSiO3 geometry, nominal 0.88 Mg + 0.12 Fe on A sites. Explicit topology proxy, not an Fe-bearing specimen refinement. |
| phase_d | 54 | CC0 COD 9001848, reported partial occupancies. No invented Fe occupancy; the Fe-bearing EOS specimen and structural surrogate differ. |
| sno2_cubic_27gpa | 29 | Reported Pa-3 Sn 4a / O 8c u=0.346 numerical model; ordered approximation to a source reporting oxygen disorder. No separate source-publication grant established. |
| sno2_pa_3_at_48gpa | 29 | Same positional model at the existing reference cell. The structural reference does not extend the separate EOS's 16.09–28.85 GPa range. |

## What is licensed and what is assessed

The six unchanged CIF inputs retain COD's [CC0 data/database declaration](https://www.crystallography.net/cod/),
with original headers and citations intact. Checksums and source URLs are in
[the input manifest](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/calculated-pattern-inputs.json).
The ideal rocksalt model is independently authored. Contributor-held rights in
the generated numerical arrays and input transcription/arrangement are
dedicated to CC0; calculation code remains MIT. This does not license a source
publication, source prose or artwork, author software or an inherited card.

FeH2 and the two SnO2 inputs use a **separate numerical-content assessment**,
not a newly authenticated open source-publication grant. Only the small,
attributed numerical structural models are represented in an independently
authored format. The calculation does not use the 2018 FeH2 supplement, whose
repository metadata declares CC BY-NC 4.0. No publisher PDF/CIF, figure, source
layout or source intensity collection is added. Unspecified source terms remain
unspecified; output CC0 applies only to rights held by the contributors.
This follows the project's [numerical-content/collection distinction](eos-structure-compilation-rights.md#numerical-parameter-and-structural-content),
without asserting worldwide public-domain status, absence of contracts or
Debian approval.

The FeH2 parameters are directly stated in [Pépin et al. (2014)](https://doi.org/10.1103/PhysRevLett.113.265504),
Fig. 2 and pp. 265504-2/-4. Fe coordinates are experimentally refined;
hydrogen positions are from the original DFT model, with possible vacancies
remaining a limitation. No claim is made that it supersedes later H models.
The SnO2 citation is corrected to [Sneed et al. (2019)](https://doi.org/10.1107/S1600577519003904),
J. Synchrotron Rad. 26, 1245–1252, the Introduction. The previous Smith/page/DOI
combination was incorrect. The original numerical u=0.346 is unchanged; it is the introductory
reference value, not the separately refined u=0.3343(5) at 44.5 GPa.
The ice CIF's archival year field says 1984; the corresponding
[Besson paper](https://doi.org/10.1103/PhysRevB.49.12540) is from 1994.
The source file is preserved verbatim and the metadata discrepancy documented.

## Calculation and reproducibility

The repository-only maintenance tool is
[`scripts/calculate_reference_patterns.py`](https://github.com/CPrescher/peritheos/blob/main/scripts/calculate_reference_patterns.py).
Installed packages ship `peritheos/data/calculated-pattern-inputs.json`, the six
CIFs, all calculated arrays and the [scoped calculation notice](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/calculated-patterns.LICENSE.md).
Gemmi is needed only to regenerate patterns; it is not a runtime dependency.

```bash
uv run --with gemmi==0.7.4 python scripts/calculate_reference_patterns.py --check
```

The model expansion uses explicit CIF symmetry operations, preserving
majorite's nonstandard origin. Decimal special fractions within 5e-6 of a
twelfth are snapped to avoid splitting special-position orbits. The independent
reciprocal sphere covers 0.8 <= d <= 20 angstrom. Full-sphere |F|² weights
include multiplicity and Friedel partners; coincident d values rounded to
10 decimals are merged and normalized to maximum 100. Numerical-zero weights
are filtered as specified in the input manifest. Gemmi 0.7.4 supplies neutral-atom
IT92 X-ray form factors. No source peak indices, cutoffs or intensities are
used to select the new lines.

Lorentz-polarization, absorption, texture, anomalous dispersion, displacement
factors and instrumental broadening are omitted. In particular, these line
weights are not a simulation of an arbitrary beamline experiment or source
Rietveld profile. The calculation is structurally/numerically checked, not
experimentally validated.

## Preservation and remaining issues

The [replacement receipt](../data/calculated-pattern-replacement.json) records
before/after document and scientific hashes, per-array original/calculated
hashes, unchanged-content projections and the private preservation-manifest
hash. Original arrays are not embedded in a public replacement fixture.
All EOS records, observation declarations/resource bytes, lattice metrics,
uncertainty conventions and scientific classifications remain unchanged.
Only ten peak arrays, their root provenance/notes and the explicitly identified
Fe3S proxy sites change. The other 224 material files remain byte-identical.

Historical receipts and frozen source pins are preserved. Current scientific
hashes change; tests compare exact unchanged projections across the authorized
replacement instead of pretending that the old full hash still describes the
current catalog. Compatibility accepts only the recorded exact source pairs.

The ten inherited collections are no longer current distribution inputs.
**Seven Holland/THERMOCALC parameter extracts remain unresolved**, and the three
numerical structural models retain the source-term caveats above. The optional
observation package remains 451 material plus two Barker entries, 392 CSVs,
with the same five held original Sakai/Campbell entries. No published archive
or Git history is rewritten.

Repository checks include exact recalculation, independently derived rocksalt
intensity ratios, cell contents, analytical spacings, systematic absences,
normalization, unchanged scientific projections and installed-package exports:

```bash
uv run --with gemmi==0.7.4 python -m pytest -q -W error tests/test_calculated_reference_patterns.py
uv run --with gemmi==0.7.4 python -m pytest --rootdir=. -c .github/pytest-compatibility.ini -q -W error
uv run python -m scripts.review_eos_structure_rights --check
uv run python -m scripts.generate_material_summaries --check
```

Two additional historical checks remain stale: the older frozen-release check
fails first on the unchanged `fe_fcc` document, and the coverage snapshot
predates seven already-added dataset declarations across five other materials.
The coverage's whole-catalog checksum also changes with these new patterns.
Neither historical baseline is repinned here. The current scientific-preservation
receipt, exact ten source pairs and current package exports are checked
separately; these checks do not constitute experimental validation.
