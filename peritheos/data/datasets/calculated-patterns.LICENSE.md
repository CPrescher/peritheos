# Independently calculated reference powder patterns

The ten arrays named in `../calculated-pattern-inputs.json` were independently
calculated by Peritheos contributors on 5 November 2026. To the extent they
hold rights in the generated arrays, input transcription/arrangement and the
independently authored ideal model, the contributors dedicate those rights to
CC0 1.0 (full text: `LICENSES/CC0-1.0.txt` in the distribution). This is not a
license for a third-party publication, database, source artwork or software.
Contributor calculation code remains MIT.

Six unchanged COD CIFs have the COD CC0 data/database declaration:

| COD | Material/model | Source |
| --- | --- | --- |
| 9007338 | Fe3S proxy: all Fe/Ni -> Fe, P -> S | Skala and Cisarova, meteoritic schreibersite, 2005 |
| 1566661 | Deuterated ice VIII | Besson et al., 1994 |
| 9010029 | Majorite, including Mg/Si disorder | Angel et al., 1989 |
| 9006822 | Phase A, including H sites | Kagi et al., 2000 |
| 9005757 | MgSiO3 topology; nominal 0.88 Mg / 0.12 Fe substitution | Yagi et al., 1978 |
| 9001848 | Phase D structural surrogate; no Fe occupancy invented | Yang, Prewitt and Frost, 1997 |

The ice CIF retains an apparent source-metadata error: its journal year is 1984,
while Physical Review B volume 49, page 12540 is the 1994 Besson paper
(DOI 10.1103/PhysRevB.49.12540). The archival CIF is left unchanged.

Their exact source references, original headers and symmetry operations remain
in `pattern-source-cod-*.cif`; checksums are in the input manifest and
`cif-license-review.json`. COD declaration: https://www.crystallography.net/cod/
The original research publications are separate from the COD data dedication.
Reference-cell overrides, substitutions and structural proxies are explicit.
The Mg0.4Fe0.6O model is an independently authored ideal rocksalt model; its
reference lattice remains attributed to Richet et al. (1989).

Three models use a separate numerical-content assessment, **not a newly found
open source-publication license**:

* FeH2: Pépin et al. (2014), DOI 10.1103/PhysRevLett.113.265504, Fig. 2 and
  pp. 265504-2/-4: Fe 4e z=0.853 and DFT H 4c/4d. The retained reference
  cell differs from the 67 GPa experimental cell. Nominal H occupancy is a
  model assumption, not an X-ray determination of exact stoichiometry.
* Two Pa-3 SnO2 reference cells: Sneed et al. (2019), DOI
  10.1107/S1600577519003904, the Introduction: Sn 4a and O 8c u=0.346. The
  source reports oxygen disorder, which these ordered models do not represent.
  Here u=0.346 is the introductory reference value, not the separately refined
  u=0.3343(5) at 44.5 GPa.

Only these attributed numerical structural parameters are represented in an
independently authored input format; no publisher CIF, prose, figure, table
layout, measured peak intensities or NC supplement is distributed here.
Unspecified publication terms stay unspecified. This assessment distinguishes
factual numerical inputs and independently computed results from protected
expression, collection selection and acquisition contracts; it is not a
worldwide legal clearance or a downstream-distribution approval. No source
license is inferred from a calculation or from the output CC0 dedication.

## Calculation and limits

The arrays are kinematic **X-ray** calculations using Gemmi 0.7.4 neutral-atom
IT92 form factors; D uses the H electron form factor. Each model is expanded
using the source CIF's explicit symmetry operations (including majorite's
nonstandard origin) or the named space group for the minimal inline models.
Decimal special-position fractions within 5e-6 of a twelfth are snapped before
expansion. All nonzero reciprocal vectors with 0.8 <= d <= 20 angstrom are
considered; coincident spacings rounded to 10 decimal places are merged.
Integrated line weights sum |F|^2 over the full sphere, including multiplicity
and Friedel pairs, and are normalized to a maximum of 100. Unmerged weights
at or below 1e-12 of the maximum are discarded as numerical zero; reported
relative intensities are rounded to 8 decimal places and rounded zeros omitted.

No Lorentz-polarization, absorption, texture, anomalous scattering, thermal
motion or instrumental broadening is modeled. These are reference sticks,
not observed peak strengths, a Rietveld model or experimental validation.
The inherited arrays/full original records are preserved privately; their
checksums, replacement scope and scientific preservation projections are in
`docs/data/calculated-pattern-replacement.json` in the source repository.
Historical released archives and Git history remain unchanged.

Repository-only reproduction:

```
uv run --with gemmi==0.7.4 python scripts/calculate_reference_patterns.py --check
```

Gemmi is an execution-only maintenance dependency under MPL-2.0; its code and
form-factor coefficient tables are not copied into these packages. See
https://github.com/project-gemmi/gemmi/tree/v0.7.4 and
https://gemmi.readthedocs.io/en/latest/scattering.html . Independently generated
outputs do not imply an open license for every input publication or the EOS.
