# Independent Sakai lineage and Campbell pressure integration

Reviewed 5 November 2026 against `cf30ca1b`. This follows the
[Sakai/Ismailova/iridium numerical review](restricted-source-numerical-review.md).
The purpose is general redistribution of useful, independently generated data.
The original source records and their numerical values remain unchanged.

**A separate Sakai evidence graph can be included. Full Campbell absolute tables
cannot yet be regenerated from the recovered inputs.** The new Campbell output
is explicitly limited to the independently calculated NNO pressure contribution.
It is not relabeled as the missing absolute table.

| New resource | Content | Optional material entries |
| --- | --- | ---: |
| `sakai-2025-independent-calibration-lineage.csv` | 37 calibration-evidence edges from 23 primary resources | 1, under gold |
| `campbell-2009-independent-nno-pressure-integral.csv` | 1717 calculated relative pressure contributions | 2, under nickel and nickel oxide |

The optional selection becomes **451 of 456 material entries**, with **390 CSV
resources**. All five original holds remain: the S8 compilation and four views
of the shared Campbell integration collection. This adds separate reconstructions
rather than silently substituting them under the original identifiers. Core
`peritheos` remains independently installable without observations.

## Sakai: independently assembled evidence, not a copied index

The generator reads only its explicit, checksummed primary-resource inputs.
It does **not** read `sakai-2025-table-s8-calibration-lineage.csv`, the raw
workbook archive containing S8, or old reconstruction reports. Tests reject
any unlisted input. The generator's sources and independent selection rules
are embedded in the installed declaration, with each edge carrying its DOI,
table locator, resource hash, observation count and selected-row digest.

This reuses prior primary-table transcriptions and their source provenance;
it is not represented as a fresh independent PDF transcription of every study.
The graph is broader than the author-selected S8 routes in some respects and
incomplete in others. It has a new functional ordering and contributor-written
annotations. Having inspected S8 previously is acknowledged: input independence
is not a legal clean-room claim.

Three relations remain distinct:

- **Paired volumes:** the seven present-study tables contain 393 paired states.
  The overlapping 326-row Au ancestry alias is not counted again.
- **Reported gauge/pressure:** original ruby coordinates, pressure-gauge labels
  and Fei's marker identities are retained. A ruby pressure is not silently
  relabeled as an IPPS-Ruby2020 reduction.
- **Co-measurement:** simultaneous primary table entries establish a shared
  observation, not necessarily a directed calibration or author fit selection.
  Dorfman pairs must share both run and observation index.

For example, Fei (2016) Table S1 has **49 Ne, 18 MgO and 7 NaCl-B2** marker
rows. These three groups are preserved, rather than assigning the whole table
to NaCl. Sakai (2014)'s pure-iron rows supply a NaCl partner but no populated
MgO-volume pairs; no missing Fe/MgO edge is fabricated. The graph also preserves
Anzellini's helium-gauge context and literature co-measurements outside the
selected S8 routes. Direction in a co-measured pair is for display only.

The exact S8 author selection, excluded masks, fitting sequence, unpublished
Ishimatsu inputs and complete recalibration choices are **not reconstructed**.
Links lacking direct numerical/primary-field support are omitted, not inferred
from the held index. The new graph is a useful partial lineage resource, not
a replacement claim for the original 35-row compilation or the author workflow.

## Campbell: independently integrate the recoverable pressure term

The final [Campbell paper](https://doi.org/10.1016/j.epsl.2009.07.022),
p.560, separates the absolute standard-state fugacity from a pressure integral:

```text
delta_log10_fO2(P,T) = 0.8686 * 1000 / (R*T)
                     * integral_0^P [V_NiO(p,T) - V_Ni(p,T)] dp
```

Volumes are cm³/mol and pressure is GPa, giving 1000 J/mol per product unit.
The printed 0.8686 factor is retained; R is 8.31446261815324 J/(mol K).
The source BM3 plus variable-exponent Debye formulation and printed Table 1
parameters are used at 295 K reference. Eight-point Gauss–Legendre quadrature
is applied to each 1 GPa interval. All coefficients and numerical conventions
ship inside each declaration. No new EOS parameters are fitted.

The generator reads only that explicit parameter configuration. The held
original output collection, source polynomial approximation and its outputs
are not generation inputs. It calculates pressure contributions for
1000–2600 K in 100 K steps and 0–100 GPa in 1 GPa steps. Every zero-pressure
contribution is zero by definition. This mathematical lower limit is explicitly
0 GPa; it is not represented as an independently recovered 1-bar Gibbs baseline.
The grid is a calculation domain, not a phase-stability or experimental-coverage
claim. The output remains diagnostic, not a recommendation for predictions.

Comparison is separate: the original Table S5 zero-pressure value is subtracted
from each original isotherm, then compared with the newly integrated increments.
Those baseline values are **never added to the generated output**.

| Comparison over 1717 coordinates | Difference in delta log10 fO2 |
| --- | ---: |
| RMS | 0.00528258402 |
| Maximum absolute | 0.01425912762 |

This is not exact reproduction of the original numerical pressure contribution.
The discrepancy remains unexplained; no alternative thermal law, parameter fit
or offset is used to make it disappear. Numerical quadrature convergence is
checked separately from source agreement. Existing published coefficients,
scientific dispositions and the original table bytes are retained.

## Why full absolute regeneration is still unresolved

Campbell specifies **Chase (1998)** for stoichiometric FeO and **Knacke et al.
(1991)** for NiO, uses only subsolidus Gibbs data, and extrapolates to higher
T without publishing the exact numerical implementation. IW also switches
from fcc Fe to the adopted Dewaele hcp model at the Shen (1998) boundary.

The [official JANAF crystalline FeO table](https://janaf.nist.gov/tables/Fe-018.txt)
was located and its cached numerical text inspected. It distinguishes the
crystalline data and melting annotation; liquid values are not a substitute for
Campbell's subsolidus extrapolation. A fresh direct text download returned 403;
no unavailable binary hash or exact Campbell extrapolation is claimed.

The [Shen primary phase-diagram paper](https://duffy.princeton.edu/sites/g/files/toruqf616/files/shen_etal_grl_1998.pdf),
DOI 10.1029/97GL03776, was inspected through accessible primary PDF text.
It does not recover Campbell's exact switching implementation. A fresh direct
PDF fetch returned 403. A fitted switch inferred from held Campbell outputs
would not be an independent reconstruction, so none is introduced.

No authenticated Knacke (1991) NiO entry was found in the checked Zotero and
Dropbox filenames or publicly accessible primary sources. The maintainer also
confirmed that no copy is available to supply. A different published
NiO compilation was located but not substituted. The exact Knacke entry and
source extrapolation/integration implementation remain the route to full
absolute reconstruction. An author calculation file would also resolve these
choices. A missing input is not replaced by an original table's baseline,
a polynomial fit, a modern thermodynamic database or an arbitrary Gibbs law.
The already included [polynomial approximation](morard-campbell-redistribution-followup.md)
remains separately identified; it does not resolve these gaps.

## Rights, preservation and verification

Contributor-held rights in the new graph selection, annotations, calculations
and arrangement are dedicated CC0-1.0; scripts are MIT. Primary numerical
resources retain their own source notices and scoped applicability assessments.
No original source workbook, article text, figure, author code, S8 commentary
or original Campbell collection is admitted by this work. The original source
works are not relicensed. The independent-output assessment concerns the actual
new content; it does not guarantee absence of every third-party right or imply
that a reformatted source compilation would be unrestricted.

The [bound receipt](../data/independent-lineage-buffers.json) is embedded as
`independent_lineage_buffers` in the installed review. Removing exactly the
three appended declarations restores the previous complete scientific snapshot.
Removing their serialized insertions restores every original material byte.
All original dataset arrays, coefficients, scientific classifications, held
CSV files and prior immutable review receipts remain intact. Original source
PDF/workbook evidence and pre-change copies are preserved privately.

Repository maintenance commands:

```bash
uv run python scripts/generate_sakai_2025_lineage.py --check
uv run python scripts/generate_campbell_2009_pressure_integral.py --check
uv run pytest -q -W error tests/test_independent_lineage_buffers.py tests/test_domain_review_preservation.py
```

Tests cover input independence, exact source-byte preservation, gauge identities,
missing evidence, closed-form quadrature/units, numerical convergence, source
misfit retention and installed exclusion of both held original resources.
