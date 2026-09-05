# Zhang and Weidner (1999): aluminous silicate perovskite

## Primary source and audit boundary

This audit covers all eight LitCurate rows grouped under Jianzhong Zhang and
Donald J. Weidner, *Thermal Equation of State of Aluminum-Enriched Silicate
Perovskite*, **Science 284**, 782--784 (1999),
[doi:10.1126/science.284.5415.782](https://doi.org/10.1126/science.284.5415.782).

The publisher and PubMed abstracts establish that the study reports new
pressure-volume-temperature measurements for silicate perovskite containing
5 mol% Al2O3. The current Science landing page and PDF endpoint did not expose
the article body. The publisher's archived 2012--2015 abstract/full-text URLs
expose only the abstract, citation metadata, and a registration/login page.
OpenAlex classifies the article as closed and lists no repository full text,
and exact title/DOI searches found no author manuscript, official supplement,
or primary data deposit. No local paper copy was found.

Consequently this paper contributes **no production EOS record** in this
batch. The three purported run parameterizations cannot meet Peritheos's
primary-authority, exact-equation-mapping, and numerical-reproduction
requirements. The other five rows are values copied from earlier MgSiO3
perovskite publications and remain citation traces.

## Material identity correction

LitCurate labels all eight rows `MgSiO3`, but that is not the identity of the
three source-result rows. The primary abstract explicitly identifies the
sample as MgSiO3 perovskite containing 5 mol% Al2O3. A 95:5 molar mixture of
MgSiO3 and Al2O3 contains Mg95 Al10 Si95 O300, which normalizes to
`Mg0.95Al0.10Si0.95O3`. This is the charge-balanced coupled-substitution
formula already represented by the bundled
`mg095al010si095o3_bridgmanite` material. The three rows therefore must not be
attached to pure `bridgmanite`.

The nominal formula is a composition normalization, not evidence that Zhang
and Weidner refined a unique Al site distribution. A future implementation
must preserve the paper's own sample characterization and substitution
language once the body is available.

## Purported run-specific fits

The discovery ledger attributes the following values to three experimental
runs:

| run | LitCurate `V0` (A^3, conventional cell) | `K0` (GPa) | `K0'` |
|---|---:|---:|---:|
| 1 | 163.051 | 234 | 4 fixed |
| 2 | 163.184 | 232 | 4 fixed |
| 3 | 163.220 | 236 | 4 fixed |

Later primary papers summarize Zhang and Weidner's overall room-pressure
result as `K0T=234+/-2 GPa` with `K0'=4`, which is numerically consistent with
the center and spread of the three ledger values. That secondary consistency
does not establish that Science published three independent standalone BM3
parameterizations. Without the body or data it is unresolved whether the
numbers are separate fits, run normalizations feeding one joint thermal fit,
or intermediate results whose covariance is shared.

The inaccessible source is also needed to resolve:

- the exact thermal EOS expression and whether each listed `K0` is a
  separately executable 300 K isotherm;
- which P-V-T rows belong to each run, their ranges, uncertainties, exclusions,
  and weighting;
- the pressure standard and temperature measurement/calibration;
- coefficient uncertainties and covariance beyond the later summary; and
- an independent off-reference source value or curve for numerical
  reproduction.

A standard BM3 calculation from the ledger coefficients would only reproduce
the ledger itself. It would not independently reproduce the primary paper and
therefore was not used to promote these candidates.

## Disposition of every same-DOI LitCurate row

| LitCurate identifier | reported row | disposition | reason |
|---|---|---|---|
| `litcurate_0de9bb03bd0be90b` | Al-Pv run 1: `V0=163.051 A^3`, `K0=234 GPa`, `K0'=4` | **held / primary data unavailable** | The primary body, equation definition, run data, range, calibration, and coefficient uncertainties could not be recovered; independent reproduction is impossible. |
| `litcurate_061ebc9e4a0a8335` | Al-Pv run 2: `V0=163.184 A^3`, `K0=232 GPa`, `K0'=4` | **held / primary data unavailable** | Same evidence gap; it is also unresolved whether this is an independent EOS or one run normalization in a joint thermal analysis. |
| `litcurate_e6b35bc95210b026` | Al-Pv run 3: `V0=163.220 A^3`, `K0=236 GPa`, `K0'=4` | **held / primary data unavailable** | Same evidence gap; no source P-V-T table or digitizable primary plot was recovered. |
| `litcurate_8ffed77f7af8d194` | MgSiO3 Pv, reference 13: `K0=261 GPa`, `K0'=4` | **citation trace only** | Earlier publication, missing `V0`; source lineage must remain with the cited primary paper. |
| `litcurate_ed509c83d0db4e9b` | MgSiO3 Pv, reference 14: `K0=261 GPa`, `K0'=4` | **citation trace only** | Earlier publication, missing `V0`; not a Zhang-Weidner source result. |
| `litcurate_e2ddb60982ab5753` | MgSiO3 Pv, reference 16: `K0=261 GPa`, `K0'=4` | **citation trace only** | Earlier publication, missing `V0`; not a complete executable EOS. |
| `litcurate_5c1f8dc63c75c909` | MgSiO3 Pv, reference 17: `K0=261 GPa`, `K0'=4` | **citation trace only** | Earlier publication, missing `V0`; equation identity is unresolved. |
| `litcurate_f7f00f45db766f48` | MgSiO3 Pv, reference 18: `K0=261 GPa`, `K0'=4` | **citation trace only** | Earlier publication, missing `V0`; audit under the underlying publication. |

## Duplicate and future-work check

No bundled record uses DOI `10.1126/science.284.5415.782`. The existing
`mg095al010si095o3_bridgmanite` material is the correct nominal-composition
target and already contains independent EOS records from later studies; none
duplicates these three held coefficient sets.

Reconsider the three source candidates only if the final Science body, an
author-provided copy, or an official data file becomes available. The future
audit should first determine whether there is one joint thermal EOS or three
scientifically independent run records, then transcribe or digitize every
recoverable source observation and reproduce at least one off-reference
published value before adding an executable record.
