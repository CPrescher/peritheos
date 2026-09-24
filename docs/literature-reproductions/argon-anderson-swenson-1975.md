# Anderson and Swenson (1975): argon source-access audit

M. S. Anderson and C. A. Swenson, *Experimental equations of state for the
rare gas solids*, Journal of Physics and Chemistry of Solids **36**(3),
145-162 (1975), DOI
[10.1016/0022-3697(75)90004-9](https://doi.org/10.1016/0022-3697(75)90004-9).

## Outcome

**Source access blocked; not reproduced.** The requested study is retained as
`supporting_studies/argon_anderson_swenson_1975` in `argon_fcc.eosmat`.
This is non-executable provenance metadata, preserved by `.eosmat` round trips.
There is no new EOS record, model implementation, or measured-point dataset.
The failure to recover the primary text does **not** establish that the paper
lacks an independently reproducible EOS.

The publisher's indexed abstract identifies piston-displacement measurements
to 20 kbar (2 GPa), from 4.2 K to each solid's triple point, with estimated
volume accuracy about 0.1%. It describes individual isotherms using relations
inspired by a Lennard-Jones potential. These abstract-level facts do not supply
the equation or coefficients. The abstract is indexed in the recommended
article section of this [publisher page](https://www.sciencedirect.com/science/article/abs/pii/0375960173905689);
that page's host article is a **different 1973 publication** and must not be
downloaded or imported as the requested paper.

## Scientific boundaries

- The fcc material association is provisional context, not a primary-source
  phase determination. Wittlinger's separate hcp entry is untouched.
- No source-defined reference pressure or temperature was recovered. In
  particular, neither 0 K nor 300 K is silently assigned to a missing curve.
- 20 kbar is converted exactly to 2 GPa. The full argon-specific sampling and
  fit ranges, lower pressure limit and temperature grid remain unknown.
- Original volume units and normalization are unverified. No cell, atomic or
  molar volume has been synthesized. The existing material's four-atom fcc
  cell must not be mistaken for this study's source volume convention.
- The approximate 0.1% overall volume accuracy is not a one-sigma point error
  or parameter uncertainty. No coefficient covariance or residual statistic
  is available.
- Pressure calibration and its temperature dependence are unresolved. A
  piston-displacement experiment is not evidence for any ruby pressure scale.
- No observations, theory curves or digitized markers were recovered from
  the original. Thus no numerical residual, fit validation or curve comparison
  is claimed. Later fits cannot be presented as the 1975 equation.

The coordinator supplied Xiao's later preprint
([Research Square v1](https://www.researchsquare.com/article/rs-5370077/v1.pdf)).
Its printed p. 22 (PDF p. 24), inspected visually, attributes bulk-modulus
data to this study. That is a useful recovery lead, not the original paper or
an authority for its published equation; no values were imported from it here.

## Retrieval attempts (2026-09-24)

| Route | Observed result |
|---|---|
| Zotero local API, title and `Swenson` searches | Healthy API; no matching item or attachment found. No writes performed. |
| Crossref DOI record | Verified title, author initials, journal, volume 36, issue 3 and pages 145-162. Links target Elsevier PII `0022369775900049`. |
| ScienceDirect article and abstract pages | Web fetch failed; browser stopped at a CAPTCHA. Challenge was not solved. |
| Elsevier Article Retrieval XML | HTTP 200 basic metadata only; no body, tables or abstract. `view=FULL` returned HTTP 401 Unauthorized. |
| Elsevier PDF representation | HTTP 406 Not Acceptable. No PDF bytes obtained. |
| Elsevier plain-text representation | HTTP 400 Bad Request. Subsequent metadata archival request was rate limited (429); no further retries. |
| OpenAlex W2062619387 | Closed; `oa_url=null`, `any_repository_has_fulltext=false`, sole publisher location has `pdf_url=null`. |
| Semantic Scholar DOI record | `openAccessPdf.status=CLOSED`, empty PDF URL. |
| Exact-title/DOI searches and targeted Ames/Iowa State/OSTI searches | References and later papers found, no original manuscript or supplement. OSTI exact-title API query returned no record. |

The stable local directory
`/Users/clemens/Documents/Peritheos Sources/argon/anderson-swenson-1975/`
contains Crossref, OpenAlex and Semantic Scholar **metadata only**, with
checksums in `retrieval-manifest.json`. None is a PDF or a substitute for one.
No copyrighted PDF is committed. The primary PDF path, origin and SHA256 are
all null in the handoff because no matching PDF was obtained or verified.

## Reproduction and integration

The existing argon tests check that the unchanged executable fcc curves,
source datasets and historical hcp diagnostic still work. JSON Schema and
Python material round-trip checks additionally verify that this supporting
study survives serialization without becoming an EOS record.

Studio should expose this as a literature/source-access status entry, with
no selectable fitted curve and no experimental-point toggle. Generic JSON
consumers can retain `supporting_studies`; dedicated Studio presentation is
the integration task's responsibility. Metadata import to Zotero C18 can
proceed, but **the user's requested PDF attachment remains unfulfilled**.
Revisit the original equation, tables, calibration and identity checks when
a legitimate matching full-text copy becomes available.
