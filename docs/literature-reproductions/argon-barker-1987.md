# Barker (1987): theoretical argon isotherm, full-text access blocked

J. A. Barker, *High pressure equation of state for solid argon from interatomic
potentials*, **The Journal of Chemical Physics 86**(3), 1509–1511 (1 February
1987), [doi:10.1063/1.452187](https://doi.org/10.1063/1.452187).

The packaged [study record](../../peritheos/data/studies/argon-barker-1987.json)
is explicitly **theoretical**, **not_reproduced**, and **fulltext_access_blocked**.
It adds a traceable supporting-study result without creating a selectable EOS.
It is not a statement that the original model lacks an independent EOS: that
cannot be determined without its full text. No numerical reproduction or fit
validation has been claimed.

## Verified primary evidence

The [publisher record](https://pubs.aip.org/jcp/article/86/3/1509/660425/High-pressure-equation-of-state-for-solid-argon)
and [IBM Research's author-affiliation record](https://research.ibm.com/publications/high-pressure-equation-of-state-for-solid-argon-from-interatomic-potentials)
identify a 298 K Monte Carlo calculation with the Barker–Fisher–Watts (BFW)
pair potential plus the Axilrod–Teller (AT) three-body interaction. The abstract
compares it with Ross and coworkers' experimental compression to 800 kbar
(80 GPa), and with Aziz–Chen, Koide–Meath–Allnatt and Ross exp-6 potentials.
Those are separate alternatives; the exp-6 potential must not silently replace
the BFW+AT model.

The publisher's visible reference list identifies BFW as Barker, Fisher and
Watts (1971), [doi:10.1080/00268977100101821](https://doi.org/10.1080/00268977100101821),
and the experimental comparison as Ross, Mao, Bell and Xu (1986),
[doi:10.1063/1.451346](https://doi.org/10.1063/1.451346).
The bibliography includes Mao pressure-scale references, but merely citing them
does not verify which scale was applied to which comparison points.

Crossref and the publisher agree on the February 1 publication date. IBM's
January 1 date is generic yearly metadata and is not used as the article date.

## Scientific boundaries

- **Phase:** solid argon verified. The precise simulated crystal structure has
  not been verified from the full paper, so the record does not assert fcc or hcp.
- **State:** 298 K is the simulation isotherm, not an inferred zero-pressure
  reference temperature for a parameterized cold curve.
- **Range:** 800 kbar is the abstract's experimental comparison upper limit.
  Simulation volumes, lower bound, fit range and extrapolation limits are unknown.
- **Units and normalization:** kbar converts to GPa by division by ten. The source
  volume unit and whether volume is molar, atomic or cell-based remain unknown.
  No multiplication by four is applied speculatively.
- **Uncertainties and scale:** potential uncertainties, numerical sampling errors,
  comparison pressure scale and any corrections remain unresolved.
- **Data:** no printed simulation table or figure has been recovered. There are
  no bundled Barker observations or synthetic model points. Ross measurements
  belong to the Ross study and must retain that provenance.
- **Implementation:** no callable Python or Rust EOS is registered. Implementing
  BFW+AT requires the exact source potential conventions, coefficients, cutoffs,
  simulation protocol, finite-size/tail corrections and any quantum correction.
  A static lattice sum or convenient analytic EOS would not establish reproduction
  of this finite-temperature Monte Carlo calculation.

## Full-text search and precise blocker (24 September 2026)

The Crossref-deposited [publisher PDF URL](https://pubs.aip.org/aip/jcp/article-pdf/86/3/1509/18963166/1509_1_online.pdf)
returned HTTP 403 to a direct request. Browser navigation reached the article
successfully, then displayed that current access was unavailable. Following the
page's PDF access link returned to the abstract with `redirectedFrom=PDF`.
The publisher says the full content is available only as a PDF. No purchase,
credential changes or author contact were attempted.

IBM's primary publication page contains the abstract and a DOI link, with no
manuscript or supplementary download. Exact-title and DOI searches, including
institutional/repository searches, found citations and later papers rather than
this full text. Those PDFs are not substitutes for the requested Barker PDF.
OpenAlex (`W2001562305`) reports closed access, one publisher location and no
repository full text. Semantic Scholar (`780aba2e3893a66efd59b63f6f2f18420d02b4d4`)
reports CLOSED and an empty PDF URL. Read-only Zotero searches for `Barker` and
`1.452187` both returned no items.

**No matching full-text PDF was obtained.** Its path and SHA256 are null, not
hashes of an abstract or HTML access page. The requested PDF import remains
unfulfilled and must be reported as such. Source captures and checksums are
available outside git in
`/Users/clemens/Documents/Peritheos/argon-sources/barker-1987/`.
The integration handoff identifies them separately from a full-text PDF.

## Integration behavior

Display this as a supporting theoretical study with an access/reproduction
limitation. Do not show a pressure curve, fitted parameter table, experimental
points, successful reproduction badge or inferred reference state. The packaged
`data/studies/` JSON is an explicit handoff artifact; the existing material/EOS
loader does not automatically index this new supporting-study directory.
The integrator must wire it into Studio's study/outcome presentation if desired.
Existing fcc EOS records and Wittlinger's hcp/not_reproduced result are untouched.

## ResearchGate follow-up

The exact [Barker publication page](https://www.researchgate.net/publication/254001099_High_pressure_equation_of_state_for_solid_argon_from_interatomic_potentials)
was additionally checked on 24 September 2026 in the existing signed-in browser
session. It matches the title and DOI but offers only **Request full-text**.
The public page explicitly says no full text is available. An institutional
publishing-fee coverage notice does not establish reading access; no download
button was available for this paper. No author request was sent. The PDF and
reproduction statuses therefore remain unchanged.
