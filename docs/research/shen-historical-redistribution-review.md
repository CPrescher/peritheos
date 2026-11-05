# Shen/Smith and historical numerical collection review

Reviewed 5 November 2026 against the preserved sources and checkout at
`db93e2cf`. At this stage the unreleased `peritheos[datasets]` selection admitted **427 of
453 material dataset entries**. This adds **12 entries backed by three unchanged
CSV files**: ten views of the shared Shen/Smith results, the 25-cell copper
handbook excerpt and the 20-row Anderson gold table. These are source-specific
maintainer assessments, not new licenses, legal opinions or Debian acceptance.

The subsequent [Sakai/Ismailova/iridium review](restricted-source-numerical-review.md)
raises the current selection to **448 of 453 entries, with five held**, retaining
all original NC/ND source flags. The 427-entry count above is historical.

All original declarations, numerical resources, published parameters, source
rights assertions and scientific classifications are unchanged. The complete
[machine-readable receipt](../data/shen-historical-redistribution-review.json)
is also embedded as `shen_historical_review` in the installed redistribution
ledger. It binds each declaration, original resource and inspected source copy.

## Shen and Smith: complete experimental numerical collection

The final [Shen and Smith (2026) article](https://journals.aps.org/prb/abstract/10.1103/fxgq-96sg),
*Physical Review B* 113, 144113, describes two purpose-run simultaneous XRD
experiments. Sections II.B–II.C explain the creation of diffraction data,
peak fitting, unit-cell volumes, errors and pressure-drift checks. Section III
and the data-availability statement identify Supplemental Table S1 as the
experimental volume record. The final article was reread, and its methods page
was visually checked. Its SHA-256 matches the earlier curation copy:
`fb101fe78237c91077586d8999103a10c7ee75813a13ada441ee81a880c83f72`.

The admitted CSV has six functional columns: experiment, run, phase, sequence,
cell volume and volume standard error. There is no source artwork, article
prose, original workbook layout or author software. It retains all numerical
pairs, including lower-quality/rejected-fit subsets, rather than rewriting the
original evidence to match the preferred fits. The ten entries share one file;
they do not represent ten independent copies or additional experiments.

| Preserved profile | Result |
| --- | --- |
| Numerical volume/error pairs | 3511 |
| DAC-1 | 1595 pairs, contiguous numbered runs 1–228 |
| DAC-2 | 1916 pairs, contiguous numbered runs 1–194 |
| Cu reference readings | 422, with original run identities |
| Pt sequence labels | 421 first, 422 last; first Pt absent at DAC-2 run 194 |
| Other sequence cells | 2668 blank; no fabricated sequence labels |
| Phase labels | Cu, Pt, Au, Ta, Mo, W, MgO, bcc/hcp Fe and B1/B2 NaCl remain distinct |

The paper states 229 DAC-1 compression steps, whereas the preserved normalized
workbook transcription has 228 numbered runs. That difference is recorded;
no extra run is invented. No pressure column is added to this measured-volume
resource. Existing paired-Cu reconstruction and its scientific limits remain
separate. In particular, the MgO/NaCl-B2 diagnostic holds and all published
coefficients remain unchanged; distribution admission does not validate an EOS.

The prior exact normalization check used the official workbook with SHA-256
`339eb0d160316ef3e32abe67c5d5e917dfc6be4b4fbf726287680cca61596565`.
This is the [preserved earlier verification](../literature-reproductions/shen-smith-2026-cu-refits.md),
not a fresh workbook comparison. The current supplement URL returns 403,
the harvest endpoint returns 401, and the official landing page identifies the
supplement as subscription-required. The final article prints APS copyright.
**No separate open dataset grant was located.** Source access is not permission.

The assessment distinguishes investment in creating measurements and their
in-measurement uncertainty/drift checks from independent database investment.
The inspected source describes the former; it does not establish a separate
substantial investment in obtaining, verifying or presenting an independently
assembled collection. The flat scientific fields and functional run order are
assessed separately from protected presentation. This is a reasoned basis to
admit the numerical record, **not proof that qualifying independent investment
or acquisition restrictions are absent**. The extraction covers the complete
collection and is not presented as insubstantial. No general exemption for
experimental data, row-count threshold, author nationality or public-domain
status is inferred from laboratory affiliation or funding.

## AIP Handbook: limited copper numerical excerpt

This is the fixed 1972 third edition, Table 4d-12, Cu column, p. 4-100,
attributed to Kennedy and Keeler. The 25 populated pressure/relative-volume
cells are normalized into scientific columns; no table image, other material
columns, caption or handbook prose is included. The existing all-rights-reserved
source statement remains verbatim. The original book's copyright is not
asserted to have expired.

The [prior Cu-column visual check](../literature-reproductions/aip-handbook-cu-vinet-refit.md)
and preserved handbook-method record identify these as calculated 25 C
(298.15 K) isotherms reduced from shock data. The method cites van Thiel's
*Compendium of Shock Wave Data*, UCRL-50108 (1966). Exact selected shock
observations and reduction inputs remain unavailable. No US-government or open
license is inferred for that report. The 15 kbar token `09900` remains explicitly
interpreted as `0.9900`; blank 5/10 kbar cells and a zero-pressure identity are
not invented as observations. Phase remains unresolved and no refit changes.

The historic book SHA-256 is
`dda86a046830f63d95933b01e5b56bce9a46572a038ba37e68ceee9f2e059dad`.
The original full book is no longer at its recorded location, so this review
does not claim a fresh page/binary verification. The prior verification record,
exact CSV and its checksum are preserved. The admission concerns a limited
functional numerical excerpt, with the historical collection-term consideration
below; it does not extend to the whole handbook or other tables.

## Anderson gold: heterogeneous scientific input, not new observations

[Anderson, Isaak and Yamamoto (1989)](https://doi.org/10.1063/1.342969),
Table I, p. 1535, contains 20 rows of heat capacity, expansivity and elastic
modes. The table/methods and references 10–16 were reread; Table I was visually
checked. The current Zotero PDF has SHA-256
`0d2e21ec06c895c2dd601a0d6004c2207b1edc5040208f1de6c3ad009805abb8`.
It differs from the historical curation PDF, SHA-256
`38a4a6276270c3777d20c0437e22103612c32f5c0b622ebba32f3532499f9024`;
a cover/download timestamp is present. These are not claimed byte-identical.

| Numerical field/regime | Attributed source and Anderson treatment |
| --- | --- |
| Low-temperature Cp, 14 rows | Hultgren et al. (1963), with interpolation to some target temperatures |
| High-temperature Cp, five rows | Barin/Knacke (1973), mathematical expression in Eq. (4) |
| Low-temperature expansivity, 11 rows | White/Collins (1972), combined with AIP/Touloukian data and hand smoothing |
| Higher-temperature expansivity, eight rows | AIP (1972) and Touloukian et al. (1975), quadratic evaluation from Eq. (5) |
| Low-temperature elastic modes, 14 rows | Neighbours/Alers (1958) |
| Higher-temperature elastic modes, six rows | Chang/Himmel (1966), including a second 300 K branch |

The receipt provides the seven full source identities. Identifying them does
not claim a fresh audit of every upstream observation or reproduce the hand
smoothing. Standalone values, scientific method outputs and functional temperature
ordering are assessed separately from expressive table design. Only normalized
fields and ordinary source identifiers are distributed; no figure or prose is added.

The two 300 K elastic rows remain distinct, including the blank Cp/expansivity
fields for the second branch. The literal printed Cp `1.1288` remains beside
the separately labeled reconstructed `1.2888`. Missing uncertainties do not become
zero. The [existing staged EOS audit](../literature-reproductions/anderson-1989-gold-thermal-eos.md)
and `not_refittable` classification are unchanged: this is heterogeneous input,
not a common P-V-T regression matrix or a new global fit.

## Collection rights and limits

[UrhG §87a](https://www.gesetze-im-internet.de/urhg/__87a.html) addresses
qualifying investment and substantially changed databases;
[§87d](https://www.gesetze-im-internet.de/urhg/__87d.html) gives database-maker
rights a fifteen-year term. [§137g](https://www.gesetze-im-internet.de/urhg/__137g.html)
sets the German transitional window at 1983–1997 with a term beginning in 1998.
The fixed 1972 edition precedes that window, and a transitional term for the
fixed 1989 collection has elapsed by 2026. The named upstream collections are
older. No later substantially reinvested collection is extracted in this review.
That supports the historical numerical admissions; **publication copyright is a
separate question and is not declared expired or globally absent**.

[UrhG §2](https://www.gesetze-im-internet.de/urhg/__2.html) and
[17 USC §103](https://www.copyright.gov/title17/92chap1.html#103) distinguish
protected authorship/compilation contributions from their underlying material.
The [CC database FAQ](https://creativecommons.org/faq/#which-components-of-databases-are-protected-by-copyright)
also separates facts, sheet design, original selection/arrangement and database
rights. The [CJEU C-203/02 judgment](https://eur-lex.europa.eu/legal-content/EN/ALL/?uri=celex:62002CJ0203),
paras. 31–36, distinguishes data creation and independently qualifying collection
investment; **only official indexed excerpts were accessible**, not fresh full
judgment text. An empty/challenge response is recorded as unusable evidence.

There is no blanket declaration that all scientific tables are unprotected,
all experimental collections are exempt, or small excerpts are automatically
safe. The decisions concern these exact resources, with scientific/functionally
constrained fields and no copied original presentation. Acquisition contracts
have not been exhaustively audited. Older database-maker terms, later qualifying
investment, continuing publication copyright and other jurisdictions remain
separate considerations. No new CC0, CC BY or author permission is created.
Downstream distributors may choose narrower exclusions; this is not Debian
acceptance or legal advice.

## Distribution and preservation

The installed scoped notice is `notices/shen-historical-collections.LICENSE.md`.
Every admitted declaration/resource is bound by the package policy and review.
The original CSVs and relevant material cards are preserved privately as exact
baseline copies, in addition to unchanged tracked originals and Git history.
The inspected PDFs are private evidence; no book, article, image or workbook is
added to either package. The complete scientific snapshot remains unchanged.

**Still held:** 22 entries associated with the known NC/ND sources and four
entries pointing to the original Campbell integration tables. The separately
calculated Campbell polynomial approximations remain available. The two Barker
supporting-study entries remain deferred for interface/source review.

Relevant checks include the source-preservation and distribution tests, the
unchanged AIP/Anderson scientific replays, package builds, all selected installed
dataset loads, held-resource exclusion and installed notice/review verification.
