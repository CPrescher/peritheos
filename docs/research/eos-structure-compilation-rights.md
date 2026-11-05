# EOS and structure compilation rights

Reviewed 5 November 2026 against `73bb757a`. This completes a component-level
rights/provenance inventory of the **685 EOS records and 234 structural/diffraction
components** in the current checkout. It does not establish universal rights
clearance. EOS values, observations and recommendation classifications are unchanged.
The later [powder-pattern replacement](calculated-reference-patterns.md) explicitly
changes ten diffraction arrays and the Fe3S structural proxy; historical receipts remain intact. A subsequent [COD replacement](post-perovskite-cod-replacement.md) changes only the two post-perovskite source/notes blocks and resolves their former ICSD acquisition cases. No EOS or observation selection changes.

| Assessed content | Components | Finding |
| --- | ---: | --- |
| Normalized EOS parameters | 678 | Scoped numerical-content assessment; no new source grant |
| Holland/THERMOCALC endmember extracts | 7 | Collection/acquisition review remains necessary |
| Inline structures and numerical fallbacks | 183 | Scoped numerical-content assessment; source/proxy qualifications retained |
| Structures tied to four reviewed COD resources | 5 | Recorded COD CC0 source scope, with numerical changes distinguished |
| Materials with no lattice, sites or peaks | 36 | No structural model values to assess; no missing values invented |
| Independently calculated reference patterns | 10 | Replaces inherited arrays; six COD CC0 inputs, one ideal model, three numerical-parameter assessments |
| Former ICSD-linked post-perovskite topology | 0 unresolved | Two inputs replaced with CC0 COD 9009217; included in the five COD components |

These are components, not independent experiments. The EOS references form
**272 source groups**; repeated fits, uncertainty variants and recalibrations
are considered together rather than counted as separate sources. Structural
and EOS provenance are assessed separately, including within the same material.

## Exact, reproducible scope

The installed
[review receipt](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/eos-structure-rights-review.json)
has an entry for every component, its source pointer, original-content and
runtime-export fingerprints, source references and assessment. It contains no
observations or copies of source publications. The source-bound
[policy](../data/eos-structure-rights-policy.json) requires all 234 material
files to match their reviewed contents exactly. The separate COD replacement receipt chains the two authorized source-only changes to their predecessors. New or changed records require a new review;
running the generator does not silently approve them.

The [component index](../data/license-inventory/eos-structure-review.csv) lists
each decision and remaining action. The
[EOS source-group index](../data/license-inventory/eos-source-groups.csv) groups
all 685 records by their actual citation. In the main inventory,
`redistribution_assessment`, `review_evidence` and `remaining_issue` are separate
from the existing source-declaration `status`. Unspecified terms do not become
open-license declarations merely because numerical redistribution was assessed.
Historical release rows receive no backdated assessment. Dedicated
`structure_reference` and `structure_provenance` fields are now also indexed;
an EOS citation alone is not a structural source.

## Numerical parameter and structural content

The numerical-content assessments concern the existing normalized scientific
coefficients, uncertainty/covariance conventions, reference states, lattice
metrics, symmetry, coordinates and occupancies. Mathematical model choices,
independent implementations and factual numerical results are distinguished
from source text, artwork, table layout, author code and database collections.
Ideal or proxy coordinates are not relabeled measured refinements. This review
is not another scientific-source validation or author-workflow reproduction.

The factual/collection distinction follows the
[Creative Commons database explanation](https://creativecommons.org/faq/#which-components-of-databases-are-protected-by-copyright),
[17 USC 102(b) and 103](https://www.copyright.gov/title17/92chap1.html), and
[CJEU C-203/02](https://curia.europa.eu/juris/showPdf.jsf?docid=64559&doclang=EN).
The latter distinguishes creating data from investment in obtaining existing
materials; it does not exempt all experimental collections. No numerical
row-count threshold, worldwide public-domain finding or absence of contracts
is asserted. Whole-source extraction, repeated extraction, original selection
and qualifying verification/presentation investment still need consideration.
The [IUCr's collection guidance](https://www.iucr.org/what-we-do/digital-standards/data/datasets)
also distinguishes attribution from collection-specific reuse conditions.

These are maintainer applicability assessments of the inspected catalog
content, not a fresh legal analysis of every publication, authentication of
every contributor's title or audit of historical acquisition agreements.
Contributor and dataset declarations remain limited to their original scopes.

## Dioptas migration: an actual notice to preserve

The material manifest names Dioptas 0.10.0 at commit
`5a8bfd81d10bfab3499039603380aae34576d60a`. Its exact Git objects were available
locally: 120 original material JSON files and the root `LICENSE.txt` were
inspected/preserved privately. There are 144 remaining EOS records carrying
that migration provenance, across 114 materials. This is not a claim that all
later records or structural refinements came from that import.

The pinned [upstream license](https://github.com/Dioptas/Dioptas/blob/5a8bfd81d10bfab3499039603380aae34576d60a/LICENSE.txt)
is MIT, with copyright lines for GSECARS/University of Chicago, University of
Cologne, DESY and Clemens Prescher. Its exact text is now retained in
`dioptas-materials.LICENSE.md`, including in both packages and their license
metadata. This covers the upstream contributors' own contributions, not an
assumption that they owned all underlying publications or inherited cards.

## Remaining cases requiring action

### Seven Holland/THERMOCALC extracts

The seven Holland et al. (2013), DOI
[10.1093/petrology/egt035](https://doi.org/10.1093/petrology/egt035), records are
the apv, mpv, cpv, fper, fpv, per and stv four-coefficient isotherms. Their exact
IDs appear in the component index. A fresh download of the
[official metabasite bundle](https://hpxeosandthermocalc.org/wp-content/uploads/2025/05/tc-thermoinput-metabasite-2022-01-30.zip)
recovered the same `tc-ds62.txt` as the scientific audit:
SHA-256 `46fa2ecb4fccdeaf44984da81342455d10e920012ab8d6db35b0bc94a7e0781c`.
It contains 256 endmembers. Neither that file nor the inspected bundle README
provides an explicit redistribution grant. The
[official site](https://hpxeosandthermocalc.org/) describes academic use, while
its [metabasite page](https://hpxeosandthermocalc.org/the-hpx-eos/the-hpx-eos-families/metabasite-set/)
supports implementations in other software. These statements do not establish
unrestricted redistribution of the collection extract; equally, they do not
prove that every published coefficient is protected.

**Action:** obtain file/extract-specific permission from the relevant rights
holders, or complete a collection/acquisition assessment establishing a basis
for the exact extract. The runtime contains converted numerical coefficients,
not the complete thermodynamic collection, axfiles or THERMOCALC executable.
The separate commercial product Thermo-Calc is unrelated; its terms do not
apply merely because the names resemble one another.

### Resolved by replacement: ten inherited peak lists

| Material | Historical rows (private) |
| --- | ---: |
| fe3s | 18 |
| feh2 | 10 |
| iceviii | 33 |
| majorite | 103 |
| mg7si2o8 | 35 |
| mgfe60o | 10 |
| perovskite_orthorhombic | 156 |
| phase_d | 65 |
| sno2_cubic_27gpa | 11 |
| sno2_pa_3_at_48gpa | 11 |

Before replacement, all **452** normalized hkl/d/intensity rows match the corresponding pinned
Dioptas arrays exactly, including the `mgsio3.json` origin of majorite. Upstream
file hashes and peak-array digests bind those comparisons. No primary
acquisition agreement or original card-maker grant was established. The
[ICDD's published policy](https://estore.icdd.com/pages/Customer-Service.html)
describes its Powder Diffraction File as licensed/copyrighted; a `.jcpds`
filename alone does not establish that these arrays came from that product.
Conversely, public Git hosting and a later MIT notice do not prove they did not.

The [calculation follow-up](calculated-reference-patterns.md) now replaces all
of these arrays with independent calculated patterns, using explicit input
models rather than the inherited reflection selection/intensities. Originals
remain privately preserved. This resolves the inherited collection as a current
distribution input; the three factual-parameter models retain source-term caveats.
Other JCPDS card-number citations without imported peak arrays are provenance
leads, not automatically additional copies of the ICDD database.

The two Nisr hydrous-silica fallbacks and phase_e_fe0072 contain a separate
total of 24 position-only markers: indexed hkl, recalculated d and uniform
display values. They are assessed separately from inherited intensity arrays;
no experimental intensity or complete atomic model is inferred.

### Resolved: two ICSD-linked structures

The [COD replacement follow-up](post-perovskite-cod-replacement.md) reacquired
both structural inputs from CC0 COD 9009217, AMCSD 0015437. Three freshly
retrieved CIF versions agree except for SVN Date/URL comments. The new source
regenerates every prior cell/site number exactly. The data identify ab initio
coordinates using the observed 121 GPa, 300 K cell; they are not experimentally
refined positions. The natural Fe-bearing material retains its explicitly
limited pure-MgSiO3 topology proxy. Original full records, previous rights
reviews and source evidence are private, checksummed preservation copies.

This resolves the current distribution input by independent open-source
acquisition, not by claiming a new license for the historical ICSD path.
Coesite-V still has no structural model and is not an additional ICSD input.

## Distribution and preservation

The ten reviewed, checksum-bound COD resources retain their
recorded source declaration. COD's [current mirror](https://cod.psdi.ac.uk/)
confirms a CC0 dedication of the data and database. The new review does not
extend it to other inline COD references without a verified resource binding.
Superhydrous phase B's separate erratum correction and all existing structural
normalizations remain explicit. The original He thesis structure remains
privately preserved and excluded; it is not reinstated.

The runtime and optional packages ship this review and the scoped compilation
and Dioptas notices. The optional package still selects 451 material plus two
Barker entries, using 392 CSV resources. Its original five holds are unchanged.
The seven Holland components needing further review above remain in the current runtime
export; the pattern replacement does not resolve their collection/acquisition scope. They need resolution or explicit downstream exclusion before claiming
an entirely cleared runtime archive. Neither the optional split nor this
assessment is Debian approval. Published 0.17.3/0.17.4 archives remain unchanged.

Original material files, observation resources, historical scientific reports,
the original inventory, pinned upstream files, and freshly retrieved policy/
THERMOCALC evidence are preserved. The later COD replacement changes only two material source/notes blocks and adds the attributed CIF; original model values and all scientific findings are preserved. Notices, review metadata, inventory outputs, documentation, tests and packaging notice inclusion record these boundaries.

Repository maintenance checks:

```bash
uv run python -m scripts.review_eos_structure_rights --check
uv run python scripts/inventory_distribution_licenses.py \
  --archives .git/pr27-domain-review/license-archives \
  --cargo-metadata .git/pr27-domain-review/eos-structure-rights/cargo-metadata.json --check
uv run pytest -q -W error tests/test_eos_structure_rights.py tests/test_license_inventory.py tests/test_distribution_split.py
```

The private cache paths are examples for this checkout, not installed tools or
distributed files. The inventory command requires verified release archives
and the release-matched Cargo metadata described in the inventory guide.
