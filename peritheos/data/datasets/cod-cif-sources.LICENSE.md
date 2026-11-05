# COD crystallographic CIF source terms

Reviewed 4 November 2026; brucite and post-perovskite rechecked 5 November 2026.
This notice covers only:

| Bundled file | COD entry and original source |
| --- | --- |
| `vanadium-james-1960-9012770.cif` | [9012770](https://www.crystallography.net/cod/9012770.html), James and Straumanis (1960) |
| `superhydrous-phase-b-koch-muller-2005-mkm0105.cif` | [9007379](https://www.crystallography.net/cod/9007379.html), Koch-Muller et al. (2005) |
| `brucite-zigan-1967-cod-1000054.cif` | [1000054](https://www.crystallography.net/cod/1000054.html), Zigan and Rothbauer (1967) |
| `mgsio3-post-perovskite-murakami-2004-cod-9009217.cif` | [9009217](https://www.crystallography.net/cod/9009217.html), Murakami et al. (2004), Science 304, 855–858, DOI 10.1126/science.1095932 |

[COD's current source declaration](https://www.crystallography.net/cod/#license)
dedicates its data and database to the public domain under
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
The full CC0 text is distributed in `LICENSES/CC0-1.0.txt`.
Retain the original authors and COD identifiers for scientific attribution.

The vanadium copy matches the [Granada COD mirror](https://qiserver.ugr.es/cod/9012770.cif)
byte-for-byte; the superhydrous-phase-B copy matches the
[primary COD file](https://www.crystallography.net/cod/9007379.cif) byte-for-byte.
Primary/mirror variants differ only in SVN date and URL comments. All other
CIF bytes match, including the legacy attribution headers. Those headers and
all structural values are preserved unchanged. Brucite now supplies the active structural reference. A fresh primary COD
download on 5 November matches the preserved archival copy except for SVN
Date/URL comments; all other bytes are identical. The original 1967 article
has not been directly inspected.

The post-perovskite CIF is a fresh, byte-identical primary COD download. The
Granada and PSDI copies differ only in SVN Date/URL comments; all remaining
bytes match. The original legacy attribution header remains unchanged. The
record's current footer explicitly declares COD data/database CC0. Its MgSiO3
Cmcm coordinates are identified as **ab initio positions using the observed
121 GPa, 300 K cell**, not experimentally refined atomic coordinates. The
original Science article has not been freshly re-audited.

COD 9009217 now supplies both the pure post-perovskite topology and the
explicitly labeled pure-MgSiO3 topology proxy for the natural Fe-bearing
post-perovskite material. Reference cells are isotropically scaled, with the
two-step rounding convention recorded in each source block. All pre-existing
cell/site numbers are regenerated exactly. This supplies no Fe/Al/Ca site
partitioning or composition-specific intensity model. Original ICSD-linked
records and earlier reviews are preserved privately; no new grant is claimed
for those earlier acquisition routes. EOS and observation source terms remain
separate from this COD data declaration.

Exact source URLs, hashes, comparison method and limitations are recorded in
`cif-license-review.json`, installed beside this notice. This is a recorded
COD data grant, not a new license for the original journal articles.

The former He (2023) thesis transcription was replaced after 0.17.4. Its CIF
and copied coordinates/refinement are no longer distribution inputs; original
evidence is preserved privately. This COD grant does not extend to that thesis.
The historical unresolved review is retained under `removed_files` in the receipt.
