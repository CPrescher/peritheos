# Post-perovskite structures reacquired from COD

Reviewed 5 November 2026 against `62f147de`. The two structures previously
linked to ICSD 158959 now use freshly acquired **CC0 COD 9009217**, AMCSD
0015437. Original full EOSMAT records, rights reviews and source evidence are
preserved privately with checksums. Git history and published archives remain.
No license is invented for the original acquisition route.

| Material | Current structural source and scope |
| --- | --- |
| `mgsio3_post_perovskite` | Pure MgSiO3 Cmcm topology from COD 9009217; source cell isotropically expanded to the existing 162.2 Å³ reference |
| `mg090fe009al0005ca0005sio3_post_perovskite` | Same pure-MgSiO3 positional topology, explicitly a proxy for the natural Fe-bearing material; reference cell scaled to 164.9 Å³ |

## Source and license

The [COD record](https://qiserver.ugr.es/cod/9009217.html) attributes the data to
Murakami, Hirose, Kawamura, Sata and Ohishi (2004), *Science* 304, 855–858,
DOI [10.1126/science.1095932](https://doi.org/10.1126/science.1095932). The
record and [current COD declaration](https://www.crystallography.net/cod/)
explicitly dedicate the data and database under CC0. Its four-site Cmcm model
has Z=4 and a reference state of 121 GPa, 300 K. **The atomic coordinates are
ab initio results using the observed cell**, not experimentally refined
positions. This corrects the previous source wording without changing values.
The original Science article was not freshly re-audited.

The bundled CIF is the unchanged primary COD download:
SHA-256 `b034f666553e0c4c9454139a17f7e44e024d362c9a1a3faedf09f7e2b55e76ce`.
Fresh [Granada](https://qiserver.ugr.es/cod/9009217.cif) and
[PSDI](https://cod.psdi.ac.uk/cod/9009217.cif) copies differ only in complete
SVN Date/URL comment lines. Every other byte, including the legacy attribution
header, agrees. The full CIF is preserved; the source header is not edited to
insert a license. The installed `cif-license-review.json` and COD notice record
the current declaration, exact server hashes and comparison method. This is a
source-declared data license, not a grant for the associated publication or
independent authentication of every underlying contributor right.

## Regeneration and scientific limits

The source cell is a=2.456 Å, b=8.042 Å, c=6.093 Å, all angles 90°. Its metric
volume is 120.343769136 Å³ (the CIF prints 120.344 Å³). Reference axes follow:

```text
pure_axes = round(source_axes * (162.2 / product(source_axes))**(1/3), 10)
proxy_axes = round(pure_axes * (164.9 / product(pure_axes))**(1/3), 10)
```

Rounding is elementwise, to 10 decimal places. The second operation deliberately
starts from the rounded pure reference to reproduce the existing compatibility
convention exactly. Ten displayed decimals are not measurement precision.

| Reference | a (Å) | b (Å) | c (Å) |
| --- | ---: | ---: | ---: |
| Pure reference | 2.7129228447 | 8.8832758621 | 6.7303904287 |
| Natural-material proxy | 2.7278932922 | 8.9322955441 | 6.7675300609 |

All existing lattice, fractional-coordinate, occupancy and Wyckoff values are
reproduced exactly. Expansion of the CIF's 16 symmetry operations gives four
Mg, four Si and twelve O positions per cell. Stoichiometric full occupancies
remain the adopted model; no experimentally refined occupancy uncertainties
are inferred. No coordinate covariance or pressure-dependent coordinates are
supplied. The scaled references are not measured ambient structures or new
EOS fits. The Fe-bearing record does not acquire Fe/Al/Ca site partitioning or
a composition-specific scattering model from a pure-MgSiO3 CIF.

The source and installed runtime metadata retain these distinctions. No new
peak list, diffraction-intensity measurement, EOS parameter or recommended
scientific model is introduced. All EOS records and dataset declarations are
unchanged, including uncertainty conventions and recommendation statuses.

## Preservation, bindings and remaining cases

The [replacement receipt](../data/post-perovskite-cod-replacement.json) records
the original/new document checksums, private preservation-manifest digest,
source retrievals and exact release-document compatibility pairs. Projecting
out only root `notes` and `source` recovers the previous complete content hash
of each record. The catalog scientific fingerprint remains
`cfa0fff415cf0b01eb3fb2f59edfe177aba3b0153d0a78b38cf474f4a6167240`.
Historical scientific receipts and frozen source pins remain unchanged;
compatibility accepts only the two exact reviewed metadata pairs. Unrelated
material files remain byte-identical.

The [component rights inventory](eos-structure-compilation-rights.md) now has
five structural components tied to four named COD resources. The two ICSD
cases are resolved for the current distribution by replacement, while the
**seven Holland/THERMOCALC extracts and ten inherited peak collections remain
unresolved**. The optional observation selection remains 451 material entries
plus two Barker entries, with the same five held original Sakai/Campbell entries.
Core archives include the new CIF and existing full CC0 text; optional data
archives retain the updated scoped notices. No old source grant, publication
license or Debian acceptance is asserted.

Repository maintenance checks:

```bash
uv run pytest -q -W error tests/test_post_perovskite_cod.py tests/test_cif_licensing.py tests/test_eos_structure_rights.py
uv run python -m scripts.review_eos_structure_rights --check
uv run python -m scripts.generate_material_summaries --check
```

An additional frozen-release check currently fails first on the unchanged
`fe_fcc` document; the prior commit's binding helper rejects the same hash.
The stored catalog-coverage snapshot also predates seven already-added dataset
declarations in five other materials. These baseline differences are recorded
separately and are not repinned or regenerated as part of this structure task.
The exact two post-perovskite source pairs, unchanged global scientific digest
and current package exports are independently checked above.
