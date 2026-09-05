# Holland et al. (2013): NCFMAS modified-Tait endmembers

## Audit result

Seven of eight LitCurate rows are accepted as genuine composition/phase-distinct modified-Tait reference isotherms. The Na-bearing `npv` row is held because its complete four-coefficient EOS cannot be established from the primary parameter artifact. The official artifact also reveals that LitCurate mistranscribed `fpv` V0 as 2.534; the source value is 2.548 J bar-1 mol-1.

Primary article: [publisher PDF endpoint](https://academic.oup.com/petrology/article-pdf/54/9/1901/4353550/egt035.pdf), DOI [10.1093/petrology/egt035](https://doi.org/10.1093/petrology/egt035). Complete coefficient authority: the authors' [official HPx-eos/THERMOCALC ds62 bundle](https://hpxeosandthermocalc.org/wp-content/uploads/2025/05/tc-thermoinput-metabasite-2022-01-30.zip). Extracted `tc-ds62.txt` SHA-256: `46fa2ecb4fccdeaf44984da81342455d10e920012ab8d6db35b0bc94a7e0781c`.

Each ds62 block supplies `V0`, `K0`, `K0'`, and the indispensable `K0''`. The ledger retained only the first three, but Peritheos ModifiedTait requires all four. Source conversions are `1 J/bar/mol = 10 cm3/mol`, `1 cm3/mol = 1.6605390671738467 A3/formula`, `10 kbar = 1 GPa`, and `1 kbar^-1 = 10 GPa^-1`. No BM3 surrogate is used.

## Same-DOI row disposition

| LitCurate ID | Endmember | Disposition | Record / audit detail |
|---|---|---|---|
| `litcurate_2d2e9920973cac04` | `mpv`, MgSiO3 perovskite | accepted | `bridgmanite_holland_2013_mpv_modified_tait` |
| `litcurate_30a0eba3590777a5` | `fpv`, FeSiO3 perovskite | accepted with correction | `fesio3_bridgmanite_holland_2013_fpv_modified_tait`; V0 is 2.548, not 2.534 J/bar/mol. |
| `litcurate_6f93835769d546c6` | `apv`, AlAlO3 perovskite | accepted | `al2o3_perovskite_holland_2013_apv_modified_tait` |
| `litcurate_aee72f71b40279a3` | `npv`, Na-bearing perovskite endmember | held | The row is absent from the official ds62 artifact; LitCurate labels V0/K0/K0' assumed and supplies no K0''. Inferring it from `apv` would create a non-source EOS. |
| `litcurate_ec9908a767c7323f` | `cpv`, CaSiO3 perovskite | accepted | `ca_perovskite_holland_2013_cpv_modified_tait` |
| `litcurate_2d3f926a8e10b0e0` | `per`, MgO periclase | accepted | `mgo_holland_2013_per_modified_tait` |
| `litcurate_555f952c38b721f1` | `fper`, FeO periclase endmember | accepted | `feo_holland_2013_fper_modified_tait` |
| `litcurate_d9edde4e235c2ee1` | `stv`, SiO2 stishovite | accepted | `sio2_stv_andr_holland_2013_stv_modified_tait` |

These are internally consistent thermodynamic reference parameterizations, not seven independent compression experiments. The 298.15 K isothermal slices are executable; the broader THERMOCALC heat-capacity/activity models are deliberately outside these records. The reproduction checks every source-to-library unit conversion and a 30 GPa analytic inverse/pressure round trip.

## Reproduction

Run:

```bash
uv run python scripts/reproduce_holland_2013_modified_tait.py
```

## Zotero-ready metadata

```bibtex
@article{HollandEtAl2013NCFMAS,
  author = {Holland, Tim J. B. and Hudson, Neil F. C. and Powell, Roger and Harte, Ben},
  title = {New Thermodynamic Models and Calculated Phase Equilibria in NCFMAS for Basic and Ultrabasic Compositions through the Transition Zone into the Uppermost Lower Mantle},
  journal = {Journal of Petrology},
  year = {2013},
  volume = {54},
  number = {9},
  pages = {1901--1920},
  doi = {10.1093/petrology/egt035}
}
```
