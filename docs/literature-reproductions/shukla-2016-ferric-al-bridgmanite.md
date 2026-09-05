# Shukla et al. (2016): Fe3+- and Al-bearing bridgmanite

## Audit result

No production EOS record is accepted from this DOI. The three direct calculated compositions are scientifically distinct, but Table 1 reports only V0, K0, and K0' while the primary article does not state the cold-compression functional form. Treating those triples as BM3 or Vinet would be an unsupported model assignment. Three other rows are explicit comparisons to earlier publications.

Primary source: [free-access publisher article](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016GL069332), DOI [10.1002/2016GL069332](https://doi.org/10.1002/2016GL069332). The article describes LDA+Usc quasiharmonic calculations on 40-atom supercells for `x=0.125`, with 10-12 pressure points, but calls only on a high-temperature EOS for density and does not name its analytical cold curve.

## Same-DOI row disposition

| LitCurate ID | Table 1 row | Disposition | Reason |
|---|---|---|---|
| `litcurate_ae6fe425cad16fdd` | pure bdg, x=0 | rejected as citation row | Value attributed to earlier work, not newly fitted here. |
| `litcurate_07a276bc1acd6da1` | Fe2+-bdg, x=0.125 | rejected as citation row | Value attributed to Shukla et al. (2015), a separate primary DOI. |
| `litcurate_fc233ec735390ab4` | Al-bdg, x=0.125 | held | Direct source row, but equation family is unresolved. |
| `litcurate_7b772f4f43012b67` | experimental Al-bdg | rejected as citation row | Comparison from Jackson et al. (2004, 2005), not a new fit. |
| `litcurate_477fce3ed82ea93a` | Fe3+-bdg, x=0.125 | held | Direct source row, but equation family is unresolved; the spin crossover makes arbitrary extrapolation especially unsafe. |
| `litcurate_668d496ea6b360eb` | Fe3+-Al-bdg, x=0.125 | held | Direct source row, but equation family is unresolved. |

The held rows are useful targets if the authors' original EOS input/output or supporting files can later establish the equation. Their reported triples are preserved in the LitCurate ledger and are not silently discarded or converted.

## Zotero-ready metadata

```bibtex
@article{ShuklaEtAl2016FerricAlBdg,
  author = {Shukla, Gaurav and Cococcioni, Matteo and Wentzcovitch, Renata M.},
  title = {Thermoelasticity of Fe3+- and Al-Bearing Bridgmanite: Effects of Iron Spin Crossover},
  journal = {Geophysical Research Letters},
  year = {2016},
  volume = {43},
  number = {11},
  pages = {5661--5670},
  doi = {10.1002/2016GL069332}
}
```
