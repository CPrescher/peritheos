# Zhang and Wentzcovitch (2022): anharmonic MgSiO3 Pv and PPv

## Audit result

Four direct, source-reported 300 K EOS rows are accepted. The two experimental range rows in the same table are citation summaries and are not duplicated under this DOI.

Primary source: [author manuscript (arXiv 2204.13159)](https://arxiv.org/pdf/2204.13159), checked against DOI [10.1103/PhysRevB.106.054103](https://doi.org/10.1103/PhysRevB.106.054103). Downloaded manuscript SHA-256: `85bb91f0a2832c5328b0d91c55cd41e6feba7dedb633f4c7003a29bfdc80515d`.

The method section reports five volumes and six temperatures between 300 and 5000 K for each exchange-correlation functional. Section III says the Helmholtz free energies were fitted to a third-order Birch-Murnaghan finite-strain expansion. Table I gives separate 300 K Pv/PPv and LDA/PBE coefficient triples in A3 per MgSiO3. The stored conventional-cell volumes multiply these values by `Z=4`; no molar conversion is involved.

## Same-DOI row disposition

| LitCurate ID | Table I row | Disposition | Peritheos record / reason |
|---|---|---|---|
| `litcurate_de4c3c7e3aa017f4` | Pv LDA: 40.15, 255.72, 3.91 | accepted | `bridgmanite_zhang_wentzcovitch_2022_phq_lda_300k_bm3` |
| `litcurate_720fd80b9e9294bd` | Pv PBE: 42.44, 226.47, 3.89 | accepted | `bridgmanite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3` |
| `litcurate_081bf45b6736075c` | Pv experimental ranges | rejected as citation row | Composite minima/maxima from Refs. 33-35 and 44-49, not a fit performed by this paper. |
| `litcurate_faa7c12c60bab5df` | PPv LDA: 40.02, 227.43, 4.24 | accepted | `mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_lda_300k_bm3` |
| `litcurate_e591fafb49a5023c` | PPv PBE: 42.56, 194.91, 4.24 | accepted | `mgsio3_post_perovskite_zhang_wentzcovitch_2022_phq_pbe_300k_bm3` |
| `litcurate_5c4acb917edcc312` | PPv experimental ranges | rejected as citation row | Composite ranges from Refs. 35, 42, and 43, with K0' fixed in the cited experiments. |

No coefficient uncertainty or covariance is published for the four computed rows. Figure 2 supplies plotted curves but not numerical F(V) observations, so the reproduction script checks the exact table coefficients, zero-pressure condition, 100 GPa inversion, and pressure-volume round trip without claiming an independent refit.

## Reproduction

Run:

```bash
uv run python scripts/reproduce_zhang_driver_thermoelastic_eos.py
```

## Zotero-ready metadata

```bibtex
@article{ZhangWentzcovitch2022Anharmonic,
  author = {Zhang, Zhen and Wentzcovitch, Renata M.},
  title = {Anharmonic Thermodynamic Properties and Phase Boundary across the Postperovskite Transition in MgSiO3},
  journal = {Physical Review B},
  year = {2022},
  volume = {106},
  pages = {054103},
  doi = {10.1103/PhysRevB.106.054103}
}
```
