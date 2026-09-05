# Akber-Knutson et al. (2005): Al-bearing pv and ppv

Primary DOI: <https://doi.org/10.1029/2005GL023192>
Author-deposited final article: <https://authors.library.caltech.edu/records/ac43h-2q991>

## Outcome

All ten source rows in Table 1 were accepted as ten distinct static 0 K GGA BM3 parameterizations. They represent two phases and five substitution states, not repeated runs. The paper explicitly defines a third-order Birch-Murnaghan fit to fixed-volume `E(V)` calculations, reports the parameter triplets, and states that the calculated grids span approximately 105–65% of `V0` and reach just above 200 GPa.

Table 1 normalizes every value to an approximately five-atom basis. Peritheos stores conventional `Z=4` volumes, so each printed `V0` is multiplied by four. The exact OVM and AlHM supercell stoichiometries are reconstructed from the source's 79- and 81-atom counts; LitCurate's generic MgSiO3 labels are not retained for those rows. No numerical `E(V)` observations, weights, residuals, coefficient errors, or covariance are published. Consequently these are parameterization-only records with an independent high-compression envelope check, not direct coefficient refits.

## Candidate-row disposition

| LitCurate row | Source identity | Production record |
|---|---|---|
| `litcurate_8f807ef87788f85f` | MgSiO3 pv, no Al | `bridgmanite_akber_knutson_2005_gga_bm3` |
| `litcurate_2bf137ddeb5434cd` | CCM 6.25 mol.% Al2O3 pv | `mg09375al0125si09375o3_bridgmanite_akber_knutson_2005_gga_bm3` |
| `litcurate_5272b35caf831df6` | CCM 100 mol.% Al2O3 pv | `al2o3_perovskite_akber_knutson_2005_gga_bm3` |
| `litcurate_bf1450f777d5bbe8` | OVM 6.25 mol.% Al2O3 pv | `mgal0125si0875o29375_bridgmanite_akber_knutson_2005_gga_bm3` |
| `litcurate_29261d2e7516b73b` | AlHM 3.125 mol.% Al2O3 pv | `mgal00625h00625si09375o3_bridgmanite_akber_knutson_2005_gga_bm3` |
| `litcurate_02a65cdb4bb78819` | MgSiO3 ppv, no Al | `mgsio3_post_perovskite_akber_knutson_2005_gga_bm3` |
| `litcurate_028263f4d2c0e0fc` | CCM 6.25 mol.% Al2O3 ppv | `mg09375al0125si09375o3_post_perovskite_akber_knutson_2005_gga_bm3` |
| `litcurate_6e8b79e9dce41899` | CCM 100 mol.% Al2O3 ppv | `al2o3_post_perovskite_akber_knutson_2005_gga_bm3` |
| `litcurate_2dc493b3352b5541` | OVM 6.25 mol.% Al2O3 ppv | `mgal0125si0875o29375_post_perovskite_akber_knutson_2005_gga_bm3` |
| `litcurate_0a5eeb36d988137b` | AlHM 3.125 mol.% Al2O3 ppv | `mgal00625h00625si09375o3_post_perovskite_akber_knutson_2005_gga_bm3` |

## Numerical check

`scripts/reproduce_akber_knutson_2005_al_perovskites.py` independently evaluates the standard BM3 expression. At `V/V0=0.65`, all ten records produce 205.91–230.83 GPa, reproducing the source's independently stated endpoint “just above 200 GPa.” At `V/V0=1.05` they give −8.55 to −10.44 GPa. This checks the equation convention, units, and fourfold volume conversion, while honestly remaining weaker than a direct `E(V)` refit.

## Zotero-ready citation

Akber-Knutson, S., G. Steinle-Neumann, and P. D. Asimow (2005). “Effect of Al on the sharpness of the MgSiO3 perovskite to post-perovskite phase transition.” *Geophysical Research Letters* 32, L14303. DOI: `10.1029/2005GL023192`.
