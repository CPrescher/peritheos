# Dewaele 2008 PAW research EOSMAT collection

The 54 Table VI PAW parameterizations are preserved as a separate research
collection. They are excluded from the main EOS database, catalog searches,
default selection and catalog ledgers. The 13 Table IV experimental EOS
remain in the main catalog. Moving files does not change scientific findings.

The collection contains nine native EOSMAT files, their Figure 4 checkpoint
resources, a source-audit ledger and the complete numerical checks.
Fifty sets have graphical source-equation parity; two bcc-Fe LDA sets lack
an independent benchmark and two hcp-Fe LDA sets retain the unresolved
table/semicore-curve identity discrepancy. Graph checks are not PAW refits.

[Primary publication](https://doi.org/10.1103/PhysRevB.78.104102) ·
[Source investigation](../literature-reproductions/dewaele-2008-metals.md) ·
[Collection files](https://github.com/CPrescher/peritheos/tree/main/peritheos/data/research/dewaele-2008-paw)

The full collection is shipped as package data and can be opened explicitly.
These are conventional-cell EOSMAT volumes; the table below also shows the
atomic normalization used in the publication. Parameter errors remain
unreported. Static-lattice and corrected 298 K fits are distinct isotherms,
not a continuous thermal model. The 298 K corrections include zero-point motion.

```python
from importlib.resources import as_file, files
from peritheos import Dataset, Material, load_eosmat

collection = files("peritheos.data").joinpath("research/dewaele-2008-paw")
with as_file(collection.joinpath("silver.eosmat")) as path:
    document = load_eosmat(path)
    material = Material.from_eosmat(document)
checkpoints = Dataset.from_mapping(document["datasets"][0], resource_root=collection)
print(len(material.eos_records), len(checkpoints))
```

For an external copy, keep the nine files and their `datasets/` directory
together and pass that collection directory as `resource_root`.
Deferred and diagnostic Fe restrictions remain in force when opening files.

## Fe

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/fe.eosmat)

<a id="record-fe_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`fe_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.586 | 21.172 | 239 | 1.88 | 0 |

**Outcome:** not_refittable; **scope:** source_equation.

No separately identified with-semicore alpha-Fe LDA finite-compression curve in Figure 4; raw PAW points are not tabulated. Coefficients are preserved for interchange, pending an independent source benchmark.

- No separately identified with-semicore alpha-Fe LDA finite-compression curve in Figure 4; raw PAW points are not tabulated. Coefficients are preserved for interchange, pending an independent source benchmark.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-fe_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`fe_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.631 | 21.262 | 237 | 1.88 | 298 |

**Outcome:** not_refittable; **scope:** source_equation.

No separately identified with-semicore alpha-Fe LDA finite-compression curve in Figure 4; raw PAW points are not tabulated. Coefficients are preserved for interchange, pending an independent source benchmark.

- No separately identified with-semicore alpha-Fe LDA finite-compression curve in Figure 4; raw PAW points are not tabulated. Coefficients are preserved for interchange, pending an independent source benchmark.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-fe_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`fe_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 11.483 | 22.966 | 186 | 5.27 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-fe_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`fe_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 11.546 | 23.092 | 181 | 5.31 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-fe_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`fe_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.881 | 21.762 | 213 | 6.5 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-fe_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`fe_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.931 | 21.862 | 206 | 6.56 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Iron

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/iron.eosmat)

<a id="record-iron_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`iron_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 9.697 | 19.394 | 334 | 4.74 | 0 |

**Outcome:** not_reproduced; **scope:** source_equation.

The hcp-Fe LDA Table VI parameterization disagrees with the thick LDA Figure 4 path assigned to with-semicore PAW beyond digitization and printed-precision bounds. The follow-up finds agreement with the thin path instead (static maximum difference 0.00244 A3/atom; 298 K 0.00138 A3/atom). The source's table/semicore identity remains unresolved; diagnostic only, not recommended for predictions or pressure calibration. See docs/literature-reproductions/dewaele-2008-metals.md#follow-up-hcp-fe-lda-curve-identity.

- Usage recommendation: not_recommended
- The hcp-Fe LDA Table VI parameterization disagrees with the thick LDA Figure 4 path assigned to with-semicore PAW beyond digitization and printed-precision bounds. The follow-up finds agreement with the thin path instead (static maximum difference 0.00244 A3/atom; 298 K 0.00138 A3/atom). The source's table/semicore identity remains unresolved; diagnostic only, not recommended for predictions or pressure calibration. See docs/literature-reproductions/dewaele-2008-metals.md#follow-up-hcp-fe-lda-curve-identity.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-iron_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`iron_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 9.726 | 19.452 | 329 | 4.75 | 298 |

**Outcome:** not_reproduced; **scope:** source_equation.

The hcp-Fe LDA Table VI parameterization disagrees with the thick LDA Figure 4 path assigned to with-semicore PAW beyond digitization and printed-precision bounds. The follow-up finds agreement with the thin path instead (static maximum difference 0.00244 A3/atom; 298 K 0.00138 A3/atom). The source's table/semicore identity remains unresolved; diagnostic only, not recommended for predictions or pressure calibration. See docs/literature-reproductions/dewaele-2008-metals.md#follow-up-hcp-fe-lda-curve-identity.

- Usage recommendation: not_recommended
- The hcp-Fe LDA Table VI parameterization disagrees with the thick LDA Figure 4 path assigned to with-semicore PAW beyond digitization and printed-precision bounds. The follow-up finds agreement with the thin path instead (static maximum difference 0.00244 A3/atom; 298 K 0.00138 A3/atom). The source's table/semicore identity remains unresolved; diagnostic only, not recommended for predictions or pressure calibration. See docs/literature-reproductions/dewaele-2008-metals.md#follow-up-hcp-fe-lda-curve-identity.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-iron_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`iron_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.486 | 20.972 | 229 | 5.44 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-iron_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`iron_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.532 | 21.064 | 224 | 5.46 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-iron_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`iron_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 9.88 | 19.76 | 301 | 4.93 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-iron_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`iron_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 9.913 | 19.826 | 297 | 4.94 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Cobalt

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/cobalt_hcp.eosmat)

<a id="record-cobalt_hcp_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`cobalt_hcp_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.071 | 20.142 | 286 | 2.84 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-cobalt_hcp_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`cobalt_hcp_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.124 | 20.248 | 281 | 2.85 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-cobalt_hcp_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`cobalt_hcp_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.914 | 21.828 | 215 | 4.81 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-cobalt_hcp_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`cobalt_hcp_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.992 | 21.984 | 208 | 4.85 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-cobalt_hcp_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`cobalt_hcp_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.478 | 20.956 | 249 | 4.18 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-cobalt_hcp_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`cobalt_hcp_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.542 | 21.084 | 243 | 4.21 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Nickel

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/nickel.eosmat)

<a id="record-nickel_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`nickel_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.018 | 40.072 | 254 | 4.99 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-nickel_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`nickel_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.099 | 40.396 | 244 | 5.04 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-nickel_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`nickel_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.893 | 43.572 | 195 | 5.17 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-nickel_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`nickel_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 11.008 | 44.032 | 184 | 5.24 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-nickel_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`nickel_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.383 | 41.532 | 227 | 5.08 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-nickel_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`nickel_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 10.477 | 41.908 | 217 | 5.13 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Zinc (hcp)

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/zinc_hcp.eosmat)

<a id="record-zinc_hcp_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`zinc_hcp_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 13.585 | 27.17 | 89 | 5.8 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-zinc_hcp_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`zinc_hcp_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 13.805 | 27.61 | 82 | 5.9 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-zinc_hcp_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`zinc_hcp_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.365 | 30.73 | 57 | 6.04 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-zinc_hcp_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`zinc_hcp_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.76 | 31.52 | 50 | 6.19 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-zinc_hcp_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`zinc_hcp_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 14.205 | 28.41 | 79 | 5.8 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-zinc_hcp_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`zinc_hcp_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 14.468 | 28.936 | 71 | 5.91 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Molybdenum

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/molybdenum.eosmat)

<a id="record-molybdenum_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`molybdenum_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.063 | 30.126 | 301 | 3.99 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-molybdenum_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`molybdenum_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.145 | 30.29 | 295 | 4.01 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-molybdenum_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`molybdenum_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.849 | 31.698 | 262 | 4.14 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-molybdenum_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`molybdenum_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.947 | 31.894 | 255 | 4.17 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-molybdenum_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`molybdenum_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.323 | 30.646 | 287 | 4.02 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-molybdenum_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`molybdenum_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.41 | 30.82 | 280 | 4.05 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Silver

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/silver.eosmat)

<a id="record-silver_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`silver_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.046 | 64.184 | 137 | 5.87 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-silver_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`silver_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.178 | 64.712 | 131 | 5.93 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-silver_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`silver_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 17.907 | 71.628 | 88 | 6.13 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-silver_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`silver_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 18.137 | 72.548 | 82 | 6.23 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-silver_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`silver_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.694 | 66.776 | 116 | 5.98 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-silver_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`silver_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.856 | 67.424 | 110 | 6.05 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Platinum

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/platinum.eosmat)

<a id="record-platinum_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`platinum_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 14.823 | 59.292 | 302 | 5.61 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-platinum_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`platinum_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 14.898 | 59.592 | 293 | 5.65 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-platinum_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`platinum_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.793 | 63.172 | 228 | 6.14 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-platinum_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`platinum_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.899 | 63.596 | 219 | 6.2 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-platinum_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`platinum_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.123 | 60.492 | 282 | 5.6 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-platinum_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`platinum_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 15.204 | 60.816 | 274 | 5.64 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

## Gold

[Native EOSMAT file](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/research/dewaele-2008-paw/gold.eosmat)

<a id="record-gold_dewaele_2008_paw_lda_static_vinet"></a>

### Dewaele et al. (2008), PAW LDA, static lattice, with semicore

`gold_dewaele_2008_paw_lda_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.603 | 66.412 | 195 | 5.85 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-gold_dewaele_2008_paw_lda_298k_vinet"></a>

### Dewaele et al. (2008), PAW LDA, 298 K, with semicore

`gold_dewaele_2008_paw_lda_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 16.724 | 66.896 | 187 | 5.92 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-gold_dewaele_2008_paw_pbe_static_vinet"></a>

### Dewaele et al. (2008), PAW PBE, static lattice, with semicore

`gold_dewaele_2008_paw_pbe_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 17.988 | 71.952 | 140 | 6.07 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-gold_dewaele_2008_paw_pbe_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBE, 298 K, with semicore

`gold_dewaele_2008_paw_pbe_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 18.173 | 72.692 | 132 | 6.15 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-gold_dewaele_2008_paw_pbesol_static_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, static lattice, with semicore

`gold_dewaele_2008_paw_pbesol_static_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 17.122 | 68.488 | 176 | 5.62 | 0 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.

<a id="record-gold_dewaele_2008_paw_pbesol_298k_vinet"></a>

### Dewaele et al. (2008), PAW PBESOL, 298 K, with semicore

`gold_dewaele_2008_paw_pbesol_298k_vinet`

| V0, Å³/atom | V0, Å³/cell | K0, GPa | K0 prime | Reference T, K |
| ---: | ---: | ---: | ---: | ---: |
| 17.261 | 69.044 | 169 | 5.69 | 298 |

**Outcome:** parity; **scope:** source_equation.

Source-equation parity against identified graphical checkpoints within declared reading and printed-coefficient bounds; no independent parameter refit or author-workflow reproduction.

- Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction.
- This named numerical assessment alone does not establish author-workflow reproduction, physical accuracy or recommendation for use.
- No reconstructed coefficient table is recorded for this assessment; no missing values have been inferred.
