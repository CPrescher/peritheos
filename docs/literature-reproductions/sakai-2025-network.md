# Sakai (2025): measured-pair sequential calibration

Audit date: 2026-10-19. **Eight conditional subset calibration stages are executable.** Published catalog coefficients and classifications are preserved. A coefficient mismatch on an incomplete subset does not establish failure of the full published fit.

[Source article](https://doi.org/10.1038/s43246-025-00792-5) · [Nine-material audit](sakai-2025-nine-materials.md) · [Full five-family Au check](sakai-2025-gold.md) · [Numerical results, hashes and every residual](../data/sakai-2025-network-reconstruction.json)

## Inputs and method

The checksum-pinned Supplementary Data 1 cell archive supplies **393 measured volume pairs** on S1–S7, including all 48 Mo–W and 19 Fe–W pairs. S6 stores W volume in column J; D/F are Fe lattice a/c. Original numeric tokens, source cells, run markers and missing errors are retained in `peritheos/data/datasets/sakai-2025-network-paired-volumes.csv`. These rows overlap the earlier Au ancestor archive and represent the same experiments.

The NaCl control adds 103 existing Sakai2011/Dorfman2012 measured pairs and the 32.3 GPa / 27.25 Å³ Matsui reference point identified in Table S8. This is one equal-residual reference constraint, not a new paired observation or hard equality. Sata2002 peak-combination choices and the other Table S8 earlier-study constraints are omitted from these subset controls. These omissions describe this audit's scope, not exclusions attributed to the original authors.

Cu stays at the published Table 2 anchor. Published-upstream controls hold the other Table 2 calibrants fixed. Sequential fits carry fitted Re into Pt, Pt into MgO, W into Au/Mo/Fe, and fitted Pt/MgO/Au into NaCl. V0/K0/K∞′ remain fixed except that Fe and NaCl refine V0/K0/K0′, as in Table 2. **No Table S9 calculated output is fitted.** Equal pressure residuals are an explicit analyst assumption. Numerical bounds keep positive K0 and reference volume above observed volumes; the implemented finite K′/V0/K0 bounds are diagnostic safeguards, not source constraints.

**Errors:** parameter errors below are conditional residual-scaled regression standard errors. They omit upstream uncertainty, shared calibrant covariance and fixed-input errors. Source `err.` fields do not establish this statistical convention. RMS is an unweighted descriptive pressure error. Residual = Pmodel(Vtarget) − Pcalibrant(Vcalibrant). Both curves use identical inputs. Pressure is horizontal and volume vertical. The NaCl reference constraint enters its RMS.

## Results

Re from S1 gives K0′=4.410969, rounding to the published 4.411 under the stated Cu anchor and objective. This is conditional central-value agreement. Other subsets show the effect of restricted observations and sequential calibration; the earlier full five-family Au central-value parity result remains separate. Original weights, covariance, complete historical constraints and full-network workflow parity remain unresolved.

### Recovered historical constraints

The third mode adds traceable earlier inputs: Dewaele2004 Cu–W rows paired by their shared Table I coordinates, Anzellini2014 ruby run 1 and measured W run 3, Dewaele2004 Pt/ruby, Dewaele2000 300 K MgO/Pt rows with a reported Pt lattice parameter, Dorfman2012 MgO/Pt and Mo/MgO pairs, Dewaele2008 Mo/ruby, the previously audited five-family 210-state Au selection, Dewaele2006 hcp-Fe/W, and Sakai2014 300 K Fe/NaCl-B2. Blank calibrant coordinates and Fe/NaCl-B1 states are excluded; they are not replaced by back-calculated marker volumes. Ruby pressures are rounding-limited coordinate reconstructions using the source Mao1986 or Dorogokupets2007 equations and Shen2020, not recovered optical spectra.

This historical control fits W before Re because Re's older input uses W. It then carries Re → Pt → MgO, W/Re/MgO → Au, W/MgO → Mo, Pt/MgO/Au → NaCl and W/NaCl → Fe. Its W, Re and Au coefficients move close to their source values. The upstream-fitted Au K0′ is approximately 5.78075; the earlier 5.78010 central-value parity result conditions on published calibrants and remains distinct. Some original constraints, row masks and weights remain unavailable or ambiguous, so this is still a conditional reconstruction rather than full author-workflow reproduction. The source hashes and every contributing source/row identity are in the numerical report.

### Published upstream

| Material | Inputs | Fitted parameters ± conditional errors | Published RMS (GPa) | Refit RMS (GPa) |
|---|---:|---|---:|---:|
| rhenium | 139 | K0_prime = 4.411 ± 0.009 | 5.84 | 5.84 |
| tungsten | 54 | K0_prime = 4.418 ± 0.022 | 7.88 | 7.88 |
| platinum | 33 | K0_prime = 5.141 ± 0.016 | 4.37 | 4.32 |
| gold | 19 | K0_prime = 5.713 ± 0.018 | 5.72 | 4.29 |
| molybdenum | 48 | K0_prime = 4.152 ± 0.011 | 4.92 | 4.88 |
| iron | 19 | V0 = 22.0 ± 0.9 Å³/cell; K0 = 170 ± 60 GPa; K0_prime = 5.7 ± 0.7 | 4.22 | 3.84 |
| mgo | 81 | K0_prime = 4.278 ± 0.007 | 2.97 | 2.35 |
| nacl_b2 | 104 | V0 = 38.7 ± 2.3 Å³/cell; K0 = 43 ± 13 GPa; K0_prime = 4.5 ± 0.3 | 3.93 | 3.89 |

<details markdown="1">
<summary>rhenium: curves and residuals</summary>

![rhenium conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-rhenium.svg#only-light)
![rhenium conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-rhenium-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>tungsten: curves and residuals</summary>

![tungsten conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-tungsten.svg#only-light)
![tungsten conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-tungsten-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>platinum: curves and residuals</summary>

![platinum conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-platinum.svg#only-light)
![platinum conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-platinum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>gold: curves and residuals</summary>

![gold conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-gold.svg#only-light)
![gold conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-gold-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>molybdenum: curves and residuals</summary>

![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-molybdenum.svg#only-light)
![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-molybdenum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>iron: curves and residuals</summary>

![iron conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-iron.svg#only-light)
![iron conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-iron-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>mgo: curves and residuals</summary>

![mgo conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-mgo.svg#only-light)
![mgo conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-mgo-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>nacl_b2: curves and residuals</summary>

![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-nacl_b2.svg#only-light)
![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/published_upstream-nacl_b2-dark.svg#only-dark)

</details>

### Sequential measured subsets

| Material | Inputs | Fitted parameters ± conditional errors | Published RMS (GPa) | Refit RMS (GPa) |
|---|---:|---|---:|---:|
| rhenium | 139 | K0_prime = 4.411 ± 0.009 | 5.84 | 5.84 |
| tungsten | 54 | K0_prime = 4.418 ± 0.022 | 7.88 | 7.88 |
| platinum | 33 | K0_prime = 5.141 ± 0.016 | 4.37 | 4.32 |
| gold | 19 | K0_prime = 5.716 ± 0.018 | 5.63 | 4.29 |
| molybdenum | 48 | K0_prime = 4.154 ± 0.011 | 4.95 | 4.88 |
| iron | 19 | V0 = 22.0 ± 0.9 Å³/cell; K0 = 170 ± 60 GPa; K0_prime = 5.7 ± 0.7 | 4.23 | 3.84 |
| mgo | 81 | K0_prime = 4.287 ± 0.007 | 3.19 | 2.36 |
| nacl_b2 | 104 | V0 = 38.6 ± 2.5 Å³/cell; K0 = 43 ± 15 GPa; K0_prime = 4.5 ± 0.3 | 4.39 | 4.33 |

<details markdown="1">
<summary>rhenium: curves and residuals</summary>

![rhenium conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-rhenium.svg#only-light)
![rhenium conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-rhenium-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>tungsten: curves and residuals</summary>

![tungsten conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-tungsten.svg#only-light)
![tungsten conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-tungsten-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>platinum: curves and residuals</summary>

![platinum conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-platinum.svg#only-light)
![platinum conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-platinum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>gold: curves and residuals</summary>

![gold conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-gold.svg#only-light)
![gold conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-gold-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>molybdenum: curves and residuals</summary>

![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-molybdenum.svg#only-light)
![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-molybdenum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>iron: curves and residuals</summary>

![iron conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-iron.svg#only-light)
![iron conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-iron-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>mgo: curves and residuals</summary>

![mgo conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-mgo.svg#only-light)
![mgo conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-mgo-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>nacl_b2: curves and residuals</summary>

![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-nacl_b2.svg#only-light)
![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/sequential_measured_subsets-nacl_b2-dark.svg#only-dark)

</details>

### Sequential with recovered history

| Material | Inputs | Fitted parameters ± conditional errors | Published RMS (GPa) | Refit RMS (GPa) |
|---|---:|---|---:|---:|
| tungsten | 96 | K0_prime = 4.416 ± 0.016 | 5.93 | 5.93 |
| rhenium | 187 | K0_prime = 4.412 ± 0.008 | 5.08 | 5.08 |
| platinum | 69 | K0_prime = 5.142 ± 0.011 | 3.04 | 3.01 |
| mgo | 166 | K0_prime = 4.237 ± 0.009 | 3.34 | 3.32 |
| gold | 210 | K0_prime = 5.781 ± 0.004 | 3.15 | 3.15 |
| molybdenum | 94 | K0_prime = 4.151 ± 0.008 | 3.58 | 3.55 |
| nacl_b2 | 104 | V0 = 38.8 ± 2.3 Å³/cell; K0 = 42 ± 13 GPa; K0_prime = 4.5 ± 0.3 | 3.92 | 3.87 |
| iron | 72 | V0 = 22.71 ± 0.28 Å³/cell; K0 = 145 ± 15 GPa; K0_prime = 5.73 ± 0.23 | 3.14 | 3.11 |

<details markdown="1">
<summary>tungsten: curves and residuals</summary>

![tungsten conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-tungsten.svg#only-light)
![tungsten conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-tungsten-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>rhenium: curves and residuals</summary>

![rhenium conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-rhenium.svg#only-light)
![rhenium conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-rhenium-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>platinum: curves and residuals</summary>

![platinum conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-platinum.svg#only-light)
![platinum conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-platinum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>mgo: curves and residuals</summary>

![mgo conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-mgo.svg#only-light)
![mgo conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-mgo-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>gold: curves and residuals</summary>

![gold conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-gold.svg#only-light)
![gold conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-gold-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>molybdenum: curves and residuals</summary>

![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-molybdenum.svg#only-light)
![molybdenum conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-molybdenum-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>nacl_b2: curves and residuals</summary>

![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-nacl_b2.svg#only-light)
![nacl_b2 conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-nacl_b2-dark.svg#only-dark)

</details>

<details markdown="1">
<summary>iron: curves and residuals</summary>

![iron conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-iron.svg#only-light)
![iron conditional subset curves and residuals](../data/sakai-2025-network/sequential_with_recovered_history-iron-dark.svg#only-dark)

</details>

## Reproduction

Repository-only maintenance commands:

```sh
uv run python -m scripts.reproduce_sakai_2025_network
uv run python -m scripts.report_sakai_2025_network
# Both commands also support --check
```

Source rights and the contributor-only CC0 scope remain those of [the original recovery](sakai-2025-gold.md#rights-and-integration-scope). No author outreach or replacement EOS is part of this audit.

## Public dataset example

The paired observations are accessible through the installed package API:

```python
from peritheos import Material, get_material_document

mgo = Material.from_eosmat(
    get_material_document("mgo"), record_identifiers=["mgo_sakai_2025_rydberg_stacey_1"]
)
pairs = mgo.get_dataset("mgo_sakai_2025_table_s7_measured_pairs")
print(len(pairs))  # 81 measured MgO–Pt volume pairs
print(pairs["target_volume_a3"], pairs["calibrant_volume_a3"])
```
