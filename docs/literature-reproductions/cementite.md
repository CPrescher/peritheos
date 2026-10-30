# Takahashi (2019): nonmagnetic Fe3C source EOSs

## Takahashi: source surface versus conditional coefficient refit

The catalog file is `peritheos/data/materials/cementite_nonmagnetic.eosmat`.
The earlier standalone research snapshot remains available for traceability.
It represents the source-assigned nonmagnetic branch separately from the
existing cementite card. Its zero-pressure V0 values are extrapolated reference
intercepts, not measurements of stable ambient nonmagnetic Fe3C. The KBr and
MgO variants describe the same observations; the latter points to the former
with `pressure_scale_variant_of`. That link does not authorize recalibration.

Equations 1–6 map to BM3 plus `MieGruneisenDebye`, with Tr=300 K,
theta0 fixed to 314 K, n=4 atoms per Fe3C formula, and Z=4 formula units per
cell. Gamma follows gamma0(V/V0)^q; theta follows the integrated Grüneisen law.
Cell volumes convert to molar J/bar with N_A × 10⁻²⁵ / Z.

| Pressure-scale variant | V0 (Å³/cell) | K0 (GPa) | K0′ | gamma0 | q |
| --- | ---: | ---: | ---: | ---: | ---: |
| KBr, Model 1 | 148.8 ± 1.0 | 311 ± 17 | 3.40 ± 0.10 | 1.1 ± 0.4 | 1.9 ± 1.7 |
| MgO, Model 2 | 147.3 ± 1.0 | 323 ± 17 | 3.43 ± 0.09 | 1.4 ± 0.3 | 1.0 ± 1.0 |

All 60 worksheet rows 4–63 are preserved, including the last ten 300 K
SiO2-marker rows, 50 paired KBr/MgO rows, all marker cells and error tokens.
Raw-cell JSON retains row IDs and Excel number formats. For example,
P5 is displayed as 121.30, so Q5 `(8)` means 0.08 Å³, not 0.8 Å³.
The worksheet has a 2600 K row despite the article's 2450 K maximum. It is
retained; no exclusion is inferred. Source error confidence is unspecified,
and missing temperature errors at 300 K remain missing.

The source surfaces reproduce the independently specified 329 GPa, 5000 K
density checks (using Fe=55.845 and C=12.011 g/mol): Model 1 gives 12.8355 g/cm³ against the printed 12.83 ± 0.67;
Model 2 gives 12.6406 against a manual Figure 3a reading of 12.64 with a
0.04 g/cm³ graphical reading bound. These are **extrapolated calculated source
benchmarks**, not measured core states or coefficient-refit parity. Figure 3's
caption reverses the dotted/dashed naming relative to its legend; the visible
5000 K branches have lower density than their 300 K counterparts.

One conditional simultaneous five-parameter refit uses all 60 rows and equal
pressure weights. The paper does not give its exact weights or row mask.
No alternate EOS family, temperature law or optimized exclusion was tried.

| Variant | Published-curve pressure RMS (GPa) | Conditional-refit RMS (GPa) | Conditional K0 (GPa) | Conditional K0′ |
| --- | ---: | ---: | ---: | ---: |
| KBr | 2.4753 | 2.2226 | 266.35 | 3.6567 |
| MgO | 1.9654 | 1.9064 | 350.97 | 3.2894 |

Under this equal-pressure objective, both are **similar**, with central-coefficient differences beyond
printed source widths. Full source/reconstructed coefficients, conditional
errors, covariance, selection and density-check bounds are preserved in
[`docs/data/cementite-reproduction.json`](../data/cementite-reproduction.json).
No refit replaces a source coefficient. These numerical findings do not establish
that the publication is wrong: exact source weighting and selections are missing.

The [30 October Takahashi-specific investigation](takahashi-2019-parity.md)
now demonstrates that both variants reach **conditional central-coefficient
parity** when the printed coordinate-error widths are used, retaining all 60
rows and the source model. Effective-variance P/V/T weighting recovers the KBr
central values at printed rounding. This explains the objective sensitivity of
the earlier equal-weight result; it does not identify the authors' exact
workflow or establish the source uncertainty confidence convention. The historical equal-weight
reconstruction is retained separately; both current catalog assessments and
ledger entries are now **parity**, conditional on that documented objective.
KBr coefficients and uncertainties reproduce published rounding; MgO
uncertainties and exact printed values remain unreproduced. Published
coefficients remain intact and original author workflow remains unconfirmed.

The diffraction structure comes from Prescher Table S1 at 47.4 GPa and
293 K, where the main paper assigns the nonmagnetic branch. Its tabulated axes
are cyclically mapped to standard Pnma:
(a,b,c)=(b,c,a)_table and (x,y,z)=(y,z,x)_table, with no origin shift.
The transformed cell is 4.820 × 6.357 × 4.226 Å (129.4878 Å³), consistent
with reported 129.50(5) Å³. Fe sites 8d and 4c plus C 4c give Z=4 under
explicitly assumed ideal full occupancy; the table does not report refined
occupancies. All original coordinate/U tokens remain in the raw transcription.
This measured finite-pressure structure is not rescaled to either Takahashi V0.

The installed package can load the added material directly:

```python
from peritheos import get_material

material = get_material("cementite_nonmagnetic")
record = material.get_eos_record("cementite_nonmagnetic_takahashi_2019_kbr_bm3_mgd")
pressure_gpa = record.pressure(120.0, 2000.0)
```

KBr is the default within this high-pressure material; that selection is not
an ambient-stability or pressure-calibrant recommendation. Exact marker registry
links and observation-level pressure re-reduction remain unaudited.

The repository-only command `uv run python -m scripts.reproduce_cementite`
replays the audit. It is not an installed-package authoring requirement.
Source checksums, access results, exact locations and literal unmerged values
are in `peritheos/data/research/cementite/sources.json`. The dataset license
note distinguishes source rights from transcription rights.
