# Dewaele (2019) static-DAC metal EOS audit

## Outcome

The 30 catalog records are the 15 metal rows in Table 1 of Dewaele (2019), each reduced once with the Mao et al. (1986) ruby calibration and once with the Dorogokupets--Oganov (2007) calibration. They are alternative reductions of the same observations, not 30 independent experiments.

This audit independently refits **28 of the 30 records** from tabulated source rows. Every unweighted pressure-residual refit is inside the 95% coefficient errors printed in Dewaele (2019). The recovered Fe EPAPS source resolves both epsilon-Fe records: all 63 rows are bundled with raw gauge observations, the 53 hcp rows reproduce the two 2019 coefficient sets nearly exactly, and those two records are now linked to bundled primary data. The family is nevertheless **not fully reproducible** because the exact hcp-Pb rows from Kuznetsov et al. (2002) and the corrected/unrounded Be row at 24.6 GPa remain unavailable. Twenty-one hcp-Pb square markers are stored as an explicit Figure 3 digitization, including pixel centers, fitted-axis calibration, symbol fill, confidence, and digitization uncertainties. The paper also does not state the least-squares dependent variable, row weights, or covariance treatment, so the broader family must not be presented as fully reproduced or committed as a completed 30-record audit yet.

Executable audit: `scripts/audit_dewaele_2019_static_dac.py`. It emits the 30 identifiers, source-data checksums, fitted coefficients for pressure- and volume-residual objectives, residual norms, and the unresolved source gaps.

## Primary-source chain and redistribution status

| Role | Primary source and locations | Local PDF SHA-256 | Data/redistribution status |
|---|---|---|---|
| Published 30 coefficients and common method | Dewaele, *Minerals* **9** (2019) 684, [doi:10.3390/min9110684](https://doi.org/10.3390/min9110684), Equation 1, Table 1, Sections 2--4, and Pb Table 2 | `8ae1bb7611c6925a6a66ac39ed3b6c5bd435df304e7939f2de0362346cb9466d` | Publisher article is CC BY 4.0. The PDF is not copied into the repository. The 17 Pb Table 2 rows may be redistributed with attribution. |
| Al, Cu, Pt, Ta, W observations | Dewaele, Loubeyre, and Mezouar, *Phys. Rev. B* **70** (2004) 094112, [doi:10.1103/PhysRevB.70.094112](https://doi.org/10.1103/PhysRevB.70.094112), Table I and Section II | `1f9e02e4a77f49a868b1a6261290465bdb00597b1660691b1790c7502ec9aeeb` | APS article; no explicit reusable dataset license located. Only factual numerical table transcriptions are stored, not the article or table image. This is not a license grant for article content. |
| Ag, Co, Mo, Ni, Zn observations and Fe cross-check | Dewaele et al., *Phys. Rev. B* **78** (2008) 104102, [doi:10.1103/PhysRevB.78.104102](https://doi.org/10.1103/PhysRevB.78.104102), Tables II--III | `a06fd5cb9433b99f35b20582be5214b1e805014c150a2907fdf8933776cf6daf` | APS article; same no-explicit-data-license limitation. Factual rows only are transcribed. |
| Au observations | Takemura and Dewaele, *Phys. Rev. B* **78** (2008) 104119, [doi:10.1103/PhysRevB.78.104119](https://doi.org/10.1103/PhysRevB.78.104119), Table III | `f8456c574ec4ec25eaf9561c84160de4182b477366349d816fe927d310e732f0` | APS article; no explicit reusable dataset license located. Existing factual table transcription only. |
| Be observations | Lazicki et al., *Phys. Rev. B* **86** (2012) 174118, [doi:10.1103/PhysRevB.86.174118](https://doi.org/10.1103/PhysRevB.86.174118), Table I; [DOE/OSTI accepted manuscript](https://www.osti.gov/servlets/purl/1226986) | `97383f7156fe51e4037a084ed192a25259831d02f0e29706f804847085be3d28` | Public-access manuscript, but no explicit reusable table-data license located. No Be table transcription is added while the erroneous row remains unresolved. |
| Re observations | Anzellini et al., *J. Appl. Phys.* **115** (2014) 043511, [doi:10.1063/1.4863300](https://doi.org/10.1063/1.4863300), Sections II.C and Table III | inherited source checksum `4975ba6e1e17c0ab7bd27938d8ca3be59fb282283c981e62e99be782518514d5` | AIP article; no explicit reusable dataset license located. Existing factual table transcription only. |
| Complete Fe observations | Dewaele et al., *Phys. Rev. Lett.* **97** (2006) 215504, [doi:10.1103/PhysRevLett.97.215504](https://doi.org/10.1103/PhysRevLett.97.215504), EPAPS `E-PRLTAO-97-048648`; [current APS Supplemental Material landing page](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.97.215504) | article PDF `75a61d6625b075f5a0a3a6bf03c8f3ab4ad1594e3ff064f3a5e6f066f8839c16`; supplied TeX `311b0191c843bca7b2232dfb1faa8e94ee7ceb8d6cae70ce4e630bf50cc9e4fa` | Complete 63-row supplementary table recovered: ten bcc and 53 hcp observations from five helium/neon runs, including reported pressure, ruby wavelength, W lattice parameter, atomic volume, and hcp c/a. The factual table transcription is bundled; no explicit reusable article license was located. |
| Older hcp-Pb observations | Kuznetsov et al., *Solid State Communications* **122** (2002) 125--127, [doi:10.1016/S0038-1098(02)00112-6](https://doi.org/10.1016/S0038-1098(02)00112-6), Figure 3 and Table 2; [author-hosted article copy](https://acce-research.fr/MyRecentPublications/FCC_HCP_Boundary_in%20Lead.pdf) | `c4bfc71c5d4186f7f3159015bafa51552349f9e0f4ba305d13604673d0b3a61d` | Copyright Elsevier 2002; no reusable dataset license located. The 21 visible hcp square markers are digitized as factual coordinates, not as an article image. Open symbols are source-author reductions of high-temperature observations; solid symbols are room-temperature observations. Overlaps are explicitly marked with wider uncertainties and lower confidence. |

The Re source-PDF checksum above is pre-existing catalog metadata; this audit did not obtain a stable official AIP PDF download from which to recompute it. That unresolved checksum provenance is another reason not to make a clean provenance commit.

## Units, cell normalization, and equation

Dewaele (2019) reports `V0` in Å³ per atom for the metals, `P` and `K0` in GPa, and dimensionless `K0'`, all at 300 K. Catalog records use conventional-cell volumes, so the published/refitted atomic `V0` is multiplied by `Z=4` for fcc Au, Ag, Al, Cu, Ni, and Pt, and by `Z=2` for bcc Mo, Ta, and W and hcp Be, Co, Fe, Pb, Re, and Zn. The audit reverses that normalization before comparing coefficients. Au volumes are derived as `a111³/4`; Re conventional hcp volumes are divided by two; the other source tables already print atomic volume.

All 30 records use the Rydberg--Vinet form

`P = 3 K0 (1-x) x^-2 exp[1.5 (K0' - 1) (1-x)]`, with `x=(V/V0)^(1/3)`.

For these 15 metal rows, `V0`, `K0`, and `K0'` are all fitted; none is fixed. Parenthetical Table 1 uncertainties are 95% fit intervals. Dewaele (2019) says the Dorogokupets reductions have the same errors as the Mao reductions, but does not print a covariance matrix.

## Pressure recalibration

The executable calibrations are:

- Mao et al. (1986): `P = (1904/7.665) [(lambda/lambda0)^7.665 - 1]` GPa.
- Dorogokupets--Oganov (2007): `P = 1884 s (1 + 5.5 s)` GPa, where `s=lambda/lambda0-1`.

Conversion is performed through the common wavelength ratio, never by interpolating pressure. For source rows printed on the Mao scale, the audit inverts Mao and evaluates Dorogokupets--Oganov; for Au and Re rows printed on the Dorogokupets scale it performs the inverse operation. Dewaele (2019) explicitly says the older data were recalibrated to these alternatives and that W gauges used for Re and Fe were themselves calibrated against ruby.

There is an important limit: some high-pressure Ni/Zn and Re rows use helium or tungsten EOS gauges rather than a recorded ruby wavelength. Treating their reported source-scale pressure as an invertible effective ruby pressure reproduces the published coefficients, but it is an inference, not a row-wise recalculation from raw calibrant observations. Exact recalibration would require the original gauge observations and the exact historical He/W EOS evaluation path.

## Row selection and independent refits

The audit uses every metal-specific row printed in the cited source table: 86 Au rows; 36 Pt; 42 Cu; 36 Ta; 40 Al; 42 W; 25 Co; 34 Ag; 29 Mo; 42 Ni; 43 Zn; 59 Re; and 53 hcp Fe rows from the recovered EPAPS table. It also fits the 17 open Dewaele (2019) Pb rows. The separately bundled Kuznetsov digitization contains 21 hcp square markers from Figure 3 over 10.527--35.121 GPa: 14 open reduced-high-temperature symbols and seven solid room-temperature symbols. It is not silently substituted for exact rows in the coefficient-parity count because the crowded low-pressure markers overlap fcc circles and fitted curves, and Dewaele (2019) does not identify the exact selected subset, weights, or Brown-NaCl scale handling.

For Fe, the EPAPS table's reported pressure is the ruby-linked effective pressure used by the source for both ruby- and W-gauge rows. Treating it as the Dorogokupets reduction and converting it through the common ruby wavelength ratio for the Mao alternative reproduces Dewaele (2019). This also confirms the effective-pressure conversion used elsewhere in the audit. The raw ruby wavelengths and tungsten lattice parameters remain bundled so alternative gauge recalculations can be tested without revisiting the supplement.

The following are unweighted pressure-residual fits. Values are atomic `(V0, K0, K0')`; all differences from the Table 1 coefficients are smaller than the printed 95% errors.

| Metal | Mao refit | Dorogokupets refit | Row-level disposition |
|---|---|---|---|
| Au | `(16.98316, 166.3667, 5.4661)` | `(16.98578, 163.3770, 6.0298)` | complete table reproduction |
| Pt | `(15.09899, 273.44, 4.8327)` | `(15.09798, 270.80, 5.4976)` | complete table reproduction |
| Cu | `(11.81000, 135.285, 4.9086)` | `(11.8104, 133.127, 5.3786)` | complete table reproduction |
| Ta | `(18.020, 197.84, 3.1734)` | `(18.019, 196.06, 3.6417)` | complete table reproduction |
| Al | `(16.5727, 76.324, 4.1551)` | `(16.5842, 74.188, 4.5244)` | complete table reproduction |
| W | `(15.8619, 298.283, 3.8152)` | `(15.8577, 298.571, 4.3673)` | complete table reproduction |
| Co | `(11.077, 197.03, 3.8491)` | `(11.077, 194.85, 4.3584)` | complete table reproduction |
| Ag | `(17.0697, 100.388, 5.6897)` | `(17.0876, 96.822, 6.2085)` | complete table reproduction |
| Mo | `(15.569, 270.31, 3.3348)` | `(15.567, 269.33, 3.8728)` | complete table reproduction |
| Ni | `(10.9569, 176.049, 4.8675)` | `(10.9546, 174.519, 5.3659)` | complete table reproduction |
| Zn | `(15.146, 64.334, 5.2974)` | `(15.155, 62.255, 5.7036)` | complete table reproduction |
| Re | `(14.7383, 349.467, 4.0128)` | `(14.7356, 349.324, 4.6510)` | rows complete; mixed-gauge recalibration inference |
| Fe | `(11.20801, 164.6379, 4.9621)` | `(11.17599, 168.6495, 5.3281)` | complete 53-row EPAPS hcp selection; both fits reproduced |
| Pb (2019 rows only) | `(28.0643, 71.727, 4.3970)` | `(28.0678, 69.790, 4.7722)` | incomplete selection; older rows absent |

The close pressure-residual reproduction is strong evidence that the published fit minimized pressure residuals without a dominating nonuniform weight. It is not proof: Dewaele (2019) merely says the P--V data were fitted. The audit therefore also reports unweighted volume-residual fits. Those often fall outside the small published coefficient intervals (for example Au, Ag, Al, Cu, Pt, Ta, W, and Zn), establishing that the unspecified objective materially matters.

## Dataset checksums

| Dataset | SHA-256 |
|---|---|
| `aluminum-dewaele-2004-table1-compression.csv` | `61fa04908b3a4827fb6566e4ceea29706ac52cbaae1cee5bb0590d1af44423b7` |
| `cobalt-dewaele-2008-table2-compression.csv` | `7c3645e11f62bf9ea41a1cfafa68bc962e0a1e165756bae2f1b1618fc54cf296` |
| `copper-dewaele-2004-table1-compression.csv` | `0f4535c7e29960690c098408711f63ae559be1d73bba18da73047d75f266120e` |
| `gold-takemura-2008-table3-compression.csv` | `04e2015f65f5203d89966ebf4b08eb3000d9a3fb7bfc7e3fd032e253e72b3858` |
| `iron-dewaele-2006-epaps-compression.csv` | `5babd4420a8d34ba6e5c19aeafb29d7b861b971a2ca9347e65f601bb942b38f9` |
| `lead-dewaele-2019-table2-compression.csv` | `37d1dbcbd19c4e750785552bde3af89c044504c6c5b551a63b7ff9ee0405172f` |
| `lead-hcp-kuznetsov-2002-figure3-digitized.csv` | `d89597df2dcb6736d178652bb831633fe99546a36ebdbe4de9f5dcad8f42de23` |
| `molybdenum-dewaele-2008-table2-compression.csv` | `a518f3f0a53785b2c0d875bb2aee8b05db94daef6e06cdd2e96800a341684e2f` |
| `nickel-dewaele-2008-table2-compression.csv` | `ca114872c320ef2364e5073ac1f05e8bb47d479f3995c6b5dc32f1134b1958be` |
| `platinum-dewaele-2004-table1-compression.csv` | `b0230e9fdd1f8a36cecc35cb79a6d56bd4ee170007a465f1842ef76544c47273` |
| `rhenium-anzellini-2014-table3-compression.csv` | `6ea5f5036862f9131d489d0589ce075848a1a140234d3a3b49239064d8ffc11f` |
| `silver-dewaele-2008-table2-compression.csv` | `b369fc32298a0b4d43d54ff75af524227c9d43ccf66e813a215a4967a1fb157c` |
| `tantalum-dewaele-2004-table1-compression.csv` | `3c69e1a9943dc1f882d60f63d01cf6708c114a4569a20e04a8b2278e5c09bd78` |
| `tungsten-dewaele-2004-table1-compression.csv` | `d73f1d9114c974be80175c6ade256bf17cb4a529e911480e6cd1e337394ef523` |
| `zinc-dewaele-2008-table2-compression.csv` | `9edf4625e61d2fdf4b20fd51c1134b7c45495d99d8070f25b555cebb0b476f11` |

## Exact remaining work

1. Obtain the exact numerical P--V rows behind Kuznetsov et al. (2002), preferably from an author or repository, and identify the subset/weighting and Brown-NaCl scale treatment used by Dewaele (2019). Until then, retain the current CSV's explicit `digitized_from_figure` status and uncertainty/confidence fields.
2. Obtain the corrected or unrounded Be Table I row at 24.6 GPa. The printed `a=2.262 Å, c=3.405 Å` is inconsistent with neighboring compression points and cannot be silently changed. Excluding it reproduces the published coefficients, but exclusion is not documented by the source.
3. Confirm with the authors or fitting code that the dependent variable was pressure, whether rows were weighted, and how 95% intervals were calculated.
4. Resolve explicit permission/redistribution status for the newly transcribed APS tables and replace the abbreviated inherited Re PDF checksum with a verified complete digest from a stable primary download.

Only after those items are resolved should the 30 material records be updated to `primary_data_check: bundled`, linked to exact dataset identifiers/checksums, and marked with independent-refit parity. At that point run the full suite from a clean environment, stage only this audit's files and the 30 narrowly changed records, create a task branch, commit, and push.
