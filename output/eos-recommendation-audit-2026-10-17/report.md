# Catalog recommendation audit — 17 October 2026

> Historical screening before the authorized changes. The seven designations,
> wadsleyite hold and fifteen audit repairs have now been implemented. See the
> [current dispositions](../../docs/research/eos-recommendation-dispositions-2026-10-17.md)
> and `implementation-verification.json`. Counts and assessments below preserve
> the pre-change evidence; `audit.py` expects that pre-change catalog state.

This is a read-only assessment under the new EOS curation stopping rule. It
screens all **663 stored EOS records in 232 material documents**, using the
current scientific ledgers, record metadata and existing source audits.
It includes sixteen direct evaluations of existing coefficient sets against
bundled observations. **No EOS was refitted, corrected, removed or reclassified.**
This is a catalog-wide screening, not a fresh independent primary-paper audit
of every record. Recommendation decisions below are curation judgments based
on the stated evidence, not findings that an author's experiment is wrong.

## Recommended new designations

**Seven records should be explicitly marked not recommended under the new rule.**
All seven currently have `scientific_validation.status=primary_source_validated`.
None is its material's current default. Source transcription checks remain
useful, but that status alone conceals unresolved agreement with the source data.

| Source parameterization | Records | Evidence for the designation |
|---|---:|---|
| Fu (2023), published cubic CaSiO3 BM3-MGD | 1 | The completed normalization audit gives 5.917 GPa pressure RMS on the 140 selected Sun observations at physical n=5. Source units, formula/molar conversion and numerical evaluator have been cross-checked. The high-temperature discrepancy is unresolved. |
| Baty (2024), experimental Pd BM3 | 1 | All 78 official rows, volume conversion and BM3 equation checked. Published-curve RMS is 1.2731 GPa; diagnostic RMS is 0.7043 GPa. Existing pressure/volume/error diagnostics do not recover the stated parameter triplet. |
| Sakai (2014), five Fe0.9Ni0.1 thermal variants | 5 | Direct replay of all 32 P7 hot rows gives 7.085–7.214 GPa RMS and maxima 11.109–12.164 GPa, versus approximately 1.1 GPa reported in Table 7. Correct Table 8 extrapolation checkpoints do not resolve the data discrepancy. |

The affected identifiers are:

- `ca_perovskite_fu_2023_bm3_mgd_refit`
- `fe09ni01_hcp_sakai_2014_type1_thermal`
- `fe09ni01_hcp_sakai_2014_type2_thermal`
- `fe09ni01_hcp_sakai_2014_type4_1_thermal`
- `fe09ni01_hcp_sakai_2014_type4_2_thermal`
- `fe09ni01_hcp_sakai_2014_type4_3_thermal`
- `palladium_baty_2024_bm3_1`

The Fu identifier's historical `_refit` suffix refers to **Fu's published
literature-data parameterization**, not the independent Peritheos candidate-data
refit. The Pd finding concerns the experimental BM3 only; Baty's separate 0 K
DFT record and the current Pd default are outside this finding. The Sakai
finding concerns the five thermal variants, not an automatic rejection of
all source isotherms or NaCl calibrations.

Evidence: [Fu source and normalization audit](../../docs/literature-reproductions/fu-2023-eos.md),
[Pd source audit](../../docs/literature-reproductions.md#palladium-baty-2024),
[Sakai source audit](../../docs/literature-reproductions/iron-fei-input-papers.md).
The current implementations should not be recommended until their discrepancies
are reconciled. No replacement model is required by this audit.

## Priority hold on a current default

`wadsleyite_katsura_2009_bm3_1` is the **current wadsleyite default**. Direct replay
of all 85 source-table states confirms 1.89443 GPa pressure RMS and 3.46380 GPa
maximum absolute residual. The existing audit identifies unresolved pressure
calibration and a possible MGD reference-state/Debye convention mismatch.

I recommend withholding a recommendation of the **current implementation**
until that mapping check is complete. This is provisional: the present evidence
must not be presented as a demonstrated error in Katsura's publication. This
case is distinct from the seven source discrepancies above and should be the
first targeted mapping review because it affects a material default.

## Existing flags and holds

Eleven other records already have `not_reproduced`, `deferred`, or an explicit
not-recommended warning. Preserve those restrictions and make the recommendation
wording consistent. A reproduction flag or a nondefault position is not itself
a clear user-facing recommendation statement.

| Source / record group | Count | Current treatment |
|---|---:|---|
| Maltby (2024), fcc Ar | 1 | `not_reproduced`; unresolved published sample calculation; no structured usage recommendation |
| Chen (2010), Ar printed-constant BM3 reconstruction | 1 | `not_reproduced`, explicit DO NOT USE warning; diagnostic reconstruction, not a verified author equation |
| Miozzi (2020), hcp Fe BM3 / Vinet / thermal BM3-MGD | 3 | `not_reproduced`; cold records already prohibit pressure calibration/core prediction; thermal wording should be equally explicit |
| Yokoo (2009), Au/Pt table-output reconstructions | 2 | `not_reproduced`, qualified diagnostics; no structured usage recommendation; analytical source PVT remains unreproduced |
| Campbell (2009), fcc Fe and FeO thermal | 2 | `deferred`, excluded from executable discovery |
| Morard (2026), FeS | 1 | `deferred`, explicitly prohibited for prediction/calibration; coefficient similarity does not resolve the published ±3 GPa residual-envelope discrepancy |
| Gleason (2008), goethite | 1 | Already `not_recommended_for_quantitative_use`, but scientific status remains `primary_source_validated`; transcription and recommendation semantics need to stay distinct |

These eleven exclude brucite, whose unresolved thermal curves/refits were already
withdrawn. They also exclude two Funamori (1998) MgAl2O4 endpoint records that
are already not recommended for a separate historical-data-quality reason.

## Do not use the 75 ledger failures as a blacklist

The existing coefficient-parity ledger has 75 `parity_not_achieved` entries.
Their evidence does not support 75 automatic recommendation withdrawals:

| Assessment of those 75 entries | Count |
|---|---:|
| Recommend a new not-recommended designation | 7 |
| Already flagged or held (subset of the eleven above) | 6 |
| Provisional implementation hold pending mapping review | 1 |
| Generic audit comparison defects / superseded comparisons | 15 |
| Insufficient evidence for an automatic downgrade | 46 |

The last 46 are **not certified safe by this assessment**. They require the
source-specific distinctions in `audit.json`; missing joint-fit inputs,
unknown weights, or poor coefficient identifiability must not be mistaken for
an established failure of a source calculation.

### Confirmed comparison defects

| Generic diagnostic | Defect / stronger existing evidence |
|---|---|
| Jacobsen (2002), three compositions | Uses quartz calibrant volume as sample volume. Correct sample a^3 gives published-curve RMS 0.00506–0.01287 GPa, versus the ledger's 23.3–25.5 GPa. |
| Jacobsen (2005), two B1 compositions | Uses the Fe0.93O rhombohedral volume for both cubic samples. Correct sample columns and source masks must be used. |
| Finkelstein (2017), helium/neon cubic | Uses `hexagonal_c_over_a` as volume and pools phase/medium branches. Correct replay gives 0.27792 / 0.11825 GPa RMS, agreeing with the dedicated source audit. |
| Dobrosavljevic (2019), B1/rhombohedral | Pools all 28 observations for each phase. Explicit source fit flags select 14 / 12 rows; their direct RMS values are 0.24056 / 0.72815 GPa. These corrected residuals alone do not establish coefficient parity. |
| Thompson (2017), HC epsilon-FeOOH | Includes five lower-pressure HOC/transition rows in the HC fit. The twelve >=30 GPa HC rows give 0.10819 GPa source-curve RMS, versus the all-row 1.56406 GPa. |
| Gong (2004), bridgmanite | Compares a shock-reduced reference isentrope directly to Hugoniot states, without the required thermal reduction. |
| Dewaele (2019), two Fe scale variants | Generic result uses all 63 upstream rows and lacks the source scale conversion. The existing dedicated 53-row hcp/two-scale audit reports both fits reproduced. |
| Solomatova (2016), two Fp48 reference branches | Fits individual HS/LS branch curves to mixed-spin observations. The dedicated coupled calculation is `similar`; the generic standalone branch failures do not test that model. |

These are ledger/validation issues. **Do not label the corresponding source EOSs
unreproducible on the basis of those comparisons.** Repairing that validation
selection is a separate implementation task; this assessment has not edited it.

### Examples where coefficient non-parity is insufficient

- B4C (Somayazulu 2023): the published heated curve has chi-square/N=0.961;
  q-only non-parity does not show failed data agreement.
- Fei (2007) NaCl-B2: the calculated source curves are reproduced within line
  width, while q is very sensitive to digitization. Do not interpret the
  diagnostic q fit alone as an unreproduced source EOS.
- Fei (2007) Pt: the thermal curve spacing agrees within source readout
  precision; the source's complete fitted/fixed list is not recovered. The
  existing paused fit investigation was not resumed.
- Dorfman (2012): fixed-metal coefficients fall within the authors' approximately
  2% realistic uncertainty and printed table pressures reproduce within rounding.
  Exact free-fit coefficient parity remains unresolved; no blanket six-record
  rejection is supported by that fact alone.
- Dorogokupets (2015) joint thermoelastic fits and the two diamond Helmholtz
  models: available P-V/P-V-T checks are not their complete upstream fit
  objectives. Coesite and hcp Pb also have incomplete source selections.
- Si-VII and the two SnO2 aliases: high-pressure-only refits hit bounds and
  poorly constrain an ambient reference state. This does not identify a source
  error without recovering the source constraints.

Two additional Ar records need recommendation-language review:
`argon_hcp_wittlinger_1997_bm2_1` (the hcp material default) and
`argon_fcc_errandonea_2006_bm3`. Both retain `primary_source_validated` alongside
explicit `reproduction_status=not_reproduced`, largely reflecting incomplete
figure/raw-data recovery. That inconsistency needs a source-benchmark decision;
missing the exact original fit alone is not enough to assert a numerical source
contradiction. Neither was automatically downgraded here.

## Record-level review and reproducibility

[Machine-readable assessment](audit.json) includes every one of the 75 flagged
identifiers, their live default status, existing recommendation/status, decision,
reason, evidence pointer and generic diagnostic column/selection. It also records
all eleven existing restrictions, sixteen direct evaluations and a catalog
fingerprint. The local `audit.py` only reads catalog/source data and evaluates
existing coefficients; it performs no regression or catalog mutation.

| Record | Live material default | Assessment |
|---|---|---|
| `akimotoite_dorogokupets_2015_298k_rydberg_stacey` | No | Source-specific review required |
| `b4c_somayazulu_2023_bm3_1` | Yes | Source-specific review required |
| `bridgmanite_dorogokupets_2015_298k_rydberg_stacey` | No | Source-specific review required |
| `ca_perovskite_fu_2023_bm3_mgd_refit` | No | New not recommended |
| `coesite_i_iii_bykova_2018_300k_bm3` | No | Source-specific review required |
| `coo_clendenen_1966_murnaghan_1` | Yes | Source-specific review required |
| `diamond_benedict_2014_double_debye_4` | No | Source-specific review required |
| `diamond_correa_2008_double_debye_log_moment_5` | No | Source-specific review required |
| `e_feooh_hc_low_spin_thompson_2017_bm3_1` | Yes | Repair audit comparison |
| `fe093o_b1_jacobsen_2005_bm3_1` | No | Repair audit comparison |
| `fe09ni01_hcp_sakai_2014_p1_bm3` | Yes | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p1_v0_dewaele_bm3` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p1_v0_dewaele_vinet` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p1_v0_mao_bm3` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p1_v0_mao_vinet` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p1_vinet` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p4_bm3` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p4_vinet` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p5_bm3` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p6_bm3` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_p6_vinet` | No | Source-specific review required |
| `fe09ni01_hcp_sakai_2014_type1_thermal` | No | New not recommended |
| `fe09ni01_hcp_sakai_2014_type2_thermal` | No | New not recommended |
| `fe09ni01_hcp_sakai_2014_type4_1_thermal` | No | New not recommended |
| `fe09ni01_hcp_sakai_2014_type4_2_thermal` | No | New not recommended |
| `fe09ni01_hcp_sakai_2014_type4_3_thermal` | No | New not recommended |
| `fe_fcc_campbell_2009_bm3_mgd` | No | Existing flag / hold |
| `feo_campbell_2009_bm3_mgd` | No | Existing flag / hold |
| `ferropericlase_fe27_jacobsen_2002_bm3_1` | Yes | Repair audit comparison |
| `goethite_gleason_2008_bm3_1` | Yes | Existing flag / hold |
| `gold_dorfman_2012_tange_mgo_k0_fixed_vinet` | No | Source-specific review required |
| `gold_dorfman_2012_tange_mgo_k0_free_vinet` | No | Source-specific review required |
| `iron_dewaele_2019_dor_vinet` | No | Repair audit comparison |
| `iron_dewaele_2019_mao_vinet` | No | Repair audit comparison |
| `iron_miozzi_2020_bm3` | No | Existing flag / hold |
| `iron_miozzi_2020_bm3_mgd` | No | Existing flag / hold |
| `iron_miozzi_2020_vinet` | No | Existing flag / hold |
| `iron_sakai_2014_p1_bm3` | No | Source-specific review required |
| `iron_sakai_2014_p1_v0_dewaele_bm3` | No | Source-specific review required |
| `iron_sakai_2014_p1_v0_dewaele_vinet` | No | Source-specific review required |
| `iron_sakai_2014_p1_v0_mao_bm3` | No | Source-specific review required |
| `iron_sakai_2014_p1_v0_mao_vinet` | No | Source-specific review required |
| `iron_sakai_2014_p4_bm3` | No | Source-specific review required |
| `iron_sakai_2014_p6_bm3` | No | Source-specific review required |
| `iron_yamazaki_2012_bm3_thermal` | No | Source-specific review required |
| `iron_yamazaki_2012_vinet_thermal` | No | Source-specific review required |
| `lead_hcp_dewaele_2019_dor_vinet` | No | Source-specific review required |
| `lead_hcp_dewaele_2019_mao_vinet` | Yes | Source-specific review required |
| `magnesiowustite_fe56_jacobsen_2002_bm3_1` | Yes | Repair audit comparison |
| `magnesiowustite_fe75_jacobsen_2002_bm3_1` | Yes | Repair audit comparison |
| `mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1` | Yes | Repair audit comparison |
| `mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2` | No | Repair audit comparison |
| `mg0490fe0483ti0027o_solomatova_2016_hs_bm3_reference_1` | Yes | Repair audit comparison |
| `mg0490fe0483ti0027o_solomatova_2016_ls_bm3_reference_2` | No | Repair audit comparison |
| `mg073fe027o_jacobsen_2005_bm3_1` | Yes | Repair audit comparison |
| `mg075fe025o_solomatova_2016_hs_bm3_reference_1` | No | Source-specific review required |
| `mg075fe025o_solomatova_2016_ls_bm3_reference_2` | No | Source-specific review required |
| `mg080fe020o_fei_2007_ls_b1_bm3` | No | Source-specific review required |
| `mg092fe008sio3_bridgmanite_gong_2004_adiabatic_bm3` | Yes | Repair audit comparison |
| `mgfe94o_b1_dobrosavljevic_2019_bm3_1` | Yes | Repair audit comparison |
| `mgfe94o_rhombohedral_dobrosavljevic_2019_bm3_1` | Yes | Repair audit comparison |
| `mgsio3_post_perovskite_dorogokupets_2015_298k_rydberg_stacey` | No | Source-specific review required |
| `molybdenum_dorfman_2012_tange_mgo_k0_fixed_vinet` | No | Source-specific review required |
| `molybdenum_dorfman_2012_tange_mgo_k0_free_vinet` | No | Source-specific review required |
| `nacl_b2_fei_2007` | No | Source-specific review required |
| `nickel_oxide_b1_campbell_2009_bm3_mgd` | Yes | Source-specific review required |
| `palladium_baty_2024_bm3_1` | No | New not recommended |
| `platinum_dorfman_2012_tange_mgo_k0_fixed_vinet` | No | Source-specific review required |
| `platinum_dorfman_2012_tange_mgo_k0_free_vinet` | No | Source-specific review required |
| `platinum_fei_2004_bm3_mgd` | No | Source-specific review required |
| `platinum_fei_2007_vinet_mgd` | No | Source-specific review required |
| `silicon_vii_anzellini_2019_vinet_1` | Yes | Source-specific review required |
| `sno2_cubic_27gpa_ono_2000_bm3_1` | Yes | Source-specific review required |
| `sno2_pa_3_at_48gpa_ono_2000_bm3_1` | Yes | Source-specific review required |
| `wadsleyite_katsura_2009_bm3_1` | Yes | Provisional implementation hold |
