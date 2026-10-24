# Catalog coverage overview

Generated from the current `.eosmat` documents and the public catalog API.
This describes the checkout, including local material edits; it is not a release manifest.
Regenerate with `python -m scripts.generate_catalog_coverage`; use `--check`
to detect drift. The [machine-readable snapshot](data/catalog-coverage.json)
includes every material and record, plus a fingerprint of the material files.

Use the [active backlog](material-eos-candidates.md) for priorities and next actions.
Scientific validation, fit reproduction, and availability answer different questions.
A missing thermal EOS is a coverage gap, not a defect in a static record.

## Current inventory

| Measure | Count |
|---|---:|
| Stored material documents | 233 |
| Stored EOS records | 678 |
| Executable materials | 230 |
| Executable EOS records, including references and qualified explicit-selection records | 675 |
| Source EOS choices in normal discovery | 624 |
| Published pressure-scale reference variants | 51 |
| Materials with an executable thermal EOS | 63 |
| Executable records with parameter covariance | 17 |
| Executable records with known parameter-error confidence | 87 |

## Scientific validation of stored records

| Status | Records |
|---|---:|
| `deferred` | 3 |
| `not_reproduced` | 18 |
| `primary_source_validated` | 657 |

## Executable equation kinds

| Status | Records |
|---|---:|
| `hugoniot` | 3 |
| `isothermal` | 538 |
| `thermal` | 134 |

## Pressure-calibration provenance of executable records

| Status | Records |
|---|---:|
| `not_applicable` | 163 |
| `partially_resolved` | 187 |
| `resolved` | 238 |
| `unresolved` | 87 |

## Observation pressure-recalculation readiness

| Status | Records |
|---|---:|
| `missing_calibrant_observations` | 195 |
| `not_applicable` | 160 |
| `not_possible` | 123 |
| `ready` | 136 |
| `reference_eos_not_bundled` | 54 |
| `reference_model_not_supported` | 7 |

Readiness values are the recorded assessments. A missing reference can mean a
genuinely absent scale or a stale/unresolved link to an existing scale. `ready`
does not establish original-fit parity or complete calibration-error propagation.
P–T envelopes are marginal coverage, not measured rectangles or phase-stability fields.

## Reproduction ledger and reconciliation

Outcomes below describe the checked-in refit ledger, not an inferred assessment
of records missing from that ledger. No scientific classifications are changed here.

| Outcome | Records |
|---|---:|
| `not_refittable` | 177 |
| `parity` | 262 |
| `parity_not_achieved` | 62 |
| `similar` | 160 |
| `source_reconstruction` | 17 |

**Catalog records absent from the refit ledger:** none.

**Refit entries absent from the catalog:** none.

**Catalog records absent from the source-audit index:** none.

**Source-audit entries absent from the catalog:** none.

## Material coverage

Counts of isothermal, thermal, and Hugoniot records include executable records only.
Qualified means `not_reproduced`; deferred records remain stored source evidence.
Dataset counts are declarations, not unique independent experiments or completeness claims.
Multiple temperature slices do not imply a continuous thermal EOS.

| Material | Formula | Stored | Executable | Isothermal | Thermal | Hugoniot | Qualified | Deferred | Datasets |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `akimotoite` | MgSiO3 | 11 | 11 | 11 | 0 | 0 | 0 | 0 | 5 |
| `al2o3_perovskite` | Al2O3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| `al2o3_post_perovskite` | Al2O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `al_bearing_phase_d` | Mg0.90Al0.64Si1.29H3.10O6 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `almandine` | Fe3Al2Si3O12 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `alpha_quartz` | SiO2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `alumina` | Al2O3 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 2 |
| `alumina_rh2o3_ii` | Al2O3 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `aluminum` | Al | 5 | 5 | 4 | 1 | 0 | 0 | 0 | 1 |
| `aragonite` | CaCO3 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `argon_fcc` | Ar | 15 | 15 | 11 | 4 | 0 | 2 | 0 | 17 |
| `argon_hcp` | Ar | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `b4c` | B4C | 3 | 3 | 0 | 3 | 0 | 0 | 0 | 1 |
| `beryllium_hcp` | Be | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `boron_nitride` | BN | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `boron_nitride_hexagonal` | BN | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `boron_phosphide` | BP | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `breyite` | CaSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `bridgmanite` | MgSiO3 | 24 | 24 | 21 | 3 | 0 | 0 | 0 | 9 |
| `brucite` | Mg(OH)2 | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 3 |
| `ca0988mg0918fe0078mn0016c2o6_dolomite` | Ca0.988Mg0.918Fe0.078Mn0.016C2O6 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `ca0988mg0918fe0078mn0016c2o6_dolomite_iii` | Ca0.988Mg0.918Fe0.078Mn0.016C2O6 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `ca_perovskite` | CaSiO3 | 29 | 29 | 21 | 8 | 0 | 1 | 0 | 12 |
| `ca_perovskite_tetragonal` | CaSiO3 | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 2 |
| `cairo3_perovskite` | CaIrO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cairo3_post_perovskite` | CaIrO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 2 |
| `calcite` | CaCO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `calcium_carbonate_post_aragonite` | CaCO3 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `calcium_germanate_perovskite` | CaGeO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `calcium_titanate_perovskite` | CaTiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cao` | CaO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cao_b2` | CaO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `casio2o5_titanite` | CaSi2O5 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `casio3_perovskite_tetragonal` | CaSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cementite` | Fe3C | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cerium_dioxide` | CeO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `chromium` | Cr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `cobalt_hcp` | Co | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 3 |
| `coesite` | SiO2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 3 |
| `coesite_iv` | SiO2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 3 |
| `coesite_v` | SiO2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `coo` | CoO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `copper` | Cu | 7 | 7 | 6 | 1 | 0 | 0 | 0 | 2 |
| `cscl` | CsCl | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 2 |
| `cu_handbook_1972_legacy` | Cu | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `delta_alooh` | AlOOH | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `diamond` | C | 7 | 7 | 2 | 5 | 0 | 0 | 0 | 4 |
| `diopside` | CaMgSi2O6 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 3 |
| `e_feooh` | FeOOH | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 2 |
| `e_feooh_hc_low_spin` | FeOOH | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fayalite` | Fe2SiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe` | Fe | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 2 |
| `fe073si027_d03` | Fe0.73Si0.27 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `fe084si016_d03` | Fe0.84Si0.16 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `fe084si016_hcp` | Fe0.84Si0.16 | 2 | 2 | 0 | 2 | 0 | 0 | 0 | 1 |
| `fe088sio3_bridgmanite` | Fe0.88SiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe093o_rhombohedral` | Fe0.93O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe09ni01_hcp` | Fe0.9Ni0.1 | 21 | 21 | 16 | 5 | 0 | 5 | 0 | 1 |
| `fe2o3` | Fe2O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe3o4` | Fe3O4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe3s` | Fe3S | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fe_bearing_al_phase_d` | Al1.53Fe0.22Si0.86O6H3.33 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `fe_fcc` | Fe | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 3 |
| `feh2` | FeH2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `feh3` | FeH3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `feo` | FeO | 5 | 4 | 3 | 1 | 0 | 0 | 1 | 6 |
| `feo_b8_2` | FeO | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `ferropericlase_fe27` | Mg0.725Fe0.268O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `ferrous_ringwoodite` | Fe2SiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fes` | FeS | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `fes_vi` | FeS | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 5 |
| `fesi_b2` | FeSi | 2 | 2 | 0 | 2 | 0 | 0 | 0 | 1 |
| `fesi_b20` | FeSi | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `fesio3_bridgmanite` | FeSiO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| `fesio3_liquid` | FeSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 2 |
| `fesio3_post_perovskite` | FeSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `fesio3_post_perovskite_ii` | FeSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `forsterite` | Mg2SiO4 | 3 | 3 | 1 | 2 | 0 | 0 | 0 | 3 |
| `geo2_rutile` | GeO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `goethite` | FeOOH | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `gold` | Au | 23 | 23 | 14 | 9 | 0 | 1 | 0 | 29 |
| `graphite` | C | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `hydrous_silica_cacl2_nisr_2017` | Si0.954O2H0.184 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 2 |
| `hydrous_stishovite_nisr_2017` | Si0.954O2H0.184 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 2 |
| `ice_vi` | H2O | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `ice_vii` | H2O | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 2 |
| `iceviii` | D2O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `indium_nitride` | InN | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `iridium` | Ir | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `iron` | Fe | 38 | 38 | 26 | 11 | 1 | 3 | 0 | 15 |
| `iron_carbide_fe7c3` | Fe7C3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `kbr_b1` | KBr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `kbr_b2` | KBr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `kcl` | KCl | 10 | 10 | 1 | 9 | 0 | 0 | 0 | 10 |
| `kcl_b1` | KCl | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 2 |
| `klb1_ca_perovskite` | CaSiO3 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `klb1_ferropericlase` | (Mg,Fe)O | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `klb1_mg_perovskite` | (Mg,Fe,Al)(Si,Al)O3 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `larnite` | Ca2SiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `lead_fcc` | Pb | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `lead_hcp` | Pb | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 3 |
| `li_bcc` | Li | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `lif_b1` | LiF | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `magnesiowustite_fe56` | Mg0.423Fe0.541O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `magnesiowustite_fe75` | Mg0.240Fe0.717O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `magnesite` | MgCO3 | 3 | 3 | 2 | 1 | 0 | 0 | 0 | 3 |
| `magnesium_bcc` | Mg | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `magnesium_hcp` | Mg | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `majorite` | MgSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `manganese_alpha` | Mn | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg0125fe0875o` | Mg0.125Fe0.875O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg0215fe0762vac0023o_b1` | Mg0.215Fe0.762□0.023O | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg025fe075o` | Mg0.25Fe0.75O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg0375fe0625o` | Mg0.375Fe0.625O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg042fe058o` | Mg0.42Fe0.58O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg0490fe0483ti0027o` | Mg0.490Fe0.483Ti0.027O | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 2 |
| `mg05fe05al05si05o3_bridgmanite` | Mg0.5Fe0.5Al0.5Si0.5O3 | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 2 |
| `mg05fe05o` | Mg0.5Fe0.5O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg061fe039o` | Mg0.61Fe0.39O | 7 | 7 | 7 | 0 | 0 | 0 | 0 | 2 |
| `mg0625fe0375o` | Mg0.625Fe0.375O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg065fe035o` | Mg0.65Fe0.35O | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg073fe027o` | Mg0.73Fe0.27O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg075fe025o` | Mg0.75Fe0.25O | 6 | 6 | 5 | 1 | 0 | 0 | 0 | 2 |
| `mg075fe025sio3_bridgmanite` | Mg0.75Fe0.25SiO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| `mg080fe020o` | Mg0.80Fe0.20O | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 2 |
| `mg080fe020sio3_bridgmanite` | Mg0.8Fe0.2SiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg083fe017o` | Mg0.83Fe0.17O | 3 | 3 | 2 | 1 | 0 | 0 | 0 | 2 |
| `mg0875fe0125o` | Mg0.875Fe0.125O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg087fe013sio3_bridgmanite` | Mg0.87Fe0.13SiO3 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `mg088fe010al014si090o3_bridgmanite` | Mg0.88Fe0.10Al0.14Si0.90O3 | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 2 |
| `mg090fe009al0005ca0005sio3_bridgmanite` | Mg0.90Fe0.09Al0.005Ca0.005SiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg090fe009al0005ca0005sio3_post_perovskite` | Mg0.90Fe0.09Al0.005Ca0.005SiO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg090fe010o` | Mg0.90Fe0.10O | 15 | 15 | 15 | 0 | 0 | 0 | 0 | 3 |
| `mg092fe008o` | Mg0.92Fe0.08O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg092fe008sio3_bridgmanite` | Mg0.92Fe0.08SiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mg09375al0125si09375o3_bridgmanite` | Mg0.9375Al0.125Si0.9375O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg09375al0125si09375o3_post_perovskite` | Mg0.9375Al0.125Si0.9375O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mg09375fe00625o` | Mg0.9375Fe0.0625O | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 |
| `mg09375fe00625sio3_bridgmanite` | Mg0.9375Fe0.0625SiO3 | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 |
| `mg095al010si095o3_bridgmanite` | Mg0.95Al0.10Si0.95O3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg096875fe003125o` | Mg0.96875Fe0.03125O | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| `mg0991fe0008mn0001co3_magnesite` | Mg0.991Fe0.008Mn0.001CO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg09al02si09o3_bridgmanite` | Mg0.9Al0.2Si0.9O3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `mg09fe01sio3_bridgmanite` | Mg0.9Fe0.1SiO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 2 |
| `mg7si2o8` | Mg7Si2O14H6 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mgal00625h00625si09375o3_bridgmanite` | MgAl0.0625H0.0625Si0.9375O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgal00625h00625si09375o3_post_perovskite` | MgAl0.0625H0.0625Si0.9375O3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgal0125si0875o29375_bridgmanite` | MgAl0.125Si0.875O2.9375 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgal0125si0875o29375_post_perovskite` | MgAl0.125Si0.875O2.9375 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgal2o4_cafe2o4` | MgAl2O4 | 4 | 4 | 2 | 2 | 0 | 0 | 0 | 3 |
| `mgal2o4_cati2o4` | MgAl2O4 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 2 |
| `mgfe60o` | Mg0.4Fe0.6O | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 1 |
| `mgfe94o_b1` | Mg0.058Fe0.942O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mgfe94o_rhombohedral` | Mg0.058Fe0.942O | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `mggeo3_post_perovskite` | MgGeO3 | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 2 |
| `mgo` | MgO | 37 | 37 | 30 | 6 | 1 | 1 | 0 | 19 |
| `mgo_b2` | MgO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgo_liquid` | MgO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `mgsio3_liquid` | MgSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 2 |
| `mgsio3_post_perovskite` | MgSiO3 | 14 | 14 | 12 | 2 | 0 | 0 | 0 | 9 |
| `mgsioh6_365a_phase` | MgSiH6O6 | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 1 |
| `molybdenum` | Mo | 9 | 9 | 8 | 1 | 0 | 0 | 0 | 6 |
| `molybenum_carbide_mo2c` | Mo2C | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `na088al099fe013si094o4_calcium_ferrite` | Na0.88Al0.99Fe0.13Si0.94O4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `na093al102si100o4_calcium_ferrite` | Na0.93Al1.02Si1.00O4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `naalsi2o6` | NaAlSi2O6 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `naalsio4_calcium_ferrite` | NaAlSiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `nacl_b1` | NaCl | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `nacl_b2` | NaCl | 14 | 14 | 12 | 2 | 0 | 1 | 0 | 7 |
| `namg2al5sio12_cf` | NaMg2Al5SiO12 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `namg2al5sio12_nal` | NaMg2Al5SiO12 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `neon_fcc` | Ne | 3 | 3 | 1 | 2 | 0 | 0 | 0 | 5 |
| `nickel` | Ni | 6 | 6 | 5 | 1 | 0 | 0 | 0 | 4 |
| `nickel_oxide` | NiO | 2 | 2 | 1 | 0 | 1 | 0 | 0 | 1 |
| `nickel_oxide_b1` | NiO | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `niobium` | Nb | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `nis` | NiS | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `orthoenstatite` | MgSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `osmium` | Os | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `palladium` | Pd | 7 | 7 | 7 | 0 | 0 | 3 | 0 | 3 |
| `pbs_b1` | PbS | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `perovskite_orthorhombic` | Mg0.88Fe0.12SiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `phase_d` | Mg1.11Si1.89H2.22O6 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `phase_e_fe0072` | Mg1.96Fe0.072Si1.04H3.7O6 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `phase_egg` | Al0.98Si0.92H1.39O4 | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 3 |
| `phase_h` | MgSiO4H2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `platinum` | Pt | 21 | 21 | 13 | 8 | 0 | 1 | 0 | 17 |
| `praseodymium_dioxide` | PrO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `pseudowollastonite` | CaSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `pyrope` | Mg3Al2Si3O12 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `qandilite` | Mg2.00Ti1.00O4 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `rbbr` | RbBr | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `rbcl` | RbCl | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 3 |
| `rbi` | RbI | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `rhenium` | Re | 6 | 6 | 5 | 1 | 0 | 0 | 0 | 8 |
| `rhodium` | Rh | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `ringwoodite` | Mg2SiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `ruthenium` | Ru | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `seifertite` | SiO2 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `siderite_fe095mn005` | Fe0.95Mn0.05CO3 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 1 |
| `silica_cacl2` | SiO2 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `silicon` | Si | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silicon_carbide_b1` | SiC | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silicon_carbide_b3` | SiC | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silicon_v` | Si | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silicon_vii` | Si | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silicon_x` | Si | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `silver` | Ag | 6 | 6 | 5 | 1 | 0 | 0 | 0 | 2 |
| `sio2_stv_andr` | SiO2 | 7 | 7 | 6 | 1 | 0 | 0 | 0 | 4 |
| `sno2` | SnO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `sno2_cubic_27gpa` | SnO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `sno2_pa_3_at_48gpa` | SnO2 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `sro` | SrO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `sro_b2` | SrO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `srsio3_6h_perovskite` | SrSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `srsio3_cubic_perovskite` | SrSiO3 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 1 |
| `superhydrous_phase_b_lt` | Mg10Si3O18H4 | 7 | 7 | 2 | 5 | 0 | 0 | 0 | 2 |
| `tantalum` | Ta | 7 | 7 | 6 | 1 | 0 | 0 | 0 | 3 |
| `tho2_cotunnite` | ThO2 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `thorianite` | ThO2 | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 1 |
| `titanium_alpha` | Ti | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `titanium_omega` | Ti | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `tungsten` | W | 7 | 7 | 6 | 1 | 0 | 0 | 0 | 4 |
| `vanadium_bcc` | V | 1 | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| `wadsleyite` | Mg2SiO4 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| `wollastonite` | CaSiO3 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| `zinc_hcp` | Zn | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 2 |
| `zinc_oxide_rocksalt` | ZnO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `zinc_oxide_wurtzite` | ZnO | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `zircon` | ZrSiO4 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `zirconium_alpha` | Zr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `zirconium_beta` | Zr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| `zirconium_omega` | Zr | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
