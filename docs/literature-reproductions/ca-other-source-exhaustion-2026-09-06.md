# CaSiO3, FeSiO3-liquid, and residual-source exhaustion audit

Audit date: 2026-09-06.

This audit resolves the ten LitCurate extracting-paper DOI families assigned to the Ca/other batch after the two Mg-only families were rerouted to their material owners. LitCurate remains discovery evidence only. A record was accepted only when the primary paper defines a complete executable pressure-volume equation, its volume basis, reference state, phase identity, and coefficient attribution.

## Outcome

| Extracting paper | LitCurate rows | Outcome | Production action |
| --- | ---: | --- | --- |
| Matsui (1993), [10.5940/jcrsj.35.190](https://doi.org/10.5940/jcrsj.35.190) | 82-83 | No new record | The source MD row gives `V0=24.56 cm3/mol` and `K0=250 GPa` but no analytical EOS family or `K0'`; the observed row is a cited comparison. |
| Morishima et al. (1994), [10.1029/94GL00844](https://doi.org/10.1029/94GL00844) | 85 | No new record | This is a fixed-pressure thermal-expansion experiment at about 20.5 GPa. The candidate's `K0=263 GPa`, `K0'=3.9` are assumed inputs and `V0` is absent, so they do not define a source EOS. |
| Karki and Crain (1998), [10.1029/98GL51952](https://doi.org/10.1029/98GL51952) | 186-191 | One accepted | Source-owned 0 K cubic CaSiO3 BM3. Five comparison rows remain attributed to Wang, Chizmeshya, Wentzcovitch, Sherman, and Mao. |
| Shim et al. (2002), [10.1029/2002GL016148](https://doi.org/10.1029/2002GL016148) | 407-409 | One accepted | Source-owned room-temperature tetragonal CaSiO3 BM3; `V0` and `K0'` are explicitly adopted/fixed. Two incomplete comparison moduli are rejected under this DOI. |
| Ferre et al. (2009), [10.2138/am.2009.3003](https://doi.org/10.2138/am.2009.3003) | 716-726 | No new record | A dislocation/Peierls-Nabarro study. The source row and ten comparison rows omit `V0` and do not specify an independently executable pressure-volume fit. |
| Ono (2013), [10.3390/e15104300](https://doi.org/10.3390/e15104300) | 831 | One accepted | The candidate static BM3 is adopted from Shim et al. (2000), but Ono's Equation (9) and Table 1 provide a source-owned logarithmic-volume thermal-pressure extension. Both attribution layers are encoded. |
| Kawai and Tsuchiya (2015), [10.1002/2015GL063446](https://doi.org/10.1002/2015GL063446) | 857 | No new record | The candidate is citation-reported and summarizes the authors' earlier thermoelastic parameterization. The primary underlying Kawai-Tsuchiya EOS is already represented under DOI 10.1002/2013JB010905. |
| Mookherjee et al. accepted manuscript (2018), [10.2138/am-2018-6694](https://doi.org/10.2138/am-2018-6694) | 1012-1014 | Alias; no new record | The accepted manuscript instructs citation to final DOI 10.2138/am-2019-6694. Both LP and HP phase-Egg EOS records already use that canonical DOI; the Vanpeteghem comparison is not duplicated. |
| Sun et al. (2019), [10.1029/2018GL081421](https://doi.org/10.1029/2018GL081421) | 1018 | One accepted | The explicit 2500 K BM4 reference isotherm is executable. The source's larger specialized P-V-T thermal term is outside the current schema and is not approximated. |
| Miyajima et al. (2025), [10.1029/2025GL115280](https://doi.org/10.1029/2025GL115280) | 1277-1281 | No new record | The source reports electron-diffraction volumes and a residual-pressure estimate, not an EOS fit. Three volumes are literature comparisons, and the coexisting bridgmanite volume is another single-state observation. |

Net production yield: four new EOS records from four papers, plus six paper-family dispositions with no new record.

## Candidate-level dispositions

- `litcurate_c6b889d4337dad0f` (row 82): rejected as a citation-reported observed comparison without an equation or `K0'`.
- `litcurate_b97f1de5f96da283` (row 83): held as incomplete source MD output (`V0`, `K0`, but no EOS family or `K0'`).
- `litcurate_8c532429033be821` (row 85): rejected; assumed compression coefficients with missing `V0` in a fixed-pressure thermal-expansion paper.
- `litcurate_3b5b744aaef477d8` (row 186): accepted as `ca_perovskite_karki_crain_1998_static_bm3`.
- `litcurate_cbe7970fb4ccc774`, `litcurate_21405bcfd736220e`, `litcurate_7b5d7de321a5f603`, `litcurate_9c7fb407f008bf6a`, and `litcurate_de1c4345673472f3` (rows 187-191): rejected under the extracting DOI because they are explicit literature comparisons.
- `litcurate_5f609165a848d8ec` (row 407): accepted as `casio3_perovskite_tetragonal_shim_2002_bm3_1`.
- `litcurate_e4596ed619dd7a98` and `litcurate_8df6ff440c929ae1` (rows 408-409): rejected as incomplete citation comparisons.
- `litcurate_b1d7c88529d6e84a` (row 716): held because the source dislocation paper does not publish `V0` or an unambiguous pressure-volume equation for the `K0=220 GPa`, `K0'=3.69` pair.
- `litcurate_96ea4a7f394919a1`, `litcurate_aedd1bae376c5729`, `litcurate_973dd53c0bf56a44`, `litcurate_cc7d653ed19b6fcb`, `litcurate_2b020b7325f5cbcf`, `litcurate_f0f2f2828b9aff24`, `litcurate_fd1c797b942620ed`, `litcurate_c34b60a10a461199`, `litcurate_368c46ab97a7f8bb`, and `litcurate_6be15fc1b98a2063` (rows 717-726): rejected as citation-reported comparison pairs with no `V0`.
- `litcurate_85d4d9fe2fc09159` (row 831): the static coefficients are not reassigned to Ono, but the audited paper yields `ca_perovskite_ono_2013_bm3_log_thermal` with explicit source lineage.
- `litcurate_5a79a8f7aaa2a668` (row 857): rejected as a citation-reported earlier parameterization already represented from its primary source.
- `litcurate_dfe076f60c19b6f9` and `litcurate_2d449f80a6565739` (rows 1012-1013): already represented as the canonical 2019 LP/HP phase-Egg records.
- `litcurate_2d02ffd4982dbbbc` (row 1014): rejected under the alias DOI as a Vanpeteghem comparison, already documented in the canonical phase-Egg audit.
- `litcurate_3383b2e20e0f10fc` (row 1018): accepted as `fesio3_liquid_sun_2019_2500k_bm4_1`.
- `litcurate_f5ea2881984b7f1c`, `litcurate_b1c9615e2d5f691d`, `litcurate_a51afdbe647da0f9`, `litcurate_ae7e0c081458274f`, and `litcurate_e1e140ba90fbceeb` (rows 1277-1281): rejected as single-state crystallographic measurements or literature comparisons, not fitted EOS parameterizations.

## Reproduction evidence

`peritheos/data/datasets/fesio3-liquid-sun-2019-table1-pvt.csv` transcribes all 46 numerical P-V-T states and printed pressure standard errors in Sun et al. Table 1. `Vx=40.72 cm3/mol` is source-defined; the formula-unit-volume column is an exact Avogadro conversion. Table 1 does not label which simulations Figure 1 classifies as liquid, so the table is not falsely declared to be the exact regression subset.

Run:

```bash
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen python scripts/reproduce_ca_other_source_exhaustion.py
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen pytest -q tests/test_ca_other_source_exhaustion.py
```

The executable Sun BM4 predicts 1.084850 GPa at the Table 1 2500 K, `V/Vx=1` state, versus 1.07(12) GPa in the paper. Ono's thermal implementation is checked directly against Equation (9). The Karki and Shim records are checked by pressure-volume inversion at their reported maximum pressures.

## Zotero-ready production metadata

- Karki, Bijaya B., and Jason Crain (1998), “First-principles determination of elastic properties of CaSiO3 perovskite at lower mantle pressures,” *Geophysical Research Letters* 25, 2741-2744, DOI 10.1029/98GL51952.
- Shim, Sang-Heon, Raymond Jeanloz, and Thomas S. Duffy (2002), “Tetragonal structure of CaSiO3 perovskite above 20 GPa,” *Geophysical Research Letters* 29, 2166, DOI 10.1029/2002GL016148.
- Ono, Shigeaki (2013), “Elastic Properties of CaSiO3 Perovskite from ab initio Molecular Dynamics,” *Entropy* 15, 4300-4309, DOI 10.3390/e15104300.
- Sun, Yicheng, Huiqun Zhou, Kun Yin, and Xiancai Lu (2019), “First-Principles Study of Thermodynamics and Spin Transition in FeSiO3 Liquid at High Pressure,” *Geophysical Research Letters* 46, 3706-3716, DOI 10.1029/2018GL081421.

The phase-Egg discovery DOI 10.2138/am-2018-6694 must remain an alias of canonical DOI 10.2138/am-2019-6694 rather than a second Zotero item or production paper.
