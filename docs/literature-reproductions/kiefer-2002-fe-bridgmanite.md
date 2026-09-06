# Kiefer, Stixrude, and Wentzcovitch (2002) bridgmanite audit

## Primary source

- Primary article: B. Kiefer, L. Stixrude, and R. M. Wentzcovitch, “Elasticity of (Mg,Fe)SiO3-perovskite at high pressures,” *Geophysical Research Letters* **29** (2002), 1539, DOI [10.1029/2002GL014683](https://doi.org/10.1029/2002GL014683).
- Open institutional copy: [UCL Discovery record and PDF](https://discovery.ucl.ac.uk/id/eprint/112632/).
- Checked locations: computational method, Table 1 and its notes, and the paragraph immediately following Table 1.

The source computed static GGA Pbnm perovskite for `x=0` and an ordered `x=0.25` Fe substitution. It explicitly reports complete EOS triplets after Table 1: 159.6 A3, 264.3 GPa, 3.94 for MgSiO3 and 160.5 A3, 270.0 GPa, 3.96 for (Mg0.75Fe0.25)SiO3. These volumes are 20-atom conventional cells (`Z=4`).

Table 1 separately reports Voigt-Reuss-Hill aggregate elastic moduli and their pressure derivatives: 263/4.06 and 268/4.10. Those are said to be *consistent with* the EOS values, not alternative BM3 fits. LitCurate split them into pseudo-records; they must not be accepted.

## All nine candidate rows

| LitCurate row | Disposition | Evidence |
|---|---|---|
| `litcurate_6fc87b3bb31010b7` (398), MgSiO3 EOS 159.6/264.3/3.94 | accepted | Complete source-generated EOS triplet. |
| `litcurate_2add53b06987f13d` (399), Fe25 EOS 160.5/270.0/3.96 | accepted | Complete source-generated composition branch. |
| `litcurate_f948de0084755234` (400), MgSiO3 aggregate 263/4.06 | rejected pseudo-split | Table 1 elastic aggregate; no V0 and not an EOS fit. |
| `litcurate_ff5b474ba6734f85` (401), Fe25 aggregate 268/4.10 | rejected pseudo-split | Table 1 elastic aggregate; no V0 and not an EOS fit. |
| `litcurate_8774ba014956f30e` (402), prior theory 259 GPa | rejected citation row | External comparison `c` in Table 1. |
| `litcurate_75b3b031cfd6ccd6` (403), prior theory 258 GPa | rejected citation row | External comparison `d`. |
| `litcurate_aa7e18d94a2e00f1` (404), prior theory 231 GPa | rejected citation row | External comparison `e`. |
| `litcurate_2c92df1a000e5fdc` (405), Brillouin 264 GPa | rejected citation row | External experimental comparison `f`. |
| `litcurate_4e5c23611ff4a5ae` (406), measured molar volume 24.46 cm3/mol | rejected citation row | External experimental ambient-volume comparison and no EOS coefficients. |

## Reproduction and result

`scripts/reproduce_kiefer_2002_fe_bridgmanite.py` verifies both serialized source triplets, evaluates pressure at `V/V0=0.8`, and performs inverse pressure-volume round trips. The computed checkpoints are analytical reproductions because the source does not publish its energy-volume grid.

- Accepted: **2** net-new records.
- Rejected pseudo-splits: **2**.
- Rejected citation rows: **5**.
- Zotero action: not performed; metadata above are import-ready.
