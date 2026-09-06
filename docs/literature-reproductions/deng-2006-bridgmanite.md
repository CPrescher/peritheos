# Deng et al. (2006) MgSiO3 bridgmanite audit

## Primary source

- L.-W. Deng, J.-J. Zhao, G.-F. Ji, Z.-Z. Gong, and D.-Q. Wei, “First-Principles Study of Orthorhombic Perovskites MgSiO3 up to 120 GPa and Its Geophysical Implications,” *Chinese Physics Letters* **23** (2006), 2334–2337, DOI [10.1088/0256-307X/23/8/101](https://doi.org/10.1088/0256-307X/23/8/101).
- Open primary authority: [publisher record](https://cpl.iphy.ac.cn/article/id/40292) and [publisher PDF](http://cpl.iphy.ac.cn/en/article/pdf/preview/40292.pdf), SHA-256 `36401c07cb6603b949713c2660ca38ca249e6074bae8841f7a0297560fba4625`.
- Checked locations: method and EOS equation on page 2335; Table 1 and its caption on page 2336.

The paper reports static LDA calculations for a 20-atom (`Z=4`) Pbnm cell from 0 to 120 GPa. Table 1 contains two source fit branches. Its caption unambiguously says superscript `b` is the third-order Birch fit and superscript `a` is fourth order. The `b` branch is complete and executable as BM3. The `a` branch lists only `V0`, `K0`, and `K0'`; it omits the independent `K0''` needed to reproduce a fourth-order Birch equation and the numerical energy-volume grid is not published. It is therefore held rather than silently downgraded to BM3.

## All seven candidates

| LitCurate row | Disposition | Reason |
|---|---|---|
| `litcurate_e35eb6e4756f8473` (574), 162.253/264/3.83 | held incomplete BM4 | Primary Table 1 identifies this as fourth order, but provides no `K0''`; three numbers cannot reproduce BM4. |
| `litcurate_d61fdcb2a93f7b6f` (575), 161.658/255/4.00 | accepted | Complete third-order Birch branch `b`; `K0'=4` is fixed by the truncation. |
| `litcurate_1402a7cf80b58a03` (576), experiment ref. 1 | rejected citation row | External comparison and no K0. |
| `litcurate_dbecb6f0878cef2b` (577), experiment ref. 9 | rejected citation row | External comparison and no K0/K0'. |
| `litcurate_df4f6d2ea24654be` (578), experiment ref. 3 | rejected citation row | External Fiquet comparison and no K0. |
| `litcurate_fc6bb3f5f0b1a72c` (579), experiment ref. 22 | rejected citation row | External modulus-only comparison. |
| `litcurate_d97945919d0deabc` (580), experiment ref. 23 | rejected citation row | External modulus-only comparison. |

## Reproduction and result

`scripts/reproduce_deng_2006_bridgmanite.py` loads the serialized BM3, verifies the Table 1 triplet, evaluates pressure at `V/V0=0.8`, and round-trips the resulting pressure through the inverse solver. This is an analytical execution check; the source's 13 calculated structures are not numerically tabulated and were not fabricated from the plot.

- Accepted: **1**.
- Held incomplete source fit: **1**.
- Rejected citation rows: **5**.
- Zotero action: not performed; metadata above are ready for import.
