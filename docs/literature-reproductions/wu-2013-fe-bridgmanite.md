# Wu, Zhao, and Tian (2013) Fe-bearing bridgmanite audit

## Primary authority

- D. Wu, J.-J. Zhao, and H. Tian, “Effect of substitution Fe2+ on physical properties of MgSiO3 perovskite at high temperature and high pressure,” *Acta Physica Sinica* **62** (2013), 049101, DOI [10.7498/aps.62.049101](https://doi.org/10.7498/aps.62.049101).
- Open primary [publisher PDF](https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.62.049101.pdf), SHA-256 `765f2c28d5c76c557d61048877b49624e3cd8e3a15440ae5a61900e5163104fd`.
- Checked locations: Section 2, Equation (10), Section 3.1, and Table 1.

The paper defines its third-order Birch-Murnaghan pressure equation explicitly and reports source-generated static coefficient triplets for Pbnm MgSiO3 and an ordered 20-atom high-spin (Mg0.75Fe0.25)SiO3 cell. Those are distinct composition branches. The paper also applies a quasiharmonic Debye calculation through 2000 K, but it does not publish a complete supported thermal-EOS coefficient set; no thermal pseudo-record is created.

## All four LitCurate rows

| Candidate | Disposition | Reason |
|---|---|---|
| `litcurate_1fd4855251635f0a` (832), MgSiO3 162.378/244/4.0 | accepted | Complete source-generated BM3. |
| `litcurate_eb3b372892561702` (833), Fe25 164.977/248/4.03 | accepted | Complete high-spin composition branch. |
| `litcurate_2edf2b54fb35f7a2` (834), experimental 162.345/264 | rejected citation row | Table 1 explicitly labels an external 300 K experimental comparison; incomplete. |
| `litcurate_8d4961bd7435eb6b` (835), literature 246–272/3.9 | rejected citation row | External literature range, not this calculation. |

`scripts/reproduce_wu_2013_fe_bridgmanite.py` verifies both Table 1 triplets, evaluates Equation (10) at `V/V0=0.8`, and round-trips through the inverse solver. The numerical P-V grid shown in figures is not printed and was not invented.

- Accepted: **2**.
- Rejected citation rows: **2**.
- Zotero action: not performed; citation metadata above are ready for import.
