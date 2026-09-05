# Matsui et al. (2012) ferropericlase audit

## Source

- Primary paper: Masanori Matsui et al., “Static compression of (Mg0.83,Fe0.17)O and (Mg0.75,Fe0.25)O ferropericlase up to 58 GPa at 300, 700, and 1100 K,” *American Mineralogist* 97, 176–183 (2012), [doi:10.2138/am.2012.3937](https://doi.org/10.2138/am.2012.3937).
- Authoritative archival PDF: [RRUFF American Mineralogist archive](https://www.rruff.net/odr/view/downloadfile/81129), SHA-256 `c1f4465e866cd1cef5bb5ecff1704aedb89cb6e44404058d965314385450b24f` (audited 2026-09-05).
- Pressure scale: Masanori Matsui, “High temperature and high pressure equation of state of gold,” *Journal of Physics: Conference Series* 215, 012197 (2010), [doi:10.1088/1742-6596/215/1/012197](https://doi.org/10.1088/1742-6596/215/1/012197).

## Primary-source result

Equations (1)–(6) define a third-order Birch–Murnaghan 300 K reference isotherm plus a Mie–Grüneisen–Debye thermal-pressure increment:

`P(V,T) = P_BM3(V,300 K) + gamma(V)/V * [E_D(V,T) - E_D(V,300 K)]`.

The subtraction at 300 K is explicit and therefore maps to `thermal_pressure_reference = "reference_temperature"`; interpreting the thermal term as absolute Debye pressure would change the published reference state. The source uses `gamma = gamma0 (V/V0)^q`, the integrated Grüneisen Debye-temperature law `theta = theta0 exp[(gamma0-gamma)/q]`, and `n=2` atoms per formula unit.

The two synthesized, stoichiometric samples are cubic B1 rocksalt (`Fm-3m`, conventional `Z=4`). Table 3 reports:

| Composition | V0 (A^3/cell) | K0T (GPa) | K0T' | theta0 (K) | gamma0 | q |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| (Mg0.83Fe0.17)O | 75.849(11) | 160 fixed | 4.08(2) | 500 fixed | 1.53(4) | 0.7(2) |
| (Mg0.75Fe0.25)O | 76.372(16) | 160 fixed | 4.22(3) | 500 fixed | 1.64(4) | 0.7(2) |

The 300 K BM3 parameters were fitted first; `gamma0` and `q` were then least-squares fitted to the 700 and 1100 K measurements. Only measurements below 47 GPa were used for the high-spin equations. The source does not report fit weights, parameter covariance, or an uncertainty confidence convention.

Tables 1 and 2 are complete primary P–V–T tables and require no digitization. They also contain Au `V/V0`, observed-pressure errors, calculated pressures, and observed-minus-calculated residuals. Peritheos retains every row, while `fit_included` marks the strict below-47-GPa selection. The excluded rows document the growing volume anomaly through the high-spin to low-spin crossover rather than silently extending the EOS through it.

## LitCurate disposition: all seven same-DOI candidates

| LitCurate identifier | Candidate | Disposition | Evidence and mapping |
| --- | --- | --- | --- |
| `litcurate_1306566ad2bcaa9c` | (Mg0.83Fe0.17)O | **Accepted** | Source-reported thermal EOS in Equations (1)–(6) and Table 3. Implemented as `mg083fe017o_matsui_2012_bm3_mgd_1`; all 34 Table 1 rows are bundled, 23 below-47-GPa rows are fitted. |
| `litcurate_51fedd92173da0ac` | (Mg0.75Fe0.25)O | **Accepted** | Source-reported thermal EOS in Equations (1)–(6) and Table 3. Implemented as `mg075fe025o_matsui_2012_bm3_mgd_1`; all 39 Table 2 rows are bundled, 30 below-47-GPa rows are fitted. |
| `litcurate_e1ea89b341cbad7b` | (Mg0.64Fe0.36)O, van Westrenen et al. (2005) | **Hold / citation trace** | Matsui et al. quote the earlier fit (`V0=77.369(8) A^3`, `K0T=159(1) GPa`, `K0T'=3.8(1)`, `theta0=576 K`, `gamma0=1.63(2)`, `q=0.54(3)`) but did not originate it. Audit and reproduce it from the underlying primary paper, [doi:10.1016/j.pepi.2005.03.003](https://doi.org/10.1016/j.pepi.2005.03.003), not from this secondary citation. |
| `litcurate_103b130359c9b477` | (Mg0.64Fe0.36)O recalibrated to Matsui (2010) Au | **Hold / not reproducible from this paper alone** | Page 182 reports a source recalculation with `V0=77.369(8) A^3` fixed, `K0T=160 GPa` fixed, and BM3 `K0T'=4.35(7)`. It combines Fei et al. (1992) and van Westrenen et al. (2005) 300 K observations after pressure recalibration, but publishes neither the recalculated row-wise pressures nor its fit selection/weights. A production record should wait until both underlying tables and the Matsui (2010) Au calibration are reconstructed. |
| `litcurate_e429c9c80d6bff5e` | MgO, Jackson and Niesler (1982) | **Rejected from EOS queue** | The paper cites the ambient adiabatic modulus `K0S=162.5(2) GPa` as an elasticity comparison. No `V0`, `K0'`, or EOS fit is reported here. |
| `litcurate_ae5b5611e53f2a1b` | MgO, Chang and Barsch (1969) and related literature | **Rejected from EOS queue** | The cited `K0T=159.7 GPa` is one ambient elastic-modulus value used to justify fixing `K0T=160 GPa`. It is not a Matsui et al. EOS fit and has no complete EOS coefficient set in this source. |
| `litcurate_699aa6b747fdb9b7` | FeO, Jackson et al. (1990) | **Rejected from EOS queue** | The discussion compares an approximately 174 GPa ambient adiabatic modulus for corrected “stoichiometric” FeO. It is an elasticity citation, not an EOS fit, and the LitCurate row has no usable coefficient set. |

This resolves the two directly actionable candidates, preserves the two Fe36 candidates as explicitly scoped follow-up work, and removes the three incomplete elastic-modulus snippets from the EOS implementation queue.

## Reproduction

Run:

```bash
python scripts/reproduce_matsui_2012_ferropericlase.py
```

Using the rounded published parameters, the implemented curves reproduce the paper’s `Pcalc` column over all rows to 0.013918 GPa maximum for Fe17 and 0.005240 GPa maximum for Fe25. Against the included observed pressures, the RMSE values are 0.123303 and 0.142610 GPa, respectively. If the excluded high-pressure observations are included in that diagnostic, residuals grow systematically because those states cross into the mixed-spin regime; they are data provenance, not part of the high-spin fit.

An exact independent coefficient refit is intentionally not claimed: the paper states least squares but omits its weighting and covariance. The published-curve-to-published-table comparison is sufficient to verify equation order, cell-volume basis, the 300 K subtraction convention, and correct coefficient transcription.

The Au observations are preserved for future pressure-scale recalculation. That recalculation is currently marked `reference_eos_not_bundled` because the cited Matsui (2010) Au thermal EOS is not yet an executable Peritheos calibration record.
