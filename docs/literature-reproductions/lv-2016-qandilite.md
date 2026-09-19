# Lv et al. (2016): synthetic qandilite at 300 K

Audit date: 2026-09-19. Two experimental records accepted for the measured
**Mg2.00(1)Ti1.00(1)O4** qandilite phase. No computational, empirical,
ultrasonic, or cited comparison EOS is imported.

## Primary authority and retrieval

Lv, M., Liu, X., Shieh, S. R., Xie, T., Wang, F., Prescher, C., and
Prakapenka, V. B. (2016), *Physics and Chemistry of Minerals* **43**, 301–306,
[doi:10.1007/s00269-015-0794-1](https://doi.org/10.1007/s00269-015-0794-1).
The six-page Springer publication of record was recovered from existing Zotero
item `268TSW4J`, attachment `3CBMV38U`. PDF SHA-256:
`2e61e8d2b32f573c3f992d762970527ac2bf25535113250fd60c6be0d31e46a5`.
The publisher landing page and Crossref record were checked: no linked
correction or official supplementary dataset was found. The article itself
has no electronic-supplement reference. Exact-DOI correction searching found
no relevant correction. The existing Zotero reference was read without
changing its collection membership.

The full paper, not its abstract or a coefficient database, supplies both
reference volumes. Table 1 and the equation were checked visually in the PDF.
The source PDF remains in Zotero; this repository bundles its factual table
transcription and checksum, not the copyrighted article.

## Candidate inventory and exact mapping

| Candidate | V0 (Å³/cell) | KT0 (GPa) | K′T0 | Disposition |
|---|---:|---:|---:|---|
| `qandilite_lv_2016_bm2` | 603.21(18) | 172(1) | 4 fixed | Accepted, default; the reduced fit plotted in Figure 3 and supported by the strain plot |
| `qandilite_lv_2016_bm3` | 603.10(25) | 175(5) | 3.5(7) | Accepted; independently meaningful free-derivative alternative with Figure 5 confidence ellipses |
| Figure 4 weighted F–f line | — | — | — | Diagnostic of these same fits; no separately tabulated complete EOS |
| Other Table 2 rows | — | — | — | Cited Fe2TiO4 fits, an earlier ultrasonic Mg2TiO4 modulus, and an empirical estimate; not source-owned fits |

The unnumbered equation on p. 303 is exactly

\[
P={3K_0\over2}(\eta^7-\eta^5)
[1+{3\over4}(K'_0-4)(\eta^2-1)],\qquad
\eta=(V_0/V)^{1/3}.
\]

The existing BM3 and BM2 implementations represent it directly. BM2 fixes
K′0=4 by model definition; it does not carry a redundant fitted parameter.
V0 and K0 are fitted in both records. Reference state: 0 GPa, 300 K;
conventional cubic cell with Z=8. No molar or primitive-cell conversion occurs.
The separately measured ambient value, **602.58(49) Å³**, is an observation,
not either fitted V0. Parentheses denote one standard deviation (Table 1/2
footnotes; Figure 5 caption). Missing numerical covariance is retained as
missing; the illustrated negative correlation is not converted into an
invented covariance matrix.

## Composition, crystallography, and actual scope

The major phase has ten EMPA analyses with the composition above (p. 302).
The synthetic product contains less than approximately 5% geikielite,
Mg1.00(2)Ti1.00(1)O3. This card describes the separately indexed qandilite
phase alone. Synthesis: 1673 K, 52 h; cooling at 5 K/min. Seventeen DAC
patterns show no decomposition, new peaks, splitting, or apparent broadening
through 14.9 GPa. The retained cubic phase is potentially metastable; neither
this pressure interval nor the single 300 K isotherm establishes equilibrium
phase stability or a thermal EOS.

The diffraction structure is separately sourced from O’Neill, Redfern,
Kesson, and Short (2003), *American Mineralogist* **88**, 860–865,
[doi:10.2138/am-2003-5-615](https://doi.org/10.2138/am-2003-5-615),
[primary PDF](https://rruff.net/doclib/am/vol88/AM88_860.pdf), SHA-256
`e488f6205e8ca980d3598154215969fbcaf1cf85a1263fe1029b9a042611671d`.
Table 1 first column (Si standard present, 111 included) gives
**a=8.44192(5) Å, uO=0.2594(1)**. Table 3 gives Fd-3m; pp. 861–862
establish the fully inverse room-temperature distribution, x=1.000(3), and
reversal of high-temperature disorder during cooling. The conventional
origin-2 Wyckoff realization is Mg on 8a (1/8,1/8,1/8), equal Mg/Ti on
16d (1/2,1/2,1/2), and O on 32e (u,u,u). Multiplicity × occupancy gives
Mg16Ti8O32, hence Z=8. This is an explicitly documented same-composition
structural fallback, not a claim to have refined Lv's site occupancies.
Its independently measured lattice is retained without rescaling it to an
EOS reference volume. The approximately 4% disorder mentioned by Lv at the
synthesis temperature is not asserted to persist at room temperature.

## Observations and calibration

All 18 Table 1 rows (p. 303) are in
`peritheos/data/datasets/qandilite-lv-2016-table1-pv.csv`,
SHA-256 `4a21ad66ea1fa6dffbd2f0b92d5456023165d3540575815b9f65933bd7fb5f40`.
They preserve original pressure, lattice parameter, conventional cell volume,
and every printed uncertainty. The ambient pressure is 0.0001 GPa; its
pressure uncertainty is unreported. The 17 DAC pressures are averages of
before/after fluorescence readings with assumed 0.1 GPa uncertainty.

The calibration is explicitly **Mao et al. (1978)** ruby fluorescence,
linked to `ruby_mao_1978`; the medium is neon. Neither the individual readings
nor row-wise ruby R1/reference wavelengths are published. Pressure-scale
re-reduction is therefore unavailable even though the cited scale is
executable. No calibrant wavelengths or neon volumes are synthesized.

## Numerical reproduction and independent fitting

Run `uv run python scripts/reproduce_lv_2016_qandilite.py` to regenerate
[`lv-2016-qandilite.json`](../data/lv-2016-qandilite.json).
The standalone evaluator implements the printed equation directly, independent
of Peritheos. It agrees with native evaluation to better than 2e-12 GPa.
The following independent measured checkpoint is Table 1's final state,
**14.9(1) GPa, 559.91(34) Å³**:

| Published model | P(559.91) (GPa) | V(14.9) (Å³) | All-row pressure RMSE (GPa) |
|---|---:|---:|---:|
| BM2 | 14.871946 | 559.841558 | 0.128393 |
| BM3 | 14.801423 | 559.664781 | 0.127302 |

Both calculated volumes agree within the measured 0.34 Å³ one-sigma error.
A pressure tolerance of 0.18 GPa includes the propagated volume error and the
0.1 GPa pressure error. Maximum all-row residuals are 0.262 and 0.274 GPa,
respectively; those residuals are disclosed rather than requiring each noisy
observation to match a smooth fit exactly.

The paper specifies least squares and EoS fit 5.2 (Figure 5), but does not
provide the regression input, exact numerical weights, residual direction,
ambient-row treatment, or covariance matrix. Thus an exact replay of its
regression and confidence ellipses is unavailable. The registered independent
check simultaneously fits all 18 printed rows using equal-weight pressure
residuals, retaining the source model order. Starting values are 600 Å³,
180 GPa, and (for BM3) K′0=4, not the reported optimum. No staged sequence is
stated or invented.

| Audit selection/objective | BM2 (V0, K0) | BM3 (V0, K0, K′0) |
|---|---|---|
| All 18, unweighted pressure | 603.192506, 172.110517 | 603.069635, 175.511883, 3.494719 |
| All 18, weighted volume | 603.194609, 172.187305 | 603.086351, 175.241780, 3.538544 |
| All 18, effective variance | 603.155954, 172.402969 | 603.027904, 176.212347, 3.421501 |
| 17 DAC only, unweighted pressure | 603.313380, 171.401637 | 603.270493, 172.421887, 3.852837 |

The weighted-volume sensitivity minimizes Σ[(Vmodel(P)-V)/σV]². Effective
variance minimizes Σ[(Pmodel(V)-P)²/(σP²+(dP/dV σV)²)]; only this sensitivity
assumes an exact ambient pressure, explicitly separate from the missing
source error. Both weighted alternatives include all 18 rows. The DAC-only
sensitivity omits only the independent ambient measurement. These are checks
of unresolved source choices, not alternative published coefficients.

Every coefficient from the registered all-row fit agrees within its published
one-sigma uncertainty. The generic primary-refit validator independently fits
the same selection through Peritheos. Published parameters remain untouched;
no opt-in refit record is warranted by these small diagnostic differences.
