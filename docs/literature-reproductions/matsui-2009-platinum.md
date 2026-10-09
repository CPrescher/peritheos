# Matsui et al. (2009): Pt pressure scale with recovered electronic table

## Source and disposition

The supplied primary paper, Matsui et al., *The temperature-pressure-volume
equation of state of platinum*, J. Appl. Phys. 105, 013505 (2009),
[DOI 10.1063/1.3054331](https://doi.org/10.1063/1.3054331), was read and visually
checked on 2026-10-09: Equations (1)-(13), Tables I-IV, and the surrounding
experimental and fitting descriptions. The supplied PDF's SHA-256 is retained
in `peritheos/data/datasets/platinum-matsui-2009-source.json`.

`platinum_matsui_2009_vinet_300k` remains the canonical reference isotherm,
unchanged numerically. The complete pressure parameterization is registered as
`platinum_matsui_2009_vinet_mgd_electronic`, using all 51 Pt electronic nodes
from the supplied Tsuchiya-Kawamura paper. The erratum changes Au only.
**Linear interpolation between 100 K nodes is an explicit Peritheos numerical
choice**, since neither source specifies interpolation. The record reproduces
source nodes and supplies a qualified continuous implementation; the original
author optimization is still not reproduced. Its determination method is
`hybrid`: experimental shock/assessed expansion constraints plus theoretical
electronic pressures.

## Equations, normalization and parameter provenance

With r = V/V0, the printed Equations (5)-(10) specify

\[
P(V,T)=P_{300}(V)+\frac{\gamma(V)}{V_m}
       [E_D(V,T)-E_D(V,300)]+\Delta P_{el}(T),
\]

\[
P_{300}=3K_0\frac{1-r^{1/3}}{r^{2/3}}
          \exp[\tfrac32(K'_0-1)(1-r^{1/3})],\quad
\gamma=\gamma_0 r^q,\quad
\theta=\theta_0\exp[(\gamma_0-\gamma)/q],
\]

\[
E_D=\frac{9nRT}{(\theta/T)^3}
       \int_0^{\theta/T}\frac{z^3}{e^z-1}\,dz.
\]

Table II uses conventional fcc-cell V0 = 60.38 A^3, four Pt atoms per cell.
Equation (9) uses n = 1 and molar R: energy is J/mol of Pt atoms, so
`V_m = V_cell * N_A * 1e-30 / 4` in m^3/mol. Divide the energy/volume
pressure by 1e9 for GPa. Peritheos thermal kernels instead use J/bar/mol:
`V_m_native = V_cell * N_A * 1e-25 / 4`, with pressure converted from bar
to GPa by 1e4. The four-atom conversion applies to both V and V0.

| Quantity | Value | Source role | Printed parameter error |
|---|---:|---|---|
| V0 (A^3/cell) | 60.38 | adopted reference volume | absent |
| K0 (GPa) | 273 | optimized | absent |
| K0' | 5.20 | optimized | absent |
| gamma0 | 2.70 | optimized | absent |
| q | 1.10 | optimized, constant with volume | absent |
| theta0 (K) | 230 | fixed from earlier work | absent |
| Tr (K) | 300 | pressure reference | absent |
| n | 1 | Pt atoms per formula unit | not statistical |

No confidence convention, parameter covariance or electronic parameter errors
are supplied. Parentheses in Table II's Fei and Zha comparison columns do not
belong to the Matsui parameters. Matsui's integrated-q temperature law differs
from the printed variable-exponent Fei law; substituting the latter changes
the pressure. The integrated law obeys `-d ln(theta)/d ln(V) = gamma(V)`.

The raw [Tsuchiya and Kawamura (2002)](https://doi.org/10.1103/PhysRevB.66.094115)
Table I pressures are referenced to 0 K. All 51 nodes from 0 to 5000 K in
100 K steps are transcribed in
`tsuchiya-kawamura-2002-table1-electronic-pressure.csv`. Matsui instead uses
`Delta P_el(T) = P_el_table(T) - P_el_table(300)` with the raw 300 K value
0.04 GPa. Thus raw 0.25, 1.64 and 3.82 GPa give increments 0.21, 1.60 and
3.78 GPa at 1000, 3000 and 5000 K. Subtract the reference once; Matsui's
printed examples are already subtracted.

The supplied [2003 erratum](https://doi.org/10.1103/PhysRevB.67.019902)
corrects nine **Au** entries from 3100 through 3900 K to 0.13, 0.14, 0.15,
0.17, 0.18, 0.19, 0.21, 0.22 and 0.24 GPa. **Pt is unchanged.** Original
and corrected Au columns remain separate; no Au EOS is changed in this audit.
The apparently irregular printed Pt value 0.95 GPa at 2100 K is retained.
Both supplied PDFs were read and visually checked, including the equations,
Table I, Figures 4-5 and the complete erratum. Their hashes are retained in
the source manifest.

Tsuchiya-Kawamura Equations (1)-(6) define
`P_el = -d[F_el(V,T)-F_el(V,0)]/dV`, with `F_el=U_el-T*S_el` and
Fermi-Dirac DOS integration. The calculation uses full-relativistic FP-LMTO
in LSDA. Linear fits of electronic free energy versus volume yield the
adopted volume-independent pressure in Figure 4's volume range (approximately
7.5-9 cm^3/mol for Pt). Matsui adopts this approximation over its broader
model envelope. A numerical DOS, fitted analytic coefficients and a
sub-grid interpolation rule are not published. A Sommerfeld/T-squared
approximation is explicitly unsuitable at high temperature in this source.
The table is theoretical input, not measured pressure.

For `T_i <= T <= T_(i+1)`, Peritheos implements
`P_el(T) = P_i + (P_(i+1)-P_i)*(T-T_i)/(T_(i+1)-T_i)` and subtracts the
same interpolated pressure at Tr. Nodes are fixed serialized configuration,
excluded from scalar fit parameters. Temperatures outside 0-5000 K are
rejected; the thermal API requires positive T. The source's recommended
300-3000 K envelope remains separate from this numerical table domain.
Pressure/volume and bounded temperature inversions are available. Electronic
heat capacity, entropy and free energy are unsupported because this pressure
table does not recover the complete caloric potential.

## Independent numerical checks

Run `.venv/bin/python -m scripts.audit_matsui_2009_platinum --check` and
`.venv/bin/python -m pytest tests/test_matsui_2009_platinum.py`.
The independent SI implementation uses explicit Vinet and Debye quadrature;
the production EOS is only the comparator. The audit output is
[`matsui-2009-platinum-audit.json`](../data/matsui-2009-platinum-audit.json).

- Full pressure agrees below 5e-13 GPa with an independent SI implementation;
  cold pressure agrees below 5e-13 GPa and the lattice increment below
  4e-14 GPa over the seven Table III volume ratios and five temperatures.
- The seven printed 300 K states agree within 0.00545 GPa. This is a
  rounded-output check, not exact equality of all printed decimal places.
- All 35 Table III pressure states are reproduced within 0.007066 GPa using
  the printed coefficients and electronic table. The small mismatch exceeds
  half a 0.01 GPa digit at some states: hidden precision is not recovered.
- All 12 Table I calculated pressures are reproduced within 0.008583 GPa.
  Both checks use source electronic nodes and are independent of interpolation.
  Omitting electronics misses about 0.21 and 1.60 GPa at 1000 and 3000 K.
- Table I's printed observed-minus-calculated residuals give RMS
  0.09229 GPa and maximum magnitude 0.19 GPa. These are checks of the
  paper's printed validation table. Reconstructed predictions and observed-minus-
  reconstructed residuals are also retained separately in the audit.

For a worked state at V/V0 = 1 and 3000 K, the cold pressure is zero,
the lattice increment is 19.944840 GPa and the electronic increment is
1.60 GPa, giving 21.544840 GPa versus the printed 21.54 GPa.

Subtracting the independently calculated cold and lattice pressures from
Table III constrains electronic **derived output** at its grid temperatures.
The common intersections of +/-0.005 GPa pressure-rounding intervals are:

| T (K) | Electronic increment (GPa), conditional interval |
|---:|---:|
| 500 | 0.037383-0.037934 |
| 1000 | 0.206394-0.207969 |
| 2000 | 0.810471-0.811673 |
| 3000 | 1.599229-1.600160 |

These diagnostic bounds assume the rounded Table II coefficients are exact; they are
neither statistical confidence bounds nor replacements for the recovered
upstream electronic table.
Even the cold table exceeds half a last decimal at two states, so hidden
precision remains unresolved. The 1000 K interval is consistent with the
two-decimal electronic example's rounding interval, although fixing its
printed value to exactly 0.21 misses two pressures by more than 0.005 GPa.
At 300 K the reference subtraction is exactly zero, regardless of table
rounding residuals.

A proposed `A(T^2-300^2)` term would require A = 2.30769e-7, 1.79574e-7,
and 1.51746e-7 GPa/K^2 from the three printed examples. They disagree well
beyond two-decimal pressure rounding. A quadratic term cannot be credited
as Matsui's equation. The implementation uses the upstream Table I nodes;
sub-grid interpolation remains explicitly qualified.

## Recovered input types and original-fit limits

The material document links four recovered evidence datasets, excluded from
`fit_datasets`:

| Evidence | Rows/states | Interpretation |
|---|---:|---|
| Matsui Table I | 12 | Measured static V/V0 and T; Pobs calibrated with Matsui, Parker and Leslie (2000) MgO; Pcalc retained separately |
| Matsui Table III | 7 rows, 35 states | Calculated pressure output; Dewaele comparison values retained separately |
| Arblaster (1997) Table II, 100-2000 K | 32 | Assessed crystallographic expansion from selected dilatometry with vacancy corrections |
| Tsuchiya-Kawamura (2002) Table I | 51 | First-principles raw electronic pressure input, with original/corrected Au preserved separately |

Table I preserves both parenthetical volume/pressure errors, with confidence
unspecified. No temperature error, covariance or raw paired MgO volume is
provided. Its Pobs is not a pressure-scale-free measurement; the fitted Pt
model is developed independently of that MgO scale and Table I tests it.
The static envelope is 21.04-41.86 GPa and 300-1600 K. Twenty external Zha
states are discussed as further validation after MgO re-reduction; their
recalculated row pressures are not printed here and have not been recovered.

The primary [Arblaster assessment](https://doi.org/10.1595/003214097X4111221)
was downloaded from the publisher and Table II visually checked. It tabulates
atomic volume and a(T) on an ITS-90 basis; the audit normalizes the rounded
lattice parameters as `[a(T)/a(300)]^3` rather than silently identifying
Arblaster's ambient cell volume with Matsui's V0. Its
[2006 methodology and erratum](https://doi.org/10.1595/147106706X129088)
were also read and visually checked. They replace the high-T
dilatometry/vacancy combination with lattice-expansion polynomials and
correct notation/heat-capacity equations; the original 1997 table cited by
Matsui is preserved separately. It is an assessment of observations, not 32
new direct measurements. Matsui does not publish its exact selected states,
unrounded expansion values or normalization procedure.

The existing seven primary Holmes (1989) shots remain unchanged. Only TaPt5,
TaPt6 and TaPt8 are below 290 GPa. Their measured Us/Up, density and reduced
shock pressure are checked against the momentum and mass conservation
identities; the four higher-pressure shots are excluded from this audit's
Matsui subset. Matsui fits the combined Holmes/Marsh **velocity data** over
28-290 GPa, obtaining C = 3.604(19) km/s and S = 1.543(13). The confidence
convention is absent. Those reduced shock states have no measured T column.

Equations (1)-(4) give Us = C + S Up, PH = rho0 Us Up,
VH/V0 = 1-Up/Us, and EH-E0 = PH(V0-VH)/2. The pressure product is in GPa
when density is g/cm^3 and velocities are km/s; the energy relation requires
a consistent mass-specific or molar volume basis and SI conversion.
Matsui subsequently samples the **derived** linear-fit Hugoniot between
16 and 290 GPa, iterates Equations (5) and (11)-(13), and uses expansion
to constrain K0, K0', gamma0 and q with theta0 fixed. The lower limit of
16 GPa is distinct from the observed velocity fit's 28 GPa limit.

The Marsh rows, exact derived-Hugoniot sampling, weights, optimization
objective, convergence procedure and covariance remain unavailable. The
later Zhu workbook contains a broader compiled shock grid, but using it
as Matsui's original selected velocity data would lose source identity and
the initial-density/velocity provenance. It is not substituted here.
Fitting Table III would merely recover calculated output; fitting Table I
would replace the authors' shock/expansion objective with a different one.
The author's joint fit is therefore **not reproduced**, while cold and
phonon/electronic pressure equations and table-node reproduction are independently
verified. Linear interpolation is an implementation choice, not author-fit parity.

The proposed 0-300 GPa, 300-3000 K pressure-scale envelope is an extrapolated
source recommendation, separate from shock constraints, zero-pressure
expansion and static validation. It is not a measured P-T rectangle or a
claim that Pt remains solid at every point in that rectangle.
