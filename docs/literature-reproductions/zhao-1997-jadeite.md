# Zhao et al. (1997) jadeite thermal EOS audit

## Source and experimental reduction

The primary source is Y. Zhao, R. B. Von Dreele, T. J. Shankland, D. J.
Weidner, J. Zhang, Y. Wang, and T. Gasparik, “Thermoelastic equation of state
of jadeite NaAlSi2O6: An energy-dispersive Reitveld Refinement Study of low
symmetry and multiple phases diffraction,” *Geophysical Research Letters*
**24**, 5–8 (1997),
[doi:10.1029/96GL03769](https://doi.org/10.1029/96GL03769). The audit used the
official four-page publisher article, especially Experimental Aspects and
Table 1 on journal pages 5–6, Equation (1) on page 6, and the Discussion and
Table 2 on page 7.

The experiment used a DIA-6 multi-anvil press (SAM-85) at NSLS beamline X-17B.
Energy-dispersive spectra were collected with a Ge detector at fixed
`2theta=5.847 degrees`; jadeite and hBN peak positions and lattice parameters
were refined simultaneously with a modified GSAS whole-pattern Rietveld
method. Temperature came from a Pt–Pt/10%Rh thermocouple. Pressure came from
NaCl in the other half of the hBN chamber using Decker's 1971 NaCl EOS,
[doi:10.1063/1.1660714](https://doi.org/10.1063/1.1660714).

Cold compression produced significant deviatoric stress. The authors therefore
compressed to the maximum pressure, heated to about 800–1000 K or above to
remove it, and collected the reported cooling/decompression paths after the
sample became pseudo-hydrostatic. Table 1 is already the source-selected set:
all 31 of its rows are included in the EOS fit, and no bundled row is excluded.
The rejected cold-compression measurements are not tabulated, so they cannot be
reconstructed or counted independently.

The exact pressure scale is now identified in the material record, but the
paper does not print the simultaneous NaCl lattice parameters. Consequently the
reported pressures cannot be recalculated from the calibrant. It gives only an
aggregate pressure precision of about 0.08–0.12 GPa, not row-wise pressure
errors.

## Source equation and preferred fit

Equation (1) is a high-temperature third-order Birch–Murnaghan relation:

\[
P=3K_Tf(1+2f)^{5/2}\left[1-\frac{3}{2}(4-K'_0)f\right],
\qquad
f=\frac{1}{2}\left[\left(\frac{V_T}{V(P,T)}\right)^{2/3}-1\right],
\]

with

\[
K_T=K_{T0}+\dot K(T-300),
\qquad
V_T=V_0\exp\left[\int_{300}^{T}\alpha(0,T')\,dT'\right],
\qquad
\alpha(0,T)=a+bT.
\]

The source writes the more general expansion as `a + bT - c/T^2`, but explicitly
neglects the `c/T^2` term for this data range. It also neglects higher-order
temperature dependences of `K`, `K'`, and the mixed derivative. Peritheos's
`AlphaKT` model with `thermal_expansion_law="linear_temperature"` and
`reference_volume_law="integrated_expansivity"` is therefore the source
relation exactly, with `alpha0=a`, `alpha1=b`, and `dK_dT=dot(K)`.

For every tradeoff calculation Zhao et al. fixed `V0=403 A^3`. This is a model
coefficient, not the independently printed ambient Table 1 observation
`403.32(8) A^3`; the earlier Peritheos audit had incorrectly substituted the
latter. The preferred fit also fixes `K0'=5` and varies four coefficients:

| Quantity | Stored value | Published error | Fit role |
|---|---:|---:|---|
| `V0` (A^3) | 403 | not reported | fixed |
| `K0` (GPa) | 124.5 | 4.0 | varied |
| `K0'` | 5.0 | not applicable | fixed |
| `dK_dT` (GPa/K) | -0.0165 | 0.0049 | varied |
| `alpha0 = a` (K^-1) | 2.56e-5 | 2.2e-6 | varied |
| `alpha1 = b` (K^-2) | 2.6e-9 | 1.8e-9 | varied |

The prose gives the unrounded `K0=124.5(4.0) GPa`; Table 2 rounds it to
`125(4) GPa`. The complete tradeoff table is retained here because it is the
source's warning against interpreting the preferred column as a unique
unconstrained solution:

| Constraint | `a` (1e-5 K^-1) | `b` (1e-8 K^-2) | `K0` (GPa) | `-dK_dT` (1e-2 GPa/K) |
|---|---:|---:|---:|---:|
| `K0'=4` | 2.51(23) | 0.24(19) | 127(5) | 1.40(51) |
| `K0'=5` (preferred) | 2.56(22) | 0.26(18) | 125(4) | 1.65(49) |
| `K0'=6` | 2.62(21) | 0.29(18) | 122(4) | 1.89(48) |
| `K0'=9.7(9)` fitted | 2.91(17) | 0.25(15) | 114(5) | 2.53(37) |

The paper says its pressure range is insufficient to resolve `K0'` confidently
and prefers a constrained value compatible with earlier measurements.

## Refit and remaining ambiguity

The source does not state whether the EOS regression minimized pressure,
volume, or orthogonal residuals. It also omits numerical row weights,
temperature errors, parameter covariance, covariance scaling, and the
confidence convention for Table 2's parenthetical errors. Those details cannot
be recovered from the article and are not inferred.

The committed primary-refit ledger therefore uses the narrowest reproducible
choice: ordinary pressure residuals for all 31 Table 1 rows, with `V0=403 A^3`
and `K0'=5` held fixed. The four parameters span eleven orders of magnitude, so
the validator optimizes dimensionless scaled coordinates before transforming
the covariance back to physical units.

| Varied coefficient | Published | Peritheos refit |
|---|---:|---:|
| `K0` (GPa) | 124.5(4.0) | 123.7653(1.3079) |
| `alpha0` (K^-1) | 2.56(22)e-5 | 2.69078(22)e-5 |
| `alpha1` (K^-2) | 2.6(1.8)e-9 | 1.412(3.336)e-9 |
| `dK_dT` (GPa/K) | -0.0165(49) | -0.0152745(2798) |

The published surface has a pressure RMSE of `0.11001 GPa` on the rounded
table; the refit lowers it to `0.10058 GPa`. Every varied coefficient is within
its published uncertainty. The ledger nevertheless labels the result
`similar`, not strict `parity`, because `alpha1` shifts by about 46% relative to
its very small central value. This is consistent with the source's broad
`alpha1` error and strong tradeoffs, but the missing objective and covariance
prevent a claim of exact protocol reproduction.
