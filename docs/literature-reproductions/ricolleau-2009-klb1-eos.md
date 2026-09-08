# Ricolleau et al. (2009) KLB-1 EOS audit

Audit date: 2026-09-08.

## Primary evidence and redistribution boundary

The authoritative source is Ricolleau et al., “Density profile of pyrolite
under the lower mantle conditions,” *Geophysical Research Letters* **36**,
L06302 (2009), [doi:10.1029/2008GL036759](https://doi.org/10.1029/2008GL036759).
The publisher identifies `grl25568-sup-0002-ts01.txt` as Table S1,
“Experimental conditions and unit cell volumes of phases.” Its accompanying
readme says that pressure is calculated from the Au volume with Fei et al.
(2007), that Mg-perovskite, ferropericlase, Ca-perovskite, pressure-medium, and
Au volumes come from diffraction refinements, and that parenthetical values are
one sigma.

The repository retains the factual table transcription, the small amount of
context needed to interpret its run labels and footnotes, and complete source
attribution. It does not redistribute the article text or figures. The
publisher download was 15,129 bytes with SHA-256
`2900b2d436207ae4a62289b8fe9464a714f1a538fb7dc08f76651b7bd5d9d4d9`;
the repository text normalizes line endings and has SHA-256
`e3f504dcb8953dbd6ef0cc0d341090a05116d32c865b748b9c95173fb94e7aa1`.

The raw TSV preserves all labels, footnotes, values, and the eight continuation
lines carrying a second pressure-medium volume. The normalized CSV contains all
153 numerical observations: 17 at 300 K and 136 heated observations. It expands
parenthetical uncertainties by printed decimal precision and never fills an
unprinted value.

## Published protocol that can be reconstructed

Section 3 and Table 1 specify a high-temperature Birch–Murnaghan EOS with
`K0_prime=4`, hence BM2, and

`K(T) = K0 + dK_dT * (T - 300 K)`

with `alpha(T) = alpha0 + alpha1*T`. Peritheos integrates that expansivity as

`V0(T) = V0 * exp[alpha0*(T-300) + 0.5*alpha1*(T^2-300^2)]`.

The article fixes Mg-perovskite and Ca-perovskite `V0`, fits their remaining
thermal-BM2 coefficients, fixes ferropericlase high-spin and low-spin `K0` to
158 and 170 GPa respectively, obtains separate 300 K limiting-branch `V0`
values, and fits heated ferropericlase from the high-spin reference branch.

The paper does **not** report an exact row mask, residual objective, weighting,
temperature uncertainties, parameter bounds, or covariance scaling. The audit
therefore uses a single transparent reconstruction: orthogonal-distance
regression with the printed one-sigma pressure and phase-volume errors,
temperature held fixed, and unscaled covariance. This is source-constrained,
but it is not represented as the authors’ undocumented implementation.

Table S1 contains Au volume and temperature for every row, and the exact cited
Fei et al. (2007) Au thermal EOS is executable as `gold_fei_2007_vinet_2`.
Re-reducing all 153 rows gives recalculated-minus-published pressure residuals
with a +0.669 GPa mean, 0.701 GPa RMS, and 0.941 GPa maximum absolute value.
Thus recalculation is operationally ready and closely reproduces the scale, but
the printed pressures are not claimed to be bit-for-bit regenerated. No hidden
correction or alternate Fei implementation is inferred.

## Ferropericlase spin selection

The seven 300 K points below 50 GPa form the limiting high-spin branch. Three
300 K points at 52.22, 53.00, and 54.89 GPa occupy the observed volume collapse
and remain in the dataset without assignment to either branch. The next seven
300 K points, from 74.21 to 98.36 GPa, form the separated limiting low-spin
branch. The audit’s 70 GPa selection boundary merely names the gap in the
observed pressure sequence; it is not claimed as a published transition
cutoff.

This selection is strongly checked by the published result. With `K0` fixed,
the high-spin fit gives `V0 = 76.43828 ± 0.01771 A^3` versus `76.44(2) A^3`, and
the low-spin fit gives `V0 = 74.04288 ± 0.02538 A^3` versus `74.04(2) A^3`.
Both central values and their unscaled errors reproduce the printed rounding.

## Thermal refit results

| Phase | Coefficient | Published | Refit ± unscaled standard error |
|---|---|---:|---:|
| Mg-perovskite | K0 (GPa) | 245(1) | 245.1917 ± 0.2242 |
|  | dK/dT (GPa/K) | -0.036(1) | -0.036371 ± 0.000981 |
|  | alpha0 (K^-1) | 3.19(17)e-5 | 3.4591e-5 ± 1.3225e-6 |
|  | alpha1 (K^-2) | 0.88(16)e-8 | 6.5642e-9 ± 1.2051e-9 |
| Ca-perovskite | K0 (GPa) | 244(1) | 243.6059 ± 1.0851 |
|  | dK/dT (GPa/K) | -0.035(2) | -0.035130 ± 0.001394 |
|  | alpha0 (K^-1) | 3.06(19)e-5 | 3.4593e-5 ± 1.7573e-6 |
|  | alpha1 (K^-2) | 0.87(18)e-8 | 5.6142e-9 ± 1.5534e-9 |
| High-spin ferropericlase, heated rows | dK/dT (GPa/K) | -0.034(1) | -0.034019 ± 0.000707 |
|  | alpha0 (K^-1) | 2.20(20)e-5 | 2.7597e-5 ± 1.6322e-6 |
|  | alpha1 (K^-2) | 3.61(27)e-8 | 3.0875e-8 ± 1.9537e-9 |

Every Mg- and Ca-perovskite coefficient agrees within combined two sigma. The
ferropericlase `alpha0` reconstruction lies just outside combined two sigma but
within the ledger’s 30% thermal-coefficient similarity threshold; the other
ferropericlase thermal coefficients agree within combined two sigma. Vertical
pressure RMSEs are 1.306 GPa for Mg-perovskite, 1.710 GPa for Ca-perovskite,
and 1.479 GPa for heated ferropericlase.

Exact thermal-fit parity is irreducibly blocked by the source’s missing fit
mask, objective/weights, temperature-error treatment, and covariance
convention. No alternative exclusions, weights, or uncertainties are inferred
to close the remaining differences.

## Reproduction

Run `python scripts/audit_ricolleau_2009_klb1_eos.py --write` to regenerate the
normalized CSV and structured audit. The four records are integrated into the
common primary-EOS ledger by `scripts/validate_primary_eos_refits.py`.
