# Suzuki (2016) low-pressure epsilon-FeOOH thermal EOS

## Scope and primary evidence

This record audits Akio Suzuki (2016), *Pressure-volume-temperature equation
of state of epsilon-FeOOH to 11 GPa and 700 K*, *Journal of Mineralogical and
Petrological Sciences* 111, 420--424,
[doi:10.2465/jmps.160719c](https://doi.org/10.2465/jmps.160719c). The final
five-page J-STAGE PDF was checked directly; its SHA-256 is
`0066dfbd7b3ab45ce7aae866aef24978c63571a433651ae470912cd0cc8b0289`.
J-STAGE lists no supplementary files.

The executable record is `e_feooh_suzuki_2016_bm3_thermal_2` in the existing
`e_feooh` material. It represents InOOH-type epsilon-FeOOH over the observed
0.0001--11.07 GPa and 300--700 K range. It does not represent the natural
goethite starting material, which was transformed to epsilon-FeOOH before the
measurements. It is also deliberately separate from
`e_feooh_hc_low_spin_thompson_2017_bm3_1`, the static hydrogen-centered Pnnm
low-spin branch fitted over 30--140 GPa, and from the approximately 53 GPa
spin transition discussed by Gleason et al. (2013).

The existing diffraction model is the hydrogen-off-center P21nm (#31)
epsilon-FeOOH structure with two FeOOH formula units per conventional cell.
Its coordinates come from the separately cited 20 GPa static structure of
Insixiengmay and Stixrude (2023); they are not treated as Suzuki observations.
Suzuki's Table 1 volumes are likewise conventional orthorhombic-cell volumes:
the ambient lattice product
`4.9544 * 4.4594 * 2.9999 = 66.2779 A^3` agrees with the printed
`66.278(6) A^3`.

## Exact equation and parameters

Equation 1 is the standard third-order Birch-Murnaghan pressure equation,
evaluated at each temperature:

\[
P(V,T)=\frac{3K_{T0}(T)}{2}
\left[\left(\frac{V_{0T}}{V}\right)^{7/3}
-\left(\frac{V_{0T}}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_T-4)
\left[\left(\frac{V_{0T}}{V}\right)^{2/3}-1\right]\right\}.
\]

The three following displayed relations define

\[
K_{T0}(T)=K_{T0}(300)
+\left(\frac{\partial K_T}{\partial T}\right)_P(T-300),
\]

\[
V_{0T}=V_{0,300}\exp\left[\int_{300}^{T}\alpha(T')dT'\right],
\qquad
\alpha(T)=a_0+a_1(T-300).
\]

Therefore the exact integrated reference volume is

\[
V_{0T}=V_{0,300}\exp\left[a_0(T-300)
+\frac{a_1}{2}(T-300)^2\right].
\]

This last convention matters. Peritheos previously supported a linear law
whose intercept is at absolute zero, `alpha0 + alpha1*T`. The new
`linear_reference_temperature` configuration keeps Suzuki's `alpha0` at
`Tr = 300 K` and stores the published coefficients without an opaque
intercept transformation.

| Parameter | Published and stored value | Fit status |
|---|---:|---|
| `V0` | 66.278(6) A^3 | fixed |
| `K0` | 135(3) GPa | fitted |
| `K0_prime` | 6.1(9) | fitted; temperature-independent |
| `dK_dT` | -0.05(2) GPa/K | fitted |
| `alpha0` | 2.6(7)e-5 K^-1 | fitted; defined at 300 K |
| `alpha1` | 1.0(3)e-7 K^-2 | fitted |
| `Tr` | 300 K | fixed |

The Table 2 footnote says parenthetical values are standard deviations. The
paper does not provide the parameter covariance, numerical weights, scalar
fit statistic, or EosFit7c version. It explicitly notes that `K0` and
`K0_prime` are inversely correlated.

## Primary observations and pressure scale

All 33 rows of article Table 1 are transcribed in
`epsilon-feooh-suzuki-2016-table1-pvt.csv`. The dataset retains pressure,
temperature, three lattice parameters, conventional-cell volume, source row
order, and every printed parenthetical one-standard-deviation uncertainty.
The ambient `0.0001 GPa` value and all temperatures have no printed errors.
No figure digitization was needed.

The experiment used NaCl as the pressure standard and explicitly cites Brown
(1999), [doi:10.1063/1.371596](https://doi.org/10.1063/1.371596). Table 1 does
not report the row-wise NaCl lattice parameters, so pressures cannot be
recalculated observation by observation. The source page states copyright
2016 Japan Association of Mineralogical Sciences and gives no open data reuse
license. The factual table transcription is retained for scientific
validation, with that licensing limitation recorded rather than asserting an
open license.

## Numerical reproduction and independent refit

A direct implementation of the printed equations agrees with the native
Peritheos model to machine precision. Evaluated at all 33 printed Table 1
volumes and temperatures, the published parameters give a pressure RMSE of
`0.1020197851 GPa` and a maximum absolute residual of `0.2432053548 GPa`.
At the high-temperature endpoint (`62.63 A^3`, `700 K`), the curve gives
`10.8335628033 GPa`, compared with the printed `11.07(2) GPa` observation.

The independent Peritheos fit uses the 32 finite-pressure rows with their
reported pressure and volume standard deviations. The ambient row fixes
`V0,300`; it has no printed pressure uncertainty and is not added as an
invented weighted residual. Temperatures are not adjusted because the source
reports no temperature errors. The joint pressure-volume errors-in-variables
result is:

| Parameter | Published | Peritheos refit | Difference / published sigma |
|---|---:|---:|---:|
| `K0` (GPa) | 135 | 133.81657 | -0.394 |
| `K0_prime` | 6.1 | 6.39080 | +0.323 |
| `alpha0` (K^-1) | 2.6e-5 | 2.70888e-5 | +0.156 |
| `alpha1` (K^-2) | 1.0e-7 | 1.04091e-7 | +0.136 |
| `dK_dT` (GPa/K) | -0.05 | -0.0508147 | -0.041 |

The refit has reduced chi-square `0.49305845` and a raw-coordinate pressure
RMSE of `0.09607800 GPa`. Every coefficient is within the corresponding
published one-standard-deviation interval. Because the source's precise
EosFit weighting, unrounded inputs, and covariance are unavailable, this is a
parity diagnostic rather than a new replacement EOS record; the published
parameterization remains the only Suzuki record.

Run the reproduction from the repository root with:

```bash
uv run python scripts/reproduce_suzuki_2016_epsilon_feooh.py
```

Suzuki also notes a possible change in compression slope between 3 and 7 GPa,
attributing it either to deviatoric stress or hydrogen-bond symmetrization.
Peritheos preserves the source's single smooth HTBM fit and records this
caveat; it does not infer a transition model or extend the fitted surface
beyond the observations.
