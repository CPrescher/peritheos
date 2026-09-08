# Anderson et al. (1989) gold thermal EOS audit

## Outcome

The Anderson, Isaak, and Yamamoto gold EOS is reproducible as a published
parameterization, but **not as a global coefficient fit**. The paper does not
fit a common set of pressure-volume-temperature observations. It constructs
Equation (29) in stages from heterogeneous literature properties, derived
coefficient tables, graphical smoothing, separate linear regressions, and four
thermodynamic-consistency trials.

Peritheos therefore keeps `gold_anderson_1989_bm3_1` as `not_refittable`, while
bundling every numerical table that constrains or verifies the executable
equation. Fitting the 126 Table V values would be circular: they are calculated
outputs of Equation (29), not observations.

Primary source: Anderson, Isaak, and Yamamoto, *Anharmonicity and the equation
of state for gold*, Journal of Applied Physics **65**, 1534-1543 (1989),
[doi:10.1063/1.342969](https://doi.org/10.1063/1.342969).

## What the paper actually optimized

The source protocol is a sequence rather than one regression:

| Stage | Numerical input | Operation | EOS consequence |
|---|---|---|---|
| Heat capacity | Hultgren et al. below 300 K; Barin and Knacke Equation (4) above 300 K | interpolation plus analytical evaluation | enters \(\gamma=\alpha K_S/(\rho C_p)\) |
| Expansivity | White and Collins, AIP Handbook, and Touloukian et al. | hand-smoothed below 250 K; quadratic least squares from 250-1000 K | supplies \(\alpha(T)\) |
| Elasticity | Neighbours and Alers through 300 K; Chang and Himmel from 300-550 K | separate linear fits to ultrasonic elastic modes | supplies \(K_S(T)\), \(K_T(T)\), and their slopes |
| Constant-volume modulus | Table III and Equation (8) | numerical integration for four trial \(K'_T\) values: 6.39, 6.12, 5.5, and 5.21 | gives four \((\partial K_T/\partial T)_V\) slopes |
| Consistency selection | Table IV and Equation (18) | compare two sides of a thermodynamic identity | favors \(K'_T\) near 5.2-5.5; no scalar objective is reported |
| Thermal slope | Table IV \(P_{TH}\) above \(1.3\theta_D\) | linear least squares | \(\alpha K_T=7.14\times10^{-3}\) GPa K\(^{-1}\) |
| Final static curve | revised Heinz-Jeanloz EOS | adopt \(K_{T0}=166.65\) GPa and \(K'_T=5.4823\) | supplies the 300 K BM3 term in Equation (29) |

The reconstructed Table III slopes are -5.14, -7.07, -11.49, and
-13.55 MPa K\(^{-1}\), which round to the four reported values. The Table IV
thermal-pressure regression gives 7.144 MPa K\(^{-1}\) and a -0.403 GPa
intercept, reproducing Equation (20)'s 7.14 MPa K\(^{-1}\) and -0.40 GPa.
Separate fits to the rounded low- and high-temperature Table III ambient
\(K_T\) branches give -50.46 and -50.29 MPa K\(^{-1}\), consistent at table
precision with Table II's -51.0 and -50.0 MPa K\(^{-1}\) coefficients.

## Equation represented by Peritheos

With \(V_a\) the ambient four-atom fcc-cell volume and \(T\) in kelvin,

\[
P(V,T)=P_{BM3}(V;V_a,166.65,5.4823)
+\left[0.00714-0.0115\ln\left(\frac{V_a}{V}\right)\right](T-300).
\]

The sign follows the stored `LogVolumeThermalPressure` convention:
`dK_dT_V=-0.0115 GPa/K` multiplies `ln(V0/V)`. Evaluating this model at all
126 top-row Table V states gives an RMS difference of 0.00270 GPa and a maximum
absolute difference of 0.00496 GPa, exactly the expected 0.01 GPa table-rounding
limit.

The ambient volume is independently recoverable from the paper's density and
atomic mass:

\[
V_a=\frac{196.967}{19.30}\frac{4\times10^{24}}{N_A}
=67.7868\ \text{A}^3,
\]

which rounds to the stored 67.79 A\(^3\).

## Table I 300 K heat-capacity inconsistency

Table I visibly prints \(C_p=1.1288\) in units of 0.1 J g\(^{-1}\) K\(^{-1}\)
at 300 K. Used literally in Equation (12), that gives \(\gamma=3.395\), not the
Table IV value 2.974. Two independent reconstructions agree instead:

- Equation (4) gives 1.28820 in the Table I units.
- Inverting Equation (12) using the printed \(\alpha\), \(K_S\), density, and
  Table IV \(\gamma\) gives 1.28843.

The dataset preserves `cp_printed_0p1_j_g_k=1.1288` and separately records the
evident intended value `cp_reconstructed_0p1_j_g_k=1.2888`. No source value is
silently overwritten.

## Why a global refit is unavailable

A source-faithful global optimization cannot be reconstructed because the
article does not publish:

- the points underlying the hand-smoothed low-temperature expansivity curve;
- a common observation matrix or cross-property residual definition;
- weights relating calorimetric, expansivity, and ultrasonic constraints;
- the exact interpolation and numerical-integration algorithm for Equations
  (8), (13), and (19); or
- a joint parameter covariance matrix.

Moreover, \(K_{T0}\) and the final \(K'_T=5.4823\) are explicitly adopted from
Heinz and Jeanloz rather than optimized by Anderson et al. A newly designed
P-V regression would answer a different scientific question and is therefore
not introduced.

## Reproducible artifacts

The audit adds machine-readable transcriptions of Tables I-V in
`peritheos/data/datasets/`, an executable reconstruction in
`scripts/reproduce_anderson_1989_gold_thermal_eos.py`, and a committed result at
`docs/data/anderson-1989-gold-thermal-eos-reproduction.json`.

Run:

```bash
uv run python scripts/reproduce_anderson_1989_gold_thermal_eos.py --check
```

This checks the staged regressions, the separate Table II coefficient branches,
the 300 K heat-capacity inconsistency, the ambient-volume conversion, all
dataset hashes, and all 126 Equation (29) Table V values without pretending the
derived grid is a fit dataset.
