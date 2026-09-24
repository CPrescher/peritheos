# Xiao et al. (2025): solid argon Helmholtz equation

**Source:** Xiao, Sriskandaruban, Maynard-Casely, Thol, Falloon, Span and May,
*International Journal of Thermophysics* **46**, 14 (2025),
[doi:10.1007/s10765-024-03469-2](https://doi.org/10.1007/s10765-024-03469-2).
The full [author preprint](https://doi.org/10.21203/rs.3.rs-5370077/v1)
was verified against the publication's official calculation workbook. The
publisher's typeset article was not downloaded. See the source manifest for
version, origin, local paths and SHA-256 hashes.

## Exact model and reference state

The fcc solid has four atoms per conventional cell. The source uses molar
volume in cm³/mol, MPa, K and J/mol. The bundled material uses conventional
cell Å³ and GPa; the molar Python thermal API uses J/bar/mol. Consequently
V00 = 22.555 cm³/mol = 149.81383464042446 Å³/cell = 2.2555 J/bar/mol.
V00 is a fixed extrapolated 0 K, zero-pressure volume (section 5), not a
measured ambient room-temperature solid volume.

With z=V00/Vm and L=ln(z), the cold Helmholtz energy (Eq. 18) is
V00[c1 L²/2+c2 L³/3+c3 L⁴/4]. It includes zero-point energy; adding a
separate zero-point Debye term would double count it. The exact NaturalStrain4
mapping is K0=c1, K0'=2+2c2/c1, and
K0''=[6c3/c1−1−(K0'−2)−(K0'−2)²]/K0. Table 4 gives c1=2656.5 MPa,
c2=7298 MPa, c3=10 MPa. The transformed derivatives are algebraic parameter
conversions, not a new refit.

The Debye term is RT[3ln(1−exp(−theta/T))−D3(theta/T)].
Gamma=2.68(Vm/V00)^0.0024 and
 theta=86.44 exp{(2.68/0.0024)[1−(Vm/V00)^0.0024]} K (Eqs. 19–24).
The anharmonic term (Eq. 25) is
b1 R theta0 (T/theta0)^4/[1+b2(T/theta0)^2] exp[b3(Vm/V00−1)],
with b1=0.0128, b2=0.388 and b3=7.85.
The author workbook uses R=8.31451 J/mol/K; this value is explicitly retained
for exact published-model reproduction. Using the current SI R shifts the
workbook pressure by about 0.00072 MPa at Vm=23 cm³/mol, 70 K.

The thermal contribution is absolute relative to 0 K. Tr=70 K only anchors
reported thermal increments and matches the workbook checkpoint; it does not
subtract a 70 K isotherm from total pressure. Gas reference entropy 130.37
J/mol/K is used in the coupled solid-fluid phase-equilibrium calculation and
does not enter solid P(V,T). The implementation does not claim fluid fugacity
or phase-equilibrium calculations.

## Source discrepancies and coverage

The official supplementary DOCX labels Table S1 “Optimized parameter values”
but gives c1=2897, c2=8000, c3=941 and most Debye/anharmonic parameters equal
to one. Section 5 identifies these as the initial guesses. They are not the
fitted coefficients. Main Table 4 agrees with the executable VBA and cached
numerical results in the official XLSM; these are the coefficient authority.

The advertised marginal envelope is up to 760 K and 6300 MPa, with a
zero-temperature limit. It is not a rectangular guarantee that argon is a
stable solid. Calculations must follow the solid branch below melting. The
multiproperty fit uses heterogeneous legacy pressure scales rather than a
single recoverable common calibration. Neither coefficient uncertainties nor
covariance are supplied. Source estimated k=1 volume uncertainties are 0.1%
along sublimation and 0.5% along melting and for compressed solid.

The 22 new neutron observations in Table 5 cover 7.96–50.02 K. All rows are
bundled with u(T)=0.02 K and relative u(V)=0.001 at k=1. Their printed quantity
is molar volume, despite the table's “cell volume” description. The sample
pressure was slightly above sublimation and was not tabulated (section 3).
No measured pressures have been fabricated. The source calls their volume
effect negligible; diagnostics use a clearly stated P=0 approximation.

## Independent diagnostics

Run `python -m scripts.reproduce_argon_xiao_2025` to reproduce
[the report](../data/argon-xiao-2025-reproduction.json). It uses direct
Helmholtz derivatives and independent adaptive SciPy Debye quadrature.
At 70 K and 23 cm³/mol, the official workbook gives P=78.3517375843629 MPa,
Cv=22.854917560882143 J/mol/K, A=−998.2296606304144 J/mol and
K=2247.940903627391 MPa. A, U, S, P, Cv and K agree to floating-point
precision. Three independent workbook volume inversions agree within
6e−15 relative. These checks are calculated theory values, not measurements.

The Table 5 P≈0 comparison gives 0.08073% AAD and 0.09518% RMS relative
volume discrepancy, consistent with the reported 0.1% experimental uncertainty.
No uncertainty was inferred from these residuals. The original global
multiproperty refit remains **not reproduced**: the consolidated selected
legacy rows, penalty functions and iterative fitting protocol are unavailable.
An artificial fit to only the 22 new rows would not reproduce that regression.

## Equal-weight comparison requested during review

The original weights are **published**, not unknown. Section 5, Eq. 62 and
Table 2 define a relative-residual objective with empirical property weights
(from 1 for compressed-solid volume to 210 for heat capacity above 12 K),
reduced weights for selected data, exclusions, and physical penalty terms.
The numerical weights are transcribed in the
[fitting-protocol report](../data/argon-xiao-2025-fitting-protocol.json).

The requested alternative uses equal weight per recovered observation in
squared relative residuals. Relative normalization follows the source and
avoids combining numerical values with incompatible units. All recovered
rows are retained and their provenance and coverage are reported explicitly.
This is an alternate regression assumption; agreement with it would not
establish that the authors used equal weights.

The official workbook contains sample calculations and fluid-EOS extrapolations,
not the consolidated experimental inputs summarized in Tables 6–8. Those
calculated values are not substituted for missing observations. Solid-only
data also cannot identify the adjustable gas-reference entropy parameter.

The [equal-weight diagnostic](../data/argon-xiao-2025-equal-weight-refit.json)
now records five starting points for each full nine-parameter and restricted
(two thermal parameters) fit. For all 22 neutron volumes, the published
volume RMS is 0.09518%; the best converged full fit gives 0.03224%, while
fitting only theta0 and gamma0 gives 0.03399%. The full fit's nearly singular
sensitivity matrix and parameter-bound hits prevent unique coefficient recovery.

A separate 27-entry case includes all five Anderson–Swenson Table 1 K0
estimates, clearly labeled extrapolated isotherm-fit summaries rather than
independent zero-pressure measurements. Its best converged full fit has about
0.06482% combined relative RMS, but remains poorly constrained and bound-dependent.
Xiao originally selected two of these five estimates; the alternate comparison
retains all five as requested. Every numerical candidate, convergence flag,
chosen bound and source hash is retained. Published catalog coefficients remain
unchanged; the missing full experimental collection is not replaced with theory.

The lowest-residual full nine-parameter 27-entry candidate also gives negative
isochoric heat capacity at some sampled states. Its small residual is therefore
not evidence of a physically acceptable replacement EOS. Every solution now
includes sampled heat-capacity and bulk-modulus admissibility flags.
