# Speziale et al. (2001): MgO thermal equations and the later Au calibration

**The existing integrated MgO calculation implements the printed variable-q
gamma equation correctly. It does not reproduce the pressures labeled Speziale
MgO in Fei (2004) and Hirose (2006).** A different gamma convention explains
almost all of that gap numerically. This is evidence for a pressure-reduction
convention difference, not recovery of the authors' executable calibration.

The cold record `mgo_speziale_2001_bm3_2` remains a cold-only constrained BM3
record. No thermal model is silently attached to it. The published Au
coefficients and qualified staged parity remain unchanged. This audit does not
refit Au, Pt, or the independent Miozzi iron reconstruction.

## Primary authority and reference state

The complete relevant primary source is
[Speziale et al., JGR 106, 515–528 (2001)](https://duffy.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf),
DOI [10.1029/2000JB900318](https://doi.org/10.1029/2000JB900318).
The newly retrieved author PDF has SHA-256
`c0431bdddf8f4bc8269a7e0d1b79d30d92a38b4a778d9edddad9cb4c45b1f94b`,
matching the earlier cached copy. Equations 1–4, Sections 4.2 and 5.1–5.4,
Tables 1–3, and Figures 6, 9 and 10 were checked, including the rendered
equations and tables. Fei (2004) Tables 1 and 3 were checked against the cached
original article. The Hirose inputs retain the verified publisher transcription.

With x=V/V0, the source equations are

\[
q(V)=q_0x^{q_1},\qquad
\gamma(V)=\gamma_0\exp\left[\frac{q_0}{q_1}(x^{q_1}-1)\right],
\]

\[
\ln\frac{\theta(V)}{\theta_0}=-\int_1^x\gamma(u)\,d\ln u,
\qquad
P(V,T)=P_{BM3}(V,300)+\frac{\gamma(V)}{V_m}
[E_D(\theta(V),T)-E_D(\theta(V),300)].
\]

Here V is the conventional cell in angstrom^3,
V_m=V*N_A*10^-30/4 in m^3/mol MgO, and the Debye energy uses n=2 atoms
per MgO formula unit. The equivalent cell basis is n=8 and Z=1. All pressures
are converted from Pa to GPa after dividing energy by molar volume.

V0=74.71, K0=160.2 GPa, K0-prime=3.99, gamma0=1.524, q0=1.65,
q1=11.8, theta0=773 K and Tr=300 K are the relevant final model coefficients.
The abstract's gamma0=1.49(3) is the preliminary static thermal fit, not the
fixed final gamma0. Table 3 is a set of **constant-q** fits: it is not a fit of
q0 and q1 simultaneously. q0 is constrained by thermodynamic derivatives,
while Section 5.4 constrains q1 from q near 0.01 at 200 GPa and the shock
comparison. This is not one recoverable flat unconstrained PVT regression.

Equation 4, theta=theta0*exp[(gamma0-gamma)/q], follows for constant q.
Substituting instantaneous q(V) into that expression is not its integrated
variable-q extension. The paper does not give an explicit executable variable-q
theta prescription. Integration is the thermodynamically consistent extension,
not proof of the original software's choices.

## What is independently validated

The new implementation in
[`reproduce_speziale_2001_mgo.py`](https://github.com/CPrescher/peritheos/blob/main/scripts/reproduce_speziale_2001_mgo.py)
checks both differential identities, the constant-q and q=0 limits, and
equivalent formula-unit/cell energy normalizations. Independently coded adaptive
integrals agree with the 64-point Gauss-Legendre evaluator on 35 states from
x=0.60 to 1.02 and 300–3663 K to **3.70e-13 GPa**. The pre-existing Au
MgO diagnostic agrees to **1.60e-13 GPa** on all 28 calibration rows.

At x=0.64, gamma=1.326084; the infinite-compression limit is 1.325127.
This agrees with the approximately 1.33 high-compression limit and monotonic
decrease described by Speziale and shown in Figure 9. It independently rules
out silently interpreting Eq. 11 as a local-power gamma law.

| Recoverable original constraint | Result | Qualification |
| --- | --- | --- |
| Table 1, 32 cold rows, K0 fixed | Equal-pressure-residual replay: V0=74.712174, K0-prime=3.995497, within printed widths | The existing library fit gives V0=74.711431 and K0-prime=3.996280; both are cold parity, not a thermal validation |
| Dewaele (2000) Table 2, 41 hot rows | Constant-q fit: q=1.237317 versus Speziale's 1.3(5) | Fixed Speziale anchors; equal pressure weights; original thermal objective/weights unknown |
| Same 41 hot rows, final Eq. 11 | Mean absolute relative volume discrepancy 0.382299%; pressure RMS 1.108894 GPa | Close to the paper's 0.4% aggregate volume statement, but this is only its Dewaele subset |
| Fiquet (1999) COD, 34 rows above 300 K | Constant-q fit: q=0.607759 versus 1.0(4) | Near the lower printed width; 36-of-37 total COD recovery, rounded volumes and unknown original selection/weights |
| Svendsen–Ahrens (1987) Table 6, four shots | Equation-consistent model leaves roughly 6–9 GPa pressure discrepancies for absolute-density-derived volumes | This conditional reconstruction does not establish the original shock tolerance |

The new shock table is recovered from the
[original Caltech-hosted article](https://web.gps.caltech.edu/~asimow/TJA_LindhurstLabWebsite/ListPublications/Papers_pdf/Seismo_1409.pdf),
Table 6 on page 685. Its greybody temperatures are 3081, 3071, 3281 and
3663 K. **The adjacent 2913, 3028, 3257 and 3667 K values are calculated model
temperatures, and are not used as measured thermal targets.** Pressures and shock
velocities also belong to the model-estimate block. Deriving compression from
them uses correlated reduction outputs, not independent measured P/V inputs.

The archive preserves two explicit shock-volume normalizations. The absolute
cell volume derives from each initial density, M=40.304 g/mol and mass/momentum
conservation. Its residuals are -8.79, -5.97, -8.72 and -8.09 GPa; predicted
volume discrepancies are -1.29%, -0.85%, -1.21% and -1.07%. Multiplying each
shot's relative compression by adopted V0 instead gives -4.82, -8.22, -7.49
and -8.86 GPa. Neither route recovers all four central points within the stated
6 GPa tolerance. Source/model errors, their correlations and the original
shock normalization prevent assigning this discrepancy to a particular author
calculation. The full published thermal analysis is therefore **not reproduced**.

## Diagnosis of Fei/Hirose pressure discrepancies

The audit keeps all reported pressure targets intact. The alternative convention
is explicitly a diagnostic substitution:

\[
\widetilde\gamma(V)=\gamma_0x^{q(V)}.
\]

It uses the printed q0 and q1 without optimizing them. Its Debye temperature
is integrated from this *alternative* gamma. It differs from Eq. 11 because

\[
\frac{d\ln\widetilde\gamma}{d\ln V}
=q(V)[1+q_1\ln x],
\]

not q(V). The derivative turns negative for x<exp(-1/11.8)=0.91875,
and gamma eventually returns to gamma0 with compression, in disagreement
with Speziale's Figure 9. Numerical agreement with reported pressures does
not make this the source-consistent physical model.

| Reduction convention | Fei 26-row RMS / max (GPa) | Hirose pressures at 1340 / 2330 K (GPa) |
| --- | --- | --- |
| Printed Eq. 11, integrated theta, 300 K | 0.254647 / 0.361844 | 105.609013 / 115.947872 |
| Local-power gamma, integrated theta, 300 K | 0.008997 / 0.014237 | 106.493943 / 117.880818 |
| Local-power gamma, integrated theta, 298 K | 0.003264 / 0.004900 | 106.502789 / 117.889630 |
| Reported targets | — | 106.5 / 117.9 |

The 300 K substitution removes about 96.5% of the Fei RMS discrepancy and
brings both Hirose rows within 0.02 GPa. Changing Tr to 298 K puts every Fei
row within half a printed pressure digit, but **the 298 K choice is inferred,
not established by Speziale's 300 K equations**. Several alternative theta
prescriptions also agree at the printed input/pressure precision, so the tables
do not uniquely identify the complete reducer.

Other plausible causes were tested separately: constant q, theta held fixed,
Eq. 4 using local q or q0, power-law theta, 298/298.15 K reference subtraction,
R=8.31451, Table 1 measured ambient lattice parameter, and formula-unit/cell
normalizations. With temperature fixed at its printed value, lattice/pressure
last-digit rounding cannot explain the Eq. 11 discrepancy for any of the 28
rows. The half-digit bounds are about 0.011–0.013 GPa for Fei and 0.071 GPa
for Hirose. These deterministic bounds are not measurement uncertainties or
parameter-confidence intervals.

Fei Table 3 prints MgO “3R (J/g K)” = 0.12664, which cannot be used literally
as the MgO formula-unit Dulong-Petit heat capacity: 6R/M=1.237762 J/g K.
The audit treats the literal entry and the 51.041 J/mol K normalization
discussed by [Dorogokupets (2010)](https://doi.org/10.1007/s00269-010-0367-2)
as separate controls. Neither accounts for the full pressure gap. That primary
reanalysis independently documents several differing implementations called
“Speziale MgO”; it corroborates a historical convention issue, not the exact
local-power reducer inferred here.

## Reproduction and remaining evidence

```bash
python -m scripts.reproduce_speziale_2001_mgo
python -m scripts.reproduce_speziale_2001_mgo --check
python -m pytest tests/test_speziale_2001_mgo.py tests/test_fei_2007_gold.py
```

[`speziale-2001-mgo-reproduction.json`](../data/speziale-2001-mgo-reproduction.json)
archives all constraints, sensitivities, hashes and explicit reproduction flags.
[`speziale-2001-mgo-calibration-residuals.csv`](../data/speziale-2001-mgo-calibration-residuals.csv)
retains reported pressures alongside every diagnostic prediction.
[`mgo-svendsen-1987-source.json`](../../peritheos/data/datasets/mgo-svendsen-1987-source.json)
records the original shock recovery and separates optical measurements from
model outputs. No third-party source-rights transfer or source CC0 claim is made.

Still needed for exact reproduction: original selected Fei (1999)/Utsumi
(1998) thermal inputs; original regression objective, weights and covariance;
the variable-q theta prescription and unrounded coefficient/input files used
in Speziale's calculations; exact shock-volume normalization; and the actual
Fei/Hirose pressure-reduction code. The available evidence establishes the
equation implementation and conditional subset agreement, while leaving the
complete author thermal calibration and reducer unidentified.
