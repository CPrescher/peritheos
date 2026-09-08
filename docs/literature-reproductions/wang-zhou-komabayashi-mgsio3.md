# Wang, Zhou and Komabayashi: original-scale MgSiO3 EOS refits

Four independently published EOS records are now available: two Wang (2004)
akimotoite BM3 fits, Zhou's (2014) isothermal akimotoite BM3, and Komabayashi's
(2008) thermal post-perovskite BM3. All recovered free coefficients agree within
the papers' printed errors. This does **not** reconstruct Dorogokupets' later
joint thermoelastic optimization, and it does not establish the original
authors' exact weighting or optimizer. The common ledger therefore labels these
results `similar`, not strict statistical `parity`.

Published coefficients remain in the catalog; refitted values are diagnostics.
Existing default EOS selections are unchanged. These are four fits, not four
new independent experimental datasets: they reuse the three bundled source tables.

## Published versus recovered

Parentheses below retain the source's last-digit error notation, without
inferring a confidence convention. Volumes are in A3 per conventional cell and
bulk moduli in GPa.

| Source and observations | Fixed values | Published free coefficients | Independent refit |
| --- | --- | --- | --- |
| Wang T0133, 13 rows at 298 K, Anderson-Au pressures | K0=210 | V0=264.2(2); K0'=4.8(5) | V0=264.3240; K0'=4.4577 |
| Wang T0150, 10 rows at 298 K, Decker-NaCl pressures | K0=210 | V0=263.9(2); K0'=5.6(8) | V0=263.9299; K0'=5.5689 |
| Zhou, 17 rows at 300 K, Tsuchiya-Au pressures | K0'=4.6 | V0=262.45(26); K0=207(3) | V0=262.4092; K0=206.7705 |
| Komabayashi PPv, 3 room-temperature then 21 P-V-T rows, Speziale-MgO pressures | V0=163.813; K0'=4; alpha(300)=1.7e-5 K^-1 | K0=223.2(2); alpha1=0.113(114)e-8 K^-2; dK/dT=-0.0085(11) GPa/K | K0=223.1706; alpha1=0.10970e-8; dK/dT=-0.0084156 |

The Wang Au derivative is the least exact recovery, but its difference of
0.3423 is smaller than the printed error of 0.5. The preferred NaCl fit and
Zhou's static fit are particularly close. Komabayashi's alpha1 has a large
published uncertainty: agreement alone does not establish precise thermal
expansion independently of the imposed constraint.

## Data and reconstruction choices

- **Wang:** retain all 134 P-V-T rows in Tables 1-2, but fit each run's complete
  298 K subset separately. Keep the final -0.05 GPa T0150 point. The authors
  prefer the NaCl result because of concerns about stress in the Au marker.
  Do not mix the hot data into the room-temperature BM3 fit or rescale to
  Sokolova. [Primary paper](https://doi.org/10.1016/j.pepi.2003.08.007),
  section 3.2 and Tables 1-2.
- **Zhou:** retain all 58 akimotoite rows, including 55 acoustic measurements.
  The static EOS refit uses the 17 printed 300 K densities, including three
  ambient density-only observations. Convert with `V=6*M/(N_A*rho)`, using
  M=100.387 g/mol and the exact Avogadro constant. These are volumes recovered
  from rounded densities, not the unavailable unrounded diffraction output.
  K0T=207 GPa is isothermal; the separately published acoustic K0S=219.4 GPa
  is not an alternative static EOS coefficient.
  [Primary paper](https://doi.org/10.1016/j.pepi.2013.06.005), section 3.1 and Table 1.
- **Komabayashi:** retain all 22 table rows; exclude only the row lacking a PPv
  volume. Fit K0 to the three 300 K points first, then fix it in the thermal
  stage using all 21 PPv observations. V0 and K0' remain fixed throughout.
  [Primary paper](https://doi.org/10.1016/j.epsl.2007.10.036), Tables 1-2,
  equations 1-4 and the PPv thermal-expansion constraint on page 522.

The primary residual is `P_model(V,T)-P_printed`, with equal weight per selected
observation. No source error bars are used as weights and no optimizer-derived
covariance is presented as a published uncertainty. Source row masks, hashes,
pressure ranges and RMSE are included in the machine-readable report.

## Komabayashi thermal equation and staging

The implemented source equation is BM3 evaluated with:

```text
K0(T) = K300 + (dK/dT)*(T-300)
alpha(T) = alpha0 + alpha1*T
V0(T) = V300*exp(alpha0*(T-300) + alpha1*(T*T-300*300)/2)
alpha0 = 1.7e-5 - 300*alpha1
```

Thus alpha0 is dependent, not an additional unconstrained fit coefficient.
The refit gives alpha0=1.6670901e-5 K^-1. The catalog stores the actual printed
alpha0=1.667e-5, including its very small rounding inconsistency with the
printed alpha1 and the exact ambient constraint.

Using the published rounded K300=223.2 instead of the independently fitted
223.1706 in stage two gives alpha1=0.11773e-8 and dK/dT=-0.0085105. Both choices
remain consistent with the publication. The staged refit's all-temperature
pressure RMSE (0.525963 GPa) is slightly higher than the printed curve's
(0.525096 GPa): its K300 minimizes the room-temperature objective, not the joint
all-temperature objective. This is expected for the stated two-stage procedure.

## Scope and reproducibility

The Wang thermal fits also use external zero-pressure expansion measurements;
Komabayashi's companion perovskite fit uses Funamori et al. data. Those are not
claimed reproduced here. The existing Zhou acoustic diagnostics are retained
as observations and regressions, not misclassified as equilibrium P-V-T records.

```sh
uv run --frozen python scripts/reproduce_wang_zhou_komabayashi_mgsio3.py
uv run --frozen python scripts/reproduce_wang_zhou_komabayashi_mgsio3.py --check
uv run --frozen pytest -q tests/test_wang_zhou_komabayashi_mgsio3.py
```

Results: [machine-readable refits](../data/wang-zhou-komabayashi-mgsio3-refit.json).
The record-by-record [common ledger](../primary-eos-refits.md) includes all four
records and their fixed/free parameter masks. See also the separate
[Dorogokupets reconstruction](dorogokupets-2015-mgsio3-refit.md).
