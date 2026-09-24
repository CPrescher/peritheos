# Finger et al. (1981): crystalline fcc argon

The executable record `argon_fcc_finger_1981_murnaghan2_debye` represents
[Finger, Hazen, Zou, Mao and Bell, Applied Physics Letters 39(11), 892–894
(1981)](https://doi.org/10.1063/1.92597). The author-hosted
[three-page primary PDF](https://hazen.carnegiescience.edu/sites/default/files/061-argon-1981.pdf)
was downloaded and all pages visually inspected. Its title, authors, journal,
year and page numbering match. The older `/sites/hazen.gl.ciw.edu/files/` URL
in the structural citation returned HTTP 404; that citation now uses the working
URL. PDF provenance and Zotero keys are in
[`argon-fcc-finger-1981-sources.json`](../../peritheos/data/datasets/argon-fcc-finger-1981-sources.json).
The copyrighted PDF is outside git.

## Published equation and reference state

The paper's “second-order Murnaghan” is **not Birch–Murnaghan BM2**. With
`x=(V0/V)^q`, Eq. (6) is

```
q = sqrt(K0_prime^2 - 2*K0*K0_double_prime)
Ps = Pz0 + 2*K0*(x-1)/(q*(x+1) - K0_prime*(x-1))
```

This integrates `K(Ps-Pz0)=K0+K0_prime*(Ps-Pz0)+K0_double_prime*(Ps-Pz0)^2/2`.
The new `SecondOrderMurnaghan` class implements this real positive-q branch in
Python/Rust with explicit static offset `P0=Pz0`. It is registered in EOSMAT,
material construction, native evaluation and fitting.

Equations (2)–(5) add the zero-point pressure `Pz` and thermal pressure `PT`:

```
gamma = 1/2 + gamma1*(V/V0)
theta = theta0*(V0/V)^(1/2)*exp(gamma1*(1-V/V0))
Pz = 9*gamma*R*theta/(8*Vm)
PT = 3*gamma*R*T*D3(theta/T)/Vm
P = Ps + Pz + PT
```

Here `Vm` is molar volume and pressure conversion must accompany its units.
The printed Eq. (1) has an apparent extra equals sign before `PT`; the additive
decomposition follows the text's definitions and Eqs. (2)–(6). It also matches
the plotted pressure scale and measured points. This typography issue is not
silently treated as a different physical model.

The existing `Dewaele2006` implementation has **exactly** this Debye/gamma law
when `gamma0=2.7`, `gamma_inf=0.5`, `beta=1`, and anharmonic/electronic coefficients
are zero. The new explicit `zero_point_pressure="included"` option, combined
with `thermal_pressure_reference="absolute_zero"`, supplies Eq. (2). This option
requires the absolute-zero baseline; the default remains `"omitted"`.
`Tr=293 K` anchors pressure **increments** only. Native fit reconstruction now
preserves both baseline and zero-point configuration, including previously
available absolute-zero models.

Table II supplies:

| Quantity | Printed value | Stored value |
|---|---|---|
| V0, fixed | 22.557 cm³/mol | 149.827118952962 Å³ / four-atom fcc cell |
| K0, fixed | 23.701 kbar | 2.3701 GPa |
| K0′, fitted | 6.97 ± 0.11 | same |
| K0″, fitted | −0.040 ± 0.010 kbar⁻¹ | −0.40 ± 0.10 GPa⁻¹ |
| Pz0, fixed | −1.0289 kbar | −0.10289 GPa |
| theta0, adopted | 93.3 K | same |
| gamma1, adopted | 2.20 | gamma0=2.70, gamma_inf=0.50 |

V0 and K0 are adopted low-temperature constraints from Peterson et al. (1966),
not 293 K zero-pressure solid properties. Theta0 comes from Finegold and
Phillips (1969), and gamma1 from Tilford and Swenson (1972), as cited in Table II.
No errors are invented for these fixed parameters. The two fitted errors are
retained; their confidence level and parameter covariance are not specified.

**The printed reference offset is internally imperfect:** Eq. (2), using the
printed parameters and modern R, gives Pz(V0)=0.104460181 GPa. Consequently the
printed total at V0 and 0 K is 0.001570181 GPa, not precisely zero. At V0 and
293 K the total is 0.776333131 GPa. All printed coefficients are preserved;
neither the offset nor the reference volume was adjusted to force cancellation.

## Observations, phase, scale and ranges

All 19 argon rows of Table I are retained in
[`argon-fcc-finger-1981-table1.csv`](../../peritheos/data/datasets/argon-fcc-finger-1981-table1.csv).
They are **measured single-crystal diffraction points**, at 293 ± 1 K and
1.28–8.17 GPa, identified as fcc Fm-3m, Z=4. No solid-solid transition was
observed in that interval. The separate crystal/liquid coexistence observation
is 11.5(5) kbar = 1.15(5) GPa at 293 K; no volume is given, so it is documented
here rather than fabricated as a P–V row. Neon measurements are outside this
argon-only addition.

The table retains source kbar pressures, lattice constants in Å, molar volumes
in cm³/mol, and all parenthetic errors. Derived public cell volumes are `a³`;
volume errors use first-order propagation `3*a²*da`. Modern Avogadro conversion
is `Vm = Vcell*NA*1e-24/4`. Source molar volumes and modern cell-derived molar
volumes differ slightly (maximum 0.00703 cm³/mol), including source rounding
and historical conversion differences; neither printed column is overwritten.
The fit audit uses the printed molar-volume column and separately reports
cell-edge-based residuals.

Pressure was determined from ruby R1 fluorescence using Piermarini et al.
(1975), reference 5, DOI 10.1063/1.321957. The text states pressure measurement
to ±0.5 kbar; individual table errors are preserved separately. Raw wavelengths
are absent, so these pressures cannot be independently recalibrated.
The [calibration paper's author-institution copy](https://physics.byu.edu/download/publication/2242)
confirms the DOI. No modern ruby scale was substituted.

The paper's Fig. 2 isotherms at 4, 77, 293, 400 and 500 K are **calculated
curves**. Only the 293 K squares are observations. The catalog's experimental
range remains 292–294 K, describing the measured temperature and its variation,
not validation of the entire theoretical temperature range or a phase boundary.
Wittlinger's separate hcp material and `not_reproduced` outcome are untouched.

## Reproduction and tests

Run `python -m scripts.reproduce_argon_finger1981` to regenerate
[`argon-finger-1981-reproduction.json`](../data/argon-finger-1981-reproduction.json).
`peritheos.eos.finger1981.Finger1981Argon` independently evaluates Eqs. (1)–(6)
in molar units with adaptive quadrature, providing component-level diagnostics.
The registered model uses the native Debye implementation. Their largest
pressure difference on the five source theory temperatures is 3.10e-11 GPa.
The independent evaluator also supplies the analytic T=0 limit; the standard
thermal material API retains the library's positive-temperature requirement.

The printed parameters give RMS residual 0.095355 GPa using source molar volumes
(0.097220 GPa using lattice-derived cell volumes), and maximum absolute residual
0.180378 GPa. Illustrative observed/calculated pressures are 1.28/1.29898,
3.33/3.28318, 4.77/4.89032 and 8.17/8.03295 GPa.

Diagnostic refits keep every source fixed parameter fixed:

| Objective | K0′ | K0″ (GPa⁻¹) | RMS (GPa) |
|---|---:|---:|---:|
| Equal pressure weight | 6.903805 | −0.336633 | 0.094127 |
| Quoted pressure errors | 6.829505 | −0.253168 | 0.097617 |
| Effective pressure + volume errors | 7.106669 | −0.514831 | 0.103521 |

These are objective-sensitivity diagnostics, **not a recovery of the unpublished
original weights or covariance**, and do not replace the published curve.
Equation/parameter transcription and independent equation evaluation are
validated; exact recovery of the original regression is not claimed.

Tests cover native/reference/fallback equation agreement, zero-point and
baseline semantics, a pressure derivative check against the analytical bulk
modulus, pressure-volume inversion, source checksums/normalization, the
first-order-Murnaghan limit, EOSMAT round trips and preservation of thermal
configuration during native fitting. Studio should evaluate the registered EOS
and draw the 19 observation rows separately.
