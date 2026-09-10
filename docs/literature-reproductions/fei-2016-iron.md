# Fei et al. (2016): hcp-iron

Fei, Y., Murphy, C., Shibazaki, Y., Shahar, A., and Huang, H. (2016),
*Thermal equation of state of hcp-iron: Constraint on the density deficit of
Earth's solid inner core*, Geophysical Research Letters **43**, 6837–6843,
[doi:10.1002/2016GL069456](https://doi.org/10.1002/2016GL069456).

The published article was read from Zotero item `K7V7ZTEW`
(BibTeX key `fei2016ThermalEquationState`). The official supplementary
workbooks and supporting-information document were supplied by the user.
Their SHA-256 hashes and the article hash are recorded in
`peritheos/data/datasets/iron-fei-2016-source.json`.

## Published records

Three published parameterizations extend the existing `iron` hcp material.
Existing catalog defaults and the separate bcc `fe` material are preserved.

| Record | Equation | Source coefficients |
|---|---|---|
| `iron_fei_2016_bm3` | Preferred 300 K BM3 | Fixed density 8.2695 g/cm³, K0 = 172.7(1.4) GPa, K0′ = 4.79(0.05) |
| `iron_fei_2016_free_density_bm3` | Alternative 300 K BM3 | Density 8.3602(0.0241) g/cm³, K0 = 191.44(5.3) GPa, K0′ = 4.52(0.08) |
| `iron_fei_2016_bm3_debye_quadratic` | Preferred BM3 plus Equation (2) | θ0 = 422 K, γ0 = 1.74, q = 0.78, β0 = 0.07 J kg⁻¹ K⁻², k = 1.34, γe = 2 |

The two static fits are printed in Section 3.1, page 6840. Equation (2) and
the thermal coefficients appear in Section 3.2, pages 6840–6841. Density is
converted to conventional hcp cell volume using two Fe atoms per cell,
M = 55.845 g/mol, and the exact Avogadro constant. The preferred V0 is
22.4276689537 Å³. The free-density volume error is propagated linearly from
the reported density error. No confidence level or covariance is inferred.
The existing P63/mmc structure, Z = 2, and occupied 2c site remain independent
of these extrapolated zero-pressure EOS volumes.

## Exact thermal mapping

For x = V/V0 = ρ0/ρ, the vibrational part uses
γ = γ0 x^q and the integrated Debye-temperature law
θ = θ0 exp[γ0(1 − x^q)/q]. The Debye energy difference is taken from 300 K.
The additional term follows directly by integrating the source's
mass-specific Cv,e = β0 x^k T:

\[
P_e = \frac{\gamma_e\beta_0\rho_{0,SI}}{2\times10^9}
      x^{k-1}(T^2-300^2).
\]

The reusable `DebyeQuadraticThermalPressure` model therefore stores
A = 5.78865 × 10⁻⁷ GPa/K² and m = k − 1 = 0.34. Its internal molar volumes
use J/bar/mol, with n = 1 Fe atom per formula; `.eosmat` volumes remain
conventional-cell Å³. This normalization reproduces the mass-specific
expression without changing the source's independently adopted γe = 2.

This is an empirical pressure surface. Since γe differs from k, the electronic
pressure and stated electronic heat capacity do not arise from a single
Helmholtz potential with that volume law. The new class deliberately offers
pressure, volume, temperature, derivatives, fitting, and interchange, but
does not claim caloric properties or an adiabatic geotherm reproduction.

## Complete new observations and scope

Both workbooks contain a single `Sheet1` worksheet. All 74 S1 observations
(rows 3–76) and 22 S2 observations (rows 3–24) are bundled as CSV. Numeric
cached cell values, row order, sample names, densities, lattice parameters,
all reported errors, and pressure-marker identities are preserved. No source
workbook is modified or recalculated.

S1 preserves Ne, NaCl-B2, and MgO marker lattice readings. The 18 absent MgO
lattice errors stay blank. S1 has no row-level pressure or volume errors.
S2 reports temperature standard deviations for averages of upstream and
downstream measurements. The confidence convention of lattice errors is not
specified. Rounding of a and c explains cell-volume differences below
0.0002 Å³; printed volume values are retained.

The methods identify Fei (2007) Ne/NaCl/Pt scales and Tange (2009) MgO,
using only the MgO (200) peak. Figure 1's caption groups MgO with Fei (2007);
the explicit methods take precedence. Exact scale variants are not specified
unambiguously, so no invented executable calibration link is supplied.
Thermal pressures use separate Fe and Pt heating cycles matched through the
Ne peak, as described in the SI introduction. S2 does not contain the matched
Pt lattice observations needed to repeat this reduction.

The global static fit combines earlier studies after pressure intercalibration
and excludes Mao (1990). Those complete re-reduced rows and weights are not
provided by these supplements. The thermal fit additionally uses Brown (2000)
shock data; its row-level temperature reduction and weights are absent here.
The general refit ledger classifies these combined fits as `not_refittable`,
while the independent partial-data diagnostics below remain executable.

Catalog validity stores the measured subset envelopes: S1 spans
18.72–204.6 GPa at 300 K; S2 spans 53.96–99.69 GPa and 1208–1795 K.
These are marginal observation ranges, not rectangular phase-stability
guarantees. The paper's combined compression curve reaches 280 GPa and its
thermal synthesis includes shock states near 200 GPa and 5000 K. Inner-core
applications are extrapolations beyond the new static measurements.

## Independent numerical checks

`scripts/reproduce_fei_2016_iron.py` implements Equation (1), numerical Debye
quadrature in SI units, and the mass-specific electronic term independently
of Peritheos. Its deterministic output is
`docs/data/fei-2016-iron-reproduction.json`.

| Diagnostic | Result |
|---|---:|
| Preferred BM3 RMS pressure residual, all 74 S1 rows | 1.37690 GPa |
| Free-density BM3 RMS pressure residual, all S1 rows | 1.35266 GPa |
| Complete thermal model RMS residual against all 22 S2 P values | 1.60934 GPa |
| Thermal-pressure RMS residual against all S2 Pth values | 1.13755 GPa |
| S2 first state: calculated / observed total pressure | 53.87281 / 53.96 GPa |
| S2 first state: calculated / observed thermal pressure | 6.56267 / 6.7 GPa |
| S2 last state: calculated / observed thermal pressure | 11.47860 / 12.2 GPa |

The S2 Pth column differs from P minus the preferred printed 300 K BM3 by
up to 2.06562 GPa. Both published quantities are retained independently;
the difference is not silently corrected or interpreted as an uncertainty.
Numerical tests compare independently evaluated equations tightly, while
measurement comparisons use explicit residual tolerances.

Unweighted S1-only pressure refits give K0 = 174.1451 GPa and K0′ = 4.74696
with fixed density, or V0 = 22.05149 Å³, K0 = 205.8664 GPa and K0′ = 4.31570
with free density. A one-parameter S2-only diagnostic, fixing k and all
vibrational parameters, gives β0 = 0.1169811 J kg⁻¹ K⁻². This is not the
paper's two-parameter combined shock/static fit. These diagnostics neither
replace published coefficients nor create catalog refit records.

Reproduce the bundled diagnostics:

```bash
uv run python scripts/reproduce_fei_2016_iron.py
```

To verify every CSV field against the supplied workbooks, pass
`--source-directory` pointing to their containing directory. `--write`
explicitly regenerates the normalized CSVs and diagnostic report.

## Newly retrieved upstream sources

The [five-paper primary-source audit](iron-fei-input-papers.md) now supplies
original Yamazaki, Sakai and Brown tables, alongside the existing Dewaele
EPAPS data. Their original pressure scales and absent shock temperatures
remain explicit. These additions do not recover the missing Dubrovinsky
rows or Fei's exact combined regression and Table S2 pressure reduction.
