# Chantel et al. (2012): Mg and Fe-bearing bridgmanite

Primary DOI: <https://doi.org/10.1029/2012GL053075>
Open final article: <https://zenodo.org/records/3653697>

## Outcome

The audit resolves the apparent contradiction between a bundled density/velocity
table and a `not_refittable` ledger entry. Chantel et al. did perform an
executable acoustic fit, but did **not** perform a simultaneous acoustic and
thermoelastic fit:

- eight 300 K MgSiO3 density--Vp--Vs rows were fitted with the third-order
  Eulerian acoustic finite-strain equations of Davies and Dziewonski (1975), as
  implemented by Li and Zhang (2005);
- the preferred room-temperature coefficients additionally include numerical
  observations from Li and Zhang (2005), which Chantel et al. do not republish;
- the two 700/1200 K rows only test a Stixrude-type thermoelastic model. Paragraph
  15 explicitly says they are insufficient to refine thermal parameters;
- theta0, gamma0, q, and eta in Table 3 are italicized terms adopted from Xu et
  al. (2008), not results fitted to Table 1.

Peritheos can now execute and fit the 300 K acoustic finite-strain model with
`EulerianFiniteStrainAcoustic` and `fit_acoustic_finite_strain`. It can also
execute the stored BM3--Mie-Grüneisen-Debye pressure surface. It cannot represent
the complete temperature-dependent velocity model because the coupled thermal
shear formulation using G, G', and eta is not part of the EOS record.

## Reconstructed acoustic objective

For reference density rho0 and compression-positive Eulerian strain

\[
f=\frac{1}{2}\left[\left(\frac{\rho}{\rho_0}\right)^{2/3}-1\right],
\]

the third-order relations are

\[
K_S(\rho)=(1+2f)^{5/2}
\left[K_{S0}+(3K_{S0}K'_{S0}-5K_{S0})f\right],
\]

\[
G(\rho)=(1+2f)^{5/2}
\left[G_0+(3K_{S0}G'_0-5G_0)f\right],
\]

with

\[
V_S=\sqrt{G/\rho},\qquad
V_P=\sqrt{(K_S+4G/3)/\rho}.
\]

The units g/cm3, km/s, and GPa are mutually consistent in these identities. The
source says only that volume and velocity data were least-squares fitted. It does
not publish whether residuals were formed in velocity, modulus, or transformed
finite-strain space; it also omits numerical weights and covariance. The
reproduction therefore reports two inspectable alternatives rather than claiming
bitwise recovery:

| Reconstruction | K_S0 (GPa) | K'_S0 | G0 (GPa) | G'_0 |
|---|---:|---:|---:|---:|
| Published, this study | 247(4) | 4.5(2) | 176(2) | 1.6(1) |
| Unweighted simultaneous Vp/Vs residuals | 246.402 | 4.54295 | 177.208 | 1.60648 |
| Diagonal density/Vp/Vs errors-in-variables sensitivity | 246.395 | 4.54358 | 177.209 | 1.60650 |

Both reconstructions recover all four current-study coefficients within combined
two-sigma uncertainty. The unweighted velocity RMSE values are 0.01925 km/s for
Vp and 0.01447 km/s for Vs. The diagonal errors-in-variables fit has reduced
chi-square 0.154 and adjusts density by at most 0.000303 g/cm3, but its formal
errors must not be mistaken for the authors' errors because correlations and the
confidence convention are unavailable.

The preferred combined Table 2 row is K_S0=252(1) GPa, K'_S0=4.1(1), G0=175(1)
GPa, and G'_0=1.7(1). Evaluated only against the eight Chantel rows, it gives Vp
and Vs RMSE values of 0.02250 and 0.02104 km/s. That is a useful curve check, not
a refit of the preferred coefficients, because the Li--Zhang input rows are
missing.

## Table 3 semantic boundary

Table 3 labels the bulk entry K_T=252 GPa, while Table 2 labels the identical
combined acoustic result K_S=252(1) GPa. The article does not document an
adiabatic-to-isothermal conversion. Peritheos retains the printed Table 3 pressure
parameter but records this unresolved label transition explicitly. The 300 K BM3
curve evaluated at all nine room-temperature density states has pressure RMSE
0.628 GPa and maximum absolute residual 1.166 GPa. Those Au-scale pressures are
independent curve checks, not inputs to the pressure-independent acoustic fit.

The dataset flags each role directly: one ambient rho0 anchor, eight acoustic-fit
rows, and two high-temperature validation-only rows. No row-level pressure errors,
within-row density/Vp/Vs correlations, four-parameter acoustic covariance,
combined-fit numerical Li--Zhang inputs, or exact least-squares weights are
available.

The two external parameter sources are Li and Zhang (2005), DOI
`10.1016/j.pepi.2005.02.004`, for the added acoustic observations and Xu et al.
(2008), DOI `10.1016/j.epsl.2008.08.012`, for the adopted thermal terms. Their
published coefficients are traceable, but Chantel's exact combined input matrix
and weighting remain unavailable from the Chantel article.

## Seven source-row dispositions

| LitCurate row | Source result | Disposition |
|---|---|---|
| `litcurate_31bc34e8292357ee` | Mg-Pv density BM2, KT=257(2) GPa | **Withheld:** unweighted Table 1 refit gives 243.1 GPa; published curve exceeds the high-P densities by several printed density errors. |
| `litcurate_079361f39d81f154` | Fe-Pv density BM2, KT=246(2) GPa | **Withheld:** Table 1 refit gives 229.3 GPa; coefficient parity is not achieved. |
| `litcurate_542dc7a1f73d9d73` | Mg-Pv acoustic KS=247(4), K'=4.5(2) | **Reproduced as acoustic finite strain:** not added as a volumetric pressure EOS. |
| `litcurate_259fceb6dd4123fc` | Mg-Pv combined acoustic KS=252(1), K'=4.1(1) | **Curve checked but not exactly refittable:** numerical Li--Zhang observations are external. |
| `litcurate_59accf80318ce824` | Fe-Pv acoustic KS=236(2), K'=4.7(1) | **Withheld from this pure-Mg record:** a separate composition and acoustic fit. |
| `litcurate_4c8341508b0d1a6a` | Mg-Pv Table 3 composite thermoelastic model | **Accepted:** `bridgmanite_chantel_2012_bm3_mgd`, with composite provenance and model limits explicit. |
| `litcurate_4de85b1c154b864a` | Fe-Pv Table 3 model | **Withheld:** printed `V0=25.50 cm3/mol` implies about 3.999 g/cm3, irreconcilable with the independently printed 4.161(1) g/cm3 reference density. |

## Fourteen comparison rows

These are all citation-reported literature comparisons and are rejected from
this paper's production count: `litcurate_646cbed7a410ad5d`,
`litcurate_67cb8f3a4139db77`, `litcurate_c6e72c7743e81ed6`,
`litcurate_13f91e0f265496af`, `litcurate_82fa467c469361c3`,
`litcurate_32212ca6614abcc8`, `litcurate_ed458727faf26675`,
`litcurate_25bcc1466760438f`, `litcurate_94fc588faf4261ea`,
`litcurate_e28b7d2629a0bd84`, `litcurate_b0e9cd7bc1f10123`,
`litcurate_781d39e2b88b3db0`, `litcurate_adaf5cd6f8d92989`, and
`litcurate_b2e503f20b0da628`.

## Zotero-ready citation

Chantel, J., D. J. Frost, C. A. McCammon, Z. Jing, and Y. Wang (2012).
“Acoustic velocities of pure and iron-bearing magnesium silicate perovskite
measured to 25 GPa and 1200 K.” *Geophysical Research Letters* 39, L19307.
DOI: `10.1029/2012GL053075`.
