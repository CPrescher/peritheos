# Hirose et al. (2005): MgGeO3 post-perovskite

## Outcome

Both same-paper LitCurate candidates are accepted.  They are the source's
preferred fixed-`K0'=4` BM3 fit and its explicit free three-parameter BM3
alternative to the same nine observations.  Neither duplicates the later Kubo
et al. (2006) experiment, which used a different high-pressure selection and
fit strategy.

Kei Hirose, Katsuyuki Kawamura, Yasuo Ohishi, Shigehiko Tateno, and Nagayoshi
Sata, “Stability and equation of state of MgGeO3 post-perovskite phase,”
*American Mineralogist* **90**, 262-265 (2005),
[doi:10.2138/am.2005.1702](https://doi.org/10.2138/am.2005.1702).
The audit used the RRUFF archival final PDF, SHA-256
`a0c0584d4d82a731b7912af429b37a93cf6c84b88f21b72181f21bd0be609ae9`.

## Data and fits

Table 3 prints all nine 300 K decompression observations from 78.55 to 6.25
GPa, including pressure, a, b, c, conventional-cell volume, and parenthetical
uncertainties.  They are transcribed without digitization.  Table 1 and the
structure discussion identify Cmcm with Z=4.

| LitCurate ID | Source result | Disposition |
|---|---|---|
| `litcurate_ff36511f92d6a719` | V0=182.2(11) A3, K0=210(20) GPa, K0'=3.5(5) | ACCEPT: free BM3 sensitivity fit. |
| `litcurate_c8a702f01054104b` | V0=183.1(8) A3, K0=192(5) GPa, K0'=4 fixed | ACCEPT: abstract-selected constrained BM3. |

The article explicitly prints the Birch-Murnaghan pressure expression.  The
fixed-derivative fit is presented in the abstract and is therefore listed
first among the Hirose records.  The free fit remains executable because it is
a complete source result and exposes the substantial V0-K0-K0' tradeoff.

Unweighted fixed-derivative refitting gives V0=183.0161 A3 and K0=192.6917
GPa.  A pressure-uncertainty-weighted free refit gives V0=182.1746 A3,
K0=212.6214 GPa, and K0'=3.4632.  Both lie within the broad source uncertainty
regions.  Published-curve pressure RMSE is 0.872 and 0.897 GPa respectively;
the largest residual is the isolated 72.27 GPa row, and no point was silently
removed.

Pressures use Pt and the Holmes et al. (1989) EOS,
[doi:10.1063/1.344177](https://doi.org/10.1063/1.344177).  The paper does not
print the row-wise Pt volumes, so the pressure reduction itself cannot be
recalculated.  The source reports partial back-transformation during the final
6 GPa step; the records are conservatively valid only over the direct
6.25-78.55 GPa range.
