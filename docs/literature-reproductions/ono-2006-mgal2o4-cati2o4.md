# Ono et al. (2006): CaTi2O4-type MgAl2O4

Audit date: 2026-09-08. **Accepted:** one published 300 K BM2 record,
`mgal2o4_cati2o4_ono_2006_bm2_2`, in the existing Cmcm material.

## Primary authority and identity

Ono, S., Kikegawa, T., and Ohishi, Y. (2006), “The stability and
compressibility of MgAl2O4 high-pressure polymorphs,” *Physics and Chemistry
of Minerals* 33, 200–206, [doi:10.1007/s00269-006-0068-z](https://doi.org/10.1007/s00269-006-0068-z).
The [publisher landing page](https://link.springer.com/article/10.1007/s00269-006-0068-z)
confirms the citation and abstract. The complete seven-page version of record
was read from this [repository PDF](https://repository.geologyscience.ru/server/api/core/bitstreams/55584c23-e366-4fc8-ae97-f94d81827cbf/content),
with visual verification of the equation and Tables 2–3 on pp. 203–204.
No supplement is linked by the publisher preview; Table 3 supplies the primary
compression observations. The article PDF is not redistributed in the bundle.

The target is pure MgAl2O4 with conventional C-centered orthorhombic Cmcm
structure, Z=4, not the natural MORB aluminous phase. The existing Funamori
(1998) structure and default EOS remain separately attributed. Ono's fitted
reference volume does not replace that recovered structural cell. No Mg/Al/O
coordinates are refined in this paper; none are invented.

## Equation, parameters, and scope

The unnumbered equation on p. 203 is
`P = 1.5 K0 (x^-7 - x^-5) [1 + 0.75 (K0' - 4) (x^-2 - 1)]`,
where `x = (V/V0)^(1/3)`. Fixing `K0'=4` gives exactly BM2.

| Quantity | Published value | Treatment |
|---|---|---|
| V0 | 238.9(9) Å³/conventional cell | Fitted; error 0.9 Å³ |
| K0 | 219(6) GPa | Fitted; error 6 GPa |
| K0′ | 4 | Fixed by BM2, not a free coefficient |
| Temperature | 300 K | Table 3 and Figure 4 |

The abstract prints `238.9(±9)`, while the body uses `238.9(9)`; the
last-digit body notation is the basis for the stored 0.9 Å³ error. Confidence
level, covariance, objective, and fitting weights are not specified.

Table 3 has two recovered 0 GPa cells and twelve measurements at
42.7–91.2 GPa. Its 14 rows, all axes, volumes, and printed errors are retained
in the checksummed CSV. Parenthetical errors are called only “error” in the
footnote, so the typed columns use `uncertainty`, not an asserted standard
deviation. Printed `(0)` entries remain zero-rounded values; they do not
supply infinite statistical weights. The separately printed volumes are used
without replacing them by products of rounded axes.

The high-temperature phase-stability range (approximately 45–117 GPa) and
Table 1's heating temperatures are not the validity limits of this isotherm.
The ambient cells are recovered metastable material, not evidence of ambient
equilibrium stability. No thermal EOS or epsilon-phase EOS is inferred.

Pressures use Holmes et al. (1989) Pt, linked to the existing executable
`platinum_holmes_1989_vinet_1`. The electronic thermal-pressure correction
mentioned on p. 201 concerns heated experiments. Table 3 supplies no paired
Pt lattice measurements or row pressure errors; the cited typical 2–3 GPa
high-P/T error is not assigned to every 300 K row.

## Independent reproduction

Run `python scripts/reproduce_ono_2006_mgal2o4.py`. This implements the printed
equation independently of Peritheos and minimizes unweighted pressure
residuals, with both V0 and K0 free and K0′ fixed implicitly.

| Diagnostic | V0 (Å³) | K0 (GPa) | Pressure RMSE (GPa) |
|---|---:|---:|---:|
| Published curve, all 14 rows | 238.9 | 219 | 1.326248 |
| Independent all-row fit | 238.880747 | 219.046970 | 1.326064 |
| Sensitivity: twelve nonzero-pressure rows | 235.630226 | 241.533861 | 1.124259 |

The all-row coefficients recover the published fit well within its reported
errors. This supports including both recovered cells, but does not prove the
authors' undocumented weighting or selection. Excluding them substantially
changes the extrapolated coefficients, despite a smaller pressure residual.
Neither diagnostic replaces the published record or creates a new default.

At the independently printed Table 3 endpoint `V=186.4 Å³, P=91.2 GPa`, the
published equation gives 89.370518 GPa (residual −1.829482 GPa), the largest
absolute table residual. This discrepancy is retained; the point is not
claimed to lie within an unavailable pressure confidence interval. Tests
compare the public evaluator against the independent equation across the
entire table, verify coefficient recovery, and check pressure-volume inversion.

## Candidate and citation dispositions

The discovery ledger `docs/data/litcurate-eos-candidates.json` contains no row
with this publication DOI or title as of this audit. This is a direct primary
literature addition; no LitCurate identifier or source row is fabricated.
The accepted direct lead is registered in `docs/material-eos-candidates.md`.

| Source item | Disposition |
|---|---|
| Results p. 203, Table 3, Figure 4: this-study Cmcm MgAl2O4 | **Accept** one 300 K BM2 |
| Epsilon phase, pp. 202–204 | **No EOS**; structure unresolved and no fit coefficients |
| Catti (2001) and Gracia et al. (2002) theoretical moduli | **Citation traces only**; not Ono-owned fits |
| Figure 5 spinel, CaFe2O4-type, and oxide curves | **Citation traces only** to the named primary papers |
| Figure 5 MORB `K0=184(8), V0=243.0(13), K0'=4` | **Exclude from this record**; different composition and attribution to Ono et al. (2005a) |
