# Grossular: recovery of the remaining nominal channels

Follow-up audit, 12 November 2026. This extends the
[observation audit and conditional 95-channel refit](grossular-angel-2022.md).
No author was contacted, no additional fit was run, and no catalog record,
coefficient, covariance, default or scientific classification was changed.
The existing pyrope email remains **waiting for Ross Angel's reply**.

## What was recovered

**All 52 Gwanmesia modulus rows are recovered from the final publisher table.**
**Seven Thiéblot expansion points are recovered by digitizing Angel's figure.**
These are different evidence levels: the latter are already processed plotted
points, not the original Thiéblot table or experimental uncertainties.

| Evidence | Nominal S1 channels covered | Remaining qualification |
| --- | ---: | --- |
| Earlier tabulated or table-derived observations | 161 | Original selections, anchors and some errors incomplete |
| Gwanmesia 2014, final Table 1 | 52 | All 12 numerical columns retained; error confidence/covariance unspecified |
| Thiéblot 1998, Angel 2022 Figure 4 | 7 | Digitized processed points; original selected table values/errors unavailable |
| **Total at mixed evidence levels** | **220** | **Not 220 recovered original fit inputs or a reproduced full fit** |

The [follow-up receipt](../data/garnet-angel-2022/grossular/followup-recovery.json)
records numerical comparisons and input hashes. The earlier 161-channel
recovery, observation audit and 95-channel fit receipts remain intact as
historical results. Their missing-input statements describe that earlier stage.

## Gwanmesia: complete primary table and fixed-source check

Source: Gwanmesia et al. (2014), *Elasticity and sound velocities of
polycrystalline grossular garnet (Ca3Al2Si3O12) at simultaneous high pressures
and high temperatures*, PEPI 228, 80–87,
[final primary publication](https://doi.org/10.1016/j.pepi.2013.09.010).
Institutional access exposed the complete publisher HTML Table 1. The PDF
download did not complete; this recovery uses the HTML table.

The [transcription](../data/garnet-angel-2022/grossular/gwanmesia-2014-table-1.csv)
preserves all **624 numerical cells**, including trailing zeros. Its body
fingerprint matches a separate extraction from the publisher's displayed table.
The [source receipt](../data/garnet-angel-2022/grossular/gwanmesia-2014-source.json)
retains conversions, calibration, error statements and dependencies.

The ten room-temperature rows are **30 °C = 303.15 K**. Seven further groups
at 127–727 °C contain six rows each. Pressures span **0.73–9.76 GPa**, using
NaCl with the Decker (1971) scale. No numerical temperature error is supplied.
The source reports typical pressure precision of 0.7% and a Kₛ error estimate
of 0.8%; these are not automatically independent one-sigma uncertainties.
Caption and later text give different density error estimates (0.2% and 0.1%);
both are retained. The largest volume uncertainty is not assigned to every row.

Using the unchanged Angel coefficients, independent Debye quadrature and the
published **calculated-modulus multiplier 1.00551**, the 52 Kₛ residuals have:

| Observed minus calculated | RMS, GPa | Mean, GPa | Maximum absolute, GPa |
| --- | ---: | ---: | ---: |
| With published multiplier | 0.71176 | −0.00241 | 1.62394 |
| Without multiplier, diagnostic comparison | 1.25892 | 1.04717 | 2.03574 |

The multiplier is taken directly from source parameter 54, **not fitted here**.
This supports consistency of the recovered table with the published scaled
modulus curve at this comparison scope. It is not a full goodness-of-fit test:
pressure errors and shared calibration covariance are not propagated here.
The printed identity Kₛ = L − 4G/3 agrees within 0.13334 GPa, compatible with
rounding all three modulus columns to 0.1 GPa. No table value is repaired.

Lengths, density, velocities and moduli share volume and transit-time inputs.
The 52 XRD volumes remain excluded, following S1; they are not additional
independent constraints alongside these 52 Kₛ values.

## Isaak: specimen identity traced to the plotted series

[Isaak et al. (1992)](https://doi.org/10.1007/BF00198608), Table 2, measures
expansion of **GR-2B: 76.1 mol% grossular, 22.1 mol% andradite**.
The **GGD, about 96.53 mol% grossular**, composition in Angel's S1 belongs to
Isaak's elasticity specimen instead. Isaak p.109 explicitly defines the strain
reference as **293 K**, resolving the earlier audit's reference-temperature
assumption. The original data were interpolated onto the published grid.

We extracted all 15 red markers from Angel Figure 4. With **no fitted scale**,
the expression

\[
V_{\mathrm{plot}}=125.35\,[1+0.01\,Y_{V,\mathrm{GR\text{-}2B}}]
\]

matches the observed GR-2B Table 2 strain column with RMS **0.0007085 cm³/mol**
and maximum difference **0.001239 cm³/mol**. Both are below the declared
digitization bound of 0.003 cm³/mol. Temperature differences are below 0.735 K,
within the 2 K digitization bound.

**The plotted series is therefore consistent with GR-2B strain directly scaled
to the adopted volume, despite S1's GGD composition label.** At this resolution,
the nearly identical calculated GR-2B strain column also agrees; digitization
cannot uniquely distinguish those two columns. No resolved compositional
correction or 293-to-298 K re-referencing is evident from this comparison.
The last Table 2 column estimates end-member expansivity through a compositional
assumption; it is not measured pure-grossular volume.

This resolves the origin of the plotted expansion series at the stated
resolution. It does **not** establish that the fit-project input was identical
to the plotted series. The 15 mixed-composition expansion rows remain outside
the diagnostic end-member refit.

## Thiéblot: seven plotted points, original table still missing

The [1998 primary article](https://doi.org/10.1127/ejm/10/1/0007) remained
inaccessible: GeoScienceWorld recognized the institutional connection but
explicitly denied full-text access; Schweizerbart supplied only a one-page
preview. Local library and repository checks did not recover the full table.

Angel Figure 4 supplies **seven filled magenta markers below 1000 K**. Their
[digitized coordinates](../data/garnet-angel-2022/grossular/figure4-thieblot-selected.csv)
cover approximately 300–975 K and 125.356–127.586 cm³/mol. Six visible open
markers above 1000 K are retained in a
[separate excluded series](../data/garnet-angel-2022/grossular/figure4-thieblot-excluded-visible.csv).
These 13 visible markers do not recover all 15 original observations.

The [digitization receipt](../data/garnet-angel-2022/grossular/figure4-recovery.json)
records the exact source PDF hash, vector paths, axis calibration and conservative
bounds: **2 K and 0.003 cm³/mol**. Bounds describe digitization resolution,
not experimental errors or statistical confidence. Extra decimal places in
the CSV preserve the extraction calculation, not experimental precision.
The source figure has already applied volume processing; it cannot recover
the original absolute measurements, reference anchor or covariance.

There is also a **source selection conflict**. The main article p.8 says the
Skinner, Thiéblot and Rodehorst (2004) expansion sets are excluded. S1 and the
official EOS comments retain seven Thiéblot points below 1000 K, consistent
with the seven filled markers. The exact fit mask remains unconfirmed.

## What still prevents a full author-workflow reproduction

The remaining gaps concern original fit inputs and reduction, rather than
52 missing Gwanmesia observations: original Thiéblot table values/errors;
confirmation of the Isaak input and compositional treatment; exact row masks,
normalization anchors and induced covariance; assigned errors where absent
from the primary tables; paired quartz readings; and fit controls. S1's
220-channel sum still conflicts with its ambient-anchor counting footnote.

We have not created errors to force a nominal 220-channel fit. The existing
95-channel refit remains a conditional diagnostic, with its original receipt.
The follow-up adds a fixed-source check of all 52 recovered Gwanmesia moduli.
No grossular request has been sent to Ross.

These are **repository-only maintenance commands**:

```bash
# Requires optional PyMuPDF and a privately obtained final source PDF:
python -m scripts.recover_grossular_angel_2022_figure4 --pdf /path/to/article.pdf
python -m scripts.audit_grossular_angel_2022_recovery
pytest -q tests/test_grossular_followup_recovery.py
```

The [follow-up verification receipt](../data/garnet-angel-2022/grossular/followup-verification.json)
records the completed checks. Numerical observations are attributed; no
publisher figures or full text are redistributed.
