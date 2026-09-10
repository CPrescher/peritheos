# Sakai et al. (2018): rhenium on the Yokoo platinum scale

Sakai, T., Yagi, T., Irifune, T., Kadobayashi, H., Hirao, N., Kunimoto, T.,
Ohfuji, H., Kawaguchi-Imada, S., Ohishi, Y., Tateno, S., and Hirose, K.
(2018), *High pressure generation using double-stage diamond anvil technique:
problems and equations of state of rhenium*, **High Pressure Research 38**(2),
107-119. [DOI: 10.1080/08957959.2018.1448082](https://doi.org/10.1080/08957959.2018.1448082).

The source is the user-supplied publisher PDF `08957959.2018.1448082.pdf`,
containing a cover sheet and 13 article pages. The `.eosmat` source check
records its SHA-256. Page numbers below refer to article pagination; article
page 10 corresponds to final journal page 116. The PDF is retained in the
user's Zotero library rather than redistributed with Peritheos.

## Published record

`rhenium_sakai_2018_yokoo_pt_vinet` extends the existing hcp-Re material
(conventional two-atom cell, space group P63/mmc). It preserves the fit
reported in Section 3, immediately below Figure 10:

| Parameter | Value | Source treatment |
|---|---:|---|
| V0 | 29.47 Å³/conventional cell | Fixed |
| K0 | 358(10) GPa | Fitted |
| K0' | 4.8(2) | Fitted |

The source names the standard Vinet EOS. With `x = (V/V0)^(1/3)`, its mapping is
`P = 3 K0 (1-x) x^-2 exp[1.5 (K0'-1) (1-x)]`. The two quoted errors are
preserved, while confidence level and covariance remain unknown. A fixed V0
does not imply a measured error of zero. The room-temperature experiment is
represented at 300 K by catalog convention; the article does not report an
exact experimental temperature or a thermal EOS.

The fit uses the **platinum Vinet pressure scale of Yokoo et al. (2009)**,
[Physical Review B 80, 104114](https://doi.org/10.1103/PhysRevB.80.104114).
Figure 10 also shows a Dewaele-Pt reduction and previous Re curves; these are
comparison data and do not define additional Sakai records. The exact Yokoo-Pt
record is not currently bundled, so no approximate calibration link is supplied.

## Experimental coverage and observations

Section 2 describes conventional-DAC run RP01 (silica glass, tungsten gasket)
and double-stage run Micro17 (glycerin, direct loading of Re and Pt). Figure 10
shows their distinct symbols. The article does not specify exact membership of
the fitted rows or the residual coordinate and weights.

The approximate pressure envelope **130-310 GPa** is read from the Yokoo-Pt
symbols in Figure 10. It is marked qualitative. The text's 280 GPa RP01 limit
uses the lower pressure scale; the Yokoo-Pt reduction reaches roughly 310 GPa
in the figure. These are experimental coverage estimates, not phase boundaries.

All pages of the attached article were extracted, and the parameter/benchmark
pages were visually inspected. Its only table describes anvil and sample setup;
there is no numerical P-V or paired Pt-Re observation table. The publisher
listing exposes no supplementary dataset; direct full-text web access returned
403. No other author-supplied files were available. The original numerical table
remains unavailable, but Figure 10 supports the qualified graphical refit below.
Its data status remains `plot_only`; its refit outcome is now **`similar`**.

## Independent coefficient comparison

The PDF contains vector marker paths. The extraction retains all **58 marker
centers**: 26 conventional-DAC RP01 observations on each of two platinum scales,
plus six double-stage Micro17 observations on the Yokoo-Pt scale. The two RP01
series are alternative reductions of the same observations, so they must not
be pooled as independent measurements. Separate fill and outline paths for a
single marker are counted once; distinct overlapping markers are retained.

The extractor calibrates PDF coordinates against all nine labeled pressure
ticks and seven labeled volume ticks. Maximum tick residuals are 0.065 GPa and
0.000036 Å³. The checksum-protected CSV stores the original graphic coordinates,
path indices, series labels and transformed P-V values. Legend symbols and
literature curves do not enter the dataset or refits.

With **V0 = 29.47 Å³ fixed**, an unweighted pressure-residual Vinet fit to the
**26 RP01 Yokoo-Pt gray diamonds** gives:

| Parameter | Published | Independent graphical refit | Difference |
|---|---:|---:|---:|
| K0 (GPa) | 358(10) | 357.192 +/- 10.227 | -0.808 (-0.23%) |
| K0' | 4.8(2) | 4.83483 +/- 0.18015 | +0.03483 (+0.73%) |

Both coefficients lie inside the paper's quoted error widths. Pressure RMSE is
**3.106 GPa**, compared with **3.213 GPa** for the published coefficients on
these same digitized rows. Conditional least-squares standard errors are
10.23 GPa and 0.180, with correlation -0.996. These describe residual scatter
under this assumed objective; they do not include pressure-scale or source
measurement uncertainties.

The surrounding Section 3 text discusses RP01 before reporting the fit, and
this selection recovers the published coefficients and their approximate error
widths. It is a supported inference about source-row membership, not an
explicit row list supplied by the authors. All selection sensitivities remain
visible:

| Diagnostic variant | Rows | K0 (GPa) | K0' | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|
| RP01 Yokoo-Pt; pressure residuals | 26 | 357.192 | 4.83483 | 3.106 |
| RP01 + Micro17 Yokoo-Pt; pressure residuals | 32 | 327.075 | 5.37348 | 4.463 |
| RP01 Yokoo-Pt; volume residuals | 26 | 346.771 | 5.02503 | 3.180 |
| RP01 Dewaele-Pt; alternate pressure scale | 26 | 344.408 | 4.64529 | 2.883 |

Adding Micro17 or switching to a volume-residual objective moves the coefficients
outside the individual quoted error widths. Systematic marker-coordinate
offsets of +/-0.1 PDF point give K0 = 355.07-359.33 GPa and K0' = 4.8163-4.8534;
these deliberately assumed perturbations are sensitivity bounds, not measured
uncertainties. The exact PDF paths avoid manual pixel-center choices but cannot
recover precision lost when the authors produced the graph.

The ledger outcome is **`parity`**, based on agreement of both coefficients
within the reported parameter error widths and comparable refit standard errors.
The refit errors above are conditional 1-sigma standard errors. The explicit
`parity_basis = within_reported_parameter_errors` records this numerical criterion;
the unknown confidence convention of the published errors leaves the formal
combined-two-sigma comparison unset. Digitization, inferred row selection,
assumed weighting and the selection sensitivities above remain qualifications.
The published EOS record remains unchanged; these refits are diagnostics.

## Independent pressure checkpoint

Article page 10 reports 464 GPa for the maximum compression of run Micro16.
Page 11 separately gives approximately V = 18.65 Å³ and V/V0 = 0.633.
The printed Vinet coefficients yield **461.5558 GPa**, a difference of
**-2.4442 GPa (0.53%)** from the reported pressure. This difference is retained.

The reproduction evaluates the equation independently and verifies the catalog
against it. It also evaluates half-last-digit rounding bounds for V, V0, K0,
and K0'. The source's 464 GPa lies within those bounds. This establishes
compatibility with the printed precision; it does not prove which unrounded
numbers the authors used, and the bounds are not a statistical error interval.

This maximum-pressure evaluation is an **extrapolation** beyond the Re-Pt
co-compression fit. Reproducing it does not validate experimental calibration
at 464 GPa or independently recover the fitted coefficients.

Run:

```bash
python -m scripts.reproduce_sakai_2018_rhenium
python -m scripts.reproduce_sakai_2018_rhenium --check
pytest -q tests/test_sakai_2018_rhenium.py
```

The [machine-readable result](../data/sakai-2018-rhenium-reproduction.json)
contains the fitted coefficients, residuals, selection/objective sensitivities,
checkpoint, rounding bounds, and scope. To repeat the vector extraction with
`pdfplumber` installed, run `python scripts/extract_sakai_2018_figure10.py
/path/to/08957959.2018.1448082.pdf`. The input checksum is verified before
extracting, and numerical refits subsequently use only the bundled CSV.
