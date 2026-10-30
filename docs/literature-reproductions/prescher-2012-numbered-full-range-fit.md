# Prescher p1-p19: full-range BM3 refit

Source: [Prescher et al. (2012), PRB 85, 140402(R)](https://doi.org/10.1103/PhysRevB.85.140402), Figure 1 and final supplementary Table S1.

All 19 numbered points use their final Table S1 volumes and published ruby pressures (1.8-47.4 GPa, nominal 293 K). The source p16 volume is 133.25 A^3. The ambient candidate is outside this selection.

The primary fit uses both printed pressure and volume widths, conditionally as independent absolute 1-sigma inputs. V0, K0 and K0-prime are free. Confidence/covariance of the source widths is unknown; ruby widths are before/after reading differences.

| Fit | V0 (A^3) | K0 (GPa) | K0-prime |
| --- | ---: | ---: | ---: |
| Published full range | not printed | 161 +/- 2 | 5.9 +/- 0.2 |
| Numbered points, coordinate widths, V0 free (primary) | 154.48775 +/- 0.11678 | 160.53550 +/- 3.46985 | 6.03174 +/- 0.22887 |
| Numbered points, coordinate widths, adopted V0 fixed (control) | 154.56788637 fixed | 158.35935 +/- 1.42698 | 6.15883 +/- 0.14135 |
| Numbered points, equal pressure weights, V0 free (control) | 154.19436 +/- 0.24776 | 170.27611 +/- 6.77416 | 5.44801 +/- 0.37894 |
| Numbered points, equal pressure weights, adopted V0 fixed (control) | 154.56788637 fixed | 160.99401 +/- 2.53219 | 5.91948 +/- 0.24087 |

At the source's printed digits, the primary fit is K0=161 +/- 3 GPa and K0-prime=6.0 +/- 0.2. The fixed-V0 equal-pressure control gives 161 +/- 3 GPa and 5.9 +/- 0.2: both central coefficients match the paper's digits, but error(K0)=2.53219 GPa rounds directly to 3, not the source's 2. Do not double-round intermediate error displays.

Primary weighted SSE/(19-3) = 3.34023. Residual-scaled formal errors are V0 +/- 0.21344 A^3, K0 +/- 6.34160 GPa and K0-prime +/- 0.41828; these use a different error convention from the primary comparison.

Fixed-V0 errors condition on exact adopted curve V0 and exclude its uncertainty. The source does not print V0 or recoverably specify whether it was fixed. The old 22-marker parity comparison uses digitized graph coordinates and equal pressure weights; it remains a separate conditional comparison.

Independent ODRPACK and linear-pressure checks pass. Numerical comparison does not reproduce the author's input set, weighting, covariance or confidence convention. Published coefficients and observations remain unchanged. This Studio snapshot adopts the fixed-V0 equal-pressure control as its recorded coefficient comparison, retaining the K0 uncertainty mismatch.

Complete numeric evidence: `docs/data/prescher-2012-numbered-full-range-fit.json`.

Repository-only replay: `.venv/bin/python -m scripts.check_prescher_2012_numbered_full_range`.
