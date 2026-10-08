# Tateno-Holmes and Campbell-Heinz: primary-source gap audit

Audit date: 2026-10-08. Run `python -m scripts.audit_kcl_source_gaps --check`.
The [numerical report](../data/kcl-source-gaps-audit.json) keeps the conditional
regressions separate from published EOS records. Published coefficients,
defaults, all original CSV bytes and the 110 earlier Pt pairing corrections
are unchanged.

## Campbell and Heinz (1991): calibration identified

The complete final [Campbell-Heinz article](https://doi.org/10.1016/0022-3697(91)90181-X)
was recovered from an existing local publisher PDF. Both local copies have
SHA256 `bba2946d8398c43dc5342f48073426d446a61508e0211fb72de1cae776eec2c5`.
The author publications page links to the DOI; the publisher web route exposed
only an abstract/access challenge during this audit. The PDF is not redistributed.

The methods on pages 495-496 explicitly assign pressure to ruby fluorescence
using reference 5, [Mao et al. (1978)](https://doi.org/10.1063/1.325277), identified
on page 499. This resolves the scale as `ruby_mao_1978`, whose bundled power law
has A=1904 GPa and B=5. No pressure-transmitting medium was used. Gold at zero
pressure calibrated the camera distance; it was not the pressure standard.

Each pressure is the mean of five spatial readings. The quoted uncertainty is
their standard deviation, not the standard error of the mean. Volume errors are
standard deviations across diffraction-line volumes. The dataset metadata now
identifies these meanings and declares all 14 pressures `as_reported` on Mao
(1978); `resolve_dataset_pressure` returns them with verified row identities
and checksum. Parameter-error confidence remains unknown.

A later nonlinear ruby conversion can normalize the printed mean pressure,
but cannot recover the mean of the five separately converted readings. The
original spatial readings and raw/reference wavelengths are absent. Source
standard deviations are retained without assigning them to another scale.
The stored 298 K reference is unchanged; the article says room temperature.

The source Table 1 basis is B1 V01=37.521 cm3/mol, and the methods adopt B1
lattice a0=6.2931 A. These printed source quantities remain separate from the
executable record's documented Campbell-Dewaele composite, V01=62.36 A3 per
formula unit. The composite V02 and its uncertainty are unchanged.

## Campbell source regression: conditional numerical recovery

Page 496 specifies weighted linear regression of normalized stress
against Jeanloz effective strain:

```text
r = V2/V01
 g = (r^(-2/3)-1)/2
 G = P/[3*(1+2*g)^(5/2)]
 G = b + m*g
r0 = V02/V01 = (1-2*b/m)^(-3/2)
K02 = m/r0^(7/3), with K02'=4 fixed
```

The report compares three explicit response-weight choices on all 14 printed
rows. Pressure-only response errors use sigmaG=sigmaP/[3*(1+2g)^(5/2)]. The
combined choice adds the volume contribution (5/3)*G*sigmaR/r in quadrature.
Neither choice accounts for error in g or its shared-volume covariance with G.

| Diagnostic weights | V02/V01 | K02 (GPa) | Within both printed error widths |
|---|---:|---:|---|
| Equal G weights | 0.83458188 | 31.06987 | No |
| Pressure-only G response errors | 0.84832774 | 28.75442 | Yes |
| Pressure plus volume G response errors | 0.84833045 | 28.72280 | Yes |
| Published | 0.8483 +/- 0.0057 | 28.7 +/- 0.6 | Source values |

This recovers the reported coefficients under plausible stated assumptions.
It does not establish the exact author's weights, covariance, or uncertainty
calculation. No new EOS record is generated, and the catalog's generic
pressure-residual refit must not be mistaken for this original regression.

## Tateno (2019): final methods and full official archive checked

The final [Tateno article](https://doi.org/10.2138/am-2019-6779) has the previously
recorded SHA256 `b53a139b2c068d31e7869464e9afec10ebb6a709945131ee10a4913da2997db9`.
Pages 720-722 specify simultaneous elastic and thermal fitting, fixed V0=54.5 A3,
fixed theta0=235 K for MGD, and the two Table 1 Holmes alternatives. They cite
Holmes (1989) without specifying which thermodynamic implementation was used.
The accepted manuscript at Edinburgh remains a discovery route; its superseded
thermal coefficients are not used.

The complete [official MSA AM-19-56779 archive](http://www.minsocam.org/MSA/AmMin/TOC/2019/May2019_data/AM-19-56779.zip)
was retrieved and verified against the existing archive hash
`652fa6bcbd9c00c76d8bb62199750a49f281edf3dc445144ff435ff1afce4210`.
Its substantive contents are:

- `6779TableS1 revised.xlsx`: the unchanged bundled 39-row source workbook.
- `6779_supp.pdf`: one page containing Figure S1, the unrolled diffraction-image
  comparison after annealing and following compression. It has no additional
  thermal equations or regression data.

The remaining archive members are filesystem metadata and an Excel lock file.
The source manifest at `peritheos/data/datasets/kcl_variant_sources/manifest.json`
now records each member's size and checksum. No article text or figures are bundled.
The workbook has one sheet, no formula cells and no external links. It supplies
paired Sokolova P, T, Pt/KCl volumes and source errors, with no Holmes-pressure
column, implementation, numerical fit weights or covariance.

## Tateno coordinate precision and thermal limitation

The [existing reproduction](kcl-published-variants.md) retains all 39 raw
Sokolova-pressure rows and separate conditional Holmes coordinates calculated
from the paired Pt volumes with Holmes Equations 11-12. The full
[Holmes audit](holmes-1989-platinum.md) distinguishes the exact Equation 11 curve
from its approximate Equation 12 thermal extension and the incompletely
deposited ionic/electronic thermodynamic construction.

Using the corrected two-decimal Pt-volume transcription instead of the full
workbook Pt decimals changes the derived Holmes coordinate by at most
0.12374833 GPa. This check holds T and KCl volumes at workbook precision and
uses equal pressure weights in both fits. The MGD diagnostic changes K0 from
17.49251 to 17.47604 GPa; the linear diagnostic changes it from 17.55493 to
17.53781 GPa. The result demonstrates sensitivity to available precision; it
cannot identify the author's actual regression precision.

Seven rows, Excel rows 25, 27, 30, 32, 34, 47 and 49, have T>=2000 K and exceed
Holmes's stated domain for the approximate thermal term. Tateno reaches 2560 K.
Neither the final methods nor the complete inspected deposit identifies a
higher-temperature correction or full Pt thermodynamic implementation.
No extension of Holmes's thermodynamic accuracy is inferred from coefficient
agreement. No Sokolova pressure uncertainty is assigned to Holmes, and no
complete derived-pressure uncertainty is claimed.

## Concrete follow-ups requiring source evidence

No authors were contacted. The remaining questions are recorded for a future
source request:

1. Tateno: provide the actual Holmes P(V,T) code/worksheet and coefficients;
   identify whether Equation 12, the full ionic/electronic model, or another
   thermal correction was used, particularly for the seven T>=2000 K states.
2. Tateno: provide the row-aligned Holmes pressures, the exact Pt/KCl/T numerical
   precision used in regression, residual variable, numerical weights,
   coordinate-covariance treatment, fit covariance and parameter-error meaning.
3. Campbell-Heinz: provide the five original spatial ruby readings/reference
   wavelengths per state and the numerical G-versus-g weights, shared-volume
   covariance treatment and parameter-error calculation.

The Campbell scale identity is closed. Exact statistical reproduction and raw
ruby recovery remain unavailable. Tateno's exact Holmes reduction and original
regression remain conditional after checking all identified final-source files.
