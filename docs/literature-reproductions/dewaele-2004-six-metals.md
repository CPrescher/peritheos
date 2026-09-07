# Dewaele et al. (2004): six-metal compression and ruby-scale audit

## Scope and primary evidence

This audit covers the eight Dewaele et al. records added by the experimental
metal batch and the four older revised-ruby records that share the same source
observations. The scientific authority is Dewaele, Loubeyre, and Mezouar,
“Equations of state of six metals above 94 GPa,” *Physical Review B* **70**,
094112 (2004), [doi:10.1103/PhysRevB.70.094112](https://doi.org/10.1103/PhysRevB.70.094112).
The [publisher PDF](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.70.094112/fulltext)
used for the audit has SHA-256
`1f9e02e4a77f49a868b1a6261290465bdb00597b1660691b1790c7502ec9aeeb`.

- Table I (journal page 094112-4) prints every pressure-volume observation,
  split into the first two Cu/W/Al runs and last three Au/Pt/Ta runs.
- Equation (1) defines the revised ruby scale.
- The Figure 2 caption defines the Vinet pressure equation and
  `x=(V/V0)^(1/3)` convention.
- Table II prints the two pressure-scale fits, fixed coefficients in bold, and
  95% confidence intervals on fitted coefficients.

The article and APS reuse pages do not attach an open data license to this 2004
table. APS copyright applies to the article. The six CSVs are direct factual
numerical transcriptions for scientific validation; their metadata records that
they are not asserted to be openly licensed. No table image, prose, or layout is
redistributed.

## Observation transcriptions and normalization

The CSVs retain source row order and printed decimal precision. Blank Table I
cells are omitted rather than inferred: Pt is blank at the `15.04/15.2 GPa`
state and Ta is blank at `78.0/80.7 GPa`. Every retained row carries the common
printed atomic-volume uncertainty of `0.01 A^3/atom`. The paper gives only a
global pressure-uncertainty description—`0.05 GPa` at `1 GPa` increasing to
`2 GPa` at `150 GPa`—so no row-wise pressure-error column is invented.

| Material | Rows | CSV SHA-256 | Atomic-to-card volume factor |
| --- | ---: | --- | ---: |
| Al | 40 | `61fa04908b3a4827fb6566e4ceea29706ac52cbaae1cee5bb0590d1af44423b7` | 4 (fcc) |
| Cu | 42 | `0f4535c7e29960690c098408711f63ae559be1d73bba18da73047d75f266120e` | 4 (fcc) |
| Au | 37 | `c3a52a211b94cef2cd46198577e18efe4fc66befd7d38aa4f20063d290e97438` | 4 (fcc) |
| Pt | 36 | `b0230e9fdd1f8a36cecc35cb79a6d56bd4ee170007a465f1842ef76544c47273` | 4 (fcc) |
| Ta | 36 | `3c69e1a9943dc1f882d60f63d01cf6708c114a4569a20e04a8b2278e5c09bd78` | 2 (bcc) |
| W | 42 | `d73f1d9114c974be80175c6ade256bf17cb4a529e911480e6cd1e337394ef523` | 2 (bcc) |

Table II gives atomic `V0`. Each `.eosmat` record stores conventional-cell
volume, so the exact fcc or bcc factor above is applied. The Vinet compression
ratio is unchanged by this normalization.

## Pressure reductions

Both Table I pressure columns are preserved as source values. They use

`P = (A/B) * ((lambda/lambda0)^B - 1)`,

with `A=1904 GPa`, `B=7.665` for the classical Mao et al. (1986) `P_R` scale,
and `B=9.5` for the revised Dewaele et al. `P'_R` scale. The source does not
print raw R1 wavelengths. Inverting a rounded `P_R` value and applying the
revised equation reproduces the paired `P'_R` column within `0.54 GPa`; most of
that bound is the explicitly retained `37.0 -> 37.1 GPa` printed pair, for which
the equations give `37.6363 GPa`. This is an internal source-table inconsistency,
not a transcription error. The audit therefore fits each printed pressure
column directly and does not silently replace any value with an equation-derived
one.

The exact observed ranges are `0-144.3/153.0 GPa` on the classical/revised
scales for Cu, W, and Al, and `0-90.0/93.6 GPa` for Au, Pt, and Ta. This corrects
the older Au record's unsupported `154 GPa` range. All records use the paper's
`298 K` reference temperature.

## Fitting choices and independent refits

The source first determines `V0` from `0-5 GPa` measurements and then fixes it
for the whole-data Vinet fits. The six classical-scale records vary `K0` and
`K0'`. The added revised-scale Pt and Ta records, and the older revised-scale Au
record, reproduce the Table II alternatives with both `V0` and the bold `K0`
fixed; the older revised Al, Cu, and W records vary `K0` and `K0'`. The source
does not report numerical regression weights or covariance.

The independent check selects every available row, converts atomic volume to
the card basis, and uses errors-in-variables least squares with the common
volume uncertainty and no invented row-wise pressure sigma. `absolute_sigma`
is true. All twelve represented fits achieve uncertainty parity: every refitted
coefficient lies within the combined two-sigma interval of the published and
refit uncertainties.

| Record | Scale; free coefficients | Rows | Published -> refit | Published/refit RMSE (GPa) |
| --- | --- | ---: | --- | ---: |
| `aluminum_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 40 | `76.3,4.16 -> 76.3832,4.15118` | `0.416615/0.373463` |
| `copper_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 42 | `135.1,4.91 -> 135.4938,4.89854` | `0.384285/0.284969` |
| `gold_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 37 | `172.5,5.40 -> 172.3775,5.41110` | `0.226794/0.203359` |
| `platinum_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 36 | `275.3,4.78 -> 275.0509,4.78955` | `0.216788/0.183973` |
| `tantalum_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 36 | `198.2,3.07 -> 195.2784,3.24406` | `0.387667/0.264718` |
| `tungsten_dewaele_2004_mao_ruby_vinet` | classical; `K0,K0'` | 42 | `298.3,3.81 -> 298.6920,3.80193` | `0.578538/0.476182` |
| `platinum_dewaele_2004_revised_ruby_vinet` | revised; `K0'` | 36 | `5.08 -> 5.08448` | `0.261721/0.220957` |
| `tantalum_dewaele_2004_revised_ruby_vinet` | revised; `K0'` | 36 | `3.52 -> 3.57078` | `0.394586/0.282194` |

The four older revised-scale records remain independently checked in the global
primary-EOS refit ledger: Al (`74.3472,4.46406`), Cu
(`132.7865,5.30404`), fixed-`K0` Au (`K0'=5.99904`), and W
(`295.6709,4.30545`). Their observation counts are 40, 42, 37, and 42,
respectively.

The refits validate the represented coefficients despite the source's omitted
weights and rounded observations; they do not claim to recover an unpublished
bit-for-bit regression implementation.
