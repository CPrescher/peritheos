# Ono and Mibe (2010): original experimental Fe3C EOS

Audited and integrated **1 November 2026** as
`cementite_ono_mibe_2010_bm3`, a **nondefault published** record on the existing
`cementite` material. Use is limited to **300 K and 0–34.3 GPa**. The original
experimental coefficients are retained; the ambient volume is explicitly adopted
from the measured cell, and the original fitting constraints and weights are
unknown. The 0 K calculations and later Litasov recalibration remain separate.

## Source and backlog decision

The authority is Ono and Mibe, *Magnetic transition of iron carbide at high
pressures*, *Physics of the Earth and Planetary Interiors* **180**, 1–6,
[DOI 10.1016/j.pepi.2010.03.008](https://doi.org/10.1016/j.pepi.2010.03.008).
The six-page final Elsevier PDF was recovered from the local Iron Carbides
project. Tables 1 and 3, Figure 4 and the pressure equation were rendered and
visually checked. Its version/checksum and source paths are retained in
`peritheos/data/datasets/cementite-ono-mibe-2010-source.json`. Publisher files
and figures remain local evidence and are not redistributed with the dataset.

The live candidate queue and the cementite audit were checked against the
installed catalog, reproduction scripts, scientific ledgers and completed
investigations. No original Ono–Mibe record or completed original-source audit
was present. Litasov's completed investigation contained later comparison
values only. The repository issue searches for Ono/Mibe and open catalog work
returned no matching active investigation. This addition therefore did not
reopen the completed Scott, Litasov, Sata or Prescher investigations.

## Identity, equation and coefficients

The experiment concerns stoichiometric orthorhombic Pnma Fe3C, with **Z=4**,
Fe12C4 and 16 atoms per conventional cell (p. 2). Table 1 uses conventional-cell
Å³ directly. The existing material structure and default are preserved; neither
is relabelled as an Ono–Mibe structural refinement.

Figure 4 explicitly identifies the experimental fit as **third-order
Birch–Murnaghan below 35 GPa**. The pressure expression on p. 4 is

\[
P=\frac32K_0(x^{-7}-x^{-5})
\left[1+\frac34(K'_0-4)(x^{-2}-1)\right],\qquad
x=(V/V_0)^{1/3},
\]

which maps exactly to native BM3. The Table 3 experimental row, footnote e,
gives **K0=167 GPa, K′0=6.7**. Neither has a reported error. The stored
**V0=155.4(2) Å³** is adopted from the measured ambient Table 1 cell, not from
Litasov, DFT or an independent refit. The source does not explicitly state
whether V0 was fixed in fitting. Consequently, the record does **not** claim
`fixed=true`; the separately declared diagnostic holds that measured anchor
fixed. Its 0.2 Å³ width remains a measured-cell error with unspecified confidence.
The Figure 4 curve check independently supports this reference-volume choice
at graphical resolution.

| Branch | Original source coefficients | Disposition |
| --- | --- | --- |
| Experimental, 300 K, below 35 GPa | V0 adopted 155.4(2) Å³; K0=167 GPa; K′0=6.7 | Installed original BM3, measured limit 34.3 GPa |
| DFT P–V, ferromagnetic, 0 K, ≤35 GPa | K0=216.5 GPa; K′0=4.15 | Separate theoretical result; not installed by this experimental addition |
| DFT E–V, ferromagnetic, 0 K, ≤35 GPa | K0=227.7 GPa; K′0=3.36 | Separate fit objective and coefficients; not merged with P–V |
| DFT P–V, nonmagnetic, 0 K, ≤400 GPa | K0=315.5 GPa; K′0=4.37 | Separate theoretical branch; not a 300 K experimental EOS |

The original theoretical reference volumes are not numerically specified in
Table 3. This addition does not infer precise V0 values from later comparison
tables or invent experimental rows from theoretical curves. The E–V/P–V
coefficient difference is preserved as a source distinction.

Litasov (2013) Table 2 later gives an experimental **Vinet** comparison with
V0=155.4 Å³, K0=166(6) GPa and K′0=6.7(3). That is a different fit/pressure-scale
provenance; its values and errors do not replace the original **BM3 167/6.7**.
See the completed [Litasov audit](litasov-2013-fe3c.md).

## Observations, limits and calibration

All **14 measured Table 1 rows**, 0–66.5 GPa, are retained in the native
`cementite_ono_mibe_2010_table1` dataset, with original tokens in the packaged
CSV. This includes every a/b/c/V value and parenthetical width. The printed
volumes are preserved rather than replaced by products of rounded axes.
Pressure errors are not reported: missing widths remain **null**, not zero.
The footnote calls the parenthetical values errors but specifies no confidence
level or covariance. No 1σ/2σ convention is invented.

The source low-pressure fit contains the ambient state and five nonambient
states: **14.6, 19.7, 24.7, 28.0 and 34.3 GPa**. The ambient row supplies the
reference anchor in the conditional diagnostic. All eight rows above 35 GPa
remain available as observations and are explicitly outside this installed
curve's validity. Table 3's nominal ≤35 GPa limit is recorded alongside the
stricter actual measured maximum, 34.3 GPa.

Figure 4 identifies the experiment as **300 K**. Laser annealing at 1500–2000 K
relaxed differential stress before room-temperature measurement; it is not hot
P–V–T coverage. The reported discontinuity near 55 GPa and its magnetic
interpretation do not extend the low-pressure fit or constitute an implemented
magnetic phase diagram. XRD alone is not a direct magnetic-state measurement.

The primary pressure calibrant and medium were **NaCl**, citing Brown (1999),
*The NaCl pressure standard*, and Ono, Kikegawa and Ohishi (2006), *Structural
property of CsCl-type sodium chloride under pressure*. Ruby with the
Dorogokupets–Oganov (2007) calibration supplied preliminary pressures. The
per-row NaCl cells and detailed phase/scale assignment are not supplied, so
the original fit-domain pressures cannot be independently re-reduced. Figure 2
does give **B2 NaCl a=2.872 Å at 58.5 GPa**; that isolated, out-of-fit-domain
checkpoint is preserved in dataset provenance. It does not recover the missing
low-pressure calibrant readings. No later NaCl calibration is substituted.

Source numerical-data licensing is unspecified. CC0 applies only to
contributor-held transcription, normalization, metadata and arrangement rights;
it does not relicense the article, text, figures or PDF.

## Numerical checks and separate diagnostic

Run the repository maintenance script:

```bash
uv run python -m scripts.ono_mibe_2010_fe3c
```

The [numerical report](../data/ono-mibe-2010-fe3c.json) keeps three comparisons
separate:

1. **Source-equation check.** Three manual readings of the dashed experimental
   Figure 4 curve, at 10, 20 and 30 GPa, agree with the installed expression
   within **0.41 Å³**. Pixel coordinates, axes calibration and source-image
   checksum are retained. The tolerance covers four original-image pixels,
   line thickness and coordinate-reading ambiguity; it is a graphical bound,
   not an experimental error or parameter uncertainty. These are calculated
   curve checkpoints, never additional observations.
2. **Observation agreement.** For the five nonambient fit-domain observations,
   the published curve has pressure RMS **0.41015 GPa**, maximum **0.67863 GPa**,
   and maximum volume residual **0.37737 Å³**. Four volume residuals are within
   their printed widths; the 14.6 GPa row is about 1.26 widths away. With missing
   pressure errors and unspecified volume-error confidence/covariance, this is
   a descriptive comparison, not a rejection or a formal goodness-of-fit test.
3. **Conditional independent refit.** Equal pressure weights on those five
   nonambient rows, with measured V0=155.4 Å³ fixed, give **K0=165.92506 GPa,
   K′0=6.70657**, pressure RMS **0.37986 GPa**. Differences from the source are
   −1.07494 GPa and +0.00657. The original objective, weights and V0 constraint
   are unknown, and no source coefficient errors exist for an uncertainty
   comparison. This is a separately labelled diagnostic, not parameter parity
   or recovery of the authors' workflow. It is not a new installed EOS record.

The common ledger therefore reports **source reconstruction / source-equation
scope**, rather than asserting coefficient parity. Source verification,
graphical equation agreement, observation residuals, conditional refitting and
recommendation remain distinct.

## Integration and verification

The original experimental curve is an ordinary named, nondefault source option,
with qualified use inside 300 K and 0–34.3 GPa. Existing Scott/Li/Prescher
records, datasets, structures and defaults are preserved. Sata's original
published EOS and its separate recalculated-pressure diagnostic are unchanged
apart from the user-requested simple scatter caveat.

Focused checks cover literal transcription, missing errors, source/DFT/
recalibration separation, source checkpoints, diagnostic selection, public
evaluation/inversion, validity, Python/Rust loading and metadata round trips.
Catalog/source/refit ledgers and generated summaries are updated together.

Verification results: **546 Python checks** across the focused scientific,
ledger, EOSMAT and material suites, and **two Rust source-record tests**, pass.
Scoped Ruff/formatting and the maintained mypy check pass. Coverage and paper
register generation are current. The new entry matches a fresh full numerical
ledger calculation; pre-existing global ledger ordering/metadata drift is
preserved instead of rewriting unrelated findings. Strict MkDocs completes
notebook processing but aborts on the existing Litasov link to a repository
script; it reports no new Ono–Mibe or Sata warning.
Native reconstruction assessments match the ledger semantically. The byte-level
standardization check passes for `cementite`; for `cementite_nonmagnetic` it
still flags an assessment-key ordering difference already present in the
pre-change snapshot. No scientific value or assessment differs there.

Snapshot comparisons preserve all 683 earlier ledger entries and their order,
every earlier cementite record/dataset, both material defaults and all other
material files. The Sata material differs only in its short warning; its
coefficients, observations, limits, conditional parity and all recalculated-
pressure/confidence diagnostics are unchanged. No commit, push, release or
deployment was made.
