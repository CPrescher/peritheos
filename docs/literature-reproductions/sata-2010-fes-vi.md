# Sata and Ohfuji: cold FeS-VI EOS audit

Audit date: 2026-10-03. Scope: the cold FeS-VI curve re-reported in
[Sata et al. (2010), JGR 115, B09204](https://doi.org/10.1029/2009JB006975),
with earlier-study provenance to Ohfuji (2007) and Sata (2008).

**The recovered cold data can be refitted directly in EosFit.** The equal-weight
refit approximately reproduces Sata's coefficients and their quoted errors
when expressed at the same fixed reference volume. The exact original
weighting and calibration protocol remain unresolved. This audit does not
promote a replacement catalog record or validate the Morard thermal EOS.

## Observations and provenance

The publisher's rendered Table 1 contains 13 FeS-VI and six FeS-VII rows.
The factual transcription
[`fes-sata-2010-table1.csv`](../data/fes-sata-2010-table1.csv)
retains both phases, source order, run identifiers, Z, and original numerical
tokens with parenthesized errors. All 13 VI rows are fitted; VII is excluded
by its phase label. VI measurements at 185.5, 198.2 and 214.5 GPa are retained
even though the phases coexist above the transition pressure.

The VI measurements are at 300 K after annealing, covering 39–214.5 GPa and
61.31–82.92 Å³/cell. The cell contains four FeS formula units, hence eight
atoms. Parenthesized errors refer to the last digits; their confidence level
and cross-row covariance are unavailable. Absent VII volume errors remain
unknown rather than zero. Molar volumes are derived using Z and Avogadro's
constant, without altering the observed cell volumes.

Table 1 footnote e attributes the FeS data to
[Ohfuji et al. (2007)](https://doi.org/10.1007/s00269-007-0151-0) and
[Sata et al. (2008)](https://doi.org/10.2138/am.2008.2762).
Figure 5 identifies VI with Ohfuji and VII with Sata. The complete original
2007/2008 articles, supplements and author fit files were not recovered.
This is a recovery of their numerical re-report, not proof of complete
original-input recovery.

The FeS subtable identifies MgO/Ar calibration collectively but supplies
neither a per-row marker assignment nor marker volumes. Independent pressure
recalibration is therefore unavailable. The
[`source manifest`](../data/fes-sata-2010-source-manifest.json)
records the cited scale coefficients, source hash, retrieval and uncertainty
limitations. No full copyrighted article is bundled.

## Reference convention and equation

The visually inspected Equation 1 is modified BM3 at a finite pressure.
Writing \(f=[(V_r/V)^{2/3}-1]/2\),

\[
 P=(1+2f)^{5/2}\left[P_r+(3K_r-5P_r)f+
 \left(\frac92K_r(K'_r-4)+\frac{35}{2}P_r\right)f^2\right].
\]

It satisfies \(P(V_r)=P_r\), \(K(V_r)=K_r\), and
\(dK/dP|_{V_r}=K'_r\). Simply adding \(P_r\) to ordinary BM3 would give a
different equation. At \(P_r=0\) it reduces to ordinary BM3.

Sata selects \(V_r=12.37\) Å³/atom = 98.96 Å³/cell =
14.89877624024 cm³/mol FeS, and reports

| Parameter at fixed Vr | Published | Independent equal-weight fit | Direct EosFit, re-expressed at Vr |
|---|---:|---:|---:|
| Pr (GPa) | 0.0 ± 4.2 | −0.04423 ± 4.21354 | −0.04402 ± 4.21361 |
| Kr (GPa) | 148 ± 16 | 148.36883 ± 15.69412 | 148.36813 ± 15.69216 |
| K′r | 4.53 ± 0.34 | 4.525118 ± 0.334080 | 4.525131 ± 0.334083 |

The selected reference volume has no quoted uncertainty. Its printed precision
must not be turned into a measured V0 error. The close reproduction of all
three central coefficients and errors is evidence consistent with equal
pressure weights; it does not establish the authors' actual objective.

The earlier Ohfuji parameterization quoted in the 2010 HTML has
Pr=36 GPa, Kr=306 GPa and K′r=3.81. **306 GPa is a modulus at 36 GPa, not K0.**
The HTML prints its reference volume as “12.615 cm³/atom”, a dimensionally
implausible unit. A diagnostic interpretation as cm³/mol FeS gives
83.79080133 Å³/cell. At that volume the Sata curve evaluates to
P=35.96989 GPa, K=298.51207 GPa and K′=3.964278, within the earlier quoted
parameter errors. Re-expressing the same Sata curve at this reference preserves
pressures to 1.6×10⁻¹³ GPa. The unit interpretation is conditional until the
original Ohfuji source is inspected.

## Independent fit and residuals

The audit uses a linear SVD benchmark in (Pr, Kr, Kr×K′r) at fixed Vr;
the physical-parameter covariance is transformed from that solution.
Its equal-weight pressure RMS is 1.980935 GPa. Replay of the unchanged
published coefficients has RMS 1.986037 GPa and residuals
Pcal−Pobs from −3.774724 to +4.234706 GPa.

Pressure-error weighting and effective variance including volume errors are
separate diagnostics, with their objectives, covariance, rank, degrees of
freedom and convergence in
[`sata-2010-fes-vi-audit.json`](../data/sata-2010-fes-vi-audit.json).
Their RMS values remain approximately 2 GPa. The scatter substantially exceeds
many quoted point errors. Morard's ±3 GPa claim is not used as an acceptance
criterion for this different source.

## Actual EosFit console refits

The unmodified official Mac **EosFit7c 7.60**, version date 13-May-2021,
converged with all 13 VI points and ordinary BM3, refining V0, K0 and K′0.
This is a zero-pressure representation of the same mBM3 curve family.
The source-reference values above are a subsequent analytic curve
re-expression, with local covariance propagation, not another EosFit fit.

| EosFit objective | V0 (Å³/cell) | K0 (GPa) | K′0 | Pressure RMS | Reported weighted χ²/dof |
|---|---:|---:|---:|---:|---:|
| Equal pressure weights | 98.93066 ± 2.80286 | 148.56732 ± 34.40257 | 4.52374 ± 0.45851 | 1.980916 GPa | — |
| Supplied P and V errors | 99.76350 ± 1.31052 | 138.01448 ± 16.68440 | 4.67266 ± 0.28835 | 2.003175 GPa | 12.2004 |

The free-zero-pressure modulus error differs from the fixed-Vr source modulus
error because these are different reference coordinates. Comparing them
without transforming their correlated covariance would be misleading.
All errors in this table are EosFit's reported, conditionally scaled esds;
the source's error confidence remains unspecified.

The actual inputs, executable hash, macros, stdout, logs, saved EOS files
and full saved covariance are retained in the
[direct-run manifest](../data/fes-eosfit7c/manifest.json).
See [direct EosFit FeS refits](fes-eosfit-direct.md) for the thermal and joint
checks and rerun instructions.

## Validation and disposition

`python -m scripts.audit_sata_2010_fes_vi --check` verifies the independent
report. Tests check raw error tokens and phase selection, the finite-pressure
reference identities and numerical derivatives, reduction to native BM3,
reference transformation invariance, conditional covariance, and actual
EosFit input/output hashes and convergence.

Disposition: **approximate coefficient reproduction, original protocol
unresolved**. The recovered cold table is directly refittable in EosFit.
Pressure-scale verification and independent physical validation remain open.
No author outreach was performed.
