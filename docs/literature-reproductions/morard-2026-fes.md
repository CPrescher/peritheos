# Morard et al. (2026): FeS phase diagram and thermal EOS

Audit date: 2026-10-02. Source: Morard, Antonangeli, Miozzi, Baron, Edmund,
Cerantola and Mezouar, *Physical Review Materials* **10**, 063605,
[DOI 10.1103/h4pj-rvxx](https://doi.org/10.1103/h4pj-rvxx).

The published BM3–MGD coefficients and all 146 supplied Table S1 observations
are preserved in `fes_vi_morard_2026_bm3_mgd`. The record is **deferred** from
the normal executable catalog. Standard MGD replay does not reproduce the
paper's claim that every calculated-minus-observed pressure lies within
±3 GPa. This finding concerns the available source model and data versions;
it does not establish an error in the experiments. No diagnostic fit replaces
the published coefficients.

## Sources and material identity

The user supplied the final nine-page article and `Table S1.xlsx`, with one
worksheet (`Sheet1`), headers in row 2 and 146 complete rows in A3:M148.
The factual CSV transcription retains every supplied digit, all thirteen
columns, original order and Excel row numbers. Source files and transcription
are checksummed in
`peritheos/data/datasets/morard_2026_sources/manifest.json`.
The supplied files' source-data license is unspecified.

The article's data-availability reference resolves to
[FeS XRD Data, Figshare version 1](https://doi.org/10.6084/m9.figshare.28143683.v1).
Its CC BY 4.0 numerical GSAS export and PT workbook are preserved separately,
with their hashes and download links. The workbook contains calibrant and
pressure-reduction information, but it is a separate source version. For
example, public `SI_FeS1` row 4, run 9, reports 29.214913231545133 GPa at
1450 K, whereas supplied Table S1 row 3 reports 29.59643816464321 GPa.
No public-deposit pressure is silently substituted for the final table.
Matching the GSAS export's full FeS lattice triple and temperature recovers
the KCl marker paired with every final Table S1 row. Identical duplicate
GSAS blocks have identical marker values; all their source-line locations are
retained in the audit report.
The deposit's diffraction patterns remain available at Figshare.

APS's supplemental landing page required subscription during the initial audit.
The initial supplied files lacked Table S2, Figure S5 and an EosFit file.
On 2026-10-03 the user additionally supplied the figures supplement,
`FeS6.eos` and two original author EosFit inputs. These inputs recover 21
literature cold observations and, in one version, 11 additional 300 K
observations. User-reported personal communication from Guillaume Morard confirms
that these 11 quenched observations were unsuitable and excluded from the final
fit. The database exposes only the selected 167-row input and literature cold
subset. Cold-then-thermal EosFit fits holding the cold coefficients fixed are
documented in the [selected-data staged audit](fes-eosfit-direct.md#selected-data-cold-then-thermal-fits-2026-10-03).
The historical two-input replay and refits are documented in the
[author-input follow-up](fes-eosfit-direct.md#author-input-follow-up-2026-10-03).
The earlier reconstructed-data diagnostics remain available separately.

The new `fes_vi` card represents nonmagnetic MnP-type FeS-VI. Its independent
experimental diffraction model is Ono et al. (2008),
[DOI 10.1016/j.epsl.2008.05.017](https://doi.org/10.1016/j.epsl.2008.05.017),
Table 1, page 483: Pnma, Z=4, measured at 56 GPa and 300 K, with
a=4.9371(9), b=3.0844(6), c=5.1576(11) Å and V=78.540(15) Å³.
The fully occupied Fe and S sites are both 4c:
(0.0051, 1/4, 0.2001) and (0.2085, 1/4, 0.5800), respectively.
These are the experimental Rietveld values, not the calculated columns in
that paper's Table 2. The cell contains four Fe and four S atoms.
Its measured cell is distinct from the extrapolated zero-pressure EOS volume.

Morard fits one common volume relation across FeS-VI, IV and IX because the
transitions show no volume discontinuity. This card does not supply IV or IX
atomic coordinates or calculate phase boundaries. Its fixed FeS-VI reference
structure must not be interpreted as the diffraction structure of all hot
polymorphs. No additional EOS parameterization from the structure paper was
added.

## Printed coefficients and declared audit interpretation

Section III.D, pages 063605-6 and 063605-7, reports:

| Coefficient | Published value | Treatment |
|---|---:|---|
| V0 | 15.40 ±0.45 cm³/mol FeS | Cold fit |
| K0 | 115.5 ±27.91 GPa | Cold fit |
| K0′ | 4.99 ±0.51 | Cold fit |
| θ0 | 417 K | Fixed, adopted hcp-Fe value |
| γ0 | 2.42 ±0.03 | Sole refined thermal coefficient |
| q | 1 | Assumed fixed |
| Tr | 300 K | Figure 6 reference isotherm |
| n | 2 atoms per FeS | Stoichiometric thermal normalization |

The public conventional-cell volume is Vcell=Vmolar ×4×10²⁴/NA, so
V0=102.28920653790897 Å³ and its reported error is
0.45 ×4×10²⁴/NA Å³. The exact SI Avogadro constant is used.
The thermal evaluator internally uses J/bar/mol volumes and n=2;
Z=4 is not an extra thermal-energy multiplier. Published coefficient-error
confidence and covariance are unspecified. No missing error is replaced by zero.

The paper names MGD and EOSFIT but does not print the volume-dependent Debye
law or provide the final configuration in the supplied files. The deferred
implementation declares the standard integrated interpretation:
γ=γ0(V/V0)^q, θ=θ0 exp[−γ0((V/V0)^q−1)/q], and
Pthermal=γ[E(V,T)−E(V,300 K)]/V. A variable-exponent
θ=θ0(V/V0)^−γ(V) sensitivity is also audited. Neither interpretation is
claimed as confirmed recovery of the authors' exact implementation.

## Source-table arithmetic and scope

Every supplied V equals a×b×c exactly. Every supplied dV equals
V×(err_a/a + err_b/b + err_c/c). This is linear error addition,
not independent-error quadrature; its statistical confidence is unspecified.
The CSV's uncertainty-column names do not establish one-sigma confidence.
All supplied errors are retained unchanged.

Every supplied C equals c/(√3 b), whereas the article's Figure 5 caption
and prose define C=b/(c/√3). The table therefore uses the reciprocal.
The largest absolute difference from the printed definition is 0.0591434.
This column is retained as a source order parameter and is not used to
assign phases automatically.

Actual Table S1 coverage is 29.596438–149.948070 GPa, 1220–2590 K and
66.899197–92.292831 Å³. This differs from the approximate prose ranges.
The paper describes the common model as applying at 40–150 GPa up to melting.
The 300 K reference comes from a separate literature fit. These are marginal
coverage statements, not a rectangular stability guarantee.

Pressure calibration is KCl using Dewaele et al. (2012),
DOI 10.1103/PhysRevB.85.214105. Section II averages marker temperature between
the cold diamond contact and hot FeS surface, following the cited earlier
protocols. The bundled `kcl_b2_dewaele_2012_vinet_3` record has the required
Vinet reference (54.5 Å³, 17.2 GPa, K0′=5.89) and linear thermal increment
0.00224(T−300 K) GPa. Independent replay of this expression with GSAS marker
volumes and T_KCl=0.75 T_FeS+0.25×295 K recovers all 146 final pressures with
maximum error 1.3×10⁻¹³ GPa. The temperature law is numerically reconstructed;
the article does not explicitly print these weights. Observed marker edges
and errors are bundled in a separate typed
`peritheos/data/datasets/fes-morard-2026-paired-kcl.csv`;
its calculated volumes, temperatures and replay pressures are declared as
derived audit quantities. Original cold-fit scale mapping remains unresolved.

## Independent replay and conditional fits

`scripts/audit_morard_2026_fes.py` evaluates BM3
and the Debye energy integral independently using 48-point Gauss–Legendre
quadrature, with physical n=2 and unchanged source coefficients.
The native Peritheos implementation is checked against that calculation away
from both the volume and temperature reference states. This is an equation
implementation check; it does not establish source-result reproduction.

| Debye-temperature law | Pressure RMS | Residual range | Rows exceeding ±3 GPa |
|---|---:|---:|---:|
| Integrated | 1.911067 GPa | −3.777566 to +6.149320 GPa | 15/146 |
| Variable exponent | 1.981292 GPa | −3.771041 to +6.315895 GPa | 16/146 |

Restricting the residual check to P≥40 GPa leaves 13 failing rows for the
integrated law and 14 for the variable-exponent law. Thus low-pressure
coverage alone does not explain the mismatch. For the integrated law,
Table S1 row 126 (116.740331 GPa, 2000 K) has the largest positive residual.

Gamma-only fits hold V0, K0, K0′, θ0, q and Tr fixed. Equal pressure weights,
supplied dP weights, and effective variance including supplied dP/dV/dT are
reported separately. Integrated-law fitted γ0 values are approximately
2.28198, 2.30425 and 2.31809, respectively. The equal-pressure fit has
conditional standard error 0.03231 and pressure RMS 1.80311 GPa.
These local Jacobian errors condition on the printed cold curve and do not
include cold-parameter uncertainty or cross-row systematics. Weight choices
are audit diagnostics because the source does not specify its full objective.
They do not recover the authors' exact protocol or establish coefficient parity.

The deterministic report
[`morard-2026-fes-audit.json`](../data/morard-2026-fes-audit.json)
retains every residual, row identifier, alternative fit, covariance for the
one free coefficient, convergence status, table checks and source extrema.
Run `python -m scripts.audit_morard_2026_fes --check` to verify it.

## Disposition

The source is retained in the database with `scientific_validation.status:
deferred`, `default: false`, and its failed replay documented. The normal
catalog excludes it. For an explicit audit calculation only, load the
document with `Material.from_eosmat(document, require_primary_validation=False)`;
this retains the deferred status. No validated default or replacement refit
is added.

The decisive remaining inputs are the final EosFit model/input/output files,
cold-fit observation selection and source pressure scales, exact thermal
convention, Table S2 and Figure S5. The final Table S1 pressure reduction has
been recovered independently; the exact source regression remains unavailable.
See the
[reviewed data request](../dataset-requests.md#morard-et-al-2026-fes).
No author outreach was performed.

## EosFit manual and direct-console follow-up (2026-10-03)

The official EosFit 7.60 manual confirms that full MGD uses the integrated
Grüneisen law already evaluated above; its q-compromise input prompt defaults
to N. The distinct q-compromise approximation has no refinable q and keeps
theta and gamma/V constant. This resolves the software's standard convention,
without claiming access to the authors' actual selected settings. Manual
source hashes are retained in
[eosfit-mgd-manual-audit.json](../data/eosfit-mgd-manual-audit.json).

The **actual official EosFit7c executable** has now been run on Table S1.
The unchanged published coefficients give RMS=1.910838 GPa and maximum
positive Pcal−Pobs residual=6.14860 GPa, with 15 rows outside ±3 GPa.
An equal-weight gamma-only refit converges to 2.28211 ± 0.03231. Thus the
residual discrepancy persists in the software named by the paper.
See [direct EosFit runs](fes-eosfit-direct.md) for all six fits, covariance,
warnings, macros and raw outputs, and the
[separate cold FeS-VI audit](sata-2010-fes-vi.md) for the 13 numerical
observations re-reported in Sata (2010). These newly recovered observations
improve cold-data availability; they do not recover the authors' final
Morard regression file or change the deferred catalog disposition.
