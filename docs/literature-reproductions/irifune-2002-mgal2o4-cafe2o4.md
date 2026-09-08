# Irifune et al. (2002) CaFe2O4-type MgAl2O4

Audit date: 2026-09-08. Primary authority is the typeset article,
[DOI 10.1007/s00269-002-0275-1](https://doi.org/10.1007/s00269-002-0275-1),
checked through the publisher metadata and
[Toru Inoue's author-uploaded full text](https://www.researchgate.net/publication/225470056_In_situ_X-ray_observations_of_phase_transitions_in_MgAl2O4_spinel_to_40_GPa_using_multianvil_apparatus_with_sintered_diamond_anvils).
No numerical supplement was identified on those pages.

## Accepted representation

The p. 652 fit uses V0=240.6(3) A³ and K0=213(3) GPa with K0′=4.
Its BM3 reduces exactly to BM2. The existing four-formula-unit CF material
receives one published record; its independently sourced structure and default
remain intact. Figure 6 explicitly labels T=27°C, represented as 300.15 K.
The three open circles are this study; small filled dots are Yutani (1997)
and the large filled dot is Funamori (1998). Those comparison observations
are not additional inputs to the Irifune fit.

The bundled CSV preserves three compressed volumes and three paired recovered
volumes. The averaged V0 is not an extra independent observation. Only compressed
rows enter the diagnostic fit. Pressure uncertainties and paired calibrant
volumes are unavailable. The experimental section names Anderson (1989) Au and
fallback Jamieson (1982) MgO, without resolving each compressed row's assignment.

## Sparse-data decision and reproduction

Run `python scripts/reproduce_irifune_2002_mgal2o4.py`.
The independent BM2 calculation holds V0 and K0′ fixed. Its pressure-residual
least squares gives K0=213.865416 GPa, within the published error. The script
also gives K0=213.659926 GPa for normalized-volume residuals, matching the stated dependent
variable while leaving the unspecified source weighting unresolved.

Three compressed observations cannot justify an unconstrained three-coefficient
fit with residual degrees of freedom. No additional executable refit, free
pressure derivative, thermal model, covariance, or inferred confidence level is
added. The catalog refit is explicitly an unweighted pressure diagnostic;
coefficient agreement does not establish exact recovery of the author's fit.
The largest published-curve pressure residual is about 1.74 GPa, so these data
must not be described as precise pressure standards. The validity interval
records coverage, not phase stability throughout that interval.

## Candidate disposition

The discovery-only LitCurate ledger contains no row with this DOI or Irifune
as a publication author. This is a direct primary-source addition, not a new
invented LitCurate identifier. Accept exactly one CF EOS under this DOI.
The comparison values attributed to earlier papers do not create additional
Irifune records. The existing Funamori record remains separately attributed.
