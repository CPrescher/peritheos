# Peritheos runtime package: code, structures and EOS parameters

Peritheos code is MIT. Crystal structures and compiled scientific parameters
retain their own provenance; the code license is not a blanket grant of
third-party rights. Numerical facts, independently implemented equations,
publication text and source compilations are distinct questions.

This package contains no experimental observation tables, original author
artifacts, workbook-cell extractions or detailed numerical research archives.
EOSMAT records preserve coefficients, units, uncertainties, reference states,
validity limits, source citations, calibration conventions and scientific
qualifications. The source repository retains the complete scientific evidence;
runtime exports do not change published values or recommendation status.

The ten archived COD CIFs retain their public-domain/CC0 source declarations.
See the [COD notice](datasets/cod-cif-sources.LICENSE.md), `datasets/cif-license-review.json`
and the full CC0 text in the distribution's license metadata. Other inline
structures and parameter compilations retain their individual provenance;
this packaging split alone does not establish Debian clearance for every entry.

The installed `eos-structure-rights-review.json` separately assesses all 685 EOS
records and 234 structural/diffraction components. See
[the compilation notice](datasets/eos-structure-compilations.LICENSE.md) for the numerical-content basis
and the seven unresolved Holland/THERMOCALC parameter extracts, which remain
present. The ten inherited peak arrays are replaced by independently calculated
X-ray patterns; see [their scoped notice](datasets/calculated-patterns.LICENSE.md)
and `calculated-pattern-inputs.json`. Six inputs use COD CC0 CIFs, one is an
independently authored ideal model, and three use attributed numerical structural
parameters without a separate source-publication grant. Generated-output CC0
covers contributor-held rights only. No blanket clearance is asserted.
The two formerly ICSD-linked post-perovskite topologies now use freshly acquired
CC0 COD 9009217 data. Their numerical values are unchanged, with ab initio
coordinate/observed-cell and pure-MgSiO3 proxy qualifications retained.
The ten named COD resource scopes remain separate from EOS/observation rights.
The pinned Dioptas migration's complete upstream MIT attribution is retained
in [the Dioptas notice](datasets/dioptas-materials.LICENSE.md); it does not clear third-party rights.

Compiled Rust third-party licenses and source attributions are in
`third_party/rust/NOTICE.md` and the complete adjacent notice bundle. A rebuild
with another compiler requires reviewing its matching standard-library notices.

Experimental data are optional: `python -m pip install "peritheos[datasets]"` adds
the matching `peritheos-datasets` distribution. That package admits only an
explicit source-scoped selection and includes its own dataset notices and
licenses. Its installed review distinguishes recorded source grants from
maintainer factual-content assessments. Unspecified source terms remain
unspecified, and admission is not Debian approval. The 21 admitted NC/ND
numerical extractions retain source-specific qualifications; five original
Sakai/Campbell dataset entries remain held as listed in that review.
No installation downloads observations at runtime. Other evidence and terms
are linked from the material's `data_availability` metadata.

Detailed repository license inventory:
https://github.com/CPrescher/peritheos/blob/main/docs/research/license-inventory.md
