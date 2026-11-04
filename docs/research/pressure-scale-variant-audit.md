# Published pressure-scale variants

Normal discovery presents one source EOS per study, equation family and
experimental/constraint selection. Published fits using alternative pressure
scales remain executable references, identified by `pressure_scale_variant_of`.
No numerical coefficients, uncertainty values, observations, validation statuses
or existing explicit lookups are changed by this grouping.

The [machine-readable inventory](../data/pressure-scale-variants.json) contains
51 reference variants in 27 source/model groups across 19 materials. Regenerate
with `python scripts/audit_pressure_scale_variants.py`; `--check` detects drift.

| Study | Materials | Source selection | Reference variants |
|---|---|---|---:|
| [Sakai (2011)](../literature-reproductions/sakai-2011-nacl-b2.md), Table II, plus the explicit Sakai (2014) recalibration | NaCl-B2 | Matsui-Pt BM3 and Vinet; includes the later Yokoo-Pt BM3 re-reduction of the same 2011 observations | 7 |
| [Dewaele (2004)](../literature-reproductions/dewaele-2004-six-metals.md), Table II | Al, Cu, W | Revised-ruby Vinet; historical Mao fits retained as references | 3 |
| [Dewaele (2019)](../literature-reproductions/dewaele-2019-static-dac-metals.md), Table 1 | Au, Pt, Cu, Ta, Al, W, Co, Mo, Zn, Re, Fe, Be, Pb | First-listed Mao Vinet for a consistent display baseline; D&O fits retained as references | 13 |
| Dewaele (2008), Table III, and its Dewaele (2019) pressure recalibrations | Ni, Ag | Original 2008 Vinet sources; both 2019 pressure reductions share the same observations, equation and free-parameter choices | 4 |
| Sakai (2014), Tables 3–4, DOI 10.1016/j.pepi.2013.12.010 | hcp Fe and Fe0.9Ni0.1 | p1 BM3 and Vinet; p2–p6 use the alternative NaCl-B2 pressure reductions of the same data | 20 |
| [Tateno (2019)](../literature-reproductions/kcl-tateno-campbell-source-gaps.md), Table 1 | B2 KCl | Sokolova-Pt for each of the MGD and linear thermal forms | 2 |
| Ono (2006), DOI 10.2138/am.2006.2118, Tables 2–3 | MgSiO3 post-perovskite | Existing first-listed Anderson-Au BM2; Jamieson and Dewaele alternatives retained as references | 2 |

These source selections organize literature results; they do not rank the
accuracy of the pressure standards or change the material's recommended default.
BM3 versus Vinet, and MGD versus linear thermal pressure, remain different model
choices. The source temperatures and validity envelopes retain their original
conventions, including nominal 298 versus 300 K reference isotherms.

## Cases retained as distinct choices

The audit compared same-material, same-publication, same-model records and
cross-publication records explicitly sharing the same observation tables, and
inspected their calibration metadata, observations and source qualifications.
A common citation and a different calibrant do not by themselves establish
pressure-scale duplication.

- Jacobsen (2008) MgO: KCl-medium and helium-medium records use different
  measurements and stress conditions. Both remain source choices.
- Dewaele (2004) Au, Pt and Ta: the revised-ruby fits fix K0, while their Mao
  alternatives fit K0. These pairs change the fitting constraints as well as
  the pressure scale, so both remain source choices.
- Wang (2004) akimotoite: T0133 Au and T0150 NaCl are separate runs.
- Reynard (1996) akimotoite: the ruby and ice-VII assignments share the experiment,
  but their fits use different available row selections (16 versus 14), and the
  exact ruby equation and executable ice-VII calibration are unresolved.
- Sakai (2014) Fe and Fe–Ni p1 variants with externally fixed Mao/Dewaele V0:
  these change the reference-volume constraint rather than the pressure scale.
- Hirose (2008) Au: DAC-only versus combined observations, and the independent
  refit, remain separate. Miozzi (2020) Fe refits also differ in thermal/model
  choices, including atom normalization, rather than solely output scale.
- Calculated temperature isotherms, different computational approximations,
  separate phases/compositions, annealing subsets and equation families remain
  distinct even when their literature citations match.
- Anzellini (2014) and Dewaele (2019) Re remain separate analyses: the original
  fit mixes ruby, W and He gauges; the later fits use a qualified effective
  ruby-coordinate reanalysis. Even the same nominal D&O pressure label does not
  establish an identical gauge reduction or fitting workflow.

## Conversion and unresolved routes

`source.on_pressure_scale(target_node)` evaluates the source curve through the
audited calibration graph. It never reads a reference variant's sample-fit
coefficients. The converted function need not belong to the source's BM3/Vinet
family; its direct inverse uses the same calibration path in reverse.

Grouping is not a graph edge. In particular:

- Sakai (2011)'s precise Holmes implementation remains unresolved. The Holmes
  source reference cannot be automatically normalized through an invented link.
- Ono (2006)'s historical Jamieson-Au model is not bundled.
- Sakai (2014)'s mixed Brown-B1/alternative-B2 reductions lack a single resolved
  executable pressure coordinate in the existing metadata. The variants remain
  explicit references; grouping does not manufacture a conversion for them.
- Dewaele (2019) Re retains its incompletely specified W-linked calibration.
  The source audit's effective-ruby-coordinate recalculation is qualified;
  a blanket exact W-volume route is not inferred here.

Missing routes raise errors. Existing executable routes may still be limited
by their pressure, temperature and phase domains. Use `check_validity=True`;
the existing Sakai comparison documents where a full experimental-range curve
comparison extrapolates the linked standards. Conversion accuracy is separate
from experimental uncertainty and from agreement with separately rounded and
refitted published alternatives.
