# Pressure-scale normalization

Pressure-scale normalization changes the pressure coordinate of a published
sample EOS without pretending that a different calibrant was measured in the
original experiment. Treat it as a graph:

```text
sample EOS -> published source scale -> documented cross-calibration edge -> target scale
```

The target is an **internally consistent pressure-scale family**, not
necessarily one physical calibrant. An Au-calibrated experiment can finish on
the Au member of that family, a Pt-calibrated experiment on its Pt member, and
a ruby experiment through a ruby-linked Au or Pt member. The pressures are
comparable because all paths terminate on the same target family.

## Valid edge types

### The same XRD standard on two scales

If two EOS records describe the same physical standard, such as Au on scale A
and Au on scale B, transform pressure through a virtual standard volume:

\[
V_{\mathrm{Au}}^*(P,T)=
\operatorname{EOS}_{\mathrm{Au,A}}^{-1}(P,T),
\qquad
P_{\mathrm{B}}=
\operatorname{EOS}_{\mathrm{Au,B}}(V_{\mathrm{Au}}^*,T).
\]

`V_Au*` is neither the sample volume nor a measured Au volume. It is only the
coordinate at which the source Au EOS has the supplied pressure. This route is
useful even when only a published sample EOS remains available.

```python
from peritheos import recalculate_xrd_pressure_scale

result = recalculate_xrd_pressure_scale(
    source_pressure_gpa=[10.0, 50.0, 80.0],
    source_standard_eos_record="gold_dewaele_2004_vinet_5",
    target_standard_eos_record="gold_fei_2007_vinet_2",
    temperature_k=300.0,
    check_validity=True,
)

print(result.target_pressure_gpa)
print(result.implied_standard_volume)
```

Peritheos requires both records to belong to the same material document and
use the same conventional unit cell. It deliberately rejects an Au-to-Pt
call.

### Direct volume-volume cross-calibration

A co-compression experiment can measure two standards at the same physical
state and publish paired observations `(V_Au, V_Pt, T)`. Those data define a
genuine cross-material edge:

\[
V_{\mathrm{Au}}
\xrightarrow{\text{paired cross-calibration}}
V_{\mathrm{Pt}}.
\]

The paper's pressure anchor and fit assumptions remain part of the edge. The
numerical Au volume must **never** be inserted directly into a Pt EOS. Au and
Pt have different reference volumes and compression, and may use different
cell conventions.

Peritheos now exposes an explicit edge registry through
`list_cross_calibration_edges()` and `get_cross_calibration_edge()`. An edge
states whether it is executable and records its primary source, table or
method location, direction, and numerical transformation. The registry does
not pretend that every simultaneous experiment supplies a volume-volume
function: Chidester's bundled table, for example, contains KCl volume and
Pt-derived pressure but not raw Pt cell volume.

### An internally consistent multi-standard family

Some studies jointly optimize several calibrants on one anchor. Their Au, Pt,
MgO, and other EOSs are separate material equations but members of one scale:

```text
Au-A pressure -> virtual Au volume -> target-family Au EOS -> target pressure
Pt-B pressure -> virtual Pt volume -> target-family Pt EOS -> target pressure
```

This is not an Au-to-Pt volume transformation. Each source stays on its own
physical marker; consistency comes from the shared target-family construction.

The complete Sokolova et al. (2013) family is present in the shared library,
including Au `gold_sokolova_2013_holzapfel_4` and Pt
`platinum_sokolova_2013_holzapfel_3`. The compact convenience catalog also
contains all four Fei et al. (2007) members, while the shared `.eosmat` library
currently contains its Au and Ne records. Check availability before choosing a
target family.

## Ruby paths

Ruby calibrations use the corrected R1 wavelength ratio as their common
coordinate. Ruby-to-ruby conversion is

\[
R_1^*=f_C^{-1}(P_C), \qquad P_D=f_D(R_1^*).
\]

```python
from peritheos import recalculate_ruby_pressure

p_target = recalculate_ruby_pressure(
    [10.0, 50.0, 100.0],
    "ruby_mao_1986",
    "ruby_dewaele_2004",
)
```

To reach an XRD family, the path also needs an XRD-standard EOS calibrated
against a known ruby scale:

```text
source ruby C pressure
    -> R1 ratio
    -> bridge ruby scale
    -> virtual volume of the bridge XRD standard
    -> same-standard member of the target family
    -> target pressure
```

```python
from peritheos import recalculate_ruby_to_xrd_pressure

result = recalculate_ruby_to_xrd_pressure(
    source_pressure_gpa=[10.0, 50.0, 80.0],
    source_calibration_record="ruby_mao_1986",
    bridge_eos_record="gold_dewaele_2004_vinet_5",
    target_eos_record="gold_fei_2007_vinet_2",
    temperature_k=300.0,
    check_validity=True,
)
```

This converts Mao (1986) pressure to the Dewaele (2004) ruby scale through the
R1 ratio, inverts the ruby-linked Dewaele Au EOS, and evaluates the Fei Au EOS.
No measured Au volume is required.

Use `list_ruby_xrd_bridges(target_eos_record)` to find selectable direct
bridges. A
Fei Au target can currently be reached through the Dewaele et al. (2004) or
Takemura and Dewaele (2008) Au records. A paper merely comparing ruby with an
XRD standard is insufficient unless its ruby equation and linked XRD EOS are
numerically recoverable.

## Normalizing complete sample EOS records

`recalculate_eos_pressure_scale()` starts at sample volume. It evaluates the
published sample EOS, reads its `pressure_calibration`, and recursively finds
an executable path to the requested target. The returned
`calibration_path`, `edge_identifiers`, and `intermediate_states` make that
route auditable.

An Au-calibrated sample needs no explicit bridge:

```python
from peritheos import recalculate_eos_pressure_scale

normalized = recalculate_eos_pressure_scale(
    source_eos_record="forsterite_finkelstein_2014_bm3_1",
    sample_volume=[290.0, 280.0, 270.0],
    target_standard_eos_record="gold_sokolova_2013_holzapfel_4",
    sample_temperature_k=298.0,
    standard_temperature_k=298.0,
    check_validity=True,
)
```

A ruby-calibrated sample can now use automatic path discovery:

```python
normalized = recalculate_eos_pressure_scale(
    source_eos_record="tungsten_dewaele_2004_vinet_2",
    sample_volume=[31.0, 29.0, 27.0],
    target_standard_eos_record="gold_sokolova_2013_holzapfel_4",
    sample_temperature_k=300.0,
    standard_temperature_k=300.0,
    check_validity=True,
)
```

Pass `bridge_eos_record="gold_dewaele_2004_vinet_5"` only when the route must
use that particular ruby--Au bridge. Without it, the graph chooses the
shortest executable route. Inspect the returned path rather than assuming
which equally short bridge was selected.

### The Au, Pt, and ruby example

Suppose MgO was reduced with Au-A, MgSiO3 with Pt-B, and CaSiO3 with ruby-C.
To compare them on one family F:

1. Recalculate MgO through `Au-A -> Au-F`.
2. Recalculate MgSiO3 through `Pt-B -> Pt-F`.
3. Choose a ruby-linked Au or Pt bridge and transform
   `ruby-C -> bridge ruby -> bridge Au/Pt -> Au-F/Pt-F` for CaSiO3.
4. Restrict all results to their common pressure, temperature, phase, and
   stress-condition envelope.

The endpoints must be members of the same family F. Choosing Au from Fei et
al. and Pt from Sokolova et al. produces two named scales, not one normalized
dataset.

## KCl as pressure marker and bridge

B2 KCl needs special treatment because it is commonly both the thermal
insulator and the pressure marker in laser-heated DAC experiments. The
preferred record currently bundled as `kcl_b2_dewaele_2012_vinet_3` is **not
Au-calibrated**. Dewaele et al. (2012) measured its 298 K compression in helium
using the Dorogokupets--Oganov (2007) ruby scale. Its high-temperature term is
an ab initio molecular-dynamics result rather than a simultaneous Au, Pt, or
ruby measurement.

Consequently, a pressure assigned with that KCl EOS inherits the
Dorogokupets--Oganov ruby pressure basis. At 300 K it can be normalized through
the existing ruby-linked Au path:

```text
KCl Dewaele pressure
    -> Dorogokupets--Oganov ruby pressure basis
    -> Takemura--Dewaele Au bridge
    -> target-family Au EOS
```

This ancestry path does not insert a KCl volume into an Au EOS. It transforms
the pressure value through the ruby calibration on which the KCl EOS was
based. For a hot KCl measurement, the scale transformation uses the bridge
standards at their calibration reference temperature, normally 300 K; the
measured KCl temperature is used only to evaluate the KCl thermal EOS. The
model-derived KCl thermal pressure and its limitations remain part of the
uncertainty.

Two experimental records now provide stronger KCl-to-metal lineage:

- [Tateno et al. (2019)](https://doi.org/10.2138/am-2019-6779) measured B2-KCl
  and Pt simultaneously to 229 GPa at 300 K and to about 60 GPa and 2560 K.
  Pressure was assigned with the Sokolova Pt EOS. The executable record
  `kcl_b2_tateno_2019_vinet_4` is therefore a direct descendant of
  `platinum_sokolova_2013_holzapfel_3`.
- [Chidester et al. (2021)](https://doi.org/10.1103/PhysRevB.104.094107)
  measured simultaneous Pt--KCl P-V-T data to 167 GPa and 2400 K and publishes
  the observations in its supplement. The record
  `kcl_b2_chidester_2021_bm3_5` and all 155 author-deposited rows are bundled.
  It used the Dorogokupets--Oganov Pt EOS and explicitly models the lower
  effective temperature of the KCl insulation relative to the Pt surface.

The complete Dorogokupets--Oganov Pt Vinet plus four-oscillator Helmholtz EOS
is bundled as `platinum_dorogokupets_oganov_2007_vinet_4`. The direct
Chidester-to-Pt edge is therefore executable. Temperature must not simply be
copied across this edge: Chidester defines the reported KCl temperature by

\[
T_{\mathrm{KCl}}=\frac{3T_{\mathrm{surface}}+295\ \mathrm{K}}{4},
\qquad
T_{\mathrm{surface}}=\frac{4T_{\mathrm{KCl}}-295\ \mathrm{K}}{3}.
\]

KCl is the insulating layer between the hot Pt surface and the approximately
295 K diamond. The paper additionally reduces the measured surface
temperature by 3% for the axial gradient through the metal foil, so the
nominal EOS temperature is
$T_{\mathrm{Pt,avg}}=0.97T_{\mathrm{surface}}$. The graph combines both steps
before evaluating another Pt or metal EOS. Thus a deposited KCl temperature
of 1500 K implies 1901.67 K at the surface and 1844.62 K for the average Pt
foil. This is still a derived EOS transformation, not
recovery of an unpublished Pt volume: Chidester's table contains KCl volume
and Pt-derived pressure but no raw Pt cell volume, so exact observation-level
re-reduction and refitting remains unavailable.

This qualification does not make the Chidester KCl EOS fictitious. It is an
experimental pressure-marker calibration: KCl diffraction supplies volume and
simultaneous Pt diffraction supplies pressure. But its high-temperature
coordinate is modeled rather than a homogeneous KCl temperature measured
directly. Treat it as an effective-temperature-calibrated P-V-T EOS for a
comparable laser-heated DAC geometry, not as geometry-independent calorimetric
evidence for intrinsic KCl thermal properties. The room-temperature compression
data do not carry this thermal-gradient qualification.

The `.eosmat` record preserves this setup in the structured
`experimental_configuration` field: KCl was itself the pressure-transmitting
medium and thermal insulator, arranged as two dried layers around the Pt foil;
it was not a separate sample embedded in some other pressure medium.

```python
from peritheos import recalculate_pressure_calibration_path

normalized = recalculate_pressure_calibration_path(
    100.0,
    "kcl_b2_chidester_2021_bm3_5",
    "gold_fei_2007_vinet_2",
    temperature_k=1500.0,
)
print(normalized.calibration_path)
print(normalized.intermediate_states[0]["target_temperature_k"])
# 1844.616... K average Pt-foil temperature
```

## Diamond-edge paths

Two room-temperature diamond-anvil high-frequency-edge scales are executable:

- `diamond_raman_akahama_2006`, calibrated to 310 GPa against Holmes Pt;
- `diamond_raman_eremets_2023`, calibrated to about 500 GPa against the
  Fratanduono et al. ramp-derived Au isotherm, which is bundled as
  `gold_fratanduono_2021_vinet_7`.

Use `recalculate_diamond_raman_pressure()` to preserve the measured
`omega/omega0` ratio between diamond scales, or use
`recalculate_pressure_calibration_path()` to continue from either diamond
scale to an XRD target. Diamond Raman means the high-frequency edge of the
stressed anvil, not the Raman mode of a separate hydrostatic diamond sample.

## Literature that can support an edge

The useful literature is small when “cross-calibration” is interpreted
strictly. A qualifying source must publish simultaneous paired observations, a
numerical volume-volume relation, or an explicitly joint model family with a
traceable anchor. Merely plotting independent pressure scales together does
not qualify.

| Study | Standards and domain | Role | Library status |
|---|---|---|---|
| [Fei et al. (2004)](https://doi.org/10.1016/j.pepi.2003.09.018) | Au, Pt, MgO and other metals; 300–2173 K, to about 28 GPa | Simultaneous high-temperature cross-calibration | Not yet encoded as a cross-material edge |
| [Dewaele et al. (2004)](https://doi.org/10.1103/PhysRevB.70.094112) | Ruby with Cu, W, Al, Au, Pt, and Ta; 300 K to 153 GPa | Ruby-to-XRD bridge family | Ruby equation and several member EOSs executable |
| [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) | Au, Pt, NaCl-B2, and Ne; thermal EOSs | Internally consistent model family | Complete in convenience catalog; partial in shared `.eosmat` records |
| [Dorogokupets and Dewaele (2007)](https://doi.org/10.1080/08957950701659700) | MgO, Au, Pt, and NaCl | Joint model family | Not yet complete in shared library |
| [Dorogokupets and Oganov (2007)](https://doi.org/10.1103/PhysRevB.75.024115) | Ruby, metals, MgO, and diamond | Semiempirical family and ruby scale | Ruby equation and complete Pt thermal EOS executable |
| [Takemura and Dewaele (2008)](https://doi.org/10.1103/PhysRevB.78.104119) | Ruby and Au; 300 K to 123 GPa | Ruby-to-Au bridge | Executable bridge EOS |
| [Hirose et al. (2008)](https://doi.org/10.1016/j.pepi.2008.03.002) | Au and MgO; high temperature to 140 GPa | Thermal cross-calibration | Not yet encoded as an edge |
| [Tateno et al. (2019)](https://doi.org/10.2138/am-2019-6779) | B2 KCl and Pt; 300 K to 229 GPa and high temperature to about 60 GPa | KCl bridge tied to Sokolova Pt | EOS and executable edge bundled |
| [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) | Au, Pt, Mo, MgO, NaCl-B2, and Ne; 300 K to about 265 GPa | Relative family anchored to Tange MgO | Five convenience records; paired edge not encoded |
| [Sokolova et al. (2013)](https://doi.org/10.1016/j.rgg.2013.01.005) | Diamond, MgO, Ag, Al, Au, Cu, Mo, Nb, Pt, Ta, and W; model to 4 Mbar and 3000 K | Internally consistent model family | All eleven members executable |
| [Ye et al. (2017)](https://doi.org/10.1002/2016JB013811) | Au, Pt, and MgO; to about 140 GPa and 2500 K | Simultaneous thermal cross-calibration | Not yet encoded as family or edge |
| [Shen et al. (2020)](https://doi.org/10.1080/08957959.2020.1791107) | IPPS Ruby2020; to 150 GPa | Modern ruby scale | Executable |
| [Chidester et al. (2021)](https://doi.org/10.1103/PhysRevB.104.094107) | B2 KCl and Pt; to 167 GPa and 2400 K | Simultaneous thermal KCl--Pt calibration on Dorogokupets--Oganov Pt | EOS, full KCl P-V-T table, temperature mapping, and direct Pt edge executable |
| [Eremets et al. (2023)](https://doi.org/10.1038/s41467-023-36429-9) | Diamond Raman edge and Au; to about 500 GPa | Simultaneous optical--XRD calibration on ramp-derived Au | Scale, Au anchor, and edge executable |
| [Sakai et al. (2025)](https://doi.org/10.1038/s43246-025-00792-5) | Cu, Re, Pt, W, Au, Mo, Fe, MgO, and NaCl; to about 431 GPa | Direct multi-standard volume relations, Cu anchored | Not yet bundled |
| [Shen and Smith (2026)](https://doi.org/10.1103/fxgq-96sg) | Pt, Cu, Au, Mo, Ta, W, Fe, MgO, and NaCl; 300 K to 140 GPa | Simultaneous paired calibration, Cu anchored | Member fits bundled; paired data and anchor edge not executable |

A simultaneous volume-volume dataset is the strongest cross-material edge
because the relationship is observed directly. A joint model family is usable
as a common target, but its members share an anchor and assumptions; it is not
independent experimental confirmation.

## Choosing the target family

- Sokolova et al. (2013) provides the broadest target set currently executable
  from shared records. It is a joint model family, not paired raw data.
- Fei et al. (2007) is appropriate within its thermal domain when every needed
  family member is available. Do not substitute a member from another family.
- Dorfman et al. (2012), Shen and Smith (2026), and Sakai et al. (2025) form
  attractive room-temperature experimental networks. Exact cross-material use
  requires their paired relations and anchor, not only fitted coefficients.
- At high temperature, prefer a thermal cross-calibration such as Fei et al.
  (2004), Hirose et al. (2008), or Ye et al. (2017) over extrapolating a 300 K
  edge.
- To reproduce historical results, retain the authors' scale and store the
  normalization separately rather than replacing published parameters.

## Temperature, validity, and uncertainty

Every thermal calculation in an edge must represent the same physical
temperature. A 300 K isotherm cannot supply a high-temperature
cross-calibration merely because the sample EOS is thermal. When two published
room-temperature isotherms use fixed reference temperatures such as 298 and
300 K, Peritheos evaluates each at its own stated reference temperature and
treats the link as a nominal room-temperature comparison; it does not imply
that the 2 K difference was corrected. The usable domain is the intersection
of the sample EOS, source calibration, edge, target-member EOS, phase fields,
and relevant stress conditions.

Use `check_validity=True` while developing a conversion. Record pressure
medium, hydrostatic or deviatoric-stress qualifications, and calibrant phase
transitions. Marginal P and T ranges do not imply that every combination in
their rectangle was measured.

Where available, uncertainty should include sample volume and temperature,
source-scale covariance, paired-edge scatter, common-anchor uncertainty,
target-family covariance, and systematic stress and temperature corrections.
Member EOSs from a joint fit are correlated through the anchor; treating their
printed errors as independent usually underestimates uncertainty. Current
transformations preserve provenance and validity, but do not propagate a
covariance matrix through a multi-edge graph.

## Derived normalization versus exact re-reduction

There are two different scientific products:

1. **Derived EOS transformation.** Evaluate a published sample EOS on a
   volume grid, transform its pressures, and optionally refit the curve. This
   is what `recalculate_eos_pressure_scale()` enables without raw data.
2. **Observation-level re-reduction.** Recompute each experimental row from
   its measured sample volume, calibrant volume or ruby wavelength,
   temperature, corrections, and covariance, then refit. This needs the
   original row-wise data.

A refit of a transformed curve must use `record_kind: refit`, link to the
published record, identify the target scale and complete path, and document
its volume grid and fit objective. It must not overwrite the published EOS.

## Provenance required for a result

Store enough information to replay the path:

```json
{
  "source_sample_eos_record": "forsterite_finkelstein_2014_bm3_1",
  "source_pressure_calibration": "gold_fei_2007_vinet_2",
  "edge_kind": "same_standard_eos",
  "target_scale_family": "sokolova_2013",
  "target_standard_eos_record": "gold_sokolova_2013_holzapfel_4",
  "standard_temperature_k": 300.0,
  "validity_policy": "intersection; no extrapolation",
  "software": "peritheos",
  "result_kind": "derived_eos_transformation"
}
```

For a direct cross-material edge, also store the paired dataset or fit ID,
primary citation and table location, anchor EOS, phases and cell conventions,
P-T range, residual scatter, and covariance. The bundled edge registry records
these facts explicitly and never infers family membership from author names.

## Intentionally rejected operations

- Passing an Au volume to a Pt EOS or assuming equal compression ratios.
- Connecting different materials without paired data or a joint target family.
- Combining target members from unrelated families and calling it one scale.
- Using a ruby comparison without a recoverable equation and linked XRD EOS.
- Extending a room-temperature edge to laser-heating conditions without a
  documented thermal model.
- Calling a curve transformed from published parameters an exact re-reduction.

These restrictions prevent numerical transformations from asserting more
experimental information than the source papers contain.
