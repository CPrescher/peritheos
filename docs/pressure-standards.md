# Pressure standards

Any cataloged material EOS can turn a measured crystallographic unit-cell
volume and temperature into pressure in GPa. The EOS records also invert pressure to volume,
preserve scalar or NumPy array shapes, carry the primary-source provenance and
validity envelope with the calculation, and propagate measurement and published
parameter uncertainty.

```python
from peritheos import get_material

gold_material = get_material("au_fcc")
gold = gold_material.get_eos_record("au_fcc_dorfman_2012")
pressure = gold.pressure(volume=50.0, temperature=300.0)
volume = gold.volume(pressure, temperature=300.0)
```

The public volume unit is always
`angstrom^3/conventional_unit_cell`. The conventional cell and its contents are
part of the enclosing `Material`; this avoids the common factor-of-$Z$ ambiguity
between atomic, formula-unit, primitive-cell, conventional-cell, and molar
volumes. Temperature is in K. Pressure is in GPa.

## Curated pressure-scale convenience catalog

The table below documents the compact set exposed directly by
`get_eos_record()` and `list_eos_records()`. The cross-compatible `.eosmat`
library is broader: it contains 115 material documents and 146 independently
audited EOS records, accessed with `list_material_documents()` and
`Material.from_eosmat()`.

| EOS record identifier | Material and phase | Model | Published envelope | Primary source |
|---|---|---|---|---|
| `mgo_b1_tange_2009_vinet` | MgO B1, 4 formula units/cell | Fit3-Vinet + Mie-Gruneisen-Debye | 0–196 GPa, 300–3700 K, $0.652 \leq V/V_0 \leq 1.150$ | [Tange et al. (2009)](https://doi.org/10.1029/2008JB005813) |
| `mgo_b1_sokolova_2016` | MgO B1, 4 formula units/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `diamond_sokolova_2016` | diamond, 8 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `al_fcc_sokolova_2016` | Al fcc, 4 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `cu_fcc_sokolova_2016` | Cu fcc, 4 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `ag_fcc_sokolova_2016` | Ag fcc, 4 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `au_fcc_sokolova_2016` | Au fcc, 4 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `pt_fcc_sokolova_2016` | Pt fcc, 4 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `nb_bcc_sokolova_2016` | Nb bcc, 2 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `ta_bcc_sokolova_2016` | Ta bcc, 2 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `mo_bcc_sokolova_2016` | Mo bcc, 2 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `w_bcc_sokolova_2016` | W bcc, 2 atoms/cell | Holzapfel + Sokolova thermal | 0–400 GPa, 298.15–3000 K | [Sokolova et al. (2016)](https://doi.org/10.1016/j.cageo.2016.06.002) |
| `au_fcc_fei_2007` | Au fcc, 4 atoms/cell | Vinet + Fei Mie-Gruneisen-Debye | 0–125 GPa, 300–2330 K | [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) |
| `pt_fcc_fei_2007` | Pt fcc, 4 atoms/cell | Vinet + Fei Mie-Gruneisen-Debye | 0–94 GPa, 300–1873 K | [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) |
| `nacl_b2_fei_2007` | NaCl B2, 1 formula unit/cell | Vinet + Fei Mie-Gruneisen-Debye | 34–107 GPa, 300–1000 K | [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) |
| `ne_fcc_fei_2007` | Ne fcc, 4 atoms/cell | Vinet + Fei Mie-Gruneisen-Debye | 5–115 GPa, 300–1000 K | [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104) |
| `au_fcc_dorfman_2012` | Au fcc, 4 atoms/cell | 300 K Vinet | 1–250 GPa | [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) |
| `pt_fcc_dorfman_2012` | Pt fcc, 4 atoms/cell | 300 K Vinet | 5–228 GPa | [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) |
| `mo_bcc_dorfman_2012` | Mo bcc, 2 atoms/cell | 300 K Vinet | 43–213 GPa | [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) |
| `nacl_b2_dorfman_2012` | NaCl B2, 1 formula unit/cell | 300 K Vinet | 34–250 GPa | [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) |
| `ne_fcc_dorfman_2012` | Ne fcc, 4 atoms/cell | 300 K Vinet | 5–250 GPa | [Dorfman et al. (2012)](https://doi.org/10.1029/2012JB009292) |
| `lif_b1_dewaele_2019` | LiF B1, 4 formula units/cell | 300 K Vinet | 0–109 GPa | [Dewaele (2019)](https://doi.org/10.3390/min9110684) |
| `nacl_b1_dewaele_2019` | NaCl B1, 4 formula units/cell | 300 K Vinet | 0–35 GPa | [Dewaele (2019)](https://doi.org/10.3390/min9110684) |
| `nacl_b2_dewaele_2019` | NaCl B2, 1 formula unit/cell | 300 K Vinet | 37–155 GPa | [Dewaele (2019)](https://doi.org/10.3390/min9110684) |
| `kcl_b1_dewaele_2012` | KCl B1, 4 formula units/cell | 298 K Vinet | 0–2.6 GPa | [Dewaele et al. (2012)](https://doi.org/10.1103/PhysRevB.85.214105) |
| `kcl_b2_dewaele_2012` | KCl B2, 1 formula unit/cell | Vinet + linear thermal pressure | 2.6–200 GPa, 300–7000 K | [Dewaele et al. (2012)](https://doi.org/10.1103/PhysRevB.85.214105) |
| `kbr_b1_dewaele_2012` | KBr B1, 4 formula units/cell | 298 K Vinet | 0–2.3 GPa | [Dewaele et al. (2012)](https://doi.org/10.1103/PhysRevB.85.214105) |
| `kbr_b2_dewaele_2012` | KBr B2, 1 formula unit/cell | Vinet + linear thermal pressure | 2.3–200 GPa, 300–7000 K | [Dewaele et al. (2012)](https://doi.org/10.1103/PhysRevB.85.214105) |
| `cbn_zincblende_datchi_2007` | c-BN, 4 BN/cell | 295 K Vinet | 0–162 GPa | [Datchi et al. (2007)](https://doi.org/10.1103/PhysRevB.75.214104) |
| `diamond_datchi_dewaele_2008` | diamond, 8 atoms/cell | Vinet + Mie-Gruneisen-Debye | 0–80 GPa, 298–900 K | [Dewaele et al. (2008)](https://doi.org/10.1103/PhysRevB.77.094106) |
| `ni_fcc_dewaele_2008` | Ni fcc, 4 atoms/cell | 300 K Vinet | 0–157 GPa | [Dewaele et al. (2008)](https://doi.org/10.1103/PhysRevB.78.104102) |
| `ag_fcc_dewaele_2008` | Ag fcc, 4 atoms/cell | 300 K Vinet | 0–124 GPa | [Dewaele et al. (2008)](https://doi.org/10.1103/PhysRevB.78.104102) |
| `re_hcp_anzellini_2014` | Re hcp, 2 atoms/cell | 300 K Vinet | 0.64–144 GPa | [Anzellini et al. (2014)](https://doi.org/10.1063/1.4863300) |
| `aragonite_martinez_1996_bm2_2` | CaCO3 aragonite, 4 formula units/cell | staged BM2 with linear reference state | 0–8.18 GPa, 298–973 K | [Martinez et al. (1996)](https://doi.org/10.2138/am-1996-5-608) |

The bundled material library additionally includes the ten independently
validated, Cu-anchored 300 K Vinet fits from [Shen and Smith
(2026)](https://doi.org/10.1103/fxgq-96sg):

| `.eosmat` record identifier | Material and phase | Fitted pressure interval |
|---|---|---|
| `platinum_shen_2026_vinet_2` | Pt fcc | 0–140 GPa |
| `gold_shen_2026_vinet_3` | Au fcc | 0–140 GPa |
| `tantalum_shen_2026_vinet_2` | Ta bcc | 0–140 GPa |
| `tungsten_shen_2026_vinet_3` | W bcc | 0–140 GPa |
| `molybdenum_shen_2026_vinet_1` | Mo bcc | 0–140 GPa |
| `mgo_shen_2026_vinet_3` | MgO B1 | 0–140 GPa |
| `nacl_b1_shen_2026_vinet_1` | NaCl B1 | 0–31 GPa |
| `nacl_b2_shen_2026_vinet_2` | NaCl B2 | 33–140 GPa |
| `fe_shen_2026_vinet_1` | Fe bcc | 0–15.5 GPa |
| `iron_shen_2026_vinet_2` | Fe hcp | 16–140 GPa |

They are separate material/phase records rather than objects assigned a
special “calibrant” role. Load one through the ordinary material interface:

```python
from peritheos import Material, get_material_document

document = get_material_document("gold")
gold = Material.from_eosmat(
    document,
    record_identifiers=["gold_shen_2026_vinet_3"],
).eos_records[0]

pressure = gold.pressure(volume=60.0, temperature=300.0)
result = gold.pressure_with_uncertainty(volume=60.0, volume_sigma=0.01)
```

The paper reports fit errors for `K0` and `K0'`, but neither their confidence
level nor covariance. Peritheos preserves the printed errors without guessing
either quantity and states the independent-parameter assumption in propagated
results. These are 300 K isotherms, not thermal EOS records; supplying a
temperature other than 300 K or a temperature uncertainty is therefore
rejected.

The `gold_anderson_1989_bm3_1` material record is a thermal Au
parameterization over the calculation grid printed in Anderson et al. Table V:
`300-3000 K` and `0.66 <= V/V0 <= 1`. It uses the same ordinary material
loading interface shown above. Its correction is not a constant `alpha*K_T`
term and not a shifted reference state; the temperature slope varies as
`alpha_KT_ref + dK_dT_V*ln(V0/V)`. The source reports only a partial numerical
uncertainty for `dK_dT_V` and notes an additional unquantified contribution
from `K0'`, which remains visible in the record notes.

The Tange domain is the marginal envelope of several pressure-scale-free data
sets, not a rectangular guarantee that every combination of its extrema was
measured. Table 5 in that paper prints calculations to 4000 K, but values beyond
the 3700 K data limit are flagged as extrapolations and are not admitted by the
default range check.

The five Dorfman entries are one internally consistent, relative 300 K scale.
They are anchored to the Tange MgO scale; co-compression by itself does not make
an independent absolute pressure scale. The paper reports internal agreement to
3% up to 2.5 Mbar and estimates a 2–3% contribution from non-hydrostatic stress.
The per-entry ranges above follow the experimental run envelope in Table 1 and
cap the common claim at 250 GPa.

The eleven Sokolova records use every material column in Table 1. Dashes in
the anharmonic/electronic rows mean that the corresponding term is inactive;
the catalog encodes the disabled coefficient as zero and preserves that fact
in parameter provenance. The paper states a calculation range of at least
400 GPa and 3000 K, but that is not a claim that every point in the resulting
rectangle was measured. Its stated no-more-than-3–4% high-P-T marker
uncertainty has no confidence convention and is therefore recorded, not
silently converted into a one-sigma parameter covariance.
The individual Appendix A workbook files were not present beside the supplied
PDF during this audit. Table 1 and the calculation output printed in Figure 2
therefore provide the catalog's primary numerical provenance; no gaps were
filled from another pressure-standard library.

Fei et al. publish four thermal scales: Au, Pt, NaCl-B2, and Ne. MgO is a
consistency anchor in that paper, not a new Fei MgO parameter set. The Fei
records explicitly select `MieGruneisenDebye`'s `variable_exponent`
Debye-temperature law rather than its `integrated_gruneisen` default, because
the two formulas differ when $q\ne0$.

## Common DAC workflows

### B2 KCl as pressure medium and marker

The same primary scale is available through the compact EOS catalog and the
cross-compatible material library. The latter also carries the ideal B2
structure needed by Dioptas:

```python
from peritheos import Material, get_eos_record, get_material_document

volume = 30.0  # A^3 for the one-formula-unit B2 conventional cell
temperature = 2000.0  # K

scale = get_eos_record("kcl_b2_dewaele_2012")
pressure = scale.pressure(volume, temperature)

document = get_material_document("kcl")
material_scale = Material.from_eosmat(
    document,
    record_identifiers=["kcl_b2_dewaele_2012_vinet_3"],
).eos_records[0]
assert abs(material_scale.pressure(volume, temperature) - pressure) < 1e-12
```

The shared material record is the default in `kcl.eosmat`. Its 298 K Vinet
fit is experimental to 165 GPa; Equation (2)'s linear thermal-pressure term is
derived from molecular dynamics and is published for 300-7000 K, 0-200 GPa,
and `0.4 <= V/V0 <= 1`. The B2 lower bound remains 2.6 GPa. Because the paper
prints no coefficient errors or covariance, Peritheos propagates measured
volume and temperature errors but does not invent parameter uncertainty.

The same material also exposes
`kcl_campbell_1991_bm2_1` for workflows that intentionally use the older
Campbell--Heinz B2 isotherm. It is labeled as a composite, not as a verbatim
single-paper parameter set: Campbell and Heinz publish
`V0(B2)/V0(B1)=0.8483(57)`, `K0=28.7(6) GPa`, and fixed `K0'=4`; Peritheos
combines that ratio with Dewaele et al.'s experimental B1
`V0=62.36 angstrom^3` to obtain `V0(B2)=52.899988 angstrom^3`. The propagated
`0.355452 angstrom^3` uncertainty covers the ratio only, because the Dewaele
table gives no B1-`V0` error. The Dewaele P-V-T record remains the default.

### Sokolova 2016 markers

```python
from peritheos import get_eos_record

gold = get_eos_record("au_fcc_sokolova_2016")
pressure = gold.pressure(volume=55.0, temperature=2000.0)
recovered = gold.volume(pressure, temperature=2000.0)

assert abs(recovered - 55.0) < 1e-8
print(gold.reference.doi)
print(gold.parameter_provenance["QE1o"])
```

The public input is the conventional fcc cell volume in
$\mathring{\mathrm{A}}^3$, not the Table 1 molar volume. Peritheos performs the conversion at the
catalog boundary. Use `mgo_b1_sokolova_2016`, `diamond_sokolova_2016`,
`al_fcc_sokolova_2016`, `cu_fcc_sokolova_2016`,
`ag_fcc_sokolova_2016`, `au_fcc_sokolova_2016`,
`pt_fcc_sokolova_2016`, `nb_bcc_sokolova_2016`,
`ta_bcc_sokolova_2016`, `mo_bcc_sokolova_2016`, or
`w_bcc_sokolova_2016` to select the material explicitly.

### Compare named Au scales without conflating them

```python
from peritheos import get_eos_record

volume = 55.0
temperature = 1800.0

p_sokolova = get_eos_record("au_fcc_sokolova_2016").pressure(volume, temperature)
p_fei = get_eos_record("au_fcc_fei_2007").pressure(volume, temperature)
p_dorfman_300k = get_eos_record("au_fcc_dorfman_2012").pressure(volume, 300.0)
```

These are separate literature scales with different reference isotherms,
thermal terms, fit inputs, and validity envelopes. The 300 K-only Dorfman
entry cannot be evaluated at 1800 K.

### Hot MgO marker

```python
import numpy as np

from peritheos import get_eos_record

mgo = get_eos_record("mgo_b1_tange_2009_vinet")
volumes = np.array([69.0, 65.0, 60.0])  # A^3/conventional cell
temperatures = np.array([1200.0, 1800.0, 2500.0])  # K
pressures = mgo.pressure(volumes, temperatures)  # GPa
recovered = mgo.volume(pressures, temperatures)
```

The thermal term is a genuine Tange Mie-Gruneisen-Debye pressure relative to
the 300 K Vinet isotherm. It is not the empirical fractional DAC confinement
term described in [Advanced DAC analysis](dac-thermal-pressure.md), and it is
not a shock Hugoniot.

### Measurement and parameter uncertainty

```python
result = mgo.pressure_with_uncertainty(
    volume=67.23,
    temperature=3000.0,
    volume_sigma=0.01,
    temperature_sigma=50.0,
)

print(result.value, result.standard_error)
print(result.assumptions)
```

Tange et al. publish standard errors but no covariance matrix. Peritheos
therefore propagates those errors as independent and states that assumption in
the result. It does not add the paper's 0.8 GPa total RMS model residual as if it
were parameter variance.

Dorfman et al. warn that formal errors from their constrained fit are
artificially low and recommend approximately 2% realistic uncertainty for
$K_0$ and $K'_0$. The catalog propagates that recommendation independently
because no covariance was published. It records, but does not fold in, the
separate 2–3% non-hydrostatic-stress estimate.

Fei Table 1 gives parenthetical or plus-minus uncertainties but does not state
their confidence convention or publish covariance. The catalog propagates
them as independent one-standard-deviation errors and reports that assumption.
Sokolova Table 1 gives no individual parameter errors, so those entries
propagate measured volume and temperature uncertainty only; the paper's
aggregate 3–4% high-P-T estimate remains a separately reported scale caveat.

Temperature uncertainty is rejected for the 300 K-only entries. A zero
temperature derivative would falsely suggest those scales were insensitive to
temperature; they simply do not contain a thermal model.

When a paper supplies no parameter errors (the KCl/KBr fits), Peritheos still
propagates measured volume and temperature uncertainty and explicitly labels
the missing parameter-uncertainty block. It never substitutes zero errors.

### Explicit extrapolation

Catalog methods check their published envelope by default:

```python
gold = get_eos_record("au_fcc_dorfman_2012")
gold.pressure(gold.reference_volume)  # raises: below the calibrated range

# The underlying equation can still be evaluated for a deliberate diagnostic.
p0 = gold.pressure(gold.reference_volume, check_validity=False)
assert p0 == 0.0
```

`within_validity(volume, temperature)` checks states without raising. Disabling
the check does not expand the scientific validation, stabilize an absent phase,
or create a thermal correction.

## Executable `.eosmat` materials

`Material.to_dict()` and `Material.to_eosmat()` return the same canonical,
JSON-safe `peritheos.material` format-3 document used for file exchange with
Dioptas. There is no second public material format:

```python
import json

from peritheos import (
    Material,
    get_material,
)

gold = get_material("au_fcc")
with open("gold.eosmat", "w", encoding="utf-8") as stream:
    json.dump(gold.to_dict(), stream, indent=2)

with open("gold.eosmat", encoding="utf-8") as stream:
    loaded_gold = Material.from_dict(json.load(stream))
```

The document stores identity and public units once. Optional `symmetry`,
`lattice`, `space_group`, `atom_sites`, and `peaks` fields survive a Peritheos
round trip even though Peritheos does not use them numerically. Each EOS record
keeps:

- its stable identifier and display label;
- separate `eos` reference-isotherm and optional `thermal` components;
- explicit reference volume, public and internal model volume units, and the
  public-to-model conversion factor;
- primary reference, DOI, equation/table/page or supplement locations;
- component-specific parameter provenance, standard errors, and covariance;
- the published validity envelope, assumptions, ambiguities, and caveats.

The loader accepts only model identifiers in Peritheos's fixed registry. It
does not dynamically import an implementation path. `Material.from_eosmat()`
also fails closed when a record is `pending_primary_source_check` or `deferred`.
An explicit `require_primary_validation=False` is needed to construct such a
record, and `record_identifiers=(...)` can select a supported record from a
larger file. That opt-in does not promote its scientific status.

`to_snapshot_dict()` and the old `peritheos.material-snapshot` version-2 reader
remain temporarily available only for backward compatibility. New files must
use `.eosmat` format 3. See the [schema reference](eosmat-schema.md).

JCPDS remains a deliberately lossy legacy export. Its conventional BM3 and
simple thermal fields cannot represent DOI-level parameter provenance,
covariance, arbitrary thermal-model composition, or the validity semantics
above.

## Catalog inventory and selection

The inventory pass was a discovery exercise only. No BurnMan or Pytheos source,
tests, equations, parameters, metadata, or implementation structure were used.

- BurnMan's public catalog documentation names Anderson Au; Armentrout Co;
  Campbell, Dewaele, Pigott, and other Ni scales; Decker and Matsui NaCl B1;
  Dewaele KCl/KBr B2, corundum, and CsCl; Dorogokupets Ag, Al, Au, Cu,
  diamond, MgO, Pt, Ta, W, and Fe polymorphs; Fei Au, Pt, and hcp Fe; Holmes
  Pt; Huang/Sokolova/Zeng and Litasov Mo/W; Le Godec and Zhao hBN; Ono and Zha
  Re; Speziale and Tange MgO; and multiple Au/Pt/KCl alternatives. Some public
  entries explicitly combine sources, so each still needs a primary-paper
  audit before adoption.
- Pytheos's public documentation and example notebooks name Au scales from
  Jamieson, Tsuchiya, Fei, Dorogokupets, Yokoo, and Ye; a Yokoo Pt scale; and
  Jamieson and Tange MgO scales. Its public API also documents constant-$q$,
  Speziale, Dorogokupets, and Tange thermal families.
- Dioptas 0.10.0 supplies 120 material documents and 147 publication-attributed
  EOS records through an offline JSON database and `.eosmat` exchange files.
  Its phase-library workflow makes Au, Pt, MgO, NaCl, Ne, Re, Ar, Ag, Mo, and
  common calibrants especially valuable, but its parameters and regression
  values are not treated as primary evidence here.

The expanded Peritheos catalog chooses Tange MgO because the paper supplies an
independent P–V–T analysis, errors, a data envelope, and a printed numerical
table. Dorfman then adds five high-value, phase-explicit 300 K standards in one
primary co-compression framework. Direct primary-source additions cover LiF,
both NaCl phases, both KCl and KBr phases, c-BN, diamond, Ag, and Ni. The
complete Sokolova Table 1 contributes eleven thermal markers, and Fei et al.
contribute four internally consistent thermal scales.

The Dorogokupets and Sokolova names require special care.
`MultiOscillatorGruneisenThermalEOS` is the reusable equation implementation
used for the 2016 spreadsheet model, which is based on a *modified*
Dorogokupets–Oganov formalism. Author/year remains in material record IDs such
as `au_fcc_sokolova_2016`, not in the equation model name. A 2007 and a 2016
pressure scale may share an equation family but must retain separate scale
identifiers, parameters, citations, and provenance; they are not blind aliases.

## Deferred entries

The quasi-hydrostatic 300 K Re scale is now included. A high-temperature Re
extension and other high-value candidates such as hcp Fe, corundum, hBN, and
newer Au/Pt scales remain inventory items until the same primary-evidence gate
is completed. `DEFERRED_EOS_RECORDS` records the specific unresolved source
issue for candidates that have reached implementation review.

## Numerical validation

The Tange implementation is regressed directly against printed Table 5 values
over 300–3000 K and $0.65 \leq V/V_0 \leq 1.00$. Values match the table's
0.01 GPa precision. Reference-state identities, analytic Grüneisen/Debye
consistency, scalar and array behavior, inverse round trips, invalid states,
validity enforcement, and V/T/parameter uncertainty are tested separately.
Additional regressions use LiF Table 4, diamond Table I, and Ag/Ni Table II;
the tolerances retain the published experimental residuals rather than forcing
measured data to lie exactly on a fitted curve. KCl/KBr tests reproduce the
printed equation (2) thermal increment exactly.

The Re Vinet record is regressed against four lattice-parameter rows from
Anzellini et al. Table III. Its reported parameter-error bars are explicitly
stored as 95% confidence half-widths and converted to normal-equivalent
standard errors only at uncertainty-propagation time.

The same confidence metadata is retained for the primary-source fits of Al,
Cu, W, Ni, Ag, diamond, alpha- and omega-Ti, Si-V, Si-VII, Si-X, corundum,
and LiF whose tables also define their quoted errors as 95% intervals. A
catalog-wide regression loads every such `.eosmat` record and checks the
conversion rather than relying only on the raw JSON value.

All eleven Sokolova records are checked against the parameters printed in
Table 1. MgO is additionally regressed at eight compressions against the
pressure-calculation output printed in Figure 2, including the 300 K versus
298.15 K reference-temperature difference. Existing diamond cases reproduce
the accompanying workbook calculation. The Fei tests verify all four Table 1
parameter sets and equation 3's publication-specific Debye-temperature
convention, then exercise reference states and P-V-T round trips.

No BurnMan or Pytheos numerical output is used as a validation baseline. A
versioned public-API-only comparison is recorded under
[External black-box comparisons](validation.md#external-black-box-comparisons)
to expose convention differences; agreement remains secondary evidence rather
than scientific validation.

The paper's 2010 correction, [doi:10.1029/2010JB007959](https://doi.org/10.1029/2010JB007959),
was also checked. It corrects only a phase-boundary plot in Figure 11 and does
not alter the EOS equations, Table 4 parameters, or Table 5 regression values
used here.
