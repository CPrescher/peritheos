# Jacobsen et al. (2008) MgO EOS audit

## Outcome

[Jacobsen et al. (2008)](https://doi.org/10.2138/am.2008.2988) supports two
executable, pressure-scale-specific room-temperature BM3 records for B1 MgO:

| Record | Medium and pressure coordinate | V0 (A^3) | K0 (GPa) | K0' | Fitted range (GPa) |
|---|---|---:|---:|---:|---:|
| `mgo_jacobsen_2008_bm3_kcl_mao1978` | KCl; Mao et al. (1978) ruby | 74.698(7) | 164.1(9) | 4.05(4) | 0-86.6 |
| `mgo_jacobsen_2008_bm3_helium_mao1986` | He; Mao et al. (1986) ruby plus six Mao-linked diamond-Raman pressures | 74.687(6) | 159.6(6) | 3.74(3) | 0-111.0 |

They are not interchangeable fits. The paper uses their difference to argue
that either MgO compression changes above about 40 GPa in helium or the
Mao-1986 ruby scale requires adjustment for the softer helium medium.

No third Jacobsen MgO EOS was added for the `MgO-scale P` column in Table 1.
Those pressures were calculated from the primary Zha et al. (2000) MgO scale
and are therefore derived values on an existing external parameterization,
not a new fit reported by Jacobsen et al.

## Source and material identity

The canonical DOI is `10.2138/am.2008.2988`, American Mineralogist 93,
1823-1828. The final six-page article and its complete Tables 1-2 were checked;
the journal and author article records expose no separate official supplement.

The samples are commercially sourced, greater-than-99.95% MgO single-crystal
plates. The paper states that MgO remains in the B1 halite structure over the
represented range. Up to eight `hk0` reflections were used in cubic-cell
refinement. The tabulated volumes are conventional B1-cell volumes in A^3;
Peritheos therefore uses the existing `Fm-3m`, Z=4 MgO material identity
without conversion.

The EOS is the standard third-order Birch-Murnaghan form. The paper supplies
the defining finite-strain convention directly:

```text
f = ((V0/V)^(2/3) - 1) / 2
F = P / (3 f (1 + 2 f)^(5/2))
F = K0 + 3 K0 (K0' - 4) f / 2
```

This is exactly Peritheos' `birch_murnaghan_3` convention. Both fits vary V0,
K0, and K0'; neither preferred record has a fixed EOS coefficient. The nominal
reference temperature is 300 K because these are unheated, room-temperature
static-compression experiments. It is not a row-wise measured temperature.

The source reports parenthesized parameter errors but does not state their
confidence convention or provide parameter covariance. Peritheos preserves
the errors, sets `parameter_error_confidence` and `parameter_covariance` to
null, and does not manufacture correlations.

## Primary observations and calibration

The bundled Table 1 resource contains 52 rows: one ambient-volume anchor and
51 helium observations from runs 1-3. It preserves the normalized ruby R1
shift, pressure and volume standard deviations, and the calculated Zha-2000
pressure column. Forty-five high-pressure rows use the executable Mao-1986
ruby calibration. The last six use Equation 4a of Sun et al. (2005), calibrated
against Mao-1986; their Raman shifts are not printed, so those six pressures
cannot be recalculated exactly from primary observables.

The bundled Table 2 resource contains all 26 KCl-medium observations and their
printed pressure and volume standard deviations. Every pressure uses the
executable Mao-1978 ruby scale. The reported source-scale pressures can be
inverted to normalized wavelength ratios, so this series is ready for
pressure-scale conversion.

The article also proposes a helium ruby calibration. It reports both a
two-parameter fit, A=1930(6) GPa and B=9.44(23), and its preferred Mao-form
fixed-A variant, A=1904 GPa and B=10.32(7), over 23-140 GPa. No calibration
record was added in this focused EOS change because the current calibration
interchange requires an absolute default ambient R1 wavelength, while this
paper defines and tabulates only the normalized shift. Inventing that missing
default would violate the primary-source rule; both calibration variants
remain documented here.

## Numerical reproduction

The stored published curves reproduce independent high-pressure table rows:

| Record | Table state | Published pressure (GPa) | Peritheos pressure (GPa) |
|---|---:|---:|---:|
| KCl/Mao-1978 | V=55.792 A^3 | 86.6 | 86.6703 |
| He/Mao-1986 | V=52.239 A^3 | 111.0 | 110.7638 |

An independent errors-in-variables refit used every non-ambient row, both
printed coordinate standard deviations, and the independently measured
V0=74.698 A^3 as a fixed anchor:

| Series | Refit K0 (GPa) | Refit K0' | Published K0 (GPa) | Published K0' |
|---|---:|---:|---:|---:|
| KCl | 164.1036 | 4.05239 | 164.1(9) | 4.05(4) |
| He | 159.3536 | 3.74517 | 159.6(6) | 3.74(3) |

Both independent results agree with the published values inside the reported
errors. They are validation refits, not replacement EOS records. A bit-for-bit
free-V0 regeneration is impossible because the source does not identify its
fitting software, exact objective, ambient-anchor weighting, covariance
convention, or parameter correlations.

## Resolved source inconsistency and omitted alternatives

The abstract gives helium `V0=74.697(6) A^3`, but the Results section explicitly
assigns `V0=74.687(6) A^3` to the three-parameter fit of all three runs. The
same Results paragraph separately gives `74.697(7) A^3` after excluding the six
diamond-Raman rows. The executable all-data record therefore uses the Results
value `74.687(6) A^3` and stores the abstract discrepancy explicitly.

The same paragraph prints `GPa` after the fixed-K0 sensitivity fit's
`V0=74.695(6)` value. This is retained as a documented source typo and resolved
to A^3 because V0 is the volume parameter and both neighboring fits use A^3.
The ruby-calibration discussion likewise attaches `GPa` to its dimensionless
free-fit B coefficient; the documented value remains dimensionless as required
by the printed power-law equation.

The paper's remaining parameter sets are sensitivity analyses, not additional
preferred records: KCl without the ambient anchor, KCl with K0' fixed to the
Zha value, helium without the six Raman rows, and helium with K0 fixed to the
Zha value. Their coefficients and fixed/free status remain in the record's
`scientific_validation.reported_parameterizations` metadata.
