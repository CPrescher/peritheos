# Phase E database integration — G06/C11

Integrated 2026-10-16 as `phase_e_fe0072`, with published/default record
`phase_e_fe0072_crichton_ross_2000_bm3`. The bundled card is
`peritheos/data/materials/phase_e_fe0072.eosmat`; normal material and record
lookup load it directly. The source audit, observations, independent diagnostic
fits and [residual plot](../data/phase-e-crichton-ross-2000-parity.png) are
retained. See the [complete scientific audit](phase-e-crichton-ross-2000.md).

## Accepted scientific assessment

**The published physical volume EOS is successfully recovered: parity within
reported uncertainty.** The independent weighted refit shifts each coefficient
by less than 0.123 times its published error width and changes the pressure
curve by at most about 0.0096 GPa over the measured volume interval.
Published coefficients remain unchanged. All eleven primary P,a,c,V rows and
reported errors are bundled; the differing thesis version is retained separately.

The common record-level refit ledger classifies this result as `parity`, with
`parity_basis: within_reported_parameter_errors`. The paper-level ledger
classifies Crichton and Ross (2000) as `Reproduced`. Exact original workflow
availability remains separate provenance, not a failed-fit classification.
Diagnostic refits are not additional catalog records.

## Identity, diffraction and calibration

The material is the measured Mg1.96(7)Fe0.072(5)Si1.04(5)H3.7(8)O6 sample,
not an idealized Mg endmember. Its R-3m, No. 166 structure uses conventional
hexagonal axes and Z=1 on the empirical six-O basis. Hydrogen is inferred
from the EMPA deficit; target-sample sites and Fe partition are unresolved.

The documented diffraction fallback transfers eight observed phase E hkl
identities from Shieh (2000), with d spacings calculated from the target cell.
Uniform markers support peak positions and indexing, without asserting
quantitative intensities, hydrogen positions or local disorder. No occupied
sites from a different composition are assigned to this sample. Axial fits
remain supporting source numbers, separate from the executable volume EOS.

The pressure scale links the existing `alpha_quartz_angel_1997_bm3_1` standard;
observation-level pressure
recalculation is `missing_calibrant_observations`. This absence is recorded
as optional workflow/calibration evidence in the dataset-request page and
does not negate the accepted fit parity. No calibration registry changes
are needed. Source reuse terms remain unspecified, with any contributor CC0
limited to transcription, normalization, metadata and arrangement rights.
Thermal data for the distinct Wang (2022) Fe-free sample remain separate.

## Catalog and documentation updates

Phase E adds one stored material document, one executable primary-source-validated
record, one distinct dataset, eleven retained rows, one common-ledger parity
result and one reproduced primary paper. The validated local integration snapshot (including other accepted additions) has:

| Measure | Count |
|---|---:|
| Stored material documents | 231 |
| Stored EOS records | 659 |
| Executable materials | 228 |
| Executable EOS records | 656 |
| Primary-source-validated records | 649 |

The local manifest, primary-source audit, generated record and paper ledgers,
coverage snapshot, candidate backlog, source-audit index, navigation,
catalog documentation, dataset requests and changelog reflect the addition.
The scoped phase E commit preserves unrelated shared work: catalog-wide count
reconciliation and the new coverage infrastructure remain separate checkout
changes. Its committed scientific ledgers add only the phase E record and
leave all pre-existing record assessments unchanged.
The common-ledger adapter calls the same dedicated reproduction used by
the material tests. Other scientific records and their classifications are
preserved.

## Reproducible checks

Verification on 2026-10-16: all **585 material, catalog, calibration, ledger
and documentation-consistency tests passed** after integration. Source-PDF
hashes, numerical report freshness, source-index/paper-ledger/coverage checks,
Ruff, formatting and Python type checks passed. The Rust workspace tests,
formatting and Clippy passed; the installed phase E Rust test also passed
after final metadata adjustments. Combined Python coverage was 90.72%.

The latest global strict documentation build reports one unrelated link in
`docs/research/ferrous-ringwoodite-armentrout-2011-integration.md` pointing
outside the documentation tree to a ringwoodite staging image. Phase E's
local audit/integration/image links resolve, and its pages and residual plot
are present in the built documentation. The separate ringwoodite work was
not changed by this integration.

```sh
uv run python scripts/reproduce_phase_e_crichton_ross_2000.py --check
uv run --group docs python scripts/reproduce_phase_e_crichton_ross_2000.py --check --plot
uv run pytest -q tests/test_phase_e_crichton_ross_2000.py -W error
cargo test -p peritheos --test phase_e_crichton_ross_2000 --locked
uv run python scripts/validate_primary_eos_refits.py --check
uv run python scripts/generate_paper_investigation_ledger.py --check
uv run python -m scripts.generate_catalog_coverage --check
uv run python scripts/generate_source_audit_index.py --check
```

The Python and Rust tests load the installed card rather than a staging copy.
Optional `--source-dir PATH` verifies all six retrieved source PDFs against
their recorded hashes. PDFs are not redistributed. No authors were contacted;
integration does not commit, push or publish the checkout.
