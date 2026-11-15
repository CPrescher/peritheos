use peritheos::{export_eosmat_v3, load_eosmat_str, load_eosmat_str_with_audit, Reproduction};

/// A fictional EOSMAT v5 record file; no experimental validation is claimed.
const RECORD_FILE: &str = r#"{
  "format_version": 5,
  "identifier": "example_v5",
  "formula": "X",
  "eos_records": [{
    "identifier": "example_v5_bm3",
    "model": "BM3",
    "temperature_ref": 300,
    "parameters": {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0},
    "reference": "Synthetic documentation example",
    "check": {
      "reproduction": "close",
      "warning": null,
      "source": "deferred",
      "audit": "../audits/example_v5.audit.json#example_v5_bm3"
    }
  }]
}"#;

/// The matching fictional audit sidecar.
const AUDIT_SIDECAR: &str = r#"{
  "format": "peritheos.eosmat-audit",
  "format_version": 5,
  "material": "example_v5",
  "records": {
    "example_v5_bm3": {
      "evidence": {
        "scientific_validation": {"status": "deferred"},
        "reproduction": {
          "assessment": {
            "status": "similar",
            "scope": "parameter_agreement",
            "basis": "Synthetic refit used only to illustrate the format."
          }
        }
      }
    }
  }
}"#;

const MATERIAL: &str = r#"{
  "format": "peritheos.material",
  "format_version": 3,
  "identifier": "example_material",
  "name": "Example material",
  "formula": "X",
  "units": {
    "pressure": "GPa",
    "temperature": "K",
    "volume": "angstrom^3/conventional_unit_cell"
  },
  "eos_records": [{
    "identifier": "example_bm3",
    "label": "Example BM3",
    "default": true,
    "reference": "Synthetic documentation example",
    "eos": {
      "type": "BM3",
      "model": "birch_murnaghan_3",
      "parameters": {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0}
    },
    "parameter_errors": {},
    "fixed_parameters": [],
    "scientific_validation": {"status": "deferred"}
  }]
}"#;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let material = load_eosmat_str(MATERIAL)?;
    let record = material
        .default_equilibrium_record()
        .expect("example has a default record");
    let pressure = record.pressure(9.0, 300.0)?;
    let recovered = record.volume(pressure, 300.0)?;

    println!("{} / {}: {pressure:.3} GPa", material.name, record.label);
    assert!((recovered - 9.0).abs() < 1.0e-9);

    // Serialization validates the retained document and preserves extensions.
    let serialized = material.to_json()?;
    let round_tripped = load_eosmat_str(&serialized)?;
    assert_eq!(round_tripped.document, material.document);
    // Without any assessment, the headline is untested rather than positive.
    assert_eq!(record.reproduction(), Reproduction::Untested);

    // An EOSMAT v5 record file keeps a short check headline; its evidence may
    // live in an audit sidecar. A string load has no location to resolve the
    // sidecar from, so the headline is available and the audit is not.
    let record_file = load_eosmat_str(RECORD_FILE)?;
    let record = &record_file.eos_records[0];
    println!(
        "{}: reproduction {}, warning {:?}",
        record.identifier,
        record.reproduction(),
        record.warning()
    );
    assert!(record.audit_is_external() && record.audit().is_none());

    // Supplying the sidecar restores the complete record and its evidence.
    let complete = load_eosmat_str_with_audit(RECORD_FILE, AUDIT_SIDECAR)?;
    let record = &complete.eos_records[0];
    assert_eq!(record.reproduction(), Reproduction::Close);
    let audit = record.audit().expect("sidecar was merged");
    println!(
        "assessment basis: {}",
        audit["evidence"]["reproduction"]["assessment"]["basis"]
    );
    // v3 exports carry the evidence inline and omit the v5-only check block.
    let legacy = export_eosmat_v3(&complete.document)?;
    assert!(legacy["eos_records"][0].get("check").is_none());
    Ok(())
}
