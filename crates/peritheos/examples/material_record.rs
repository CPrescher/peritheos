use peritheos::load_eosmat_str;

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
    "eos": {
      "type": "BM3",
      "model": "birch_murnaghan_3",
      "parameters": {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0}
    },
    "parameter_errors": {},
    "fixed_parameters": [],
    "scientific_validation": {"status": "example_only"}
  }]
}"#;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let material = load_eosmat_str(MATERIAL)?;
    let record = material
        .default_record()
        .expect("example has a default record");
    let pressure = record.pressure(9.0, 300.0)?;
    let recovered = record.volume(pressure, 300.0)?;

    println!("{} / {}: {pressure:.3} GPa", material.name, record.label);
    assert!((recovered - 9.0).abs() < 1.0e-9);
    Ok(())
}
