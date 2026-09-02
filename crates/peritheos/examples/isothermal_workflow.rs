use peritheos::{batch::IsothermalEosBatch, isothermal::BM3, IsothermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // V0 and every evaluated volume share a unit; K0 and pressure use GPa.
    let eos = BM3::new(10.0, 160.0, 4.2)?;
    let volumes = [10.0, 9.5, 9.0, 8.5];
    let pressures = eos.pressures(&volumes)?;

    println!("volume\tpressure (GPa)\tbulk modulus (GPa)");
    for (&volume, &pressure) in volumes.iter().zip(&pressures) {
        println!(
            "{volume:.3}\t{pressure:.3}\t\t{:.3}",
            eos.bulk_modulus(volume)?
        );
    }

    let recovered = eos.volumes(&pressures)?;
    assert!(recovered
        .iter()
        .zip(volumes)
        .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9));
    Ok(())
}
