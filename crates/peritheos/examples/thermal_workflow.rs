use peritheos::{isothermal::BM3, thermal::MieGruneisenDebye, CaloricEos, ThermalEos};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Energy-based thermal models use molar volume in J bar^-1 mol^-1.
    let reference = BM3::new(1.02, 165.0, 5.0)?;
    let eos = MieGruneisenDebye::new(reference, 300.0, 170.0, 2.9, 1.0, 1.0)?;
    let volume = 0.95;
    let temperature = 1_500.0;

    let pressure = eos.pressure(volume, temperature)?;
    let recovered_volume = eos.volume(pressure, temperature)?;
    let recovered_temperature = eos.temperature(pressure, volume)?;

    println!("P({volume:.3}, {temperature:.0} K) = {pressure:.3} GPa");
    println!(
        "Cv = {:.3} J mol^-1 K^-1",
        eos.molar_heat_capacity_v(volume, temperature)?
    );
    println!(
        "Cp = {:.3} J mol^-1 K^-1",
        eos.molar_heat_capacity_p(volume, temperature)?
    );

    assert!((recovered_volume - volume).abs() < 1.0e-9);
    assert!((recovered_temperature - temperature).abs() < 1.0e-6);
    Ok(())
}
