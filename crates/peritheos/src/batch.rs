//! Dependency-free ordered batch evaluation for native EOS models.

use crate::{CaloricEos, EosError, EosResult, IsothermalEos, ThermalEos};

fn matching_lengths(left: usize, right: usize) -> EosResult<()> {
    if left == right {
        Ok(())
    } else {
        Err(EosError::InvalidState {
            name: "batch",
            reason: "input slices must have matching lengths",
        })
    }
}

fn map_values<F>(values: &[f64], mut evaluate: F) -> EosResult<Vec<f64>>
where
    F: FnMut(f64) -> EosResult<f64>,
{
    values.iter().map(|&value| evaluate(value)).collect()
}

fn map_pairs<F>(left: &[f64], right: &[f64], mut evaluate: F) -> EosResult<Vec<f64>>
where
    F: FnMut(f64, f64) -> EosResult<f64>,
{
    matching_lengths(left.len(), right.len())?;
    left.iter()
        .zip(right)
        .map(|(&left, &right)| evaluate(left, right))
        .collect()
}

/// Ordered batch operations available for every isothermal EOS.
pub trait IsothermalEosBatch: IsothermalEos {
    /// Evaluate pressure for every volume in input order.
    ///
    /// # Errors
    ///
    /// Returns the first scalar model error in input order.
    fn pressures(&self, volumes: &[f64]) -> EosResult<Vec<f64>> {
        map_values(volumes, |volume| self.pressure(volume))
    }

    /// Evaluate isothermal bulk modulus for every volume in input order.
    ///
    /// # Errors
    ///
    /// Returns the first scalar model error in input order.
    fn bulk_moduli(&self, volumes: &[f64]) -> EosResult<Vec<f64>> {
        map_values(volumes, |volume| self.bulk_modulus(volume))
    }

    /// Invert every pressure on the model's documented branch.
    ///
    /// # Errors
    ///
    /// Returns the first scalar inversion error in input order.
    fn volumes(&self, pressures: &[f64]) -> EosResult<Vec<f64>> {
        map_values(pressures, |pressure| self.volume(pressure))
    }
}

impl<T: IsothermalEos + ?Sized> IsothermalEosBatch for T {}

/// Ordered batch operations available for every thermal EOS.
pub trait ThermalEosBatch: ThermalEos {
    /// Evaluate thermal pressure for paired volume-temperature states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn thermal_pressures(&self, volumes: &[f64], temperatures: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.thermal_pressure(volume, temperature)
        })
    }

    /// Evaluate reference-relative thermal-pressure increments for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn thermal_pressure_increments(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.thermal_pressure_increment(volume, temperature)
        })
    }

    /// Evaluate total pressure for paired volume-temperature states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn pressures(&self, volumes: &[f64], temperatures: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.pressure(volume, temperature)
        })
    }

    /// Invert paired pressure-temperature states to volume.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar
    /// inversion error in input order.
    fn volumes(&self, pressures: &[f64], temperatures: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(pressures, temperatures, |pressure, temperature| {
            self.volume(pressure, temperature)
        })
    }

    /// Predict heated volumes from paired cold pressures and temperatures.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid confinement fraction, mismatched slice
    /// lengths, or the first scalar inversion error in input order.
    fn volumes_with_dac_confinement(
        &self,
        cold_pressures: &[f64],
        temperatures: &[f64],
        f_dac: f64,
    ) -> EosResult<Vec<f64>> {
        map_pairs(
            cold_pressures,
            temperatures,
            |cold_pressure, temperature| {
                self.volume_with_dac_confinement(cold_pressure, temperature, f_dac)
            },
        )
    }

    /// Invert paired pressure-volume states to temperature.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar
    /// inversion error in input order.
    fn temperatures(&self, pressures: &[f64], volumes: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(pressures, volumes, |pressure, volume| {
            self.temperature(pressure, volume)
        })
    }

    /// Evaluate isothermal bulk modulus for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid step size, mismatched slice lengths, or
    /// the first scalar model error in input order.
    fn bulk_moduli(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
        relative_step: f64,
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.bulk_modulus(volume, temperature, relative_step)
        })
    }

    /// Evaluate isothermal compressibility for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn isothermal_compressibilities(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.isothermal_compressibility(volume, temperature)
        })
    }

    /// Evaluate volumetric thermal expansivity for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid step size, mismatched slice lengths, or
    /// the first scalar model error in input order.
    fn thermal_expansivities(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
        relative_step: f64,
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.thermal_expansivity(volume, temperature, relative_step)
        })
    }

    /// Infer temperature from paired ambient and heated volumes.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid confinement fraction, mismatched slice
    /// lengths, or the first scalar inversion error in input order.
    fn temperatures_from_volumes(
        &self,
        ambient_volumes: &[f64],
        heated_volumes: &[f64],
        f_dac: f64,
    ) -> EosResult<Vec<f64>> {
        map_pairs(
            ambient_volumes,
            heated_volumes,
            |ambient_volume, heated_volume| {
                self.temperature_from_volumes(ambient_volume, heated_volume, f_dac)
            },
        )
    }
}

impl<T: ThermalEos + ?Sized> ThermalEosBatch for T {}

/// Ordered caloric batch operations for models with a defined potential.
pub trait CaloricEosBatch: CaloricEos {
    /// Evaluate constant-volume molar heat capacity for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn molar_heat_capacities_v(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.molar_heat_capacity_v(volume, temperature)
        })
    }

    /// Evaluate constant-pressure molar heat capacity for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn molar_heat_capacities_p(
        &self,
        volumes: &[f64],
        temperatures: &[f64],
    ) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.molar_heat_capacity_p(volume, temperature)
        })
    }

    /// Evaluate thermodynamic Gruneisen parameter for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn gruneisen_parameters(&self, volumes: &[f64], temperatures: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.gruneisen_parameter(volume, temperature)
        })
    }

    /// Evaluate adiabatic bulk modulus for paired states.
    ///
    /// # Errors
    ///
    /// Returns an error for mismatched slice lengths or the first scalar model
    /// error in input order.
    fn adiabatic_bulk_moduli(&self, volumes: &[f64], temperatures: &[f64]) -> EosResult<Vec<f64>> {
        map_pairs(volumes, temperatures, |volume, temperature| {
            self.adiabatic_bulk_modulus(volume, temperature)
        })
    }
}

impl<T: CaloricEos + ?Sized> CaloricEosBatch for T {}
