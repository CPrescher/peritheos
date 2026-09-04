//! Shock Hugoniot equations of state.

use crate::validation::{finite_parameter, finite_result, positive_parameter, positive_state};
use crate::{EosError, EosResult};

/// Common behavior of a pressure-volume shock Hugoniot.
pub trait Hugoniot {
    /// Initial/reference volume in any consistent volume unit.
    fn reference_volume(&self) -> f64;

    /// Initial density in g cm^-3.
    fn initial_density(&self) -> f64;

    /// Initial pressure in `GPa`.
    fn initial_pressure(&self) -> f64;

    /// Hugoniot pressure in `GPa` at a compressed volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn pressure(&self, volume: f64) -> EosResult<f64>;

    /// Compressed volume on the Hugoniot at a pressure in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error when pressure lies outside the invertible branch.
    fn volume(&self, pressure: f64) -> EosResult<f64>;

    /// Particle velocity in km s^-1 at a compressed volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn particle_velocity(&self, volume: f64) -> EosResult<f64>;

    /// Shock velocity in km s^-1 at a compressed volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn shock_velocity(&self, volume: f64) -> EosResult<f64>;

    /// Density in g cm^-3 at a compressed volume.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn density(&self, volume: f64) -> EosResult<f64>;

    /// Specific internal-energy increase in MJ kg^-1.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn specific_internal_energy_change(&self, volume: f64) -> EosResult<f64>;

    /// Tangent stiffness `-V dP_H/dV` along the Hugoniot in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error when the volume is outside the compressive branch.
    fn tangent_modulus(&self, volume: f64) -> EosResult<f64>;
}

/// Principal or branch Hugoniot defined by `U_s = c0 + s u_p`.
///
/// Density is in g cm^-3 and velocities are in km s^-1, so the momentum
/// relation yields pressure in `GPa` without an additional unit conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearUsUpHugoniot {
    /// Initial/reference volume.
    pub v0: f64,
    /// Initial density in g cm^-3.
    pub rho0: f64,
    /// Zero-particle-velocity intercept in km s^-1.
    pub c0: f64,
    /// Dimensionless linear slope.
    pub s: f64,
    /// Initial pressure in `GPa`.
    pub p0: f64,
}

impl LinearUsUpHugoniot {
    /// Construct a linear shock-velocity--particle-velocity Hugoniot.
    ///
    /// # Errors
    ///
    /// Returns an error unless `v0`, `rho0`, `c0`, and `s` are positive and
    /// finite and `p0` is finite.
    pub fn new(v0: f64, rho0: f64, c0: f64, s: f64, p0: f64) -> EosResult<Self> {
        Ok(Self {
            v0: positive_parameter(v0, "V0")?,
            rho0: positive_parameter(rho0, "rho0")?,
            c0: positive_parameter(c0, "c0")?,
            s: positive_parameter(s, "s")?,
            p0: finite_parameter(p0, "P0")?,
        })
    }

    fn compression(&self, volume: f64) -> EosResult<(f64, f64)> {
        let volume = positive_state(volume, "volume")?;
        if volume > self.v0 {
            return Err(EosError::InvalidState {
                name: "volume",
                reason: "must not exceed V0 on a compression Hugoniot",
            });
        }
        let mu = 1.0 - volume / self.v0;
        let denominator = 1.0 - self.s * mu;
        if denominator <= 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        Ok((mu, denominator))
    }

    /// Shock velocity in km s^-1 for a particle velocity in km s^-1.
    ///
    /// # Errors
    ///
    /// Returns an error unless particle velocity is finite and non-negative.
    pub fn shock_velocity_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        let particle_velocity = finite_parameter(particle_velocity, "particle_velocity")?;
        if particle_velocity < 0.0 {
            return Err(EosError::InvalidState {
                name: "particle_velocity",
                reason: "must be non-negative",
            });
        }
        finite_result(self.c0 + self.s * particle_velocity)
    }

    /// Pressure in `GPa` for a particle velocity in km s^-1.
    ///
    /// # Errors
    ///
    /// Returns an error unless particle velocity is finite and non-negative.
    pub fn pressure_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        let shock_velocity = self.shock_velocity_from_particle_velocity(particle_velocity)?;
        finite_result(self.p0 + self.rho0 * shock_velocity * particle_velocity)
    }

    /// Volume for a particle velocity in km s^-1.
    ///
    /// # Errors
    ///
    /// Returns an error unless particle velocity is finite, non-negative, and
    /// maps to a positive volume.
    pub fn volume_from_particle_velocity(&self, particle_velocity: f64) -> EosResult<f64> {
        let shock_velocity = self.shock_velocity_from_particle_velocity(particle_velocity)?;
        let compression = particle_velocity / shock_velocity;
        let volume = self.v0 * (1.0 - compression);
        if volume <= 0.0 || !volume.is_finite() {
            return Err(EosError::OutsideInvertibleRange);
        }
        Ok(volume)
    }
}

impl Hugoniot for LinearUsUpHugoniot {
    fn reference_volume(&self) -> f64 {
        self.v0
    }

    fn initial_density(&self) -> f64 {
        self.rho0
    }

    fn initial_pressure(&self) -> f64 {
        self.p0
    }

    fn pressure(&self, volume: f64) -> EosResult<f64> {
        let (mu, denominator) = self.compression(volume)?;
        finite_result(self.p0 + self.rho0 * self.c0 * self.c0 * mu / denominator.powi(2))
    }

    fn volume(&self, pressure: f64) -> EosResult<f64> {
        let pressure = finite_parameter(pressure, "pressure")?;
        if pressure < self.p0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        let y = (pressure - self.p0) / (self.rho0 * self.c0 * self.c0);
        let mu = if y == 0.0 {
            0.0
        } else {
            2.0 * y / (1.0 + 2.0 * self.s * y + (1.0 + 4.0 * self.s * y).sqrt())
        };
        let volume = self.v0 * (1.0 - mu);
        if volume <= 0.0 || !volume.is_finite() || 1.0 - self.s * mu <= 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        Ok(volume)
    }

    fn particle_velocity(&self, volume: f64) -> EosResult<f64> {
        let (mu, denominator) = self.compression(volume)?;
        finite_result(self.c0 * mu / denominator)
    }

    fn shock_velocity(&self, volume: f64) -> EosResult<f64> {
        let (_, denominator) = self.compression(volume)?;
        finite_result(self.c0 / denominator)
    }

    fn density(&self, volume: f64) -> EosResult<f64> {
        self.compression(volume)?;
        finite_result(self.rho0 * self.v0 / volume)
    }

    fn specific_internal_energy_change(&self, volume: f64) -> EosResult<f64> {
        let density = self.density(volume)?;
        let pressure = self.pressure(volume)?;
        finite_result(0.5 * (pressure + self.p0) * (1.0 / self.rho0 - 1.0 / density))
    }

    fn tangent_modulus(&self, volume: f64) -> EosResult<f64> {
        let (mu, denominator) = self.compression(volume)?;
        finite_result(
            (1.0 - mu) * self.rho0 * self.c0 * self.c0 * (1.0 + self.s * mu) / denominator.powi(3),
        )
    }
}
