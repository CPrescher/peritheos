//! Literal Maltby 2024 equations (DOI 10.1063/5.0237497), NOT validated.
//!
//! The source does not specify the CSM cutoff; callers must choose it explicitly.
//! The Table 8 volume is not reproduced. See the study reproduction documentation.
//! Volume uses J/bar/mol, pressure `GPa`, and unshifted molar energy J/mol.

use crate::thermal::debye_function_3;
use crate::validation::{finite_result, finite_state, positive_state};
use crate::{EosError, EosResult};

const R: f64 = 8.314_462_618;
const NA: f64 = 0.602_214_076;
const VREF: f64 = 22.56;

/// Experimental literal coefficient implementation, excluded from EOSMAT.
#[derive(Clone, Debug, PartialEq)]
pub struct Maltby2024Published {
    /// Squared cutoff radius in units of the nearest-neighbor distance.
    shell_cutoff_squared: u32,
    /// Outer zero of the Buckingham potential, in angstroms.
    sigma_angstrom: f64,
    shells: Vec<(f64, f64)>,
}

impl Maltby2024Published {
    /// Construct with an explicit occupied fcc shell cutoff, between 1 and 256.
    ///
    /// # Errors
    /// Returns an error for an invalid cutoff. No default cutoff is implied.
    // The validated cutoff bounds every integer/float conversion below by 512.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    pub fn new(shell_cutoff_squared: u32) -> EosResult<Self> {
        if !(1..=256).contains(&shell_cutoff_squared) {
            return Err(EosError::InvalidParameter {
                name: "shell_cutoff_squared",
                reason: "must lie in [1, 256]",
            });
        }
        let extent = (2.0 * f64::from(shell_cutoff_squared)).sqrt().ceil() as i32;
        let mut counts = vec![0_u32; shell_cutoff_squared as usize + 1];
        for i in -extent..=extent {
            for j in -extent..=extent {
                for k in -extent..=extent {
                    let squared = i * i + j * j + k * k;
                    if (i + j + k) % 2 == 0
                        && squared > 0
                        && squared <= 2 * shell_cutoff_squared as i32
                    {
                        counts[(squared / 2) as usize] += 1;
                    }
                }
            }
        }
        if counts[shell_cutoff_squared as usize] == 0 {
            return Err(EosError::InvalidParameter {
                name: "shell_cutoff_squared",
                reason: "must be an occupied fcc shell",
            });
        }
        let shells = counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(m, n)| (m as f64, f64::from(*n)))
            .collect();
        let (mut lo, mut hi): (f64, f64) = (0.8, 1.0);
        for _ in 0..64 {
            let mid = (lo + hi) / 2.0;
            if 6.0 * (14.19 * (1.0 - mid)).exp() - 14.19 / mid.powi(6) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(Self {
            shell_cutoff_squared,
            sigma_angstrom: 3.802 * (lo + hi) / 2.0,
            shells,
        })
    }

    #[allow(clippy::many_single_char_names)] // Symbols follow the source equations.
    fn energy_and_dv(&self, volume: f64, temperature: f64) -> EosResult<(f64, f64)> {
        let v = positive_state(volume, "volume")? * 10.0;
        let t = finite_state(temperature, "temperature")?;
        if t < 0.0 {
            return Err(EosError::InvalidState {
                name: "temperature",
                reason: "must be nonnegative",
            });
        }
        let rnn = (2.0_f64.sqrt() * v / NA).cbrt();
        let rc = f64::from(self.shell_cutoff_squared).sqrt() * rnn;
        let (eps, steepness, rmin) = (134.7, 14.19, 3.802_f64);
        let decay = steepness / rmin;
        let (mut pair, mut pair_dv) = (0.0, 0.0);
        for (m, n) in &self.shells {
            let radius = m.sqrt() * rnn;
            let rep = eps * 6.0 / (steepness - 6.0) * (steepness - decay * radius).exp();
            let att = eps * steepness / (steepness - 6.0) * (rmin / radius).powi(6);
            pair += n * (rep - att) / 2.0;
            pair_dv += n * (-decay * radius * rep / 3.0 + 2.0 * att) / (2.0 * v);
        }
        let rep_rc = eps * 6.0 / (steepness - 6.0) * (steepness - decay * rc).exp();
        let b = eps * steepness / (steepness - 6.0) * rmin.powi(6);
        let prefactor = 2.0 * std::f64::consts::PI * NA / v;
        let tail = prefactor
            * (rep_rc * (rc.powi(2) / decay + 2.0 * rc / decay.powi(2) + 2.0 / decay.powi(3))
                - b / (3.0 * rc.powi(3)));
        let tail_dv = -tail / v - prefactor / v * (rep_rc - b / rc.powi(6)) * rc.powi(3) / 3.0;
        pair += tail;
        pair_dv += tail_dv;
        let correction = 3.202e5 * NA / (eps * self.sigma_angstrom.powi(6));
        let mut energy = pair * (1.0 - correction / v);
        let mut derivative = pair_dv * (1.0 - correction / v) + pair * correction / v.powi(2);
        let zpv = 140.0 * (2.34 / 0.683 * (1.0 - (v / 19.7).powf(0.683))).exp();
        energy += zpv;
        derivative -= zpv * 2.34 * (v / 19.7).powf(0.683) / v;
        let cold1 = 0.7204 * (-1.614 * (1.0 - v / VREF)).exp();
        let cold2 = -0.01943 * VREF / v * (-27.64 * (1.0 - v / VREF)).exp();
        energy += cold1 + cold2;
        derivative += cold1 * 1.614 / VREF + cold2 * (-1.0 / v + 27.64 / VREF);
        if t > 0.0 {
            let theta_t = 92.0 + 10.4 * (-0.02503 * t.powi(2) - 0.000_156_8 * t.powi(3)).exp_m1();
            let gamma_d = 2.563 * (v / VREF).powf(0.2874);
            let theta = theta_t * (2.563 / 0.2874 * (1.0 - (v / VREF).powf(0.2874))).exp();
            let x = theta / t;
            // The paper's D3 is 1/3 of the convention used by Peritheos.
            let d3 = debye_function_3(x)? / 3.0;
            let weight = 1.0 - 0.0261 - 0.03784 - 0.04512;
            energy += 3.0 * weight * t * ((-(-x).exp_m1()).ln() - d3);
            derivative -= 9.0 * weight * t * gamma_d * d3 / v;
            for (w, theta0, gamma0) in [
                (0.0261, 77.81, 6.221),
                (0.03784, 550.0, 1.617e-6),
                (0.04512, 45.36, 3.127),
            ] {
                let theta_i = theta0 * (gamma0 * (1.0 - v / VREF)).exp();
                let y = theta_i / t;
                energy += 3.0 * w * t * (-(-y).exp_m1()).ln();
                let occupation = if y > 700.0 { 0.0 } else { 1.0 / y.exp_m1() };
                derivative -= 3.0 * w * theta_i * occupation * gamma0 / VREF;
            }
            let anh = -0.000_447_5 * theta_t * (t / 92.0).powi(4)
                / (1.0 + 2.041e-6 * (t / 92.0).powi(2))
                * (5.75e-7 * (v / VREF - 1.0)).exp();
            energy += anh;
            derivative += anh * 5.75e-7 / VREF;
        }
        Ok((finite_result(R * energy)?, finite_result(R * derivative)?))
    }

    /// Unshifted molar Helmholtz energy (no fluid-reference offsets).
    /// # Errors
    /// Returns an error for invalid states or nonfinite results.
    pub fn molar_helmholtz_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        Ok(self.energy_and_dv(volume, temperature)?.0)
    }

    /// Pressure, including explicit density dependence of the effective potential.
    /// # Errors
    /// Returns an error for invalid states or nonfinite results.
    pub fn pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        Ok(-self.energy_and_dv(volume, temperature)?.1 * 1e-3)
    }

    /// Pressure increment relative to 300 K; this is not the total pressure.
    /// # Errors
    /// Returns errors from pressure evaluation.
    pub fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        Ok(self.pressure(volume, temperature)? - self.pressure(volume, 300.0)?)
    }

    /// Isothermal bulk modulus from the pressure derivative.
    /// # Errors
    /// Returns an error for invalid states or derivative step.
    pub fn bulk_modulus(
        &self,
        volume: f64,
        temperature: f64,
        relative_step: f64,
    ) -> EosResult<f64> {
        positive_state(volume, "volume")?;
        if !relative_step.is_finite() || relative_step <= 0.0 || relative_step >= 1.0 {
            return Err(EosError::InvalidState {
                name: "relative_step",
                reason: "must lie in (0, 1)",
            });
        }
        let h = volume * relative_step;
        finite_result(
            -volume
                * (self.pressure(volume + h, temperature)?
                    - self.pressure(volume - h, temperature)?)
                / (2.0 * h),
        )
    }

    /// Stable root in the numerical 8--26 cm3/mol search interval.
    /// This interval is not a phase-stability claim.
    /// # Errors
    /// Returns an error when the state is invalid or no stable root is present.
    pub fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        finite_state(pressure, "pressure")?;
        let mut lo = 0.8;
        let mut flo = self.pressure(lo, temperature)? - pressure;
        for i in 1..=180 {
            let mut hi = 0.8 + f64::from(i) * 0.01;
            let fhi = self.pressure(hi, temperature)? - pressure;
            if flo >= 0.0 && fhi <= 0.0 {
                for _ in 0..48 {
                    let mid = (lo + hi) / 2.0;
                    if self.pressure(mid, temperature)? > pressure {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let result = (lo + hi) / 2.0;
                if self.bulk_modulus(result, temperature, 1e-5)? > 0.0 {
                    return Ok(result);
                }
            }
            lo = hi;
            flo = fhi;
        }
        Err(EosError::OutsideInvertibleRange)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_quadrature_and_thermodynamic_identity() {
        let model = Maltby2024Published::new(64).unwrap();
        // SciPy direct quadrature of source Helmholtz equations, not Table 8.
        assert!((model.pressure(2.397, 70.0).unwrap() + 0.00422597965).abs() < 2e-9);
        assert!((model.pressure(1.8, 300.0).unwrap() - 2.21438377236).abs() < 8e-9);
        assert!((model.pressure(1.2, 0.0).unwrap() - 16.50973820028).abs() < 8e-9);
        assert!((model.molar_helmholtz_energy(2.397, 70.0).unwrap() + 9059.06647911).abs() < 1e-6);
        for v in [1.2, 1.8, 2.397] {
            for t in [0.0, 5.0, 70.0, 300.0] {
                let h = v * 1e-5;
                let derivative = -(model.molar_helmholtz_energy(v + h, t).unwrap()
                    - model.molar_helmholtz_energy(v - h, t).unwrap())
                    / (2.0 * h)
                    * 1e-4;
                assert!((model.pressure(v, t).unwrap() - derivative).abs() < 1e-7);
                let p = model.pressure(v, t).unwrap();
                assert!((model.volume(p, t).unwrap() - v).abs() < 1e-10);
            }
        }
        assert_eq!(model.thermal_pressure_increment(2.0, 300.0).unwrap(), 0.0);
    }

    #[test]
    fn validates_inputs_and_retains_source_discrepancy() {
        let model = Maltby2024Published::new(64).unwrap();
        assert!(Maltby2024Published::new(0).is_err());
        assert!(Maltby2024Published::new(30).is_err());
        assert!(model.pressure(0.0, 70.0).is_err());
        assert!(model.pressure(2.0, -1.0).is_err());
        assert!(model.pressure(2.0, f64::NAN).is_err());
        assert!(model.bulk_modulus(2.0, 70.0, 0.0).is_err());
        // Explicitly guard against silently tuning published coefficients.
        assert!((model.volume(0.001, 70.0).unwrap() * 10.0 - 23.97).abs() > 0.07);
    }
}
