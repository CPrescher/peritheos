//! Bounded pressure tables for electronic and derived residual corrections.

use super::{debye_function_3, MieGruneisenDebye, GAS_CONSTANT};
use crate::validation::{finite_result, finite_state, positive_parameter, positive_state};
use crate::{EosError, EosResult, IsothermalEos, ThermalEos};

/// Finite piecewise-linear pressure table with no temperature extrapolation.
#[derive(Clone, Debug, PartialEq)]
pub struct PressureTable {
    temperatures: Vec<f64>,
    pressures: Vec<f64>,
}

impl PressureTable {
    /// Construct a table, optionally requiring nondecreasing pressures.
    ///
    /// # Errors
    /// Rejects mismatched, nonfinite, or unordered coordinates.
    pub fn new(temperatures: Vec<f64>, pressures: Vec<f64>, monotone: bool) -> EosResult<Self> {
        if temperatures.len() < 2
            || temperatures.len() != pressures.len()
            || temperatures
                .iter()
                .chain(&pressures)
                .any(|x| !x.is_finite())
            || temperatures[0] < 0.0
            || temperatures.windows(2).any(|x| x[0] >= x[1])
            || (monotone && pressures.windows(2).any(|x| x[0] > x[1]))
        {
            return Err(EosError::InvalidParameter {
                name: "pressure table",
                reason: "requires matching finite arrays and increasing nonnegative temperatures",
            });
        }
        Ok(Self {
            temperatures,
            pressures,
        })
    }

    fn bounds(&self) -> [f64; 2] {
        [
            self.temperatures[0],
            self.temperatures[self.temperatures.len() - 1],
        ]
    }

    fn pressure(&self, temperature: f64) -> EosResult<f64> {
        let [lo, hi] = self.bounds();
        if !temperature.is_finite() || temperature < lo || temperature > hi {
            return Err(EosError::InvalidState {
                name: "temperature",
                reason: "outside pressure table",
            });
        }
        let i = self
            .temperatures
            .partition_point(|t| *t < temperature)
            .clamp(1, self.temperatures.len() - 1);
        let fraction = (temperature - self.temperatures[i - 1])
            / (self.temperatures[i] - self.temperatures[i - 1]);
        finite_result(
            self.pressures[i - 1] + fraction * (self.pressures[i] - self.pressures[i - 1]),
        )
    }
}

/// Reference-temperature Debye pressure plus a bounded electronic table.
#[derive(Clone, Debug, PartialEq)]
pub struct DebyeTabulatedThermalPressure<R> {
    /// Debye model, including its reference isotherm.
    pub debye: MieGruneisenDebye<R>,
    electronic: PressureTable,
}

impl<R: IsothermalEos> DebyeTabulatedThermalPressure<R> {
    /// Compose an integrated Debye model and electronic pressure table.
    ///
    /// # Errors
    /// Rejects a reference temperature outside the table.
    pub fn new(debye: MieGruneisenDebye<R>, electronic: PressureTable) -> EosResult<Self> {
        electronic.pressure(debye.tr)?;
        Ok(Self { debye, electronic })
    }
}

impl<R: IsothermalEos> ThermalEos for DebyeTabulatedThermalPressure<R> {
    type Reference = R;
    fn reference_eos(&self) -> &R {
        &self.debye.rt_eos
    }
    fn reference_temperature(&self) -> f64 {
        self.debye.tr
    }
    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.debye.thermal_pressure(volume, temperature)?
                + self.electronic.pressure(temperature)?
                - self.electronic.pressure(self.debye.tr)?,
        )
    }
}

/// Absolute phonon pressure with separate cold/ambient volumes and bounded tables.
/// This pressure reconstruction supplies no electronic caloric potential.
#[derive(Clone, Debug, PartialEq)]
pub struct AsymptoticDebyeTabulatedPressure<R> {
    /// Cold 0 K reference curve.
    pub rt_eos: R,
    tr: f64,
    theta0: f64,
    gamma0: f64,
    a: f64,
    b: f64,
    n: f64,
    cold_volume_ratio: f64,
    electronic: PressureTable,
    residual: PressureTable,
    volume_ratio_range: [f64; 2],
    temperature_range: [f64; 2],
}

impl<R: IsothermalEos> AsymptoticDebyeTabulatedPressure<R> {
    /// Construct the bounded absolute pressure surface, in molar J/bar/mol units.
    ///
    /// # Errors
    /// Rejects invalid coefficients, bounds, or incomplete temperature tables.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rt_eos: R,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        a: f64,
        b: f64,
        n: f64,
        cold_volume_ratio: f64,
        electronic: PressureTable,
        residual: PressureTable,
        volume_ratio_range: [f64; 2],
        temperature_range: [f64; 2],
    ) -> EosResult<Self> {
        for (bounds, positive) in [(volume_ratio_range, true), (temperature_range, false)] {
            if bounds.iter().any(|x| !x.is_finite())
                || bounds[0] < 0.0
                || (positive && bounds[0] <= 0.0)
                || bounds[0] >= bounds[1]
            {
                return Err(EosError::InvalidParameter {
                    name: "reconstruction bounds",
                    reason: "requires ordered finite nonnegative bounds",
                });
            }
        }
        if !a.is_finite()
            || !(0.0..=1.0).contains(&a)
            || tr < temperature_range[0]
            || tr > temperature_range[1]
        {
            return Err(EosError::InvalidParameter {
                name: "a or Tr",
                reason: "outside reconstruction bounds",
            });
        }
        for t in temperature_range {
            electronic.pressure(t)?;
            residual.pressure(t)?;
        }
        Ok(Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            theta0: positive_parameter(theta0, "theta0")?,
            gamma0: positive_parameter(gamma0, "gamma0")?,
            a,
            b: positive_parameter(b, "b")?,
            n: positive_parameter(n, "n")?,
            cold_volume_ratio: positive_parameter(cold_volume_ratio, "cold_volume_ratio")?,
            electronic,
            residual,
            volume_ratio_range,
            temperature_range,
        })
    }

    fn ambient_volume(&self) -> f64 {
        self.rt_eos.reference_volume() / self.cold_volume_ratio
    }

    fn volume_root(&self, pressure: f64, temperature: f64, f_dac: f64) -> EosResult<f64> {
        let target = finite_state(pressure, "pressure")?;
        bounded_root(
            |v| {
                Ok(self.pressure(v, temperature)?
                    - f_dac * self.thermal_pressure_increment(v, temperature)?
                    - target)
            },
            self.volume_ratio_range.map(|x| x * self.ambient_volume()),
        )
    }
}

impl<R: IsothermalEos> ThermalEos for AsymptoticDebyeTabulatedPressure<R> {
    type Reference = R;
    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }
    fn reference_temperature(&self) -> f64 {
        self.tr
    }
    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let ratio = positive_state(volume, "volume")? / self.ambient_volume();
        let [lo, hi] = self.volume_ratio_range;
        if ratio < lo - 16.0 * f64::EPSILON
            || ratio > hi + 16.0 * f64::EPSILON
            || !temperature.is_finite()
            || temperature < self.temperature_range[0]
            || temperature > self.temperature_range[1]
        {
            return Err(EosError::InvalidState {
                name: "volume or temperature",
                reason: "outside reconstruction domain",
            });
        }
        let gamma = self.gamma0 * (1.0 + self.a * (ratio.powf(self.b) - 1.0));
        let theta = self.theta0
            * ratio.powf(-(1.0 - self.a) * self.gamma0)
            * (-(gamma - self.gamma0) / self.b).exp();
        let energy = if temperature == 0.0 {
            0.0
        } else {
            3.0 * self.n * GAS_CONSTANT * temperature * debye_function_3(theta / temperature)?
        };
        finite_result(
            gamma * energy / volume / 1.0e4
                + self.electronic.pressure(temperature)?
                + self.residual.pressure(temperature)?,
        )
    }
    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        Ok(self.thermal_pressure(volume, temperature)? - self.thermal_pressure(volume, self.tr)?)
    }
    fn volume(&self, pressure: f64, temperature: f64) -> EosResult<f64> {
        self.volume_root(pressure, temperature, 0.0)
    }
    fn volume_with_dac_confinement(
        &self,
        pressure: f64,
        temperature: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        if !f_dac.is_finite() || !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        self.volume_root(pressure, temperature, f_dac)
    }
}

fn bounded_root<F: Fn(f64) -> EosResult<f64>>(
    function: F,
    [mut lo, mut hi]: [f64; 2],
) -> EosResult<f64> {
    let mut left = function(lo)?;
    let right = function(hi)?;
    if left.abs() <= 1e-10 {
        return Ok(lo);
    }
    if right.abs() <= 1e-10 {
        return Ok(hi);
    }
    if left.signum() == right.signum() {
        return Err(EosError::OutsideInvertibleRange);
    }
    for _ in 0..256 {
        let mid = (lo + hi) / 2.0;
        let value = function(mid)?;
        if value.abs() <= 1e-10 || (hi - lo).abs() <= 1e-13 * mid.abs().max(1.0) {
            return Ok(mid);
        }
        if value.signum() == left.signum() {
            lo = mid;
            left = value;
        } else {
            hi = mid;
        }
    }
    Err(EosError::OutsideInvertibleRange)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isothermal::BM3;

    #[test]
    fn table_validation_and_interpolation_are_bounded() {
        for (t, p, monotone) in [
            (vec![0.0], vec![0.0], true),
            (vec![0.0, 300.0], vec![0.0], true),
            (vec![300.0, 300.0], vec![0.0, 1.0], true),
            (vec![-1.0, 300.0], vec![0.0, 1.0], true),
            (vec![0.0, f64::NAN], vec![0.0, 1.0], true),
            (vec![0.0, 300.0], vec![1.0, 0.0], true),
        ] {
            assert!(PressureTable::new(t, p, monotone).is_err());
        }
        let table =
            PressureTable::new(vec![0.0, 300.0, 1000.0], vec![0.0, 0.3, 1.0], true).unwrap();
        for (t, p) in [
            (0.0, 0.0),
            (150.0, 0.15),
            (300.0, 0.3),
            (650.0, 0.65),
            (1000.0, 1.0),
        ] {
            assert!((table.pressure(t).unwrap() - p).abs() < 1e-12);
        }
        assert!(table.pressure(-1.0).is_err());
        assert!(table.pressure(1001.0).is_err());
        assert!(table.pressure(f64::NAN).is_err());
        let reference = BM3::new(1.0, 160.0, 4.0).unwrap();
        let debye = MieGruneisenDebye::new(reference, 300.0, 800.0, 1.5, 1.0, 2.0).unwrap();
        let model = DebyeTabulatedThermalPressure::new(debye, table).unwrap();
        assert!(model.thermal_pressure(0.8, 300.0).unwrap().abs() < 1e-12);
        assert!(model.thermal_pressure(0.8, 1001.0).is_err());
    }
}
