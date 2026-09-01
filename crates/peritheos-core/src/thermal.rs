//! Built-in thermal equations of state and caloric models.

use crate::isothermal::{Holzapfel, ModifiedTait};
use crate::quadrature::integrate;
use crate::root::solve_temperature_function;
use crate::validation::{
    finite_parameter, finite_result, finite_state, positive_parameter, positive_state,
};
use crate::{CaloricEos, EosError, EosResult, IsothermalEos, ThermalEos};

/// Universal molar gas constant frozen to the `SciPy` 1.18.1 baseline value.
pub const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Third-order Debye function with stable small- and large-argument limits.
///
/// # Errors
///
/// Returns an error unless `argument` is positive and finite or if numerical
/// quadrature fails to converge.
pub fn debye_function_3(argument: f64) -> EosResult<f64> {
    let argument = positive_state(argument, "Debye-function argument")?;
    if argument < 1.0e-3 {
        return Ok(1.0 - 3.0 * argument / 8.0 + argument.powi(2) / 20.0 - argument.powi(4) / 1680.0);
    }
    if argument > 150.0 {
        return Ok(std::f64::consts::PI.powi(4) / (5.0 * argument.powi(3)));
    }
    let integral = integrate(
        |value| {
            if value.abs() < f64::EPSILON {
                Ok(0.0)
            } else {
                finite_result(value.powi(3) / value.exp_m1())
            }
        },
        0.0,
        argument,
    )?;
    finite_result(3.0 * integral / argument.powi(3))
}

/// Shared representation underlying the public Debye and Einstein aliases.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MieGruneisen<R, const DEBYE: bool> {
    /// Reference isothermal EOS.
    pub rt_eos: R,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Reference characteristic temperature in kelvin.
    pub theta0: f64,
    /// Reference Gruneisen parameter.
    pub gamma0: f64,
    /// Volume exponent of the Gruneisen parameter.
    pub q: f64,
    /// Number of atoms in the formula unit.
    pub n: f64,
}

/// Mie--Gruneisen--Debye thermal equation of state.
pub type MieGruneisenDebye<R> = MieGruneisen<R, true>;

/// Mie--Gruneisen--Einstein thermal equation of state.
pub type MieGruneisenEinstein<R> = MieGruneisen<R, false>;

fn new_mie_gruneisen<R, const DEBYE: bool>(
    rt_eos: R,
    tr: f64,
    theta0: f64,
    gamma0: f64,
    q: f64,
    n: f64,
) -> EosResult<MieGruneisen<R, DEBYE>> {
    Ok(MieGruneisen {
        rt_eos,
        tr: positive_parameter(tr, "Tr")?,
        theta0: positive_parameter(theta0, "theta0")?,
        gamma0: finite_parameter(gamma0, "gamma0")?,
        q: finite_parameter(q, "q")?,
        n: positive_parameter(n, "n")?,
    })
}

impl<R> MieGruneisen<R, true>
where
    R: IsothermalEos,
{
    /// Construct a Mie--Gruneisen--Debye model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite thermal parameters.
    pub fn new(rt_eos: R, tr: f64, theta0: f64, gamma0: f64, q: f64, n: f64) -> EosResult<Self> {
        new_mie_gruneisen(rt_eos, tr, theta0, gamma0, q, n)
    }
}

impl<R> MieGruneisen<R, false>
where
    R: IsothermalEos,
{
    /// Construct a Mie--Gruneisen--Einstein model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite thermal parameters.
    pub fn new(rt_eos: R, tr: f64, theta0: f64, gamma0: f64, q: f64, n: f64) -> EosResult<Self> {
        new_mie_gruneisen(rt_eos, tr, theta0, gamma0, q, n)
    }
}

impl<R, const DEBYE: bool> MieGruneisen<R, DEBYE>
where
    R: IsothermalEos,
{
    /// Volume-dependent Gruneisen parameter.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a non-finite result.
    pub fn volume_gruneisen_parameter(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        finite_result(self.gamma0 * (self.q * (volume / self.rt_eos.reference_volume()).ln()).exp())
    }

    /// Volume-dependent characteristic lattice temperature.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume or a non-finite result.
    #[allow(clippy::float_cmp)]
    pub fn characteristic_temperature(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let logarithmic_volume = (volume / self.rt_eos.reference_volume()).ln();
        let exponent = if self.q == 0.0 {
            -self.gamma0 * logarithmic_volume
        } else {
            -self.gamma0 * (self.q * logarithmic_volume).exp_m1() / self.q
        };
        finite_result(self.theta0 * exponent.exp())
    }

    /// Vibrational thermal energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or quadrature failure.
    pub fn thermal_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let theta = self.characteristic_temperature(volume)?;
        if DEBYE {
            finite_result(
                3.0 * self.n * GAS_CONSTANT * temperature * debye_function_3(theta / temperature)?,
            )
        } else {
            let ratio = theta / temperature;
            let decay = (-ratio).exp();
            finite_result(3.0 * self.n * GAS_CONSTANT * theta * decay / (-(-ratio).exp_m1()))
        }
    }

    /// Vibrational entropy in J mol^-1 K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or quadrature failure.
    pub fn thermal_entropy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let ratio = self.characteristic_temperature(volume)? / temperature;
        let log_term = (-(-ratio).exp_m1()).ln();
        if DEBYE {
            finite_result(self.n * GAS_CONSTANT * (4.0 * debye_function_3(ratio)? - 3.0 * log_term))
        } else {
            let occupation = (-ratio).exp() / (-(-ratio).exp_m1());
            finite_result(3.0 * self.n * GAS_CONSTANT * (ratio * occupation - log_term))
        }
    }

    /// Unreferenced vibrational pressure in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or model evaluation.
    pub fn vibrational_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.volume_gruneisen_parameter(volume)? * self.thermal_energy(volume, temperature)?
                / volume
                / 1.0e4,
        )
    }

    /// Vibrational Helmholtz free energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or model evaluation.
    pub fn thermal_helmholtz_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.thermal_energy(volume, temperature)?
                - temperature * self.thermal_entropy(volume, temperature)?,
        )
    }

    /// Vibrational enthalpy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or model evaluation.
    pub fn thermal_enthalpy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.thermal_energy(volume, temperature)?
                + self.vibrational_pressure(volume, temperature)? * volume * 1.0e4,
        )
    }

    /// Vibrational Gibbs free energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid state or model evaluation.
    pub fn thermal_gibbs_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.thermal_helmholtz_free_energy(volume, temperature)?
                + self.vibrational_pressure(volume, temperature)? * volume * 1.0e4,
        )
    }
}

impl<R, const DEBYE: bool> ThermalEos for MieGruneisen<R, DEBYE>
where
    R: IsothermalEos,
{
    type Reference = R;

    fn reference_eos(&self) -> &Self::Reference {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let energy_difference =
            self.thermal_energy(volume, temperature)? - self.thermal_energy(volume, self.tr)?;
        finite_result(self.volume_gruneisen_parameter(volume)? * energy_difference / volume / 1.0e4)
    }
}

impl<R, const DEBYE: bool> CaloricEos for MieGruneisen<R, DEBYE>
where
    R: IsothermalEos,
{
    fn molar_heat_capacity_v(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let temperature = positive_state(temperature, "temperature")?;
        let step = 1.0e-5 * temperature;
        finite_result(
            (self.thermal_energy(volume, temperature + step)?
                - self.thermal_energy(volume, temperature - step)?)
                / (2.0 * step),
        )
    }

    fn gruneisen_parameter(&self, volume: f64, _temperature: f64) -> EosResult<f64> {
        self.volume_gruneisen_parameter(volume)
    }
}

/// Holland--Powell thermal modified Tait equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalModifiedTait {
    /// Modified Tait reference isotherm.
    pub rt_eos: ModifiedTait,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Einstein characteristic temperature in kelvin.
    pub theta: f64,
    /// Reference volumetric thermal expansivity in `K^-1`.
    pub alpha0: f64,
    /// Number of atoms in the formula unit.
    pub n: f64,
    cv0: f64,
    pressure_factor: f64,
}

/// Compatibility alias for [`ThermalModifiedTait`].
pub type HollandPowell2011 = ThermalModifiedTait;

impl ThermalModifiedTait {
    /// Construct a Holland--Powell thermal modified Tait model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite thermal parameters.
    pub fn new(rt_eos: ModifiedTait, tr: f64, theta: f64, alpha0: f64, n: f64) -> EosResult<Self> {
        let tr = positive_parameter(tr, "Tr")?;
        let theta = positive_parameter(theta, "theta")?;
        let alpha0 = finite_parameter(alpha0, "alpha0")?;
        let n = positive_parameter(n, "n")?;
        let cv0 = einstein_heat_capacity(theta, n, tr)?;
        let pressure_factor = alpha0 * rt_eos.k0 / cv0;
        finite_result(pressure_factor)?;
        Ok(Self {
            rt_eos,
            tr,
            theta,
            alpha0,
            n,
            cv0,
            pressure_factor,
        })
    }

    /// Reference constant-volume heat capacity used by the pressure model.
    #[must_use]
    pub fn reference_heat_capacity_v(&self) -> f64 {
        self.cv0
    }

    fn einstein_energy(&self, temperature: f64) -> EosResult<f64> {
        einstein_energy(self.theta, self.n, temperature)
    }
}

impl ThermalEos for ThermalModifiedTait {
    type Reference = ModifiedTait;

    fn reference_eos(&self) -> &Self::Reference {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        finite_result(
            self.pressure_factor
                * (self.einstein_energy(temperature)? - self.einstein_energy(self.tr)?),
        )
    }
}

impl CaloricEos for ThermalModifiedTait {
    fn molar_heat_capacity_v(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        positive_state(volume, "volume")?;
        einstein_heat_capacity(self.theta, self.n, temperature)
    }

    fn gruneisen_parameter(&self, volume: f64, _temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        finite_result(volume * self.pressure_factor * 1.0e4)
    }
}

fn einstein_energy(theta: f64, n: f64, temperature: f64) -> EosResult<f64> {
    let temperature = positive_state(temperature, "temperature")?;
    let ratio = theta / temperature;
    finite_result(3.0 * n * GAS_CONSTANT * theta * (-ratio).exp() / (-(-ratio).exp_m1()))
}

fn einstein_heat_capacity(theta: f64, n: f64, temperature: f64) -> EosResult<f64> {
    let temperature = positive_state(temperature, "temperature")?;
    let ratio = theta / temperature;
    let decay = (-ratio).exp();
    finite_result(3.0 * n * GAS_CONSTANT * ratio * ratio * decay / (1.0 - decay).powi(2))
}

/// Parameters of the Sokolova et al. (2016) thermal-pressure model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SokolovaParameters {
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// First Einstein characteristic temperature in kelvin.
    pub qe1o: f64,
    /// First Einstein multiplicity.
    pub me1: f64,
    /// Second Einstein characteristic temperature in kelvin.
    pub qe2o: f64,
    /// Second Einstein multiplicity.
    pub me2: f64,
    /// Additive Gruneisen normalization.
    pub delta: f64,
    /// Generalized Gruneisen parameter.
    pub t: f64,
    /// Intrinsic anharmonicity parameter in `10^-6 K^-1`.
    pub a_0: f64,
    /// Anharmonic volume exponent.
    pub m: f64,
    /// Electronic volume exponent.
    pub g: f64,
    /// Free-electron parameter in `10^-6 K^-1`.
    pub e_0: f64,
    /// Volume-dependent correction to `t`.
    pub beta: f64,
    /// First Bose-mode reference characteristic temperature.
    pub qbo: f64,
    /// First Bose-mode dispersion.
    pub d: f64,
    /// First Bose-mode multiplicity.
    pub mb: f64,
    /// Second Bose-mode reference characteristic temperature.
    pub qb1o: f64,
    /// Second Bose-mode dispersion.
    pub d1: f64,
    /// Second Bose-mode multiplicity.
    pub mb1: f64,
}

impl SokolovaParameters {
    /// Construct the reduced model used by historical Peritheos releases.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn reduced(
        tr: f64,
        qe1o: f64,
        me1: f64,
        qe2o: f64,
        me2: f64,
        delta: f64,
        t: f64,
        a_0: f64,
        m: f64,
        g: f64,
        e_0: f64,
    ) -> Self {
        Self {
            tr,
            qe1o,
            me1,
            qe2o,
            me2,
            delta,
            t,
            a_0,
            m,
            g,
            e_0,
            beta: 0.0,
            qbo: 1.0,
            d: 1.0,
            mb: 0.0,
            qb1o: 1.0,
            d1: 1.0,
            mb1: 0.0,
        }
    }
}

/// Sokolova et al. (2016) thermal-pressure equation with a Holzapfel reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sokolova2016 {
    /// Holzapfel reference isotherm.
    pub rt_eos: Holzapfel,
    /// Validated thermal parameters.
    pub parameters: SokolovaParameters,
}

#[derive(Clone, Copy, Debug)]
struct SokolovaVolumeTerms {
    volume: f64,
    gamma: f64,
    qb: f64,
    qb1: f64,
    qe1: f64,
    qe2: f64,
    reference_oscillator_pressure: f64,
    squared_temperature_coefficient: f64,
}

impl Sokolova2016 {
    /// Construct the Sokolova thermal-pressure model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid characteristic temperatures,
    /// dispersions, multiplicities, or other non-finite parameters.
    pub fn new(rt_eos: Holzapfel, parameters: SokolovaParameters) -> EosResult<Self> {
        let parameters = SokolovaParameters {
            tr: positive_parameter(parameters.tr, "Tr")?,
            qe1o: positive_parameter(parameters.qe1o, "QE1o")?,
            me1: nonnegative_parameter(parameters.me1, "mE1")?,
            qe2o: positive_parameter(parameters.qe2o, "QE2o")?,
            me2: nonnegative_parameter(parameters.me2, "mE2")?,
            delta: finite_parameter(parameters.delta, "delta")?,
            t: finite_parameter(parameters.t, "t")?,
            a_0: finite_parameter(parameters.a_0, "a_0")?,
            m: finite_parameter(parameters.m, "m")?,
            g: finite_parameter(parameters.g, "g")?,
            e_0: finite_parameter(parameters.e_0, "e_0")?,
            beta: finite_parameter(parameters.beta, "beta")?,
            qbo: positive_parameter(parameters.qbo, "QBo")?,
            d: positive_parameter(parameters.d, "d")?,
            mb: nonnegative_parameter(parameters.mb, "mb")?,
            qb1o: positive_parameter(parameters.qb1o, "QB1o")?,
            d1: positive_parameter(parameters.d1, "d1")?,
            mb1: nonnegative_parameter(parameters.mb1, "mb1")?,
        };
        Ok(Self { rt_eos, parameters })
    }

    /// Integrate the model Gruneisen parameter from `x` to the reference ratio.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid volume ratio, model evaluation, or
    /// quadrature failure.
    pub fn gruneisen_integral(&self, x: f64) -> EosResult<f64> {
        let x = positive_state(x, "volume ratio")?;
        integrate(
            |ratio| {
                let volume = ratio * self.rt_eos.v0;
                let pressure = self.rt_eos.pressure(volume)?;
                let bulk_modulus = self.rt_eos.bulk_modulus(volume)?;
                let derivative = self.rt_eos.bulk_modulus_derivative(volume, 1.0e-6)?;
                sokolova_gruneisen_integrand(
                    ratio,
                    pressure,
                    bulk_modulus,
                    derivative,
                    self.parameters.delta,
                    self.parameters.t,
                    self.parameters.beta,
                )
            },
            x,
            1.0,
        )
    }

    fn volume_terms(&self, volume: f64) -> EosResult<SokolovaVolumeTerms> {
        let volume = positive_state(volume, "volume")?;
        let x = volume / self.rt_eos.v0;
        let pressure = self.rt_eos.pressure(volume)?;
        let bulk_modulus = self.rt_eos.bulk_modulus(volume)?;
        let derivative = self.rt_eos.bulk_modulus_derivative(volume, 1.0e-6)?;
        let generalized_t = self.parameters.t - self.parameters.beta * x.cbrt();
        let gamma = (-3.0 * bulk_modulus
            + 2.0 * pressure * generalized_t
            + 9.0 * bulk_modulus * derivative
            - 6.0 * generalized_t * bulk_modulus)
            / (6.0 * (3.0 * bulk_modulus - 2.0 * pressure * generalized_t))
            + self.parameters.delta;
        let scale = self.gruneisen_integral(x)?.exp();
        let qb = self.parameters.qbo * scale;
        let qb1 = self.parameters.qb1o * scale;
        let qe1 = self.parameters.qe1o * scale;
        let qe2 = self.parameters.qe2o * scale;
        let reference_oscillator_pressure = self.parameters.mb
            * GAS_CONSTANT
            * bose_mode_energy(qb, self.parameters.tr, self.parameters.d)?
            * gamma
            / volume
            + self.parameters.mb1
                * GAS_CONSTANT
                * bose_mode_energy(qb1, self.parameters.tr, self.parameters.d1)?
                * gamma
                / volume
            + self.parameters.me1
                * GAS_CONSTANT
                * sokolova_einstein_energy(qe1, self.parameters.tr)?
                * gamma
                / volume
            + self.parameters.me2
                * GAS_CONSTANT
                * sokolova_einstein_energy(qe2, self.parameters.tr)?
                * gamma
                / volume;
        let squared_temperature_coefficient =
            1.5 * self.rt_eos.n * GAS_CONSTANT / 1_000_000.0 / volume
                * (self.parameters.a_0 * x.powf(self.parameters.m) * self.parameters.m
                    + self.parameters.e_0 * x.powf(self.parameters.g) * self.parameters.g);
        finite_result(gamma)?;
        finite_result(scale)?;
        finite_result(reference_oscillator_pressure)?;
        finite_result(squared_temperature_coefficient)?;
        Ok(SokolovaVolumeTerms {
            volume,
            gamma,
            qb,
            qb1,
            qe1,
            qe2,
            reference_oscillator_pressure,
            squared_temperature_coefficient,
        })
    }

    #[allow(clippy::similar_names)]
    fn thermal_pressure_from_terms(
        &self,
        terms: SokolovaVolumeTerms,
        temperature: f64,
    ) -> EosResult<f64> {
        let temperature = positive_state(temperature, "temperature")?;
        let pressure_b = self.parameters.mb
            * GAS_CONSTANT
            * bose_mode_energy(terms.qb, temperature, self.parameters.d)?
            * terms.gamma
            / terms.volume;
        let pressure_b1 = self.parameters.mb1
            * GAS_CONSTANT
            * bose_mode_energy(terms.qb1, temperature, self.parameters.d1)?
            * terms.gamma
            / terms.volume;
        let pressure_e1 = self.parameters.me1
            * GAS_CONSTANT
            * sokolova_einstein_energy(terms.qe1, temperature)?
            * terms.gamma
            / terms.volume;
        let pressure_e2 = self.parameters.me2
            * GAS_CONSTANT
            * sokolova_einstein_energy(terms.qe2, temperature)?
            * terms.gamma
            / terms.volume;
        let pressure_bar = pressure_b + pressure_b1 + pressure_e1 + pressure_e2
            - terms.reference_oscillator_pressure
            + terms.squared_temperature_coefficient
                * (temperature * temperature - self.parameters.tr * self.parameters.tr);
        finite_result(pressure_bar / 1.0e4)
    }
}

impl ThermalEos for Sokolova2016 {
    type Reference = Holzapfel;

    fn reference_eos(&self) -> &Self::Reference {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.parameters.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        self.thermal_pressure_from_terms(self.volume_terms(volume)?, temperature)
    }

    fn temperature(&self, pressure: f64, volume: f64) -> EosResult<f64> {
        let pressure = finite_state(pressure, "pressure")?;
        let volume = positive_state(volume, "volume")?;
        let target = pressure - self.rt_eos.pressure(volume)?;
        let terms = self.volume_terms(volume)?;
        solve_temperature_function(
            |temperature| self.thermal_pressure_from_terms(terms, temperature),
            target,
            self.parameters.tr,
        )
    }
}

fn nonnegative_parameter(value: f64, name: &'static str) -> EosResult<f64> {
    let value = finite_parameter(value, name)?;
    if value < 0.0 {
        Err(EosError::InvalidParameter {
            name,
            reason: "must not be negative",
        })
    } else {
        Ok(value)
    }
}

fn sokolova_gruneisen_integrand(
    x: f64,
    pressure: f64,
    bulk_modulus: f64,
    bulk_modulus_derivative: f64,
    delta: f64,
    t: f64,
    beta: f64,
) -> EosResult<f64> {
    let generalized_t = t - beta * x.cbrt();
    finite_result(
        (delta
            + (-3.0 * bulk_modulus
                + 2.0 * pressure * generalized_t
                + 9.0 * bulk_modulus * bulk_modulus_derivative
                - 6.0 * generalized_t * bulk_modulus)
                / (6.0 * (3.0 * bulk_modulus - 2.0 * pressure * generalized_t)))
            / x,
    )
}

fn sokolova_einstein_energy(theta: f64, temperature: f64) -> EosResult<f64> {
    let temperature = positive_state(temperature, "temperature")?;
    let ratio = theta / temperature;
    let decay = (-ratio).exp();
    finite_result(theta / 2.0 + theta * decay / (-(-ratio).exp_m1()))
}

fn bose_mode_energy(theta: f64, temperature: f64, dispersion: f64) -> EosResult<f64> {
    let temperature = positive_state(temperature, "temperature")?;
    let exponent = dispersion * (theta / (temperature * dispersion)).ln_1p();
    let decay = (-exponent).exp();
    let occupation = decay / (-(-exponent).exp_m1());
    let zero_point = theta * (dispersion - 1.0) / (2.0 * dispersion);
    let thermal_part =
        temperature * theta * dispersion * occupation / (temperature * dispersion + theta);
    finite_result(zero_point + thermal_part)
}
