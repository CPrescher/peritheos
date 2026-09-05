//! Built-in thermal equations of state and caloric models.

use crate::isothermal::{
    Holzapfel, ModifiedTait, Murnaghan, NaturalStrain2, NaturalStrain3, NaturalStrain4, Vinet, BM2,
    BM3, BM4,
};
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

/// Supported volume laws for the Debye characteristic temperature.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DebyeTemperatureLaw {
    /// Integrate `gamma = -d ln(theta) / d ln(V)` for the power-law gamma.
    #[default]
    IntegratedGruneisen,
    /// Apply `theta = theta0 (V/V0)^(-gamma(V))` directly.
    VariableExponent,
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
    /// Debye-temperature convention; Einstein models use the integrated law.
    pub debye_temperature_law: DebyeTemperatureLaw,
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
    debye_temperature_law: DebyeTemperatureLaw,
) -> EosResult<MieGruneisen<R, DEBYE>> {
    Ok(MieGruneisen {
        rt_eos,
        tr: positive_parameter(tr, "Tr")?,
        theta0: positive_parameter(theta0, "theta0")?,
        gamma0: finite_parameter(gamma0, "gamma0")?,
        q: finite_parameter(q, "q")?,
        n: positive_parameter(n, "n")?,
        debye_temperature_law,
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
        new_mie_gruneisen(
            rt_eos,
            tr,
            theta0,
            gamma0,
            q,
            n,
            DebyeTemperatureLaw::IntegratedGruneisen,
        )
    }

    /// Construct a model with an explicit Debye-temperature convention.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite thermal parameters.
    pub fn new_with_temperature_law(
        rt_eos: R,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        q: f64,
        n: f64,
        debye_temperature_law: DebyeTemperatureLaw,
    ) -> EosResult<Self> {
        new_mie_gruneisen(rt_eos, tr, theta0, gamma0, q, n, debye_temperature_law)
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
        new_mie_gruneisen(
            rt_eos,
            tr,
            theta0,
            gamma0,
            q,
            n,
            DebyeTemperatureLaw::IntegratedGruneisen,
        )
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
        let exponent =
            if DEBYE && self.debye_temperature_law == DebyeTemperatureLaw::VariableExponent {
                -self.volume_gruneisen_parameter(volume)? * logarithmic_volume
            } else if self.q == 0.0 {
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

/// An isothermal EOS that can be reconstructed at a new `V0` and `K0`.
pub trait ReferenceStateEos: IsothermalEos + Copy {
    /// Reference bulk modulus in `GPa`.
    fn reference_bulk_modulus(&self) -> f64;

    /// Reconstruct the same equation family with new reference values.
    ///
    /// # Errors
    ///
    /// Returns an error when either reference value is invalid.
    fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self>;
}

macro_rules! impl_two_parameter_reference_state {
    ($model:ty, $constructor:path) => {
        impl ReferenceStateEos for $model {
            fn reference_bulk_modulus(&self) -> f64 {
                self.k0
            }

            fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
                $constructor(volume, bulk_modulus)
            }
        }
    };
}

macro_rules! impl_three_parameter_reference_state {
    ($model:ty, $constructor:path) => {
        impl ReferenceStateEos for $model {
            fn reference_bulk_modulus(&self) -> f64 {
                self.k0
            }

            fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
                $constructor(volume, bulk_modulus, self.k0_prime)
            }
        }
    };
}

macro_rules! impl_four_parameter_reference_state {
    ($model:ty, $constructor:path) => {
        impl ReferenceStateEos for $model {
            fn reference_bulk_modulus(&self) -> f64 {
                self.k0
            }

            fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
                $constructor(volume, bulk_modulus, self.k0_prime, self.k0_double_prime)
            }
        }
    };
}

impl_two_parameter_reference_state!(BM2, BM2::new);
impl_two_parameter_reference_state!(NaturalStrain2, NaturalStrain2::new);
impl_three_parameter_reference_state!(BM3, BM3::new);
impl_three_parameter_reference_state!(Murnaghan, Murnaghan::new);
impl_three_parameter_reference_state!(NaturalStrain3, NaturalStrain3::new);
impl_three_parameter_reference_state!(Vinet, Vinet::new);
impl_four_parameter_reference_state!(BM4, BM4::new);
impl_four_parameter_reference_state!(ModifiedTait, ModifiedTait::new);
impl_four_parameter_reference_state!(NaturalStrain4, NaturalStrain4::new);

impl ReferenceStateEos for Holzapfel {
    fn reference_bulk_modulus(&self) -> f64 {
        self.k0
    }

    fn with_reference_state(&self, volume: f64, bulk_modulus: f64) -> EosResult<Self> {
        Self::new(volume, bulk_modulus, self.k0_prime, self.n, self.z)
    }
}

/// Constant `alpha K_T` thermal-pressure correction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearThermalPressure<R> {
    /// Reference isotherm.
    pub rt_eos: R,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Constant thermal-pressure slope in `GPa K^-1`.
    pub alpha_kt: f64,
}

impl<R: IsothermalEos> LinearThermalPressure<R> {
    /// Construct a constant-slope thermal-pressure model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid thermal parameters.
    pub fn new(rt_eos: R, tr: f64, alpha_kt: f64) -> EosResult<Self> {
        Ok(Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            alpha_kt: finite_parameter(alpha_kt, "alpha_KT")?,
        })
    }
}

impl<R: IsothermalEos> ThermalEos for LinearThermalPressure<R> {
    type Reference = R;

    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        finite_result(self.alpha_kt * (temperature - self.tr))
    }
}

/// Linear-in-temperature thermal pressure with a logarithmic volume slope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LogVolumeThermalPressure<R> {
    /// Reference isotherm.
    pub rt_eos: R,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Thermal-pressure slope at the reference volume in `GPa K^-1`.
    pub alpha_kt_ref: f64,
    /// Constant-volume derivative of bulk modulus in `GPa K^-1`.
    pub dk_dt_v: f64,
}

impl<R: IsothermalEos> LogVolumeThermalPressure<R> {
    /// Construct a logarithmic-volume thermal-pressure model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid thermal parameters.
    pub fn new(rt_eos: R, tr: f64, alpha_kt_ref: f64, dk_dt_v: f64) -> EosResult<Self> {
        Ok(Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            alpha_kt_ref: finite_parameter(alpha_kt_ref, "alpha_KT_ref")?,
            dk_dt_v: finite_parameter(dk_dt_v, "dK_dT_V")?,
        })
    }
}

impl<R: IsothermalEos> ThermalEos for LogVolumeThermalPressure<R> {
    type Reference = R;

    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let slope =
            self.alpha_kt_ref + self.dk_dt_v * (self.rt_eos.reference_volume() / volume).ln();
        finite_result(slope * (temperature - self.tr))
    }
}

/// Absolute bivariate second-order Taylor thermal pressure.
///
/// The wrapped isothermal EOS is a cold curve. With
/// `eta = 1 - V/V0`, `delta_eta = eta - eta0`, and
/// `delta_temperature = T - Tr`, the additive pressure is
///
/// `c0 + c1*delta_eta + c2*delta_temperature`
/// `+ c3*delta_eta^2/2 + c4*delta_temperature^2/2`
/// `+ c5*delta_eta*delta_temperature/2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SecondOrderTaylorThermalPressure<R> {
    /// Cold compression curve.
    pub rt_eos: R,
    /// Temperature coordinate of the Taylor expansion in kelvin.
    pub tr: f64,
    /// Compression coordinate of the Taylor expansion.
    pub eta0: f64,
    /// Constant pressure coefficient in `GPa`.
    pub c0: f64,
    /// Linear compression coefficient in `GPa`.
    pub c1: f64,
    /// Linear temperature coefficient in GPa/K.
    pub c2: f64,
    /// Quadratic compression coefficient in `GPa`.
    pub c3: f64,
    /// Quadratic temperature coefficient in GPa/K^2.
    pub c4: f64,
    /// Compression-temperature cross coefficient in GPa/K.
    pub c5: f64,
}

impl<R: IsothermalEos> SecondOrderTaylorThermalPressure<R> {
    /// Construct an absolute second-order thermal-pressure model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rt_eos: R,
        tr: f64,
        eta0: f64,
        c0: f64,
        c1: f64,
        c2: f64,
        c3: f64,
        c4: f64,
        c5: f64,
    ) -> EosResult<Self> {
        Ok(Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            eta0: finite_parameter(eta0, "eta0")?,
            c0: finite_parameter(c0, "c0")?,
            c1: finite_parameter(c1, "c1")?,
            c2: finite_parameter(c2, "c2")?,
            c3: finite_parameter(c3, "c3")?,
            c4: finite_parameter(c4, "c4")?,
            c5: finite_parameter(c5, "c5")?,
        })
    }
}

impl<R: IsothermalEos> ThermalEos for SecondOrderTaylorThermalPressure<R> {
    type Reference = R;

    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let delta_eta = 1.0 - volume / self.rt_eos.reference_volume() - self.eta0;
        let delta_temperature = temperature - self.tr;
        finite_result(
            self.c0
                + self.c1 * delta_eta
                + self.c2 * delta_temperature
                + 0.5 * self.c3 * delta_eta.powi(2)
                + 0.5 * self.c4 * delta_temperature.powi(2)
                + 0.5 * self.c5 * delta_eta * delta_temperature,
        )
    }

    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.thermal_pressure(volume, temperature)? - self.thermal_pressure(volume, self.tr)?,
        )
    }
}

/// Temperature dependence of instantaneous volumetric expansivity.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ThermalExpansionLaw {
    /// `alpha(T) = alpha0`.
    #[default]
    Constant,
    /// `alpha(T) = alpha0 + alpha1 T`.
    LinearTemperature,
    /// `alpha(T) = alpha0 + alpha1 (T - Tr)`.
    LinearReferenceTemperature,
}

/// Relationship used to construct the temperature-dependent reference volume.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReferenceVolumeLaw {
    /// Exponentially integrate the instantaneous expansivity.
    #[default]
    IntegratedExpansivity,
    /// Apply `V0(T)=V0(Tr)[1+alpha0(T-Tr)]` directly.
    LinearTemperature,
    /// Apply the Berman (1988) truncated quadratic reference-volume law.
    Berman,
}

/// EOS with temperature-dependent reference volume and bulk modulus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalReferenceState<R> {
    /// Reference-temperature EOS.
    pub rt_eos: R,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Constant or intercept expansivity in K^-1.
    pub alpha0: f64,
    /// Temperature derivative of the reference bulk modulus in `GPa K^-1`.
    pub dk_dt: f64,
    /// Linear temperature coefficient of expansivity in K^-2.
    pub alpha1: f64,
    /// Instantaneous expansivity law.
    pub thermal_expansion_law: ThermalExpansionLaw,
    /// Reference-volume construction law.
    pub reference_volume_law: ReferenceVolumeLaw,
}

impl<R: ReferenceStateEos> ThermalReferenceState<R> {
    /// Construct a temperature-dependent reference-state EOS.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid parameters or an inconsistent law pair.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rt_eos: R,
        tr: f64,
        alpha0: f64,
        dk_dt: f64,
        alpha1: f64,
        thermal_expansion_law: ThermalExpansionLaw,
        reference_volume_law: ReferenceVolumeLaw,
    ) -> EosResult<Self> {
        let model = Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            alpha0: finite_parameter(alpha0, "alpha0")?,
            dk_dt: finite_parameter(dk_dt, "dK_dT")?,
            alpha1: finite_parameter(alpha1, "alpha1")?,
            thermal_expansion_law,
            reference_volume_law,
        };
        if model.thermal_expansion_law == ThermalExpansionLaw::Constant && model.alpha1 != 0.0 {
            return Err(EosError::InvalidParameter {
                name: "alpha1",
                reason: "must be zero for constant thermal expansion",
            });
        }
        if model.reference_volume_law == ReferenceVolumeLaw::LinearTemperature
            && (model.thermal_expansion_law != ThermalExpansionLaw::Constant || model.alpha1 != 0.0)
        {
            return Err(EosError::InvalidParameter {
                name: "reference_volume_law",
                reason: "linear temperature volume requires constant expansivity configuration",
            });
        }
        if model.reference_volume_law == ReferenceVolumeLaw::Berman
            && model.thermal_expansion_law != ThermalExpansionLaw::LinearTemperature
        {
            return Err(EosError::InvalidParameter {
                name: "reference_volume_law",
                reason: "Berman volume requires linear-temperature expansivity configuration",
            });
        }
        Ok(model)
    }

    fn state_eos(&self, temperature: f64) -> EosResult<R> {
        let temperature = positive_state(temperature, "temperature")?;
        let delta = temperature - self.tr;
        let reference_volume = match self.reference_volume_law {
            ReferenceVolumeLaw::LinearTemperature => {
                self.rt_eos.reference_volume() * (1.0 + self.alpha0 * delta)
            }
            ReferenceVolumeLaw::Berman => {
                self.rt_eos.reference_volume()
                    * (1.0 + self.alpha0 * delta + 0.5 * self.alpha1 * delta * delta)
            }
            ReferenceVolumeLaw::IntegratedExpansivity => {
                let mut exponent = self.alpha0 * delta;
                if self.thermal_expansion_law == ThermalExpansionLaw::LinearTemperature {
                    exponent += 0.5 * self.alpha1 * (temperature * temperature - self.tr * self.tr);
                } else if self.thermal_expansion_law
                    == ThermalExpansionLaw::LinearReferenceTemperature
                {
                    exponent += 0.5 * self.alpha1 * delta * delta;
                }
                self.rt_eos.reference_volume() * exponent.exp()
            }
        };
        let bulk_modulus = self.rt_eos.reference_bulk_modulus() + self.dk_dt * delta;
        if !reference_volume.is_finite() || reference_volume <= 0.0 {
            return Err(EosError::InvalidState {
                name: "temperature",
                reason: "produces a non-positive reference volume",
            });
        }
        if !bulk_modulus.is_finite() || bulk_modulus <= 0.0 {
            return Err(EosError::InvalidState {
                name: "temperature",
                reason: "produces a non-positive bulk modulus",
            });
        }
        self.rt_eos
            .with_reference_state(reference_volume, bulk_modulus)
    }
}

impl<R: ReferenceStateEos> ThermalEos for ThermalReferenceState<R> {
    type Reference = R;

    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        finite_result(
            self.state_eos(temperature)?.pressure(volume)? - self.rt_eos.pressure(volume)?,
        )
    }

    fn bulk_modulus(&self, volume: f64, temperature: f64, relative_step: f64) -> EosResult<f64> {
        positive_state(relative_step, "relative_step")?;
        self.state_eos(temperature)?.bulk_modulus(volume)
    }
}

/// Tange-type asymptotic-power-law Mie--Gruneisen--Debye EOS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymptoticPowerLawMieGruneisenDebye<R> {
    /// Reference isotherm.
    pub rt_eos: R,
    /// Reference temperature in kelvin.
    pub tr: f64,
    /// Reference Debye temperature in kelvin.
    pub theta0: f64,
    /// Reference Gruneisen parameter.
    pub gamma0: f64,
    /// Fractional asymptotic coefficient.
    pub a: f64,
    /// Volume exponent.
    pub b: f64,
    /// Number of atoms per formula.
    pub n: f64,
}

impl<R: IsothermalEos> AsymptoticPowerLawMieGruneisenDebye<R> {
    /// Construct the asymptotic-power-law model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid thermal parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rt_eos: R,
        tr: f64,
        theta0: f64,
        gamma0: f64,
        a: f64,
        b: f64,
        n: f64,
    ) -> EosResult<Self> {
        let a = finite_parameter(a, "a")?;
        if !(0.0..=1.0).contains(&a) {
            return Err(EosError::InvalidParameter {
                name: "a",
                reason: "must lie between zero and one",
            });
        }
        Ok(Self {
            rt_eos,
            tr: positive_parameter(tr, "Tr")?,
            theta0: positive_parameter(theta0, "theta0")?,
            gamma0: finite_parameter(gamma0, "gamma0")?,
            a,
            b: finite_parameter(b, "b")?,
            n: positive_parameter(n, "n")?,
        })
    }

    /// Volume-dependent Gruneisen parameter.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or non-finite result.
    pub fn volume_gruneisen_parameter(&self, volume: f64) -> EosResult<f64> {
        let ratio = positive_state(volume, "volume")? / self.rt_eos.reference_volume();
        finite_result(self.gamma0 * (1.0 + self.a * (ratio.powf(self.b) - 1.0)))
    }

    /// Volume-dependent Debye temperature.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or non-finite result.
    #[allow(clippy::float_cmp)]
    pub fn characteristic_temperature(&self, volume: f64) -> EosResult<f64> {
        let ratio = positive_state(volume, "volume")? / self.rt_eos.reference_volume();
        let logarithmic_ratio = ratio.ln();
        let exponent = if self.b == 0.0 {
            -self.gamma0 * logarithmic_ratio
        } else {
            -self.gamma0
                * ((1.0 - self.a) * logarithmic_ratio
                    + self.a * (self.b * logarithmic_ratio).exp_m1() / self.b)
        };
        finite_result(self.theta0 * exponent.exp())
    }

    /// Debye vibrational energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn thermal_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let temperature = positive_state(temperature, "temperature")?;
        let theta = self.characteristic_temperature(volume)?;
        finite_result(
            3.0 * self.n * GAS_CONSTANT * temperature * debye_function_3(theta / temperature)?,
        )
    }

    /// Debye vibrational entropy in J mol^-1 K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn thermal_entropy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let temperature = positive_state(temperature, "temperature")?;
        let ratio = self.characteristic_temperature(volume)? / temperature;
        let log_term = (-(-ratio).exp_m1()).ln();
        finite_result(self.n * GAS_CONSTANT * (4.0 * debye_function_3(ratio)? - 3.0 * log_term))
    }

    /// Unreferenced vibrational pressure in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or non-finite result.
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
    /// Returns an error for an invalid state or non-finite result.
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
    /// Returns an error for an invalid state or non-finite result.
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
    /// Returns an error for an invalid state or non-finite result.
    pub fn thermal_gibbs_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.thermal_helmholtz_free_energy(volume, temperature)?
                + self.vibrational_pressure(volume, temperature)? * volume * 1.0e4,
        )
    }
}

impl<R: IsothermalEos> ThermalEos for AsymptoticPowerLawMieGruneisenDebye<R> {
    type Reference = R;

    fn reference_eos(&self) -> &R {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let energy_difference =
            self.thermal_energy(volume, temperature)? - self.thermal_energy(volume, self.tr)?;
        finite_result(self.volume_gruneisen_parameter(volume)? * energy_difference / volume / 1.0e4)
    }
}

impl<R: IsothermalEos> CaloricEos for AsymptoticPowerLawMieGruneisenDebye<R> {
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

#[derive(Clone, Copy, Debug, PartialEq)]
struct DoubleDebyeModeTerms {
    theta_a: f64,
    gamma_a: f64,
    theta_b: f64,
    gamma_b: f64,
    weight_a: f64,
    weight_a_prime: f64,
}

/// Vinet curve plus a double-Debye Helmholtz contribution.
///
/// When `tr` is `None`, the Vinet member is a motionless-ion 0 K cold curve
/// and [`Self::thermal_pressure`] is an absolute contribution including
/// zero-point pressure. When `tr` is present, the Vinet member is a complete
/// reference isotherm and the non-cold contribution is rebased to zero at it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleDebyeHelmholtz {
    /// Vinet cold curve, or complete reference isotherm when `tr` is present.
    pub rt_eos: Vinet,
    /// Characteristic reference volume in J bar^-1 mol^-1.
    pub vp: f64,
    /// First Debye cutoff at `vp`, in kelvin.
    pub theta_a0: f64,
    /// Exponential volume coefficient for the first cutoff.
    pub a_a: f64,
    /// Power-law volume coefficient for the first cutoff.
    pub b_a: f64,
    /// Second Debye cutoff at `vp`, in kelvin.
    pub theta_b0: f64,
    /// Exponential volume coefficient for the second cutoff.
    pub a_b: f64,
    /// Power-law volume coefficient for the second cutoff.
    pub b_b: f64,
    /// First phonon moment at `vp`, in kelvin.
    pub theta_1_0: f64,
    /// Exponential volume coefficient for the first phonon moment.
    pub a_1: f64,
    /// Power-law volume coefficient for the first phonon moment.
    pub b_1: f64,
    /// Number of atoms per formula unit.
    pub n: f64,
    /// Reference anharmonic coefficient in K^-1.
    pub alpha0: f64,
    /// Anharmonic reference volume in J bar^-1 mol^-1.
    pub ve: f64,
    /// Anharmonic volume exponent.
    pub kappa: f64,
    /// Cold energy at the Vinet reference volume in J mol^-1.
    pub phi0: f64,
    /// Optional complete-reference-isotherm temperature in kelvin.
    pub tr: Option<f64>,
}

impl DoubleDebyeHelmholtz {
    /// Conventional temperature used to select the nearest inversion branch.
    pub const REFERENCE_TEMPERATURE: f64 = 300.0;

    /// Construct a Vinet/double-Debye Helmholtz model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite parameters.
    #[allow(clippy::too_many_arguments, clippy::similar_names)]
    pub fn new(
        rt_eos: Vinet,
        vp: f64,
        theta_a0: f64,
        a_a: f64,
        b_a: f64,
        theta_b0: f64,
        a_b: f64,
        b_b: f64,
        theta_1_0: f64,
        a_1: f64,
        b_1: f64,
        n: f64,
        alpha0: f64,
        ve: f64,
        kappa: f64,
        phi0: f64,
    ) -> EosResult<Self> {
        let alpha0 = finite_parameter(alpha0, "alpha0")?;
        if alpha0 < 0.0 {
            return Err(EosError::InvalidParameter {
                name: "alpha0",
                reason: "must not be negative",
            });
        }
        Ok(Self {
            rt_eos,
            vp: positive_parameter(vp, "Vp")?,
            theta_a0: positive_parameter(theta_a0, "theta_a0")?,
            a_a: finite_parameter(a_a, "a_a")?,
            b_a: finite_parameter(b_a, "b_a")?,
            theta_b0: positive_parameter(theta_b0, "theta_b0")?,
            a_b: finite_parameter(a_b, "a_b")?,
            b_b: finite_parameter(b_b, "b_b")?,
            theta_1_0: positive_parameter(theta_1_0, "theta_1_0")?,
            a_1: finite_parameter(a_1, "a_1")?,
            b_1: finite_parameter(b_1, "b_1")?,
            n: positive_parameter(n, "n")?,
            alpha0,
            ve: positive_parameter(ve, "Ve")?,
            kappa: finite_parameter(kappa, "kappa")?,
            phi0: finite_parameter(phi0, "phi0")?,
            tr: None,
        })
    }

    /// Rebase the non-cold contribution onto the supplied Vinet isotherm.
    ///
    /// # Errors
    ///
    /// Returns an error unless `tr` is positive and finite.
    pub fn with_reference_temperature(mut self, tr: f64) -> EosResult<Self> {
        self.tr = Some(positive_parameter(tr, "Tr")?);
        Ok(self)
    }

    fn increment_reference_temperature(&self) -> f64 {
        self.tr.unwrap_or(Self::REFERENCE_TEMPERATURE)
    }

    fn nonnegative_temperature(temperature: f64) -> EosResult<f64> {
        let temperature = finite_state(temperature, "temperature")?;
        if temperature < 0.0 {
            Err(EosError::InvalidState {
                name: "temperature",
                reason: "must not be negative",
            })
        } else {
            Ok(temperature)
        }
    }

    fn temperature_law(&self, volume: f64, theta0: f64, a: f64, b: f64) -> EosResult<(f64, f64)> {
        let volume = positive_state(volume, "volume")?;
        let theta = theta0 * (-b * (volume / self.vp).ln() + a * (self.vp - volume)).exp();
        if !theta.is_finite() || theta <= 0.0 {
            return Err(EosError::NonFiniteResult);
        }
        Ok((theta, finite_result(a * volume + b)?))
    }

    /// Return the two Debye cutoffs and first phonon moment in kelvin.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn debye_temperatures(&self, volume: f64) -> EosResult<(f64, f64, f64)> {
        Ok((
            self.temperature_law(volume, self.theta_a0, self.a_a, self.b_a)?
                .0,
            self.temperature_law(volume, self.theta_b0, self.a_b, self.b_b)?
                .0,
            self.temperature_law(volume, self.theta_1_0, self.a_1, self.b_1)?
                .0,
        ))
    }

    #[allow(clippy::similar_names)]
    fn mode_terms(&self, volume: f64) -> EosResult<DoubleDebyeModeTerms> {
        let volume = positive_state(volume, "volume")?;
        let (theta_a, gamma_a) = self.temperature_law(volume, self.theta_a0, self.a_a, self.b_a)?;
        let (theta_b, gamma_b) = self.temperature_law(volume, self.theta_b0, self.a_b, self.b_b)?;
        let (theta_1, gamma_1) =
            self.temperature_law(volume, self.theta_1_0, self.a_1, self.b_1)?;
        let denominator = theta_b - theta_a;
        let numerator = theta_b - theta_1;
        let scale = theta_a.max(theta_b).max(theta_1);

        let (weight_a, weight_a_prime) = if denominator.abs() > 1.0e-10 * scale {
            let theta_a_prime = -gamma_a * theta_a / volume;
            let theta_b_prime = -gamma_b * theta_b / volume;
            let theta_1_prime = -gamma_1 * theta_1 / volume;
            (
                numerator / denominator,
                ((theta_b_prime - theta_1_prime) * denominator
                    - numerator * (theta_b_prime - theta_a_prime))
                    / denominator.powi(2),
            )
        } else {
            let gamma_denominator = gamma_b - gamma_a;
            let limiting_weight = if gamma_denominator.abs() > 1.0e-12 {
                (gamma_b - gamma_1) / gamma_denominator
            } else {
                0.5
            };
            (limiting_weight, 0.0)
        };
        Ok(DoubleDebyeModeTerms {
            theta_a,
            gamma_a,
            theta_b,
            gamma_b,
            weight_a: finite_result(weight_a)?,
            weight_a_prime: finite_result(weight_a_prime)?,
        })
    }

    /// Return the volume-dependent weights of Debye modes A and B.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn double_debye_weights(&self, volume: f64) -> EosResult<(f64, f64)> {
        let weight_a = self.mode_terms(volume)?.weight_a;
        Ok((weight_a, 1.0 - weight_a))
    }

    fn single_debye_free_energy(theta: f64, temperature: f64) -> EosResult<f64> {
        let temperature = Self::nonnegative_temperature(temperature)?;
        let mut result = 9.0 * GAS_CONSTANT * theta / 8.0;
        if temperature > 0.0 {
            let ratio = theta / temperature;
            result += GAS_CONSTANT
                * temperature
                * (3.0 * (-(-ratio).exp_m1()).ln() - debye_function_3(ratio)?);
        }
        finite_result(result)
    }

    fn single_debye_internal_energy(theta: f64, temperature: f64) -> EosResult<f64> {
        let temperature = Self::nonnegative_temperature(temperature)?;
        let mut result = 9.0 * GAS_CONSTANT * theta / 8.0;
        if temperature > 0.0 {
            result += 3.0 * GAS_CONSTANT * temperature * debye_function_3(theta / temperature)?;
        }
        finite_result(result)
    }

    fn single_debye_heat_capacity(theta: f64, temperature: f64) -> EosResult<f64> {
        let temperature = Self::nonnegative_temperature(temperature)?;
        if temperature == 0.0 {
            return Ok(0.0);
        }
        let ratio = theta / temperature;
        let occupation = ratio * (-ratio).exp() / (-(-ratio).exp_m1());
        finite_result(3.0 * GAS_CONSTANT * (4.0 * debye_function_3(ratio)? - 3.0 * occupation))
    }

    /// Return the Vinet cold-curve energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn cold_energy(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let delta = self.rt_eos.k0_prime - 1.0;
        let x = (volume / self.rt_eos.v0).cbrt();
        let reduced = if delta.abs() < 1.0e-7 {
            let y = x - 1.0;
            1.125 * y.powi(2) - 1.125 * delta * y.powi(3)
        } else {
            let exponent = 1.5 * delta * (x - 1.0);
            (-(-exponent).exp_m1() - exponent * (-exponent).exp()) / delta.powi(2)
        };
        finite_result(self.phi0 + 4.0 * self.rt_eos.v0 * self.rt_eos.k0 * 1.0e4 * reduced)
    }

    /// Return the weighted double-Debye zero-point energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn zero_point_energy(&self, volume: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        finite_result(
            self.n * 9.0 * GAS_CONSTANT / 8.0
                * (terms.weight_a * terms.theta_a + (1.0 - terms.weight_a) * terms.theta_b),
        )
    }

    /// Return the double-Debye ionic Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn ion_helmholtz_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        let free_a = Self::single_debye_free_energy(terms.theta_a, temperature)?;
        let free_b = Self::single_debye_free_energy(terms.theta_b, temperature)?;
        finite_result(self.n * (terms.weight_a * free_a + (1.0 - terms.weight_a) * free_b))
    }

    /// Return the volume-dependent anharmonic coefficient in K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn anharmonic_coefficient(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        finite_result(self.alpha0 * (volume / self.ve).powf(self.kappa))
    }

    /// Return the anharmonic Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or a non-finite result.
    pub fn anharmonic_helmholtz_free_energy(
        &self,
        volume: f64,
        temperature: f64,
    ) -> EosResult<f64> {
        let temperature = Self::nonnegative_temperature(temperature)?;
        finite_result(
            -0.5 * self.n
                * GAS_CONSTANT
                * self.anharmonic_coefficient(volume)?
                * temperature.powi(2),
        )
    }

    /// Return the complete Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed model evaluation.
    pub fn helmholtz_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let mut non_cold = self.ion_helmholtz_free_energy(volume, temperature)?
            + self.anharmonic_helmholtz_free_energy(volume, temperature)?;
        if let Some(tr) = self.tr {
            non_cold -= self.ion_helmholtz_free_energy(volume, tr)?
                + self.anharmonic_helmholtz_free_energy(volume, tr)?;
        }
        finite_result(self.cold_energy(volume)? + non_cold)
    }

    /// Return ionic pressure, including zero-point pressure, in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn ion_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let terms = self.mode_terms(volume)?;
        let free_a = Self::single_debye_free_energy(terms.theta_a, temperature)?;
        let free_b = Self::single_debye_free_energy(terms.theta_b, temperature)?;
        let energy_a = Self::single_debye_internal_energy(terms.theta_a, temperature)?;
        let energy_b = Self::single_debye_internal_energy(terms.theta_b, temperature)?;
        finite_result(
            self.n
                * ((terms.weight_a * terms.gamma_a * energy_a
                    + (1.0 - terms.weight_a) * terms.gamma_b * energy_b)
                    / volume
                    - terms.weight_a_prime * (free_a - free_b))
                / 1.0e4,
        )
    }

    /// Return the anharmonic pressure in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or a non-finite result.
    pub fn anharmonic_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = Self::nonnegative_temperature(temperature)?;
        finite_result(
            0.5 * self.n
                * GAS_CONSTANT
                * self.kappa
                * self.anharmonic_coefficient(volume)?
                * temperature.powi(2)
                / volume
                / 1.0e4,
        )
    }
}

impl ThermalEos for DoubleDebyeHelmholtz {
    type Reference = Vinet;

    fn reference_eos(&self) -> &Vinet {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.increment_reference_temperature()
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let mut pressure = self.ion_pressure(volume, temperature)?
            + self.anharmonic_pressure(volume, temperature)?;
        if let Some(tr) = self.tr {
            pressure -= self.ion_pressure(volume, tr)? + self.anharmonic_pressure(volume, tr)?;
        }
        finite_result(pressure)
    }

    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        if self.tr.is_some() {
            return self.thermal_pressure(volume, temperature);
        }
        let reference_temperature = self.increment_reference_temperature();
        finite_result(
            self.pressure(volume, temperature)? - self.pressure(volume, reference_temperature)?,
        )
    }

    fn dac_thermal_pressure(&self, volume: f64, temperature: f64, f_dac: f64) -> EosResult<f64> {
        let f_dac = finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        finite_result(f_dac * self.thermal_pressure_increment(volume, temperature)?)
    }

    fn temperature_from_volumes(
        &self,
        ambient_volume: f64,
        heated_volume: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let ambient_volume = positive_state(ambient_volume, "volume")?;
        let heated_volume = positive_state(heated_volume, "volume")?;
        let f_dac = finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        let reference_temperature = self.increment_reference_temperature();
        let heated_reference_pressure = self.pressure(heated_volume, reference_temperature)?;
        let target = (self.pressure(ambient_volume, reference_temperature)?
            - heated_reference_pressure)
            / (1.0 - f_dac);
        if target < 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        let temperature = solve_temperature_function(
            |value| {
                self.pressure(heated_volume, value)
                    .map(|pressure| pressure - heated_reference_pressure)
            },
            target,
            reference_temperature,
        )?;
        let tolerance = 1.0e-10 * reference_temperature;
        if temperature < reference_temperature - tolerance {
            Err(EosError::OutsideInvertibleRange)
        } else {
            Ok(temperature)
        }
    }
}

impl CaloricEos for DoubleDebyeHelmholtz {
    fn molar_heat_capacity_v(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        finite_result(
            self.n
                * (terms.weight_a * Self::single_debye_heat_capacity(terms.theta_a, temperature)?
                    + (1.0 - terms.weight_a)
                        * Self::single_debye_heat_capacity(terms.theta_b, temperature)?
                    + GAS_CONSTANT
                        * self.anharmonic_coefficient(volume)?
                        * Self::nonnegative_temperature(temperature)?),
        )
    }
}

/// Vinet/double-Debye Helmholtz EOS constrained by the logarithmic phonon moment.
///
/// The mode weights conserve `theta_0` as in Correa et al. (2008), equation
/// 13. The anharmonic coefficient is volume independent, so it contributes to
/// free energy and heat capacity but not pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleDebyeLogMomentHelmholtz {
    /// Vinet cold curve, or complete reference isotherm when `tr` is present.
    pub rt_eos: Vinet,
    /// Characteristic reference volume in J bar^-1 mol^-1.
    pub vp: f64,
    /// First Debye cutoff at `vp`, in kelvin.
    pub theta_a0: f64,
    /// Exponential volume coefficient for the first cutoff.
    pub a_a: f64,
    /// Power-law volume coefficient for the first cutoff.
    pub b_a: f64,
    /// Second Debye cutoff at `vp`, in kelvin.
    pub theta_b0: f64,
    /// Exponential volume coefficient for the second cutoff.
    pub a_b: f64,
    /// Power-law volume coefficient for the second cutoff.
    pub b_b: f64,
    /// Logarithmic phonon moment at `vp`, in kelvin.
    pub theta_0_0: f64,
    /// Exponential volume coefficient for the logarithmic phonon moment.
    pub a_0: f64,
    /// Power-law volume coefficient for the logarithmic phonon moment.
    pub b_0: f64,
    /// Number of atoms per formula unit.
    pub n: f64,
    /// Volume-independent anharmonic coefficient in K^-1.
    pub anharmonic_a: f64,
    /// Cold energy at the Vinet reference volume in J mol^-1.
    pub phi0: f64,
    /// Optional complete-reference-isotherm temperature in kelvin.
    pub tr: Option<f64>,
}

impl DoubleDebyeLogMomentHelmholtz {
    /// Conventional temperature used to select the nearest inversion branch.
    pub const REFERENCE_TEMPERATURE: f64 = 300.0;

    /// Construct a logarithmic-moment double-Debye Helmholtz model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid or non-finite parameters.
    #[allow(clippy::too_many_arguments, clippy::similar_names)]
    pub fn new(
        rt_eos: Vinet,
        vp: f64,
        theta_a0: f64,
        a_a: f64,
        b_a: f64,
        theta_b0: f64,
        a_b: f64,
        b_b: f64,
        theta_0_0: f64,
        a_0: f64,
        b_0: f64,
        n: f64,
        anharmonic_a: f64,
        phi0: f64,
    ) -> EosResult<Self> {
        let anharmonic_a = finite_parameter(anharmonic_a, "anharmonic_a")?;
        if anharmonic_a < 0.0 {
            return Err(EosError::InvalidParameter {
                name: "anharmonic_a",
                reason: "must not be negative",
            });
        }
        Ok(Self {
            rt_eos,
            vp: positive_parameter(vp, "Vp")?,
            theta_a0: positive_parameter(theta_a0, "theta_a0")?,
            a_a: finite_parameter(a_a, "a_a")?,
            b_a: finite_parameter(b_a, "b_a")?,
            theta_b0: positive_parameter(theta_b0, "theta_b0")?,
            a_b: finite_parameter(a_b, "a_b")?,
            b_b: finite_parameter(b_b, "b_b")?,
            theta_0_0: positive_parameter(theta_0_0, "theta_0_0")?,
            a_0: finite_parameter(a_0, "a_0")?,
            b_0: finite_parameter(b_0, "b_0")?,
            n: positive_parameter(n, "n")?,
            anharmonic_a,
            phi0: finite_parameter(phi0, "phi0")?,
            tr: None,
        })
    }

    /// Rebase the non-cold contribution onto the supplied Vinet isotherm.
    ///
    /// # Errors
    ///
    /// Returns an error unless `tr` is positive and finite.
    pub fn with_reference_temperature(mut self, tr: f64) -> EosResult<Self> {
        self.tr = Some(positive_parameter(tr, "Tr")?);
        Ok(self)
    }

    fn increment_reference_temperature(&self) -> f64 {
        self.tr.unwrap_or(Self::REFERENCE_TEMPERATURE)
    }

    fn temperature_law(&self, volume: f64, theta0: f64, a: f64, b: f64) -> EosResult<(f64, f64)> {
        let volume = positive_state(volume, "volume")?;
        let theta = theta0 * (-b * (volume / self.vp).ln() + a * (self.vp - volume)).exp();
        if !theta.is_finite() || theta <= 0.0 {
            return Err(EosError::NonFiniteResult);
        }
        Ok((theta, finite_result(a * volume + b)?))
    }

    /// Return the two Debye cutoffs and logarithmic phonon moment in kelvin.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn debye_temperatures(&self, volume: f64) -> EosResult<(f64, f64, f64)> {
        Ok((
            self.temperature_law(volume, self.theta_a0, self.a_a, self.b_a)?
                .0,
            self.temperature_law(volume, self.theta_b0, self.a_b, self.b_b)?
                .0,
            self.temperature_law(volume, self.theta_0_0, self.a_0, self.b_0)?
                .0,
        ))
    }

    #[allow(clippy::similar_names)]
    fn mode_terms(&self, volume: f64) -> EosResult<DoubleDebyeModeTerms> {
        let volume = positive_state(volume, "volume")?;
        let (theta_a, gamma_a) = self.temperature_law(volume, self.theta_a0, self.a_a, self.b_a)?;
        let (theta_b, gamma_b) = self.temperature_law(volume, self.theta_b0, self.a_b, self.b_b)?;
        let (theta_0, gamma_0) =
            self.temperature_law(volume, self.theta_0_0, self.a_0, self.b_0)?;
        let denominator = (theta_b / theta_a).ln();
        let numerator = (theta_b / theta_0).ln();

        let (weight_a, weight_a_prime) = if denominator.abs() > 1.0e-10 {
            let numerator_prime = (gamma_0 - gamma_b) / volume;
            let denominator_prime = (gamma_a - gamma_b) / volume;
            (
                numerator / denominator,
                (numerator_prime * denominator - numerator * denominator_prime)
                    / denominator.powi(2),
            )
        } else {
            let gamma_denominator = gamma_b - gamma_a;
            let limiting_weight = if gamma_denominator.abs() > 1.0e-12 {
                (gamma_b - gamma_0) / gamma_denominator
            } else {
                0.5
            };
            (limiting_weight, 0.0)
        };
        Ok(DoubleDebyeModeTerms {
            theta_a,
            gamma_a,
            theta_b,
            gamma_b,
            weight_a: finite_result(weight_a)?,
            weight_a_prime: finite_result(weight_a_prime)?,
        })
    }

    /// Return the volume-dependent weights of Debye modes A and B.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn double_debye_weights(&self, volume: f64) -> EosResult<(f64, f64)> {
        let weight_a = self.mode_terms(volume)?.weight_a;
        Ok((weight_a, 1.0 - weight_a))
    }

    /// Return the Vinet cold-curve energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn cold_energy(&self, volume: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let delta = self.rt_eos.k0_prime - 1.0;
        let x = (volume / self.rt_eos.v0).cbrt();
        let reduced = if delta.abs() < 1.0e-7 {
            let y = x - 1.0;
            1.125 * y.powi(2) - 1.125 * delta * y.powi(3)
        } else {
            let exponent = 1.5 * delta * (x - 1.0);
            (-(-exponent).exp_m1() - exponent * (-exponent).exp()) / delta.powi(2)
        };
        finite_result(self.phi0 + 4.0 * self.rt_eos.v0 * self.rt_eos.k0 * 1.0e4 * reduced)
    }

    /// Return the weighted double-Debye zero-point energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume or a non-finite result.
    pub fn zero_point_energy(&self, volume: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        finite_result(
            self.n * 9.0 * GAS_CONSTANT / 8.0
                * (terms.weight_a * terms.theta_a + (1.0 - terms.weight_a) * terms.theta_b),
        )
    }

    /// Return the double-Debye ionic Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn ion_helmholtz_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        let free_a = DoubleDebyeHelmholtz::single_debye_free_energy(terms.theta_a, temperature)?;
        let free_b = DoubleDebyeHelmholtz::single_debye_free_energy(terms.theta_b, temperature)?;
        finite_result(self.n * (terms.weight_a * free_a + (1.0 - terms.weight_a) * free_b))
    }

    /// Return the volume-independent anharmonic coefficient in K^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid volume.
    pub fn anharmonic_coefficient(&self, volume: f64) -> EosResult<f64> {
        positive_state(volume, "volume")?;
        Ok(self.anharmonic_a)
    }

    /// Return the volume-independent anharmonic Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or a non-finite result.
    pub fn anharmonic_helmholtz_free_energy(
        &self,
        volume: f64,
        temperature: f64,
    ) -> EosResult<f64> {
        self.anharmonic_coefficient(volume)?;
        let temperature = DoubleDebyeHelmholtz::nonnegative_temperature(temperature)?;
        finite_result(-self.n * GAS_CONSTANT * self.anharmonic_a * temperature.powi(2))
    }

    /// Return zero anharmonic pressure in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state.
    pub fn anharmonic_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        self.anharmonic_coefficient(volume)?;
        DoubleDebyeHelmholtz::nonnegative_temperature(temperature)?;
        Ok(0.0)
    }

    /// Return the complete Helmholtz energy in J mol^-1.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed model evaluation.
    pub fn helmholtz_free_energy(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let mut non_cold = self.ion_helmholtz_free_energy(volume, temperature)?
            + self.anharmonic_helmholtz_free_energy(volume, temperature)?;
        if let Some(tr) = self.tr {
            non_cold -= self.ion_helmholtz_free_energy(volume, tr)?
                + self.anharmonic_helmholtz_free_energy(volume, tr)?;
        }
        finite_result(self.cold_energy(volume)? + non_cold)
    }

    /// Return ionic pressure, including zero-point pressure, in `GPa`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid state or failed Debye evaluation.
    pub fn ion_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let terms = self.mode_terms(volume)?;
        let free_a = DoubleDebyeHelmholtz::single_debye_free_energy(terms.theta_a, temperature)?;
        let free_b = DoubleDebyeHelmholtz::single_debye_free_energy(terms.theta_b, temperature)?;
        let energy_a =
            DoubleDebyeHelmholtz::single_debye_internal_energy(terms.theta_a, temperature)?;
        let energy_b =
            DoubleDebyeHelmholtz::single_debye_internal_energy(terms.theta_b, temperature)?;
        finite_result(
            self.n
                * ((terms.weight_a * terms.gamma_a * energy_a
                    + (1.0 - terms.weight_a) * terms.gamma_b * energy_b)
                    / volume
                    - terms.weight_a_prime * (free_a - free_b))
                / 1.0e4,
        )
    }
}

impl ThermalEos for DoubleDebyeLogMomentHelmholtz {
    type Reference = Vinet;

    fn reference_eos(&self) -> &Vinet {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.increment_reference_temperature()
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let mut pressure = self.ion_pressure(volume, temperature)?
            + self.anharmonic_pressure(volume, temperature)?;
        if let Some(tr) = self.tr {
            pressure -= self.ion_pressure(volume, tr)? + self.anharmonic_pressure(volume, tr)?;
        }
        finite_result(pressure)
    }

    fn thermal_pressure_increment(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        if self.tr.is_some() {
            return self.thermal_pressure(volume, temperature);
        }
        let reference_temperature = self.increment_reference_temperature();
        finite_result(
            self.pressure(volume, temperature)? - self.pressure(volume, reference_temperature)?,
        )
    }

    fn dac_thermal_pressure(&self, volume: f64, temperature: f64, f_dac: f64) -> EosResult<f64> {
        let f_dac = finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        finite_result(f_dac * self.thermal_pressure_increment(volume, temperature)?)
    }

    fn temperature_from_volumes(
        &self,
        ambient_volume: f64,
        heated_volume: f64,
        f_dac: f64,
    ) -> EosResult<f64> {
        let ambient_volume = positive_state(ambient_volume, "volume")?;
        let heated_volume = positive_state(heated_volume, "volume")?;
        let f_dac = finite_state(f_dac, "f_dac")?;
        if !(0.0..1.0).contains(&f_dac) {
            return Err(EosError::InvalidState {
                name: "f_dac",
                reason: "must lie in [0, 1)",
            });
        }
        let reference_temperature = self.increment_reference_temperature();
        let heated_reference_pressure = self.pressure(heated_volume, reference_temperature)?;
        let target = (self.pressure(ambient_volume, reference_temperature)?
            - heated_reference_pressure)
            / (1.0 - f_dac);
        if target < 0.0 {
            return Err(EosError::OutsideInvertibleRange);
        }
        let temperature = solve_temperature_function(
            |value| {
                self.pressure(heated_volume, value)
                    .map(|pressure| pressure - heated_reference_pressure)
            },
            target,
            reference_temperature,
        )?;
        let tolerance = 1.0e-10 * reference_temperature;
        if temperature < reference_temperature - tolerance {
            Err(EosError::OutsideInvertibleRange)
        } else {
            Ok(temperature)
        }
    }
}

impl CaloricEos for DoubleDebyeLogMomentHelmholtz {
    fn molar_heat_capacity_v(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let terms = self.mode_terms(volume)?;
        let temperature = DoubleDebyeHelmholtz::nonnegative_temperature(temperature)?;
        finite_result(
            self.n
                * (terms.weight_a
                    * DoubleDebyeHelmholtz::single_debye_heat_capacity(
                        terms.theta_a,
                        temperature,
                    )?
                    + (1.0 - terms.weight_a)
                        * DoubleDebyeHelmholtz::single_debye_heat_capacity(
                            terms.theta_b,
                            temperature,
                        )?
                    + 2.0 * GAS_CONSTANT * self.anharmonic_a * temperature),
        )
    }
}

/// Thermal parameters of the Dorogokupets--Oganov (2007) Helmholtz model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DorogokupetsOganov2007Parameters {
    pub tr: f64,
    pub theta_b1: f64,
    pub d_b1: f64,
    pub m_b1: f64,
    pub theta_b2: f64,
    pub d_b2: f64,
    pub m_b2: f64,
    pub theta_e1: f64,
    pub m_e1: f64,
    pub theta_e2: f64,
    pub m_e2: f64,
    pub gamma0: f64,
    pub gamma_inf: f64,
    pub beta: f64,
    /// Table-I intrinsic-anharmonicity coefficient in `10^-6 K^-1`.
    pub anharmonic_a: f64,
    pub anharmonic_m: f64,
    /// Table-I electronic coefficient in `10^-6 K^-1`.
    pub electronic_e: f64,
    pub electronic_g: f64,
    /// Vacancy formation enthalpy divided by the gas constant, in kelvin.
    pub defect_h: f64,
    pub defect_s: f64,
}

/// Dorogokupets--Oganov (2007) four-oscillator thermal equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DorogokupetsOganov2007<R> {
    pub rt_eos: R,
    pub parameters: DorogokupetsOganov2007Parameters,
    /// Number of atoms per chemical formula.
    pub n: f64,
}

impl<R: IsothermalEos> DorogokupetsOganov2007<R> {
    /// Construct the complete equations (7)--(14) model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid temperatures, dispersions,
    /// multiplicities, Gruneisen coefficients, or atom count.
    pub fn new(rt_eos: R, parameters: DorogokupetsOganov2007Parameters, n: f64) -> EosResult<Self> {
        let parameters = DorogokupetsOganov2007Parameters {
            tr: positive_parameter(parameters.tr, "Tr")?,
            theta_b1: positive_parameter(parameters.theta_b1, "theta_B1")?,
            d_b1: positive_parameter(parameters.d_b1, "d_B1")?,
            m_b1: positive_parameter(parameters.m_b1, "m_B1")?,
            theta_b2: positive_parameter(parameters.theta_b2, "theta_B2")?,
            d_b2: positive_parameter(parameters.d_b2, "d_B2")?,
            m_b2: positive_parameter(parameters.m_b2, "m_B2")?,
            theta_e1: positive_parameter(parameters.theta_e1, "theta_E1")?,
            m_e1: positive_parameter(parameters.m_e1, "m_E1")?,
            theta_e2: positive_parameter(parameters.theta_e2, "theta_E2")?,
            m_e2: positive_parameter(parameters.m_e2, "m_E2")?,
            gamma0: positive_parameter(parameters.gamma0, "gamma0")?,
            gamma_inf: positive_parameter(parameters.gamma_inf, "gamma_inf")?,
            beta: positive_parameter(parameters.beta, "beta")?,
            anharmonic_a: finite_parameter(parameters.anharmonic_a, "anharmonic_a")?,
            anharmonic_m: finite_parameter(parameters.anharmonic_m, "anharmonic_m")?,
            electronic_e: finite_parameter(parameters.electronic_e, "electronic_e")?,
            electronic_g: finite_parameter(parameters.electronic_g, "electronic_g")?,
            defect_h: positive_parameter(parameters.defect_h, "defect_H")?,
            defect_s: finite_parameter(parameters.defect_s, "defect_S")?,
        };
        let n = positive_parameter(n, "n")?;
        let multiplicity = parameters.m_b1 + parameters.m_b2 + parameters.m_e1 + parameters.m_e2;
        if (multiplicity - 3.0 * n).abs() > 5.0e-3 {
            return Err(EosError::InvalidParameter {
                name: "oscillator multiplicities",
                reason: "must sum to three times the atom count",
            });
        }
        Ok(Self {
            rt_eos,
            parameters,
            n,
        })
    }

    fn gamma(&self, ratio: f64) -> f64 {
        self.parameters.gamma_inf
            + (self.parameters.gamma0 - self.parameters.gamma_inf)
                * ratio.powf(self.parameters.beta)
    }

    fn theta(&self, theta0: f64, ratio: f64) -> f64 {
        theta0
            * ratio.powf(-self.parameters.gamma_inf)
            * (((self.parameters.gamma0 - self.parameters.gamma_inf) / self.parameters.beta)
                * (1.0 - ratio.powf(self.parameters.beta)))
            .exp()
    }

    fn anharmonic_bracket(theta: f64, temperature: f64) -> (f64, f64) {
        let exponent = theta / temperature;
        let decay = (-exponent).exp();
        let occupation = decay / (-(-exponent).exp_m1());
        let fluctuation = occupation * (occupation + 1.0);
        let energy = theta * (0.5 + occupation);
        let energy_derivative = 0.5 + occupation - theta / temperature * fluctuation;
        let bracket = energy * energy + 2.0 * theta * theta * fluctuation;
        let derivative = 2.0 * energy * energy_derivative + 4.0 * theta * fluctuation
            - 2.0 * theta * theta / temperature * fluctuation * (2.0 * occupation + 1.0);
        (bracket, derivative)
    }

    fn mode_pressure(
        &self,
        mode: (f64, f64, Option<f64>),
        ratio: f64,
        gamma: f64,
        volume: f64,
        temperature: f64,
    ) -> EosResult<(f64, f64)> {
        let (theta0, multiplicity, dispersion) = mode;
        let theta = self.theta(theta0, ratio);
        let energy = match dispersion {
            Some(value) => bose_mode_energy(theta, temperature, value)?,
            None => sokolova_einstein_energy(theta, temperature)?,
        };
        let quasiharmonic = multiplicity * GAS_CONSTANT * energy * gamma / volume;
        let (bracket, bracket_derivative) = Self::anharmonic_bracket(theta, temperature);
        let anharmonic_derivative = multiplicity
            * GAS_CONSTANT
            * self.parameters.anharmonic_a
            * 1.0e-6
            * ratio.powf(self.parameters.anharmonic_m)
            / (6.0 * volume)
            * (self.parameters.anharmonic_m * bracket - gamma * theta * bracket_derivative);
        Ok((quasiharmonic, anharmonic_derivative))
    }

    fn absolute_nonreference_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        let volume = positive_state(volume, "volume")?;
        let temperature = positive_state(temperature, "temperature")?;
        let ratio = volume / self.rt_eos.reference_volume();
        let gamma = self.gamma(ratio);
        let modes = [
            (
                self.parameters.theta_b1,
                self.parameters.m_b1,
                Some(self.parameters.d_b1),
            ),
            (
                self.parameters.theta_b2,
                self.parameters.m_b2,
                Some(self.parameters.d_b2),
            ),
            (self.parameters.theta_e1, self.parameters.m_e1, None),
            (self.parameters.theta_e2, self.parameters.m_e2, None),
        ];
        let mut quasiharmonic = 0.0;
        let mut anharmonic_derivative = 0.0;
        for (theta0, multiplicity, dispersion) in modes {
            let (mode_quasiharmonic, mode_anharmonic) = self.mode_pressure(
                (theta0, multiplicity, dispersion),
                ratio,
                gamma,
                volume,
                temperature,
            )?;
            quasiharmonic += mode_quasiharmonic;
            anharmonic_derivative += mode_anharmonic;
        }
        let electronic = 1.5
            * self.n
            * GAS_CONSTANT
            * self.parameters.electronic_e
            * 1.0e-6
            * self.parameters.electronic_g
            * ratio.powf(self.parameters.electronic_g)
            * temperature
            * temperature
            / volume;
        let defect_exponent = self.parameters.defect_s / ratio
            - self.parameters.defect_h / (temperature * ratio * ratio);
        let defect_derivative = -self.parameters.defect_s / ratio
            + 2.0 * self.parameters.defect_h / (temperature * ratio * ratio);
        let defect =
            1.5 * self.n * GAS_CONSTANT * temperature * defect_exponent.exp() * defect_derivative
                / volume;
        finite_result((quasiharmonic - anharmonic_derivative + electronic + defect) / 1.0e4)
    }
}

impl<R: IsothermalEos> ThermalEos for DorogokupetsOganov2007<R> {
    type Reference = R;

    fn reference_eos(&self) -> &Self::Reference {
        &self.rt_eos
    }

    fn reference_temperature(&self) -> f64 {
        self.parameters.tr
    }

    fn thermal_pressure(&self, volume: f64, temperature: f64) -> EosResult<f64> {
        finite_result(
            self.absolute_nonreference_pressure(volume, temperature)?
                - self.absolute_nonreference_pressure(volume, self.parameters.tr)?,
        )
    }
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

/// Generic multi-oscillator Gruneisen thermal-pressure equation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiOscillatorGruneisen<R> {
    /// Freely chosen reference isotherm.
    pub rt_eos: R,
    /// Validated thermal parameters.
    pub parameters: SokolovaParameters,
    /// Number of atoms per chemical formula.
    pub n: f64,
}

/// Compatibility alias for the historical Holzapfel-based public Rust type.
pub type Sokolova2016 = MultiOscillatorGruneisen<Holzapfel>;

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

impl MultiOscillatorGruneisen<Holzapfel> {
    /// Construct the Sokolova thermal-pressure model.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid characteristic temperatures,
    /// dispersions, multiplicities, or other non-finite parameters.
    pub fn new(rt_eos: Holzapfel, parameters: SokolovaParameters) -> EosResult<Self> {
        let n = rt_eos.n;
        Self::new_with_atom_count(rt_eos, parameters, n)
    }
}

impl<R: IsothermalEos> MultiOscillatorGruneisen<R> {
    /// Construct the model with an explicit formula atom count.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid characteristic temperatures,
    /// dispersions, multiplicities, atom count, or other non-finite parameters.
    pub fn new_with_atom_count(
        rt_eos: R,
        parameters: SokolovaParameters,
        n: f64,
    ) -> EosResult<Self> {
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
        Ok(Self {
            rt_eos,
            parameters,
            n: positive_parameter(n, "n")?,
        })
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
                let volume = ratio * self.rt_eos.reference_volume();
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
        let x = volume / self.rt_eos.reference_volume();
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
        let squared_temperature_coefficient = 1.5 * self.n * GAS_CONSTANT / 1_000_000.0 / volume
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

impl<R: IsothermalEos> ThermalEos for MultiOscillatorGruneisen<R> {
    type Reference = R;

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
