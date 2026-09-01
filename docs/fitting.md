# Fitting P-V and P-V-T data

Peritheos uses bounded nonlinear least squares and reports covariance,
correlation, residual, and information-criterion diagnostics. Measurement
uncertainties can be supplied for every observed state variable.

Fits using a named loss and an exact built-in Peritheos EOS run end-to-end in
Rust after one transfer of the input arrays. Custom Python EOS classes and
subclasses retain callback evaluation, and callable loss functions retain the
SciPy compatibility path. These backend choices do not change the public
fitting functions or result schema.

## Objective and diagnostics

With pressure as the only uncertain observation, the normalized residual for
point $i$ is

\[
r_i=\frac{P_{\mathrm{model},i}-P_{\mathrm{observed},i}}{\sigma_{P,i}}.
\]

If no pressure uncertainties are supplied, Peritheos sets the divisor to one
and the residuals retain pressure units. For a supplied within-observation
covariance $\Sigma_i$, the pressure and adjusted-coordinate residual vector
$\boldsymbol d_i$ is whitened using its Cholesky factor
$\Sigma_i=L_iL_i^{\mathsf T}$:

\[
\boldsymbol r_i=L_i^{-1}\boldsymbol d_i.
\]

At the solution, the reported least-squares statistics are

\[
\chi^2=\sum_j r_j^2,
\qquad
\nu=N_r-N_p,
\qquad
\chi^2_\nu=\frac{\chi^2}{\nu}.
\]

Here $N_r$ is the total number of residual components and $N_p$ is the
number of optimized scalar values. In an errors-in-variables fit, $N_p$
includes the adjusted latent volume and temperature values as well as EOS
parameters. Reduced chi-square is reported only when $\nu>0$.

For a linear loss, the local free-parameter covariance is obtained from the
Jacobian information matrix after profiling out latent observation coordinates:

\[
\Sigma_\theta\simeq
\left(J_\theta^{\mathsf T}J_\theta-
J_\theta^{\mathsf T}J_x
(J_x^{\mathsf T}J_x)^{-1}
J_x^{\mathsf T}J_\theta\right)^{-1}.
\]

Without latent coordinates, this reduces to
$(J_\theta^{\mathsf T}J_\theta)^{-1}$. The implementation uses a pseudoinverse
with a numerical rank tolerance where necessary. Observation-local latent
blocks remain structured during covariance profiling rather than being
assembled into a dense matrix. Unless
`absolute_sigma=True`, this covariance is multiplied by $\chi^2_\nu$. The
reported information criteria use

\[
\mathrm{AIC}=N_r\ln\left(\frac{\chi^2}{N_r}\right)+2N_p,
\qquad
\mathrm{BIC}=N_r\ln\left(\frac{\chi^2}{N_r}\right)+N_p\ln N_r.
\]

These are comparative diagnostics for fits to the same observations and error
model. They are not absolute tests that an EOS is physically adequate.

## Isothermal fitting

```python
from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos

result = fit_rt_eos(
    BM3,
    volume=volumes,
    pressure=pressures,
    initial={"V0": 10.0, "K0": 120.0, "K0_prime": 4.0},
    bounds={
        "V0": (9.0, 11.0),
        "K0": (50.0, 250.0),
        "K0_prime": (1.0, 10.0),
    },
    pressure_sigma=pressure_uncertainties,
    volume_sigma=volume_uncertainties,
    absolute_sigma=True,
)

print(result.parameters)
print(result.standard_errors)
print(result.correlation)
fitted_pressures = result.model.pressure(volumes)
```

## Reporting and export

`FitResult.summary()` produces a compact report with free and fixed parameters,
standard errors, fit diagnostics, and solver status:

```python
print(result.summary())
```

Machine-readable results use a versioned, JSON-safe schema:

```python
payload = result.to_dict()
json_text = result.to_json()
result.to_json("bm3-fit.json")
```

The export contains the model's module, class, and reconstructable parameter
values; free-parameter ordering; covariance and correlation matrices; raw and
weighted residuals; adjusted observations; diagnostics; and solver metadata.
Non-finite diagnostics, such as reduced chi-square for a fit with no degrees of
freedom, are represented as JSON `null`. The export is intended as a durable
analysis record; it does not dynamically import and execute model classes.

The fit covariance can be propagated into subsequent EOS calculations without
supplying the parameter errors again:

```python
eos_uncertainty = result.eos_uncertainty()
pressure_prediction = eos_uncertainty.pressure(new_volumes)
print(pressure_prediction.standard_error)
```

See [Uncertainty in EOS calculations](uncertainty.md) for correlated parameter
propagation, Monte Carlo intervals, and reference-EOS uncertainty.

Parameters in `fixed` are passed to the constructor but not optimized:

```python
result = fit_rt_eos(
    BM3,
    volumes,
    pressures,
    initial={"K0": 120.0, "K0_prime": 4.0},
    fixed={"V0": 10.0},
)
```

## Thermal fitting

The reference EOS is deliberately held fixed so the thermal and reference
isotherm regressions remain inspectable:

```python
from peritheos.eos.thermal import MieGruneisenEinstein
from peritheos.fitting import fit_thermal_eos

result = fit_thermal_eos(
    MieGruneisenEinstein,
    rt_eos=reference_eos,
    volume=volumes,
    temperature=temperatures,
    pressure=pressures,
    initial={"gamma0": 1.5, "q": 1.0},
    fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
    bounds={"gamma0": (0.0, 3.0), "q": (-2.0, 4.0)},
    pressure_sigma=pressure_uncertainties,
    volume_sigma=volume_uncertainties,
    temperature_sigma=temperature_uncertainties,
    absolute_sigma=True,
)
```

## Joint reference and thermal fitting

Use `fit_joint_eos` when a single P-V-T dataset should constrain both the
reference isotherm and the thermal model. Reference parameters use the same
`rt_eos.*` names returned by `ThermalEOS.parameter_values()`, while thermal
parameters retain their constructor names:

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenEinstein
from peritheos.fitting import fit_joint_eos

result = fit_joint_eos(
    MieGruneisenEinstein,
    BM3,
    volume=volumes,
    temperature=temperatures,
    pressure=pressures,
    initial={
        "rt_eos.V0": 10.0,
        "rt_eos.K0": 120.0,
        "rt_eos.K0_prime": 4.0,
        "gamma0": 1.5,
        "q": 1.0,
    },
    fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
    bounds={
        "rt_eos.K0": (50.0, 250.0),
        "gamma0": (0.0, 3.0),
        "q": (-2.0, 4.0),
    },
    pressure_sigma=pressure_uncertainties,
    absolute_sigma=True,
)
```

`result.model.rt_eos` is the fitted reference EOS. The rows and columns of
`result.covariance` follow `result.free_parameters` and include correlations
between reference and thermal parameters. Consequently,
`result.eos_uncertainty()` propagates the complete joint covariance without an
independence assumption.

## Errors in pressure, volume, and temperature

Pressure is the dependent observation. When only `pressure_sigma` is given,
the objective contains the familiar normalized pressure residuals.

When `volume_sigma` or `temperature_sigma` is supplied, Peritheos performs an
errors-in-variables fit. For each measurement it also fits an adjusted latent
volume and/or temperature. For a P-V-T observation, the minimized components
are:

\[
\begin{aligned}
\frac{P_{\mathrm{model}}(V^*,T^*)-P_{\mathrm{observed}}}
     {\sigma_P}, \qquad
\frac{V^*-V_{\mathrm{observed}}}{\sigma_V}, \qquad
\frac{T^*-T_{\mathrm{observed}}}{\sigma_T}.
\end{aligned}
\]

This is a nonlinear orthogonal-distance formulation for independent Gaussian
measurement errors. It does not merely evaluate the EOS at the measured
volume and temperature and pretend those coordinates are exact. Parameter
bounds continue to apply.

The fitted observation coordinates and corrections are available as:

```python
result.adjusted_volume
result.volume_corrections
result.adjusted_temperature
result.temperature_corrections
```

`weighted_residuals` contains the pressure components first, followed by the
volume and temperature correction components when present. `residuals`
contains the unscaled pressure residuals evaluated at the adjusted state.

When separate sigma arguments are used, errors in P, V, and T are treated as
mutually independent. Correlated errors within each observation can instead be
supplied as `observation_covariance`. For an isothermal fit, component order is
`(pressure, volume)` and the accepted shape is `(2, 2)` or
`pressure.shape + (2, 2)`. For thermal and joint fits, order is
`(pressure, volume, temperature)` and the final dimensions are `(3, 3)`:

```python
standard_deviations = np.array([pressure_error, volume_error, temperature_error])
covariance = correlation * np.outer(standard_deviations, standard_deviations)

result = fit_joint_eos(
    ...,
    observation_covariance=covariance,
    absolute_sigma=True,
)
```

A single matrix is applied to every observation; an array of matrices supplies
different errors and correlations per observation. Matrices must be finite,
symmetric, and positive definite. `observation_covariance` cannot be combined
with individual sigma arguments. Correlations between different observations
are not represented.

## Robust losses

Every fitting function exposes SciPy's robust least-squares losses:

```python
result = fit_rt_eos(
    ...,
    loss="soft_l1",
    f_scale=1.0,
    max_nfev=2000,
)
```

Available named losses are `linear` (the default), `soft_l1`, `huber`,
`cauchy`, and `arctan`; a SciPy-compatible callable is also accepted.
`f_scale` defines the soft outlier threshold in weighted-residual units, so it
should be chosen relative to the supplied measurement errors. Robust losses
reduce the influence of outliers but do not replace investigation or
documentation of anomalous measurements. `chi_square`, AIC, and BIC remain
least-squares diagnostics of the final weighted residuals even when a robust
loss determined the optimum.

## Uncertainty semantics

- With no uncertainty arguments, pressure residuals are unweighted and
  covariance is scaled by the fitted reduced chi-square.
- With uncertainty arguments and `absolute_sigma=True`, supplied standard
  deviations define the absolute covariance scale.
- With uncertainty arguments and `absolute_sigma=False`, they provide relative
  weights and covariance is rescaled by reduced chi-square.
- `sigma` remains a compatibility alias for `pressure_sigma`; do not pass both.
- Fixed parameters have a reported standard error of zero and are excluded
  from the covariance matrix.
- Joint-fit covariance includes reference/thermal cross-terms. A separate
  reference fit can only be combined under an explicit block-independence
  assumption.

All supplied uncertainties must be finite, positive, and scalar or
broadcastable to the pressure array.

Always inspect `correlation`, residuals versus pressure, and results under
reasonable changes in fitting range. A small residual does not make a strongly
correlated higher-order parameter physically meaningful.
