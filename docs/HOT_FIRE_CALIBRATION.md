# Hot-fire uncertainty calibration

## Purpose

`rocket_propulsion.propulsion.burns.calibration` converts repeated hot-fire
measurements into a traceable deterministic correction and residual propulsion
uncertainty. It replaces an unsupported “±x percent” with separately reported:

- persistent engine-to-engine process variance;
- burn-to-burn variance within an engine;
- declared measurement standard uncertainty;
- uncertainty in the estimated consensus correction.

The first implementation calibrates positive dimensionless scales used by the
tabulated-artifact ensemble: thrust, system Isp, duration, and axial efficiency.
The paired extension calibrates two or more of those scales together and turns
same-run covariance evidence into an ensemble-ready Gaussian-copula model.

## Observation contract

Every `HotFireScaleMeasurement` declares engine and run identifiers, an
operating-condition identifier, measured/model scale, measurement standard
uncertainty in the same scale units, and source.

Engine/run identities must be unique. Measurements with different condition
identifiers are rejected rather than pooled. Pressure, mixture ratio, throttle,
thermal state, hardware configuration, or instrumentation changes require a
separate stratum or an explicit regression model.

The `standard_uncertainty` field is a standard uncertainty, not an arbitrary
tolerance or an expanded interval. NIST TN 1297 distinguishes the method used
to evaluate an uncertainty component (Type A or Type B) from whether an effect
is random or systematic; this implementation preserves that distinction.

## Random-effects model

For engine `i` and run `j`, the calibration assumes

```text
y_ij = mu + b_i + e_ij + m_ij
```

where `mu` is the population consensus correction, `b_i` is a zero-mean
persistent engine effect, `e_ij` is common-variance within-engine process
variation, and `m_ij` is independent measurement error with declared standard
uncertainty. The design may be unbalanced and measurement uncertainties may
differ by run. At least two engines and one repeated observation are required.

### Within-engine component

The pooled observed within-engine mean square is

```text
MS_within = sum_i sum_j (y_ij - mean_i)^2 / (N - k)
```

For unequal measurement variances `u_ij^2`, their expected contribution to the
within sum of squares is

```text
sum_i (1 - 1/n_i) sum_j u_ij^2
```

That contribution is reported and subtracted before the nonnegative physical
within-engine variance is formed. Measurement uncertainty is therefore not
silently relabelled as engine variability.

### Between-engine component

Each engine mean receives modeled variance

```text
t_i^2 = sigma_within^2 / n_i + sum_j u_ij^2 / n_i^2
```

The nonnegative between-engine variance `tau^2` is found by the Mandel-Paule
equation

```text
w_i = 1 / (tau^2 + t_i^2)
mu  = sum_i w_i mean_i / sum_i w_i
sum_i w_i (mean_i - mu)^2 = k - 1
```

The solve is bracketed and bisected with retained iterations, residual,
tolerance, numerical variance floor, and boundary status. NIST describes this
method as a computationally simpler approximation to maximum likelihood and
notes that uncertainty intervals may be too narrow with fewer than six groups.
The implementation carries that warning when fewer than six engines are used.

## Using the result

`consensus_scale` is a deterministic correction to the base engine model. It
must be applied before residual sampling. Two residual populations are exposed:

- `CalibrationScope.POPULATION`: a future burn from a randomly selected engine,
  using between-engine plus within-engine process variance;
- `CalibrationScope.SAME_ENGINE`: another burn from an engine whose persistent
  effect is known/corrected, using within-engine process variance only.

`residual_uncertain_parameter()` returns a bounded truncated-normal parameter
centered at one, compatible with `ScenarioConditionalUncertainty`. Its source
contains the complete calibration SHA-256.

```python
calibration = calibrate_hot_fire_scale(
    "thrust_scale",
    measurements,
    source="qualification campaign Q1",
    rationale="rated-condition thrust/model ratios",
)

corrected_thrust = model_thrust * calibration.consensus_scale
residual = calibration.residual_uncertain_parameter(
    scope=CalibrationScope.POPULATION,
)
```

Do not multiply the consensus correction into every residual draw twice. The
residual parameter is centered at one because the deterministic base has
already been corrected.

## Paired multivariate calibration

`HotFireVectorMeasurement` records an ordered vector of scales from one
engine/run and the complete standard measurement-uncertainty covariance matrix
in squared scale units. Every record in a fit must have the same parameter
order and operating-condition identifier. The matrix must be symmetric and
positive semidefinite. Off-diagonal measurement covariance is not assumed zero.

For parameters `p` and `q`, the pooled within-engine cross-covariance is

```text
S_within[p,q] = sum_i sum_j
    (y_ijp - mean_ip) (y_ijq - mean_iq) / (N - k)
```

and the expected contribution of declared same-run measurement covariance is

```text
M_within[p,q] = sum_i (1 - 1/n_i) sum_j C_ij[p,q] / (N - k)
```

The physical within-engine estimate is `S_within - M_within`. Between-engine
cross-covariance uses the sample covariance of paired engine means minus the
average covariance of an engine mean. Marginal diagonal variances remain the
nonnegative Mandel-Paule results from the scalar fits. Consequently the method
is explicitly identified as a hybrid marginal/random-effects plus paired
method-of-moments estimator, not a full multivariate REML fit.

Finite data and measurement subtraction can produce a raw correlation matrix
that is not positive definite. Pairwise coefficients are first limited to the
open physical interval; if necessary, all off-diagonal terms are then shrunk
uniformly toward the identity. Clip count, shrinkage fraction, iteration count,
raw matrix, and effective matrix are retained in
`CorrelationRegularizationEvidence`. No regularization is hidden.

```python
calibration = calibrate_hot_fire_multivariate(
    paired_measurements,
    source="qualification campaign Q2",
    rationale="paired thrust, Isp, and duration residuals",
)

corrected_thrust = model_thrust * calibration.marginals[0].consensus_scale
parameters, correlation = calibration.residual_model(
    scope=CalibrationScope.POPULATION,
)
```

`rocket_propulsion_hot_fire_multivariate_calibration_v1` persistence stores the
complete measurements, covariance matrices, marginal fits, numerical evidence,
references, warnings, and SHA-256. Reopening reruns the calibration and rejects
both hash changes and derived-result tampering. See
`examples/hot_fire_multivariate_calibration.py` for a three-parameter campaign.

## Credibility boundary

This is a variance-component calibration, not a complete Bayesian engine
digital twin. It does not yet model sampling uncertainty of variance components
as a distribution, model-form discrepancy, time drift, censoring, outliers,
autocorrelation, operating-condition response surfaces, nonlinear/tail
dependence beyond a Gaussian copula, manufacturing sub-hierarchies, or
qualification-to-flight transfer. The paired method also assumes measurement
errors from different runs are independent after each run's declared covariance
is applied.

These omissions are reported in warnings and must not be hidden by interpreting
the fitted residual distribution as total uncertainty. NASA/SP-2009-569 and
NASA-STD-7009B support keeping data evidence, parameter uncertainty, model
knowledge, and intended use explicit.

## Primary references

- NIST TN 1297, *Guidelines for Evaluating and Expressing the Uncertainty of
  NIST Measurement Results*:
  <https://www.nist.gov/pml/nist-technical-note-1297>
- NIST/SEMATECH e-Handbook §7.4.4, *Variance Components*:
  <https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm>
- NIST Dataplot, *Consensus Mean*:
  <https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm>
- JCGM 100:2008, *Guide to the Expression of Uncertainty in Measurement*:
  <https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf>
- NASA/SP-2009-569, *Bayesian Inference for NASA Probabilistic Risk and
  Reliability Analysis*:
  <https://ntrs.nasa.gov/citations/20090023159>
- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>

