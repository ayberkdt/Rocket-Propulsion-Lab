# Propulsion uncertainty and ensemble contract

## Purpose

The uncertainty layer produces an ensemble of physically closed propulsion
histories for later trajectory propagation. It does not turn an uncertainty
percentage directly into orbit dispersion and it does not own reference-frame,
attitude, or guidance uncertainty. Those effects belong in Sidera when it
propagates each retained thrust/mass-flow realization.

Implementation:
`rocket_propulsion.propulsion.burns.uncertainty`.

## Credibility rules

Every uncertain input must declare:

- a unique physical parameter name;
- nominal, lower, and upper values;
- marginal distribution;
- aleatory or epistemic classification;
- an engineering source;
- a written rationale;
- standard deviation when a truncated-normal model is selected.

This follows NASA-STD-7009B's emphasis on model credibility, validation and
uncertainty qualification. A bare “plus or minus five percent” cannot enter
the propulsion ensemble without a source and rationale.

Repeated hot-fire data can replace analyst-selected bounds through
[`HOT_FIRE_CALIBRATION.md`](HOT_FIRE_CALIBRATION.md). That layer separately
estimates engine-to-engine and within-engine process variance, removes declared
measurement uncertainty from process scatter, and emits a hash-linked bounded
residual parameter. Paired hot-fire vectors additionally estimate
cross-parameter covariance, account for full same-run measurement covariance,
and emit the ordered `CorrelationModel` together with all marginal residuals.

`UncertaintyClass.ALEATORY` represents modeled realization-to-realization
variability. `EPISTEMIC` represents incomplete knowledge that may be reduced by
additional analysis or test data. The result reports their counts separately;
it does not collapse the distinction into one label.

## Supported distributions

- bounded uniform;
- bounded triangular, with the nominal value as the mode;
- bounded/truncated normal, with declared sigma;
- bounded log-uniform for positive scale-spanning quantities.

All distributions are transformed from open-interval probabilities and are
guaranteed to stay inside the declared bounds. The implementation deliberately
does not provide an unbounded normal distribution for quantities such as
pressure, thrust or specific impulse.

## Correlation

`CorrelationModel` applies an ordered Gaussian copula. The correlation matrix
must be finite, symmetric, positive definite and have a unit diagonal. A
correlation source and rationale are mandatory.

For independent marginals, Latin-hypercube strata are preserved exactly. With
a copula, the dependency structure takes priority; the transformed final ranks
are not claimed to remain an exact Latin hypercube. This limitation is reported
in the function docstring rather than hidden.

## Blowdown realization mapping

`simulate_blowdown_uncertainty` currently accepts uncertainty in:

- initial tank pressure;
- initial ullage volume;
- initial propellant inventory;
- propellant density;
- polytropic exponent;
- rated delivered thrust;
- rated system specific impulse;
- minimum usable supply pressure;
- thrust, Isp and chamber-pressure correlation exponents;
- reserve mass;
- cant efficiency.

The declared uncertain nominal must match the deterministic model input. This
prevents a distribution centered on one engine definition from being silently
applied to another.

For every draw the solver rebuilds the tank, rated operating point and pressure
law, then reruns the deterministic blowdown convergence gate. Total tank flow
and all stream flows are recomputed from delivered thrust and system Isp. A
failed realization aborts the study with its sample index; invalid samples are
never dropped to make the output distribution look better.

Each `BlowdownUncertaintySample` retains its complete
`rocket_propulsion_tabulated_burn_v1` artifact. Consequently, a later Sidera
ensemble adapter can propagate the exact sampled thrust and mass history rather
than reconstructing it from summary percentiles.

## Statistics and sensitivity

The result reports, for duration, impulse, consumed mass, thrust centroid and
final supply pressure:

- sample mean and standard deviation;
- standard error and normal-approximation 95% interval for the mean;
- minimum, P05, P50, P95 and maximum.

Pearson input/output correlations are ranked independently for every output.
They are explicitly labelled screening sensitivities: they can miss nonlinear
interactions and are not Sobol indices or evidence of causality.

The first-half/full-ensemble change in impulse mean and standard deviation is
recorded as a sampling-stability diagnostic. Failing this tolerance adds a
warning and indicates that the ensemble should be enlarged; it is not silently
reported as converged.

## Example

```python
from rocket_propulsion.propulsion.burns import (
    DistributionKind,
    UncertainParameter,
    UncertaintyClass,
    simulate_blowdown_uncertainty,
)

thrust_uncertainty = UncertainParameter(
    name="rated_delivered_thrust_n",
    nominal=100.0,
    lower=97.0,
    upper=103.0,
    distribution=DistributionKind.TRUNCATED_NORMAL,
    uncertainty_class=UncertaintyClass.EPISTEMIC,
    standard_deviation=1.0,
    source="hot-fire campaign A, calibrated load-cell result",
    rationale="95% engineering interval including calibration uncertainty",
)

ensemble = simulate_blowdown_uncertainty(
    tank,
    rated_point,
    pressure_law,
    (thrust_uncertainty,),
    sample_count=256,
    seed=42,
)
```

The names `tank`, `rated_point`, and `pressure_law` above refer to the
deterministic objects documented in [`BLOWDOWN_MODEL.md`](BLOWDOWN_MODEL.md).

## Current fidelity boundary

This module quantifies continuous parameter uncertainty; it is not itself a
failure-probability or reliability model. Discrete ignition failure, engine-out,
early-cutoff and common-cause outcomes are represented separately by the
implemented [`PROPULSION_SCENARIOS.md`](PROPULSION_SCENARIOS.md) contract.
Pointing, attitude-estimation and maneuver-epoch uncertainties must still be
composed in Sidera because their consequences depend on trajectory state and
frame.

Future higher-fidelity additions should include model-form discrepancy,
operating-condition regression, multivariate REML or Bayesian calibration when
data volume justifies it, Sobol/global sensitivity, and sequential sample-size
convergence studies.

Discrete outcomes and continuous scales can now be composed by the implemented
[`HIERARCHICAL_PROPULSION_ENSEMBLE.md`](HIERARCHICAL_PROPULSION_ENSEMBLE.md)
contract. That layer supports different conditional definitions per scenario,
retains the complete artifact for every joint realization, and keeps scenario
probability from being applied twice.

## Primary references

- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
- NASA/SP-2011-3421, *Probabilistic Risk Assessment Procedures Guide for NASA
  Managers and Practitioners*:
  <https://ntrs.nasa.gov/citations/20120001369>
- NASA/SP-2009-569, *Bayesian Inference for NASA Probabilistic Risk and
  Reliability Analysis*: <https://ntrs.nasa.gov/citations/20090023159>
- NASA SP-8112, *Pressurization Systems for Liquid Rockets*:
  <https://ntrs.nasa.gov/citations/19760015212>
- NASA CR-131400, *Reliability Model of a Monopropellant Auxiliary Propulsion
  System*: <https://ntrs.nasa.gov/citations/19730012094>
