# Model discrepancy and qualification-to-flight transfer

## Purpose

Repeated hot-fire measurements quantify engine-to-engine and run-to-run
variation inside their test domain. They do not, by themselves, quantify:

- systematic difference between a propulsion model and observed performance;
- the change from qualification hardware/environment to the intended flight
  configuration and application domain.

`rocket_propulsion.propulsion.burns.discrepancy` treats those effects as
separate epistemic components. It does not inflate the hot-fire aleatory sigma
until a desired margin appears.

## Ratio evidence

Every `DiscrepancyRatioObservation` provides an independent case, explicit
applicability identifier, positive ratio, measurement standard uncertainty,
and source.

For model-form evidence:

```text
ratio = observed performance / model prediction
```

For qualification-to-flight evidence:

```text
ratio = flight-representative performance / qualification performance
```

The numerator and denominator must use the same physical scale definition.
Thrust scale must not be mixed with axial efficiency or chamber Isp with total
system Isp. Observations with different `applicability_id` values are refused.

## Random-effects discrepancy model

For validation case `i`:

```text
r_i = mu + d_i + m_i
```

where `mu` is the consensus correction, `d_i` is between-case discrepancy with
variance `tau^2`, and `m_i` has the declared measurement variance `u_i^2`.
Mandel-Paule weights and consensus are

```text
w_i = 1 / (tau^2 + u_i^2)
mu  = sum_i w_i r_i / sum_i w_i
sum_i w_i (r_i - mu)^2 = k - 1
```

The deterministic base artifact receives `mu` once. For a future case inside
the same application domain, the epistemic predictive variance is

```text
sigma_predictive^2 = tau^2 + u_mu^2
u_mu^2             = 1 / sum_i w_i
```

This keeps uncertainty in the consensus bias as well as unexplained
case-to-case discrepancy. Referent measurement uncertainty informs the fit but
is not relabelled as engine variability.

## Layered scale components

`build_layered_calibration_plan()` composes paired hot-fire calibration and
independent discrepancy calibrations. Parameter names retain both physical
target and source class:

```text
thrust_scale::aleatory:hot-fire
thrust_scale::epistemic:model-form:validation-residual
specific_impulse_scale::epistemic:qualification-to-flight:block-transfer
```

The hierarchical ensemble multiplies all components sharing the same prefix.
For example,

```text
k_thrust = product(all thrust_scale components)
k_flow   = k_thrust / product(all specific_impulse_scale components)
```

Time scales segment boundaries; axial-efficiency scales net axial force only.
The common `scale_tabulated_artifact()` implementation is used for both
deterministic consensus corrections and stochastic realizations, preventing
the two paths from drifting apart.

The plan refuses:

- reuse of the hot-fire `evidence_id` as discrepancy evidence;
- a discrepancy outside the declared application domain;
- duplicate qualified component names;
- shared discrepancy evidence without an explicit ordered correlation model;
- an axial-efficiency envelope capable of producing axial force greater than
  delivered thrust.

The default block model preserves paired hot-fire correlation, preserves a
supplied discrepancy correlation, and explicitly assumes zero correlation
between those evidence blocks. This assumption is recorded in the plan.

## Workflow

```python
model_form = calibrate_scale_discrepancy(
    "thrust_scale",
    "validation-residual",
    DiscrepancyKind.MODEL_FORM,
    validation_ratios,
    evidence_id="MF-V1",
    source="independent model validation campaign",
    rationale="rated flight-block application domain",
)

plan = build_layered_calibration_plan(
    hot_fire_calibration,
    (model_form,),
    scenario_id="nominal",
    applicability_id="flight-block-A-rated",
    hot_fire_evidence_id="HF-Q1",
    sample_count=256,
    source="flight block A uncertainty qualification",
    rationale="separate repeatability from model-form uncertainty",
)

corrected = plan.corrected_artifact(base_artifact)
definition = plan.conditional_definition
```

Use `corrected` as the scenario artifact and `definition` as its conditional
uncertainty definition. Do not apply any consensus factor again.

Each discrepancy calibration persists as
`rocket_propulsion_discrepancy_calibration_v1`. Reopening reruns the estimator
and rejects both hash changes and falsified derived results.

## Credibility boundary

The current estimator is a constant discrepancy inside one declared
application domain. It is not a Gaussian-process discrepancy surface and does
not extrapolate across pressure, mixture ratio, throttle, thermal state, lot,
age, or environment. It assumes validation cases are independent unless their
dependence is handled before calibration.

Calibration and discrepancy can be weakly identifiable when they use the same
responses. The evidence-ID refusal reduces accidental double counting but
cannot prove statistical independence. High-consequence use still requires a
reviewed validation plan, application-domain coverage, sensitivity analysis,
and acceptance criteria.

## Primary references

- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
- NASA-HDBK-7009, *NASA Handbook for Models and Simulations*:
  <https://ntrs.nasa.gov/citations/20140002378>
- NASA, *Calibration of Predictor Models Using Multiple Validation
  Experiments*: <https://ntrs.nasa.gov/citations/20150006032>
- NASA/SP-2009-569, *Bayesian Inference for NASA Probabilistic Risk and
  Reliability Analysis*: <https://ntrs.nasa.gov/citations/20090023159>
- NIST Dataplot, *Consensus Mean*:
  <https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm>
