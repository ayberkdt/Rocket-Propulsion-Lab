# Discrete propulsion scenario contract

## Purpose

Continuous parameter uncertainty and discrete failure scenarios answer
different questions. A thrust calibration interval can be sampled as a
continuous uncertainty; an ignition failure or common-cause engine-out is a
distinct system outcome. This module keeps those concepts separate.

`rocket_propulsion.propulsion.burns.scenarios` builds an exhaustive set of
mutually exclusive propulsion histories. Every history is target-neutral and
can later be propagated by Sidera from the same orbital initial state.

## Probability policy

The library does not calculate component failure probability from engine count
and does not assume failures are independent. Every scenario probability must
come from an analyst-declared source and rationale. Scenario probabilities must
sum to one within a strict numerical tolerance; missing probability mass is an
error rather than an implicit “other” outcome.

This is particularly important for clustered engines. NASA reliability work
identifies engine count, engine-out design, start/cutoff transients, duration,
restart and health-management behavior as relevant reliability drivers. NASA
common-cause guidance also emphasizes that dependent failures can defeat
redundancy and that applicable aerospace data may be limited.

## Supported events

`ScenarioEventType` currently supports:

- `engine-disabled`: ignition failure or declared engine-out;
- `ignition-delay`: additional delay before the engine command begins;
- `early-cutoff`: command is replaced by zero after the specified relative time;
- `thrust-scale`: delivered/ideal thrust scaling with stream closure;
- `specific-impulse-scale`: system-Isp scaling with recomputed tank flow;
- `response-time-scale`: scaling of the declared first-order actuator response.

A single event may name several engines. This represents one explicitly
dependent/common-cause outcome; the implementation does not multiply separate
engine probabilities.

The event source and rationale are mandatory. Repeated event types for the same
engine, unknown engine identifiers, non-physical scale factors, and additional
modifications applied to a disabled engine are rejected.

## Physical realization

For every scenario the builder:

1. starts from the same nominal cluster members and commands;
2. applies declared engine disablement and timing changes;
3. recomputes thrust, system Isp, total tank flow and named stream flows for
   performance changes;
4. resolves every throttle discontinuity;
5. runs the cluster convergence check against analytic response integrals;
6. stores a complete `rocket_propulsion_tabulated_burn_v1` artifact.

An all-engines-disabled scenario produces an exact zero-thrust/zero-flow
artifact over the nominal command duration. It is retained, not discarded.

## Probability-weighted outputs

The ensemble reports expected value, minimum, P05, P50, P95 and maximum for:

- duration;
- delivered impulse;
- axial impulse;
- consumed propellant;
- thrust-centroid time.

It also reports exact zero-thrust probability. These are propulsion
consequences only. Orbit-state, position, velocity, element or targeting
consequences must be computed by propagating each scenario in Sidera.

## Composition with continuous uncertainty

The intended higher-fidelity workflow is hierarchical:

```text
discrete scenario
    -> continuous propulsion parameter draws conditional on that scenario
        -> one thrust/mass-flow artifact per draw
            -> Sidera trajectory propagation
                -> weighted orbit-dispersion/risk result
```

Scenario probability must be applied once. Continuous samples inside a
scenario describe the conditional distribution and must not each inherit the
full scenario probability.

This workflow is implemented by
[`build_hierarchical_propulsion_ensemble`](HIERARCHICAL_PROPULSION_ENSEMBLE.md).
For `N_s` conditional samples, each retained history receives weight
`P(s) / N_s`; deterministic scenarios use one sample. The builder checks joint
probability closure and reports each scenario's contribution to the expected
propulsion consequence.

## Example

```python
scenario_set = PropulsionScenarioSet(
    scenarios=(
        PropulsionScenario(
            "nominal",
            0.98,
            (),
            source="approved event tree",
            rationale="successful two-engine maneuver",
        ),
        PropulsionScenario(
            "single-engine-out",
            0.015,
            (
                ScenarioEvent(
                    ScenarioEventType.ENGINE_DISABLED,
                    ("engine-a",),
                    None,
                    source="approved reliability assessment",
                    rationale="engine-a fails to ignite",
                ),
            ),
            source="approved event tree",
            rationale="one-engine contingency",
        ),
        PropulsionScenario(
            "common-cause-no-ignition",
            0.005,
            (
                ScenarioEvent(
                    ScenarioEventType.ENGINE_DISABLED,
                    ("engine-a", "engine-b"),
                    None,
                    source="approved common-cause assessment",
                    rationale="shared initiation system unavailable",
                ),
            ),
            source="approved event tree",
            rationale="dependent loss of both engines",
        ),
    ),
    objective="quantify maneuver propulsion consequences",
)

ensemble = build_propulsion_scenario_ensemble(members, commands, scenario_set)
```

The numerical probabilities above are API examples only and must not be used as
generic rocket-engine reliability values.

## Primary references

- NASA/SP-2011-3421, *Probabilistic Risk Assessment Procedures Guide for NASA
  Managers and Practitioners*:
  <https://ntrs.nasa.gov/citations/20120001369>
- NASA, *Key Reliability Drivers of Liquid Propulsion Engines and a Reliability
  Model for Sensitivity Analysis*:
  <https://ntrs.nasa.gov/citations/20050207429>
- NASA, *Common Cause Failure Modeling in Space Launch Vehicles*:
  <https://ntrs.nasa.gov/citations/20160007073>
- NASA, *Transient Mathematical Modeling for Liquid Rocket Engine Systems*:
  <https://ntrs.nasa.gov/citations/20040000363>
