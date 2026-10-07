# Hierarchical propulsion ensemble

## Purpose

`rocket_propulsion.propulsion.burns.ensemble` combines two different forms of
propulsion uncertainty without confusing their probabilities:

1. mutually exclusive discrete outcomes such as nominal operation,
   single-engine loss, early cutoff, or common-cause loss;
2. continuous performance variation conditional on each outcome.

Every retained realization is a complete
`rocket_propulsion_tabulated_burn_v1` thrust/mass-flow artifact. The result is
therefore useful as a standalone propulsion-risk product and is ready for a
future Sidera trajectory ensemble consumer. It does not modify Sidera.

## Probability measure

For scenario `s`, let its declared probability be `P(s)` and let its
conditional design contain `N_s` equally weighted samples. Realization `i`
has

```text
P(i | s) = 1 / N_s
P(s, i)  = P(s) / N_s
```

The implementation verifies that the joint weights sum to one within
`1e-12`. Scenario probability is applied exactly once. A deterministic
scenario has no continuous parameters, uses one sample, and therefore retains
its full scenario probability. Zero-thrust outcomes remain in the measure.

Conditional definitions can differ between scenarios. For example, an
engine-out branch may use qualification evidence for the surviving engine,
while a common-cause no-ignition branch remains deterministic. Each definition
requires its own source and rationale.

Conditional scales may be derived from repeated test data using the implemented
[`HOT_FIRE_CALIBRATION.md`](HOT_FIRE_CALIBRATION.md) workflow. Apply the
calibrated consensus correction to the deterministic base first; the returned
residual parameter remains centered at one and is sampled conditionally here.
Independent model-form and qualification-to-flight evidence can then be added
through the layered workflow in
[`MODEL_DISCREPANCY_TRANSFER.md`](MODEL_DISCREPANCY_TRANSFER.md). Component
names retain their aleatory/epistemic source while identical physical prefixes
are multiplied into one artifact scale.

## Conditional artifact scales

The first implementation supports positive, bounded scales:

- `thrust_scale`: multiplies delivered and axial thrust;
- `specific_impulse_scale`: changes flow so system Isp changes by the declared
  scale;
- `time_scale`: dilates every segment boundary and therefore burn duration;
- `axial_efficiency_scale`: changes axial thrust without changing delivered
  thrust.

Each name may be unqualified, for backward compatibility, or use
`physical_scale::class:component`. Multiple qualified components may target the
same scale; their sampled values are multiplied. Their complete names must
remain unique and the correlation order covers components, not only physical
targets.

For a draw with thrust scale `k_F` and system-Isp scale `k_I`, total and named
tank flows scale by `k_F / k_I`. This preserves
`Isp_system = F_delivered / (g0 * mdot_tank)`. Time scaling acts on the entire
history, so impulse and propellant consumption change consistently. A proposed
axial-efficiency distribution is rejected if any allowed draw could make
axial thrust exceed delivered thrust.

The model is deliberately scale-based. It does not claim that these four
variables reproduce feed-system dynamics, combustion instability, thermal
soakback, or regulator hysteresis. Constant in-domain model/flight discrepancy
can be represented by the implemented evidence layer, but a scale is not a
substitute for a state-dependent discrepancy surface.

## Evidence and reproducibility

Sampling may use Latin hypercube or Monte Carlo. Seeds for every scenario and
for primary/audit channels are derived independently from the base seed and
scenario identifier. The primary samples become retained realizations. The
audit samples are used only to compare delivered-impulse mean and sample
standard deviation; they do not receive probability mass and are never mixed
into the reported distribution.

The replicate comparison is a diagnostic, not proof of tail convergence. A
failure adds a warning and should lead to a larger conditional sample design.
Every artifact hash includes the base artifact, scenario, complete conditional
definition, sampled inputs, seed, and method. Changing the engineering source
or rationale therefore changes provenance even when the sampled numbers happen
to match.

## Outputs

The result retains:

- all conditional definitions;
- every joint realization and its scenario, conditional, and joint weights;
- weighted mean, standard deviation, extrema, P05, P50, and P95 for duration,
  delivered impulse, axial impulse, propellant consumption, and thrust
  centroid;
- each scenario's conditional mean and contribution to the joint expectation;
- exact zero-thrust probability;
- independent-replicate evidence and warnings.

For every reported metric it also applies the law of total variance:

```text
Var(X) = E_s[Var(X | s)] + Var_s(E[X | s])
```

The first term quantifies continuous performance variation inside scenarios;
the second quantifies separation between nominal, degraded, and failed
outcomes. The result reports both fractions and a numerical closure residual.
This prevents a small within-nominal scatter from hiding a dominant discrete
failure contribution.

Weighted percentiles use the discrete empirical cumulative distribution of the
retained joint measure. They are not interpolated quantiles.

## Persistence and integrity

`rocket_propulsion_hierarchical_ensemble_v1` stores the complete conditional
definitions, every joint realization and nested tabulated artifact, statistics,
variance decomposition, scenario contributions, replicate evidence, references
and warnings. The document has its own SHA-256 identity; each nested burn keeps
its independent artifact hash.

The reader does more than compare the outer hash. It reopens and verifies every
nested artifact, then reproduces probability closure, zero-thrust probability,
weighted statistics, total-variance decomposition, scenario contributions and
sampling-evidence flags. A document with a newly calculated hash but falsified
derived statistics is therefore still rejected.

```python
text = hierarchical_ensemble_to_json(result)
reopened = hierarchical_ensemble_from_json(text)
assert hierarchical_ensemble_hash(reopened) == hierarchical_ensemble_hash(result)
```

## Sidera boundary

A future Sidera ensemble adapter should propagate each retained artifact from
the requested maneuver epoch and spacecraft state, then apply the artifact's
joint probability to trajectory outputs. Sidera owns direction, frame,
attitude, guidance, epoch dispersion, gravity, events, and orbital results.
Rocket Propulsion Lab owns thrust, axial hardware loss, tank drain, and
propulsion provenance.

The consumer must not:

- multiply a realization by the scenario probability a second time;
- replace all histories with a percentile or mean burn;
- discard failed/zero-thrust outcomes;
- interpret axial scalar thrust as a complete inertial force vector;
- mix pointing or epoch uncertainty into propulsion calibration parameters.

## Primary references

- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
- NASA/SP-2011-3421, *Probabilistic Risk Assessment Procedures Guide for NASA
  Managers and Practitioners*:
  <https://ntrs.nasa.gov/citations/20120001369>
- NASA, *Key Reliability Drivers of Liquid Propulsion Engines and a Reliability
  Model for Sensitivity Analysis*:
  <https://ntrs.nasa.gov/citations/20050207429>
- NASA, *Common Cause Failure Modeling in Space Launch Vehicles*:
  <https://ntrs.nasa.gov/citations/20160007073>

