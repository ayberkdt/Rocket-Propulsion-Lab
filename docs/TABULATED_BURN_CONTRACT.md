# Tabulated burn exchange contract

## Status and ownership

`rocket_propulsion_tabulated_burn_v1` is an isolated Rocket Propulsion Lab
prototype for a future Sidera transient-maneuver capability. No Sidera source
file or schema is modified by this implementation.

Rocket Propulsion Lab owns thrust, axial thrust after hardware cant, total tank
drain, optional per-tank drains, command/response physics, provenance and
integral closure. Sidera will continue to own epoch, maneuver ordering,
direction/frame interpretation, state propagation, gravity, events and orbital
results.

## Why the contract stores segments

A point-only table cannot represent an instantaneous throttle jump without
duplicate timestamps or unintended interpolation. This contract therefore
stores contiguous segments with independent start and end states:

```text
[t0, t1]: F_start -> F_end, mdot_start -> mdot_end
[t1, t2]: independent F_start -> F_end, independent mdot_start -> mdot_end
```

The first segment may end at full thrust while the next begins at zero thrust.
The discontinuity is exact and the consumer must split integration at `t1`.

## Required physical fields

Every segment carries SI values for:

- relative start and end time;
- delivered scalar thrust;
- axial thrust after declared cluster cant efficiency;
- total tank-depleting mass flow;
- optional named tank-flow components;
- a phase label.

The artifact also carries the source identifier, source SHA-256, explicit
interpolation policy, zero or more stable engineering `reference_ids`, warnings,
derived duration/impulse/mass/centroid summary, and its own deterministic
content hash. Changing the citations changes the artifact hash.

Reference identifiers resolve through the reviewed registry in
`rocket_propulsion.propulsion.burns.references`; they are not free-form claims.

Per-tank flows, when present, must sum to total mass flow at both segment ends.
Axial thrust cannot exceed delivered thrust. Times are contiguous, strictly
advancing within a segment and begin at zero relative time.

## Command response

The implemented transient law is

```text
tau dr/dt + r = u(t)
```

where `u(t)` is linear inside a command segment and `r(t)` is realized
throttle. The exponential solution, throttle exposure and first time moment
are analytic. A zero time constant exactly recovers the commanded linear
schedule.

Exponential histories are sampled for exchange only. Sampling is refined until
the tabulated delivered impulse, axial impulse, propellant mass and thrust
centroid agree with the analytic response integrals at the requested tolerance.
The evidence records the final step, refinement count and every residual.

## Imported traces

CSV replay requires:

```text
time_s,delivered_thrust_n,total_mass_flow_kg_s
```

Optional columns are `axial_thrust_n`, `phase`, and any number of
`stream:<tank-id>` drains. Time must begin at zero and increase strictly.
`linear` and zero-order `hold` interpolation are explicit choices. Input text
is SHA-256 hashed and samples are not smoothed.

## Pressure-fed blowdown producer

The isolated producer can generate this contract from a polytropic tank model
and either a bounded pressure correlation or a measured pressure-performance
table. Its mass-domain integration is refined against independent Simpson
integrals for mass, impulse, duration and centroid. See
[`BLOWDOWN_MODEL.md`](BLOWDOWN_MODEL.md). This does not change the current
Sidera consumer status.

## Uncertainty and discrete-scenario producers

Continuous propulsion uncertainty retains one complete artifact per sampled
realization. Discrete ignition/engine-out, common-cause, delay, early-cutoff and
performance-degradation scenarios likewise retain one probability-labelled
artifact per mutually exclusive outcome. Summary percentiles are never used to
reconstruct a representative burn.

A future Sidera ensemble consumer should propagate every retained artifact,
then apply scenario weights and conditional continuous-sample weights to the
trajectory outputs. It must not apply scenario probability once per inner
continuous sample or silently drop zero-thrust outcomes. See
[`PROPULSION_UNCERTAINTY.md`](PROPULSION_UNCERTAINTY.md) and
[`PROPULSION_SCENARIOS.md`](PROPULSION_SCENARIOS.md).

The implemented hierarchical producer formalizes that composition and emits
the joint weight beside every retained artifact. Its exact weighting,
conditional scaling and Sidera ownership rules are documented in
[`HIERARCHICAL_PROPULSION_ENSEMBLE.md`](HIERARCHICAL_PROPULSION_ENSEMBLE.md).

The dependency-free command-line path is:

```powershell
rocket-propulsion burn tabulate-trace examples/solid_thrust_trace.csv `
  --source-id static-fire-001 --interpolation linear `
  --output solid-profile.json
```

## Future Sidera consumer rules

A Sidera implementation may consume this contract only if it:

1. verifies schema and artifact hash;
2. splits propagation at every segment boundary;
3. interpolates thrust and mass flow independently;
4. uses total tank drain for spacecraft mass depletion;
5. uses axial thrust, not pre-cant delivered thrust, as the net force magnitude;
6. applies direction/frame logic inside Sidera, not this artifact;
7. reports the source artifact hash with trajectory results;
8. rejects unsupported schemas, negative mass, non-closing streams and
   overlapping maneuvers;
9. never silently reduces a transient artifact to a constant burn.

The weighted execution wrapper is specified in
[`SIDERA_ENSEMBLE_HANDOFF.md`](SIDERA_ENSEMBLE_HANDOFF.md).

The current production Sidera adapter remains unchanged and continues to
accept only exact constant burns or an explicitly requested equivalent-constant
reduction.
