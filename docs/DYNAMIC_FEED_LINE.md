# Dynamic liquid feed-line model

`rocket_propulsion_dynamic_feed_line_v1` is a two-state, mission-level
surrogate for the first hydraulic mode of a liquid-propellant feed line. It is
intended for simultaneous integration with an orbit propagator when a static
`K mdot²` pressure loss is too weak, but a distributed water-hammer solver is
too expensive or unsupported by the available qualification data.

It does not modify or import Sidera. The framework-neutral output supplies the
state derivatives, event roots, conservation evidence, and analytic partials a
future native force-model wrapper needs.

## State and equations

The propagated state order is fixed by the model manifest:

1. `mass_flow_kg_s` — liquid flow from the tank into the line, `q`;
2. `manifold_pressure_pa` — engine-side line/manifold pressure, `p_m`.

For upstream tank pressure `p_up`, engine demand `q_e`, inertance `L`,
compliance `C`, and quadratic resistance `R`,

```text
L dq/dt   = p_up - p_m - R q |q|
C dp_m/dt = q - q_e
```

The compliance state stores incremental liquid mass. Its rate is

```text
dm_stored/dt = q - q_e.
```

The tank liquid rate is `-q`, while the external vehicle mass rate produced by
the engine is `-q_e`. Therefore

```text
dm_tank/dt + dm_stored/dt - dm_vehicle/dt = 0.
```

`DynamicFeedLineEvaluation.external_mass_closure_error_kg_s` reports the
absolute residual of this identity at every evaluation.

The uncoupled local modal diagnostics are

```text
omega_n = sqrt(1 / (L C))
zeta    = R |q| sqrt(C / L).
```

They describe the linearization about the current operating point. They are
not a substitute for a stability assessment of the complete engine/feed
system, whose engine demand gradient changes the composed Jacobian.

## Engine and pressurization coupling

`evaluate_dynamic_feed_propulsion()` performs one simultaneous evaluation:

1. manifold pressure is supplied to a bounded
   `RectilinearPerformanceSurface`;
2. the `PropagatorPropulsionBridge` returns force and engine mass-flow demand;
3. line flow drains the liquid tank and expands the ullage in
   `RegulatedFeedSystem`;
4. engine demand changes the compliant line storage;
5. the engine pressure gradient is composed into the two-state line Jacobian.

There is no algebraic pressure/flow iteration in this formulation because
manifold pressure and line flow are propagated states.

The regulated feed configuration must set
`feed_line_resistance_pa_s2_kg2 = 0`. Applying a static line loss as well as the
dynamic line resistance is rejected because it would count the same loss
twice.

While the engine is off or inhibited, engine demand and its pressure gradient
are exactly zero, but pressurization and line relaxation may continue. This is
important across ignition and shutdown boundaries.

## Jacobians

The local line derivatives are

```text
d(qdot)/dq     = -2 R |q| / L
d(qdot)/dp_m   = -1 / L
d(qdot)/dp_up  =  1 / L
d(pdot)/dq     =  1 / C
d(pdot)/dq_e   = -1 / C.
```

At `q = 0`, the implementation marks the flow coordinate as non-smooth. The
composed manifold-pressure column additionally contains
`-(dq_e/dp_m) / C`. Upstream tank-pressure derivatives are chained through the
regulated-feed pressure partials. Acceleration derivatives with respect to the
line states are returned in the resolved thrust frame.

The automated tests compare these terms with central finite differences of the
fully composed model rather than checking only copied formulas.

## Event surfaces

Each event function is positive inside the accepted domain and zero at its
boundary:

| Event | Root function | Action |
|---|---|---|
| reverse flow | `q` | inhibit the burn before sustained reverse operation |
| minimum manifold pressure | `p_m - p_min` | inhibit the burn |
| maximum manifold pressure | `p_max - p_m` | inhibit the burn |
| maximum absolute flow | `q_max - |q|` | inhibit the burn |

A consuming propagator should also register the regulated-feed propellant
reserve, bottle-pressure-margin, and tank-overpressure roots. When a dynamic
line is active, its minimum-manifold-pressure event replaces the static
regulated-feed injector-pressure event: with static resistance disabled, the
latter is only the tank pressure and no longer represents engine inlet
pressure.

The burn start, profile-boundary, stop, dry-mass, tank-reserve, and performance
map-domain roots remain owned by `PropagatorPropulsionBridge`.

## Propagator state ownership

A native Sidera or Orekit-style consumer should integrate, in one state vector:

- position and velocity;
- total spacecraft mass;
- any explicitly bound propellant tank masses;
- the six regulated-feed states;
- line flow and manifold pressure.

The consumer owns epoch, frame transformation, attitude/guidance, root
location, reset/inhibition policy, numerical tolerances, and trajectory output.
The propulsion package owns force magnitude, mass/inventory derivatives,
validity domains, and propulsion-state Jacobians.

## Fidelity and calibration limits

This model represents one effective inertance/compliance mode. It omits:

- distributed acoustic modes and wave travel time;
- water hammer and valve slam;
- priming fronts, gas ingestion, cavitation, and two-phase flow;
- line-wall structural modes and fluid-structure interaction;
- temperature-dependent real-fluid properties;
- injector/chamber combustion instability coupling.

`L`, `C`, and `R` are hardware-specific parameters. They must be estimated from
geometry plus fluid properties and correlated against component, cold-flow,
hot-fire, or a validated higher-order nodal model over the declared operating
domain. The implementation must not be interpreted as flight-qualified merely
because its numerical residuals close.

## Primary references

- NASA, *Liquid Rocket Propellant Feedline Dynamics*:
  <https://ntrs.nasa.gov/citations/19740028545>
- NASA NESC Technical Bulletin 22-03, transient pressure and water-hammer
  screening: <https://ntrs.nasa.gov/citations/20220006583>
- NASA, *Nodal Modeling of Liquid Propellant Feed and Pressurization System*:
  <https://ntrs.nasa.gov/citations/20240003493>
- NASA SP-8112, *Pressurization Systems for Liquid Rockets*:
  <https://ntrs.nasa.gov/citations/19760015212>
- Orekit 13.1.5 `EventDetector` and `PropulsionModel` APIs:
  <https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html>
  and
  <https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html>

