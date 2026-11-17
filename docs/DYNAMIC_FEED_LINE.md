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
absolute residual of this identity using the pressure derivative returned to
the integrator. The residual is `|C * exported_dp_dt - (q - q_e)|` rounded to
binary64. It measures numerical closure for the supplied state and can be zero
at binary64 resolution.

The force balance retains the supplied binary64 values as exact rationals
through the pressure subtraction, signed loss and division by inertance.
The pressure derivative and the rational Jacobian coefficients likewise round
after their complete expression. An output that would become infinite or
lose a nonzero value to rounding is refused with its field name. The modal
square roots use binary64, with the damping product and quotient balanced
before export.

In the coupled single- and dual-feed models, negative line states remain
available to bracket the zero-flow cutoff. The regulated tank then uses the
same signed flow in its mass and ullage-work equations. Such a result exposes
`reverse_flow_event_extension=True`; it is a mathematical event-location
continuation, not sustained reverse hardware operation. The propagator must
terminate physical propagation at the zero-flow event. Direct regulated-feed
calls reject negative flow unless `for_event_location=True` is explicitly set,
and that extension requires zero static line resistance. The dual-feed vehicle
closure uses the actual returned tank derivatives. See
[the conservation review](COUPLED_FEED_REVIEW.md) for reproduced failures and
an integrated cutoff-root check.

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

## Verification against independent controls

- A 1 Pa head opposed by a quadratic loss with `R = 100` and binary64
  `q = 0.1` previously yielded a zero flow derivative. A 1,200-digit Decimal
  calculation gives approximately `−1.11022e−22 kg/s²` at `L = 1e6`.
  Reversing the head and flow reverses the sign. The public motor coupling
  reproduces the same failure mechanism with a 2 MPa tank and 1.8 MPa
  manifold. These are checks on the supplied numbers at near balance.
- With `R = L = 1e308` and `q = 0.5`, the flow Jacobian is `−1 s⁻¹`.
  Doubling the resistance before forming the quotient previously overflowed.
  Separate cases reject nonzero friction or state rates that round to zero.
- For `R = 0`, constant `p_up`, and a step engine demand `Q` from an initial
  state `q = 0, p_m = p_up`, the closed form is
  `q = Q(1 − cos(ωt))` and `p_m = p_up − Q sin(ωt)/(Cω)`.
  With `p_up = 2 MPa`, `Q = 0.2 kg/s` and `ω = 10 rad/s`, the public
  pressure event stops integration at 1 MPa after 0.05235988 s, with flow
  0.02679492 kg/s. Reversing the event sign causes the solver to miss that
  cutoff and the test to fail on the returned integration status.

The targeted regression covers the public line, regulated feed, active
single-feed and dual-feed motor coupling, and feed API cases. Isolated
mutations reintroduce premature friction rounding, an incorrect Jacobian
coefficient, a false zero mass residual and the reversed pressure cutoff.
Each selected mutant must fail its assertion, with a passing original test
both before and after injection.

## Experimental check on the model's usable bandwidth

Figure 33 of [NASA CR-124308](https://ntrs.nasa.gov/citations/19730020456)
plots measured terminal pressure per imposed volumetric flow for a water-filled
aluminum line with a linear terminal resistor. The line was 27 ft long, with
a 3 in bore and a 0.050 in wall. The side-branch pulser was just above the
terminal. The report identifies first and second longitudinal modes.

Five open-circle measurements, read from the scan rather than the theory
curve, show the two elevated regions and intervening valley. The coordinates
use a 2× raster of PDF page 95, printed page 79. The horizontal log axis uses
10 and 100 Hz at pixels 705 and 1145. The vertical log axis uses `10⁵` and
`10⁶ lbf·s/ft⁵` at pixels 848 and 410.

| Marker | Pixel (x, y) | Approx. frequency (Hz) | Approx. amplitude (lbf·s/ft⁵) |
|---|---:|---:|---:|
| Before first resonance | 704, 919 | 9.95 | 68,849 |
| First resonance region | 891, 385 | 26.47 | 1,140,454 |
| Valley | 1045, 1271 | 59.26 | 10,820 |
| Second resonance region | 1102, 528 | 79.85 | 537,766 |
| After second resonance | 1145, 971 | 100.00 | 52,382 |

Allowing ±8 pixels on each reading leaves each elevated region more than
96× and 45× above the valley, respectively. These scan-reading intervals
carry no claim about the experiment's instrument uncertainty or the precise
resonance frequencies. The PDF SHA-256 is
`e4c8915cee8829af4d98584d68f394b7a4a573e07abc64c2d911332f8bc25a9a`.

For the lumped model at fixed upstream pressure, let
`r = 2R|q₀|` be the local line resistance and `g` a static incremental
terminal conductance. With the imposed flow disturbance entering at the
manifold, the linearized pressure transfer function is

```text
P(s) / D(s) = −(Ls + r) / [(Cs + g)(Ls + r) + 1].
```

Writing `x = ω²` gives squared magnitude `(a + b x)/(d₀ + d₁ x + d₂ x²)`,
where `a = r²`, `b = L²`, `d₀ = (1 + rg)²` and `d₂ = (LC)²`.
The sign of its slope comes from
`b d₀ − a d₁ − 2a d₂ x − b d₂ x²`, a strictly decreasing function for
`x > 0`. Under these boundary assumptions, the transfer can have at most one
local peak and cannot have a valley followed by another rise. Changing the
constant coefficients cannot reproduce the measured shape across both
longitudinal modes.

This test constrains the scope of the surrogate. A first-mode fit needs
configuration-specific gain and phase measurements and a stated frequency
band. An integrated motor and feed acceptance test also needs synchronized
pressure and flow, the valve timing and inlet boundary conditions, and the
experimental uncertainty. Pressure peaks beyond the first mode require
separate measurement or a distributed model validated in the applicable
configuration. [NASA NESC Technical Bulletin 22-03](https://ntrs.nasa.gov/citations/20220006583)

A separate comparison used 15 measured oscillation periods from an acrylic
water line with a trapped air pocket and a closed end. With the initial
20 ft liquid length held fixed, inertance and equilibrium gas compliance
were calculated from the rig geometry and pressure, without fitting them
to the periods. The library's period errors ranged from −13.62% to +9.03%.
This comparison tests the local period only. The full valve-opening
transient, peak pressure and damping remain outside its scope. The
calculation and measured values are in
`output/assurance/entrapped_air_period_benchmark_20261110`.
[Nahm Ho Lee, Georgia Tech, Table 5.1](https://repository.gatech.edu/bitstreams/ec90ac9c-69b7-4b2b-b809-be7e2bd3e981/download)

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

