# Propagator propulsion bridge

`PropagatorPropulsionBridge` is the numerical-propagation boundary for the
isolated propulsion core. It turns a verified
`rocket_propulsion_tabulated_burn_v1` artifact into simultaneous force,
vehicle-mass, and optional tank-state derivatives without importing Sidera or
Orekit.

An optional `RectilinearPerformanceSurface` makes the same transient profile
state-dependent on supply pressure, ambient pressure, throttle, mixture ratio,
power or other explicitly named conditions. See
[`OPERATING_CONDITION_SURFACES.md`](OPERATING_CONDITION_SURFACES.md).
`RegulatedFeedSystem` can produce the supply-pressure condition from six
propagated feed states and close the pressure/engine-flow algebraic loop; see
[`REGULATED_FEED_DYNAMICS.md`](REGULATED_FEED_DYNAMICS.md).

## Why this is a separate layer

The tabulated artifact owns propulsion truth: axial force, total tank drain,
named stream drain, discontinuities, interpolation, provenance, and citations.
The astrodynamics framework owns epoch conversion, attitude, frame transforms,
position/velocity integration, and event handling. Keeping this boundary narrow
prevents an engine model from quietly becoming a second orbit propagator.

The contract follows the responsibilities exposed by Orekit's
[`PropulsionModel`](https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html):
acceleration, mass derivative, parameter drivers, initialization-time event
registration, and field/variational compatibility. It extends the practical
bookkeeping beyond a single mass scalar by returning named tank-state drains.

## Evaluation contract

For relative propagation time \(t\), ignition bias \(b\), current vehicle mass
\(m\), and a normalized direction \(\hat{u}\), the bridge samples the artifact
at

\[
  \tau = t-b.
\]

With axial force \(F_a(\tau)\), total flow \(\dot m_p(\tau)\), thrust scale
\(s_F\), and mass-flow scale \(s_m\), it returns

\[
  \mathbf F = s_F F_a\hat{u},\qquad
  \mathbf a = \frac{\mathbf F}{m},\qquad
  \dot m = -s_m\dot m_p.
\]

The force channel is already the net axial force. Hardware cant or divergence
represented by the artifact producer is not applied again. Every bound stream
has a negative additional-state derivative; several streams may drain the same
tank state and are summed. Complete stream binding is the default refusal rule.

## Discontinuities and event surfaces

The bridge exposes root functions for:

- burn start;
- every segment boundary;
- burn stop;
- protected vehicle dry mass; and
- every declared tank reserve.

Time roots use `g = t - (event_time + ignition_time_bias)`. Inventory roots use
`g = current_mass - minimum_mass`. Segment boundaries have explicit left and
right limits. Normal propagation should restart with the right-hand derivative
after locating a root. A reserve or dry-mass event stops this burn model, not
the complete orbit propagation.

Outside the firing interval, or after an inventory limit is reached, force and
all mass derivatives are exact zeros. The profile is never extrapolated.

## Estimation and variational equations

Three named, bounded drivers are provided:

| Driver | Meaning | Default |
|---|---|---:|
| `thrust_scale` | scales force only | 1 |
| `mass_flow_scale` | scales total and named tank drain only | 1 |
| `ignition_time_bias_s` | shifts profile and all time events together | 0 s |

Force and flow scale are intentionally independent at the orbit-determination
boundary. A calibrated propulsion realization should normally keep both at one;
independent estimation avoids silently forcing measured acceleration and tank
bookkeeping into a fabricated Isp identity.

The evaluation returns analytic local derivatives of acceleration with respect
to vehicle mass and all three drivers, plus total and tank mass-rate derivatives.
Inside a linear segment,

\[
  \frac{\partial\mathbf a}{\partial m}=-\frac{\mathbf a}{m},
  \quad
  \frac{\partial\mathbf a}{\partial s_F}=\frac{F_a\hat{u}}{m},
  \quad
  \frac{\partial\mathbf a}{\partial b}=-\frac{s_F\dot F_a\hat{u}}{m}.
\]

At a declared jump the timing derivative is mathematically non-smooth. The
bridge reports that fact instead of returning a misleading differentiable
value. Tests compare every smooth analytic parameter derivative against central
finite differences.

When a performance surface is attached, the same evaluation also returns
analytic acceleration, total-mass-rate and tank-rate partials with respect to
every condition. Lower and upper validity bounds are exposed as positive-inside
event surfaces, so a consuming propagator can stop the burn before any map
extrapolation.

For liquid systems, `DynamicFeedLine` can make the pressure condition a
propagated engine-manifold state instead of an algebraic input. Its coupled
evaluation drains the tank with line flow, drives compliant storage with the
difference between line and engine flow, and inserts the engine demand
gradient into the line-state Jacobian. See `DYNAMIC_FEED_LINE.md`.

## Integrity and refusal rules

Construction reopens and verifies the source artifact through its public hash
boundary. The bridge manifest binds the artifact hash, stream map, inventory
limits, driver definitions, units, semantics, and primary references into a
new deterministic `model_hash`.

The bridge refuses:

- unknown, duplicated, or incompletely bound streams;
- constraints without a corresponding propagated state;
- missing/non-finite additional states;
- zero or non-finite direction vectors;
- out-of-bounds or unknown estimator parameters; and
- unsupported schemas.

## Sidera integration shape

Sidera does not need a new astrodynamics architecture. Its numerical force-model
loop can adapt one bridge evaluation as follows:

1. convert the current epoch to burn-relative seconds;
2. resolve the guidance/attitude direction in the propagation frame;
3. pass current total mass and registered tank additional states;
4. add `acceleration_m_s2` to translational dynamics;
5. add `mass_derivative_kg_s` and named state derivatives;
6. register every `event_surfaces()` root once during force-model setup; and
7. map named analytic partials into Sidera's variational/estimation ordering.

When the dynamic line is selected, Sidera also integrates line flow, manifold
pressure, and the regulated-feed states. It registers the dynamic-line roots
and replaces the regulated feed's static injector-pressure root with the
minimum-manifold-pressure root.

No equivalent-constant reduction is used in this path.

## Primary references

- [Orekit 13.1.5 `PropulsionModel`](https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html)
- [Orekit 13.1.5 `EventDetectorsProvider`](https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html)
- [NASA transient liquid-engine modeling](https://ntrs.nasa.gov/citations/20040000363)
- [NASA SP-8112, Pressurization Systems for Liquid Rockets](https://ntrs.nasa.gov/citations/19760015212)
- [CCSDS 502.0-B-3](https://ccsds.org/Pubs/502x0b3e1.pdf)
