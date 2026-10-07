# Operating-condition performance surfaces

`RectilinearPerformanceSurface` supplies the state-dependent layer between a
time-resolved burn artifact and the numerical propagator. It maps an explicit
set of operating conditions to independent axial-force and total-tank-flow
multipliers.

Typical axes include:

- supply or injector inlet pressure;
- ambient pressure or altitude-derived back pressure;
- commanded/realized throttle;
- mixture ratio;
- thermal state;
- available electrical power; and
- thrust-vectoring deflection.

NASA's parametric aerospike model used mixture ratio, power level,
thrust-vectoring level and altitude to generate trajectory-facing engine
performance tables. The implementation here generalizes that table boundary
without claiming that any particular set of axes is universally sufficient.

## Data contract

Each `PerformanceAxis` declares a stable name, SI-oriented unit and at least two
strictly increasing coordinates. `PerformanceScalePoint` records

\[
  s_F(\mathbf{x}),\qquad s_{\dot m}(\mathbf{x})
\]

at one Cartesian grid node. The complete tensor-product grid must be present
exactly once. Missing, duplicated or off-grid points are rejected during
construction.

Force and mass-flow scales are deliberately independent. Their ratio represents
a change in effective system specific impulse; forcing the two scales to be
equal would erase measured off-design performance.

The source identifier, source SHA-256, axes, units, nodes and references are
bound into `surface_hash`. JSON reopening reconstructs the surface and rejects
hash-altered content.

## Interpolation and gradients

Within one N-dimensional cell, tensor-product multilinear interpolation is
used:

\[
  y(\mathbf{x}) =
  \sum_{\mathbf{c}\in\{0,1\}^N}
  y_{\mathbf{c}}
  \prod_{j=1}^{N}
  \left[c_j\lambda_j+(1-c_j)(1-\lambda_j)\right].
\]

The code differentiates the same expression analytically with respect to every
axis. It therefore supplies exact local partials for variational equations,
state-transition matrices and orbit determination; the propagator does not
need to finite-difference a table lookup.

At grid-cell boundaries the interpolated value is continuous, but the gradient
need not be. Consumers should treat a knot crossing as a derivative-regime
change when high-order sensitivity continuity is required.

## Validation domain and events

No extrapolation is permitted. A coordinate more than numerical endpoint
tolerance outside its axis raises a domain error. The surface exposes two
positive-inside root functions per axis:

\[
  g_\mathrm{lower}=x-x_\min,
  \qquad
  g_\mathrm{upper}=x_\max-x.
\]

These roots allow Sidera or another numerical propagator to stop the propulsion
model exactly at the validated-domain boundary. This is materially safer than
discovering an invalid map query after an integration stage has stepped beyond
the applicable data.

## Propagator bridge composition

When a surface is attached to `PropagatorPropulsionBridge`, a base artifact
sample \(F_0,\dot m_0\) becomes

\[
  F=s_F^\mathrm{est}\,s_F(\mathbf{x})F_0,
  \qquad
  \dot m=-s_{\dot m}^\mathrm{est}\,
  s_{\dot m}(\mathbf{x})\dot m_0.
\]

Every named stream receives the same map flow multiplier, preserving its base
mixture split and exact total-flow closure. The evaluation returns:

- the conditions actually consumed;
- interpolated force and flow scales;
- acceleration partials with respect to every condition;
- total mass-rate partials; and
- per-tank mass-rate partials.

If, for example, `supply_pressure_pa` is also a propagated tank/feed state,
Sidera maps the named condition partial directly into the corresponding column
of its state Jacobian. Altitude-derived ambient pressure can instead be composed
through the atmosphere model with the chain rule.

Conditions are required only while the burn is active. This prevents an
inactive maneuver from querying or extrapolating a propulsion map unnecessarily.

## Limits

This layer interpolates supplied evidence; it does not generate that evidence.
It does not replace:

- feed-system differential equations;
- regulator hysteresis or valve dynamics;
- thermal soak and two-phase effects;
- model-discrepancy calibration;
- covariance across map nodes; or
- validation outside the declared domain.

The next fidelity step is to couple pressure/thermal/feed states to these axes
and propagate their uncertainty with the existing hierarchical ensemble.

## Primary references

- [NASA, Parametric Model of an Aerospike Rocket Engine](https://ntrs.nasa.gov/citations/20000031654)
- [NASA-STD-7009B, Standard for Models and Simulations](https://standards.nasa.gov/standard/NASA/NASA-STD-7009)
- [NASA transient liquid-engine modeling](https://ntrs.nasa.gov/citations/20040000363)
- [Orekit `PropulsionModel`](https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html)
