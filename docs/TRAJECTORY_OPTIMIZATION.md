# Propulsion for trajectory optimization and estimation

This page covers what a trajectory optimizer, a targeting loop or an orbit
determination filter needs beyond replaying a burn:

1. a **control-oriented engine model** whose throttle is a decision variable;
2. **smooth** throttle characteristics, so gradients have no kinks;
3. **impulsive equivalents** of finite burns for initial guesses;
4. a verified recipe for **variational equations** across burn
   discontinuities.

## 1. Throttle-commanded engine

`ThrottleableEngineModel` (`propulsion/burns/throttleable.py`):

\[
  F(u)=F_\mathrm{rated}\,f_T(u)\,s_F-p_aA_e,\qquad
  \dot m(u)=\dot m_\mathrm{rated}\,f_m(u)\,s_{\dot m},
\]

with \(u=0\) (off) or \(u\in[u_\min,u_\max]\). Values in \((0,u_\min)\) are
refused, because a deep-throttle gap is a hard engine limit. A direct method
models it with a discrete on/off variable or a complementarity constraint
rather than by interpolation. Total flow is split into tank states by fixed
`stream_fractions` that sum to one, so total and per-tank drains close
exactly.

Each evaluation returns analytic partials with respect to the throttle (the
control), vehicle mass, the supplied direction (3×3), ambient pressure and the
`thrust_scale` / `mass_flow_scale` drivers. Each is checked against central
finite differences in `tests/test_burn_throttleable.py`. Dry-mass and reserve
roots, mass closure, and a hash-verified document follow the same conventions
as the other models.

## 2. Smooth throttle characteristics

Measured engines often lose specific impulse at deep throttle, so
\(f_T(u)\neq f_m(u)\). Both fractions can be given as `MonotoneCubicCurve`, a
shape-preserving piecewise cubic Hermite interpolant (PCHIP; Fritsch–Carlson
with Fritsch–Butland interior slopes and the standard end conditions):

- it passes through every node and is continuously differentiable, so an
  optimizer's Jacobian has no kinks at the nodes (unlike linear
  interpolation);
- it is monotone wherever the data are monotone, so it never overshoots,
  which matters for physical fractions near limits;
- it never extrapolates.

The tests check node reproduction, exact reproduction of linear data,
monotonicity on random monotone data (2 000 samples per curve), slope
continuity at every knot, and the derivative against finite differences.

The N-D `RectilinearPerformanceSurface` is multilinear (C⁰, with its kinks
exposed as restart roots) by default. With `interpolation="cubic_hermite"` it
becomes C¹ with the same PCHIP nodal derivatives (see
[`OPERATING_CONDITION_SURFACES.md`](OPERATING_CONDITION_SURFACES.md)), which
is the choice for gradient-based optimization over operating conditions.

## 3. Impulsive equivalent of a finite burn

`free_space_burn_metrics(artifact, initial_mass_kg)` integrates a tabulated
burn in free space:

\[
  \Delta v=\int\frac{F(t)}{m(t)}\,dt,\qquad
  t^\ast=\frac{\int t\,a\,dt}{\int a\,dt},\qquad
  c_\mathrm{eff}=\frac{\Delta v}{\ln(m_0/m_f)} .
\]

Thrust and flow are piecewise linear, so \(m(t)\) is piecewise quadratic and
exact mass partitions precede a normalized \(2n\)-point Gauss–Legendre rule,
with exact polynomial identities where applicable. Independent directed
Bernstein enclosures bound Δv and moment; output rounding and mass-logarithm
enclosures are included in the reported numerical bounds. `quadrature_points`
is the configured refined order, not a total evaluation count. See
`FREE_SPACE_TABULATED_NUMERICS_REVIEW.md` for limits and range refinement.
\(t^\ast\) places the equivalent impulse. For constant positive thrust and
positive drain it lies after the time midpoint; varied profiles require the
actual acceleration-weighted moment.

The tests check the quadrature rule (exact for polynomials up to degree
\(2n-1\)), the closed forms for constant thrust (\(\Delta v\) and \(t^\ast\) to
1e-12), and agreement with an event-aware RK4 integration of the bridge on
random multi-segment profiles.

## 4. Variational equations across burn discontinuities

For parameters \(p\) (the bridge drivers) a filter integrates
\(S=\partial y/\partial p\) with

\[
  \dot S=\frac{\partial f}{\partial y}S+\frac{\partial f}{\partial p},
\]

using `acceleration_wrt_mass`, `acceleration_wrt_parameters`,
`mass_rate_wrt_parameters` and `additional_state_rates_wrt_parameters`. Profile
time events (start, segment boundaries, stop) move one-for-one with
`ignition_time_bias_s`. At each one the sensitivity jumps
(saltation):

\[
  S^+=S^-+\left(f^- - f^+\right)\frac{\partial\tau}{\partial p},\qquad
  \frac{\partial\tau}{\partial b}=1 .
\]

`bridge.time_event_jump(...)` returns the propulsion part of \(f^- - f^+\)
(acceleration, mass rate and per-tank rates) at an event, with inventory
limits off. `tests/test_burn_variational.py` integrates exactly this recipe
through a profile with start, mid-burn and stop jumps. All three drivers match
central finite differences of whole propagations to six significant digits.
A guard test shows that the ignition-bias sensitivity is wrong without the
jump. For electric propulsion the same test chains
`acceleration_wrt_available_power` through a time-varying power model.

## References

- F. N. Fritsch and R. E. Carlson, *Monotone Piecewise Cubic Interpolation*,
  SIAM J. Numer. Anal. 17(2), 1980; F. N. Fritsch and J. Butland, SIAM J. Sci.
  Stat. Comput. 5(2), 1984 (cited bibliographically).
- [NIST Digital Library of Mathematical Functions §3.5(v), Gauss quadrature](https://dlmf.nist.gov/3.5#v)
- [Navigation Filter Best Practices, NTRS 20180003657](https://ntrs.nasa.gov/citations/20180003657)
- [Orekit 13.1.5 `PropulsionModel`](https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html)
