# Power-limited electric propulsion

`PowerLimitedPropulsionModel` (in `propulsion/burns/electric.py`) gives an
orbit propagator the force, propellant drain, consumed power, analytic partials
and switching roots of an electric thruster whose operating point is set by the
power available on the bus. It follows the same boundary as
[`PROPAGATOR_PROPULSION_BRIDGE.md`](PROPAGATOR_PROPULSION_BRIDGE.md): no epoch,
frame, attitude or orbit is computed here.

## Ownership

| Rocket Propulsion Lab | Astrodynamics framework (Sidera) |
|---|---|
| throttle table, throttle path, hysteresis | solar distance, array degradation, eclipse |
| thrust, flow, consumed power per mode | bus loads and the resulting available power |
| mode-switch, saturation, reserve and dry-mass roots | event location and mode state |
| partials with respect to mass, direction, power, drivers | chain rule into its state/parameter ordering |

The power model stays in the propagator because it depends on the trajectory.
The propulsion model only receives the scalar available power and returns
`dT/dP` and `dmdot/dP` so the chain rule can be closed there.

## Throttle laws

### Discrete throttle path

An `ElectricThrottleTable` is evidence: each `ElectricThrottleLevel` carries
the bus input power \(P_k\), thrust \(T_k\) and mass flow \(\dot m_k\), and the
table is bound to its source by `source_id`, `source_sha256` and a content
hash. Each level is checked for physical consistency:

\[
  I_{sp,k}=\frac{T_k}{g_0\dot m_k},\qquad
  \eta_k=\frac{T_k^2}{2\dot m_k P_k}\le 1 .
\]

A `DiscreteThrottlePath` is an ordered subset with strictly increasing
\(P_k\). `DiscreteThrottlePath.from_table(table, objective=...)` builds the
maximum-thrust or maximum-Isp path: levels are visited by increasing power and a
level is kept only if it strictly improves the objective, so more power never
selects a worse level. Qualification tables can list several operating points
at similar power with different thrust and Isp (the NEXT table summarized in the
reference below has 40 conditions); a path keeps one level per power step.

Selection uses the effective power \(P_e=P_\mathrm{avail}-P_\mathrm{margin}\).
From level \(k\) the model drops when \(P_e<P_k\) and rises to \(k+1\) only when
\(P_e\ge P_{k+1}+h\), where \(h\) is the up-switch hysteresis. Inside a level the
force does not depend on power.

### Fixed-efficiency law

For low-thrust trajectory optimization a continuous law is often preferred.
With constant total efficiency \(\eta\) and specific impulse:

\[
  T=\frac{2\eta P}{g_0 I_{sp}},\qquad \dot m=\frac{T}{g_0 I_{sp}},
  \qquad P_\min\le P\le P_\max .
\]

This is the jet-power identity \(\tfrac12\dot m v_e^2=\eta P\). Above
\(P_\max\) the thruster saturates and \(\partial T/\partial P=0\); below
\(P_\min\) it is off. Restart needs \(P_\min+h\).

## Mode state and events

The level index (`-1` = off) is a propagator-owned discrete state, as in an
event-driven maneuver model:

1. initialize with `model.select_level(P_avail)`;
2. register `model.switching_surfaces(level)`; each `PowerSwitchSurface` is a
   root \(g=P_\mathrm{avail}-P_\mathrm{threshold}\) with a crossing direction and
   a fixed `to_level`;
3. after locating a root, set `level = surface.to_level` and re-register.

The target level is fixed by the crossing direction instead of being recomputed
from the power at the located root. A root found a round-off on either side of
its threshold therefore cannot leave the mode unchanged and re-trigger the same
event. A saturation surface (`crossing="either"`, `to_level == from_level`) only
marks the kink in \(\partial T/\partial P\) and asks for an integrator restart.

`evaluate` never changes the level. A Runge-Kutta stage that probes just past a
down-switch root stays on the current level's smooth branch and reports
`power_deficit_w > 0`. A deficit after an accepted step means the consumer
missed a switching event. As for the transient bridge, a consumer that locates
the dry-mass and propellant-reserve roots (`inventory_surfaces()`) itself
should pass `apply_inventory_limits=False`.

## Partials

`ElectricPropulsionPartials` contains, for the resolved direction
\(\hat u = \mathbf d/\lVert\mathbf d\rVert\):

\[
  \frac{\partial\mathbf a}{\partial m}=-\frac{\mathbf a}{m},\quad
  \frac{\partial\mathbf a}{\partial\mathbf d}
  =\frac{T}{m\lVert\mathbf d\rVert}(I-\hat u\hat u^{\mathsf T}),\quad
  \frac{\partial\mathbf a}{\partial P}=\frac{1}{m}\frac{\partial T}{\partial P}\hat u,
\]

the mass-rate partial with respect to power, and partials with respect to the
bounded `thrust_scale` and `mass_flow_scale` drivers. All are checked against
central finite differences in `tests/test_burn_electric.py`.

## Verification

`tests/test_burn_electric.py` integrates both laws with an event-driven RK4
loop under a linearly decaying power supply. The discrete path switches
L3 → L2 → L1 → off at the analytic times; final mass and velocity match the
piecewise rocket equation to 1e-10. The fixed-efficiency law saturates, follows
power down to shutdown, consumes exactly \(2\eta E/v_e^2\) of propellant for the
delivered electrical energy \(E\), and obeys \(\Delta v=v_e\ln(m_0/m_f)\).

## Limits

- The throttle table is user-supplied evidence; no thruster data are bundled.
- `duty_cycle` represents an averaged on/off pattern, not individual pulses.
- Cathode/neutralizer flows must be included in the level's total mass flow.
- Thruster degradation with throughput, multiple thruster strings and
  power-processing-unit efficiency curves are not modelled separately.
- Available power, solar-array physics and eclipse remain the consumer's model.

## Primary references

- [NASA's Evolutionary Xenon Thruster (NEXT) Ion Propulsion System Information Summary, NTRS 20090004685](https://ntrs.nasa.gov/citations/20090004685)
- [Orekit 13.1.5 `PropulsionModel`](https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html)
- [Orekit 13.1.5 `EventDetector`](https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html)
