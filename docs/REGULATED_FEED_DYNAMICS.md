# Regulated pressurant and feed dynamics

`RegulatedFeedSystem` is a propagation-ready zero-dimensional feed-network
model. It supplies the pressure and thermal states needed by an operating-
condition performance surface while remaining independent of epochs, frames,
attitude and orbital integration.

The model contains:

- a finite-volume high-pressure pressurant bottle;
- a dynamic regulator/valve opening state;
- compressible choked or unchoked regulator flow;
- a growing ideal-gas tank ullage;
- tank and bottle wall heat transfer;
- liquid propellant inventory;
- quadratic propellant feed-line pressure loss; and
- protective cutoff event surfaces.

This is the first state producer for `rocket_propulsion_propagator_bridge_v1`.
It replaces a prescribed pressure history with six simultaneous ODE states.

## State vector

The fixed state ordering is:

1. propellant mass \(m_p\);
2. tank pressurant mass \(m_t\);
3. tank gas temperature \(T_t\);
4. bottle pressurant mass \(m_b\);
5. bottle gas temperature \(T_b\); and
6. regulator opening \(z\).

For rigid tank volume \(V_T\) and incompressible propellant density \(\rho_p\),

\[
 V_u=V_T-\frac{m_p}{\rho_p},\qquad
 p_t=\frac{m_tRT_t}{V_u},\qquad
 p_b=\frac{m_bRT_b}{V_b}.
\]

Construction refuses any state that does not leave positive ullage.
An empty liquid tank has full ullage volume and can be evaluated at zero
liquid demand. Negative liquid mass and nonzero liquid demand at depletion
raise a domain error.

## Regulator and pressurant flow

The preliminary controller uses clipped proportional pressure error,

\[
 z_c=\operatorname{clip}
 \left(\frac{p_{set}-p_t}{\Delta p_{open}},0,1\right),
 \qquad
 \dot z=\frac{z_c-z}{\tau_z}.
\]

Effective area is \(A=C_dzA_{max}\). Pressurant mass flow uses the ideal-gas
isentropic restriction equation and switches to the sonic expression when

\[
 \frac{p_t}{p_b}\le
 \left(\frac{2}{\gamma+1}\right)^{\gamma/(\gamma-1)}.
\]

Reverse flow is not modeled: regulator flow becomes zero when bottle pressure
is not greater than tank pressure.

The pressure head and actuator error are computed from the unrounded
ideal-gas node pressures. The displayed pressures or command can be equal
while the flow or valve rate is nonzero. The regulator combines area,
opening, pressure and gas scaling in squared form before exporting the
square root. For relative head no greater than `2**-54`, the unchoked
squared flow factor uses its leading term `2 (1 - p_t/p_b)`. Higher-order
terms in this factor are of order `(1 - p_t/p_b)**2`. This also preserves
representable flow when the relative head itself is below the binary64 range.

Tests compare forward flow and valve rate against Decimal calculations.
They cover equal displayed pressures with positive, zero and reverse head,
both compensating extremes of area scaling, and the square-root limit at a
relative head of `2**-1200`. These are numerical stress cases. The real
regulator's deadband, leakage, hysteresis and sensor resolution need measured
data from the selected configuration.

The time-dependent feed and the direct-acting regulator use the same
isentropic restriction law. Sonic classification compares the unrounded
pressure ratio to outward Decimal bounds of the critical ratio. The
bound is refined if the first comparison cannot separate the two values.
Exact rational equality belongs to the choked side. A comparison that
remains unresolved at 1280 decimal digits raises a domain error.

The earlier binary64 dynamic comparison disagreed with the 110-digit
reference at 16 of 120 synthetic pressure pairs bracketing the critical
ratio. The test now checks all 120 pairs in both models. It also tests
exact rational thresholds, refinement for a ratio within the first bound,
and refusal when the precision limit cannot resolve the side. A separate
test catches zero valve flow at a positive one-ULP pressure head in the
former direct-acting valve calculation.
[NASA Glenn isentropic relations](https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html)
[Python Decimal rounding](https://docs.python.org/3.12/library/decimal.html)

## Mass and energy equations

Positive engine demand \(\dot m_p\) drains the liquid:

\[
 \frac{dm_p}{dt}=-\dot m_p,\qquad
 \frac{dV_u}{dt}=\frac{\dot m_p}{\rho_p}.
\]

Pressurant mass is transferred exactly between bottle and ullage:

\[
 \frac{dm_t}{dt}=\dot m_r,\qquad
 \frac{dm_b}{dt}=-\dot m_r.
\]

For a calorically perfect pressurant, the well-mixed tank energy balance is

\[
 \frac{d(m_tc_vT_t)}{dt}=
 \dot m_rc_pT_b+UA_t(T_{w,t}-T_t)-p_t\frac{dV_u}{dt}.
\]

The fixed-volume bottle balance is

\[
 \frac{d(m_bc_vT_b)}{dt}=
 -\dot m_rc_pT_b+UA_b(T_{w,b}-T_b).
\]

The implementation solves these equations explicitly for both temperature
derivatives. Tests independently reconstruct both open-system energy balances
and the pressurant mass closure.

## Engine inlet pressure

The preliminary liquid line model is

\[
 \Delta p_f=K_f\dot m_p^2,
 \qquad
 p_{inj}=p_t-\Delta p_f.
\]

The model returns analytic derivatives of injector pressure with respect to all
six feed states and propellant-flow demand. For example,

\[
 \frac{\partial p_{inj}}{\partial \dot m_p}=-2K_f\dot m_p.
\]

These derivatives connect the feed states to the `supply_pressure_pa` axis of
the multidimensional performance surface.

## Algebraic feed–engine closure

Engine flow depends on inlet pressure, while inlet pressure depends on engine
flow. `evaluate_coupled_feed_propulsion()` solves

\[
 q=Q(p),\qquad p=P(\mathbf{x},q)
\]

with a relaxed fixed-point iteration and retains iteration count, flow residual
and pressure change. Failure to meet the requested tolerance raises an error;
the last iterate is never presented as a converged propulsion state.

The function also closes the implicit Jacobian. With
\(P_q=\partial P/\partial q\) and \(Q_p=\partial Q/\partial p\),

\[
 \frac{dp}{dx_i}=\frac{P_{x_i}}{1-P_qQ_p}.
\]

It returns closed feed-state derivatives of injector pressure, propellant flow
and spacecraft acceleration. The denominator is exposed, and a near-singular
algebraic loop fails closed. These are the derivatives a Sidera variational
equation or estimator needs; independently linearizing the engine and feed line
would miss the feedback term.

## Event surfaces

Positive-safe root functions are provided for:

- propellant reserve;
- minimum injector pressure;
- minimum bottle-to-tank pressure margin; and
- optional maximum tank pressure.

These functions provide positive values on the permitted side and negative
values on the violated side. The caller registers the events, checks the
initial margins and switches the burn off at the earliest protective root.
The native Sidera adapter propagates the seven orbital and vehicle-mass
states. It handles profile boundaries and the dry-mass cutoff. It does not
integrate the additional regulated-feed states or register their pressure
surfaces. The feed HTTP route evaluates operating points and sweeps without
propagating a burn.

Ullage capacity is tested using the full expression
`V_tank - m_p/rho` before binary64 rounding. Initial gas mass, tank pressure
and pressure partials use that volume. Each pressure cutoff forms the complete
signed margin from the state, then rounds the result once. An out-of-range
margin raises a domain error. Reserve and overpressure queries do not evaluate
the regulator dynamics or the static injector loss.

Decimal oracles exercise both sides of each pressure limit and exact equality,
including synthetic margins below the spacing of the displayed pressures.
The tank-capacity test checks a positive ullage of approximately
`1.93e-19 m³` and rejects the next larger representable propellant mass.
For a sealed adiabatic tank under constant demand `q`, an independent
root-time check uses `p V^gamma = constant`, with
`V(t) = V(0) + q t / rho` and
`p_t(cutoff) = p_inj_min + R_line q²`. It compares the resulting time and
liquid mass with a numerical integration of the six-state model. These tests
verify the equations and event functions. Integration into a mission
propagator and response of a physical shutdown chain require separate evidence.

A configured zero reserve has an exact root at an empty tank. The state
can be evaluated there in coast mode. Numerical methods may query past the
root while constructing a trial step. Physical feed queries still require
nonnegative inventory and zero demand at the empty boundary.

## Coupled trajectory reference

`rocket_propulsion.integrations.regulated_feed_trajectory.propagate_regulated_feed`
integrates the motor bridge and six feed states with Cartesian position and
velocity. It accepts a single liquid stream, a total-flow performance surface
with injector pressure in Pa and a fixed inertial thrust direction. SciPy is
loaded on demand. Supply `max_step_s` for the regulator and event timescales.
The vehicle mass equals the constant nonliquid mass plus the integrated
liquid inventory. For this integration the convergence test scales to actual
flow with zero scale floor, so a small positive motor demand cannot converge
to zero feed flow under an absolute tolerance inherited from kg/s-scale burns.

The integrator checks every configured feed margin, the pressure-map bounds
and optional vehicle dry mass. It checks the right-hand demand before starting
a new profile segment, even when the requested final epoch lies exactly on a
command boundary. A continuous crossing uses the left limit at the previous
segment end. The root finder uses endpoint pressure values only for motor-map
trial queries and evaluates the signed protection margins at physical pressure.
It returns the first stop state, final state and signed final feed margins.
At commanded shutdown the motor force and liquid demand become zero. Gravity,
regulator motion and gas heat transfer continue, so a coast margin may become
negative even after the motor has stopped.

Analytic tests use sealed-tank adiabatic expansion, an integrated
pressure-dependent motor law, a wall-heating transient, the rocket equation
with force-free coast and a circular orbit under central gravity. They also
check empty and nonzero reserves, pressure jumps, initial unsafe states, map
endpoint bracketing, dry mass and small flow. The mutation suite reintroduces
faults in the integration, cutoffs and coast to verify test sensitivity.

For each accepted DOP853 step the reference evaluates the dense solution
at both endpoints and three interior points. A signed margin changing
from positive to nonpositive brackets a crossing. A rate changing from
negative to nonnegative brackets a possible minimum, where the physical
margin is checked even if the sampled values were positive. The
injector rate includes the explicit commanded-flow slope and the implicit
pressure/flow response. The first detected root takes precedence when
several margins cross inside a scan interval. The solver resumes in coast
from the dense state at that root.

The independent sealed-tank ramp regression gives a pressure violation
from 0.360370 to 0.419999 s. The original step-end check missed it at a
1 s maximum step. The dense check detects it and the shallower 0.01 Pa
excursion, rejects a 0.01 Pa near miss, and orders a pressure trip ahead
of a later dry-mass trip. A separate closed-vessel thermal case exercises
the bottle-to-tank head margin. The time error depends on the ODE
interpolation and the slope at the crossing, especially near a grazing
limit. Multiple extrema between scan points and error near a tangency
can still hide a crossing or change its classification. The maximum
step and tolerances must be justified against the relevant transient
and checked under refinement. The adapter
models instantaneous commanded thrust removal. Shutdown latency, valve
response, thrust tailoff, attitude and post-cutoff feed safety require
mission-specific models and measured evidence.
[SciPy DOP853 dense output](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.DOP853.html)

## Fidelity limits

The model is a preliminary single-node network. By itself it intentionally
omits:

- feed-line inertance and distributed compliance;
- priming, water hammer and acoustic modes;
- ullage stratification;
- regulator hysteresis, stiction and detailed poppet mechanics;
- real-gas effects;
- cryogenic phase change and propellant heating;
- dissolved pressurant and diaphragm/bladder mechanics; and
- multi-tank bipropellant mixture-ratio control.

The optional `DynamicFeedLine` adds one calibrated inertance/compliance mode
and engine-manifold pressure state. It does not remove the remaining
distributed water-hammer, priming, cavitation, two-phase, or structural
limitations. When that model is composed, static line resistance here must be
zero so the pressure loss is not counted twice.

Fast valve events or structural pressure qualification require a validated
nodal/transient tool. NASA’s recent GFSSP work explicitly solves mass, energy,
state and branch momentum equations and provides the higher-fidelity reference
architecture; this implementation is the compact orbit-propagation layer, not
a substitute for that system model.

## Primary references

- [NASA SP-8112, Pressurization Systems for Liquid Rockets](https://ntrs.nasa.gov/citations/19760015212)
- [NASA SP-8080, Liquid Rocket Pressure Regulators and Valves](https://ntrs.nasa.gov/citations/19740002611)
- [NASA, Nodal Modeling of Liquid Propellant Feed and Pressurization System](https://ntrs.nasa.gov/citations/20240003493)
- [NASA Glenn, Mass Flow Choking](https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html)
- [NASA NESC, Treatment of Transient Pressure Events](https://ntrs.nasa.gov/citations/20220006583)
