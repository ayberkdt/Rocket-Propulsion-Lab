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

Crossing zero stops the burn model but does not terminate orbit propagation.

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
