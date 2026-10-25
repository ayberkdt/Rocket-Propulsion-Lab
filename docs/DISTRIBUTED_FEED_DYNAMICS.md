# Distributed feed dynamics: water hammer, regulator hysteresis, ullage stratification

Three self-contained modules extend the single-node regulated feed
(`feed_system.py`) and the two-state line surrogate (`feed_line.py`). They are
analysis tools: none is wired into the propagator bridge, and none changes the
behaviour of an existing module. Every equation constant below was read from
the cited page on 2026-10-25; the locators are in
[the reference index](REFERENCES.md), section "Feed dynamics sources".

## 1. Water hammer by the method of characteristics

Module: `rocket_propulsion.propulsion.burns.water_hammer`.

One uniform pipe runs from a constant-pressure reservoir to a valve. The
governing equations are the NASA/TM-20210022275 Appendix B pair (eqs. B-1, B-2),
with volume flow `Q`, area `A`, Darcy factor `f` and body acceleration
`alpha sin(theta)`:

```text
dp/dt + (rho a^2 / A) dQ/dx = 0
(1/A) dQ/dt + (1/rho) dp/dx + alpha sin(theta) + f Q|Q| / (2 D A^2) = 0
```

Along `dx/dt = +a` and `-a` they reduce to compatibility equations. On the
Courant-number-one grid `dt = dx / a`, with `B = rho a / A`,
`R = rho f dx / (2 A^2 D)` and `S alpha = rho dx alpha sin(theta)` (eqs. B-11 to
B-17):

```text
Cp = p[i-1] + (B - R|Q[i-1]|) Q[i-1] - S alpha
Cm = p[i+1] - (B - R|Q[i+1]|) Q[i+1] + S alpha
p  = (Cp + Cm) / 2,   Q = (Cp - Cm) / (2 B)
```

At the reservoir, pressure is fixed. The valve obeys
`p_valve - p_down = (K_open / tau(t)^2) Q|Q|`. This is the quadratic branch law
of the NASA network code (Majumdar and Steadman, eq. 2). `tau` is the relative
flow coefficient, prescribed by a callable such as `LinearValveLaw`. A closed
valve gives `Q = 0`. The initial state is the steady solution of the same
discrete equations, so a valve that does not move leaves the line steady to
round-off.

`pipe_wave_speed_m_s` implements eq. 2,
`a = sqrt(beta/rho) / sqrt(1 + beta D C / (E t))`. The boundary factor `C` comes
from Figure 5, for anchoring cases I to IV. The thick-walled forms apply for
`D/t <= 25`. `joukowsky_surge_pa` implements eq. 1, `dP = -rho a dV`.

Evidence (`tests/test_water_hammer.py`):

- Case Study 8-1 of the TM is reproduced in lbm-in-s units: 0.088 lb/s,
  0.62 in², 39,370 in/s, 0.036 lb/in³. Both the formula and the MOC run give
  the printed 14.5 psi Joukowsky surge, which the TM itself recovered by MOC.
- The surge converges to Joukowsky as the closure time shrinks. For any
  closure no longer than the wave round trip `2L/a`, the frictionless surge
  equals `rho a V0` to 1e-10. Slower closures give monotonically smaller
  surges, as stated on TM p. 24.
- The first step after an instantaneous closure rises by exactly `rho a V0`,
  even with friction. Line packing then lifts the peak above Joukowsky, and
  that peak converges with grid refinement.
- Mass balance holds:
  `(A/(rho a^2)) integral(p dx) - integral(Q_in - Q_out) dt -> 0`.
- Acoustic energy balances: the change in kinetic plus elastic energy equals
  boundary power minus Darcy dissipation minus body-force work. Both balance
  residuals fall at first order with the reach count. In a closed, lossless
  line, energy is conserved to round-off.
- A pressure at or below the vapour pressure raises `DomainError`, because
  column separation is outside the model.

Validity limits: a single uniform pipe, with constant density and wave speed,
steady Darcy friction (no unsteady friction), and no fluid-structure
interaction. There is no cavitation or column separation, and no gas
(two-phase) content. The 1/tau² valve law assumes the valve flow coefficient
scales with `tau`. The Joukowsky bound can be non-conservative with two-phase
flow, area changes, tees, or vapour-pressure excursions (TM p. 20).

Priming is not modelled. The NASA sources read here (TM p. 13 and the 2019
TFAWS priming presentation) describe priming of evacuated or gas-filled
lines, and bound the dead-end surge by `rho a dV` for a known impact velocity.
They also report GFSSP comparisons with priming experiments. None of them gives
a closed-form model of the advancing front, so none is implemented. A priming
analysis needs a validated network code.

## 2. Regulator characteristic with lockup, droop and hysteresis

Module: `rocket_propulsion.propulsion.burns.regulator_hysteresis`.

`DirectActingRegulator` implements the NASA SP-8080 steady force balance of a
direct-acting, pressure-unbalanced regulator (eq. 2, printed p. 72):

```text
P_A A_A + L - Y X - P_R A_A + P_R A_i - B - Z X - P_S A_i = 0
```

Its differential form is eq. 3,
`dP_R = -dP_S A_i/(A_A - A_i) - dX (Y + Z)/(A_A - A_i)`. Larger flow needs a
larger stroke `X`, so the outlet pressure droops with flow; this is the
position error of eq. 3. Lockup, where all flow stops (Glossary, p. 104), is
the pressure at which the closing valve seats, `X = 0`. Coulomb damping is a
constant-magnitude friction `F_c` opposing the motion (Glossary). Under it,
the stroke moves only when the force imbalance exceeds `F_c`. This gives an
opening branch, a closing branch shifted up by `2 F_c / (A_A - A_i)`, and a
stuck region in between. The stroke is the memory state.
`regulator_operating_point` and `sweep_regulator` carry that state.

Inlet-valve gas flow is
`Cd * min(X dA/dX, A_max) * (NASA Glenn mass-flow eq. 10)`. The Mach number is
taken from the isentropic pressure ratio (NASA Glenn isentropic eq. 6) and is
choked at `M = 1`. The default `Cd = 0.6` is the SP-8080 preliminary-sizing
value (p. 76). `feed_configuration_with_regulated_pressure` copies an
operating point's outlet pressure into a `RegulatedFeedConfiguration` set
pressure. The copy repeats all validation, so an overpressure is refused.

Evidence (`tests/test_regulator_hysteresis.py`):

- Branch pressures satisfy eq. 2 with friction `+F_c`, `-F_c` or 0.
- Finite differences reproduce both sensitivities of eq. 3.
- Choked and unchoked flows match the NASA Glenn forms to 1e-12.
- Lockup equals the closing-branch zero-stroke pressure.
- Droop is monotone and, in choked flow, linear in demand.
- With zero friction the characteristic does not depend on the path. With
  friction, the branches are separated by exactly the band at equal stroke.
- A small unchoked flow reversal stays stuck inside the band.
- Saturation at full stroke is flagged. A demand beyond capacity raises
  `DomainError`.

Validity limits: the model is quasi-static and has no poppet dynamics. The
SP-8080 "response" time of 63 % (p. 105), damping orifices, and dome-loaded,
piloted or nonmodulating designs are not modelled, and neither are
temperature and acceleration errors, seat leakage or reverse flow. Spring
rates are constant; Belleville or variable-rate springs are not represented.
SP-8080 recommends a supply- and flow-error budget of about 25 % of tolerance,
and the same limit for hysteresis and friction (pp. 11, 76). These are design
targets, not model outputs.

## 3. Stratified (n-layer) ullage

Module: `rocket_propulsion.propulsion.burns.ullage_stratification`.

The ullage of a vertical cylinder is a stack of Lagrangian gas layers at one
uniform pressure. Pressurant enters at the top as a new layer at the inlet
temperature. Older layers move down as the liquid drains. Each layer exchanges
heat with the wall strips it overlaps, and the bottom layer also exchanges
heat with the liquid. This is the one-dimensional, no-axial-conduction picture
of NASA TN D-2797: assumptions (1)-(13), eqs. 1-3. The layered discretization
follows NASA MSFC's multi-node ullage model (LeClair and Majumdar, 2025).
Wetted wall is held at the liquid temperature (TN D-2797 assumption 13). End
domes are not modelled.

Heat-transfer coefficients use the NASA MSFC form of Ring's correlation
(Majumdar and Steadman, eqs. 10-13 and p. 5; LeClair and Majumdar, eqs. 1-2):

```text
h = K_H C (k / L) (Gr Pr)^n,   Gr = g beta dT L^3 rho^2 / mu^2,   Pr = mu cp / k
ullage-wall:   C = 0.54, n = 0.25
ullage-liquid: C = 0.27, n = 0.25,   K_H = 1,   L = tank diameter
```

The `L^3` comes from the nomenclature of NASA 20190030427; it is not legible
in the scanned 1999 draft. The ideal-gas `rho = p/(RT)` and `beta = 1/T` are
evaluated at the layer temperature, which is an implementation choice.
`K_H` is a calibration factor; NASA set it to 0.5 in one model to match a
measured boil-off.

Each step:

1. Constant-pressure heat exchange. Sink temperatures are frozen over the
   step, and each layer relaxes exponentially toward them, so the step is
   stable for any `dt`. The wall receives exactly the heat the gas loses.
2. Isentropic compression of the old layers.
3. A new layer is admitted at the inlet temperature.
4. The new pressure is solved from `sum(m R T) = p V`.

In regulated mode the admitted mass restores the set pressure. If heating
alone exceeds it, nothing is admitted and the pressure rises, as at regulator
lockup. Merging every layer after each step gives the well-mixed single node,
which is the ullage assumption of `RegulatedFeedSystem`.

`compare_ullage_models` returns, for the same inputs:

- the stratified and well-mixed histories (pressure, mean, top and bottom
  temperatures, added mass, heat to wall and liquid);
- the ideal adiabatic mass `p dV / (R T_in)` (TN D-2797 eq. 24);
- the collapse factors (eq. 25);
- the SP-8112 eq. 14 first-cut ullage mass.

Evidence (`tests/test_ullage_stratification.py`):

- Adiabatic regulated outflow needs exactly the ideal mass (collapse factor 1).
- Adiabatic expansion keeps `p V^gamma` constant to 1e-12.
- Closed adiabatic filling converges to `dp/dt = gamma R T_in mdot / V`.
  Its pressure error and energy residual halve with the step.
- Mass is exact. The energy residual vanishes to round-off at constant
  pressure.
- The added mass converges with the step size.
- The final profile is hot at the top and bounded by the inlet and liquid
  temperatures. The wetted wall stays at liquid temperature.
- With heat transfer, both collapse factors exceed 1, and the final mass
  stays below the SP-8112 eq. 14 bound.
- In an interface-dominated case, the stratified ullage loses less heat to
  the liquid and needs less pressurant than the well-mixed one, as SP-8112
  p. 90 states.

In the wall-dominated default case of the tests (a development run on
2026-10-25 gave collapse factors 1.30 stratified and 1.28 well-mixed),
stratification needed slightly more pressurant than full mixing: hot top gas meets cold freshly exposed wall, and
the correlation grows as `dT^1.25`. Neither model should be assumed
conservative without hardware data.

Validity limits: a calorically perfect ideal gas, and constant transport
properties supplied by the user. There is no condensation, evaporation,
dissolution or propellant-vapour partial pressure, no mass transfer at the
interface, no inlet diffuser jet or mixed convection, no slosh, no dome heat
transfer and no axial wall conduction. The outflow is constant.
LeClair and Majumdar report that the natural-convection correlation
under-predicts measured helium consumption, and that a mixed-convection
correlation did better. Use the results for trade studies and calibrate
against test data.

## Unverified or omitted

- Wylie and Streeter, *Fluid Transients*, was not read. The method of
  characteristics is implemented and cited from NASA/TM-20210022275
  Appendix B.
- E. Ring (1964) was not read. Its constants are used only as quoted, with
  page locators, by NASA MSFC papers.
- No priming model is included (see section 1).
- The property evaluation temperature of the convection correlation is not
  specified in the sources.

## References

- NASA/TM-20210022275, transient pressure events: <https://ntrs.nasa.gov/citations/20210022275>
- NASA SP-8080, regulators and valves: <https://ntrs.nasa.gov/citations/19740002611>
- NASA SP-8112, pressurization systems: <https://ntrs.nasa.gov/citations/19760015212>
- NASA TN D-2797, stratified outflow pressurization: <https://ntrs.nasa.gov/citations/19650012768>
- Majumdar and Steadman, propellant tank pressurization: <https://ntrs.nasa.gov/citations/19980236009>
- LeClair and Majumdar, multi-node ullage: <https://ntrs.nasa.gov/citations/20240006787>
- Bandyopadhyay et al., multi-node cryogenic tank model: <https://ntrs.nasa.gov/citations/20190030427>
- LeClair, Boehm and Majumdar, priming transients: <https://ntrs.nasa.gov/citations/20190030443>
- NASA Glenn, Mass Flow Choking: <https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html>
- NASA Glenn, Isentropic Flow: <https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html>
- NASA Glenn, Specific Heats: <https://www.grc.nasa.gov/www/k-12/airplane/specheat.html>
