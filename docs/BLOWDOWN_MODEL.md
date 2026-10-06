# Pressure-fed blowdown model

## Status and ownership

This is an isolated, preliminary L2 propulsion model. It produces thrust and
tank-drain histories for Rocket Propulsion Lab and a future Sidera adapter; it
does not propagate an orbit, select a frame, or determine thrust direction.

The implementation is in
`rocket_propulsion.propulsion.burns.blowdown`. Every public model/function has
a `References` section, and every generated tabulated artifact carries stable
reference identifiers inside its deterministic content hash.

## Tank pressure model

`BlowdownTankModel` represents a rigid tank, an incompressible liquid and a
declared polytropic ullage gas:

```text
V_u(m) = V_u0 + (m_0 - m) / rho_p
p(m)   = p_0 [V_u0 / V_u(m)]^n
```

Here `m` is remaining propellant, `rho_p` is the declared liquid density and
`n` is the user-supplied polytropic exponent. `n = 1` is the isothermal
ideal-gas limit. The exponent is not inferred from a pressurant name.

The inverse relation is used to stop a burn exactly at the declared minimum
supply pressure. A propellant reserve is an independent cutoff; whichever
leaves more propellant in the tank binds first.

Primary basis: NASA SP-8112 describes blowdown pressurization as ullage-gas
expansion with decreasing tank pressure and the resulting pressure-fed engine
performance decline. NASA CR-131400 documents time-history performance
modeling for an unregulated monopropellant blowdown system.

## Pressure-to-performance models

Two interchangeable, non-extrapolating models implement
`PressurePerformanceModel`:

- `PressurePerformanceLaw` uses declared power correlations
  `F/F_ref = (p/p_ref)^a`, `Isp/Isp_ref = (p/p_ref)^b`, and an independent
  chamber-pressure exponent. Exponents are calibration inputs, not universal
  engine constants.
- `PressurePerformanceTable` linearly interpolates absolute thrust, system Isp
  and chamber pressure between ordered hot-fire or validated-model rows. It
  refuses values outside the table.

NASA CR-122347 reports the tested device's thrust and Isp across supply
pressure and provides the empirical rationale for pressure-indexed maps. The
code does not copy that hardware's coefficients into a different engine.

For either model, total tank flow is recomputed from the system identity

```text
mdot_tank = F_delivered / (g0 Isp_system)
```

and every explicit propellant stream is scaled by the same flow ratio. Thus
stream closure, total tank drain and system Isp remain consistent even when
the rated point contains non-thrust-producing dump streams.

## Integration and convergence

The solver uses expelled propellant mass as the independent variable:

```text
dt/dm = 1 / mdot_tank(m)
dJ/dm = F / mdot_tank = g0 Isp_system(m)
```

Composite Simpson integration supplies the reference duration and impulse.
The same run independently creates a piecewise-linear
`rocket_propulsion_tabulated_burn_v1` artifact. The mass-domain mesh is halved
until the artifact's mass, impulse, duration change and thrust-centroid change
all meet the requested relative tolerance. `BlowdownIntegrationEvidence`
records the final step, refinement count and residuals.

## Example

```python
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    simulate_pressure_fed_blowdown,
)

rated = ConstantPerformanceProvider.monopropellant(
    delivered_thrust_n=100.0,
    system_specific_impulse_s=220.0,
    chamber_pressure_pa=1.5e6,
).operating_point()

tank = BlowdownTankModel(
    initial_pressure_pa=2.0e6,
    initial_ullage_volume_m3=0.01,
    initial_propellant_mass_kg=10.0,
    propellant_density_kg_m3=1000.0,
    polytropic_exponent=1.2,
)

law = PressurePerformanceLaw(
    reference_supply_pressure_pa=2.0e6,
    minimum_supply_pressure_pa=1.0e6,
    thrust_pressure_exponent=1.0,
    specific_impulse_pressure_exponent=0.02,
)

result = simulate_pressure_fed_blowdown(tank, rated, law)
print(result.cutoff_reason, result.artifact.artifact_hash)
```

For measured data, replace `law` with a `PressurePerformanceTable` of
`PressurePerformanceSample` rows. The integrator API is unchanged.

Sourced input distributions can be propagated with
`simulate_blowdown_uncertainty`. Each draw reruns this deterministic solver and
retains a complete exchange artifact; see
[`PROPULSION_UNCERTAINTY.md`](PROPULSION_UNCERTAINTY.md).

## Explicit limitations and refusal rules

The current model omits regulator dynamics, feed-line pressure loss and
compliance, pressurant heat transfer, vapor pressure, dissolved gas,
diaphragm/bladder mechanics, slosh, acceleration head, valve transients,
combustion instability and mixture-ratio migration. It must not be presented
as a flight-qualified digital twin.

The implementation refuses:

- non-physical tank, pressure, density or Isp values;
- pressure correlation/table extrapolation;
- unordered or duplicate table pressures;
- a reserve at or above initial inventory;
- an initial tank pressure outside the performance model range;
- a run that cannot meet its declared convergence tolerance.

## Primary references

- NASA SP-8112, *Pressurization Systems for Liquid Rockets*:
  <https://ntrs.nasa.gov/citations/19760015212>
- NASA CR-131400, *Reliability Model of a Monopropellant Auxiliary Propulsion
  System*: <https://ntrs.nasa.gov/citations/19730012094>
- NASA CR-122347, *Monopropellant Hydrazine Resisto Jet*:
  <https://ntrs.nasa.gov/citations/19720011120>
- NASA, *Transient Mathematical Modeling for Liquid Rocket Engine Systems*:
  <https://ntrs.nasa.gov/citations/20040000363>
- NASA Glenn, *Rocket Thrust Equation*:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/>
