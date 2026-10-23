# Engine-cycle component balance review

This is an idealized component sizing model. It does not close a real engine,
select motor/feed hardware, demonstrate a relevant environment or establish TRL.
The full review and all nine hardware evidence gates remain incomplete.

## Routed mass and power contract

The previous default gas-generator model charged the pump for chamber flow only,
while adding generator flow to tank drain. NASA TM X-64749 uses the whole pumped
oxidizer/fuel flows in its pump-power balance. The new default is a declared
single equivalent liquid, common cold-side head and shared pump. This is a
simplified specialization, not an implementation/validation of the entire report.

Let m be chamber flow, s=DeltaP/(rho*eta_p), and w be available turbine work per
mass. Shared-pump GG satisfies Ppump=(m+g)*s=g*w, hence g=m*s/(w-s), total=m*w/(w-s),
generator fraction=s/w and delivered fraction=1-s/w. Positive duty requires w>s.
For staged return flow, total pumped flow remains m, g=m*s/w and w>=s is allowed;
returning turbine flow is not counted again as tank drain. Zero head gives zero
shaft/turbine demand and unchanged chamber delivery.

`external-gg-feed` explicitly retains the old main-pump-only component basis:
g=m*s/w, pump flow=m, tank drain=m+g. Its generator pressure source and energy
cost are excluded, identified in energy_balance_scope and a critical warning.
This is not claimed to be a closed complete generator cycle.

Turbine work is eta_t*cp*T0*(1-r**(-(gamma-1)/gamma)); r is total inlet pressure
divided by total outlet pressure. Constant cp/gamma and component efficiencies
are idealizations. Hot-gas pressure compatibility, separate pump/mixture maps,
cavitation, mechanical losses and integrated hardware capability remain unclosed.
The zero power residual is an algebraic identity, not independent measured closure.

## Arithmetic and capacity comparisons

Cold-side required pressure and head retain exact binary rational ledgers before
separate export. Coefficient products are not exported before cancellation,
multiplication or division; this repairs erased head and recoverable pump/work
underflow or overflow. The existing directed Decimal/checked positive export
kernels enclose log, exp and dependent ratios. Each output rounds separately.
The finite precision schedule can fail with ConvergenceError, which is not proof
of physical infeasibility or unrepresentability.

Rational powers receive an exact path: for reduced p/q, a rational base has a
rational p/q power only if its coprime numerator and denominator are perfect
q-th powers. Integer-root searches are bounded by bit length. This recognizes
exact shaft-capacity equalities; no absolute feasibility tolerance is used.
Otherwise directed work/head bounds certify which side of capacity applies.
The same unrounded work determines flows and delivery, including close feedback
boundaries. The 8% generator warning is a heuristic on the exported coefficient,
not an operating limit, empirical validation or acceptance criterion.

## Pressure-fed sizing, records and API

Pressure-fed mode sizes the injector-feed requirement and compares the supplied
inlet pressure exactly. Insufficient supply returns a critical warning and false
pressure_fed_supply_sufficient, even when rounded displayed pressures coincide.
Pump/turbine inputs can be absent in this mode. Supplied inputs still receive
strict scalar validation. Pumped modes require their active inputs.

The additional signed inlet margin is optional when binary64 cannot represent
its nonzero value: null plus inlet_pressure_margin_representable=false preserves
that limitation without restricting otherwise representable pressure-fed sizing.
The sufficiency decision still uses the exact sign. Mandatory physical outputs
remain checked finite quantities, with typed errors for erased nonzero exports.

Results retain frozen EngineCycleInputs, pump routing, numerical policy, scope,
assumptions and warnings. Every field replays from source inputs; JSON input and
warning mappings are copied and validated. Replay is internal consistency, not
data/campaign authentication. API metadata advances to 2.0; default GG results
change because generator pump duty is included. The UI disables inactive fields,
exposes external feed and pressure shortfall, and preserves model/policy versions.

## Evidence and limits

Tests compare independent Decimal equations at two precisions, exact capacity
fixtures, near-unity gas factors, recoverable products, routing, malformed input,
forged records and JSON ownership. The first repair preserved a JSON warning
reopen failure before fixing the reader. Instrumented Node executions cover the
actual control/render functions, not browser layout or all frontend flows.
An additional mapping protocol attack exposed bare KeyError in both retained
input and warning conversion; these failures are archived. Known malformed
mapping failures now become typed domain errors while unrelated RuntimeError
callbacks still propagate. The reference index also retains evidence of an
encoding/index failure and an exact recovery against the Git source plus its
known added citation; this was documentation repair, not a numerical model fix.
Source-bound failed baselines and controlled-fault reports remain locally under
ignored output/assurance. These cannot replace selected requirements, independent
operating data, measurement uncertainty or relevant-environment demonstration.

## References

- NASA TM X-64749, Schmucker, 1973, pump-power balance and whole pumped flows:
  <https://ntrs.nasa.gov/citations/19730016059>.
- NASA Glenn, EngineSim theory, ideal turbine temperature/pressure and shaft work:
  <https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf>.
- Python 3.12, correctly rounded Decimal functions and directed contexts:
  <https://docs.python.org/3.12/library/decimal.html>.
