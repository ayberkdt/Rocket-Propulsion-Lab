# Multi-burn propulsion timelines

`propulsion/burns/timeline.py` places several propagator bridges on one clock
and adds discrete mass jettison (stage separation, fairing, adapters). It
answers the questions a propagator would otherwise have to hard-code: which
burn is firing now, what all its roots are, how much mass leaves at a
separation, and how to persist the whole sequence as one verified document.

## Building a timeline

```python
timeline = PropulsionTimeline(
    burns=(
        ScheduledBurn("stage-1", stage_1_bridge, ignition_time_s=0.0),
        ScheduledBurn("stage-2", stage_2_bridge, ignition_time_s=165.0),
        ScheduledBurn("stage-2-restart", stage_2_restart_bridge, ignition_time_s=3600.0),
    ),
    jettisons=(
        MassJettison("stage-1-dry", 160.0, 4_200.0, "first-stage separation"),
        MassJettison("fairing", 190.0, 900.0, "payload fairing"),
    ),
)
```

Construction checks:

- burn and jettison identifiers are unique, and burn identifiers contain no `/`;
- burns are sorted by ignition time;
- burns may not overlap, though touching at one instant is allowed, unless
  `allow_overlapping_burns=True`. With that flag, concurrent forces and drains
  are summed, for example attitude thrusters during a main burn;
- burns may share tank states, so a restartable stage drains one inventory.

## Evaluation

```python
value = timeline.evaluate(
    time_s=t,
    vehicle_mass_kg=m,
    directions={"stage-2": u_stage_2},     # only active burns need a direction
    additional_states=tank_masses,
    stopped_burns=stopped,                 # burns the consumer stopped at a root
    ambient_pressure_pa=p_atm,             # forwarded only to bridges with an exit area
    operating_conditions={"stage-2": {...}},
    parameter_overrides={"stage-2": {"thrust_scale": 1.01}},
    side="right",
    apply_inventory_limits=False,
)
```

The result carries the summed acceleration, mass rate, per-tank drains and
\(\partial\mathbf a/\partial m\), plus every individual
`PropulsionDynamicsEvaluation` with its full partials. Estimator drivers stay
per burn, so one burn's `thrust_scale` is never confused with another's. A
non-overlapping timeline that finds two burns active raises an error rather
than summing them.

## Events

`timeline.events()` returns:

- **time events**, sorted: each burn's start, boundaries and stop (ids
  `"<burn>/burn:start"` …) and each jettison (`"jettison/<id>"`);
- **state events**: each burn's dry-mass and reserve roots. When one fires,
  the consumer adds the burn to `stopped_burns`.

At a jettison event the consumer calls `timeline.apply_jettison(event_id, m)`,
which refuses a non-positive remaining mass. `jettisoned_mass_kg(t)` supports
mass-closure checks across separations. Time events use nominal ignition
times; a non-zero `ignition_time_bias_s` driver shifts that burn's events by
the same amount.

## Persistence

`to_json()` / `from_json()` store the schedule, the jettisons and every
bridge's complete document. Reopening verifies each artifact hash, each
surface hash, each bridge's `model_hash` and finally the `timeline_hash`.
Changing any ignition time, jettison mass or profile sample is detected.

## Verification

`tests/test_burn_timeline.py` integrates a two-burn sequence with an
interstage jettison between the burns. The velocity matches
\(c\ln(m_0/m_1)+c\ln\big((m_1-m_J)/m_2\big)\) to 1e-10, and tank and vehicle
masses close exactly. Further tests cover touching burns (left and right
limits select the correct burn), summed concurrent burns, per-burn ambient
pressure forwarding, stopped burns and tamper detection.

## References

- Concurrent forces and drains add linearly (Newton's second law with summed
  thrust; mass rate is the sum of propellant flows): NASA Glenn, *Rocket
  Thrust Equation*,
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/>.
- The jettison oracle in the tests (rocket equation across burns with a
  discrete mass drop): NASA Glenn, *Ideal Rocket Equation*,
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>.
- Event roots and restart semantics: Orekit 13.1.5 `EventDetectorsProvider`
  API,
  <https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html>.
