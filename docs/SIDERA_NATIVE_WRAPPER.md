# Native Sidera wrapper for the propagation bridge

`rocket_propulsion.integrations.sidera_native` flies a
`PropagatorPropulsionBridge` burn through Sidera's public numerical propagator.
It does not modify Sidera and does not reduce the burn to an equivalent
constant Isp. Sidera stays optional: nothing is imported until a
`SideraNativeBurn` is built, and an absent or incompatible installation raises
`FeatureUnavailableError`.

## Sidera public API used

Checked against the Sidera checkout at `D:\sidera` (version `0.2.1.dev0`).

| Sidera member | Location | Use |
|---|---|---|
| `propagate(dynamics, y0, cfg, *, time_cfg, extra_events, ...)` | `src/sidera/core/propagation/propagator.py:175` | one call per coast or burn segment |
| `PropagationDynamics` protocol | `src/sidera/core/propagation/dynamics_protocol.py:45` | the segment dynamics satisfy it; the base engine must too |
| `RestartPlan` | `src/sidera/core/propagation/dynamics_protocol.py:27` | forwarded unchanged from the base engine |
| `maneuver_frame_basis(state, frame)` | `src/sidera/core/propagation/plans.py:681` | resolves `inertial` / `ric` / `vnb` directions from the state |
| extra-event attributes `terminal`, `direction`, `uses_full_state` | `src/sidera/core/propagation/event_wiring.py:236-254`, `src/sidera/core/events.py:214` | dry-mass terminal event on `y[6]` |

The tests also use `build_body_dynamics` (`src/sidera/core/body_engine/orbits.py:714`)
and `load_body_config` (`src/sidera/core/config.py:798`) to build a point-mass
Earth engine.

## How the bridge maps onto Sidera

- **Thrust.** Each profile segment is propagated with a
  `SideraBridgeSegmentDynamics`. It forwards every protocol member to the base
  engine and adds the bridge's net axial `F / m` to `dy[3:6]` and the bridge's
  total mass rate to `dy[6]`. The direction is resolved from the current state
  on Sidera's maneuver axes at every evaluation, and the bridge is evaluated at
  the absolute Sidera time with `evaluate_at_epoch`.
- **Start, boundary and stop.** These become segment boundaries at the
  canonical dates `ignition + bias + t_k`, each rounded once. This is the same
  way Sidera flies its own `ProfiledBurn` (`propagator.py:872-905`). The right
  limit is taken at a segment start and the left limit at a segment end, so no
  accepted step crosses a jump or kink.
- **Dry mass.** `protected_dry_mass_kg` becomes a terminal Sidera extra event
  `g = y[6] - m_dry` with `direction = -1` and `uses_full_state = True`. When
  it fires, thrust stops at the located root and the run coasts on the base
  engine. With this event registered the bridge is evaluated with
  `apply_inventory_limits=False`, so the clamp does not bias the steps.
- **Capability negotiation.** `inspect_sidera_native_capabilities()` checks the
  four required entry points: `propagate`, the dynamics protocol,
  `extra_events` and the frame basis. `require_sidera_native_capabilities()`
  also asks for `additional_states_available` when the bridge has stream
  bindings or reserves, `operating_conditions_available` when it has a
  performance surface, and `ambient_pressure_available` when it has a nozzle
  exit area.

## What cannot be done natively, and why

- **Tank states and tank reserves.** Sidera's `propagate` accepts only the
  6-state or the 7-state with mass (`src/sidera/core/propagation/result.py:12-26`).
  The body engine accepts nothing else either (`orbits.py:92-99`), and Sidera
  has no public `PropulsionProvider` with additional states. Bridges with
  stream bindings or `AdditionalStateConstraint` reserves are refused. The
  vehicle-level reserve (`protected_dry_mass_kg`) is supported.
- **Performance surfaces and ambient pressure.** Sidera has no public channel
  that passes operating conditions or ambient pressure to a thrust provider, so
  these bridges are refused.
- **A restart-plan event instead of segmentation.** Sidera's restart solver
  (`scipy_runner.py:236`) ends a segment at a root and keeps the same
  right-hand side across segments. A right-hand side evaluated exactly at a
  burn boundary could then be either the end of one step or the start of the
  next, and it cannot tell which one-sided limit to use without keeping state.
  Separate segments remove that ambiguity. This is reasoning from the source;
  the restart-plan path was not implemented or tested.
- **Analytic partials.** The bridge's partials are not passed on, because
  Sidera's `propagate` has no variational channel for a user force.

## Tests and tolerances

`tests/test_sidera_native.py` puts Sidera on `sys.path` only when
`RPLAB_SIDERA_SRC` is set. It also sets `sys.dont_write_bytecode` and
redirects `NUMBA_CACHE_DIR` to the system temp directory, so that neither
bytecode nor Numba JIT caches are written into the checkout. The refusal tests run in both modes. The
end-to-end tests fly the two-segment example burn: 6 kg of a 100 kg vehicle,
constant `F/q = 50 m/s`, ignition at t = 10 s. The orbit is circular at
350 000 km altitude, which keeps the gravity-gradient difference against a
coast reference negligible.

| Check | Tolerance | Measured | Meaning |
|---|---:|---:|---|
| final mass vs `m0 - consumed_propellant_kg` | 1e-10 kg | 1.4e-14 kg | Sidera integrates the bridge's exact drain |
| `|v_burn - v_coast|` vs `free_space_burn_metrics` Δv | 1e-8 m/s | 4.4e-10 m/s | delivered free-space Δv. The residual mixes integrator error (`rtol=1e-10`) and the gravity-gradient term; the two are not separated |
| Δv along +x vs `50 ln(m0/mf)` | 1e-8 m/s | 4.4e-10 m/s | rocket equation, exact here because `F/q` is constant |
| dry-mass root time vs `t_ign + 1 + 0.5/3` | 1e-9 s | 3.6e-15 s | terminal-event location |

Burn start, boundary and stop are checked at their exact dates. The mass is
also checked to stay constant before ignition and after stop, and the ignition
bias is checked to shift all three time events.

## References

- MIT-522-L5, thrust acceleration F/m (equation 19):
  https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf
- NASA-7009B, model credibility and verification reporting:
  https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf
- Sidera public API (`propagate`, `PropagationDynamics`, `maneuver_frame_basis`)
  read from the local checkout; file and line locators are listed above.
