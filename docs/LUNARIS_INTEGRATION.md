# Lunaris integration

[Lunaris](https://github.com/ayberkdt/lunaris) is a lunar-orbit propagation
framework. It propagates a 6- or 7-element state (mass at `y[6]`) and supports
impulsive maneuvers: `ImpulsiveManeuver(t_burn_s, dv_mps, frame, isp_s)` in a
strictly time-ordered `ManeuverPlan`, with frames `inertial` and `ric` and the
ideal-rocket mass update \(m\,e^{-|\Delta v|/(g_0 I_{sp})}\). It accepts any
dynamics object that provides `build_rhs()`.

`rocket_propulsion.integrations.lunaris` connects Rocket Propulsion Lab
timelines to Lunaris **without modifying Lunaris**. Lunaris is an optional
import; without it, the impulsive equivalents are still produced as data. The
adapter was written against Lunaris `0.1.0rc3.dev0` (`lunaris.core.propagation`).

## Path 1 — exact-mass impulsive equivalents

```python
from rocket_propulsion.integrations import (
    lunaris_impulsive_equivalents, to_lunaris_maneuver_plan,
)

impulses = lunaris_impulsive_equivalents(
    timeline,
    initial_mass_kg=1000.0,
    directions={"lower": ((0.0, -1.0, 0.0), "ric"), "raise": ((0.0, 1.0, 0.0), "ric")},
)
plan = to_lunaris_maneuver_plan(impulses)   # lunaris ManeuverPlan
```

For each burn, taken in time order and starting from the previous burn's final
mass:

- \(\Delta v\) is the burn's free-space velocity change (Gauss–Legendre, with
  an error estimate);
- the impulse is placed at the acceleration-weighted centroid \(t^\ast\);
- \(I_{sp}=c_\mathrm{eff}/g_0\) with \(c_\mathrm{eff}=\Delta v/\ln(m_0/m_f)\).
  Lunaris' mass law then reproduces the finite burn's final mass **exactly**
  (both sides use \(g_0 = 9.80665\) m/s²).

What is approximated: the burn's time distribution, i.e. the finite-duration /
gravity loss. Each `LunarisImpulsiveEquivalent` reports
`finite_burn_duration_s`, so the size of that approximation is visible.

## Path 2 — finite burns inside Lunaris

```python
from rocket_propulsion.integrations import (
    LunarisFiniteBurnDynamics, lunaris_segmentation_plan, velocity_direction,
)

dynamics = LunarisFiniteBurnDynamics(
    lunaris_dynamics, timeline, velocity_direction,
    parameter_overrides={"raise": {"ignition_time_bias_s": 2.0}},  # optional
    initial_mass_kg=1000.0,           # required when a bridge has protected_dry_mass_kg
)
result = lunaris.core.propagation.propagate(
    dynamics, y0_with_mass, cfg, time_cfg=time_cfg,
    maneuver_plan=dynamics.segmentation_plan(),
)
```

Build the segmentation plan from the wrapper (or pass the same
`parameter_overrides` and `initial_mass_kg` to `lunaris_segmentation_plan`).
A plan built from nominal times while a burn is biased would leave its
discontinuities inside integration steps; a test shows the resulting loss of
accuracy.

- The wrapper adds the timeline's thrust acceleration and mass rate to the
  wrapped right-hand side. It overrides both `build_rhs` and Lunaris' trusted
  `_build_solver_rhs`; without the latter, Lunaris would bypass the
  propulsion. All other attributes (`grav`, `flags`, …) are delegated.
- `lunaris_segmentation_plan` places zero-Δv impulses (no `isp_s`, so no mass
  change) at every burn discontinuity. Lunaris restarts integration at
  maneuver times, so each segment is smooth.
- Each segment builds a fresh right-hand side. Its first evaluation takes the
  right-hand limit at a discontinuity; later evaluations on a boundary take
  the left-hand limit.
- `direction_law(burn_id, t, y)` returns the inertial thrust direction
  (`velocity_direction` = prograde).
- **Ignition bias.** Every time event of a burn is shifted by its resolved
  `ignition_time_bias_s`. A bias that makes two burns overlap is refused.
- **Exact dry-mass cut.** Lunaris event functions only see `y[:6]`, so a
  mass-threshold event cannot be located by Lunaris. For a mass-only burn the
  mass flow is a function of time alone, so with `initial_mass_kg` the time at
  which the vehicle reaches `protected_dry_mass_kg` is solved exactly
  (piecewise-quadratic mass, stable quadratic root) and becomes a segment
  boundary; the burn is off from that time on (`lunaris_burn_windows` reports
  it). The bridge's internal clamp is then switched off by default
  (`apply_inventory_limits=None`), because a clamp evaluated inside
  Runge–Kutta stages next to the cut biases the final mass; a test shows that
  bias.
- **Mass consistency.** With `initial_mass_kg`, every segment's first
  evaluation compares `y[6]` with the mass predicted from the timeline
  (`predicted_mass_kg`, relative tolerance `mass_tolerance=1e-9`) and fails
  closed on a mismatch: a different start mass or another mass-changing
  maneuver would make the precomputed cut wrong.

Fail-closed rules: Lunaris carries no tank states, so bridges must be
mass-only (no stream bindings or reserve constraints; `protected_dry_mass_kg`
still works, with `initial_mass_kg`). Performance surfaces, ambient pressure
and jettison events have no Lunaris counterpart and are refused, and
6-element states are rejected. Impulsive equivalents refuse a burn that its
protected dry mass would truncate; use the finite-burn path for it.

## Verification

`tests/test_lunaris_integration.py`:

- impulsive equivalents reproduce the finite burn's \(\Delta v\) and, through
  Lunaris' mass law, its final mass to 1e-14 relative;
- the finite-burn wrapper, integrated with Lunaris-style segmentation in free
  space, reproduces the rocket equation to 1e-10;
- **end-to-end with the real Lunaris propagator** (skipped when Lunaris is not
  installed): a two-burn sequence around a point-mass Moon, propagated once as
  finite burns and once as impulses. The two paths agree on final mass, and
  their semi-major-axis changes agree within 0.1 % of the change.
- bias-shifted segmentation reproduces the unbiased free-space result to
  1e-10, while nominal segmentation of a biased burn is measurably worse;
- the dry-mass cut time matches closed forms inside a ramp and in the steady
  phase (1e-12 s); the finite burn ends exactly at the dry mass with the
  rocket-equation velocity (1e-10);
- **end-to-end with the real Lunaris propagator**: a biased 600 s burn cut by
  its dry mass after 250 s matches an explicitly truncated 250 s burn
  (position to 1 mm, velocity to 1 µm/s), and a wrong start mass is refused.

## Not covered

- Lunaris' GPU / batch backends (the wrapper targets the CPU `propagate`).
- Tank-level inventory inside Lunaris (use the propagator bridge directly in a
  propagator that carries extra states).
- Burn targeting: computing impulses or finite burns that reach a target orbit
  is the mission layer's job.

## References

- Ideal rocket equation \(\Delta v = c\ln(m_0/m_f)\), used for the exact-mass
  impulse (\(I_{sp}=c_\mathrm{eff}/g_0\)) and the rocket-equation oracles:
  NASA Glenn, *Ideal Rocket Equation*,
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>.
- \(g_0 = 9.80665\) m/s² (exact, conventional): BIPM SI Brochure, 9th ed.,
  2019, <https://www.bipm.org/en/publications/si-brochure>; registered in
  `rocket_propulsion.core.constants`.
- Gauss–Legendre quadrature for \(\Delta v\) and \(t^\ast\): NIST DLMF §3.5(v),
  <https://dlmf.nist.gov/3.5#v>.
- Event and discontinuity handling (restart at each discontinuity, left/right
  limits): Orekit 13.1.5 `EventDetector` API,
  <https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html>.
- Dry-mass cut time: piecewise-linear flow gives a piecewise-quadratic mass;
  the root \(\tau = 2c/(q_0+\sqrt{q_0^2+4ac})\) is the cancellation-free form
  of the quadratic formula (D. Goldberg, *What Every Computer Scientist Should
  Know About Floating-Point Arithmetic*, ACM Computing Surveys 23(1), 1991,
  §1.4; bibliographic citation).
- Lunaris maneuver and event contracts: Lunaris source,
  `lunaris.core.propagation` (version `0.1.0rc3.dev0`),
  <https://github.com/ayberkdt/lunaris>.
