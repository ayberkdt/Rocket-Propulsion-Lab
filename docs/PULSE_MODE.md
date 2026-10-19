# Pulse-mode thrusters

Attitude-control and small correction thrusters fire in short electrical
pulses. For short pulses the delivered impulse is not `F_steady · t_on`:
valve opening, chamber fill and tail-off make the **impulse bit** and the
pulse-averaged specific impulse nonlinear functions of the on-time. The
valve also sets a smallest pulse, the **minimum impulse bit** (MIB). A
trajectory or station-keeping tool therefore cannot command an arbitrary
small Δv. It gets an integer number of pulses, and their total generally
differs from the request.

`PulseModeThruster` (`propulsion/burns/pulse_mode.py`) represents exactly
that.

## Model

Inputs are user-supplied tables intended to characterize hardware. The class
does not authenticate their measurement provenance or empirical validity.
The test tables are synthetic. The library provides no generic physical pulse
coefficients.

| Field | Meaning |
| --- | --- |
| `on_times_s` | electrical pulse widths, strictly increasing; the first one is the minimum on-time |
| `impulse_bits_n_s` | impulse per pulse, strictly increasing; the first one is the MIB |
| `propellant_per_pulse_kg` | propellant per pulse, positive and non-decreasing |
| `steady_thrust_n`, `steady_mass_flow_kg_s` | fully developed rates beyond the last node |
| `maximum_on_time_s` | longest permitted pulse (default: last node) |
| `minimum_off_time_s` | shortest gap between pulses in a train |

Between nodes, `I(t)` and `m(t)` use PCHIP (monotone, no overshoot).
`pulse(t)` returns their exact derivatives. Beyond the last node, extra
on-time adds steady impulse and flow linearly. The value is continuous
there; the slope may change.

## Requests and quantization

A train of `n` identical pulses has one exhaust velocity `c = I/m`, so the
rocket equation holds exactly:

\[
  \Delta v_n(t) = \frac{I(t)}{m(t)}\,\ln\frac{m_0}{m_0 - n\,m(t)} .
\]

In ideal arithmetic, \(\Delta v_n\) increases strictly with `t` and with `n`.

Budget-limited width inversion evaluates only the propellant curve. A model
owns at most 512 cached scalar inverse solutions, keyed by the complete
per-pulse mass budget; a short lock protects cache access without retaining
models globally. Full pulse values are reused only within one planning call.
Caches do not enter the model manifest, hash, equality or JSON exchange.
Tests interleave budgets, initial masses, policies and distinct models against
enumerated feasible intervals, and preserve normal copying and persistence.
These optimizations do not change search bounds or request tolerances.

- `plan_impulse(I)` and `plan_delta_v(Δv, m0, propellant_budget_kg=...)`
  return the **fewest** equal
  pulses that meet the request. The pulse count comes from bisection on `n`,
  and the on-time from bisection on `t`.
- A request is met exactly when some `n` brackets it between `n` minimum
  bits and `n` longest pulses. Otherwise it lies in a gap between `n − 1`
  longest pulses and `n` minimum bits (always the case below one MIB). The
  plan is then `quantized=True` and reports `residual = requested −
  delivered`.
  - `quantization="nearest"` picks the side with the smaller residual.
  - `quantization="no_overshoot"` never delivers more than requested.
- `plan_delta_v` needs a propellant budget below the initial mass: the dry
  structure is not propellant. For `n` pulses the feasible widths are
  `[t_min, t_max(n)]` with `t_max(n) = min(t_longest, m⁻¹(budget/n))`. Beyond
  the count the budget allows at the longest width, shorter pulses can still
  fit, and those budget-limited counts are scanned as well.
- A plan whose delivered value is within 1e-9 (relative) of the request
  counts as exact (`quantized=False`), so a request one rounding step outside
  an achievable interval is met by that interval's train rather than by one
  more pulse.
- Requests above everything the budget can deliver, or above the pulse limit
  (applied to the returned plan), fail closed.

The tests compare plans with brute-force enumeration: no integer train with
equal pulses gets closer to the request, and exact plans use the minimum
count.

## Propagator handoff

`pulse_train_artifact(plan, period_s)` returns a `TabulatedBurnArtifact`.
Each pulse is a rectangle over its on-time, with thrust `I/t_on` and flow
`m/t_on`. Impulse, propellant and exhaust velocity per pulse are therefore
exact; the pulse's internal shape is not represented. The period must leave
the minimum off-time.

The artifact drives the propagator bridge, timelines and the Lunaris
finite-burn wrapper like any other burn. A test shows the free-space Δv
through the Lunaris wrapper equal to the plan's Δv within 1e-12 m/s.

## Not modelled

- Pulse-to-pulse variation and thermal history (a cold first pulse differs
  from a warm train). Use measured tables for the relevant condition, and the
  execution-error models for scatter.
- Off-pulsing and thruster-pair torque coupling. Use `ThrusterLayout` for
  force and torque allocation.

## Reference

- J. M. Parker and M. J. Wilson, *The Minimum Impulse Thruster*, 2005,
  NTRS 20070032005: <https://ntrs.nasa.gov/citations/20070032005>. The record
  was located through the NTRS search index; its full text could not be
  opened from the development environment.
- F. N. Fritsch and R. E. Carlson, *Monotone Piecewise Cubic Interpolation*,
  SIAM J. Numer. Anal. 17(2), 1980 (cited bibliographically).
