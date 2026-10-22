# Electric mode, event-coordinate and frontier review

This bounded review follows the electric document/arithmetic review. It covers
the two built-in laws' declared power budgets, hysteresis, model bus margin,
switching roots, and discrete thrust/Isp frontier ordering. These are software
policy contracts for supplied values. Mission requirements, measured switching
behavior, PPU response, validated motor maps and relevant-environment hardware
demonstrations remain absent. Passing this review cannot establish TRL 5.

## Reproduced failures and declared policy

An absolute 1e-9 W comparison floor previously allowed accepted positive-demand
synthetic levels to operate at zero power. Rounding margin or hysteresis sums
before comparisons also admitted underfunded modes. Required power, hysteresis
and bus-margin operands now retain the exact rational values of supplied binary
scalars. No absolute power allowance is added; measurement uncertainty and
mission-derived engineering margins require separately justified inputs.

Initialization with current_level=None selects the highest bare feasible mode.
A known OFF state requires its startup hysteresis. A held mode drops to the
highest bare feasible mode only when power is strictly below its demand;
upward eligibility includes equality at demand plus hysteresis. This deliberate
distinction between initialization and memory is exercised across both laws.
OFF evaluation need not export an unrepresentable unused effective power.
Held continuous evaluation retains exact margin through saturation decisions.

## Root contract and consumer obligations

PowerSwitchSurface.threshold_ratio stores a normalized pair of integer operands;
threshold_exact_w exposes its Fraction value. threshold_w is a float display.
value(power) subtracts the true threshold before checked float export. For a
1 W demand and 1e16 W margin, the display is 1e16 but the residual there is -1 W.
Consumers must call value(), rather than subtracting the displayed threshold.

Upward targets use the first representable power at or above the exact root;
downward targets use the last representable power strictly below it. Targets
may skip intermediate levels that cannot be retained at that input coordinate.
Unreachable finite-input startup/saturation roots are omitted. Requesting an
ON mode whose sustaining demand exceeds every finite bus input fails typed.
Residual overflow fails typed rather than exporting infinity to a solver.

The consumer must apply crossing direction and to_level, handle initialization
at zero, reset detectors, and avoid interpreting a constant-zero interval as a
new crossing. At downward equality the mode is still feasible; the event target
describes the strictly lower side. The supplied sampled-controller trace checks
both directions, equality, repeated samples, mode sequence, zero deficits and
propellant accounting. It does not qualify a generic solver's root clock,
backward propagation, simultaneous hardware events or controller firmware.
SciPy documents that sign-change detection can miss multiple roots per step;
power-coordinate coalescing does not remove that consumer responsibility.

## Frontier and persisted model identity

Isp ordering compares exact thrust/flow ratios; the common positive standard
gravity cancels. Distinct ratios with identical displayed Isp remain ordered,
and ordering need not export an unrepresentable Isp. At equal power, the better
objective wins. This ranks declared points without authenticating their source.

The model manifest requires power_mode_numerics with value
exact_binary_budget_hysteresis_and_power_grid_targets_v1, included in its content
hash. Old payloads without the marker are rejected explicitly; reopening them
requires a deliberate migration/reconstruction and changes their model hash.
No legacy payload is silently assigned the new switching semantics. Subclass
callbacks retain public boundaries; this is not an adequacy review of custom laws.

## Evidence and oracle correction

The independent Decimal oracle enumerates feasible destinations with 1700-digit
precision. Cases include neighboring binary inputs, the smallest subnormal,
fractional hysteresis, tiny and large margins, every held state, finite-range
limits, actual residual signs, invalid directions/ratios, metadata identity and
two exact Isp counterexamples. Each test has a stated public behavior to falsify.
Curated mutations exercise selected faults; they are not exhaustive coverage.

The initial targeted run had one incorrect fixture: zero hysteresis made the
first upward threshold exactly representable, so expecting a two-level jump
was unjustified. Its original source and failure are archived under
output/assurance/history/electric_mode_initial_oracle_failure. The fixture now
uses 0.25 W hysteresis so both roots share the same first visible power guard.
The jump assertion remains; an extra preceding-input assertion checks that the
earlier guard cannot activate the higher mode. No numerical limits or test
execution budgets were relaxed. Final source-bound counts belong to the
verification reports, not this prospective review document.

## References

- Python 3.12 fractions: exact values of float inputs and normalized rational
  arithmetic. <https://docs.python.org/3.12/library/fractions.html>
- Python 3.12 math.nextafter: adjacent representable inputs for visible guards.
  <https://docs.python.org/3.12/library/math.html#math.nextafter>
- Python 3.12 decimal: independent precision contexts and exact float conversion.
  <https://docs.python.org/3.12/library/decimal.html>
- SciPy 1.14.1 solve_ivp: directed sign-change event detection and its within-step
  limitation. This source does not endorse the local relay policy.
  <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.solve_ivp.html>
