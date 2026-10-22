"""Exact declared-power relay decisions and observable power-grid root targets.

Required powers are positive, strictly increasing binary-input fractions supplied
by validated built-in laws. Hysteresis and bus margin are non-negative. These
are numerical policy rules, not empirical power acceptance margins.
"""

import sys
from fractions import Fraction
from math import nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer

_POWER_MODE_POLICY = "exact_binary_budget_hysteresis_and_power_grid_targets_v1"
_MAX_POWER = Fraction(sys.float_info.max)


def _mode_index(current, count, *, allow_none=False):
    if current is None and allow_none:
        return None
    current = integer(current, "Current throttle level")
    if not -1 <= current < count:
        raise DomainError(
            "Current throttle level is outside the law.", field="Current throttle level"
        )
    return current


def _select_mode(power, required, hysteresis, current):
    """None initializes without memory; known modes follow the declared relay."""
    current = _mode_index(current, len(required), allow_none=True)
    feasible = -1
    for index, demand in enumerate(required):
        if demand <= power:
            feasible = index
    if current is None or current >= 0 and power < required[current]:
        return feasible
    level = current
    while level + 1 < len(required) and power >= required[level + 1] + hysteresis:
        level += 1
    return level


def _ceil_power(threshold):
    """First finite input power at or above a reachable positive threshold."""
    if not 0 < threshold <= _MAX_POWER:
        raise DomainError(
            "Power threshold has no finite upward guard.", field="Power-switch threshold"
        )
    rounded = float(threshold)
    return nextafter(rounded, float("inf")) if Fraction(rounded) < threshold else rounded


def _below_power(threshold):
    """Last finite input strictly below a reachable positive down threshold."""
    if not 0 < threshold <= _MAX_POWER:
        raise DomainError(
            "Power threshold has no finite downward guard.", field="Power-switch threshold"
        )
    rounded = float(threshold)
    return nextafter(rounded, -float("inf")) if Fraction(rounded) >= threshold else rounded


def _mode_roots(required, hysteresis, current, identifiers, margin, *, saturation=None):
    """Retain exact sums; coalesce transitions crossed by one visible power step.

    Impossible startup/saturation roots are omitted. A requested on-mode whose
    sustaining threshold exceeds every finite bus power is an invalid event state.
    """
    current = _mode_index(current, len(required))
    result = []
    if current >= 0:
        down_threshold = required[current] + margin
        if down_threshold > _MAX_POWER:
            raise DomainError(
                "Current mode cannot be sustained by finite bus power.",
                field="Current throttle level",
            )
        guard = _below_power(down_threshold)
        down_target = _select_mode(Fraction(guard) - margin, required, hysteresis, current)
        result.append(
            (
                f"power:{identifiers[current]}:down",
                down_threshold,
                "decreasing",
                current,
                down_target,
            )
        )
    if current + 1 < len(required):
        up_threshold = required[current + 1] + hysteresis + margin
        if up_threshold <= _MAX_POWER:
            guard = _ceil_power(up_threshold)
            up_target = _select_mode(Fraction(guard) - margin, required, hysteresis, current)
            result.append(
                (
                    f"power:{identifiers[current + 1]}:up",
                    up_threshold,
                    "increasing",
                    current,
                    up_target,
                )
            )
    if saturation is not None and current == 0:
        kink = saturation + margin
        if kink <= _MAX_POWER:
            result.append((f"power:{identifiers[0]}:saturation", kink, "either", current, current))
    return tuple(result)
