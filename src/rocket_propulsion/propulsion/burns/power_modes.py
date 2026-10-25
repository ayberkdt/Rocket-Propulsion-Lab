"""Exact power-relay mode selection and power-switch root targets.

Required powers are positive, strictly increasing exact rationals of the binary
inputs supplied by validated built-in throttle laws [1]_; hysteresis and bus
margin are non-negative.  Mode decisions compare exact rationals, and each
event target is the mode selected at the first representable input power
beyond the root [2]_, so levels crossed within one float step coalesce into a
single directed event [3]_.  Hysteresis and margin are numerical policy rules,
not empirical power acceptance margins.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
       their exact binary value as a rational.
.. [2] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter": adjacent
       representable floating-point values.
.. [3] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
       EventDetector" API.
"""

import sys
from fractions import Fraction
from math import nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer

_POWER_MODE_POLICY = "exact_binary_budget_hysteresis_and_power_grid_targets_v1"
_MAX_POWER = Fraction(sys.float_info.max)


def _mode_index(current, count, *, allow_none=False):
    """Return a validated mode index in ``[-1, count)``, or ``None`` when allowed."""
    if current is None and allow_none:
        return None
    current = integer(current, "Current throttle level")
    if not -1 <= current < count:
        raise DomainError(
            "Current throttle level is outside the law.", field="Current throttle level"
        )
    return current


def _select_mode(power, required, hysteresis, current):
    """Return the relay mode for an exact effective power.

    Parameters
    ----------
    power : Fraction
        Effective power, in W.
    required : tuple of Fraction
        Strictly increasing required powers of the modes, in W.
    hysteresis : Fraction
        Non-negative up-switch hysteresis, in W.
    current : int or None
        Current mode (``-1`` is off).  ``None`` initializes without memory.

    Returns
    -------
    int
        Selected mode index, or ``-1`` (off).

    Raises
    ------
    DomainError
        If ``current`` is not ``None`` and lies outside ``[-1, len(required))``.

    Notes
    -----
    Without memory, or when ``power`` falls strictly below the current mode's
    requirement, the highest mode with ``required[k] <= power`` is returned.
    Otherwise the mode rises while ``power >= required[k + 1] + hysteresis``.
    All comparisons are exact [1]_.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational.
    """
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
    """Return the first float at or above a positive exact threshold.

    Raises
    ------
    DomainError
        If the threshold is not positive or exceeds the largest finite float.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter": adjacent
           representable floating-point values.
    """
    if not 0 < threshold <= _MAX_POWER:
        raise DomainError(
            "Power threshold has no finite upward guard.", field="Power-switch threshold"
        )
    rounded = float(threshold)
    return nextafter(rounded, float("inf")) if Fraction(rounded) < threshold else rounded


def _below_power(threshold):
    """Return the last float strictly below a positive exact threshold.

    Raises
    ------
    DomainError
        If the threshold is not positive or exceeds the largest finite float.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter": adjacent
           representable floating-point values.
    """
    if not 0 < threshold <= _MAX_POWER:
        raise DomainError(
            "Power threshold has no finite downward guard.", field="Power-switch threshold"
        )
    rounded = float(threshold)
    return nextafter(rounded, -float("inf")) if Fraction(rounded) >= threshold else rounded


def _mode_roots(required, hysteresis, current, identifiers, margin, *, saturation=None):
    """Return the directed power-switch roots of mode ``current``.

    Parameters
    ----------
    required : tuple of Fraction
        Strictly increasing required powers, in W.
    hysteresis : Fraction
        Up-switch hysteresis, in W.
    current : int
        Current mode; ``-1`` is off.
    identifiers : tuple of str
        Mode identifiers used in event names.
    margin : Fraction
        Bus margin added to every threshold, in W.
    saturation : Fraction, optional
        Saturation power of a continuous law, which adds a mode-preserving
        ``either`` kink root while mode ``0`` is active.

    Returns
    -------
    tuple of tuple
        ``(event_id, exact_threshold, crossing, from_level, to_level)`` items.

    Raises
    ------
    DomainError
        If ``current`` is invalid or an active mode's sustaining threshold
        exceeds every finite bus power.

    Notes
    -----
    Thresholds are exact sums.  Each target is the mode selected at the
    first representable power beyond the root [1]_ [2]_, so transitions
    crossed within one float step coalesce into one event [3]_.  Unreachable
    start-up and saturation roots are omitted.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational.
    .. [2] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter": adjacent
           representable floating-point values.
    .. [3] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
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
