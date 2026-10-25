"""Private numerical helpers shared by the advanced compressible-flow modules.

The helpers are generic root bracketing and unimodal maximisation routines.
They carry no physics; every physical relation lives in the public modules
with its own source.
"""

from __future__ import annotations

from collections.abc import Callable
from math import isfinite, sqrt

from rocket_propulsion.core.errors import ConvergenceError, DomainError

_GOLDEN = (sqrt(5.0) - 1.0) / 2.0


def bracketed_root(
    function: Callable[[float], float],
    lower: float,
    upper: float,
    *,
    x_tolerance: float = 1.0e-14,
    max_iterations: int = 200,
    label: str = "root",
) -> float:
    """Return a sign-change root of ``function`` inside ``[lower, upper]``.

    Parameters
    ----------
    function : callable
        Continuous scalar function of one float.
    lower, upper : float
        Bracket with ``lower < upper``; the function must change sign on it.
    x_tolerance : float, optional
        Bracket-width tolerance, scaled by ``max(1, |midpoint|)``.
    max_iterations : int, optional
        Maximum number of halvings.
    label : str, optional
        Name used in error messages.

    Returns
    -------
    float
        Midpoint of the final bracket.

    Raises
    ------
    DomainError
        If the end values are not finite or do not bracket a sign change.
    ConvergenceError
        If a non-finite value is met inside the bracket.

    Notes
    -----
    Plain bisection on the sign of the function. The acceptance test is on
    the bracket width, so the result does not depend on the magnitude of the
    residual.
    """

    f_lower = function(lower)
    f_upper = function(upper)
    if not (isfinite(f_lower) and isfinite(f_upper)):
        raise DomainError(f"{label}: non-finite value at the bracket ends.")
    if f_lower == 0.0:
        return lower
    if f_upper == 0.0:
        return upper
    if (f_lower > 0.0) == (f_upper > 0.0):
        raise DomainError(
            f"{label}: the bracket does not contain a sign change.",
            code="root_not_bracketed",
            details={"lower": lower, "upper": upper, "f_lower": f_lower, "f_upper": f_upper},
        )
    for _ in range(max_iterations):
        midpoint = 0.5 * (lower + upper)
        if upper - lower <= x_tolerance * max(1.0, abs(midpoint)):
            return midpoint
        value = function(midpoint)
        if not isfinite(value):
            raise ConvergenceError(f"{label}: non-finite value inside the bracket.")
        if value == 0.0:
            return midpoint
        if (value > 0.0) == (f_lower > 0.0):
            lower, f_lower = midpoint, value
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def maximize_unimodal(
    function: Callable[[float], float],
    lower: float,
    upper: float,
    *,
    iterations: int = 90,
) -> tuple[float, float]:
    """Golden-section search for the maximum of a unimodal function.

    Parameters
    ----------
    function : callable
        Scalar function with a single interior maximum on the interval.
    lower, upper : float
        Search interval.
    iterations : int, optional
        Number of golden-section reductions (each shrinks the interval by
        about 0.618).

    Returns
    -------
    tuple of float
        ``(x_max, f_max)``.
    """

    a, b = lower, upper
    c = b - _GOLDEN * (b - a)
    d = a + _GOLDEN * (b - a)
    fc, fd = function(c), function(d)
    for _ in range(iterations):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - _GOLDEN * (b - a)
            fc = function(c)
        else:
            a, c, fc = c, d, fd
            d = a + _GOLDEN * (b - a)
            fd = function(d)
    x = 0.5 * (a + b)
    return x, function(x)
