"""Small numerical methods that report explicit convergence evidence."""

from collections.abc import Callable
from dataclasses import dataclass

from .errors import ConvergenceError, DomainError, InputError


@dataclass(frozen=True, slots=True)
class RootResult:
    """Accepted scalar root together with the evidence used to accept it.

    Attributes
    ----------
    value : float
        Accepted abscissa (the last bisection midpoint).
    iterations : int
        Number of function evaluations at midpoints.
    residual : float
        ``abs(function(value) - target)`` at the accepted point.
    lower_bound, upper_bound : float
        Bracket that contained ``value`` when it was accepted.
    tolerance : float
        Tolerance used in the acceptance test.
    method : str
        Name of the root-finding method.
    """

    value: float
    iterations: int
    residual: float
    lower_bound: float
    upper_bound: float
    tolerance: float
    method: str = "bisection"


def solve_bisect_monotonic(
    function: Callable[[float], float],
    target: float,
    lower: float,
    upper: float,
    *,
    tolerance: float = 1.0e-12,
    max_iterations: int = 200,
) -> RootResult:
    """Solve a bracketed monotonic scalar equation by bisection.

    Parameters
    ----------
    function : callable
        Continuous, monotonic (increasing or decreasing) function of one float.
    target : float
        Function value to reach.
    lower, upper : float
        Initial bracket; ``lower`` must be strictly smaller than ``upper``.
    tolerance : float, optional
        Acceptance tolerance, applied to the function residual and scaled by
        ``max(1, abs(target))``.
    max_iterations : int, optional
        Maximum number of midpoint evaluations.

    Returns
    -------
    RootResult
        Accepted root, residual, final bracket and iteration count.

    Raises
    ------
    InputError
        If ``tolerance`` is not positive, ``max_iterations`` is below one or
        the bracket is empty or reversed.
    DomainError
        If ``target`` lies outside the function values at the bracket ends,
        widened by ``tolerance``.
    ConvergenceError
        If no midpoint satisfies the acceptance test within
        ``max_iterations`` evaluations.

    Notes
    -----
    At each step the bracket is halved and the half that still contains the
    target is kept, using the sign of ``function(midpoint) - target`` and the
    monotonic direction inferred from the end values. A midpoint is accepted
    when ``abs(function(midpoint) - target) <= tolerance * max(1, abs(target))``.
    The test is on the residual, not on the bracket width, so a flat function
    can be accepted over a wide bracket. The method needs no derivatives and
    stays predictable near sonic singularities, at the cost of linear
    convergence.
    """

    if tolerance <= 0.0:
        raise InputError("Tolerance must be greater than zero.", field="tolerance")
    if max_iterations < 1:
        raise InputError(
            "Maximum iterations must be at least one.", field="max_iterations"
        )
    if lower >= upper:
        raise InputError(
            "Lower bound must be smaller than upper bound.",
            code="invalid_bracket",
            details={"lower": lower, "upper": upper},
        )

    low_value = function(lower)
    high_value = function(upper)
    minimum, maximum = sorted((low_value, high_value))
    if not minimum - tolerance <= target <= maximum + tolerance:
        raise DomainError(
            f"Target {target:g} is outside the bracketed range "
            f"[{minimum:g}, {maximum:g}].",
            code="target_not_bracketed",
            details={
                "target": target,
                "range_minimum": minimum,
                "range_maximum": maximum,
                "lower_bound": lower,
                "upper_bound": upper,
            },
        )

    increasing = high_value >= low_value
    residual = min(abs(low_value - target), abs(high_value - target))
    for iteration in range(1, max_iterations + 1):
        midpoint = (lower + upper) / 2.0
        value = function(midpoint)
        residual = abs(value - target)
        if residual <= tolerance * max(1.0, abs(target)):
            return RootResult(
                value=midpoint,
                iterations=iteration,
                residual=residual,
                lower_bound=lower,
                upper_bound=upper,
                tolerance=tolerance,
            )

        if (value < target) == increasing:
            lower = midpoint
        else:
            upper = midpoint

    raise ConvergenceError(
        "Bisection did not converge within the iteration limit.",
        details={
            "method": "bisection",
            "iterations": max_iterations,
            "residual": residual,
            "lower_bound": lower,
            "upper_bound": upper,
            "tolerance": tolerance,
        },
    )


def bisect_monotonic(
    function: Callable[[float], float],
    target: float,
    lower: float,
    upper: float,
    *,
    tolerance: float = 1.0e-12,
    max_iterations: int = 200,
) -> float:
    """Return ``x`` such that a monotonic ``function(x)`` equals ``target``.

    Convenience wrapper around :func:`solve_bisect_monotonic` that returns
    only the accepted abscissa.

    Parameters
    ----------
    function : callable
        Continuous, monotonic (increasing or decreasing) scalar function.
    target : float
        Desired function value.
    lower, upper : float
        Initial bracket; it must contain the target.
    tolerance : float, optional
        Residual tolerance, scaled by ``max(1, abs(target))``.
    max_iterations : int, optional
        Maximum number of midpoint evaluations.

    Returns
    -------
    float
        The accepted root.

    Raises
    ------
    DomainError
        If the arguments are invalid (:class:`InputError`), the target is not
        bracketed, or bisection does not converge (:class:`ConvergenceError`);
        all of these derive from :class:`DomainError`.

    See Also
    --------
    solve_bisect_monotonic : Same solve, with convergence evidence.
    """

    return solve_bisect_monotonic(
        function,
        target,
        lower,
        upper,
        tolerance=tolerance,
        max_iterations=max_iterations,
    ).value
