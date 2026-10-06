"""Small numerical methods with explicit convergence evidence."""

from collections.abc import Callable
from dataclasses import dataclass

from .errors import ConvergenceError, DomainError, InputError


@dataclass(frozen=True, slots=True)
class RootResult:
    """Converged scalar root and the evidence used to accept it."""

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
    """Solve a bracketed monotonic scalar relation and report convergence."""

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

    The function may increase or decrease. The initial bounds must bracket the
    target. Bisection is slower than Newton iteration but predictable near sonic
    singularities and does not require analytical derivatives.

    Args:
        function: Continuous monotonic scalar function.
        target: Desired function value.
        lower: Lower search bound.
        upper: Upper search bound.
        tolerance: Relative/absolute convergence tolerance.
        max_iterations: Hard iteration limit.

    Raises:
        DomainError: If the target is not bracketed or convergence fails.
    """

    return solve_bisect_monotonic(
        function,
        target,
        lower,
        upper,
        tolerance=tolerance,
        max_iterations=max_iterations,
    ).value

