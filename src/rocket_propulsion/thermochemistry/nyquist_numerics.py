"""Resolve a Nyquist contour using supplied global amplitude and local slope bounds.

The positive-frequency contour starts at zero. Each accepted interval is
contained in a disk about its midpoint which excludes the origin: the radius
is an upper bound on |dF/domega| times the half width. The curve and its chord
are then homotopic without crossing zero, so the chord has the same phase
increment. The omitted tail lies in |F-1| < 1 and is closed there explicitly.

This is floating-point numerical enclosure, conditional on the supplied
bounds and no open-loop RHP poles. It is not an interval-arithmetic proof or
a validation of those assumptions for a measured transfer function.
"""

from __future__ import annotations

import cmath
from collections.abc import Callable
from itertools import pairwise
from math import fsum, hypot, isfinite, pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real


def resolve_contour(*, response: Callable[[complex], complex],
                    slope_bound: Callable[[float, float], float],
                    amplitude_bound: float, residence_time_s: float,
                    maximum_frequency_hz: float | None,
                    maximum_evaluations: int) -> tuple[int, tuple[float, ...], tuple[float, ...], float, tuple[float, ...]]:
    """Count closed-loop RHP zeros or refuse an unresolved/invalid contour.

    ``amplitude_bound / hypot(1, omega*residence_time_s)`` bounds |f| for
    *all* nonnegative frequencies, including the uncomputed tail.
    ``slope_bound(a,b)`` bounds |df/domega| throughout the interval [a,b].
    """
    if (isinstance(maximum_evaluations, bool) or not isinstance(maximum_evaluations, int)
            or not 1000 <= maximum_evaluations <= 1_000_000):
        raise DomainError("points must be an integer from 1000 to 1000000.", field="points")
    scale = amplitude_bound
    # Above this frequency the whole remaining contour lies in |F-1| <= 1/2.
    twice = 2 * scale
    half_tail = twice * sqrt(max(0.0, 1 - (1 / twice) ** 2)) / residence_time_s if scale >= 1 else 0.0
    requested = 2 * pi * maximum_frequency_hz if maximum_frequency_hz is not None else max(
        100 * pi / residence_time_s, half_tail)
    end = min(requested, half_tail) if half_tail else requested
    if not isfinite(end) or end <= 0:
        raise DomainError("Nyquist frequency range cannot be represented finitely.", field="maximum_frequency_hz")
    tail_radius = scale / hypot(1.0, end * residence_time_s)
    if tail_radius >= 1 - 1e-12:
        raise DomainError("The global admittance envelope does not close the Nyquist tail; "
                          "raise maximum_frequency_hz.", code="nyquist_tail_unresolved",
                          field="maximum_frequency_hz", details={"tail_radius_bound": tail_radius})
    cache: dict[float, complex] = {}

    def evaluate(omega: float) -> complex:
        if omega not in cache:
            if len(cache) >= maximum_evaluations:
                raise DomainError("The Nyquist contour is unresolved within the evaluation budget; "
                                  "raise points. No stability verdict is available.",
                                  code="nyquist_resolution", field="points",
                                  details={"evaluations": len(cache), "budget": maximum_evaluations})
            try:
                value = complex(response(1j * omega))
            except DomainError:
                raise
            except (TypeError, ValueError, ArithmeticError) as exc:
                raise DomainError("Admittance evaluation failed on the Nyquist contour.",
                                  code="nyquist_response", details={"omega_rad_s": omega}) from exc
            if not (isfinite(value.real) and isfinite(value.imag)):
                raise DomainError("Admittance must be finite on the Nyquist contour.", code="nyquist_response")
            envelope = scale / hypot(1.0, omega * residence_time_s)
            if abs(value) > envelope * (1 + 1e-12) + 1e-14:
                raise DomainError("The response exceeds its declared amplitude envelope.",
                                  code="nyquist_invalid_bound", details={"omega_rad_s": omega})
            cache[omega] = value
        return cache[omega]

    start_value = 1 + evaluate(0.0)
    if abs(start_value.imag) > 1e-12 * max(1.0, abs(start_value)):
        raise DomainError("A real-system Nyquist count needs a real DC response.", code="nyquist_response")
    if abs(start_value) <= 1e-12 * max(1.0, scale):
        raise DomainError("A zero at or too close to the contour prevents a stability verdict.",
                          code="nyquist_boundary_unresolved")
    end_value = 1 + evaluate(end)
    if scale < 1:
        # Small-gain disk excludes zero at every frequency, including infinity.
        return 0, (), (), 1 - scale, (0.0, end / (2 * pi))

    stack = [(0.0, end)]
    phase_steps: list[float] = []
    distance = 1 - tail_radius
    while stack:
        left, right = stack.pop()
        middle = left + (right - left) / 2
        if middle in (left, right):
            raise DomainError("Nyquist contour cannot be resolved at floating-point precision.",
                              code="nyquist_resolution", field="points")
        a, b, m = 1 + evaluate(left), 1 + evaluate(right), 1 + evaluate(middle)
        try:
            slope = real(slope_bound(left, right), "admittance_derivative_bound")
        except DomainError:
            raise
        except (TypeError, ValueError, ArithmeticError) as exc:
            raise DomainError("The admittance derivative bound cannot be evaluated finitely.",
                              code="nyquist_invalid_bound") from exc
        if slope < 0:
            raise DomainError("Derivative bounds must be nonnegative.", code="nyquist_invalid_bound")
        radius = slope * ((right - left) / 2)
        # Check observations against the claimed bound, but do not pretend
        # finite samples can validate an arbitrary callback's global contract.
        slack = 1e-12 * max(1.0, abs(a), abs(b), abs(m))
        if max(abs(a - m), abs(b - m)) > radius * (1 + 1e-12) + slack:
            raise DomainError("The response violates its derivative bound.", code="nyquist_invalid_bound")
        if isfinite(radius) and radius + slack < 0.75 * abs(m):
            step = cmath.phase(b) - cmath.phase(a)
            step = (step + pi) % (2 * pi) - pi
            phase_steps.append(step)
            distance = min(distance, max(0.0, abs(m) - radius - slack))
        else:
            stack.append((middle, right))
            stack.append((left, middle))
    # The terminal disk is in Re(F)>0. Close to F=1 within this disk;
    # reflect the contour across the real axis, allowing negative F(0).
    count = -(fsum(phase_steps) - cmath.phase(end_value)) / pi
    zeros = round(count)
    if zeros < 0 or abs(count - zeros) > 1e-7:
        raise DomainError("Nyquist winding is inconsistent with the model assumptions.",
                          code="nyquist_resolution")

    crossings: list[float] = []
    margins: list[float] = []
    grid = sorted(cache)
    for left, right in pairwise(grid):
        va, vb = abs(cache[left]) - 1, abs(cache[right]) - 1
        if va * vb >= 0:
            continue
        # Resolve each bracketed gain crossover instead of assigning a grid
        # endpoint's phase to its geometric-mean frequency.
        for _ in range(64):
            middle = left + (right - left) / 2
            vm = abs(evaluate(middle)) - 1
            if vm == 0 or right - left <= 1e-8 * max(1.0, middle):
                break
            if va * vm < 0:
                right = middle
            else:
                left, va = middle, vm
        crossings.append(middle / (2 * pi))
        margin = cmath.phase(evaluate(middle)) * 180 / pi + 180
        margins.append(margin - 360 if margin > 180 else margin)
    return zeros, tuple(crossings), tuple(margins), distance, tuple(w / (2 * pi) for w in sorted(cache))
