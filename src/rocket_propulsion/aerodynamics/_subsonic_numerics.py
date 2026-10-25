"""Private numerical helpers for the inviscid subsonic aerodynamics modules.

Dense Gaussian elimination and Gauss-Legendre quadrature, written with the
standard library only so that the package keeps no third-party dependency.
These are generic numerical-analysis tools (any numerical-methods text); the
aerodynamic equations they serve are cited in the calling modules.

References
----------
.. [1] DLMF-3.5. NIST Digital Library of Mathematical Functions, section 3.5(v),
       Gauss quadrature; Gauss-Legendre nodes are the zeros of the Legendre
       polynomials, tabulated in Table 3.5.1.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from functools import lru_cache
from math import cos, isfinite, pi

from rocket_propulsion.core.errors import ConvergenceError, DomainError


def solve_dense(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Solve a dense square linear system by Gaussian elimination.

    Parameters
    ----------
    matrix : sequence of sequence of float
        Square coefficient matrix, row major.
    rhs : sequence of float
        Right-hand side with one entry per row.

    Returns
    -------
    list of float
        Solution vector.

    Raises
    ------
    DomainError
        If the shapes do not match or an entry is not finite.
    ConvergenceError
        If a pivot is zero to working precision (singular system).

    Notes
    -----
    Partial (row) pivoting; the singularity test compares each pivot with
    ``1e-14`` times the largest absolute entry of the original matrix. Generic
    numerical linear algebra; no aerodynamic source is implied.
    """

    size = len(rhs)
    if size == 0 or len(matrix) != size or any(len(row) != size for row in matrix):
        raise DomainError("Linear system must be square and non-empty.")
    work = [[float(value) for value in row] + [float(rhs[index])] for index, row in enumerate(matrix)]
    scale = max(abs(value) for row in work for value in row[:size])
    if not all(isfinite(value) for row in work for value in row) or scale == 0.0:
        raise DomainError("Linear system entries must be finite and not all zero.")
    for column in range(size):
        pivot_row = max(range(column, size), key=lambda row: abs(work[row][column]))
        if abs(work[pivot_row][column]) <= 1.0e-14 * scale:
            raise ConvergenceError("Linear system is singular to working precision.")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
        pivot = work[column]
        pivot_tail = pivot[column:]
        for row in range(column + 1, size):
            factor = work[row][column] / pivot[column]
            if factor != 0.0:
                target = work[row]
                target[column:] = [a - factor * b for a, b in zip(target[column:], pivot_tail, strict=True)]
    solution = [0.0] * size
    for row in range(size - 1, -1, -1):
        total = work[row][size] - sum(work[row][k] * solution[k] for k in range(row + 1, size))
        solution[row] = total / work[row][row]
    return solution


@lru_cache(maxsize=8)
def gauss_legendre(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return Gauss-Legendre nodes and weights on ``[-1, 1]``.

    Parameters
    ----------
    order : int
        Number of nodes, at least 2.

    Returns
    -------
    nodes, weights : tuple of float
        Abscissae (zeros of the Legendre polynomial ``P_order``) and weights.

    Raises
    ------
    ConvergenceError
        If the Newton iteration for a node does not converge.

    Notes
    -----
    Newton iteration on ``P_n`` from the Chebyshev-type initial guess
    ``cos(pi (i - 1/4) / (n + 1/2))``; weights ``2 / ((1 - x^2) P_n'(x)^2)``.
    The nodes are the zeros of ``P_n`` [1]_; the 5-point nodes and weights are
    checked against DLMF Table 3.5.1 in the test suite.

    References
    ----------
    .. [1] DLMF-3.5. NIST Digital Library of Mathematical Functions, section
           3.5(v), Gauss quadrature; Gauss-Legendre nodes are the zeros of the
           Legendre polynomials, tabulated in Table 3.5.1.
    """

    if order < 2:
        raise DomainError("Gauss-Legendre order must be at least 2.")
    nodes: list[float] = []
    weights: list[float] = []
    for index in range(1, order + 1):
        x = cos(pi * (index - 0.25) / (order + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for degree in range(2, order + 1):
                p0, p1 = p1, ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree
            derivative = order * (x * p1 - p0) / (x * x - 1.0)
            step = p1 / derivative
            x -= step
            if abs(step) < 1.0e-15:
                break
        else:  # pragma: no cover - the iteration converges for every order used here
            raise ConvergenceError("Gauss-Legendre node iteration did not converge.")
        p0, p1 = 1.0, x
        for degree in range(2, order + 1):
            p0, p1 = p1, ((2 * degree - 1) * x * p1 - (degree - 1) * p0) / degree
        derivative = order * (x * p1 - p0) / (x * x - 1.0)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return tuple(nodes), tuple(weights)


def integrate(
    function: Callable[[float], float],
    breakpoints: Sequence[float],
    *,
    panels_per_interval: int = 16,
    order: int = 16,
) -> float:
    """Integrate a function piecewise with composite Gauss-Legendre rules.

    Parameters
    ----------
    function : callable
        Integrand ``f(t)``; must be smooth between consecutive breakpoints.
    breakpoints : sequence of float
        Increasing interval ends, including both outer limits.
    panels_per_interval : int, default 16
        Equal sub-panels per breakpoint interval.
    order : int, default 16
        Gauss-Legendre nodes per sub-panel.

    Returns
    -------
    float
        The integral.

    Notes
    -----
    Splitting at breakpoints keeps the rule spectrally accurate when the
    integrand has a derivative discontinuity there (for example the junction
    of the two NACA mean-line arcs) [1]_.

    References
    ----------
    .. [1] DLMF-3.5. NIST Digital Library of Mathematical Functions, section
           3.5(v), Gauss quadrature; Gauss-Legendre nodes are the zeros of the
           Legendre polynomials, tabulated in Table 3.5.1.
    """

    nodes, weights = gauss_legendre(order)
    total = 0.0
    for start, end in zip(breakpoints[:-1], breakpoints[1:], strict=True):
        if end < start:
            raise DomainError("Quadrature breakpoints must be increasing.")
        width = (end - start) / panels_per_interval
        for panel in range(panels_per_interval):
            left = start + panel * width
            half = 0.5 * width
            middle = left + half
            total += half * sum(
                weight * function(middle + half * node)
                for node, weight in zip(nodes, weights, strict=True)
            )
    return total
