"""Shared shape-preserving interpolation kernels.

References
----------
F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
SIAM Journal on Numerical Analysis 17(2), 1980.
F. N. Fritsch and J. Butland, A Method for Constructing Local Monotone
Piecewise Cubic Interpolants, SIAM J. Sci. Stat. Comput. 5(2), 1984.
"""

from __future__ import annotations

from itertools import pairwise


def pchip_slopes(x: tuple[float, ...], y: tuple[float, ...]) -> tuple[float, ...]:
    """Return Fritsch-Carlson / Fritsch-Butland monotone nodal slopes.

    Interior slopes are the weighted harmonic mean of adjacent secants (zero
    at a local extremum); end slopes use the shape-preserving three-point
    formula.  The resulting cubic Hermite interpolant is C1 and monotone
    wherever the data are.

    References
    ----------
    F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis 17(2), 1980.
    """

    count = len(x)
    widths = [b - a for a, b in pairwise(x)]
    secants = [(y[i + 1] - y[i]) / widths[i] for i in range(count - 1)]
    if count == 2:
        return (secants[0], secants[0])
    slopes = [0.0] * count
    for i in range(1, count - 1):
        left, right = secants[i - 1], secants[i]
        if left * right <= 0.0:
            slopes[i] = 0.0
        else:
            w1 = 2.0 * widths[i] + widths[i - 1]
            w2 = widths[i] + 2.0 * widths[i - 1]
            slopes[i] = (w1 + w2) / (w1 / left + w2 / right)

    def end_slope(h0: float, h1: float, d0: float, d1: float) -> float:
        slope = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if slope * d0 <= 0.0:
            return 0.0
        if d0 * d1 <= 0.0 and abs(slope) > abs(3.0 * d0):
            return 3.0 * d0
        return slope

    slopes[0] = end_slope(widths[0], widths[1], secants[0], secants[1])
    slopes[-1] = end_slope(widths[-1], widths[-2], secants[-1], secants[-2])
    return tuple(slopes)

