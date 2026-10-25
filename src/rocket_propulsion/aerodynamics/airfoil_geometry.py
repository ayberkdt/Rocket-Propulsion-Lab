"""NACA four- and five-digit airfoil section geometry.

Thickness form, mean lines and the combination of the two follow the defining
NACA reports [1]_ [2]_ [3]_ as summarised in NACA Report 824 [4]_. All lengths
are fractions of the chord ``c`` (multiply by the chord for metres); ``x = 0``
is the leading edge and ``x = 1`` the trailing edge of the mean line.

Thickness (TR-460, p. 4)::

    y_t = (t / 0.20) (0.29690 sqrt(x) - 0.12600 x - 0.35160 x^2
                      + 0.28430 x^3 - 0.10150 x^4)

with leading-edge radius ``r_t = 1.10 t^2`` (TR-460, p. 5). The formula
leaves a finite trailing-edge thickness ``2 y_t(1) = 0.0021 t / 0.10``; the
tabulated ordinates of TR-824 show the same value in parentheses.

Four-digit mean line (TR-460, p. 5), camber ``m`` at ``x = p``::

    y_c = m / p^2 (2 p x - x^2)                 (x <= p)
    y_c = m / (1 - p)^2 ((1 - 2p) + 2 p x - x^2)  (x >= p)

Five-digit mean lines (TR-537, p. 521-522, Table I), for the "simple" lines
210-250 with ``m`` and ``k1`` from Table I::

    y_c = k1 / 6 (x^3 - 3 m x^2 + m^2 (3 - m) x)   (x <= m)
    y_c = k1 m^3 / 6 (1 - x)                        (x >= m)

and for the reflexed lines 221-251 (also Table I), with
``k2/k1 = (3 (m - p)^2 - m^3) / (1 - m)^3`` (TR-537 p. 522; the Table I
column of this ratio is printed as 1e4 times that value)::

    y_c = k1 / 6 ((x - m)^3 - k2/k1 (1 - m)^3 x - m^3 x + m^3)        (x <= m)
    y_c = k1 / 6 (k2/k1 (x - m)^3 - k2/k1 (1 - m)^3 x - m^3 x + m^3)  (x >= m)

Table I is for a design lift coefficient of 0.3 (first digit 2). Per TR-824
(p. 5) the first digit ``L`` gives a design lift coefficient ``0.15 L`` and
all mean-line values scale linearly with it, so ``k1`` is multiplied by
``L / 2``. Upper and lower surfaces (TR-824 Eqs. 1-4)::

    x_U = x - y_t sin(theta),  y_U = y_c + y_t cos(theta)
    x_L = x + y_t sin(theta),  y_L = y_c - y_t cos(theta),  theta = atan(dy_c/dx)

References
----------
.. [1] NACA-TR-460. E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The
       Characteristics of 78 Related Airfoil Sections from Tests in the
       Variable-Density Wind Tunnel", NACA Report 460, 1933, pp. 4-5.
.. [2] NACA-TR-537. E. N. Jacobs, R. M. Pinkerton, "Tests in the
       Variable-Density Wind Tunnel of Related Airfoils Having the Maximum
       Camber Unusually Far Forward", NACA Report 537, 1936, pp. 521-522.
.. [3] NACA-TR-610. E. N. Jacobs, R. M. Pinkerton, H. Greenberg, "Tests of
       Related Forward-Camber Airfoils in the Variable-Density Wind Tunnel",
       NACA Report 610, 1937, p. 2.
.. [4] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
       "Summary of Airfoil Data", NACA Report 824, 1945, pp. 4-5, 100-101.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan, cos, isfinite, pi, sin, sqrt

from rocket_propulsion.core.errors import DomainError, InputError

# TR-537 Table I: mean-line designation -> (position of camber p, m, k1, tabulated k2/k1 column).
# The tabulated k2/k1 entries equal 1e4 times the ratio given by the TR-537 formula
# k2/k1 = (3 (m - p)^2 - m^3) / (1 - m)^3; the formula value is used for the geometry.
_FIVE_DIGIT_TABLE_I: dict[str, tuple[float, float, float, float | None]] = {
    "210": (0.05, 0.0580, 361.4, None),
    "220": (0.10, 0.1260, 51.64, None),
    "230": (0.15, 0.2025, 15.957, None),
    "240": (0.20, 0.2900, 6.643, None),
    "250": (0.25, 0.3910, 3.230, None),
    "221": (0.10, 0.1300, 51.99, 7.64),
    "231": (0.15, 0.2170, 15.793, 67.70),
    "241": (0.20, 0.3180, 6.520, 303.0),
    "251": (0.25, 0.4410, 3.191, 1355.0),
}

# TR-460 p. 4 basic thickness form coefficients for t = 0.20.
_THICKNESS_COEFFICIENTS = (0.29690, -0.12600, -0.35160, 0.28430, -0.10150)
_BASIC_THICKNESS = 0.20


@dataclass(frozen=True, slots=True)
class NacaSection:
    """Geometry of a NACA four- or five-digit section on a unit chord [1]_ [2]_ [3]_.

    Attributes
    ----------
    designation : str
        Four- or five-digit designation, e.g. ``"2412"`` or ``"23012"``.
    family : {"four-digit", "five-digit"}
        Series.
    thickness_ratio : float
        Maximum thickness ``t`` as a fraction of chord.
    max_camber : float
        Mean-line maximum ordinate: four-digit ``m``; for five-digit lines
        the ordinate ``y_c(p)`` at the tabulated position of maximum camber.
    camber_position : float
        Four-digit ``p``; five-digit position of maximum camber ``p``.
    junction : float
        Abscissa where the two mean-line arcs meet (``p`` for four-digit,
        ``m`` for five-digit).
    k1, k2_over_k1 : float
        Five-digit constants (``k1`` already scaled with the first digit);
        zero and ``None`` for four-digit sections.
    design_lift_coefficient : float or None
        Five-digit design lift coefficient ``0.15 L``; ``None`` otherwise.

    Notes
    -----
    Build instances with :func:`naca_section`; see the module docstring for
    the equations.

    References
    ----------
    .. [1] NACA-TR-460. E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The
           Characteristics of 78 Related Airfoil Sections from Tests in the
           Variable-Density Wind Tunnel", NACA Report 460, 1933, pp. 4-5.
    .. [2] NACA-TR-537. E. N. Jacobs, R. M. Pinkerton, "Tests in the
           Variable-Density Wind Tunnel of Related Airfoils Having the Maximum
           Camber Unusually Far Forward", NACA Report 537, 1936, Table I.
    .. [3] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
           "Summary of Airfoil Data", NACA Report 824, 1945, pp. 4-5.
    """

    designation: str
    family: str
    thickness_ratio: float
    max_camber: float
    camber_position: float
    junction: float
    k1: float = 0.0
    k2_over_k1: float | None = None
    design_lift_coefficient: float | None = None

    @property
    def leading_edge_radius(self) -> float:
        """Leading-edge radius ``1.10 t^2`` as a fraction of chord (TR-460 p. 5)."""

        return 1.10 * self.thickness_ratio**2

    @property
    def is_symmetric(self) -> bool:
        """True when the mean line is the chord line."""

        return self.max_camber == 0.0 and self.k1 == 0.0

    def thickness(self, x: float) -> float:
        """Half-thickness ``y_t(x)`` (TR-460 p. 4) for ``0 <= x <= 1``."""

        position = _station(x)
        a0, a1, a2, a3, a4 = _THICKNESS_COEFFICIENTS
        basic = a0 * sqrt(position) + position * (a1 + position * (a2 + position * (a3 + position * a4)))
        return self.thickness_ratio / _BASIC_THICKNESS * basic

    def mean_line(self, x: float) -> tuple[float, float]:
        """Return ``(y_c, dy_c/dx)`` of the mean line at ``x`` (``0 <= x <= 1``)."""

        position = _station(x)
        if self.family == "four-digit":
            m, p = self.max_camber, self.camber_position
            if m == 0.0:
                return 0.0, 0.0
            if position <= p:
                return m / p**2 * (2.0 * p * position - position**2), 2.0 * m / p**2 * (p - position)
            return (
                m / (1.0 - p) ** 2 * ((1.0 - 2.0 * p) + 2.0 * p * position - position**2),
                2.0 * m / (1.0 - p) ** 2 * (p - position),
            )
        m, k1, ratio = self.junction, self.k1, self.k2_over_k1
        if ratio is None:
            if position <= m:
                return (
                    k1 / 6.0 * (position**3 - 3.0 * m * position**2 + m * m * (3.0 - m) * position),
                    k1 / 6.0 * (3.0 * position**2 - 6.0 * m * position + m * m * (3.0 - m)),
                )
            return k1 * m**3 / 6.0 * (1.0 - position), -k1 * m**3 / 6.0
        tail = ratio * (1.0 - m) ** 3
        if position <= m:
            return (
                k1 / 6.0 * ((position - m) ** 3 - tail * position - m**3 * position + m**3),
                k1 / 6.0 * (3.0 * (position - m) ** 2 - tail - m**3),
            )
        return (
            k1 / 6.0 * (ratio * (position - m) ** 3 - tail * position - m**3 * position + m**3),
            k1 / 6.0 * (3.0 * ratio * (position - m) ** 2 - tail - m**3),
        )

    def upper_surface(self, x: float) -> tuple[float, float]:
        """Return ``(x_U, y_U)`` for mean-line station ``x`` (TR-824 Eqs. 1-2)."""

        y_t = self.thickness(x)
        y_c, slope = self.mean_line(x)
        theta = atan(slope)
        return x - y_t * sin(theta), y_c + y_t * cos(theta)

    def lower_surface(self, x: float) -> tuple[float, float]:
        """Return ``(x_L, y_L)`` for mean-line station ``x`` (TR-824 Eqs. 3-4)."""

        y_t = self.thickness(x)
        y_c, slope = self.mean_line(x)
        theta = atan(slope)
        return x + y_t * sin(theta), y_c - y_t * cos(theta)

    def panel_nodes(
        self, panels_per_surface: int, *, close_trailing_edge: bool = True
    ) -> tuple[tuple[float, float], ...]:
        """Return surface nodes clockwise from the lower trailing edge.

        Parameters
        ----------
        panels_per_surface : int
            Panels on each surface (total ``2 n``), at least 4.
        close_trailing_edge : bool, default True
            Close the finite TR-460 trailing-edge thickness by subtracting
            the linear ramp ``x y_t(1)`` from the half-thickness, so that the
            contour is closed while its trailing-edge angle is preserved.
            This is a library modification, not part of the NACA definition;
            it changes the ordinates by at most ``y_t(1) = 0.0105 t`` of
            chord (at the trailing edge). With ``False`` the TR-460
            ordinates are returned unchanged and the contour stays open.

        Returns
        -------
        tuple of (float, float)
            ``2 n + 1`` points: lower trailing edge, along the lower surface
            to the leading edge, along the upper surface to the upper trailing
            edge. Stations follow the cosine spacing
            ``x = (1 - cos(beta)) / 2`` with equal steps in ``beta``, which
            concentrates panels at both edges as recommended by Anderson
            (Section 4.10).
        """

        if panels_per_surface < 4:
            raise DomainError("At least 4 panels per surface are required.", field="panels_per_surface")
        stations = [0.5 * (1.0 - cos(pi * k / panels_per_surface)) for k in range(panels_per_surface + 1)]
        stations[-1] = 1.0
        ramp = self.thickness(1.0) if close_trailing_edge else 0.0

        def surface(x: float, side: float) -> tuple[float, float]:
            y_t = self.thickness(x) - ramp * x
            y_c, slope = self.mean_line(x)
            theta = atan(slope)
            return x - side * y_t * sin(theta), y_c + side * y_t * cos(theta)

        lower = [surface(x, -1.0) for x in reversed(stations)]
        upper = [surface(x, 1.0) for x in stations[1:]]
        return tuple(lower + upper)


def _station(x: float) -> float:
    position = float(x)
    if not isfinite(position) or position < 0.0 or position > 1.0:
        raise DomainError("Chord station must lie in [0, 1].", field="x")
    return position


def naca_section(designation: str) -> NacaSection:
    """Build a NACA four- or five-digit section from its designation.

    Parameters
    ----------
    designation : str
        Four digits ``MPTT`` (camber ``M`` percent at ``P`` tenths, thickness
        ``TT`` percent) or five digits ``LPQTT`` (design lift coefficient
        ``0.15 L``, mean line ``LPQ`` with ``Q = 0`` simple and ``Q = 1``
        reflexed, thickness ``TT`` percent). An optional ``"NACA"`` prefix and
        spaces are ignored.

    Returns
    -------
    NacaSection
        Section geometry.

    Raises
    ------
    InputError
        If the designation is malformed, has zero thickness, a four-digit
        camber without position, a mean line not listed in TR-537 Table I
        (including the 211 line, whose data TR-537 withholds), or a reflexed
        line with a first digit other than 2 (TR-537 tabulates reflexed lines
        only for a design lift coefficient of 0.3).

    Notes
    -----
    Numbering systems: TR-460 p. 4 (four-digit), TR-824 p. 5 and TR-537 p. 521
    (five-digit) [1]_ [2]_ [3]_.

    References
    ----------
    .. [1] NACA-TR-460. E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The
           Characteristics of 78 Related Airfoil Sections from Tests in the
           Variable-Density Wind Tunnel", NACA Report 460, 1933, p. 4.
    .. [2] NACA-TR-537. E. N. Jacobs, R. M. Pinkerton, "Tests in the
           Variable-Density Wind Tunnel of Related Airfoils Having the Maximum
           Camber Unusually Far Forward", NACA Report 537, 1936, Table I.
    .. [3] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
           "Summary of Airfoil Data", NACA Report 824, 1945, pp. 4-5.
    """

    if not isinstance(designation, str):
        raise InputError("NACA designation must be a string.", field="designation")
    digits = designation.upper().replace("NACA", "").replace(" ", "").replace("-", "")
    if not digits.isdigit() or len(digits) not in (4, 5):
        raise InputError("NACA designation must have four or five digits.", field="designation")
    thickness = int(digits[-2:]) / 100.0
    if thickness <= 0.0:
        raise InputError("NACA section thickness must be positive.", field="designation")
    if len(digits) == 4:
        camber = int(digits[0]) / 100.0
        position = int(digits[1]) / 10.0
        if (camber == 0.0) != (position == 0.0):
            raise InputError("Four-digit camber and camber position must both be zero or nonzero.",
                             field="designation")
        return NacaSection(digits, "four-digit", thickness, camber, position, position)
    lift_digit = int(digits[0])
    line = "2" + digits[1:3]
    if lift_digit == 0 or line not in _FIVE_DIGIT_TABLE_I:
        raise InputError(f"Mean line {digits[:3]} is not available from TR-537 Table I.",
                         field="designation")
    camber_position, junction, k1, ratio = _FIVE_DIGIT_TABLE_I[line]
    if ratio is not None and lift_digit != 2:
        raise InputError("Reflexed five-digit mean lines are tabulated only for a first digit of 2.",
                         field="designation")
    scaled_k1 = k1 * lift_digit / 2.0
    if ratio is not None:
        ratio = (3.0 * (junction - camber_position) ** 2 - junction**3) / (1.0 - junction) ** 3
    section = NacaSection(
        digits,
        "five-digit",
        thickness,
        0.0,
        camber_position,
        junction,
        scaled_k1,
        ratio,
        0.15 * lift_digit,
    )
    camber_at_position = section.mean_line(camber_position)[0]
    return NacaSection(
        digits, "five-digit", thickness, camber_at_position, camber_position, junction,
        scaled_k1, ratio, 0.15 * lift_digit,
    )
