"""Two-dimensional inviscid, incompressible potential flow.

Elementary flows (uniform stream, source, doublet, vortex) and their
superposition; the nonlifting and lifting flow over a circular cylinder; the
Kutta-Joukowski theorem; the incompressible pressure coefficient; and the exact
flow over a Joukowski section obtained by conformal mapping of the lifting
cylinder flow.

Sign conventions follow Anderson [1]_: a source of strength ``Lambda`` has
``V_r = Lambda / (2 pi r)``; a vortex of strength ``Gamma`` is positive
**clockwise** (``V_theta = -Gamma / (2 pi r)``), so a positive circulation on a
body in a stream from left to right produces positive (upward) lift
``L' = rho V Gamma``; the doublet of strength ``kappa`` is the one of Eq.
(3.87), whose arrow points in the negative x direction. All quantities are SI;
angles are in radians.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
       6th ed., McGraw-Hill, 2017, Sections 3.5 and 3.9-3.16.
.. [2] NACA-TR-411. T. Theodorsen, "Theory of Wing Sections of Arbitrary
       Shape", NACA Report 411, 1931.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import asin, atan2, cos, exp, hypot, isfinite, log, pi, sin, sinh, sqrt
from typing import Protocol

from rocket_propulsion.core.errors import DomainError


def _finite(value: float, label: str, *, positive: bool = False) -> float:
    """Return ``value`` as float after a finite (and optional positive) check.

    Raises
    ------
    DomainError
        If the value is not finite, or not strictly positive when required.
    """

    number = float(value)
    if not isfinite(number):
        raise DomainError(f"{label} must be finite.", field=label)
    if positive and number <= 0.0:
        raise DomainError(f"{label} must be greater than zero.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class PlanarVelocity:
    """Cartesian velocity of a planar flow at one point [1]_.

    Attributes
    ----------
    u, v : float
        Components along x and y in m/s.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.9.
    """

    u: float
    v: float

    @property
    def speed(self) -> float:
        """Magnitude ``sqrt(u^2 + v^2)`` in m/s."""

        return hypot(self.u, self.v)


class ElementaryFlow(Protocol):
    """Interface shared by the elementary flows (velocity, potential, stream function)."""

    def velocity(self, x: float, y: float) -> PlanarVelocity: ...

    def potential(self, x: float, y: float) -> float: ...

    def stream_function(self, x: float, y: float) -> float: ...


def _offset(x: float, y: float, x0: float, y0: float) -> tuple[float, float, float]:
    """Return ``(dx, dy, r^2)`` relative to a singular point, rejecting the point itself."""

    dx = _finite(x, "x") - x0
    dy = _finite(y, "y") - y0
    radius_squared = dx * dx + dy * dy
    if radius_squared == 0.0:
        raise DomainError("The flow is singular at the location of the elementary flow.")
    return dx, dy, radius_squared


@dataclass(frozen=True, slots=True)
class UniformFlow:
    """Uniform stream of speed ``V`` inclined at ``angle_rad`` to the x axis [1]_.

    Attributes
    ----------
    speed : float
        Freestream speed ``V`` in m/s, finite and positive.
    angle_rad : float, default 0
        Flow direction measured counterclockwise from the x axis.

    Notes
    -----
    For ``angle_rad = 0`` the potential and stream function are Eqs. (3.53)
    and (3.55), ``phi = V x`` and ``psi = V y``; an inclined stream is the same
    flow in a rotated frame: ``phi = V (x cos a + y sin a)``,
    ``psi = V (y cos a - x sin a)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.9, Eqs. (3.53)-(3.57).
    """

    speed: float
    angle_rad: float = 0.0

    def __post_init__(self) -> None:
        _finite(self.speed, "speed", positive=True)
        _finite(self.angle_rad, "angle_rad")

    def velocity(self, x: float, y: float) -> PlanarVelocity:
        """Return the (uniform) velocity at ``(x, y)``."""

        _finite(x, "x")
        _finite(y, "y")
        return PlanarVelocity(self.speed * cos(self.angle_rad), self.speed * sin(self.angle_rad))

    def potential(self, x: float, y: float) -> float:
        """Return the velocity potential in m^2/s."""

        return self.speed * (_finite(x, "x") * cos(self.angle_rad) + _finite(y, "y") * sin(self.angle_rad))

    def stream_function(self, x: float, y: float) -> float:
        """Return the stream function in m^2/s."""

        return self.speed * (_finite(y, "y") * cos(self.angle_rad) - _finite(x, "x") * sin(self.angle_rad))


@dataclass(frozen=True, slots=True)
class SourceFlow:
    """Line source (``strength > 0``) or sink (``strength < 0``) at ``(x0, y0)`` [1]_.

    Attributes
    ----------
    strength : float
        Volume flow per unit depth ``Lambda`` in m^2/s.
    x0, y0 : float, default 0
        Location of the source in m.

    Notes
    -----
    ``V_r = Lambda / (2 pi r)``, ``V_theta = 0`` (3.62);
    ``phi = Lambda ln(r) / (2 pi)`` (3.67); ``psi = Lambda theta / (2 pi)``
    (3.72), with ``theta`` in ``(-pi, pi]``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.10, Eqs. (3.62), (3.67), (3.72).
    """

    strength: float
    x0: float = 0.0
    y0: float = 0.0

    def __post_init__(self) -> None:
        _finite(self.strength, "strength")
        _finite(self.x0, "x0")
        _finite(self.y0, "y0")

    def velocity(self, x: float, y: float) -> PlanarVelocity:
        """Return the induced velocity at ``(x, y)``; singular at the source."""

        dx, dy, r2 = _offset(x, y, self.x0, self.y0)
        factor = self.strength / (2.0 * pi * r2)
        return PlanarVelocity(factor * dx, factor * dy)

    def potential(self, x: float, y: float) -> float:
        """Return ``Lambda ln(r) / (2 pi)``."""

        _, _, r2 = _offset(x, y, self.x0, self.y0)
        return self.strength * 0.5 * log(r2) / (2.0 * pi)

    def stream_function(self, x: float, y: float) -> float:
        """Return ``Lambda theta / (2 pi)``."""

        dx, dy, _ = _offset(x, y, self.x0, self.y0)
        return self.strength * atan2(dy, dx) / (2.0 * pi)


@dataclass(frozen=True, slots=True)
class DoubletFlow:
    """Doublet of strength ``kappa`` at ``(x0, y0)``, Anderson's orientation [1]_.

    Attributes
    ----------
    strength : float
        Doublet strength ``kappa`` in m^3/s.
    x0, y0 : float, default 0
        Location in m.

    Notes
    -----
    ``psi = -kappa sin(theta) / (2 pi r)`` (3.87) and
    ``phi = kappa cos(theta) / (2 pi r)`` (3.88). The velocity components
    follow by differentiating (3.87):
    ``V_r = -kappa cos(theta) / (2 pi r^2)``,
    ``V_theta = -kappa sin(theta) / (2 pi r^2)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.12, Eqs. (3.87), (3.88).
    """

    strength: float
    x0: float = 0.0
    y0: float = 0.0

    def __post_init__(self) -> None:
        _finite(self.strength, "strength")
        _finite(self.x0, "x0")
        _finite(self.y0, "y0")

    def velocity(self, x: float, y: float) -> PlanarVelocity:
        """Return the induced velocity at ``(x, y)``; singular at the doublet."""

        dx, dy, r2 = _offset(x, y, self.x0, self.y0)
        factor = self.strength / (2.0 * pi * r2 * r2)
        return PlanarVelocity(factor * (dy * dy - dx * dx), -2.0 * factor * dx * dy)

    def potential(self, x: float, y: float) -> float:
        """Return ``kappa cos(theta) / (2 pi r)``."""

        dx, _, r2 = _offset(x, y, self.x0, self.y0)
        return self.strength * dx / (2.0 * pi * r2)

    def stream_function(self, x: float, y: float) -> float:
        """Return ``-kappa sin(theta) / (2 pi r)``."""

        _, dy, r2 = _offset(x, y, self.x0, self.y0)
        return -self.strength * dy / (2.0 * pi * r2)


@dataclass(frozen=True, slots=True)
class VortexFlow:
    """Point vortex of circulation ``Gamma`` (clockwise positive) at ``(x0, y0)`` [1]_.

    Attributes
    ----------
    circulation : float
        ``Gamma`` in m^2/s; positive values turn clockwise.
    x0, y0 : float, default 0
        Location in m.

    Notes
    -----
    ``V_theta = -Gamma / (2 pi r)``, ``V_r = 0`` (3.105);
    ``phi = -Gamma theta / (2 pi)`` (3.112); ``psi = Gamma ln(r) / (2 pi)``
    (3.114).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.14, Eqs. (3.105), (3.112),
           (3.114).
    """

    circulation: float
    x0: float = 0.0
    y0: float = 0.0

    def __post_init__(self) -> None:
        _finite(self.circulation, "circulation")
        _finite(self.x0, "x0")
        _finite(self.y0, "y0")

    def velocity(self, x: float, y: float) -> PlanarVelocity:
        """Return the induced velocity at ``(x, y)``; singular at the vortex."""

        dx, dy, r2 = _offset(x, y, self.x0, self.y0)
        factor = self.circulation / (2.0 * pi * r2)
        return PlanarVelocity(factor * dy, -factor * dx)

    def potential(self, x: float, y: float) -> float:
        """Return ``-Gamma theta / (2 pi)`` with ``theta`` in ``(-pi, pi]``."""

        dx, dy, _ = _offset(x, y, self.x0, self.y0)
        return -self.circulation * atan2(dy, dx) / (2.0 * pi)

    def stream_function(self, x: float, y: float) -> float:
        """Return ``Gamma ln(r) / (2 pi)``."""

        _, _, r2 = _offset(x, y, self.x0, self.y0)
        return self.circulation * 0.5 * log(r2) / (2.0 * pi)


def superpose_velocity(flows: Iterable[ElementaryFlow], x: float, y: float) -> PlanarVelocity:
    """Sum the velocities of several elementary flows at one point.

    Parameters
    ----------
    flows : iterable of elementary flows
        Objects with a ``velocity(x, y)`` method.
    x, y : float
        Field point in m.

    Returns
    -------
    PlanarVelocity
        Superposed velocity.

    Notes
    -----
    Laplace's equation is linear, so velocities of superposed elementary
    flows add directly [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Sections 3.7 and 3.11, Eqs. (3.40),
           (3.46).
    """

    u = v = 0.0
    count = 0
    for flow in flows:
        velocity = flow.velocity(x, y)
        u += velocity.u
        v += velocity.v
        count += 1
    if count == 0:
        raise DomainError("At least one elementary flow is required.")
    return PlanarVelocity(u, v)


def superpose_stream_function(flows: Iterable[ElementaryFlow], x: float, y: float) -> float:
    """Sum the stream functions of several elementary flows at one point.

    Parameters
    ----------
    flows : iterable of elementary flows
        Objects with a ``stream_function(x, y)`` method.
    x, y : float
        Field point in m.

    Returns
    -------
    float
        Superposed stream function in m^2/s.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.11, Eq. (3.46).
    """

    values = [flow.stream_function(x, y) for flow in flows]
    if not values:
        raise DomainError("At least one elementary flow is required.")
    return sum(values)


def superpose_potential(flows: Iterable[ElementaryFlow], x: float, y: float) -> float:
    """Sum the velocity potentials of several elementary flows at one point.

    Parameters
    ----------
    flows : iterable of elementary flows
        Objects with a ``potential(x, y)`` method.
    x, y : float
        Field point in m.

    Returns
    -------
    float
        Superposed velocity potential in m^2/s.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.7, Eq. (3.40).
    """

    values = [flow.potential(x, y) for flow in flows]
    if not values:
        raise DomainError("At least one elementary flow is required.")
    return sum(values)


def pressure_coefficient(local_speed: float, freestream_speed: float) -> float:
    """Incompressible pressure coefficient ``Cp = 1 - (V / V_inf)^2``.

    Parameters
    ----------
    local_speed : float
        Local flow speed in m/s, finite and non-negative.
    freestream_speed : float
        Freestream speed in m/s, finite and positive.

    Returns
    -------
    float
        Pressure coefficient; 1 at a stagnation point.

    Notes
    -----
    Eq. (3.38), valid for steady incompressible flow without friction [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.5, Eq. (3.38).
    """

    speed = _finite(local_speed, "local_speed")
    if speed < 0.0:
        raise DomainError("Local speed must be non-negative.", field="local_speed")
    reference = _finite(freestream_speed, "freestream_speed", positive=True)
    return 1.0 - (speed / reference) ** 2


def pressure_coefficient_from_pressure(
    pressure: float, freestream_pressure: float, freestream_density: float, freestream_speed: float
) -> float:
    """Pressure coefficient from its definition ``(p - p_inf) / q_inf``.

    Parameters
    ----------
    pressure, freestream_pressure : float
        Local and freestream static pressure in Pa.
    freestream_density : float
        Freestream density in kg/m^3, positive.
    freestream_speed : float
        Freestream speed in m/s, positive.

    Returns
    -------
    float
        Pressure coefficient.

    Notes
    -----
    Eq. (3.36) with ``q_inf = rho_inf V_inf^2 / 2`` [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.5, Eq. (3.36).
    """

    density = _finite(freestream_density, "freestream_density", positive=True)
    speed = _finite(freestream_speed, "freestream_speed", positive=True)
    difference = _finite(pressure, "pressure") - _finite(freestream_pressure, "freestream_pressure")
    return difference / (0.5 * density * speed * speed)


def kutta_joukowski_lift_per_span(density: float, speed: float, circulation: float) -> float:
    """Lift per unit span ``L' = rho V Gamma`` of a closed two-dimensional body.

    Parameters
    ----------
    density : float
        Freestream density in kg/m^3, positive.
    speed : float
        Freestream speed in m/s, positive.
    circulation : float
        Circulation about the body in m^2/s (clockwise positive).

    Returns
    -------
    float
        Lift per unit span in N/m, perpendicular to the freestream.

    Notes
    -----
    Kutta-Joukowski theorem, Eq. (3.140) [1]_; also Eq. (VII) context of [2]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Sections 3.15-3.16, Eq. (3.140).
    .. [2] NACA-TR-411. T. Theodorsen, "Theory of Wing Sections of Arbitrary
           Shape", NACA Report 411, 1931, p. 231.
    """

    return (
        _finite(density, "density", positive=True)
        * _finite(speed, "speed", positive=True)
        * _finite(circulation, "circulation")
    )


@dataclass(frozen=True, slots=True)
class CircularCylinderFlow:
    """Lifting (or nonlifting) potential flow over a circular cylinder [1]_.

    Attributes
    ----------
    freestream_speed : float
        ``V_inf`` in m/s along +x, positive.
    radius : float
        Cylinder radius ``R`` in m, positive.
    circulation : float, default 0
        Circulation ``Gamma`` in m^2/s, clockwise positive.

    Notes
    -----
    Uniform flow + doublet of strength ``kappa = 2 pi V_inf R^2`` (3.98) +
    vortex. Field velocities (3.119)-(3.120); surface velocity (3.125);
    surface pressure coefficient (3.126); ``c_l = Gamma / (R V_inf)`` (3.138)
    based on the diameter ``c = 2R``; drag is zero (3.132).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Sections 3.13 and 3.15, Eqs. (3.92)-(3.138).
    """

    freestream_speed: float
    radius: float
    circulation: float = 0.0

    def __post_init__(self) -> None:
        _finite(self.freestream_speed, "freestream_speed", positive=True)
        _finite(self.radius, "radius", positive=True)
        _finite(self.circulation, "circulation")

    @property
    def doublet_strength(self) -> float:
        """Doublet strength ``kappa = 2 pi V_inf R^2`` in m^3/s (Eq. 3.98)."""

        return 2.0 * pi * self.freestream_speed * self.radius**2

    @property
    def lift_coefficient(self) -> float:
        """Lift coefficient ``Gamma / (R V_inf)`` on the diameter (Eq. 3.138)."""

        return self.circulation / (self.radius * self.freestream_speed)

    def velocity_polar(self, r: float, theta: float) -> tuple[float, float]:
        """Return ``(V_r, V_theta)`` at polar position ``(r, theta)``, ``r >= R``.

        Notes
        -----
        Eqs. (3.119) and (3.120).
        """

        radius = _finite(r, "r", positive=True)
        if radius < self.radius * (1.0 - 1.0e-12):
            raise DomainError("Field point lies inside the cylinder.", field="r")
        angle = _finite(theta, "theta")
        ratio = (self.radius / radius) ** 2
        radial = (1.0 - ratio) * self.freestream_speed * cos(angle)
        tangential = -(1.0 + ratio) * self.freestream_speed * sin(angle) - self.circulation / (
            2.0 * pi * radius
        )
        return radial, tangential

    def surface_velocity(self, theta: float) -> float:
        """Return the signed surface velocity ``V_theta`` at ``r = R`` (Eq. 3.125)."""

        return -2.0 * self.freestream_speed * sin(_finite(theta, "theta")) - self.circulation / (
            2.0 * pi * self.radius
        )

    def surface_pressure_coefficient(self, theta: float) -> float:
        """Return the surface pressure coefficient (Eq. 3.126, i.e. 3.38 with 3.125)."""

        return 1.0 - (self.surface_velocity(theta) / self.freestream_speed) ** 2

    def stagnation_points(self) -> tuple[tuple[float, float], ...]:
        """Return the stagnation points as ``(r, theta)`` pairs in the external flow.

        Returns
        -------
        tuple of (float, float)
            Two surface points when ``Gamma < 4 pi V R`` (Eq. 3.123), one at
            ``(R, -pi/2)`` when equal, and the single external point on the
            negative y axis when greater (outer root of Eq. 3.124). Angles are
            in ``[0, 2 pi)``.
        """

        g = self.circulation / (4.0 * pi * self.freestream_speed * self.radius)
        if abs(g) < 1.0:
            base = asin(-g)
            first = (pi - base) % (2.0 * pi)
            second = base % (2.0 * pi)
            return tuple(sorted(((self.radius, first), (self.radius, second)), key=lambda p: p[1]))
        quarter = 1.5 * pi if g > 0.0 else 0.5 * pi
        if abs(g) == 1.0:
            return ((self.radius, quarter),)
        half = abs(self.circulation) / (4.0 * pi * self.freestream_speed)
        outer = half + sqrt(half * half - self.radius**2)
        return ((outer, quarter),)

    def freestream_pressure_angles(self) -> tuple[float, ...]:
        """Return the surface angles where ``p = p_inf`` (``Cp = 0``), in ``[0, 2 pi)``.

        Notes
        -----
        ``Cp = 0`` in (3.126) means ``|V_theta| = V_inf``, i.e.
        ``sin(theta) = (+-1 - Gamma / (2 pi R V_inf)) / 2``; each admissible
        root gives two angles, as in Anderson's Example 3.17.
        """

        offset = self.circulation / (2.0 * pi * self.radius * self.freestream_speed)
        angles: set[float] = set()
        for sign in (1.0, -1.0):
            value = (sign - offset) / 2.0
            if -1.0 <= value <= 1.0:
                base = asin(value)
                angles.add(round(base % (2.0 * pi), 15))
                angles.add(round((pi - base) % (2.0 * pi), 15))
        return tuple(sorted(angles))


@dataclass(frozen=True, slots=True)
class JoukowskiAirfoil:
    """Exact potential flow about a Joukowski section with the Kutta condition [1]_ [2]_.

    The section is the image of a circle of radius ``R`` centred at
    ``(center_x, center_y)`` in the ``z'`` plane under the Joukowski map
    ``zeta = z' + a^2 / z'`` (TR-411 Eq. I). When ``R`` equals the distance
    from the centre to ``z' = a`` the circle passes through the critical point
    and the section has a cusped trailing edge at ``zeta = 2a``; a larger
    radius gives a smooth section with a rounded trailing edge (TR-411 curve
    ``B`` with ``psi > 0``).

    Attributes
    ----------
    scale : float
        Map constant ``a`` in m, positive.
    center_x, center_y : float
        Circle centre in m; ``center_x < 0`` adds thickness, ``center_y > 0``
        adds positive camber.
    radius : float or None
        Circle radius in m; ``None`` means the cusped (classical) section.
    angle_of_attack_rad : float
        Freestream direction measured from the x axis of the ``zeta`` plane.
    freestream_speed : float, default 1
        ``V_inf`` in m/s.

    Notes
    -----
    The circle flow is the lifting cylinder flow of Anderson with the
    freestream rotated by ``alpha``: on the circle,
    ``V_theta = -2 V sin(theta - alpha) - Gamma / (2 pi R)`` (Eq. 3.125 in the
    rotated frame). The circulation puts the rear stagnation point at the
    trailing-edge angle ``theta_TE``: ``Gamma = -4 pi R V sin(theta_TE - alpha)``,
    which for the cusped section (``theta_TE = -beta``,
    ``beta = asin(center_y / R)``) is TR-411 Eq. (VII),
    ``Gamma = 4 pi V R sin(alpha + beta)``. For a rounded trailing edge the
    trailing edge is taken as the rearmost point of the section, the
    "equivalent assumption" allowed in TR-411 footnote 4. Surface speeds on the
    section are the circle speeds divided by ``|d zeta / d z'|`` (TR-411
    p. 232), with ``d zeta / d z' = 1 - a^2 / z'^2``. Lift per span is
    ``rho V Gamma``; ``c_l`` uses the x extent of the section as chord.

    References
    ----------
    .. [1] NACA-TR-411. T. Theodorsen, "Theory of Wing Sections of Arbitrary
           Shape", NACA Report 411, 1931, pp. 230-233, Eqs. (I), (V)-(VII).
    .. [2] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.15, Eqs. (3.122)-(3.125),
           (3.140).
    """

    scale: float
    center_x: float
    center_y: float
    angle_of_attack_rad: float
    radius: float | None = None
    freestream_speed: float = 1.0

    def __post_init__(self) -> None:
        _finite(self.scale, "scale", positive=True)
        _finite(self.center_x, "center_x")
        _finite(self.center_y, "center_y")
        _finite(self.angle_of_attack_rad, "angle_of_attack_rad")
        _finite(self.freestream_speed, "freestream_speed", positive=True)
        critical = hypot(self.scale - self.center_x, self.center_y)
        if self.radius is not None:
            radius = _finite(self.radius, "radius", positive=True)
            if radius < critical * (1.0 - 1.0e-12):
                raise DomainError(
                    "Circle must enclose the critical point z' = -a and pass through or "
                    "enclose z' = a.",
                    field="radius",
                )
        if hypot(-self.scale - self.center_x, self.center_y) > self.circle_radius * (1.0 + 1.0e-12):
            raise DomainError("Circle must enclose (or pass through) the critical point z' = -a.",
                              field="center_x")

    @property
    def circle_radius(self) -> float:
        """Radius ``R`` of the generating circle in m."""

        if self.radius is None:
            return hypot(self.scale - self.center_x, self.center_y)
        return float(self.radius)

    @property
    def is_cusped(self) -> bool:
        """True when the circle passes through ``z' = a`` (cusped trailing edge)."""

        critical = hypot(self.scale - self.center_x, self.center_y)
        return self.radius is None or abs(self.circle_radius - critical) <= 1.0e-12 * critical

    def _circle_point(self, theta: float) -> complex:
        return complex(self.center_x, self.center_y) + self.circle_radius * complex(cos(theta), sin(theta))

    def _map(self, z: complex) -> complex:
        return z + self.scale**2 / z

    def surface_point(self, theta: float) -> tuple[float, float]:
        """Return the section point ``(x, y)`` that is the image of circle angle ``theta``."""

        zeta = self._map(self._circle_point(_finite(theta, "theta")))
        return zeta.real, zeta.imag

    def _x_of(self, theta: float) -> float:
        return self._map(self._circle_point(theta)).real

    def _x_slope(self, theta: float) -> float:
        """Return ``d Re(zeta) / d theta`` on the circle (chain rule, TR-411 p. 232)."""

        z = self._circle_point(theta)
        tangent = 1j * self.circle_radius * complex(cos(theta), sin(theta))
        return ((1.0 - self.scale**2 / (z * z)) * tangent).real

    def _extremum(self, sign: float, guess: float) -> float:
        """Locate the extremum of ``sign * x(theta)`` near ``guess``.

        A 720-point scan brackets the extremum; bisection on the sign of the
        analytic derivative ``d x / d theta`` then refines it to round-off.
        """

        samples = 720
        best = max(
            (guess + 2.0 * pi * (k / samples - 0.5) for k in range(samples + 1)),
            key=lambda t: sign * self._x_of(t),
        )
        lower, upper = best - 2.0 * pi / samples, best + 2.0 * pi / samples
        if sign * self._x_slope(lower) <= 0.0 or sign * self._x_slope(upper) >= 0.0:
            return best
        for _ in range(200):
            middle = 0.5 * (lower + upper)
            if sign * self._x_slope(middle) > 0.0:
                lower = middle
            else:
                upper = middle
            if upper - lower <= 4.0e-16 * max(1.0, abs(middle)):
                break
        return 0.5 * (lower + upper)

    @property
    def trailing_edge_angle(self) -> float:
        """Circle angle ``theta_TE`` of the trailing edge (rear stagnation point)."""

        if self.is_cusped:
            return atan2(-self.center_y, self.scale - self.center_x)
        return self._extremum(1.0, atan2(-self.center_y, self.scale - self.center_x))

    @property
    def leading_edge_angle(self) -> float:
        """Circle angle of the forward-most section point."""

        return self._extremum(-1.0, self.trailing_edge_angle + pi)

    @property
    def chord(self) -> float:
        """x extent of the section in m (reference chord for ``c_l``)."""

        return self._x_of(self.trailing_edge_angle) - self._x_of(self.leading_edge_angle)

    @property
    def circulation(self) -> float:
        """Kutta circulation ``Gamma`` in m^2/s (clockwise positive)."""

        return -4.0 * pi * self.circle_radius * self.freestream_speed * sin(
            self.trailing_edge_angle - self.angle_of_attack_rad
        )

    @property
    def lift_coefficient(self) -> float:
        """Exact ``c_l = 2 Gamma / (V_inf c)`` from ``L' = rho V Gamma``."""

        return 2.0 * self.circulation / (self.freestream_speed * self.chord)

    def surface_speed(self, theta: float) -> float:
        """Return the exact surface speed in m/s at the image of circle angle ``theta``.

        Raises
        ------
        DomainError
            At the cusped trailing edge itself, where the map is singular.
        """

        angle = _finite(theta, "theta")
        z = self._circle_point(angle)
        derivative = abs(1.0 - self.scale**2 / (z * z))
        if derivative <= 1.0e-12:
            raise DomainError("Surface speed is undefined at the cusped trailing edge point.")
        circle_speed = -2.0 * self.freestream_speed * sin(angle - self.angle_of_attack_rad) - (
            self.circulation / (2.0 * pi * self.circle_radius)
        )
        return abs(circle_speed) / derivative

    def surface_pressure_coefficient(self, theta: float) -> float:
        """Return the exact surface pressure coefficient (Eq. 3.38) at circle angle ``theta``."""

        return 1.0 - (self.surface_speed(theta) / self.freestream_speed) ** 2

    def surface_nodes(self, panel_count: int) -> tuple[tuple[float, float], ...]:
        """Return ``panel_count + 1`` section points, clockwise from the trailing edge.

        Points are equally spaced in circle angle, starting and ending at the
        trailing edge, going first along the lower surface (circle angle
        decreasing).
        """

        if panel_count < 8:
            raise DomainError("At least 8 panels are required.", field="panel_count")
        start = self.trailing_edge_angle
        points = [
            self.surface_point(start - 2.0 * pi * k / panel_count) for k in range(panel_count)
        ]
        points.append(points[0])
        return tuple(points)


def ellipse_surface_speed_ratio(
    thickness_parameter: float, angle_of_attack_rad: float, theta: float
) -> float:
    """Exact surface-speed ratio ``v / V`` on an elliptic section with the Kutta condition.

    Parameters
    ----------
    thickness_parameter : float
        ``psi_0 > 0``; the ellipse ``x = 2a cosh(psi_0) cos(theta)``,
        ``y = 2a sinh(psi_0) sin(theta)`` has ``t/c = tanh(psi_0)``.
    angle_of_attack_rad : float
        Incidence from the major axis.
    theta : float
        Ellipse parameter; the rear stagnation point is at ``theta = 0``
        (``x = +2a cosh psi_0``) and the flow is from ``-x`` to ``+x``.

    Returns
    -------
    float
        Surface speed ratio.

    Notes
    -----
    The ellipse is the Joukowski image (TR-411 Eqs. I-II) of the circle
    ``z' = a e^(psi0 + i theta)``. With the Kutta circulation
    ``Gamma = 4 pi V a e^psi0 sin(alpha)`` (TR-411 Eq. VII), the circle speed
    is ``2 V |sin(theta - alpha) + sin(alpha)|`` (Anderson Eq. 3.125 in a
    frame rotated by ``alpha``) and ``|d zeta / d z'| =
    2 e^(-psi0) sqrt(sinh^2 psi0 + sin^2 theta)`` (TR-411 p. 232-233), so
    ``v/V = e^psi0 |sin(theta - alpha) + sin(alpha)| /
    sqrt(sinh^2 psi0 + sin^2 theta)``, the form of TR-411 Eq. (XII) with
    ``psi' = epsilon = epsilon_T = 0``. The lift coefficient on the chord
    ``4 a cosh psi0`` is ``2 pi (1 + t/c) sin(alpha)``.

    References
    ----------
    .. [1] NACA-TR-411. T. Theodorsen, "Theory of Wing Sections of Arbitrary
           Shape", NACA Report 411, 1931, pp. 230-233, Eqs. (I), (II), (VII),
           (XII).
    .. [2] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
           6th ed., McGraw-Hill, 2017, Section 3.15, Eq. (3.125).
    """

    psi = _finite(thickness_parameter, "thickness_parameter", positive=True)
    alpha = _finite(angle_of_attack_rad, "angle_of_attack_rad")
    angle = _finite(theta, "theta")
    numerator = exp(psi) * abs(sin(angle - alpha) + sin(alpha))
    return numerator / sqrt(sinh(psi) ** 2 + sin(angle) ** 2)
