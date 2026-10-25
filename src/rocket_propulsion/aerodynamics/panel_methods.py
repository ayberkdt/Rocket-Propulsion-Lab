"""Constant-strength panel methods for two-dimensional incompressible flow.

Two methods are provided on a closed polygon of straight panels whose nodes
run **clockwise** (for an airfoil: trailing edge, lower surface, leading edge,
upper surface, trailing edge), so that the outward normal of panel ``j`` is at
``beta_j = Phi_j + pi/2`` (Anderson Fig. 3.42, Eqs. 3.160):

* :func:`solve_source_panels` - the nonlifting source panel method of
  Anderson Section 3.17: Eq. (3.153) with the geometric integrals ``I_ij``
  (3.163) and the tangential integrals ``J_ij`` (3.165); surface velocity
  (3.156), pressure coefficient (3.38) and the closure check (3.157).
* :func:`solve_lifting_panels` - a lifting method that combines those source
  panels with a vortex sheet of the same (unknown) strength ``gamma`` on every
  panel, closed by the Kutta condition in the form of Anderson Eq. (4.81),
  ``gamma_1 = -gamma_N`` applied to the surface velocities of the two
  trailing-edge panels (``V_1 + V_N = 0`` along the clockwise tangent, i.e.
  equal speeds leaving the trailing edge on both surfaces). Anderson
  (Section 4.10, p. 373) notes that such source-vortex combinations are the
  usual practical form of the vortex panel method. The vortex-panel influence
  integrals follow from Eqs. (4.72)-(4.77): with the substitution
  (3.160)-(3.161), the integrand of ``d theta_ij / d n_i`` is identical to that
  of ``d ln r_ij / d s_i`` in (3.165), and the integrand of
  ``d theta_ij / d s_i`` is the negative of that of ``d ln r_ij / d n_i``
  (3.159). Hence ``int_j d theta_ij / d n_i ds_j = J_ij`` and
  ``int_j d theta_ij / d s_i ds_j = -I_ij``; with the potential (4.72) a
  panel of strength ``gamma`` induces ``-gamma J_ij / (2 pi)`` normal and
  ``+gamma I_ij / (2 pi)`` tangential velocity, plus the self term
  ``gamma / 2`` on the outer side (the velocity jump across the sheet is
  ``gamma``, Anderson p. 372-373). Lift follows from Kutta-Joukowski,
  ``L' = rho V Gamma`` with ``Gamma = gamma sum S_j`` (4.82)-(4.83), and also
  from integrating the panel pressures with Eqs. (1.15)-(1.19).

Anderson's own first-order formulation with one unknown ``gamma_j`` per panel
and a dropped control point (Eqs. 4.80-4.81) was evaluated during development
and is not provided: on sections with a thin trailing edge it gave lift
coefficients that did not converge with panel refinement.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
       6th ed., McGraw-Hill, 2017, Section 1.5, Eqs. (1.15)-(1.19); Section
       3.17, Eqs. (3.141)-(3.165); Section 4.10, Eqs. (4.72)-(4.83).
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass
from math import atan, atan2, cos, degrees, hypot, isfinite, log, pi, sin, sqrt

from rocket_propulsion.core.errors import DomainError

from ._subsonic_numerics import gauss_legendre, solve_dense

# Library threshold (not a source value): below it the contour is treated as cusped.
_MINIMUM_TRAILING_EDGE_ANGLE_DEG = 2.5


@dataclass(frozen=True, slots=True)
class PanelGeometry:
    """Straight-panel discretisation of a closed clockwise contour [1]_.

    Attributes
    ----------
    nodes : tuple of (float, float)
        Boundary points, first equal to last.
    lengths : tuple of float
        Panel lengths ``S_j``.
    angles : tuple of float
        Panel angles ``Phi_j`` from the x axis to the panel, counterclockwise.
    control_points : tuple of (float, float)
        Panel midpoints ``(x_i, y_i)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 3.17, Figs.
           3.40 and 3.42, Eqs. (3.160)-(3.161).
    """

    nodes: tuple[tuple[float, float], ...]
    lengths: tuple[float, ...]
    angles: tuple[float, ...]
    control_points: tuple[tuple[float, float], ...]

    @property
    def panel_count(self) -> int:
        """Number of panels."""

        return len(self.lengths)

    @property
    def normal_angles(self) -> tuple[float, ...]:
        """Outward-normal angles ``beta_j = Phi_j + pi/2`` (Anderson Eq. 3.160)."""

        return tuple(angle + 0.5 * pi for angle in self.angles)


def panel_geometry(nodes: Sequence[Sequence[float]]) -> PanelGeometry:
    """Build and validate a closed clockwise panel contour.

    Parameters
    ----------
    nodes : sequence of (x, y)
        At least 4 points in m (or chord fractions); the first and last must
        coincide; the contour must be traversed clockwise.

    Returns
    -------
    PanelGeometry
        Panel lengths, angles and control points.

    Raises
    ------
    DomainError
        If a coordinate is not finite, the contour is open, a panel has zero
        length, or the contour runs counterclockwise (positive signed area).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 3.17.
    """

    points = tuple((float(point[0]), float(point[1])) for point in nodes)
    if len(points) < 4:
        raise DomainError("At least three panels (four nodes) are required.", field="nodes")
    if not all(isfinite(x) and isfinite(y) for x, y in points):
        raise DomainError("Panel nodes must be finite.", field="nodes")
    if points[0] != points[-1]:
        raise DomainError("The panel contour must be closed (first node equal to last).", field="nodes")
    lengths, angles, controls = [], [], []
    signed_area = 0.0
    for (x0, y0), (x1, y1) in itertools.pairwise(points):
        length = hypot(x1 - x0, y1 - y0)
        if length == 0.0:
            raise DomainError("Panels must have nonzero length.", field="nodes")
        lengths.append(length)
        angles.append(atan2(y1 - y0, x1 - x0))
        controls.append((0.5 * (x0 + x1), 0.5 * (y0 + y1)))
        signed_area += 0.5 * (x0 * y1 - x1 * y0)
    if signed_area >= 0.0:
        raise DomainError("Panel nodes must run clockwise around the body.", field="nodes")
    return PanelGeometry(points, tuple(lengths), tuple(angles), tuple(controls))


def _panel_integrals(geometry: PanelGeometry, i: int, j: int) -> tuple[float, float]:
    """Return ``(I_ij, J_ij)`` of Anderson Eqs. (3.163) and (3.165) for ``i != j``.

    Notes
    -----
    When the control point lies on the extension of panel ``j``
    (``E = sqrt(B - A^2)`` vanishes) the closed forms are 0/0; the integrals
    (3.162) and (3.165) are then evaluated by 32-point Gauss-Legendre
    quadrature of their (smooth) integrands.
    """

    xi, yi = geometry.control_points[i]
    xj, yj = geometry.nodes[j]
    phi_i, phi_j = geometry.angles[i], geometry.angles[j]
    length = geometry.lengths[j]
    a = -(xi - xj) * cos(phi_j) - (yi - yj) * sin(phi_j)
    b = (xi - xj) ** 2 + (yi - yj) ** 2
    c = sin(phi_i - phi_j)
    d = (yi - yj) * cos(phi_i) - (xi - xj) * sin(phi_i)
    e_squared = b - a * a
    if e_squared > 1.0e-20 * b:
        e = sqrt(e_squared)
        log_term = log((length * length + 2.0 * a * length + b) / b)
        atan_term = atan((length + a) / e) - atan(a / e)
        normal = 0.5 * c * log_term + (d - a * c) / e * atan_term
        tangential = (d - a * c) / (2.0 * e) * log_term - c * atan_term
        return normal, tangential
    p = (xi - xj) * cos(phi_i) + (yi - yj) * sin(phi_i)
    q = cos(phi_i - phi_j)
    nodes, weights = gauss_legendre(32)
    normal = tangential = 0.0
    for node, weight in zip(nodes, weights, strict=True):
        s = 0.5 * length * (node + 1.0)
        denominator = s * s + 2.0 * a * s + b
        normal += weight * (c * s + d) / denominator
        tangential += weight * (p - q * s) / denominator
    return 0.5 * length * normal, 0.5 * length * tangential


def influence_matrices(geometry: PanelGeometry) -> tuple[tuple[tuple[float, ...], ...], tuple[tuple[float, ...], ...]]:
    """Return the source-panel geometric integrals ``I_ij`` and ``J_ij``.

    Parameters
    ----------
    geometry : PanelGeometry
        Closed clockwise contour.

    Returns
    -------
    normal, tangential : tuple of tuple of float
        ``I_ij = int_j d(ln r_ij)/d n_i ds_j`` (Eq. 3.163) and
        ``J_ij = int_j d(ln r_ij)/d s ds_j`` (Eq. 3.165); diagonal entries are
        zero (the self terms are added by the solvers).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 3.17, Eqs.
           (3.158)-(3.165).
    """

    count = geometry.panel_count
    normal = [[0.0] * count for _ in range(count)]
    tangential = [[0.0] * count for _ in range(count)]
    for i in range(count):
        for j in range(count):
            if i != j:
                normal[i][j], tangential[i][j] = _panel_integrals(geometry, i, j)
    return tuple(tuple(row) for row in normal), tuple(tuple(row) for row in tangential)


def _freestream(speed: float, angle: float) -> tuple[float, float]:
    value = float(speed)
    if not isfinite(value) or value <= 0.0:
        raise DomainError("Freestream speed must be finite and positive.", field="freestream_speed")
    alpha = float(angle)
    if not isfinite(alpha):
        raise DomainError("Angle of attack must be finite.", field="angle_of_attack_rad")
    return value, alpha


@dataclass(frozen=True, slots=True)
class SourcePanelSolution:
    """Nonlifting source-panel solution [1]_.

    Attributes
    ----------
    geometry : PanelGeometry
        Panel discretisation.
    source_strengths : tuple of float
        ``lambda_j`` per unit length in m/s.
    surface_velocities : tuple of float
        Tangential velocity ``V_i`` at each control point, positive along the
        clockwise panel direction (Eq. 3.156).
    pressure_coefficients : tuple of float
        ``Cp_i = 1 - (V_i / V_inf)^2`` (Eq. 3.38).
    source_closure : float
        ``sum lambda_j S_j / sum |lambda_j| S_j``; zero for an exact closed
        body (Eq. 3.157).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 3.17, Eqs.
           (3.153)-(3.157).
    """

    geometry: PanelGeometry
    source_strengths: tuple[float, ...]
    surface_velocities: tuple[float, ...]
    pressure_coefficients: tuple[float, ...]
    source_closure: float


def solve_source_panels(
    nodes: Sequence[Sequence[float]],
    *,
    freestream_speed: float = 1.0,
    angle_of_attack_rad: float = 0.0,
) -> SourcePanelSolution:
    """Solve the nonlifting flow over a closed body with constant source panels.

    Parameters
    ----------
    nodes : sequence of (x, y)
        Closed clockwise contour (see :func:`panel_geometry`).
    freestream_speed : float, default 1
        ``V_inf`` in m/s.
    angle_of_attack_rad : float, default 0
        Freestream direction from the x axis.

    Returns
    -------
    SourcePanelSolution
        Strengths, surface velocities and pressure coefficients.

    Notes
    -----
    Solves ``pi lambda_i + sum_j lambda_j I_ij = -2 pi V cos(beta_i - alpha)``,
    i.e. Eq. (3.153) multiplied by ``2 pi`` (cf. Eq. 3.164), then evaluates
    ``V_i = V sin(beta_i - alpha) + sum_j lambda_j J_ij / (2 pi)`` (3.156).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 3.17, Eqs.
           (3.148)-(3.157), Example 3.19.
    """

    speed, alpha = _freestream(freestream_speed, angle_of_attack_rad)
    geometry = panel_geometry(nodes)
    normal, tangential = influence_matrices(geometry)
    count = geometry.panel_count
    betas = geometry.normal_angles
    matrix = [[pi if i == j else normal[i][j] for j in range(count)] for i in range(count)]
    rhs = [-2.0 * pi * speed * cos(betas[i] - alpha) for i in range(count)]
    strengths = solve_dense(matrix, rhs)
    velocities = tuple(
        speed * sin(betas[i] - alpha)
        + sum(strengths[j] * tangential[i][j] for j in range(count)) / (2.0 * pi)
        for i in range(count)
    )
    weighted = [strength * length for strength, length in zip(strengths, geometry.lengths, strict=True)]
    magnitude = sum(abs(value) for value in weighted)
    closure = sum(weighted) / magnitude if magnitude > 0.0 else 0.0
    return SourcePanelSolution(
        geometry=geometry,
        source_strengths=tuple(strengths),
        surface_velocities=velocities,
        pressure_coefficients=tuple(1.0 - (v / speed) ** 2 for v in velocities),
        source_closure=closure,
    )


@dataclass(frozen=True, slots=True)
class LiftingPanelSolution:
    """Lifting source-plus-vortex panel solution with the Kutta condition [1]_.

    Attributes
    ----------
    geometry : PanelGeometry
        Panel discretisation.
    source_strengths : tuple of float
        ``lambda_j`` in m/s.
    vortex_strength : float
        Uniform vortex-sheet strength ``gamma`` in m/s (clockwise positive).
    circulation : float
        ``Gamma = gamma sum S_j`` in m^2/s (Eq. 4.82 with uniform ``gamma``).
    surface_velocities : tuple of float
        Tangential velocity at each control point along the clockwise panel
        direction.
    pressure_coefficients : tuple of float
        ``1 - (V_i / V_inf)^2``.
    chord : float
        Reference chord used for the coefficients.
    lift_coefficient : float
        ``2 Gamma / (V_inf c)`` from ``L' = rho V Gamma`` (Eq. 4.83).
    pressure_lift_coefficient, pressure_drag_coefficient : float
        From integrating ``Cp`` over the panels (Eqs. 1.15-1.19); the drag is
        zero in exact potential flow and measures discretisation error.
    moment_coefficient : float
        Pitching moment about ``moment_reference`` on ``c^2``, nose-up
        positive (pressure terms of Eq. 1.17).
    moment_reference : (float, float)
        Moment reference point.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Sections 1.5, 3.17 and
           4.10.
    """

    geometry: PanelGeometry
    source_strengths: tuple[float, ...]
    vortex_strength: float
    circulation: float
    surface_velocities: tuple[float, ...]
    pressure_coefficients: tuple[float, ...]
    chord: float
    lift_coefficient: float
    pressure_lift_coefficient: float
    pressure_drag_coefficient: float
    moment_coefficient: float
    moment_reference: tuple[float, float]


def trailing_edge_angle_deg(geometry: PanelGeometry) -> float:
    """Included angle in degrees between the two panels meeting at the first node.

    Parameters
    ----------
    geometry : PanelGeometry
        Contour whose first node is the trailing edge.

    Returns
    -------
    float
        Angle between the first panel and the reversed last panel.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 4.5 (finite-angle
           and cusped trailing edges).
    """

    first = geometry.angles[0]
    last_reversed = geometry.angles[-1] + pi
    difference = (first - last_reversed + pi) % (2.0 * pi) - pi
    return abs(degrees(difference))


def solve_lifting_panels(
    nodes: Sequence[Sequence[float]],
    angle_of_attack_rad: float,
    *,
    freestream_speed: float = 1.0,
    chord: float | None = None,
    moment_reference: tuple[float, float] | None = None,
) -> LiftingPanelSolution:
    """Solve the lifting flow over a section with source panels, a vortex sheet and Kutta.

    Parameters
    ----------
    nodes : sequence of (x, y)
        Closed clockwise contour starting and ending at the trailing edge
        (as returned by ``NacaSection.panel_nodes`` or
        ``JoukowskiAirfoil.surface_nodes``).
    angle_of_attack_rad : float
        Freestream direction from the x axis.
    freestream_speed : float, default 1
        ``V_inf`` in m/s.
    chord : float, optional
        Reference chord; default is the x extent of the nodes.
    moment_reference : (float, float), optional
        Moment reference point; default is a quarter chord behind the
        forward-most node, at its height.

    Returns
    -------
    LiftingPanelSolution
        Strengths, surface velocities and force coefficients.

    Raises
    ------
    DomainError
        For an invalid contour (see :func:`panel_geometry`), a trailing-edge
        included angle below 2.5 degrees (treated as a cusp, for which the
        method was found not to converge), or invalid freestream inputs.

    Notes
    -----
    Unknowns ``lambda_1 ... lambda_N`` and ``gamma``. For each control point
    (Eq. 3.153 and 4.79 multiplied by ``2 pi``)::

        pi lambda_i + sum_j lambda_j I_ij - gamma sum_j J_ij = -2 pi V cos(beta_i - alpha)

    and the Kutta condition (4.81) on the trailing-edge panels ``1`` and ``N``
    with tangential velocities::

        V_i = V sin(beta_i - alpha) + (sum_j lambda_j J_ij + gamma (pi + sum_j I_ij)) / (2 pi)

    set to ``V_1 + V_N = 0``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 1.5, Eqs.
           (1.15)-(1.19); Section 3.17, Eqs. (3.148)-(3.165); Section 4.10,
           Eqs. (4.72)-(4.83).
    """

    speed, alpha = _freestream(freestream_speed, angle_of_attack_rad)
    geometry = panel_geometry(nodes)
    if trailing_edge_angle_deg(geometry) < _MINIMUM_TRAILING_EDGE_ANGLE_DEG:
        raise DomainError(
            "Trailing-edge included angle is below 2.5 degrees (cusp); the lifting panel method "
            "does not converge for cusped trailing edges.",
            field="nodes",
        )
    normal, tangential = influence_matrices(geometry)
    count = geometry.panel_count
    betas = geometry.normal_angles
    vortex_normal = [sum(row) for row in tangential]
    vortex_tangential = [pi + sum(row) for row in normal]
    matrix = [
        [pi if i == j else normal[i][j] for j in range(count)] + [-vortex_normal[i]]
        for i in range(count)
    ]
    rhs = [-2.0 * pi * speed * cos(betas[i] - alpha) for i in range(count)]
    last = count - 1
    matrix.append(
        [tangential[0][j] + tangential[last][j] for j in range(count)]
        + [vortex_tangential[0] + vortex_tangential[last]]
    )
    rhs.append(-2.0 * pi * speed * (sin(betas[0] - alpha) + sin(betas[last] - alpha)))
    solution = solve_dense(matrix, rhs)
    strengths, gamma = solution[:count], solution[count]
    velocities = tuple(
        speed * sin(betas[i] - alpha)
        + (sum(strengths[j] * tangential[i][j] for j in range(count)) + gamma * vortex_tangential[i])
        / (2.0 * pi)
        for i in range(count)
    )
    pressures = tuple(1.0 - (v / speed) ** 2 for v in velocities)

    xs = [x for x, _ in geometry.nodes]
    reference_chord = (max(xs) - min(xs)) if chord is None else float(chord)
    if not isfinite(reference_chord) or reference_chord <= 0.0:
        raise DomainError("Reference chord must be finite and positive.", field="chord")
    if moment_reference is None:
        leading = min(geometry.nodes, key=lambda point: point[0])
        reference = (leading[0] + 0.25 * reference_chord, leading[1])
    else:
        reference = (float(moment_reference[0]), float(moment_reference[1]))
    force_x = force_y = moment = 0.0
    for cp, length, beta, (xc, yc) in zip(pressures, geometry.lengths, betas, geometry.control_points,
                                          strict=True):
        nx, ny = cos(beta), sin(beta)
        force_x -= cp * nx * length
        force_y -= cp * ny * length
        moment += cp * length * ((xc - reference[0]) * ny - (yc - reference[1]) * nx)
    force_x /= reference_chord
    force_y /= reference_chord
    circulation = gamma * sum(geometry.lengths)
    return LiftingPanelSolution(
        geometry=geometry,
        source_strengths=tuple(strengths),
        vortex_strength=gamma,
        circulation=circulation,
        surface_velocities=velocities,
        pressure_coefficients=pressures,
        chord=reference_chord,
        lift_coefficient=2.0 * circulation / (speed * reference_chord),
        pressure_lift_coefficient=force_y * cos(alpha) - force_x * sin(alpha),
        pressure_drag_coefficient=force_y * sin(alpha) + force_x * cos(alpha),
        moment_coefficient=moment / reference_chord**2,
        moment_reference=reference,
    )
