"""Staggered supersonic characteristic field between a smooth axis and one wall.

This forward initial/boundary-value process needs a supersonic inlet and an
owned axis history. It does not construct a sonic line or design an exit
contour. Interior, axis, axis-departure and wall points use the unit processes
built on [1]_, [2]_ and [3]_. Mass diagnostics integrate the isentropic mass
flux [4]_ over the represented polylines; they do not enclose discretization
error or certify shock-free flow beyond the accepted mesh.

References
----------
.. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
       compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
.. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
       averaged characteristic coefficients.
       https://ntrs.nasa.gov/citations/19980232191
.. [3] NASA-TND4960-AXIS-CONTINUITY. Anderson, 1969, Appendix C equation C1 and
       streamline definitions. https://ntrs.nasa.gov/citations/19690005234
.. [4] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""
from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass
from decimal import Context, Decimal, DecimalException, localcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import islice, pairwise

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import integer, real, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational

from .moc_axis_departure import solve_moc_axis_departure
from .moc_boundary import (
    MocAxisInputs,
    MocBoundaryResult,
    MocWallSegment,
    _wall_slope,
    solve_moc_axis_point,
    solve_moc_wall_point,
)
from .moc_unit_process import (
    _TOLERANCE,
    MocFlowPoint,
    _coefficients,
    _coordinate_cell_radius,
    _owned,
    _upper_residual,
    solve_moc_interior_point,
)

_POLICY = 'staggered_axis_wall_characteristic_field_with_polyline_mass_audit_v1'
_MASS_POLICY = 'relative_isentropic_density_decimal80_linear_flux_ring_weight_v1'


def _bounded(value, limit, label):
    """Return at most ``limit`` items of an ordered iterable as a tuple.

    Raises
    ------
    DomainError
        If ``value`` is text, bytes, a mapping or a set, cannot be iterated, or
        has more than ``limit`` items.
    """
    if isinstance(value, (str, bytes, bytearray, memoryview, Mapping, AbstractSet)):
        raise DomainError('Ordered field records required.', field=label)
    try:
        items = tuple(islice(iter(value), limit + 1))
    except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
        raise DomainError('Malformed field iterable.', field=label) from error
    if len(items) > limit:
        raise DomainError('Field record resource budget exceeded.', field=label)
    return items


@dataclass(frozen=True, slots=True)
class MocFieldInputs:
    """Validated inputs of :func:`solve_moc_field`.

    Attributes
    ----------
    axis_result : MocBoundaryResult
        Axis result for the first inlet point; it fixes gamma and the
        planar/cylindrical mode.
    inlet_points : tuple of MocFlowPoint
        Vertical inlet section of 3 to 65 points with strictly increasing
        radius, starting at the axis point; the last point lies on the wall
        and satisfies slip.
    wall_segment : MocWallSegment
        Straight wall starting at the last inlet point.
    marching_steps : int
        Number of complete fronts to march, between 1 and 128, with at most
        4096 point solves in total.
    max_iterations : int, default 128
        Corrector budget per unit process, between 1 and 256.

    Raises
    ------
    DomainError
        If a field is malformed or violates the constraints above.
    """

    axis_result: MocBoundaryResult
    inlet_points: tuple[MocFlowPoint, ...]
    wall_segment: MocWallSegment
    marching_steps: int
    max_iterations: int = 128

    def __post_init__(self):
        axis = _owned(self.axis_result, MocBoundaryResult, 'Initial axis result')
        if axis.boundary_kind != 'axis' or not isinstance(axis.inputs, MocAxisInputs):
            raise DomainError('An owned axis result is required.', field='axis_result')
        points = tuple(_owned(p, MocFlowPoint, 'Inlet point')
                       for p in _bounded(self.inlet_points, 65, 'inlet_points'))
        if len(points) < 3:
            raise DomainError('At least two inlet intervals are required.', field='inlet_points')
        wall = _owned(self.wall_segment, MocWallSegment, 'Field wall segment')
        if points[0] != axis.point or any(p.x_m != axis.point.x_m for p in points):
            raise DomainError('The vertical inlet must start at its owned axis state.', field='inlet_points')
        if any(b.radius_m <= a.radius_m for a, b in pairwise(points)):
            raise DomainError('Inlet radii must increase strictly from the axis.', field='inlet_points')
        if (wall.x0_m, wall.radius0_m) != (points[-1].x_m, points[-1].radius_m):
            raise DomainError('The wall must start at the inlet endpoint.', field='wall_segment')
        _, _, _, _, u, v = _coefficients(points[-1], axis.inputs.gamma, 1, False)
        if abs(v - _wall_slope(wall) * u) / max(abs(u), abs(v)) > _TOLERANCE:
            raise DomainError('Inlet wall state must satisfy slip.', field='inlet_points')
        steps = integer(self.marching_steps, 'marching_steps')
        budget = integer(self.max_iterations, 'max_iterations')
        if not 1 <= steps <= 128 or steps * (2 * len(points) - 1) > 4096:
            raise DomainError('Field marching resource budget exceeded.', field='marching_steps')
        if not 1 <= budget <= 256:
            raise DomainError('Iteration budget must be between 1 and 256.', field='max_iterations')
        for name, value in (('axis_result', axis), ('inlet_points', points),
                            ('wall_segment', wall), ('marching_steps', steps), ('max_iterations', budget)):
            object.__setattr__(self, name, value)


def _cross(a, b, c):
    """Return the exact z-component of ``(b - a) x (c - a)`` in the x-r plane."""
    return ((F(b.x_m) - F(a.x_m)) * (F(c.radius_m) - F(a.radius_m))
            - (F(b.radius_m) - F(a.radius_m)) * (F(c.x_m) - F(a.x_m)))


def _check_cell(points):
    """Require every vertex turn of a cell polygon to be strictly counter-clockwise.

    The orientation is evaluated exactly on the represented coordinates; no
    tolerance is applied.

    Raises
    ------
    ConvergenceError
        If the cell is folded or collapsed.
    """
    for index in range(len(points)):
        if _cross(points[index - 1], points[index], points[(index + 1) % len(points)]) <= 0:
            raise ConvergenceError('Characteristic cell folded or collapsed.', code='moc_field_folded_cell')


def _cells(old, half, new):
    """Return the cells that tile the strip between two complete fronts.

    The strip between ``old`` and ``new`` is covered by triangles on each side
    of the half front, the axis and wall end triangles, and quadrilaterals
    between them.
    """
    return (*tuple((old[j], half[j], old[j + 1]) for j in range(len(half))),
            *tuple((half[j], new[j], new[j + 1]) for j in range(len(half))),
            (old[0], new[0], half[0]),
            *tuple((old[j], half[j - 1], new[j], half[j]) for j in range(1, len(old) - 1)),
            (old[-1], half[-1], new[-1]))


def _check_front(points, inputs, *, complete):
    """Check radial ordering and that a front stays between the axis and the wall.

    On a complete front the first point may lie on the axis and the last on
    the wall (within its coordinate rounding cell); every other point must be
    strictly inside the domain.

    Raises
    ------
    ConvergenceError
        If the ordering is lost or a point crosses the axis or the wall.
    """
    if any(b.radius_m <= a.radius_m for a, b in pairwise(points)):
        raise ConvergenceError('Field front lost radial ordering.', code='moc_field_front_order')
    wall = inputs.wall_segment
    slope = _wall_slope(wall)
    for index, p in enumerate(points):
        if p.radius_m < 0 or p.radius_m == 0 and (not complete or index != 0):
            raise ConvergenceError('Field crossed the axis.', code='moc_field_outside_domain')
        distance = F(p.radius_m) - F(wall.radius0_m) - slope * (F(p.x_m) - F(wall.x0_m))
        if complete and index == len(points) - 1:
            cell = _coordinate_cell_radius(p.radius_m) + abs(slope) * _coordinate_cell_radius(p.x_m)
            outside = abs(distance) > cell
        else:
            outside = distance >= 0
        if outside:
            raise ConvergenceError('Field crossed the prescribed wall.', code='moc_field_outside_domain')


@lru_cache(maxsize=4096)
def _flux_state(point, gamma, reference_mach):
    """Return ``(rho/rho_ref) u/a*`` and ``(rho/rho_ref) v/a*`` at ``point``.

    The density ratio is isentropic [1]_,
    ``rho/rho_ref = ((2 + (gamma - 1) M_ref^2)/(2 + (gamma - 1) M^2))^(1/(gamma - 1))``;
    common factors cancel in the diagnostics. Eighty decimal digits and a
    bounded logarithmic density range are numerical policies, not interval
    enclosures. Normalization preserves uniform extreme-Mach states even when
    absolute density cannot be exported in binary64.

    Raises
    ------
    ConvergenceError
        If the density ratio exceeds the logarithmic resource bound or cannot
        be evaluated.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """
    g, m, ref = F(gamma), F(point.mach), F(reference_mach)
    ratio = (2 + (g - 1) * ref * ref) / (2 + (g - 1) * m * m)
    if ratio == 1:
        density = F(1)
    else:
        try:
            with localcontext(Context(prec=80)):
                logarithm = (Decimal(ratio.numerator) / Decimal(ratio.denominator)).ln()
                logarithm *= Decimal((g - 1).denominator) / Decimal((g - 1).numerator)
                if abs(logarithm) > 16384:
                    raise ConvergenceError('Relative density ledger resource budget exceeded.', code='moc_field_density_resource')
                density = F(logarithm.exp())
        except DecimalException as error:
            raise ConvergenceError('Relative density could not be resolved.', code='moc_field_density_unresolved') from error
    _, _, _, _, u, v = _coefficients(point, gamma, 1, False)
    return density * u, density * v


def _edge_flux(a, b, inputs):
    """Return the outward mass flux through edge ``a -> b`` of a CCW polygon.

    The mass-flux vector varies linearly along the edge. In cylindrical flow
    the integrand is weighted by the linearly varying radius, which gives the
    end weights ``(2 r_a + r_b)/6`` and ``(r_a + 2 r_b)/6`` [1]_. The common
    factor ``2`` (planar) or ``2 pi`` (cylindrical) is omitted.

    References
    ----------
    .. [1] NASA-TND4960-AXIS-CONTINUITY. Anderson, 1969, Appendix C equation C1 and
           streamline definitions. https://ntrs.nasa.gov/citations/19690005234
    """
    mode = inputs.axis_result.inputs
    ua, va = _flux_state(a, mode.gamma, inputs.axis_result.point.mach)
    ub, vb = _flux_state(b, mode.gamma, inputs.axis_result.point.mach)
    dx, dr = F(b.x_m) - F(a.x_m), F(b.radius_m) - F(a.radius_m)
    if mode.axisymmetric:
        wa, wb = 2 * F(a.radius_m) + F(b.radius_m), F(a.radius_m) + 2 * F(b.radius_m)
        return (dr * (wa * ua + wb * ub) - dx * (wa * va + wb * vb)) / 6
    return (dr * (ua + ub) - dx * (va + vb)) / 2


def _mass_metrics(inputs, fronts, half_fronts):
    """Return the mass-flux diagnostics of the marched field.

    Returns
    -------
    tuple
        ``(relative, drift, imbalance, leakage)``: mass flux through each
        complete front relative to the inlet; the largest relative drift; the
        largest cell flux imbalance relative to the sum of absolute edge
        fluxes; and the total absolute flux through wall edges relative to the
        inlet flux.

    Raises
    ------
    ConvergenceError
        If a front carries non-positive mass flux or a cell carries none.
    """
    fluxes = tuple(sum((_edge_flux(a, b, inputs) for a, b in pairwise(row)), F(0)) for row in fronts)
    if min(fluxes) <= 0:
        raise ConvergenceError('Field section has nonpositive mass flux.', code='moc_field_nonpositive_flux')
    relative = tuple(_finite_rational(value / fluxes[0], 'Relative field mass flux') for value in fluxes)
    drift = max(abs(value / fluxes[0] - 1) for value in fluxes)
    imbalance = F(0)
    for old, half, new in zip(fronts, half_fronts, fronts[1:]):
        for cell in _cells(old, half, new):
            edges = tuple(_edge_flux(cell[j - 1], cell[j], inputs) for j in range(len(cell)))
            denominator = sum(map(abs, edges))
            if denominator == 0:
                raise ConvergenceError('Cell mass flux is unresolved.', code='moc_field_zero_cell_flux')
            imbalance = max(imbalance, abs(sum(edges)) / denominator)
    leakage = sum((abs(_edge_flux(b[-1], a[-1], inputs)) for a, b in pairwise(fronts)), F(0)) / fluxes[0]
    return relative, _upper_residual(drift), _upper_residual(imbalance), _upper_residual(leakage)


@lru_cache(maxsize=4)
def _solve(inputs):
    """March the field; return the fields of :class:`MocFieldResult` after ``inputs``.

    Each step builds a half front from adjacent parents (axis departure plus
    interior points), then a complete front closed by an axis reflection, new
    interior points and a wall point. Fronts and cells are checked for
    ordering, domain and folding, and residual diagnostics are accumulated.
    """
    axis = inputs.axis_result
    mode = axis.inputs
    row = inputs.inlet_points
    fronts = [row]
    halves = []
    diagnostics = [axis.equation_residual_max, axis.rounding_adjusted_residual_max,
                   axis.coordinate_rounding_allowance_max, 0.]
    for _ in range(inputs.marching_steps):
        departure = solve_moc_axis_departure(axis_result=axis, minus_parent=row[1], max_iterations=inputs.max_iterations)
        half_results = tuple(solve_moc_interior_point(plus_parent=a, minus_parent=b, gamma=mode.gamma,
                       axisymmetric=mode.axisymmetric, max_iterations=inputs.max_iterations)
                       for a, b in pairwise(row[1:]))
        half = tuple(result.point for result in (departure, *half_results))
        _check_front(half, inputs, complete=False)
        axis = solve_moc_axis_point(incoming=half[0], axis_history=(axis.inputs.axis_history[-1], axis.point),
                    gamma=mode.gamma, axisymmetric=mode.axisymmetric, max_iterations=inputs.max_iterations)
        new_results = tuple(solve_moc_interior_point(plus_parent=a, minus_parent=b,
                    gamma=mode.gamma, axisymmetric=mode.axisymmetric, max_iterations=inputs.max_iterations)
                    for a, b in pairwise(half))
        wall_result = solve_moc_wall_point(incoming=half[-1], wall_segment=inputs.wall_segment, gamma=mode.gamma,
                    axisymmetric=mode.axisymmetric, max_iterations=inputs.max_iterations)
        new = tuple(result.point for result in (axis, *new_results, wall_result))
        _check_front(new, inputs, complete=True)
        for cell in _cells(row, half, new):
            _check_cell(cell)
        fronts.append(new)
        halves.append(half)
        for result in (departure, *half_results, axis, *new_results, wall_result):
            diagnostics[0] = max(diagnostics[0], result.equation_residual_max)
            diagnostics[1] = max(diagnostics[1], result.rounding_adjusted_residual_max)
            diagnostics[2] = max(diagnostics[2], result.coordinate_rounding_allowance_max)
        diagnostics[3] = max(diagnostics[3], departure.averaged_plus_compatibility_defect)
        row = new
    fronts, halves = tuple(fronts), tuple(halves)
    return fronts, halves, *_mass_metrics(inputs, fronts, halves), *diagnostics


@dataclass(frozen=True, slots=True)
class MocFieldResult:
    """Marched characteristic field with mass and equation diagnostics.

    Construction re-marches the field from ``inputs`` and rejects any field
    that does not match.

    Attributes
    ----------
    inputs : MocFieldInputs
        Owned inputs.
    fronts : tuple of tuple of MocFlowPoint
        ``marching_steps + 1`` complete fronts, the first being the inlet.
    half_fronts : tuple of tuple of MocFlowPoint
        ``marching_steps`` staggered half fronts.
    relative_mass_flux : tuple of float
        Mass flux through each complete front relative to the inlet.
    max_relative_mass_drift : float
        Largest ``|mdot_front / mdot_inlet - 1|``.
    max_cell_relative_mass_imbalance : float
        Largest cell net flux relative to its sum of absolute edge fluxes.
    wall_mass_leakage_relative : float
        Absolute flux through wall edges relative to the inlet flux.
    max_equation_residual : float
        Largest unit-process residual.
    max_rounding_adjusted_equation_residual : float
        Largest residual after coordinate rounding allowances.
    max_coordinate_rounding_allowance : float
        Largest coordinate rounding allowance.
    max_axis_departure_plus_defect : float
        Largest unenforced ``C+`` defect of the axis-departure points.
    numerical_policy, mass_flux_policy : str
        Identifiers of the marching and mass-flux policies.
    physical_model_qualified, sonic_initialization_complete, full_nozzle_solution : bool
        Always ``False``.
    discretization_error_bound : None
        Always ``None``.

    Raises
    ------
    DomainError
        If a field is malformed, a qualification claim is made, or the record
        does not reproduce from ``inputs``.

    Notes
    -----
    Mass fluxes use the isentropic density ratio [1]_ and, in cylindrical
    flow, the radius-weighted flux of the continuity equation [2]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    .. [2] NASA-TND4960-AXIS-CONTINUITY. Anderson, 1969, Appendix C equation C1 and
           streamline definitions. https://ntrs.nasa.gov/citations/19690005234
    """

    inputs: MocFieldInputs
    fronts: tuple[tuple[MocFlowPoint, ...], ...]
    half_fronts: tuple[tuple[MocFlowPoint, ...], ...]
    relative_mass_flux: tuple[float, ...]
    max_relative_mass_drift: float
    max_cell_relative_mass_imbalance: float
    wall_mass_leakage_relative: float
    max_equation_residual: float
    max_rounding_adjusted_equation_residual: float
    max_coordinate_rounding_allowance: float
    max_axis_departure_plus_defect: float
    numerical_policy: str = _POLICY
    mass_flux_policy: str = _MASS_POLICY
    physical_model_qualified: bool = False
    sonic_initialization_complete: bool = False
    full_nozzle_solution: bool = False
    discretization_error_bound: None = None

    def __post_init__(self):
        inputs = _owned(self.inputs, MocFieldInputs, 'Field inputs')
        object.__setattr__(self, 'inputs', inputs)
        for name, limit, width in (('fronts', inputs.marching_steps + 1, len(inputs.inlet_points)),
                                    ('half_fronts', inputs.marching_steps, len(inputs.inlet_points) - 1)):
            rows = tuple(tuple(_owned(p, MocFlowPoint, 'Field point') for p in _bounded(row, width, name))
                         for row in _bounded(getattr(self, name), limit, name))
            if len(rows) != limit or any(len(row) != width for row in rows):
                raise DomainError('Incomplete characteristic field.', field=name)
            object.__setattr__(self, name, rows)
        flux = tuple(real(v, 'Relative mass flux') for v in _bounded(self.relative_mass_flux, inputs.marching_steps + 1, 'relative_mass_flux'))
        object.__setattr__(self, 'relative_mass_flux', flux)
        for name in ('max_relative_mass_drift', 'max_cell_relative_mass_imbalance', 'wall_mass_leakage_relative',
                     'max_equation_residual', 'max_rounding_adjusted_equation_residual',
                     'max_coordinate_rounding_allowance', 'max_axis_departure_plus_defect'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        for name in ('numerical_policy', 'mass_flux_policy'):
            object.__setattr__(self, name, text(getattr(self, name), name))
        for name in ('physical_model_qualified', 'sonic_initialization_complete', 'full_nozzle_solution'):
            if type(getattr(self, name)) is not bool or getattr(self, name):
                raise DomainError('Supersonic field cannot establish qualification or sonic initialization.', field=name)
        actual = (self.fronts, self.half_fronts, flux, self.max_relative_mass_drift,
                  self.max_cell_relative_mass_imbalance, self.wall_mass_leakage_relative,
                  self.max_equation_residual, self.max_rounding_adjusted_equation_residual,
                  self.max_coordinate_rounding_allowance, self.max_axis_departure_plus_defect)
        if (self.numerical_policy != _POLICY or self.mass_flux_policy != _MASS_POLICY
                or self.discretization_error_bound is not None or actual != _solve(inputs)):
            raise DomainError('Field does not replay from its owned inputs.', field='MOC field result')


def solve_moc_field(*, axis_result, inlet_points, wall_segment, marching_steps, max_iterations=128):
    """March a complete supersonic axis-to-wall strip, retaining every front.

    Parameters
    ----------
    axis_result : MocBoundaryResult
        Axis result for the first inlet point.
    inlet_points : sequence of MocFlowPoint
        Vertical, fully supersonic inlet section from the axis to the wall.
    wall_segment : MocWallSegment or mapping
        Straight slip wall starting at the last inlet point.
    marching_steps : int
        Number of complete fronts to march.
    max_iterations : int, default 128
        Corrector budget per unit process.

    Returns
    -------
    MocFieldResult
        All fronts and half fronts with mass and equation diagnostics.

    Raises
    ------
    DomainError
        If the inputs are invalid.
    ConvergenceError
        If a unit process fails, a cell folds, a front leaves the domain, or a
        mass flux is non-positive.

    Notes
    -----
    Initial data are a vertical, fully supersonic section, its computed axis
    state and history, and a compatible finite straight slip wall. Points are
    computed with the characteristic unit processes [1]_, [2]_, axis regularity
    from continuity [3]_ and isentropic density [4]_. No wall corner, sonic
    starting condition, nozzle optimization or shock continuation is supplied.
    Cell folds and domain crossings fail; mass drift is reported rather than
    corrected to match the inlet mass by construction.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    .. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    .. [3] NASA-TND4960-AXIS-CONTINUITY. Anderson, 1969, Appendix C equation C1 and
           streamline definitions. https://ntrs.nasa.gov/citations/19690005234
    .. [4] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """
    inputs = MocFieldInputs(axis_result, inlet_points, wall_segment, marching_steps, max_iterations)
    return MocFieldResult(inputs, *_solve(inputs))
