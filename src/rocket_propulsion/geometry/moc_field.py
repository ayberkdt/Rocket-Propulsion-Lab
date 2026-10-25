"""Staggered supersonic characteristic field between a smooth axis and one wall.

This forward initial/boundary-value process needs a supersonic inlet and an
owned axis history. It does not invent a sonic line or design an exit contour.
Mass diagnostics integrate the represented polylines; they do not enclose
discretization error or certify shock-free flow beyond the accepted mesh.
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
    return ((F(b.x_m) - F(a.x_m)) * (F(c.radius_m) - F(a.radius_m))
            - (F(b.radius_m) - F(a.radius_m)) * (F(c.x_m) - F(a.x_m)))


def _check_cell(points):
    """Exact represented orientation; no epsilon may excuse a folded cell."""
    for index in range(len(points)):
        if _cross(points[index - 1], points[index], points[(index + 1) % len(points)]) <= 0:
            raise ConvergenceError('Characteristic cell folded or collapsed.', code='moc_field_folded_cell')


def _cells(old, half, new):
    """Tile the complete strip, including both Cauchy-front triangle fans."""
    return (*tuple((old[j], half[j], old[j + 1]) for j in range(len(half))),
            *tuple((half[j], new[j], new[j + 1]) for j in range(len(half))),
            (old[0], new[0], half[0]),
            *tuple((old[j], half[j - 1], new[j], half[j]) for j in range(1, len(old) - 1)),
            (old[-1], half[-1], new[-1]))


def _check_front(points, inputs, *, complete):
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
    """rho/rho_ref times u/a*,v/a*; common factors cancel in the diagnostics.

    Eighty Decimal digits and a bounded logarithmic density range are numerical
    policies, not interval enclosures. Normalization preserves uniform extreme
    Mach states even when absolute density cannot be exported in binary64.
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
    """Outward flux along a CCW polygon edge; omit common 2 or 2*pi."""
    mode = inputs.axis_result.inputs
    ua, va = _flux_state(a, mode.gamma, inputs.axis_result.point.mach)
    ub, vb = _flux_state(b, mode.gamma, inputs.axis_result.point.mach)
    dx, dr = F(b.x_m) - F(a.x_m), F(b.radius_m) - F(a.radius_m)
    if mode.axisymmetric:
        wa, wb = 2 * F(a.radius_m) + F(b.radius_m), F(a.radius_m) + 2 * F(b.radius_m)
        return (dr * (wa * ua + wb * ub) - dx * (wa * va + wb * vb)) / 6
    return (dr * (ua + ub) - dx * (va + vb)) / 2


def _mass_metrics(inputs, fronts, half_fronts):
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
    """Adjacent-parent half front, then reflected/closed complete front."""
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

    Initial data are a vertical, fully supersonic section, its computed axis
    state/history and a compatible finite straight slip wall. No wall corner,
    sonic starting condition, nozzle optimization or shock continuation is
    supplied. Cell folds and domain crossings fail; mass drift remains visible
    instead of being corrected to match inlet mass by construction.
    """
    inputs = MocFieldInputs(axis_result, inlet_points, wall_segment, marching_steps, max_iterations)
    return MocFieldResult(inputs, *_solve(inputs))
