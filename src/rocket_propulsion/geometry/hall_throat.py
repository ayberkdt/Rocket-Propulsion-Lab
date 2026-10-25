"""First/second-order Hall throat field and its bracketed sonic curve.

Constant-gamma potential-flow expansion, not a transonic CFD solution or a
complete sonic-to-MOC initializer. The reference wall is a hyperbola with
the requested throat radius and curvature. Nonlinear PDE/slip defects remain
visible. No Mach offset is introduced at the sonic curve.

Planar coefficients: McCabe, ARC R&M 3440, Appendix A.
Axisymmetric Cartesian coefficients: Kliegel/Quan, 02874-6002-R000, eqs 4-1..4-4.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import atan2, comb, degrees, inf, isfinite, nextafter

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import integer, real, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output

from .moc_unit_process import _owned

_POLICY = 'hall_cartesian_second_order_exact_rational_stretched_coordinates_v1'
_ROOT_POLICY = 'hall_hyperbolic_wall_fraction_sonic_polynomial_bernstein_monotonic_bracket_v1'


@dataclass(frozen=True, slots=True)
class HallThroatModel:
    throat_radius_m: float
    throat_curvature_radius_m: float
    gamma: float = 1.4
    axisymmetric: bool = True

    def __post_init__(self):
        for name in ('throat_radius_m', 'throat_curvature_radius_m', 'gamma'):
            value = real(getattr(self, name), name)
            if value <= (1 if name == 'gamma' else 0):
                raise DomainError('Positive radii and gamma greater than one are required.', field=name)
            object.__setattr__(self, name, value)
        if type(self.axisymmetric) is not bool:
            raise DomainError('Geometry selector must be a boolean.', field='axisymmetric')


@lru_cache(maxsize=32)
def _coefficients(model):
    """Hall polynomials in w=A*sqrt(R)*x/rt, with A**2=(sigma+1)/(gamma+1)."""
    g = F(model.gamma)
    sigma = int(model.axisymmetric)
    eps = F(model.throat_radius_m) / F(model.throat_curvature_radius_m)
    k = F(1 + sigma) / (g + 1) / eps
    c0 = -F(1 + sigma, 2 * (3 + sigma))
    a = (3 - 2 * g) / 6
    d = ((1 + sigma) * (g + 3) + 3) / (6 * (3 + sigma))
    e = F(1, 2 * (3 + sigma)) - 2 * d
    c = e - 2 * g * c0 / 3
    p = ((1 + sigma) * (F(1, 2) + d + (g - F(1, 2)) / 4) - 2 * c0) / (5 + sigma)
    q = ((1 + sigma) * (c / 2 + e + (g + F(1, 2)) * c0) + 2 * c0) / (3 + sigma)
    t = -p - q
    f = t - c0 * c - (g - F(1, 2)) * c0 * c0
    return eps, k, c0, a, c, d, e, f, p, q, t


@dataclass(frozen=True, slots=True)
class HallThroatInputs:
    model: HallThroatModel
    axial_parameter: float
    radius_fraction: float
    radial_basis: str = 'throat'
    order: int = 2

    def __post_init__(self):
        model = _owned(self.model, HallThroatModel, 'Hall model')
        object.__setattr__(self, 'model', model)
        for name in ('axial_parameter', 'radius_fraction'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        basis = text(self.radial_basis, 'radial_basis')
        if basis not in ('throat', 'wall'):
            raise DomainError('Unknown Hall radial basis.', field='radial_basis')
        object.__setattr__(self, 'radial_basis', basis)
        order = integer(self.order, 'order')
        if order not in (1, 2):
            raise DomainError('Only first and second Hall orders are implemented.', field='order')
        object.__setattr__(self, 'order', order)
        eps, k, *_ = _coefficients(model)
        wall2 = 1 + eps * F(self.axial_parameter) ** 2 / k
        r2 = F(self.radius_fraction) ** 2 * (wall2 if basis == 'wall' else 1)
        if self.radius_fraction < 0 or r2 > wall2:
            raise DomainError('Hall point lies outside the hyperbolic wall.', field='radius_fraction')


def _upper(value):
    try:
        number = float(value)
    except OverflowError as error:
        raise DomainError('Hall residual is outside finite output range.', field='Hall residual') from error
    if not isfinite(number):
        raise DomainError('Hall residual is outside finite output range.', field='Hall residual')
    return nextafter(number, inf) if F(number) < value else number


def _signed_root(square, sign, label):
    if square == 0:
        return 0.
    value = _sqrt_output(square, label)
    return value if sign > 0 else -value


def _optional(value, label, *, squared=False, sign=1):
    try:
        return (_signed_root(value, sign, label) if squared else _finite_rational(value, label)), 'finite'
    except DomainError:
        return None, 'unrepresentable'


def _ledgers(model, w, r2, order):
    """Retain velocities squared and the common radical factor symbolically."""
    eps, k, c0, a, c, d, e, f, p, q, t = _coefficients(model)
    u = 1 + eps * (w + r2 / 2 + c0)
    h = w - c0 * (r2 - 1)
    duw = eps
    dur_over_r = eps
    hr = -c0
    hw = F(1)
    if order == 2:
        u += eps * eps * (a * w * w + (r2 + c) * w + d * r2 * r2 + e * r2 + f)
        h += eps * (w * w + (4 * d * r2 + 2 * e) * w + p * r2 * r2 + q * r2 + t)
        duw += eps * eps * (2 * a * w + r2 + c)
        dur_over_r += eps * eps * (2 * w + 4 * d * r2 + 2 * e)
        hr += eps * (4 * d * w + 2 * p * r2 + q)
        hw += eps * (2 * w + 4 * d * r2 + 2 * e)
    v2 = eps * eps * r2 * h * h / k
    speed2 = u * u + v2
    sound2 = ((F(model.gamma) + 1) - (F(model.gamma) - 1) * speed2) / 2
    if u <= 0 or sound2 <= 0:
        raise DomainError('Hall expansion left the forward positive-energy branch.', field='Hall state')
    # Divide the full potential equation by sqrt(k); all terms are rational.
    terms = ((sound2 - u * u) * duw,
             -2 * u * eps * r2 * h * dur_over_r / k,
             (sound2 - v2) * eps * (h + 2 * r2 * hr) / k,
             int(model.axisymmetric) * sound2 * eps * h / k)
    irrotational = dur_over_r - eps * hw
    return u, v2, speed2, sound2, h, terms, irrotational


@lru_cache(maxsize=64)
def _evaluate(inputs):
    model = inputs.model
    eps, k, *_ = _coefficients(model)
    rt = F(model.throat_radius_m)
    w = F(inputs.axial_parameter)
    wall2 = 1 + eps * w * w / k
    r2 = F(inputs.radius_fraction) ** 2 * (wall2 if inputs.radial_basis == 'wall' else 1)
    u, v2, speed2, sound2, h, terms, irrotational = _ledgers(model, w, r2, inputs.order)
    axial = _finite_rational(u, 'Hall axial velocity / a*')
    radial, radial_status = _optional(v2, 'Hall radial velocity / a*', squared=True, sign=h)
    perturbation, perturbation_status = _optional(u - 1, 'Hall axial perturbation / a*')
    angle = None if radial is None else degrees(atan2(radial, axial))
    if angle == 0 and h != 0 and r2 != 0:
        angle = None
    mach = _sqrt_output(speed2 / sound2, 'Hall Mach')
    denominator = sum(map(abs, terms))
    residual = _upper(abs(sum(terms)) / eps)
    relative = _upper(abs(sum(terms)) / denominator) if denominator else 0.
    slip = None
    slip_status = 'not_wall'
    if r2 == wall2:
        square = eps * eps * (wall2 * h - w * u) ** 2 / (k * wall2)
        slip, slip_status = _optional(square, 'Hall wall slip / a*', squared=True)
    return {'x_m': _signed_root(rt * rt * w * w / k, w, 'Hall x'),
        'radius_m': _signed_root(rt * rt * r2, 1, 'Hall radius'),
        'wall_radius_m': _signed_root(rt * rt * wall2, 1, 'Hall wall radius'),
        'axial_velocity_over_sonic': axial, 'radial_velocity_over_sonic': radial,
        'radial_velocity_status': radial_status, 'axial_perturbation_over_sonic': perturbation,
        'axial_perturbation_status': perturbation_status, 'mach': mach, 'flow_angle_deg': angle,
        'flow_regime': 'sonic' if speed2 == 1 else 'supersonic' if speed2 > 1 else 'subsonic',
        'representation_limited': (axial == 1 and u != 1 or mach == 1 and speed2 != 1
                                or radial is None or angle is None or slip_status == 'unrepresentable'),
        'potential_residual': residual, 'relative_term_balance': relative,
        'irrotational_residual': _upper(abs(irrotational)), 'wall_slip_velocity_over_sonic': slip, 'wall_slip_status': slip_status}


@dataclass(frozen=True, slots=True)
class HallThroatResult:
    inputs: HallThroatInputs
    x_m: float
    radius_m: float
    wall_radius_m: float
    axial_velocity_over_sonic: float
    radial_velocity_over_sonic: float | None
    radial_velocity_status: str
    axial_perturbation_over_sonic: float | None
    axial_perturbation_status: str
    mach: float
    flow_angle_deg: float | None
    flow_regime: str
    representation_limited: bool
    potential_residual: float
    relative_term_balance: float
    irrotational_residual: float
    wall_slip_velocity_over_sonic: float | None
    wall_slip_status: str
    numerical_policy: str = _POLICY
    physical_model_qualified: bool = False
    sonic_initialization_complete: bool = False
    approximation_error_bound: None = None

    def __post_init__(self):
        inputs = _owned(self.inputs, HallThroatInputs, 'Hall evaluation inputs')
        object.__setattr__(self, 'inputs', inputs)
        for name in ('x_m', 'radius_m', 'wall_radius_m', 'axial_velocity_over_sonic', 'mach',
                     'potential_residual', 'relative_term_balance', 'irrotational_residual'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        for name in ('radial_velocity_over_sonic', 'axial_perturbation_over_sonic', 'flow_angle_deg', 'wall_slip_velocity_over_sonic'):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, real(getattr(self, name), name))
        for name in ('radial_velocity_status', 'axial_perturbation_status', 'flow_regime', 'numerical_policy', 'wall_slip_status'):
            object.__setattr__(self, name, text(getattr(self, name), name))
        for name in ('representation_limited', 'physical_model_qualified', 'sonic_initialization_complete'):
            if type(getattr(self, name)) is not bool:
                raise DomainError('Hall flags must be actual booleans.', field=name)
        if self.physical_model_qualified or self.sonic_initialization_complete or self.approximation_error_bound is not None or self.numerical_policy != _POLICY:
            raise DomainError('Hall approximation cannot establish physical qualification or a complete initializer.', field='Hall evidence')
        if any(getattr(self, name) != value for name, value in _evaluate(inputs).items()):
            raise DomainError('Hall result does not replay from owned inputs.', field='Hall result')


def evaluate_hall_throat(*, model, axial_parameter, radius_fraction, radial_basis='throat', order=2):
    """Evaluate the retained Hall polynomial, exposing full-PDE and wall defects."""
    inputs = HallThroatInputs(model, axial_parameter, radius_fraction, radial_basis, order)
    return HallThroatResult(inputs=inputs, **_evaluate(inputs))


def _add(*polynomials):
    result = [F(0)] * max(map(len, polynomials))
    for poly in polynomials:
        for j, value in enumerate(poly):
            result[j] += value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def _scale(poly, value):
    return tuple(value * x for x in poly)


def _multiply(a, b):
    result = [F(0)] * (len(a) + len(b) - 1)
    for j, av in enumerate(a):
        for k, bv in enumerate(b):
            result[j + k] += av * bv
    return _add(result)


def _at(poly, value):
    answer = F(0)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def _sonic_polynomial(model, eta, order):
    eps, k, c0, a, c, d, e, f, p, q, t = _coefficients(model)
    w, one = (F(0), F(1)), (F(1),)
    r2 = (eta * eta, F(0), eta * eta * eps / k)
    u1 = _add(w, _scale(r2, F(1, 2)), (c0,))
    h = _add(w, _scale(_add(r2, (-1,)), -c0))
    u = _add(one, _scale(u1, eps))
    if order == 2:
        u2 = _add(_scale(_multiply(w, w), a), _multiply(_add(r2, (c,)), w),
                  _scale(_multiply(r2, r2), d), _scale(r2, e), (f,))
        h2 = _add(_multiply(w, w), _multiply(_add(_scale(r2, 4 * d), (2 * e,)), w),
                  _scale(_multiply(r2, r2), p), _scale(r2, q), (t,))
        u = _add(u, _scale(u2, eps * eps))
        h = _add(h, _scale(h2, eps))
    return _add(_multiply(u, u), _scale(_multiply(r2, _multiply(h, h)), eps * eps / k), (-1,))


def _positive_on_interval(poly, lo, hi, depth=8):
    """Exact sufficient positivity certificate, never a sampled sign guess."""
    if not poly:
        return False
    n = len(poly) - 1
    power = tuple(sum(poly[j] * comb(j, i) * lo ** (j - i) * (hi - lo) ** i for j in range(i, n + 1))
                  for i in range(n + 1))
    bernstein = tuple(sum(power[j] * F(comb(i, j), comb(n, j)) for j in range(i + 1)) for i in range(n + 1))
    if min(bernstein) > 0:
        return True
    if max(bernstein) <= 0 or depth == 0:
        return False
    mid = (lo + hi) / 2
    return (_positive_on_interval(poly, lo, mid, depth - 1)
            and _positive_on_interval(poly, mid, hi, depth - 1))


@dataclass(frozen=True, slots=True)
class HallSonicInputs:
    model: HallThroatModel
    wall_fraction: float
    order: int = 2
    max_iterations: int = 128

    def __post_init__(self):
        object.__setattr__(self, 'model', _owned(self.model, HallThroatModel, 'Hall sonic model'))
        value = real(self.wall_fraction, 'wall_fraction')
        if not 0 <= value <= 1:
            raise DomainError('Sonic wall fraction must lie in [0,1].', field='wall_fraction')
        object.__setattr__(self, 'wall_fraction', value)
        order = integer(self.order, 'order')
        budget = integer(self.max_iterations, 'max_iterations')
        if order not in (1, 2) or not 1 <= budget <= 256:
            raise DomainError('Unsupported Hall order or root iteration budget.', field='Hall sonic inputs')
        object.__setattr__(self, 'order', order)
        object.__setattr__(self, 'max_iterations', budget)


def _float_bounds(value):
    number = _finite_rational(value, 'Hall root parameter')
    return (nextafter(number, -inf) if F(number) > value else number,
            nextafter(number, inf) if F(number) < value else number)


def _root_bounds(square, sign):
    number = _sqrt_output(square, 'Hall root coordinate') if square else 0.
    low = nextafter(number, -inf) if F(number) ** 2 > square else number
    high = nextafter(number, inf) if F(number) ** 2 < square else number
    return (low, high) if sign >= 0 else (-high, -low)


def _radial_bounds(model, eta, lo, hi):
    eps, k, *_ = _coefficients(model)
    rt2 = F(model.throat_radius_m) ** 2
    smallest = F(0) if lo <= 0 <= hi else min(lo * lo, hi * hi)
    largest = max(lo * lo, hi * hi)
    low = _root_bounds(rt2 * eta * eta * (1 + eps * smallest / k), 1)[0]
    high = _root_bounds(rt2 * eta * eta * (1 + eps * largest / k), 1)[1]
    return low, high


@lru_cache(maxsize=32)
def _locate(inputs):
    model = inputs.model
    poly = _sonic_polynomial(model, F(inputs.wall_fraction), inputs.order)
    derivative = tuple(i * poly[i] for i in range(1, len(poly)))
    lo, hi = F(-1), F(1)
    if _at(poly, lo) > 0 or _at(poly, hi) < 0:
        raise ConvergenceError('Hall sonic root is not bracketed near the throat.', code='hall_sonic_not_bracketed')
    if not _positive_on_interval(derivative, lo, hi):
        raise ConvergenceError('Hall sonic branch has no monotonicity certificate.', code='hall_sonic_branch_unresolved')
    for count in range(1, inputs.max_iterations + 1):
        mid = (lo + hi) / 2
        value = _at(poly, mid)
        if value == 0:
            lo = hi = mid
        elif value < 0:
            lo = mid
        else:
            hi = mid
        if float(lo) == float(hi):
            break
    else:
        raise ConvergenceError('Hall sonic root iteration budget exhausted.', code='hall_sonic_not_converged')
    lower, _ = _float_bounds(lo)
    _, upper = _float_bounds(hi)
    lower_value, upper_value = _at(poly, F(lower)), _at(poly, F(upper))
    if lower_value > 0 or upper_value < 0:
        raise ConvergenceError('Exported Hall parameter endpoints do not bracket the sonic root.', code='hall_sonic_output_bracket')
    candidate = float((lo + hi) / 2)
    point = evaluate_hall_throat(model=model, axial_parameter=candidate,
              radius_fraction=inputs.wall_fraction, radial_basis='wall', order=inputs.order)
    _eps, k, *_ = _coefficients(model)
    rt = F(model.throat_radius_m)
    xl = _root_bounds(rt * rt * F(lower) ** 2 / k, lower)[0]
    xu = _root_bounds(rt * rt * F(upper) ** 2 / k, upper)[1]
    rl, ru = _radial_bounds(model, F(inputs.wall_fraction), F(lower), F(upper))
    return {'point': point, 'iterations': count, 'parameter_lower': lower, 'parameter_upper': upper,
                'x_lower_m': xl, 'x_upper_m': xu, 'radius_lower_m': rl, 'radius_upper_m': ru,
                'sonic_function_lower': -_upper(-lower_value),
                'sonic_function_upper': _upper(upper_value)}


@dataclass(frozen=True, slots=True)
class HallSonicResult:
    inputs: HallSonicInputs
    point: HallThroatResult
    iterations: int
    parameter_lower: float
    parameter_upper: float
    x_lower_m: float
    x_upper_m: float
    radius_lower_m: float
    radius_upper_m: float
    sonic_function_lower: float
    sonic_function_upper: float
    numerical_policy: str = _ROOT_POLICY
    physical_model_qualified: bool = False
    sonic_initialization_complete: bool = False
    approximation_error_bound: None = None

    def __post_init__(self):
        inputs = _owned(self.inputs, HallSonicInputs, 'Hall sonic inputs')
        object.__setattr__(self, 'inputs', inputs)
        object.__setattr__(self, 'point', _owned(self.point, HallThroatResult, 'Hall sonic point'))
        object.__setattr__(self, 'iterations', integer(self.iterations, 'iterations'))
        for name in ('parameter_lower', 'parameter_upper', 'x_lower_m', 'x_upper_m',
                     'radius_lower_m', 'radius_upper_m', 'sonic_function_lower', 'sonic_function_upper'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        object.__setattr__(self, 'numerical_policy', text(self.numerical_policy, 'numerical_policy'))
        for name in ('physical_model_qualified', 'sonic_initialization_complete'):
            if type(getattr(self, name)) is not bool or getattr(self, name):
                raise DomainError('A polynomial sonic root cannot establish physical or initializer qualification.', field=name)
        if self.approximation_error_bound is not None or self.numerical_policy != _ROOT_POLICY or any(
                getattr(self, name) != value for name, value in _locate(inputs).items()):
            raise DomainError('Hall sonic record does not replay.', field='Hall sonic result')


def locate_hall_sonic_point(*, model, wall_fraction, order=2, max_iterations=128):
    """Bracket |V|/a*=1 inside the retained Hall model and hyperbolic wall.

    wall_fraction=0 is the axis and 1 is the actual hyperbolic wall. The
    polynomial root certificate bounds numerical root error, not Hall's
    approximation error or the true transonic flow. No supersonic offset or
    downstream axis-history initialization is supplied.
    """
    inputs = HallSonicInputs(model, wall_fraction, order, max_iterations)
    return HallSonicResult(inputs=inputs, **_locate(inputs))
