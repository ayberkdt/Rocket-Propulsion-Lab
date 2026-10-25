"""First/second-order Hall throat field and its bracketed sonic curve.

Hall's method expands the steady, irrotational, constant-gamma transonic
throat flow in powers of ``eps = r_t / R_c``, the ratio of throat radius to
wall radius of curvature at the throat. It is a potential-flow expansion, not
a transonic CFD solution or a complete sonic-to-MOC initializer. The
reference wall is the hyperbola ``r_w^2 = r_t^2 + (r_t / R_c) x^2`` with the
requested throat radius and curvature. Residuals of the full nonlinear
potential equation and of wall slip are reported rather than suppressed. No
Mach offset is introduced at the sonic curve.

Planar coefficients follow McCabe, ARC R&M 3440, Appendix A [1]_.
Axisymmetric Cartesian coefficients follow 02874-6002-R000, equations 4-1 to
4-4 [2]_. Sauer's near-critical-flow theory [3]_ is cited as historical
context.

References
----------
.. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
       context, ARC R&M 3440 (report 1964, issued 1967).
       https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
.. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
       and comparison of transonic approximations, 02874-6002-R000.
       https://ntrs.nasa.gov/citations/19670009614
.. [3] SAUER-NACA-TM1147. Near-critical-flow and curved-sonic-line theory,
       NACA TM 1147, 1947.
       https://ntrs.nasa.gov/citations/20050019415
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
    """Throat geometry and gas for Hall's transonic expansion [1]_, [2]_.

    Attributes
    ----------
    throat_radius_m : float
        Throat radius (axisymmetric) or half-height (planar) ``r_t`` in m,
        positive.
    throat_curvature_radius_m : float
        Wall radius of curvature ``R_c`` at the throat in m, positive.
    gamma : float, default 1.4
        Ratio of specific heats, greater than one.
    axisymmetric : bool, default True
        ``True`` for axisymmetric flow (``sigma = 1``), ``False`` for planar
        flow (``sigma = 0``).

    Raises
    ------
    DomainError
        If a radius is not positive, ``gamma <= 1`` or ``axisymmetric`` is not
        a ``bool``.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """

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
    """Return the exact rational Hall coefficients for ``model``.

    The expansion uses the stretched axial coordinate
    ``w = A sqrt(R) x / r_t`` with ``A^2 = (sigma + 1)/(gamma + 1)`` and
    ``R = R_c / r_t``.

    Returns
    -------
    tuple of Fraction
        ``(eps, k, c0, a, c, d, e, f, p, q, t)`` with ``eps = r_t / R_c`` and
        ``k = A^2 R``, so that ``x = r_t w / sqrt(k)``; the remaining entries
        are the first- and second-order polynomial coefficients of the axial
        velocity ``u`` and of the radial-velocity factor ``h`` (``v`` is
        proportional to ``r h``) [1]_, [2]_.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
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
    """Validated evaluation point for :func:`evaluate_hall_throat`.

    Attributes
    ----------
    model : HallThroatModel
        Owned copy of the throat model.
    axial_parameter : float
        Stretched axial coordinate ``w``; ``x = r_t w / sqrt(k)``.
    radius_fraction : float
        Non-negative radial coordinate as a fraction of ``r_t``
        (``radial_basis="throat"``) or of the local wall radius
        (``radial_basis="wall"``).
    radial_basis : str, default "throat"
        ``"throat"`` or ``"wall"``.
    order : int, default 2
        Expansion order, 1 or 2.

    Raises
    ------
    DomainError
        If a field is malformed or the point lies outside the hyperbolic wall.
    """

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
    """Return the smallest binary64 value not below the exact rational ``value``.

    Raises
    ------
    DomainError
        If ``value`` overflows binary64.
    """
    try:
        number = float(value)
    except OverflowError as error:
        raise DomainError('Hall residual is outside finite output range.', field='Hall residual') from error
    if not isfinite(number):
        raise DomainError('Hall residual is outside finite output range.', field='Hall residual')
    return nextafter(number, inf) if F(number) < value else number


def _signed_root(square, sign, label):
    """Return ``sqrt(square)`` rounded to nearest-even, with the sign of ``sign``."""
    if square == 0:
        return 0.
    value = _sqrt_output(square, label)
    return value if sign > 0 else -value


def _optional(value, label, *, squared=False, sign=1):
    """Return ``(float, "finite")`` or ``(None, "unrepresentable")`` for an exact value.

    With ``squared=True`` the signed square root of ``value`` is exported.
    """
    try:
        return (_signed_root(value, sign, label) if squared else _finite_rational(value, label)), 'finite'
    except DomainError:
        return None, 'unrepresentable'


def _ledgers(model, w, r2, order):
    """Evaluate the Hall field exactly at ``(w, r^2)``.

    Velocities squared and the common radical factor are retained
    symbolically so all quantities stay rational.

    Returns
    -------
    tuple
        ``(u, v2, speed2, sound2, h, terms, irrotational)``: axial velocity
        ``u/a*``, ``(v/a*)^2``, ``(|V|/a*)^2``, ``(a/a*)^2``, the radial
        factor ``h``, the four terms of the full potential equation and the
        irrotationality defect.

    Raises
    ------
    DomainError
        If the expansion leaves the forward, positive-energy branch
        (``u <= 0`` or ``a^2 <= 0``).

    Notes
    -----
    The local speed of sound follows from the energy equation,
    ``(a/a*)^2 = ((gamma + 1) - (gamma - 1)(|V|/a*)^2)/2``. The four terms are
    those of the steady full potential equation
    ``(a^2 - u^2) u_x - 2 u v u_r + (a^2 - v^2) v_r + sigma a^2 v / r = 0``
    (using ``u_r = v_x``), each divided by a common positive factor [2]_. The
    expansion itself follows [1]_ (planar) and [2]_ (axisymmetric).

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
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
    """Return the field values of :class:`HallThroatResult` for ``inputs``.

    All quantities are evaluated exactly and rounded once on export;
    residuals are rounded upward.
    """
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
    """Hall expansion state at one point, with its equation and wall defects.

    The expansion is that of [1]_ (planar) and [2]_ (axisymmetric).

    Construction re-evaluates the expansion from ``inputs`` and rejects any
    field that does not match.

    Attributes
    ----------
    inputs : HallThroatInputs
        Owned evaluation inputs.
    x_m, radius_m, wall_radius_m : float
        Axial position (throat at 0), radius and local hyperbolic-wall radius
        in m.
    axial_velocity_over_sonic : float
        ``u/a*``.
    radial_velocity_over_sonic : float or None
        ``v/a*``; ``None`` when it is not representable as a finite binary64
        value (see ``radial_velocity_status``).
    radial_velocity_status : str
        ``"finite"`` or ``"unrepresentable"``.
    axial_perturbation_over_sonic : float or None
        ``u/a* - 1``; ``None`` when unrepresentable.
    axial_perturbation_status : str
        ``"finite"`` or ``"unrepresentable"``.
    mach : float
        Local Mach number ``|V|/a``.
    flow_angle_deg : float or None
        ``atan2(v, u)`` in degrees; ``None`` when ``v`` is unrepresentable or
        when the rounded angle is zero although the exact radial velocity is
        not.
    flow_regime : str
        ``"sonic"``, ``"supersonic"`` or ``"subsonic"`` from the exact
        ``|V|/a*``.
    representation_limited : bool
        ``True`` when binary64 export hides an exact distinction (for example
        a rounded unit velocity or Mach number) or a quantity is
        unrepresentable.
    potential_residual : float
        Upper bound of the full potential equation sum divided by ``eps``; it
        measures the truncation defect of the expansion at this point.
    relative_term_balance : float
        ``|sum of terms| / sum of |terms|`` (zero when all terms vanish).
    irrotational_residual : float
        Upper bound of the irrotationality defect ``u_r - v_x`` in expansion
        units.
    wall_slip_velocity_over_sonic : float or None
        Velocity component normal to the wall at wall points; ``None``
        elsewhere or when unrepresentable.
    wall_slip_status : str
        ``"not_wall"``, ``"finite"`` or ``"unrepresentable"``.
    numerical_policy : str
        Identifier of the numerical policy.
    physical_model_qualified, sonic_initialization_complete : bool
        Always ``False``: the expansion is neither validated against hardware
        nor a complete sonic-line initializer.
    approximation_error_bound : None
        Always ``None``; no bound on Hall's approximation error is available.

    Raises
    ------
    DomainError
        If a field is malformed, a qualification claim is made, or the record
        does not reproduce from ``inputs``.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """

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
    """Evaluate the Hall throat expansion at one point.

    Parameters
    ----------
    model : HallThroatModel or mapping
        Throat geometry and gas.
    axial_parameter : float
        Stretched axial coordinate ``w``.
    radius_fraction : float
        Radial coordinate as a fraction of ``r_t`` or of the local wall
        radius.
    radial_basis : str, default "throat"
        ``"throat"`` or ``"wall"``.
    order : int, default 2
        Expansion order, 1 or 2.

    Returns
    -------
    HallThroatResult
        Velocities, Mach number, flow angle and the residuals of the full
        potential equation, irrotationality and wall slip.

    Raises
    ------
    DomainError
        If the inputs are invalid, the point lies outside the wall, or the
        expansion leaves the positive-energy branch.

    Notes
    -----
    The first- and second-order coefficients are those of Hall's expansion
    in planar [1]_ and axisymmetric [2]_ form. The residuals are reported, not
    enforced; they indicate where the truncated expansion departs from the
    full nonlinear potential equation and the wall-tangency condition.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
    inputs = HallThroatInputs(model, axial_parameter, radius_fraction, radial_basis, order)
    return HallThroatResult(inputs=inputs, **_evaluate(inputs))


def _add(*polynomials):
    """Return the sum of coefficient tuples (ascending powers), trailing zeros removed."""
    result = [F(0)] * max(map(len, polynomials))
    for poly in polynomials:
        for j, value in enumerate(poly):
            result[j] += value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def _scale(poly, value):
    """Return ``poly`` with every coefficient multiplied by ``value``."""
    return tuple(value * x for x in poly)


def _multiply(a, b):
    """Return the exact product of two coefficient tuples (ascending powers)."""
    result = [F(0)] * (len(a) + len(b) - 1)
    for j, av in enumerate(a):
        for k, bv in enumerate(b):
            result[j + k] += av * bv
    return _add(result)


def _at(poly, value):
    """Evaluate a coefficient tuple (ascending powers) at ``value`` by Horner's rule."""
    answer = F(0)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def _sonic_polynomial(model, eta, order):
    """Return ``(|V|/a*)^2 - 1`` as an exact polynomial in ``w`` along a wall ray.

    The ray is ``r = eta r_w(x)``, so ``(r/r_t)^2 = eta^2 (1 + eps w^2 / k)``;
    the velocity components are those of the Hall expansion [1]_, [2]_ of the
    requested order. A zero of the polynomial is a sonic point (``|V| = a*``).

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
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
    """Return whether ``poly`` is certified strictly positive on ``[lo, hi]``.

    The polynomial is converted exactly to the Bernstein basis on the
    interval; if every Bernstein coefficient is positive the polynomial is
    positive there. Otherwise the interval is bisected, up to ``depth``
    levels. ``True`` is a proof of positivity; ``False`` means only that no
    certificate was found, never a sampled sign guess.
    """
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
    """Validated inputs for :func:`locate_hall_sonic_point`.

    Attributes
    ----------
    model : HallThroatModel
        Owned copy of the throat model.
    wall_fraction : float
        Ray position ``eta`` in ``[0, 1]``: 0 is the axis and 1 the hyperbolic
        wall.
    order : int, default 2
        Expansion order, 1 or 2.
    max_iterations : int, default 128
        Bisection budget, between 1 and 256.

    Raises
    ------
    DomainError
        If a field is malformed or outside its range.
    """

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
    """Return binary64 values ``(lower, upper)`` with ``lower <= value <= upper``."""
    number = _finite_rational(value, 'Hall root parameter')
    return (nextafter(number, -inf) if F(number) > value else number,
            nextafter(number, inf) if F(number) < value else number)


def _root_bounds(square, sign):
    """Return binary64 bounds on ``sqrt(square)``, negated and swapped when ``sign < 0``."""
    number = _sqrt_output(square, 'Hall root coordinate') if square else 0.
    low = nextafter(number, -inf) if F(number) ** 2 > square else number
    high = nextafter(number, inf) if F(number) ** 2 < square else number
    return (low, high) if sign >= 0 else (-high, -low)


def _radial_bounds(model, eta, lo, hi):
    """Return binary64 bounds on the ray radius ``eta r_w(x)`` for ``w`` in ``[lo, hi]``."""
    eps, k, *_ = _coefficients(model)
    rt2 = F(model.throat_radius_m) ** 2
    smallest = F(0) if lo <= 0 <= hi else min(lo * lo, hi * hi)
    largest = max(lo * lo, hi * hi)
    low = _root_bounds(rt2 * eta * eta * (1 + eps * smallest / k), 1)[0]
    high = _root_bounds(rt2 * eta * eta * (1 + eps * largest / k), 1)[1]
    return low, high


@lru_cache(maxsize=32)
def _locate(inputs):
    """Bracket the sonic root of the Hall expansion on a wall ray.

    The sonic polynomial [1]_, [2]_ must change sign from non-positive at
    ``w = -1`` to non-negative at ``w = 1`` and its derivative must carry a
    positivity certificate there, so the root is unique in the interval. Exact
    rational bisection proceeds until both ends round to the same binary64
    value; the exported endpoints are rounded outward and rechecked.

    Raises
    ------
    ConvergenceError
        With code ``hall_sonic_not_bracketed``,
        ``hall_sonic_branch_unresolved``, ``hall_sonic_not_converged`` or
        ``hall_sonic_output_bracket``.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
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
    """Bracketed sonic point of the Hall expansion on one wall ray.

    The sonic function is evaluated with the expansion of [1]_ and [2]_.

    Construction re-solves the bracket from ``inputs`` and rejects any field
    that does not match.

    Attributes
    ----------
    inputs : HallSonicInputs
        Owned solver inputs.
    point : HallThroatResult
        Expansion state at the binary64 midpoint of the final bracket.
    iterations : int
        Bisection iterations used.
    parameter_lower, parameter_upper : float
        Outward-rounded bounds on the stretched coordinate ``w`` of the root.
    x_lower_m, x_upper_m : float
        Corresponding bounds on the axial position in m.
    radius_lower_m, radius_upper_m : float
        Bounds on the ray radius over the bracket in m.
    sonic_function_lower, sonic_function_upper : float
        Outward-rounded values of ``(|V|/a*)^2 - 1`` at the bracket ends
        (non-positive and non-negative respectively).
    numerical_policy : str
        Identifier of the root-bracketing policy.
    physical_model_qualified, sonic_initialization_complete : bool
        Always ``False``.
    approximation_error_bound : None
        Always ``None``.

    Raises
    ------
    DomainError
        If a field is malformed, a qualification claim is made, or the record
        does not reproduce from ``inputs``.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """

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
    """Bracket the sonic point ``|V|/a* = 1`` of the Hall expansion on a wall ray.

    Parameters
    ----------
    model : HallThroatModel or mapping
        Throat geometry and gas.
    wall_fraction : float
        Ray position in ``[0, 1]``; 0 is the axis and 1 the hyperbolic wall.
    order : int, default 2
        Expansion order, 1 or 2.
    max_iterations : int, default 128
        Bisection budget, between 1 and 256.

    Returns
    -------
    HallSonicResult
        Outward-rounded bracket of the sonic point and the expansion state at
        its midpoint.

    Raises
    ------
    DomainError
        If the inputs are invalid.
    ConvergenceError
        If the root is not bracketed in ``-1 <= w <= 1``, monotonicity cannot
        be certified, the budget is exhausted, or the exported endpoints fail
        to bracket the root.

    Notes
    -----
    The polynomial root certificate bounds the numerical root error only, not
    the approximation error of Hall's expansion [1]_, [2]_ or the true
    transonic flow. No supersonic offset or downstream axis-history
    initialization is supplied.

    References
    ----------
    .. [1] HALL-MCCABE-RM3440. McCabe, planar Hall coefficients and nozzle-design
           context, ARC R&M 3440 (report 1964, issued 1967).
           https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023
    .. [2] KLIEGEL-QUAN-1966. Full potential equations, Cartesian Hall coefficients
           and comparison of transonic approximations, 02874-6002-R000.
           https://ntrs.nasa.gov/citations/19670009614
    """
    inputs = HallSonicInputs(model, wall_fraction, order, max_iterations)
    return HallSonicResult(inputs=inputs, **_locate(inputs))
