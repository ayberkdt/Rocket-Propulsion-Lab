"""Free-space burn metrics with directed enclosures of the numerical error.

For a tabulated burn with piecewise-linear thrust ``F(t)`` and total mass flow
``q(t)`` the module integrates the external-force-free velocity change
``dv = integral F / m dt`` with ``m(t) = m0 - integral q dt`` [1]_, the
acceleration-weighted centroid time, and the effective exhaust velocity
``c = dv / ln(m0 / mf)`` [1]_.  Integrals use Gauss-Legendre quadrature [2]_
and are enclosed by directed-rounding decimal sums [3]_ over exact rational
inputs [4]_.

The reported bounds concern the declared piecewise-linear binary inputs,
numerical integration and output rounding only.  They are not measurement or
physical-model uncertainty, and gravity, drag and attitude motion are
excluded.

References
----------
.. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
.. [2] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
.. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
       exact float conversion, https://docs.python.org/3.12/library/decimal.html.
.. [4] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
       their exact binary value as a rational.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import ROUND_CEILING, ROUND_FLOOR, Context, Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from math import ceil, cos, inf, log, nextafter, pi

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import coerce_fields, integer, real

from .numeric_outputs import _finite_rational, _fraction
from .tabulated import TabulatedBurnArtifact, _artifact_ledger, _segment_ledger


def gauss_legendre(points: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return the Gauss-Legendre nodes and weights on ``[-1, 1]``.

    Parameters
    ----------
    points : int
        Number of nodes, between 1 and 200.  The value is validated before the
        cached rule is consulted.

    Returns
    -------
    nodes, weights : tuple of float
        Binary64 nodes (the zeros of the Legendre polynomial ``P_n``) and the
        weights ``2 / ((1 - x**2) P_n'(x)**2)`` [1]_, in matching order.

    Raises
    ------
    DomainError
        If ``points`` is not an integer in ``[1, 200]``.
    ConvergenceError
        If a Newton iteration for a node does not converge; no unconverged
        rule is returned.

    References
    ----------
    .. [1] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
    """
    points = integer(points, "Quadrature points")
    if not 1 <= points <= 200:
        raise DomainError("Quadrature points must be between 1 and 200.", field="Quadrature points")
    return _gauss_legendre_rule(points)


@lru_cache(maxsize=16)
def _gauss_legendre_rule(points):
    """Compute the ``points``-node Gauss-Legendre rule by Newton iteration.

    Nodes are the zeros of ``P_n``, found from the three-term recurrence for
    ``P_n`` and its derivative; weights are ``2 / ((1 - x**2) P_n'(x)**2)``
    [1]_.  A ``ConvergenceError`` is raised rather than returning an
    unconverged rule.

    References
    ----------
    .. [1] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
    """

    def legendre(x):
        previous, current = 1.0, x
        for k in range(2, points + 1):
            previous, current = current, ((2 * k - 1) * x * current - (k - 1) * previous) / k
        return current, points * (x * current - previous) / (x * x - 1.0)

    nodes, weights = [], []
    for i in range(1, points + 1):
        x = cos(pi * (i - 0.25) / (points + 0.5))
        for _ in range(100):
            value, derivative = legendre(x)
            step = value / derivative
            x -= step
            if abs(step) < 1e-16:
                break
        else:
            raise ConvergenceError(
                "Gauss-Legendre root did not converge.", field="Quadrature nodes"
            )
        _, derivative = legendre(x)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return tuple(nodes), tuple(weights)


def _decimal(value):
    """Convert an exact rational to ``Decimal`` under the caller's rounding context."""
    return Decimal(value.numerator) / Decimal(value.denominator)


def _mass_intervals(artifact, initial, force_scale, flow_scale):
    """Yield exact sub-intervals on which the vehicle loses at most half its mass.

    Within each linear segment the remaining mass is the exact quadratic
    ``m(x) = m_start - h (q0 x + (q1 - q0) x**2 / 2)`` of the normalized time
    ``x``; a leaf is bisected until ``m(b) >= m(a) / 2``.  Each item is
    ``(start, h, f0, f1, q0, q1, mass)`` with exact rational values.

    Raises
    ------
    DomainError
        If the mass reaches zero or below within the burn.
    """
    mass_start = initial
    for segment in artifact.segments:
        duration = _segment_ledger(segment)[0]
        f0 = _fraction(segment.start_axial_thrust_n) * force_scale
        f1 = _fraction(segment.end_axial_thrust_n) * force_scale
        q0 = _fraction(segment.start_total_mass_flow_kg_s) * flow_scale
        q1 = _fraction(segment.end_total_mass_flow_kg_s) * flow_scale

        def mass_at(x, *, mass_start=mass_start, duration=duration, q0=q0, q1=q1):
            return mass_start - duration * (q0 * x + (q1 - q0) * x * x / 2)

        pending = [(Fraction(0), Fraction(1))]
        while pending:
            a, b = pending.pop()
            left, right = mass_at(a), mass_at(b)
            if right <= 0:
                raise DomainError("The burn reaches nonpositive mass.", field="Final mass")
            if 2 * right < left:
                middle = (a + b) / 2
                pending.extend(((middle, b), (a, middle)))
                continue
            yield (
                _fraction(segment.start_time_s) + duration * a,
                duration * (b - a),
                f0 + (f1 - f0) * a,
                f0 + (f1 - f0) * b,
                q0 + (q1 - q0) * a,
                q0 + (q1 - q0) * b,
                left,
            )
        mass_start = mass_at(Fraction(1))


def _quadrature(interval, points):
    """Return the quadrature estimate of ``dv`` and its time moment on one leaf.

    The weighted sum of ``F(x) / m(x)`` over the Gauss-Legendre nodes [1]_
    uses weights normalized to their exact sum.  Each term is formed exactly
    and only then converted to a 96-digit decimal [2]_; ``F / m`` is never
    evaluated in binary floating point.

    References
    ----------
    .. [1] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
           exact float conversion.
    """
    start, h, f0, f1, q0, q1, mass = interval
    nodes, weights = gauss_legendre(points)
    weight_sum = sum((_fraction(w) for w in weights), Fraction(0))
    with localcontext(Context(prec=96)):
        dv = moment = Decimal(0)
        for node, weight in zip(nodes, weights, strict=True):
            x = (_fraction(node) + 1) / 2
            force = (1 - x) * f0 + x * f1
            remaining = mass - h * (q0 * x + (q1 - q0) * x * x / 2)
            term = h * (_fraction(weight) / weight_sum) * force / remaining
            dv += _decimal(term)
            moment += _decimal(term * (start + h * x))
        return Fraction(dv), Fraction(moment)


@lru_cache(maxsize=512)
def _bernstein_enclosure(interval, precision=96):
    """Enclose ``dv`` and its time moment on one leaf with directed Bernstein sums.

    Notes
    -----
    ``u(x)`` is the consumed mass fraction of the leaf, with degree-two
    Bernstein coefficients ``(0, h q0 / (2 m), h (q0 + q1) / (2 m))``, all in
    ``[0, r]`` with ``r <= 1/2``.  The integrand is expanded as
    ``F / (1 - u) = sum_k F u**k``.  Product coefficients are non-negative, so
    each integral is the mean of its Bernstein coefficients.  The omitted
    integral after ``K`` terms is at most ``integral(F) r**(K + 1) / (1 - r)``;
    the moment tail is additionally bounded by the leaf's global end time.
    Lower and upper sums use floor and ceiling decimal rounding [1]_.

    Returns ``(dv_lower, dv_upper, moment_lower, moment_upper)`` as exact
    rationals.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
           exact float conversion.
    """
    start, h, f0, f1, q0, q1, mass = interval
    if f0 == f1 == 0:
        return (Fraction(0),) * 4
    r = h * (q0 + q1) / (2 * mass)
    r_float = float(r)
    log_r = log(r_float) if r_float else (log(r.numerator) - log(r.denominator) if r else 0.0)
    target = 90 if precision == 96 else (precision - 16) * log(10.0)
    terms = 0 if not r else max(1, ceil(target / -log_r))
    results = []
    for rounding in (ROUND_FLOOR, ROUND_CEILING):
        with localcontext(Context(prec=precision, rounding=rounding)):
            beta = _decimal(h / mass)
            t0, dt = _decimal(start), _decimal(h)
            u1, u2 = _decimal(h * q0 / (2 * mass)), _decimal(r)
            p = [_decimal(f0), _decimal(f1)]
            initial_mean = (p[0] + p[1]) / 2
            total = moment = Decimal(0)
            for k in range(terms + 1):
                degree = len(p) - 1
                integral = sum(p, Decimal(0)) / (degree + 1)
                x_integral = sum(((i + 1) * v for i, v in enumerate(p)), Decimal(0)) / (
                    (degree + 1) * (degree + 2)
                )
                total += integral
                moment += t0 * integral + dt * x_integral
                if k == terms:
                    break
                denominator = (degree + 2) * (degree + 1)
                following = []
                for j in range(degree + 3):
                    value = Decimal(0)
                    if 0 <= j - 1 <= degree:
                        value += p[j - 1] * u1 * 2 * j * (degree + 2 - j)
                    if 0 <= j - 2 <= degree:
                        value += p[j - 2] * u2 * j * (j - 1)
                    following.append(value / denominator)
                p = following
            dv, mom = beta * total, beta * moment
            if rounding == ROUND_CEILING and r:
                tail_denominator = Context(prec=precision, rounding=ROUND_FLOOR).subtract(Decimal(1), u2)
                tail_power = Decimal(1)
                for _ in range(terms + 1):
                    tail_power *= u2
                tail = beta * initial_mean * tail_power / tail_denominator
                dv += tail
                mom += (t0 + dt) * tail
            results.append((Fraction(dv), Fraction(mom)))
    return results[0][0], results[1][0], results[0][1], results[1][1]


def _mass_log_enclosure(initial, remaining, minimum_precision=96):
    """Return directed lower and upper bounds of ``ln(initial / remaining)``.

    The decimal precision grows with the digits of the relative consumption,
    so a tiny consumed fraction is retained before the logarithm.  The bounds
    are taken under floor and ceiling rounding and widened by one unit in the
    last place [1]_.  The logarithm is the mass-ratio term of the ideal rocket
    equation [2]_.

    Raises
    ------
    ConvergenceError
        If the lower bound is not positive.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
           exact float conversion.
    .. [2] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation".
    """
    loss_ratio = (initial - remaining) / initial
    precision = max(
        minimum_precision, 80 + len(str(loss_ratio.denominator)) - len(str(loss_ratio.numerator))
    )
    ratio = initial / remaining
    with localcontext(Context(prec=precision, rounding=ROUND_FLOOR)):
        lower = _decimal(ratio).ln().next_minus()
    with localcontext(Context(prec=precision, rounding=ROUND_CEILING)):
        upper = _decimal(ratio).ln().next_plus()
    if lower <= 0:
        raise ConvergenceError(
            "Mass-log enclosure lost positive consumption.", field="Mass logarithm"
        )
    return Fraction(lower), Fraction(upper)


def _upper_output(value, label):
    """Export a non-negative error bound rounded upward to the next float.

    Positive subnormal bounds are rounded outward with ``math.nextafter`` [1]_.

    Raises
    ------
    DomainError
        If the bound is negative or exceeds the finite float range.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter",
           https://docs.python.org/3.12/library/math.html#math.nextafter.
    """
    if value < 0:
        raise DomainError(f"{label} must be nonnegative.", field=label)
    if value == 0:
        return 0.0
    try:
        result = float(value)
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if _fraction(result) < value:
        result = nextafter(result, inf)
    return real(result, label)


def _reference_integrals(intervals, precision):
    """Sum the directed leaf enclosures exactly, without re-rounding them."""
    rows = tuple(_bernstein_enclosure(interval, precision) for interval in intervals)
    return tuple(sum((row[index] for row in rows), Fraction(0)) for index in range(4))


def _exact_integrals(intervals):
    """Return exact ``(dv, moment)`` when every leaf has a polynomial acceleration.

    Two cases are integrated in closed form: zero flow (constant mass, linear
    thrust) and constant flow with ``F / m`` constant.  Exact results also
    resolve ties at the edge of the float range.  ``None`` means at least one
    leaf needs the numerical path.
    """
    dv = moment = Fraction(0)
    for start, h, f0, f1, q0, q1, mass in intervals:
        if f0 == f1 == 0:
            continue
        if q0 == q1 == 0:
            part = h * (f0 + f1) / (2 * mass)
            local = h * h * (f0 + 2 * f1) / (6 * mass)
        elif q0 == q1 and f0 * (mass - h * q0) == f1 * mass:
            part = h * f0 / mass
            local = part * h / 2
        else:
            return None
        dv += part
        moment += start * part + local
    return dv, moment


def _constant_exhaust(intervals):
    """Return the exact constant ``c`` when ``F = c q`` at every leaf end, else ``None``.

    With a constant effective exhaust velocity the rocket-equation quotient
    ``dv / ln(m0 / mf)`` equals ``c`` exactly [1]_, so no approximate
    logarithmic division is needed.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation".
    """
    candidate = None
    for _, _, f0, f1, q0, q1, _ in intervals:
        for force, flow in ((f0, q0), (f1, q1)):
            if not flow:
                if force:
                    return None
                continue
            ratio = force / flow
            if candidate is not None and candidate != ratio:
                return None
            candidate = ratio
    return candidate


def _range_output(value, lower, upper, label):
    """Export ``value`` after deciding representability from its enclosure.

    Representability is decided from ``[lower, upper]``, not from an
    approximation that may have underflowed or overflowed.  ``None`` asks the
    caller for a tighter reference interval.  A true zero and exact ties are
    resolved by the exact polynomial and constant-exhaust identities.

    Raises
    ------
    DomainError
        If the whole enclosure lies outside the non-zero finite float range.
    """

    def rounded(number):
        try:
            return float(number)
        except OverflowError:
            return inf

    low, high = rounded(lower), rounded(upper)
    if lower == upper == 0:
        return 0.0, Fraction(0)
    if low == high == 0.0 or low == high == inf:
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    if low == 0.0 or high == inf:
        return None
    chosen = rounded(value)
    if chosen == 0.0 or chosen == inf:
        value = (lower + upper) / 2
    return _finite_rational(value, label), value


@dataclass(frozen=True, slots=True)
class FreeSpaceBurnMetrics:
    """Rounded free-space burn metrics with numerical error bounds.

    Parameters
    ----------
    initial_mass_kg, final_mass_kg : float
        Vehicle mass before and after the burn, in kg; positive, with
        ``final_mass_kg <= initial_mass_kg``.
    consumed_propellant_kg : float
        Integrated propellant consumption, in kg.
    delta_v_m_s : float
        External-force-free velocity change ``integral F / m dt`` [1]_, in m/s.
    impulsive_equivalent_time_s : float
        Acceleration-weighted centroid time of the burn, in s.
    effective_exhaust_velocity_m_s : float or None
        ``dv / ln(m0 / mf)`` [1]_, in m/s; ``None`` when nothing is consumed.
    quadrature_points : int
        Gauss-Legendre nodes per leaf [2]_.
    delta_v_error_estimate_m_s : float
        Upper bound of the numerical delta-v error, in m/s.
    error_estimate_kind : str, optional
        Name of the error-bound method.
    integration_subintervals : int, optional
        Number of mass-partition leaves.
    centroid_error_bound_s : float or None, optional
        Upper bound of the numerical centroid error, in s.
    effective_exhaust_velocity_error_bound_m_s : float or None, optional
        Upper bound of the numerical effective-exhaust-velocity error, in m/s.

    Raises
    ------
    DomainError
        If a mass is not positive, the final mass exceeds the initial mass,
        or a metric or bound is negative.

    Notes
    -----
    The bounds cover quadrature and output rounding of the declared
    piecewise-linear inputs only; they are enclosed with directed decimal
    arithmetic [3]_.  They exclude measurement and physical-model
    uncertainty, and the cited sources do not independently validate this
    implementation.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
    .. [2] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
    .. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
           exact float conversion, https://docs.python.org/3.12/library/decimal.html.
    """

    initial_mass_kg: float
    final_mass_kg: float
    consumed_propellant_kg: float
    delta_v_m_s: float
    impulsive_equivalent_time_s: float
    effective_exhaust_velocity_m_s: float | None
    quadrature_points: int
    delta_v_error_estimate_m_s: float
    error_estimate_kind: str = "unspecified"
    integration_subintervals: int = 0
    centroid_error_bound_s: float | None = None
    effective_exhaust_velocity_error_bound_m_s: float | None = None

    def __post_init__(self):
        coerce_fields(self, "Free-space metrics")
        if self.final_mass_kg > self.initial_mass_kg:
            raise DomainError("Final mass cannot exceed initial mass.", field="final_mass_kg")
        for name in ("initial_mass_kg", "final_mass_kg"):
            if getattr(self, name) <= 0:
                raise DomainError(f"{name} must be positive.", field=name)
        for name in (
            "consumed_propellant_kg",
            "delta_v_m_s",
            "impulsive_equivalent_time_s",
            "effective_exhaust_velocity_m_s",
            "delta_v_error_estimate_m_s",
            "centroid_error_bound_s",
            "effective_exhaust_velocity_error_bound_m_s",
            "quadrature_points",
            "integration_subintervals",
        ):
            value = getattr(self, name)
            if value is not None and value < 0:
                raise DomainError(f"{name} must be nonnegative.", field=name)


def free_space_burn_metrics(
    artifact: TabulatedBurnArtifact,
    initial_mass_kg: float,
    *,
    thrust_scale: float = 1.0,
    mass_flow_scale: float = 1.0,
    quadrature_points: int = 16,
) -> FreeSpaceBurnMetrics:
    """Compute free-space delta-v, centroid time and effective exhaust velocity.

    Parameters
    ----------
    artifact : TabulatedBurnArtifact
        Burn with piecewise-linear thrust and total mass flow.
    initial_mass_kg : float
        Vehicle mass at burn start, in kg; positive.
    thrust_scale, mass_flow_scale : float, optional
        Positive multipliers on the tabulated thrust and mass flow.
    quadrature_points : int, optional
        ``n`` in ``[1, 100]``; each leaf uses a normalized ``2n``-point
        Gauss-Legendre rule [2]_.

    Returns
    -------
    FreeSpaceBurnMetrics
        Rounded metrics with directed numerical error bounds.

    Raises
    ------
    DomainError
        If the artifact type, mass, scales or point count are invalid, if the
        burn consumes at least the whole initial mass, or if a result lies
        outside the finite float range.
    ConvergenceError
        If the output range stays unresolved after reference refinement.

    Notes
    -----
    The velocity change is ``dv = integral F / m dt`` with
    ``m(t) = m0 - integral q dt``, and the effective exhaust velocity is
    ``dv / ln(m0 / mf)`` [1]_.  Exact mass decisions and a geometric mass
    partition (each leaf loses at most half its mass) precede quadrature.
    Directed Bernstein reference sums enclose delta-v and its time moment
    independently of the quadrature, with floor/ceiling decimal rounding
    [3]_.  The public error bounds include rule error and output rounding,
    not sensor or physical-model error.  Gravity and other external forces
    are excluded.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
    .. [2] NIST DLMF §3.5(v), "Gauss quadrature", https://dlmf.nist.gov/3.5#v.
    .. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal": precision contexts and
           exact float conversion, https://docs.python.org/3.12/library/decimal.html.
    """
    if not isinstance(artifact, TabulatedBurnArtifact):
        raise DomainError("Burn metrics require a TabulatedBurnArtifact.", field="artifact")
    initial = _fraction(initial_mass_kg, "Initial mass")
    force_scale = _fraction(thrust_scale, "Thrust scale")
    flow_scale = _fraction(mass_flow_scale, "Mass-flow scale")
    if min(initial, force_scale, flow_scale) <= 0:
        raise DomainError("Initial mass and scales must be positive.", field="Mass and scales")
    points = integer(quadrature_points, "Quadrature points")
    if not 1 <= points <= 100:
        raise DomainError("Quadrature points must be between 1 and 100.", field="Quadrature points")
    consumed = _artifact_ledger(artifact.segments)[2] * flow_scale
    remaining = initial - consumed
    if remaining <= 0:
        raise DomainError("The burn consumes at least the whole initial mass.", field="Final mass")
    final = _finite_rational(remaining, "Final mass")
    exported_consumed = _finite_rational(consumed, "Consumed propellant")
    intervals = tuple(_mass_intervals(artifact, initial, force_scale, flow_scale))
    exact = _exact_integrals(intervals)
    if exact is not None:
        dv, moment = exact
    else:
        with localcontext(Context(prec=96)):
            dv_sum = moment_sum = Decimal(0)
            for interval in intervals:
                part, part_moment = _quadrature(interval, 2 * points)
                dv_sum += _decimal(part)
                moment_sum += _decimal(part_moment)
            dv, moment = Fraction(dv_sum), Fraction(moment_sum)
    constant_exhaust = _constant_exhaust(intervals) if consumed else None
    for precision in (96, 192, 384, 768, 1536):
        if exact is None:
            lo, hi, moment_lo, moment_hi = _reference_integrals(intervals, precision)
        else:
            lo = hi = dv
            moment_lo = moment_hi = moment
        dv_result = _range_output(dv, lo, hi, "Free-space delta-v")
        if dv_result is None:
            continue
        exported_dv, used_dv = dv_result
        if used_dv:
            centroid_lo, centroid_hi = moment_lo / hi, moment_hi / lo
            centroid_result = _range_output(
                moment / used_dv, centroid_lo, centroid_hi, "Free-space centroid"
            )
            if centroid_result is None:
                continue
            centroid = centroid_result[0]
        else:
            centroid = 0.0
            centroid_lo = centroid_hi = Fraction(0)
        if consumed:
            log_lo, log_hi = _mass_log_enclosure(initial, remaining, precision)
            effective_lo, effective_hi = lo / log_hi, hi / log_lo
            if constant_exhaust is not None:
                effective_lo = effective_hi = constant_exhaust
            effective_result = _range_output(
                used_dv / ((log_lo + log_hi) / 2),
                effective_lo,
                effective_hi,
                "Effective exhaust velocity",
            )
            if effective_result is None:
                continue
            effective = effective_result[0]
            effective_error = _upper_output(
                max(
                    abs(effective_lo - _fraction(effective)),
                    abs(effective_hi - _fraction(effective)),
                ),
                "Exhaust velocity numerical bound",
            )
        else:
            effective = effective_error = None
        break
    else:
        raise ConvergenceError(
            "Output range remained unresolved after reference refinement.",
            field="Numerical output range",
        )
    dv_error = _upper_output(
        max(abs(lo - _fraction(exported_dv)), abs(hi - _fraction(exported_dv))),
        "Delta-v numerical bound",
    )
    centroid_error = _upper_output(
        max(abs(centroid_lo - _fraction(centroid)), abs(centroid_hi - _fraction(centroid))),
        "Centroid numerical bound",
    )
    return FreeSpaceBurnMetrics(
        float(initial),
        final,
        exported_consumed,
        exported_dv,
        centroid,
        effective,
        2 * points,
        dv_error,
        error_estimate_kind="directed_bernstein_enclosure_including_output_rounding",
        integration_subintervals=len(intervals),
        centroid_error_bound_s=centroid_error,
        effective_exhaust_velocity_error_bound_m_s=effective_error,
    )
