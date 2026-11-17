"""Complete area and Prandtl-Meyer calculations before export to binary64.

Area inversion uses the logarithmic area relation with private Decimal context.
The PM series factors 1-a^n=(1-a)*(1+...+a^(n-1)) before its leading terms can
cancel. Area inversion uses safeguarded Newton/bisection in a branch bracket.
Computed zeros or bounds exporting to the same binary64 value terminate the
iteration; this is not an empirical validation or global root proof.

A faster first stage returns the same binary64 values for most inputs. A
native estimate proposes the nearest binary64 Mach number; lower-precision
Decimal residuals with explicit error bounds then prove that the root lies
between that value's two half-ULP boundaries. Forward angles are accepted
only when the bounded interval rounds to one value. Any uncertain sign, tie,
nonzero underflow or rejected estimate falls back to the complete algorithm
below, which also performs every admission refusal.

References
----------
.. [1] NACA-TR-1135. https://ntrs.nasa.gov/citations/19930091059
.. [2] PYTHON-3.12-DECIMAL. https://docs.python.org/3.12/library/decimal.html
"""
from decimal import (
    ROUND_HALF_EVEN,
    Context,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    getcontext,
    localcontext,
)
from decimal import Decimal as D
from functools import lru_cache
from math import atan, degrees, inf, isfinite, log, nextafter, sqrt, ulp

from rocket_propulsion.core.errors import ConvergenceError, DomainError

_AREA_INVERSE_ITERATIONS = 512
_PM_INVERSE_ITERATIONS = 512
_PM_SERIES_ITERATIONS = 512
_ATAN_SERIES_ITERATIONS = 2048
# First-stage controls: fewer digits than the complete stage, and an error
# allowance of 10**(_BOUND_MARGIN - precision) times the term magnitudes.
_FAST_PATH = True
_FAST_DIGITS_BELOW_COMPLETE = 60
_BOUND_MARGIN = 8
_ROUNDING_STEPS = 6


def _precision(gamma):
    return 90+max(0, D.from_float(gamma).adjusted())


def _context(gamma):
    """Own every context setting, including mutable decimal.DefaultContext defaults."""
    return Context(prec=_precision(gamma), rounding=ROUND_HALF_EVEN,
                   Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[],
                   traps=[InvalidOperation, DivisionByZero, Overflow])


def _atan(value):
    """Positive Decimal arctangent by half-angle reduction and alternating series."""
    multiplier = 1
    while value > D('.01'):
        value /= 1+(1+value*value).sqrt()
        multiplier *= 2
    power, squared, total = value, value*value, D(0)
    for n in range(_ATAN_SERIES_ITERATIONS):
        previous = total
        total += power/(2*n+1) if n % 2 == 0 else -power/(2*n+1)
        if total == previous:
            return multiplier*total
        power *= squared
    raise ConvergenceError('Decimal arctangent did not settle.', field='prandtl_meyer_deg')


@lru_cache(maxsize=32)
def _decimal_pi(precision):
    """Cache immutable pi by precision; even cache misses own every setting."""
    context = Context(prec=precision, rounding=ROUND_HALF_EVEN,
                      Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[],
                      traps=[InvalidOperation, DivisionByZero, Overflow])
    with localcontext(context):
        return 4*_atan(D(1))


def _pm(mach, gamma, pi_value):
    squared = (mach-1)*(mach+1)
    if not squared:
        return D(0)
    root = squared.sqrt()
    if mach > D('1.01'):
        factor = ((gamma+1)/(gamma-1)).sqrt()
        return (factor*_atan(root/factor)-_atan(root))*180/pi_value
    b = 2/(gamma+1)
    a = 1-b
    power, geometric, a_power = squared*root, D(1), D(1)
    total = D(0)
    for n in range(1, _PM_SERIES_ITERATIONS):
        term = geometric*power/(2*n+1)
        previous = total
        total += term if n % 2 else -term
        if total == previous:
            return b*total*180/pi_value
        power *= squared
        a_power *= a
        geometric += a_power
    raise ConvergenceError('Near-sonic PM series did not settle.', field='prandtl_meyer_deg')


def sonic_pm(mach, gamma):
    """Evaluate PM in degrees with independent decimal context and checked export."""
    if _FAST_PATH:
        fast = _fast_pm(mach, gamma)
        if fast is not None:
            return fast
    with localcontext(_context(gamma)):
        value = _pm(D(mach), D(gamma), _decimal_pi(_precision(gamma)))
        result = float(value)
    if value and result == 0.:
        raise DomainError('PM angle is outside the nonzero output range.',
                          code='numeric_range', field='prandtl_meyer_deg')
    return result


def sonic_area(mach, gamma):
    """Evaluate the complete normalized area expression near M=1."""
    with localcontext(_context(gamma)):
        m, g = D(mach), D(gamma)
        logarithm = ((g+1)/(2*(g-1)))*((2+(g-1)*m*m)/(g+1)).ln()-m.ln()
        return float(logarithm.exp())


def _log_area(m, gamma):
    """Complete log(A/A*) in the active Decimal context."""
    return ((gamma+1)/(2*(gamma-1)))*((2+(gamma-1)*m*m)/(gamma+1)).ln()-m.ln()


def _area_inverse(wanted, gamma, lower, upper, subsonic):
    """Solve log(A/A*) on one branch; private context is owned by sonic_inverse."""
    def value(m):
        return _log_area(m, gamma)

    log_left, log_right = value(lower), value(upper)
    if not min(log_left, log_right) <= wanted <= max(log_left, log_right):
        raise DomainError('Requested area inverse is outside the Mach bracket.',
                          code='outside_inverse_range', field='value')
    mach = (lower*upper).sqrt()
    for _ in range(_AREA_INVERSE_ITERATIONS):
        residual = value(mach)-wanted
        if not residual:
            return float(mach)
        if (residual > 0) == (not subsonic):
            upper = mach
        else:
            lower = mach
        if float(lower) == float(upper):
            return float(lower)
        derivative = (mach*mach-1)/(mach*(1+(gamma-1)*mach*mach/2))
        candidate = mach-residual/derivative if derivative else (lower*upper).sqrt()
        if not lower < candidate < upper or candidate == mach:
            candidate = (lower*upper).sqrt() if upper/lower > 4 else (lower+upper)/2
        mach = candidate
    raise ConvergenceError('Area inverse exhausted its iteration budget.',
                           code='area_inverse_not_converged', field='value')


def _pm_fast_bracket(wanted, gamma, maximum, pi_value):
    """Propose a narrow bracket and keep it only after complete sign checks."""
    target, g = float(wanted), float(gamma)
    if not (1.05 <= g <= 2.0 and 1e-4 <= target <= 100.0):
        return D(1), maximum
    factor = sqrt((g + 1.0) / (g - 1.0))
    lo, hi = 1.0, float(maximum)
    for _ in range(55):
        center = (lo + hi) / 2.0
        if center == lo or center == hi:
            break
        root = sqrt((center - 1.0) * (center + 1.0))
        value = degrees(factor * atan(root / factor) - atan(root))
        if value < target:
            lo = center
        else:
            hi = center
    width = 128.0 * ulp(center)
    proposed_lo = max(D(1), D(center - width))
    proposed_hi = min(maximum, D(center + width))
    lo_exact = proposed_lo if _pm(proposed_lo, gamma, pi_value) <= wanted else D(1)
    hi_exact = proposed_hi if _pm(proposed_hi, gamma, pi_value) >= wanted else maximum
    return lo_exact, hi_exact


def _pm_inverse(wanted, gamma, maximum, pi_value):
    """Monotone PM inversion with complete endpoint admission and safeguarded steps."""
    factor = ((gamma+1)/(gamma-1)).sqrt()
    limit = 180/((gamma-1)*(factor+1))
    if not 0 <= wanted < limit:
        raise DomainError('PM angle is outside its physical interval.', field='value')
    lower, upper = D(1), maximum
    if wanted > _pm(upper, gamma, pi_value):
        raise DomainError('PM angle corresponds to a root outside the inverse bracket.',
                          code='outside_inverse_range', field='value')
    lower, upper = _pm_fast_bracket(wanted, gamma, maximum, pi_value)
    mach = (lower*upper).sqrt()
    for _ in range(_PM_INVERSE_ITERATIONS):
        residual = _pm(mach, gamma, pi_value)-wanted
        if not residual:
            return float(mach)
        if residual > 0:
            upper = mach
        else:
            lower = mach
        if float(lower) == float(upper):
            return float(lower)
        derivative = ((mach-1)*(mach+1)).sqrt()/(mach*(1+(gamma-1)*mach*mach/2))*180/pi_value
        candidate = mach-residual/derivative if derivative else (lower+upper)/2
        if candidate == mach or not lower < candidate < upper:
            candidate = (lower+upper)/2
        mach = candidate
    raise ConvergenceError('PM inverse exhausted its iteration budget.',
                           code='prandtl_meyer_not_converged', field='value')


def sonic_inverse(target, gamma, *, area, minimum_mach, maximum_mach, subsonic=False):
    """Invert an area ratio or PM angle within the existing Mach bracket."""
    if _FAST_PATH:
        fast = (_fast_area_inverse(target, gamma, minimum_mach, maximum_mach, subsonic) if area
                else _fast_pm_inverse(target, gamma, maximum_mach))
        if fast is not None:
            return fast
    with localcontext(_context(gamma)):
        g, wanted = D(gamma), D(target)
        if area:
            lower, upper = (D(minimum_mach), D(1)) if subsonic else (D(1), D(maximum_mach))
            wanted = wanted.ln()
            return _area_inverse(wanted, g, lower, upper, subsonic)
        try:
            pi_value = _decimal_pi(_precision(gamma))
            return _pm_inverse(wanted, g, D(maximum_mach), pi_value)
        except ConvergenceError as error:
            raise ConvergenceError('PM inverse evaluation did not settle.',
                                   code='prandtl_meyer_not_converged', field='value') from error


# --- Bounded first stage ---------------------------------------------------------


def _fast_context(gamma):
    """Own every setting, as _context does, at the first-stage precision."""
    return Context(prec=_precision(gamma)-_FAST_DIGITS_BELOW_COMPLETE, rounding=ROUND_HALF_EVEN,
                   Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[],
                   traps=[InvalidOperation, DivisionByZero, Overflow])


# Midpoints of two binary64 values have finite decimal expansions well below
# this precision, so they are exact.
_EXACT = Context(prec=1200, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999,
                 capitals=1, clamp=0, flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])


def _midpoint(left, right):
    return _EXACT.divide(_EXACT.add(D(left), D(right)), 2)


def _bound(magnitude):
    """Absolute error allowance at the active first-stage precision."""
    return magnitude*D(10)**(_BOUND_MARGIN-getcontext().prec)


def _sign(residual, allowance):
    """Return the certain sign of a bounded residual, or 0 when it is uncertain."""
    if abs(residual) <= allowance:
        return 0
    return 1 if residual > 0 else -1


def _jump(point, residual, slope, value, direction, lower, upper):
    """Next binary64 candidate on the proved side: a Newton step from ``point``.

    The step only proposes a value; it is accepted later by its own boundary
    proofs. It always moves at least one binary64 step toward the root.
    """
    try:
        candidate = float(point-residual/slope(point))
    except (ArithmeticError, OverflowError):
        candidate = value
    step = nextafter(value, direction)
    if direction < 0:
        candidate = min(candidate, step) if isfinite(candidate) else step
    else:
        candidate = max(candidate, step) if isfinite(candidate) else step
    return min(max(candidate, lower), upper)


def _nearest(estimate, normalized, slope, lower, upper):
    """Return the binary64 value nearest the root, or None when it is not proved.

    ``normalized(x)`` returns ``(sign, residual)`` of an increasing residual at
    the exact decimal ``x``; the sign is 0 when the bound leaves it uncertain.
    ``slope(x)`` is its positive derivative. The root rounds to ``value`` when
    it lies strictly between that value's half-ULP boundaries; the caller has
    admitted the bracket ends. Moves follow only a proved side.
    """
    if estimate is None or not isfinite(estimate):
        return None
    value = min(max(estimate, lower), upper)
    for _ in range(_ROUNDING_STEPS):
        if value > lower:
            point = _midpoint(nextafter(value, -inf), value)
            below, residual = normalized(point)
            if below == 0:
                return None
            if below > 0:
                value = _jump(point, residual, slope, value, -inf, lower, upper)
                continue
        if value < upper:
            point = _midpoint(value, nextafter(value, inf))
            above, residual = normalized(point)
            if above == 0:
                return None
            if above < 0:
                value = _jump(point, residual, slope, value, inf, lower, upper)
                continue
        return value
    return None


def _native_pm(mach, gamma):
    """Native PM estimate in degrees; the near-sonic series avoids cancellation."""
    squared = (mach-1.0)*(mach+1.0)
    if squared <= 0.0:
        return 0.0
    root = sqrt(squared)
    if mach > 1.1:
        factor = sqrt((gamma+1.0)/(gamma-1.0))
        return degrees(factor*atan(root/factor)-atan(root))
    b = 2.0/(gamma+1.0)
    a = 1.0-b
    power, geometric, a_power, total = squared*root, 1.0, 1.0, 0.0
    for n in range(1, 200):
        term = geometric*power/(2*n+1)
        previous = total
        total += term if n % 2 else -term
        if total == previous:
            break
        power *= squared
        a_power *= a
        geometric += a_power
    return degrees(b*total)


def _native_root(function, slope, lower, upper, increasing):
    """Safeguarded native Newton/bisection estimate; None when it cannot proceed."""
    try:
        low, high = lower, upper
        mach = 0.5*(low+high) if high < 10.0 else sqrt(low*high)
        for _ in range(200):
            residual = function(mach)
            if residual == 0.0:
                return mach
            if (residual > 0.0) == increasing:
                high = mach
            else:
                low = mach
            derivative = slope(mach)
            candidate = mach-residual/derivative if derivative and isfinite(derivative) else 0.5*(low+high)
            if not low < candidate < high:
                candidate = sqrt(low*high) if high/low > 4.0 else 0.5*(low+high)
            if candidate == mach or high-low <= 4.0*ulp(high):
                return candidate
            mach = candidate
    except (ArithmeticError, ValueError):
        return None
    return mach


def _fast_pm(mach, gamma):
    """Return the binary64 PM angle when a bounded first-stage value proves it."""
    try:
        with localcontext(_fast_context(gamma)):
            m, g = D(mach), D(gamma)
            value = _pm(m, g, _decimal_pi(getcontext().prec))
            magnitude = (((g+1)/(g-1)).sqrt()+1)*90 if m > D('1.01') else 2*abs(value)
            allowance = _bound(magnitude)
            low, high = float(value-allowance), float(value+allowance)
    except (ArithmeticError, ConvergenceError):
        return None
    return low if low == high and low != 0.0 else None


def _fast_pm_inverse(target, gamma, maximum_mach):
    """Nearest binary64 Mach for a PM angle, or None to use the complete stage."""
    if not (1.05 <= gamma <= 1e6 and target > 0.0):
        return None
    try:
        limit = degrees(0.5*3.141592653589793*(sqrt((gamma+1.0)/(gamma-1.0))-1.0))
        # Conservative native admission: anything near the physical limit or
        # the inverse bracket is left to the complete stage and its refusals.
        if not (target < 0.999*limit and target < 0.999*_native_pm(maximum_mach, gamma)):
            return None
        estimate = _native_root(lambda m: _native_pm(m, gamma)-target,
                                lambda m: degrees(sqrt(max(0.0, m*m-1.0))/(m*(1.0+0.5*(gamma-1.0)*m*m))),
                                1.0, maximum_mach, True)
        with localcontext(_fast_context(gamma)):
            g, wanted = D(gamma), D(target)
            pi_value = _decimal_pi(getcontext().prec)
            factor_bound = (((g+1)/(g-1)).sqrt()+1)*90

            def normalized(m):
                value = _pm(m, g, pi_value)
                magnitude = factor_bound if m > D('1.01') else 2*abs(value)
                return _sign(value-wanted, _bound(magnitude+wanted)), value-wanted

            def slope(m):
                return ((m-1)*(m+1)).sqrt()/(m*(1+(g-1)*m*m/2))*180/pi_value

            return _nearest(estimate, normalized, slope, 1.0, maximum_mach)
    except (ArithmeticError, ConvergenceError, ValueError):
        return None


def _native_log_area(mach, gamma):
    return (gamma+1.0)/(2.0*(gamma-1.0))*log((2.0+(gamma-1.0)*mach*mach)/(gamma+1.0))-log(mach)


def _fast_area_inverse(target, gamma, minimum_mach, maximum_mach, subsonic):
    """Nearest binary64 Mach for an area ratio on one branch, or None."""
    if not (1.0 < gamma <= 1e6 and target > 1.0):
        return None
    lower, upper = (minimum_mach, 1.0) if subsonic else (1.0, maximum_mach)
    try:
        wanted_native = log(target)
        endpoint = _native_log_area(minimum_mach if subsonic else maximum_mach, gamma)
        # Conservative native admission; refusals stay with the complete stage.
        if not wanted_native < endpoint-1e-6*(abs(endpoint)+1.0):
            return None
        estimate = _native_root(lambda m: _native_log_area(m, gamma)-wanted_native,
                                lambda m: (m*m-1.0)/(m*(1.0+0.5*(gamma-1.0)*m*m)),
                                lower, upper, not subsonic)
        with localcontext(_fast_context(gamma)):
            g = D(gamma)
            wanted = D(target).ln()
            scale = abs((g+1)/(2*(g-1)))

            def normalized(m):
                # Native magnitudes overestimate generously; the bound only needs scale.
                native = float(m)
                quotient = (2.0+(gamma-1.0)*native*native)/(gamma+1.0)
                magnitude = 2*(scale*D(abs(log(quotient))+1.0)+D(abs(log(native))+1.0)+abs(wanted))
                residual = _log_area(m, g)-wanted
                certain = _sign(residual, _bound(magnitude))
                return (-certain, -residual) if subsonic else (certain, residual)

            def slope(m):
                derivative = (m*m-1)/(m*(1+(g-1)*m*m/2))
                return -derivative if subsonic else derivative

            return _nearest(estimate, normalized, slope, lower, upper)
    except (ArithmeticError, ConvergenceError, ValueError):
        return None

