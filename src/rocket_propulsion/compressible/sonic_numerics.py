"""Complete area and Prandtl-Meyer calculations before export to binary64.

Area inversion uses the logarithmic area relation with private Decimal context.
The PM series factors 1-a^n=(1-a)*(1+...+a^(n-1)) before its leading terms can
cancel. Area inversion uses safeguarded Newton/bisection in a branch bracket.
Computed zeros or bounds exporting to the same binary64 value terminate the
iteration; this is not an empirical validation or global root proof.

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
    localcontext,
)
from decimal import Decimal as D
from functools import lru_cache
from math import atan, degrees, sqrt, ulp

from rocket_propulsion.core.errors import ConvergenceError, DomainError

_AREA_INVERSE_ITERATIONS = 512
_PM_INVERSE_ITERATIONS = 512
_PM_SERIES_ITERATIONS = 512
_ATAN_SERIES_ITERATIONS = 2048


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


def _area_inverse(wanted, gamma, lower, upper, subsonic):
    """Solve log(A/A*) on one branch; private context is owned by sonic_inverse."""
    def value(m):
        return ((gamma+1)/(2*(gamma-1)))*((2+(gamma-1)*m*m)/(gamma+1)).ln()-m.ln()

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
