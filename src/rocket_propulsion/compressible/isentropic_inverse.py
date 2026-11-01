"""Complete inverses of supplied perfect-gas static/total and static/sonic ratios.

Internal callers validate finite positive ratios and gamma>1 first. Contexts
are private. Logarithmic bracket tests precede any potentially enormous power.
This numerical bracket is not an experimentally validated gas-model domain.

References
----------
.. [1] NACA-TR-1135. https://ntrs.nasa.gov/citations/19930091059
.. [2] https://docs.python.org/3.12/library/decimal.html
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

from rocket_propulsion.core.errors import DomainError


def _context(value: float, gamma: float) -> Context:
    precision = 90 + max(0, D.from_float(gamma).adjusted()) + 2*abs(D.from_float(value).adjusted())
    return Context(prec=precision, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999,
                   capitals=1, clamp=0, flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])


def ratio_inverse(kind: str, value: float, gamma: float, minimum: float, maximum: float) -> float:
    """Invert one admitted ratio, comparing the complete root to the Mach bracket."""
    with localcontext(_context(value, gamma)):
        g, r = D.from_float(gamma), D.from_float(value)
        inverse_power = ((g-1)/g if kind.startswith('pressure')
                         else g-1 if kind.startswith('density') else D(1))
        log_x = -inverse_power*r.ln()
        if '_critical_' in kind:
            log_x += ((g+1)/2).ln()
        if log_x < 0:
            raise DomainError(f'{kind} exceeds its stagnation (M=0) value.', field='value')
        lower, upper = D.from_float(minimum), D.from_float(maximum)
        log_lower = (1+(g-1)*lower*lower/2).ln()
        log_upper = (1+(g-1)*upper*upper/2).ln()
        if not log_lower <= log_x <= log_upper:
            raise DomainError(f'{kind} corresponds to a Mach number outside [{minimum:g}, {maximum:g}].',
                              code='outside_inverse_range', field='value')
        return float((2*(log_x.exp()-1)/(g-1)).sqrt())
