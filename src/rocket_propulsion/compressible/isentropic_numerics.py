"""Complete perfect-gas ratio expressions beyond safe binary64 intermediates.

Private Decimal contexts own precision, exponents, flags and traps. These are
arithmetic repairs to the ideal-gas equations, not real-gas applicability limits.

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
from math import isfinite, log

from rocket_propulsion.core.errors import DomainError


def requires_precise_ratios(mach: float, gamma: float) -> bool:
    """Conservative dispatch, not a physical admissibility restriction."""
    if gamma < 1.1 or gamma > 1e6 or mach < 1e-80 or mach > 1e3:
        return True
    energy = 1.0 + 0.5*(gamma-1.0)*mach*mach
    critical = (gamma+1.0)/(2.0*energy)
    # Rounded exponents lose accuracy in proportion to their log arguments.
    # Select complete expressions before that amplification becomes large.
    return max(abs(log(energy)), abs(log(critical)))*gamma/(gamma-1.0) > 16.0


def _context(mach: float, gamma: float) -> Context:
    # Retain 1+small kinetic term even for the smallest positive admitted Mach.
    precision = 90 + max(0, D.from_float(gamma).adjusted()) + 2*abs(D.from_float(mach).adjusted())
    return Context(prec=precision, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999,
                   capitals=1, clamp=0, flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])


def _export_exp(logarithm: D, field: str) -> float:
    # Far-out limits avoid irrelevant huge Decimal exponentiation. The actual
    # binary64 rounding boundary, including subnormal rounding, is checked next.
    if not D(-800) < logarithm < D(800):
        raise DomainError('Isentropic ratio is outside its positive finite output range.',
                          code='numeric_range', field=field)
    result = float(logarithm.exp())
    if not isfinite(result) or result <= 0:
        raise DomainError('Isentropic ratio is outside its positive finite output range.',
                          code='numeric_range', field=field)
    return result


def _logarithms(mach: float, gamma: float):
    m, g = D.from_float(mach), D.from_float(gamma)
    x = 1+(g-1)*m*m/2
    total = x.ln()
    critical = ((g+1)/(2*x)).ln()
    area = ((g+1)/(2*(g-1)))*((2*x)/(g+1)).ln()-m.ln()
    return g, total, critical, area


def precise_area(mach: float, gamma: float) -> float:
    """Evaluate area without overflowing 2*gamma or erasing gamma+1."""
    with localcontext(_context(mach, gamma)):
        return _export_exp(_logarithms(mach, gamma)[3], 'area_critical_ratio')


def precise_ratios(mach: float, gamma: float) -> dict[str, float]:
    """Return checked complete static/total, static/sonic and area ratios."""
    with localcontext(_context(mach, gamma)):
        g, total, critical, area = _logarithms(mach, gamma)
        return {
            'temperature_total_ratio': _export_exp(-total, 'temperature_total_ratio'),
            'pressure_total_ratio': _export_exp(-g/(g-1)*total, 'pressure_total_ratio'),
            'density_total_ratio': _export_exp(-total/(g-1), 'density_total_ratio'),
            'temperature_critical_ratio': _export_exp(critical, 'temperature_critical_ratio'),
            'pressure_critical_ratio': _export_exp(g/(g-1)*critical, 'pressure_critical_ratio'),
            'density_critical_ratio': _export_exp(critical/(g-1), 'density_critical_ratio'),
            'area_critical_ratio': _export_exp(area, 'area_critical_ratio'),
        }
