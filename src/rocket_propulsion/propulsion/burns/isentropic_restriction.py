"""Shared ideal-gas orifice flow with bounded sonic regime comparisons."""

from decimal import ROUND_CEILING, ROUND_FLOOR, Context, Decimal
from fractions import Fraction
from functools import lru_cache
from math import exp, expm1, log, log1p, sqrt

from rocket_propulsion.core.errors import DomainError

from .probability_numerics import _sqrt_output

_SONIC_PRECISIONS = (80, 160, 320, 640, 1280)


@lru_cache(maxsize=128)
def _critical_ratio_bounds(gamma: float, precision: int) -> tuple[Fraction, Fraction]:
    """Enclose the sonic pressure ratio ``(2/(gamma+1))**(gamma/(gamma-1))``.

    NASA Glenn isentropic pressure relation at Mach one:
    https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html

    Decimal ln and exp round to nearest. Step one representable value
    outward after each to enclose the mathematical result. The intervening
    divisions and products use directed rounding.
    https://docs.python.org/3.12/library/decimal.html
    """
    exact_gamma = Fraction(gamma)
    temperature = (exact_gamma + 1) / 2
    exponent = exact_gamma / (exact_gamma - 1)
    if exponent.denominator == 1 and exponent.numerator <= 64:
        limit = (1 / temperature)**exponent.numerator
        return limit, limit

    floor = Context(prec=precision, rounding=ROUND_FLOOR)
    ceiling = Context(prec=precision, rounding=ROUND_CEILING)
    temperature_numerator = Decimal(temperature.numerator)
    temperature_denominator = Decimal(temperature.denominator)
    exponent_numerator = Decimal(exponent.numerator)
    exponent_denominator = Decimal(exponent.denominator)

    log_lower = floor.next_minus(floor.ln(floor.divide(
        temperature_numerator, temperature_denominator,
    )))
    log_upper = ceiling.next_plus(ceiling.ln(ceiling.divide(
        temperature_numerator, temperature_denominator,
    )))
    power_lower = floor.multiply(log_lower, floor.divide(exponent_numerator, exponent_denominator))
    power_upper = ceiling.multiply(log_upper, ceiling.divide(exponent_numerator, exponent_denominator))
    lower = floor.next_minus(floor.exp(power_upper.copy_negate()))
    upper = ceiling.next_plus(ceiling.exp(power_lower.copy_negate()))
    return Fraction(lower), Fraction(upper)


def _sonic_regime(ratio: Fraction, gamma: float) -> bool:
    """Classify the unrounded pressure ratio. Exact sonic equality is choked."""
    for precision in _SONIC_PRECISIONS:
        lower, upper = _critical_ratio_bounds(gamma, precision)
        if ratio <= lower:
            return True
        if ratio >= upper:
            return False
        if precision == _SONIC_PRECISIONS[0]:
            exact_gamma = Fraction(gamma)
            exponent = exact_gamma / (exact_gamma - 1)
            # A reduced base denominator is at least 2. If its integer power
            # can have the supplied denominator, compare that power exactly.
            if exponent.denominator == 1 and exponent.numerator < ratio.denominator.bit_length():
                return ratio <= (2 / (exact_gamma + 1))**exponent.numerator
    raise DomainError("Pressure ratio could not be separated from the sonic threshold.",
                      field="Sonic pressure ratio")


def _restriction_flow(
    *,
    gamma: float,
    gas_constant: float,
    effective_area: Fraction,
    upstream: Fraction,
    temperature: float,
    downstream: Fraction,
    label: str,
) -> tuple[float, bool]:
    """Apply the isentropic mass-flow law to validated inputs.

    The effective area includes the discharge coefficient. Pressures remain
    exact through head and regime evaluation. Near zero forward head, the
    leading squared flow factor is ``2*(1-ratio)``. Squared flow is kept
    through the final square root to avoid intermediate range loss.
    https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
    """
    if upstream <= downstream or effective_area <= 0:
        return 0.0, False
    ratio = downstream / upstream
    head = 1 - ratio
    if head <= 2.0**-54:
        choked = False
        shape_squared = 2 * head
    else:
        delta_gamma = gamma - 1.0
        choked = _sonic_regime(ratio, gamma)
        if choked:
            sonic_temperature_log = log1p(0.5 * delta_gamma)
            factor = sqrt(gamma) * exp(-(.5 + 1.0 / delta_gamma) * sonic_temperature_log)
            shape_squared = Fraction(factor)**2
        else:
            log_ratio = log1p(-float(head)) if head < 0.5 else log(float(ratio))
            shape_squared = Fraction(
                2.0 * (gamma / delta_gamma) * exp(2.0 * log_ratio / gamma)
                * -expm1((delta_gamma / gamma) * log_ratio)
            )
    area_pressure = effective_area * upstream
    flow_squared = area_pressure**2 * shape_squared / (
        Fraction(gas_constant) * Fraction(temperature)
    )
    return _sqrt_output(flow_squared, label), choked
