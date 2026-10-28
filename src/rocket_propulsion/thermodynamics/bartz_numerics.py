"""Directed numerical enclosures for prescribed signed Bartz profiles.

The correlation coefficient and exponents of the Bartz relation [1]_ and the
recovery convention [2]_ are evaluated as exact decimal rationals over exact
binary64 inputs. Fractional powers use exact rational roots where they exist
and otherwise directed logarithm/exponential bounds in decimal arithmetic
[3]_, [4]_. Geometry uses the fixed binary64 value of ``math.pi``. Each
exported float is accepted only when the lower and upper bounds round to the
same value. The bounds concern this numerical model only, not its physical
error or the unsampled wall profile.

References
----------
.. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
       Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
       equations 1, 2 and 4,
       https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
.. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal - Decimal fixed-point and
       floating-point arithmetic", https://docs.python.org/3.12/library/decimal.html
.. [4] PYTHON-3.12-FRACTION. Python 3.12, "fractions - Rational numbers",
       https://docs.python.org/3.12/library/fractions.html
"""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
from math import isfinite, pi
from sys import float_info

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.propulsion.burns.marginal_numerics import (
    _PRECISIONS,
    _Intervals,
    _log_interval,
)
from rocket_propulsion.propulsion.cycle_numerics import _integer_root

_ZERO=(Decimal(0),Decimal(0))
_OVERFLOW=Fraction(float_info.max)+2**970
_UNDERFLOW=Fraction(1,2**1075)
_BARTZ_POLICY='rational_bartz_powers_directed_signed_profiles_linear_radius_flux_integral_v1'


@lru_cache(maxsize=4096)
def _power(value, exponent, precision):
    """Return a Decimal enclosure ``(lower, upper)`` of ``value ** exponent``.

    ``value`` is a positive :class:`Fraction` and ``exponent`` a rational
    :class:`Fraction` of either sign. When the integer roots of the
    numerator and denominator required by the exponent's denominator are
    exact, the power is evaluated exactly; otherwise
    ``exp(exponent * ln(value))`` is bounded outward at ``precision`` digits
    and the lower bound is clipped at zero.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12, "decimal - Decimal fixed-point and
           floating-point arithmetic", https://docs.python.org/3.12/library/decimal.html
    .. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions - Rational numbers",
           https://docs.python.org/3.12/library/fractions.html
    """
    if value == 1 or exponent == 0:return Decimal(1),Decimal(1)
    math=_Intervals(precision)
    n,d=_integer_root(value.numerator,exponent.denominator),_integer_root(value.denominator,exponent.denominator)
    if n is not None and d is not None:
        return math.rational(Fraction(n,d)**exponent.numerator)
    logarithm=math.mul(math.rational(exponent),_log_interval(value,precision))
    return (max(Decimal(0),math.near.next_minus(math.near.exp(logarithm[0]))),
            math.near.next_plus(math.near.exp(logarithm[1])))


def _signed_output(evaluate,label):
    """Return the unique binary64 value of a signed quantity from its enclosures.

    ``evaluate(precision)`` returns Decimal bounds. Precision increases until
    both bounds round to the same nonzero finite float, or both are exactly
    zero. Results that are certainly outside the nonzero finite binary64
    range raise :class:`DomainError`; an invalid enclosure or a value still
    unresolved at the largest precision raises :class:`ConvergenceError`.
    """
    for precision in _PRECISIONS:
        lower,upper=evaluate(precision)
        if not lower.is_finite() or not upper.is_finite() or lower>upper:
            raise ConvergenceError('Invalid signed thermal enclosure.',field=label)
        lo,hi=Fraction(lower),Fraction(upper)
        if lo==hi==0:return 0.
        if lo>=_OVERFLOW or hi<=-_OVERFLOW or 0<lo<=hi<=_UNDERFLOW or -_UNDERFLOW<=lo<=hi<0:
            raise DomainError(f'{label} is outside the nonzero finite output range.',field=label)
        try:a,b=float(lo),float(hi)
        except OverflowError:continue
        if a==b and a!=0 and isfinite(a):return a
    raise ConvergenceError('Thermal expression rounding unresolved within the precision budget.',field=label)


@lru_cache(maxsize=128)
def _profile(inputs,precision):
    """Return Decimal enclosures of every station quantity at one precision.

    For each station the enclosures cover the adiabatic wall temperature,
    the gas-side coefficient, the signed heat flux and the cumulative signed
    heat. With ``C = 1 + (gamma - 1) M^2 / 2`` and
    ``B = (Tw/Tc C + 1) / 2`` the coefficient is [1]_
    ``h = 0.026 cp mu^(1/5) Pr^(-3/5) D^(-1/5) (Pc/c*)^(4/5) (D/Rc)^(1/10)``
    ``Ar^(-9/10) B^(-17/25) C^(-3/25)``; the recovery temperature is
    ``Taw = Tc (1 + Pr^(1/3) (C - 1)) / C`` and the flux ``q = h (Taw - Tw)``
    [2]_. The segment heat is the exact integral of linearly reconstructed
    flux and radius, ``pi L / 3 [(2 r0 + r1) q0 + (r0 + 2 r1) q1]``.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """
    math=_Intervals(precision);f=Fraction
    tc,tw=f(inputs.chamber_temperature_k),f(inputs.wall_temperature_k)
    pr=f(inputs.prandtl_number);diameter=2*f(inputs.contour.throat_radius_m)
    # Exact rational coefficients specify the printed correlation exponents.
    base=math.rational(f(13,500)*f(inputs.specific_heat_j_kg_k))
    for value,exponent in ((f(inputs.dynamic_viscosity_pa_s),f(1,5)),(pr,f(-3,5)),
            (diameter,f(-1,5)),(f(inputs.chamber_pressure_pa)/f(inputs.characteristic_velocity_m_s),f(4,5)),
            (diameter/f(inputs.throat_radius_of_curvature_m),f(1,10))):
        base=math.mul(base,_power(value,exponent,precision))
    root=_power(pr,f(1,3),precision)
    rn,rd=_integer_root(pr.numerator,3),_integer_root(pr.denominator,3)
    exact_root=None if rn is None or rd is None else f(rn,rd)
    records:list[dict[str,tuple[Decimal,Decimal]]]=[];cumulative=_ZERO
    for index,station in enumerate(inputs.contour.stations):
        compression=1+(f(inputs.gamma)-1)*f(station.mach)**2/2
        sigma_base=(tw/tc*compression+1)/2
        h=math.mul(base,_power(f(station.area_ratio),f(-9,10),precision))
        h=math.mul(h,_power(sigma_base,f(-17,25),precision))
        h=math.mul(h,_power(compression,f(-3,25),precision))
        if exact_root is not None:
            taw_exact=tc*(1+exact_root*(compression-1))/compression
            taw=math.rational(taw_exact);difference=math.rational(taw_exact-tw)
        else:
            taw=math.add(math.rational(tc/compression),math.mul(math.rational(tc*(compression-1)/compression),root))
            # Keep subtraction in the underlying expression, before any export.
            difference=math.add(math.rational(tc/compression-tw),math.mul(math.rational(tc*(compression-1)/compression),root))
        flux=math.mul(h,difference)
        if index:
            previous=inputs.contour.stations[index-1];previous_flux=records[-1]['heat_flux_w_m2']
            length_squared=(f(station.x_m)-f(previous.x_m))**2+(f(station.radius_m)-f(previous.radius_m))**2
            length=_power(length_squared,f(1,2),precision)
            # Integral of q(u)*r(u) for linearly reconstructed q and radius.
            a=math.mul(math.rational(2*f(previous.radius_m)+f(station.radius_m)),previous_flux)
            b=math.mul(math.rational(f(previous.radius_m)+2*f(station.radius_m)),flux)
            segment=math.mul(math.rational(f(pi)/3),math.mul(length,math.add(a,b)))
            cumulative=math.add(cumulative,segment)
        records.append({'adiabatic_wall_temperature_k': taw,'gas_side_coefficient_w_m2_k': h,
                            'heat_flux_w_m2': flux,'cumulative_heat_w': cumulative})
    return tuple(records)


def _peak_index(inputs):
    """Return the index of the first station with the largest signed heat flux.

    Comparisons use the unrounded enclosures, so ties in displayed values do
    not decide the location. A station whose prescribed Mach number and area
    ratio equal those of the current peak has the same flux and does not
    replace it. An ordering
    that stays unresolved at every precision raises
    :class:`ConvergenceError`.
    """
    peak=0
    for index in range(1,len(inputs.contour.stations)):
        a,b=inputs.contour.stations[index],inputs.contour.stations[peak]
        if (a.mach,a.area_ratio)==(b.mach,b.area_ratio):continue
        for precision in _PRECISIONS:
            math=_Intervals(precision);profile=_profile(inputs,precision)
            x,y=profile[index]['heat_flux_w_m2'],profile[peak]['heat_flux_w_m2']
            if x==y==_ZERO:break
            lo,hi=math.sub(x,y)
            if lo>0:peak=index;break
            if hi<0:break
        else:raise ConvergenceError('Thermal peak ordering unresolved.',field='Peak location')
    return peak


def _profile_exports(inputs):
    """Return the exported station records, peak, total heat load and sign flag.

    Every exported float is rounded independently from its own enclosure;
    rounded values are never reused in dependent calculations.
    """
    stations=[]
    for index,station in enumerate(inputs.contour.stations):
        record={'x_m': station.x_m,'radius_m': station.radius_m,'mach': station.mach}
        for name in ('adiabatic_wall_temperature_k','gas_side_coefficient_w_m2_k','heat_flux_w_m2','cumulative_heat_w'):
            record[name]=_signed_output(lambda p, *, index=index, name=name:_profile(inputs,p)[index][name],name)
        stations.append(record)
    peak=_peak_index(inputs)
    return {'station_values': tuple(stations),'peak_heat_flux_w_m2': stations[peak]['heat_flux_w_m2'],
                'peak_location_m': stations[peak]['x_m'],'total_heat_load_w': stations[-1]['cumulative_heat_w'],
                'contains_negative_heat_flux': any(s['heat_flux_w_m2']<0 for s in stations)}
