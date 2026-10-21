"""Bounded quantiles on the binary64 grid, with directed decimal enclosures.

Transcendental comparisons must certify a sign; unresolved comparisons fail
typed. Inputs/probabilities are declared numbers, not validated physical laws.
"""

from decimal import Context, Decimal, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction
from functools import lru_cache
from math import inf, isfinite, isqrt, nextafter
from statistics import NormalDist
from struct import pack, unpack

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from .numeric_outputs import _finite_rational, _fraction

_SAMPLING_POLICY = "open_binary_strata_latent_gaussian_and_certified_bounded_quantiles_v1"
_PRECISIONS = (48, 96, 192, 384, 768, 1536, 3072)
_SIGN = 1 << 63
_MASK = (1 << 64)-1


def _sign(value):
    return (value > 0)-(value < 0)


def _key(value):
    bits = unpack('>Q', pack('>d', value))[0]
    return (~bits & _MASK) if bits & _SIGN else bits | _SIGN


def _from_key(value):
    bits = value & ~_SIGN if value & _SIGN else ~value & _MASK
    return unpack('>d', pack('>Q', bits))[0]


def _checked_inverse_result(value, compare):
    if not isfinite(value) or value == 0.0 and compare(Fraction(0)) != 0:
        raise DomainError("Marginal value is outside the nonzero finite output range.", field="Marginal value")
    return value


def _inverse_grid(lower, upper, compare, hint=None):
    """Invert a monotone sign oracle; midpoint ties use the even significand."""
    if hint is not None and isfinite(hint) and lower <= hint <= upper:
        before, after = nextafter(hint, -inf), nextafter(hint, inf)
        left = -1 if not isfinite(before) else compare((_fraction(before)+_fraction(hint))/2)
        right = 1 if not isfinite(after) else compare((_fraction(hint)+_fraction(after))/2)
        if left <= 0 <= right:
            even = unpack('>Q', pack('>d', hint))[0] % 2 == 0
            chosen = before if left == 0 and not even else after if right == 0 and not even else hint
            return _checked_inverse_result(chosen, compare)
    a, b = _key(lower), _key(upper)
    while b-a > 1:
        middle = (a+b)//2
        point = _from_key(middle)
        decision = compare(_fraction(point))
        if decision == 0:
            return point
        if decision < 0:
            a = middle
        else:
            b = middle
    lo, hi = _from_key(a), _from_key(b)
    decision = compare((_fraction(lo)+_fraction(hi))/2)
    chosen = hi if decision < 0 else lo if decision > 0 else lo if unpack('>Q',pack('>d',lo))[0] % 2 == 0 else hi
    return _checked_inverse_result(chosen, compare)


class _Intervals:
    def __init__(self, precision):
        kwargs = dict(prec=precision, Emin=-999999999, Emax=999999999)
        self.down = Context(rounding=ROUND_FLOOR, **kwargs)
        self.up = Context(rounding=ROUND_CEILING, **kwargs)
        self.near = Context(rounding=ROUND_HALF_EVEN, **kwargs)

    def rational(self, value):
        value = _fraction(value)
        n, d = Decimal(value.numerator), Decimal(value.denominator)
        return self.down.divide(n,d), self.up.divide(n,d)

    def add(self, a, b):
        return self.down.add(a[0],b[0]), self.up.add(a[1],b[1])

    def sub(self, a, b):
        return self.down.subtract(a[0],b[1]), self.up.subtract(a[1],b[0])

    def mul(self, a, b):
        return (min(self.down.multiply(x,y) for x in a for y in b),
                max(self.up.multiply(x,y) for x in a for y in b))

    def div(self, a, b):
        if b[0] <= 0:
            raise ConvergenceError("Positive interval denominator is unresolved.", field="Marginal enclosure")
        return self.mul(a, (self.down.divide(Decimal(1),b[1]), self.up.divide(Decimal(1),b[0])))


@lru_cache(maxsize=2048)
def _log_interval(value, precision):
    if value == 1:
        return Decimal(0), Decimal(0)
    math = _Intervals(precision)
    lo, hi = math.rational(value)
    return math.near.next_minus(math.near.ln(lo)), math.near.next_plus(math.near.ln(hi))


def _certified_sign(evaluate):
    for precision in _PRECISIONS:
        lower, upper = evaluate(precision)
        if lower > 0:return 1
        if upper < 0:return -1
        if lower == upper == 0:return 0
    raise ConvergenceError("Marginal comparison could not be certified within the precision budget.", field="Marginal enclosure")


def _gaussian_series_point(x, precision):
    """Enclose integral_0^x exp(-t^2/2) dt by its alternating series."""
    if x == 0:return Decimal(0), Decimal(0)
    extra = int(Fraction(x)**2/2)+80
    math = _Intervals(precision+extra)
    square = math.mul((x,x),(x,x))
    term = x,x
    total = Decimal(0),Decimal(0)
    goal = math.down.multiply(min(x,Decimal(1)),Decimal((0,(1,),-precision-12)))
    for index in range(12*(precision+extra)+2000):
        total = math.add(total,term) if index % 2 == 0 else math.sub(total,term)
        numerator = math.mul(square,(Decimal(2*index+1),)*2)
        denominator = (Decimal(2*(index+1)*(2*index+3)),)*2
        ratio = math.div(numerator,denominator)
        following = math.mul(term,ratio)
        if ratio[1] < 1 and following[1] <= goal:
            if index % 2 == 0:
                result = math.down.subtract(total[0],following[1]), total[1]
            else:
                result = total[0], math.up.add(total[1],following[1])
            return max(Decimal(0),result[0]), min(x,result[1])
        term = following
    raise ConvergenceError("Gaussian integral series did not certify its remainder.", field="Marginal enclosure")


@lru_cache(maxsize=16)
def _gaussian_limit(precision):
    # Mills: tail(R) <= exp(-R^2/2)/R < 10^(-precision-12).
    radius = isqrt(6*(precision+14))+1
    low, high = _gaussian_series_point(Decimal(radius),precision)
    upper = _Intervals(precision+40).up.add(high,Decimal((0,(1,),-precision-12)))
    return low,upper


@lru_cache(maxsize=2048)
def _gaussian_primitive(value, precision):
    if value == 0:return Decimal(0),Decimal(0)
    radius = isqrt(6*(precision+14))+1
    if value > radius:return _gaussian_limit(precision)
    # Evaluate the monotone primitive at two exact decimal bounding points.
    math = _Intervals(precision+int(value*value/2)+80)
    lo,hi = math.rational(value)
    low = _gaussian_series_point(lo,precision)[0]
    high = _gaussian_series_point(hi,precision)[1]
    return low,high


class _GaussianProbability:
    def __init__(self, value):
        self.z = _fraction(value)

    def interval(self, precision):
        if self.z == 0:return Decimal('.5'),Decimal('.5')
        math = _Intervals(precision)
        center = (Decimal('.5'),)*2
        primitive = _gaussian_primitive(abs(self.z),precision)
        denominator = math.mul((Decimal(2),)*2,_gaussian_limit(precision))
        ratio = math.div(primitive,denominator)
        result = math.add(center,ratio) if self.z > 0 else math.sub(center,ratio)
        return max(Decimal(0),result[0]), min(Decimal(1),result[1])


def _target_interval(target, precision):
    return target.interval(precision) if isinstance(target,_GaussianProbability) else _Intervals(precision).rational(target)


def _compare_rational_probability(value, target):
    if isinstance(target,Fraction):return _sign(value-target)
    return _certified_sign(lambda precision:_Intervals(precision).sub(
        _Intervals(precision).rational(value),_target_interval(target,precision)))


def _hint_float(value, lower, upper):
    try:hint=float(value)
    except OverflowError:hint=upper if value > 0 else lower
    return min(upper,max(lower,hint))


def _normal_compare(item, target, point):
    lower,upper,mean,sigma=map(_fraction,(item.lower,item.upper,item.nominal,item.standard_deviation))
    if point <= lower:return -1
    if point >= upper:return 1
    alpha=(mean-lower)/sigma;beta=(upper-mean)/sigma;z=(point-mean)/sigma
    if isinstance(target,_GaussianProbability) and z == target.z and alpha == beta:
        return _sign(z)
    def difference(precision):
        math=_Intervals(precision)
        left=_gaussian_primitive(alpha,precision);right=_gaussian_primitive(beta,precision)
        extent=_gaussian_primitive(abs(z),precision)
        mass=math.add(left,extent) if z >= 0 else math.sub(left,extent)
        return math.sub(mass,math.mul(_target_interval(target,precision),math.add(left,right)))
    return _certified_sign(difference)


def _log_compare(lower, upper, point, target):
    if point <= lower:return -1
    if point >= upper:return 1
    if isinstance(target,Fraction) and target.denominator <= 64:
        a,b=target.numerator,target.denominator
        return _sign(point**b-lower**(b-a)*upper**a)
    def difference(precision):
        math=_Intervals(precision)
        low=_log_interval(lower,precision);high=_log_interval(upper,precision)
        value=_log_interval(point,precision)
        return math.sub(math.sub(value,low),math.mul(_target_interval(target,precision),math.sub(high,low)))
    return _certified_sign(difference)


def marginal_value(item, probability):
    lower,upper,mode=map(_fraction,(item.lower,item.upper,item.nominal))
    target=probability if isinstance(probability,_GaussianProbability) else _fraction(probability)
    if isinstance(target,_GaussianProbability) and target.z == 0:target=Fraction(1,2)
    kind=item.distribution.value
    if kind=='uniform' and isinstance(target,Fraction):
        return _finite_rational((1-target)*lower+target*upper,'Marginal value')
    width=upper-lower
    approximate_p=NormalDist().cdf(float(target.z)) if isinstance(target,_GaussianProbability) else float(target)
    hint=_hint_float(lower+_fraction(approximate_p)*width,item.lower,item.upper)
    if kind=='uniform':
        compare=lambda point: -1 if point <= lower else 1 if point >= upper else _compare_rational_probability((point-lower)/width,target)
    elif kind=='triangular':
        if isinstance(target,Fraction):
            if target <= (mode-lower)/width:
                radicand=target*width*(mode-lower)
                compare=lambda point: -1 if point < lower else _sign((point-lower)**2-radicand)
            else:
                radicand=(1-target)*width*(upper-mode)
                compare=lambda point: 1 if point > upper else _sign(radicand-(upper-point)**2)
        else:
            def compare(point):
                if point <= lower:return -1
                if point >= upper:return 1
                cdf=(point-lower)**2/(width*(mode-lower)) if point <= mode else 1-(upper-point)**2/(width*(upper-mode))
                return _compare_rational_probability(cdf,target)
    elif kind=='log-uniform':
        math=_Intervals(80)
        lo=_log_interval(lower,80);hi=_log_interval(upper,80)
        logvalue=math.add(lo,math.mul(math.rational(_fraction(approximate_p)),math.sub(hi,lo)))
        hint=_hint_float(math.near.exp(logvalue[0]),item.lower,item.upper)
        compare=lambda point:_log_compare(lower,upper,point,target)
    elif kind=='truncated-normal':
        if lower+upper == 2*mode and target == Fraction(1,2):return item.nominal
        compare=lambda point:_normal_compare(item,target,point)
        # NormalDist supplies a search hint only; its rounded CDF is not the target.
        if isinstance(target,_GaussianProbability):z=float(target.z)
        else:
            sigma=_fraction(item.standard_deviation)
            try:alpha=float((lower-mode)/sigma);beta=float((upper-mode)/sigma)
            except OverflowError:alpha=-inf;beta=inf
            normal=NormalDist();a=normal.cdf(alpha);b=normal.cdf(beta)
            u=a+float(target)*(b-a)
            z=normal.inv_cdf(u) if 0 < u < 1 else 0.
        hint=_hint_float(mode+_fraction(item.standard_deviation)*_fraction(z),item.lower,item.upper)
    else:
        raise DomainError('Unsupported marginal distribution.')
    return _inverse_grid(item.lower,item.upper,compare,hint)
