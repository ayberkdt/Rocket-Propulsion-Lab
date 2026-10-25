"""Correctly rounded quantiles of bounded marginal distributions.

A quantile is returned as the binary64 value whose rounding cell contains the
exact inverse CDF. Each candidate is accepted only after the sign of
``CDF(x) - p`` is certified: rational CDFs are compared exactly, and
transcendental CDFs (normal, log-uniform) are enclosed with directed-rounding
decimal intervals whose precision grows until the sign is resolved. A comparison
that stays unresolved within the precision budget raises ``ConvergenceError``.

The standard normal CDF uses ``Phi(z) = 1/2 + G(z) / (2 G(inf))`` with
``G(x) = integral_0^x exp(-t**2 / 2) dt``, evaluated from the positive series of
[1]_; ``G(inf)`` and the normal tails are bounded with the Mills-ratio bound of
[2]_. Triangular and log-uniform inverses follow [3]_ and [4]_; interval
endpoints use directed-rounding decimal contexts [5]_.

Distribution parameters and probabilities are declared numbers; the module does
not validate that a marginal describes measured physics.

References
----------
.. [1] NIST-DLMF-7.6. NIST DLMF, positive error-function power series,
       equation 7.6.2.
.. [2] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
.. [3] NIST-TRIPPF. NIST Dataplot, triangular percent point function.
.. [4] SCIPY-1.14.1-LOGUNIFORM. SciPy, log-uniform CDF and quantile formulas.
.. [5] PYTHON-3.12-DECIMAL. Python 3.12, ``decimal``: precision contexts and
       exact float conversion.
"""

from decimal import ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, Context, Decimal
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
    """Return ``-1``, ``0`` or ``1`` according to the sign of ``value``."""
    return (value > 0)-(value < 0)


def _key(value):
    """Map a binary64 value to an unsigned integer key with the same ordering."""
    bits = unpack('>Q', pack('>d', value))[0]
    return (~bits & _MASK) if bits & _SIGN else bits | _SIGN


def _from_key(value):
    """Invert :func:`_key`, returning the binary64 value for an ordering key."""
    bits = value & ~_SIGN if value & _SIGN else ~value & _MASK
    return unpack('>d', pack('>Q', bits))[0]


def _checked_inverse_result(value, compare):
    """Return ``value`` unless it is non-finite or a zero that is not the exact root.

    Raises
    ------
    DomainError
        If ``value`` is not finite, or is ``0.0`` while the exact root is nonzero
        (the result underflowed).
    """
    if not isfinite(value) or value == 0.0 and compare(Fraction(0)) != 0:
        raise DomainError("Marginal value is outside the nonzero finite output range.", field="Marginal value")
    return value


def _inverse_grid(lower, upper, compare, hint=None):
    """Return the binary64 value nearest the root of a monotone sign oracle.

    ``compare(x)`` returns the sign of ``F(x) - target`` for a nondecreasing
    ``F`` evaluated at an exact rational ``x``. The result is the float in
    ``[lower, upper]`` whose rounding cell contains the root; a root exactly on
    a cell midpoint resolves to the float with an even significand.

    Parameters
    ----------
    lower, upper : float
        Search bounds on the binary64 grid.
    compare : callable
        Certified sign oracle taking a ``Fraction``.
    hint : float, optional
        Approximate root. Up to four neighbouring cells around the hint are
        certified first; no accuracy is assumed of the hint, and the search falls
        back to bisection over the ordered float keys of ``[lower, upper]``.

    Returns
    -------
    float
        Correctly rounded root.

    Raises
    ------
    DomainError
        If the rounded root is not finite or underflows to zero.
    """
    if hint is not None and isfinite(hint) and lower <= hint <= upper:
        # A useful hint can be off by a few ULPs. Certify adjacent cells before
        # falling back to the full grid; no accuracy bound on the hint is assumed.
        for _ in range(4):
            before, after = nextafter(hint, -inf), nextafter(hint, inf)
            left = -1 if not isfinite(before) else compare((_fraction(before)+_fraction(hint))/2)
            right = 1 if not isfinite(after) else compare((_fraction(hint)+_fraction(after))/2)
            if left <= 0 <= right:
                even = unpack('>Q', pack('>d', hint))[0] % 2 == 0
                chosen = before if left == 0 and not even else after if right == 0 and not even else hint
                return _checked_inverse_result(chosen, compare)
            hint = before if left > 0 else after
            if not isfinite(hint) or not lower <= hint <= upper:
                break
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
    """Decimal interval arithmetic with outward (floor/ceiling) rounding.

    Intervals are ``(low, high)`` pairs of ``Decimal``. Lower endpoints are
    computed in a ``ROUND_FLOOR`` context and upper endpoints in a
    ``ROUND_CEILING`` context of the requested precision [1]_, so every result
    encloses the exact value.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12, ``decimal``: precision contexts and
           exact float conversion.
           https://docs.python.org/3.12/library/decimal.html
    """

    def __init__(self, precision):
        kwargs = {"prec": precision, "Emin": -999999999, "Emax": 999999999}
        self.down = Context(rounding=ROUND_FLOOR, **kwargs)
        self.up = Context(rounding=ROUND_CEILING, **kwargs)
        self.near = Context(rounding=ROUND_HALF_EVEN, **kwargs)

    def rational(self, value):
        """Return a decimal interval enclosing the exact rational ``value``."""
        value = _fraction(value)
        n, d = Decimal(value.numerator), Decimal(value.denominator)
        return self.down.divide(n,d), self.up.divide(n,d)

    def add(self, a, b):
        """Return an enclosure of ``a + b``."""
        return self.down.add(a[0],b[0]), self.up.add(a[1],b[1])

    def sub(self, a, b):
        """Return an enclosure of ``a - b``."""
        return self.down.subtract(a[0],b[1]), self.up.subtract(a[1],b[0])

    def mul(self, a, b):
        """Return an enclosure of ``a * b`` for intervals of any sign."""
        if a[0] >= 0 and b[0] >= 0:
            return self.down.multiply(a[0],b[0]), self.up.multiply(a[1],b[1])
        return (min(self.down.multiply(x,y) for x in a for y in b),
                max(self.up.multiply(x,y) for x in a for y in b))

    def div(self, a, b):
        """Return an enclosure of ``a / b`` for a strictly positive divisor.

        Raises
        ------
        ConvergenceError
            If the divisor interval is not strictly positive at this precision.
        """
        if b[0] <= 0:
            raise ConvergenceError("Positive interval denominator is unresolved.", field="Marginal enclosure")
        return self.mul(a, (self.down.divide(Decimal(1),b[1]), self.up.divide(Decimal(1),b[0])))


@lru_cache(maxsize=2048)
def _log_interval(value, precision):
    """Return a decimal interval enclosing ``ln(value)`` for a positive rational."""
    if value == 1:
        return Decimal(0), Decimal(0)
    math = _Intervals(precision)
    lo, hi = math.rational(value)
    return math.near.next_minus(math.near.ln(lo)), math.near.next_plus(math.near.ln(hi))


def _certified_sign(evaluate):
    """Return the certified sign of an interval-valued expression.

    ``evaluate(precision)`` returns an enclosure of the expression. Precision is
    increased through ``_PRECISIONS`` until the enclosure excludes zero or
    collapses to exactly zero.

    Raises
    ------
    ConvergenceError
        If the sign is still unresolved at the largest precision.
    """
    for precision in _PRECISIONS:
        lower, upper = evaluate(precision)
        if lower > 0:return 1
        if upper < 0:return -1
        if lower == upper == 0:return 0
    raise ConvergenceError("Marginal comparison could not be certified within the precision budget.", field="Marginal enclosure")


def _gaussian_series_point(x, precision):
    """Enclose ``G(x) = integral_0^x exp(-t**2 / 2) dt`` with a positive series.

    ``G(x) = exp(-x**2 / 2) * sum(x**(2n+1) / (2n+1)!!)``, the series of [1]_
    rewritten for the unit-variance Gaussian. All terms are positive, so no large
    alternating terms cancel. Once the ratio of successive terms is below one,
    a decreasing geometric majorant bounds the entire omitted tail.

    Parameters
    ----------
    x : Decimal
        Nonnegative evaluation point.
    precision : int
        Target decimal precision of the enclosure.

    Returns
    -------
    tuple of Decimal
        Lower and upper bounds of ``G(x)``, clipped to ``[0, x]``.

    Raises
    ------
    ConvergenceError
        If the remainder bound is not certified within the term budget.

    References
    ----------
    .. [1] NIST-DLMF-7.6. NIST DLMF, positive error-function power series,
           equation 7.6.2.
    """
    if x == 0:return Decimal(0), Decimal(0)
    math = _Intervals(precision+32)
    square = math.mul((x,x),(x,x))
    term = x,x
    total = term
    epsilon = Decimal((0,(1,),-precision-12))
    exponent = math.div(math.sub((Decimal(0),)*2,square),(Decimal(2),)*2)
    exponential = (math.near.next_minus(math.near.exp(exponent[0])),
                   math.near.next_plus(math.near.exp(exponent[1])))
    for index in range(12*(precision+32)+2000):
        ratio = (math.down.divide(square[0],Decimal(2*index+3)),
                 math.up.divide(square[1],Decimal(2*index+3)))
        following = math.mul(term,ratio)
        if ratio[1] < 1:
            remaining = math.up.divide(following[1],math.down.subtract(Decimal(1),ratio[1]))
            if remaining <= math.down.multiply(total[0],epsilon):
                upper = math.up.add(total[1],remaining)
                result = math.mul((total[0],upper),exponential)
                return max(Decimal(0),result[0]), min(x,result[1])
        total = math.add(total,following)
        term = following
    raise ConvergenceError("Gaussian integral series did not certify its remainder.", field="Marginal enclosure")


@lru_cache(maxsize=16)
def _gaussian_limit(precision):
    """Enclose ``G(inf) = sqrt(pi / 2)`` from ``G(R)`` and the Mills tail bound.

    For ``R = isqrt(6 (precision + 14)) + 1`` the tail
    ``integral_R^inf exp(-t**2 / 2) dt <= exp(-R**2 / 2) / R`` [1]_ is below
    ``10**(-precision - 12)``, which is added to the upper endpoint.

    References
    ----------
    .. [1] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
    """
    # Mills: tail(R) <= exp(-R^2/2)/R < 10^(-precision-12).
    radius = isqrt(6*(precision+14))+1
    low, high = _gaussian_series_point(Decimal(radius),precision)
    upper = _Intervals(precision+40).up.add(high,Decimal((0,(1,),-precision-12)))
    return low,upper


@lru_cache(maxsize=2048)
def _gaussian_primitive(value, precision):
    """Enclose ``G(value)`` for a nonnegative rational, saturating at ``G(inf)``.

    ``G`` is increasing, so the series is evaluated at the outward decimal
    endpoints of ``value``. Beyond the Mills radius the enclosure of ``G(inf)``
    from :func:`_gaussian_limit` is returned.
    """
    if value == 0:return Decimal(0),Decimal(0)
    radius = isqrt(6*(precision+14))+1
    if value > radius:return _gaussian_limit(precision)
    # Evaluate the monotone primitive at two exact decimal bounding points.
    math = _Intervals(precision+32)
    lo,hi = math.rational(value)
    low = _gaussian_series_point(lo,precision)[0]
    high = _gaussian_series_point(hi,precision)[1]
    return low,high


class _GaussianProbability:
    """Probability ``Phi(z)`` of a standard normal variate, held symbolically.

    Correlated Gaussian-copula draws are carried as their exact latent ``z`` so
    that a marginal quantile can be certified against ``Phi(z)`` without first
    rounding the probability. ``interval`` encloses
    ``Phi(z) = 1/2 + sign(z) G(|z|) / (2 G(inf))`` [1]_ [2]_.

    References
    ----------
    .. [1] NIST-DLMF-7.6. NIST DLMF, positive error-function power series,
           equation 7.6.2.
    .. [2] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
    """

    def __init__(self, value):
        self.z = _fraction(value)

    def interval(self, precision):
        """Return a decimal enclosure of ``Phi(z)`` clipped to ``[0, 1]``."""
        if self.z == 0:return Decimal('.5'),Decimal('.5')
        math = _Intervals(precision)
        center = (Decimal('.5'),)*2
        primitive = _gaussian_primitive(abs(self.z),precision)
        denominator = math.mul((Decimal(2),)*2,_gaussian_limit(precision))
        ratio = math.div(primitive,denominator)
        result = math.add(center,ratio) if self.z > 0 else math.sub(center,ratio)
        return max(Decimal(0),result[0]), min(Decimal(1),result[1])


def _target_interval(target, precision):
    """Enclose a target probability given as a ``Fraction`` or Gaussian probability."""
    return target.interval(precision) if isinstance(target,_GaussianProbability) else _Intervals(precision).rational(target)


def _compare_rational_probability(value, target):
    """Return the certified sign of ``value - target`` for a rational CDF value."""
    if isinstance(target,Fraction):return _sign(value-target)
    return _certified_sign(lambda precision:_Intervals(precision).sub(
        _Intervals(precision).rational(value),_target_interval(target,precision)))


def _normal_quantile(probability):
    """Return the standard normal quantile ``z`` with ``Phi(z) = p``, correctly rounded.

    By the Mills-ratio bound [1]_ both normal tails beyond ``|z| = 64`` carry
    less probability than the smallest positive binary64 number, so the search
    is confined to ``[-64, 64]``. ``statistics.NormalDist`` provides only the
    starting hint; acceptance is decided by certified enclosures of ``Phi``.

    Parameters
    ----------
    probability : real or Fraction
        Probability in the open interval ``(0, 1)``.

    Returns
    -------
    float
        Correctly rounded standard normal quantile.

    Raises
    ------
    DomainError
        If ``probability`` is not in ``(0, 1)``.
    ConvergenceError
        If a CDF comparison cannot be certified within the precision budget.

    References
    ----------
    .. [1] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
    """
    target = _fraction(probability)
    if not 0 < target < 1:
        raise DomainError("Normal sampling probability must be in (0, 1).",field="Sampling probability")
    def compare(point):
        return _certified_sign(lambda precision:_Intervals(precision).sub(
            _GaussianProbability(point).interval(precision),_Intervals(precision).rational(target)))
    hint = NormalDist().inv_cdf(float(target))
    return _inverse_grid(-64.,64.,compare,hint)


def _hint_float(value, lower, upper):
    """Return ``float(value)`` clamped to ``[lower, upper]``, saturating on overflow."""
    try:hint=float(value)
    except OverflowError:hint=upper if value > 0 else lower
    return min(upper,max(lower,hint))


def _normal_compare(item, target, point):
    """Return the certified sign of ``F(point) - target`` for a truncated normal.

    With mean ``m``, standard deviation ``s``, bounds ``[a, b]``,
    ``alpha = (m - a) / s``, ``beta = (b - m) / s`` and ``z = (x - m) / s``,
    the truncated CDF is
    ``F(x) = (G(alpha) + sign(z) G(|z|)) / (G(alpha) + G(beta))``, where ``G``
    is the Gaussian integral from zero.
    """
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
    """Return the certified sign of ``F(point) - target`` for a log-uniform law.

    ``F(x) = (ln x - ln a) / (ln b - ln a)`` [1]_. For a rational target
    ``p = num / den`` with ``den <= 64`` the comparison is exact:
    ``x**den`` versus ``a**(den - num) * b**num``.

    References
    ----------
    .. [1] SCIPY-1.14.1-LOGUNIFORM. SciPy, log-uniform CDF and quantile formulas.
    """
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
    """Return the correctly rounded quantile of a bounded marginal distribution.

    The returned binary64 value is the float whose rounding cell contains the
    exact inverse CDF at ``probability``. Supported marginals, with lower bound
    ``a``, upper bound ``b``, nominal value (mode or mean) ``c`` and, for the
    truncated normal, standard deviation ``s``:

    ``uniform``
        ``x = (1 - p) a + p b``.
    ``triangular``
        ``F(x) = (x - a)**2 / ((b - a)(c - a))`` for ``x <= c`` and
        ``1 - (b - x)**2 / ((b - a)(b - c))`` otherwise [1]_.
    ``log-uniform``
        ``F(x) = (ln x - ln a) / (ln b - ln a)`` [2]_.
    ``truncated-normal``
        Normal CDF with mean ``c`` and standard deviation ``s``, renormalized
        on ``[a, b]``; the normal integral uses [3]_ and [4]_.

    Parameters
    ----------
    item : UncertainParameter
        Marginal definition supplying ``lower``, ``upper``, ``nominal``,
        ``distribution`` and, for the truncated normal, ``standard_deviation``.
    probability : real, Fraction or _GaussianProbability
        Target probability in ``[0, 1]``, or the latent standard normal value of
        a Gaussian-copula draw.

    Returns
    -------
    float
        Quantile in ``[item.lower, item.upper]``.

    Raises
    ------
    DomainError
        If the distribution is unsupported or the result is not representable.
    ConvergenceError
        If a transcendental CDF comparison cannot be certified.

    Notes
    -----
    ``statistics.NormalDist`` and decimal approximations supply only search
    hints; their rounded values are never used as the target.

    References
    ----------
    .. [1] NIST-TRIPPF. NIST Dataplot, triangular percent point function.
    .. [2] SCIPY-1.14.1-LOGUNIFORM. SciPy, log-uniform CDF and quantile formulas.
    .. [3] NIST-DLMF-7.6. NIST DLMF, positive error-function power series,
           equation 7.6.2.
    .. [4] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
    """
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
