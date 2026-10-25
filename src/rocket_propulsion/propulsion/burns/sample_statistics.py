"""Exact sample summaries of binary64 inputs with correctly rounded exports.

Sample moments [1]_, linear-interpolation percentiles [2]_ and the Pearson
coefficient [3]_ are evaluated in exact rational arithmetic on the supplied
binary values and rounded once on export. These are empirical arithmetic
contracts: they describe the retained sample and do not establish calibrated
uncertainty or interval coverage.

References
----------
.. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
       of NIST Measurement Results", 1994 (Type A evaluation: arithmetic mean and
       experimental standard deviation).
       https://www.nist.gov/pml/nist-technical-note-1297
.. [2] NUMPY-2.2-LINEAR-QUANTILE. NumPy, ``numpy.quantile``, linear
       (Hyndman-Fan type 7) sample quantile definition.
.. [3] NIST-PEARSON. NIST Dataplot, centered-product correlation definition.
"""
from dataclasses import dataclass
from fractions import Fraction
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real_tuple

from .marginal_numerics import _inverse_grid
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import _sqrt_output

_SUMMARY_POLICY = 'exact_binary_sample_moments_linear_rational_quantiles_and_pearson_v1'
_INTERVAL_KIND = 'iid_normal_formula_proxy_unqualified_coverage'
# Two-sided 95 % normal multiplier, used only in an IID-normal formula proxy.
_NORMAL_PROXY_MULTIPLIER = Fraction(49,25)


@dataclass(frozen=True,slots=True)
class _SampleLedger:
    """Exact sample values with their arithmetic mean and unbiased variance."""

    values: tuple[Fraction,...]
    mean: Fraction
    variance: Fraction


def _sample_ledger(values):
    """Return the exact mean and ``n - 1`` sample variance of finite values.

    Parameters
    ----------
    values : sequence of real
        At least two finite sample values.

    Returns
    -------
    _SampleLedger
        Exact values, mean ``m = sum(x_i) / n`` and sample variance
        ``s**2 = sum((x_i - m)**2) / (n - 1)`` [1]_.

    Raises
    ------
    DomainError
        If fewer than two values are supplied or a value is not finite.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the
           Uncertainty of NIST Measurement Results", 1994 (Type A evaluation).
           https://www.nist.gov/pml/nist-technical-note-1297
    """
    rows=tuple(Fraction(x) for x in real_tuple(values,'Sample values'))
    if len(rows)<2:
        raise DomainError('Sample summaries require at least two values.',field='Sample values')
    mean=sum(rows,Fraction(0))/len(rows)
    variance=sum(((x-mean)**2 for x in rows),Fraction(0))/(len(rows)-1)
    return _SampleLedger(rows,mean,variance)


def _linear_quantile(values,probability):
    """Return the linear-interpolation sample quantile, rounded once.

    With sorted values ``x_0 <= ... <= x_(n-1)``, position ``h = p (n - 1)``,
    ``j = floor(h)`` and ``f = h - j``, the quantile is
    ``(1 - f) x_j + f x_k`` with ``k = min(j + 1, n - 1)`` [1]_.

    Parameters
    ----------
    values : sequence of real
        Nonempty finite sample values.
    probability : real or Fraction
        Probability ``p`` in ``[0, 1]``.

    Returns
    -------
    float
        Correctly rounded binary64 quantile.

    Raises
    ------
    DomainError
        If the sample is empty or ``p`` lies outside ``[0, 1]``.

    References
    ----------
    .. [1] NUMPY-2.2-LINEAR-QUANTILE. NumPy, ``numpy.quantile``, linear
           (Hyndman-Fan type 7) sample quantile definition.
    """
    rows=tuple(sorted(Fraction(x) for x in real_tuple(values,'Percentile values')))
    p=_fraction(probability,'Percentile probability')
    if not rows or not 0<=p<=1:
        raise DomainError('Percentiles need nonempty values and a probability in [0,1].',field='Percentile probability')
    position=p*(len(rows)-1)
    low=position.numerator//position.denominator
    high=min(low+1,len(rows)-1)
    fraction=position-low
    return _finite_rational((1-fraction)*rows[low]+fraction*rows[high],'Sample quantile')


def _affine_sqrt_output(mean,radicand,direction,label):
    """Return ``mean + direction * sqrt(radicand)`` rounded once to binary64.

    Neither the mean nor the square root is rounded separately. The binary64
    result is located by bisection on the float grid, and each candidate is
    accepted or rejected by an exact squared rational comparison, so the final
    rounding cell is certified. The midpoint between the largest finite float and
    the next power of two separates the last finite result from overflow.

    Parameters
    ----------
    mean : real or Fraction
        Affine offset.
    radicand : real or Fraction
        Nonnegative value under the square root.
    direction : {-1, 1}
        Sign applied to the square root.
    label : str
        Output name used in error messages and the error ``field``.

    Returns
    -------
    float
        The correctly rounded finite value.

    Raises
    ------
    DomainError
        If the radicand is negative, ``direction`` is not ``-1`` or ``1``, or
        the result is outside the nonzero finite binary64 range.
    """
    mean,radicand=_fraction(mean),_fraction(radicand)
    if radicand<0 or direction not in (-1,1):
        raise DomainError('Invalid affine square-root expression.',field=label)
    if not radicand:
        return _finite_rational(mean,label)
    def compare(point):
        difference=point-mean
        if direction>0:
            if difference<0:return -1
            residual=difference*difference-radicand
        else:
            if difference>0:return 1
            residual=radicand-difference*difference
        return (residual>0)-(residual<0)
    maximum=float_info.max
    boundary=Fraction(maximum)+2**970
    if compare(boundary)<=0 or compare(-boundary)>=0:
        raise DomainError(f'{label} is outside the finite output range.',field=label)
    if compare(Fraction(maximum))<0:return maximum
    if compare(Fraction(-maximum))>0:return -maximum
    try:hint=float(mean)
    except OverflowError:hint=maximum if mean>0 else -maximum
    try:
        return _inverse_grid(-maximum,maximum,compare,hint)
    except DomainError as error:
        raise DomainError(f'{label} is outside the nonzero finite output range.',field=label) from error


def _pearson(left,right):
    """Return the exact-input Pearson correlation, or ``None`` when undefined.

    ``r = sum((x_i - mx)(y_i - my)) / sqrt(sum((x_i - mx)**2) sum((y_i - my)**2))``
    [1]_. The magnitude is the correctly rounded square root of the exact
    ``r**2``; its sign is the sign of the exact centered cross product.

    Parameters
    ----------
    left, right : sequence of real
        Equal-length samples with at least two values each.

    Returns
    -------
    float or None
        Correlation coefficient in ``[-1, 1]``; ``None`` if either sample is
        constant, because the coefficient is then undefined [2]_.

    Raises
    ------
    DomainError
        If the samples differ in length or have fewer than two values.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot, centered-product correlation definition.
    .. [2] SCIPY-1.14.1-PEARSON. SciPy, ``scipy.stats.pearsonr``, Pearson
           coefficient and undefined constant-input behavior.
    """
    a,b=_sample_ledger(left),_sample_ledger(right)
    if len(a.values)!=len(b.values):
        raise DomainError('Pearson vectors require equal lengths.',field='Pearson values')
    xx=sum(((x-a.mean)**2 for x in a.values),Fraction(0))
    yy=sum(((y-b.mean)**2 for y in b.values),Fraction(0))
    if not xx or not yy:return None
    covariance=sum(((x-a.mean)*(y-b.mean) for x,y in zip(a.values,b.values,strict=True)),Fraction(0))
    if not covariance:return 0.
    magnitude=_sqrt_output(covariance*covariance/(xx*yy),'Pearson correlation')
    return magnitude if covariance>0 else -magnitude


def _relative_exact(left,right):
    """Return ``|left - right| / max(|left|, |right|)``, or zero if both are zero.

    Parameters
    ----------
    left, right : real or Fraction
        Values to compare.

    Returns
    -------
    float
        Relative difference in ``[0, 2]``, rounded once to binary64.
    """
    left,right=_fraction(left),_fraction(right)
    scale=max(abs(left),abs(right))
    return _finite_rational(abs(left-right)/scale if scale else Fraction(0),'Relative sample difference')


def _deviation_relative_change(left_variance,right_variance):
    """Return the relative difference of two standard deviations from variances.

    For variances ``a`` and ``b`` the result is
    ``1 - sqrt(min(a, b) / max(a, b))``, which equals
    ``|sqrt(a) - sqrt(b)| / max(sqrt(a), sqrt(b))`` without rounding either
    square root separately.

    Parameters
    ----------
    left_variance, right_variance : Fraction
        Nonnegative exact variances.

    Returns
    -------
    float
        Relative standard-deviation difference in ``[0, 1]``; zero if both
        variances are zero.
    """
    scale=max(left_variance,right_variance)
    if not scale:return 0.
    ratio=min(left_variance,right_variance)/scale
    return _affine_sqrt_output(Fraction(1),ratio,-1,'Relative deviation difference')
