"""Exact binary-input sample summaries and checked algebraic exports.

These are empirical arithmetic contracts, not calibrated uncertainty coverage.
"""
from dataclasses import dataclass
from fractions import Fraction
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real_tuple
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import _sqrt_output
from .marginal_numerics import _inverse_grid

_SUMMARY_POLICY = 'exact_binary_sample_moments_linear_rational_quantiles_and_pearson_v1'
_INTERVAL_KIND = 'iid_normal_formula_proxy_unqualified_coverage'
_NORMAL_PROXY_MULTIPLIER = Fraction(49,25)


@dataclass(frozen=True,slots=True)
class _SampleLedger:
    values: tuple[Fraction,...]
    mean: Fraction
    variance: Fraction


def _sample_ledger(values):
    rows=tuple(Fraction(x) for x in real_tuple(values,'Sample values'))
    if len(rows)<2:
        raise DomainError('Sample summaries require at least two values.',field='Sample values')
    mean=sum(rows,Fraction(0))/len(rows)
    variance=sum(((x-mean)**2 for x in rows),Fraction(0))/(len(rows)-1)
    return _SampleLedger(rows,mean,variance)


def _linear_quantile(values,probability):
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
    """Round mean +/- sqrt(radicand) without exporting either intermediate.

    Squared rational comparisons certify the final binary rounding cell.
    Explicit overflow midpoints distinguish a last finite result from infinity.
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
    left,right=_fraction(left),_fraction(right)
    scale=max(abs(left),abs(right))
    return _finite_rational(abs(left-right)/scale if scale else Fraction(0),'Relative sample difference')


def _deviation_relative_change(left_variance,right_variance):
    scale=max(left_variance,right_variance)
    if not scale:return 0.
    ratio=min(left_variance,right_variance)/scale
    return _affine_sqrt_output(Fraction(1),ratio,-1,'Relative deviation difference')
