"""Exact binary-input weighted distributions and checked statistical exports.

Weights are normalized for moments/quantiles; supplied closure diagnostics stay
separate. This is numerical consistency, not calibrated probability evidence.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite, isqrt, ldexp

from rocket_propulsion.core.errors import DomainError

from .numeric_outputs import _finite_rational, _fraction

_PROBABILITY_POLICY = "exact_binary_normalized_joint_weights_and_rational_percentiles_v1"


def _sqrt_output(value, label):
    """Round a nonnegative rational square root to binary64, including subnormals.

    An integer floor brackets the root on the destination grid. Squared exact
    midpoint comparisons implement nearest/even without exporting the variance.
    """
    if value < 0:
        raise DomainError(f"{label} requires nonnegative variance.", field=label)
    if value == 0:
        return 0.0
    numerator, denominator = value.numerator, value.denominator
    exponent = numerator.bit_length() - denominator.bit_length()
    if (
        numerator < denominator << exponent
        if exponent >= 0
        else numerator << -exponent < denominator
    ):
        exponent -= 1
    grid = max(exponent // 2 - 52, -1074)
    if grid >= 0:
        denominator <<= 2 * grid
    else:
        numerator <<= -2 * grid
    root = isqrt(numerator // denominator)
    midpoint = 4 * numerator - (2 * root + 1) ** 2 * denominator
    if midpoint > 0 or midpoint == 0 and root % 2:
        root += 1
    try:
        result = ldexp(float(root), grid)
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if not isfinite(result) or result == 0.0:
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    return result


@dataclass(frozen=True, slots=True)
class _WeightedLedger:
    rows: tuple
    total: Fraction
    mean: Fraction
    variance: Fraction

    def quantile(self, probability):
        probability = _fraction(probability, "Quantile probability")
        if not 0 <= probability <= 1:
            raise DomainError("Quantile probability must lie in [0, 1].")
        cumulative = Fraction(0)
        for value, weight in self.rows:
            cumulative += weight
            if cumulative >= probability * self.total:
                return _finite_rational(value, "Weighted quantile")
        return _finite_rational(self.rows[-1][0], "Weighted quantile")


def _weight_ledger(values):
    rows = []
    for value, weight in values:
        value, weight = _fraction(value), _fraction(weight, "Probability weight")
        if weight < 0:
            raise DomainError(
                "Probability weights must be nonnegative.", field="Probability weight"
            )
        if weight:
            rows.append((value, weight))
    if not rows:
        raise DomainError(
            "A weighted distribution requires positive mass.", field="Probability weight"
        )
    rows = tuple(sorted(rows))
    total = sum((w for x, w in rows), Fraction(0))
    mean = sum((x * w for x, w in rows), Fraction(0)) / total
    variance = sum((w * (x - mean) ** 2 for x, w in rows), Fraction(0)) / total
    return _WeightedLedger(rows, total, mean, variance)


def _quantiles(ledger):
    return tuple(ledger.quantile(p) for p in (Fraction(1, 20), Fraction(1, 2), Fraction(19, 20)))


def _sample_moments(values):
    values = tuple(_fraction(x) for x in values)
    if not values:
        raise DomainError("Sample moments require at least one value.")
    mean = sum(values, Fraction(0)) / len(values)
    variance = (
        sum(((x - mean) ** 2 for x in values), Fraction(0)) / (len(values) - 1)
        if len(values) > 1
        else Fraction(0)
    )
    return _finite_rational(mean, "Sample mean"), _sqrt_output(
        variance, "Sample standard deviation"
    )


def _relative_change(left, right):
    left, right = _fraction(left), _fraction(right)
    scale = max(abs(left), abs(right))
    ratio = abs(left - right) / scale if scale else Fraction(0)
    return _finite_rational(ratio, "Relative sample difference")


def _joint_weight(scenario, conditional):
    return _finite_rational(
        _fraction(scenario) * _fraction(conditional), "Joint realization probability"
    )
