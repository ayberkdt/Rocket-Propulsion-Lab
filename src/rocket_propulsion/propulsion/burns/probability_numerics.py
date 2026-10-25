"""Exact weighted distributions of binary64 inputs with checked exports.

Values and probability weights are held as exact rationals. Weights are
normalized by their total for the mean [1]_, variance and quantiles; the closure
of the supplied weights to one is reported separately and is not hidden by the
normalization. Square roots are correctly rounded with integer primitives
[2]_. The results are numerically consistent summaries of the declared
weights, not calibrated probability evidence.

References
----------
.. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, *Weighted Mean*, equation 2-20,
       normalized weighted mean.
.. [2] PYTHON-3.12-MATH. Python 3.12, ``math``: ``isqrt`` and ``ldexp``
       primitives for binary64 statistical exports.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite, isqrt, ldexp

from rocket_propulsion.core.errors import DomainError

from .numeric_outputs import _finite_rational, _fraction

_PROBABILITY_POLICY = "exact_binary_normalized_joint_weights_and_rational_percentiles_v1"


def _sqrt_output(value, label):
    """Return ``sqrt(value)`` correctly rounded to binary64, including subnormals.

    The exact rational is scaled to the binary64 grid of its square root, the
    integer square root brackets the result, and an exact squared comparison
    with the midpoint applies round-half-to-even [1]_. The variance itself is
    never rounded or exported.

    Parameters
    ----------
    value : Fraction
        Nonnegative exact value, typically a variance.
    label : str
        Output name used in error messages and the error ``field``.

    Returns
    -------
    float
        Correctly rounded square root; ``0.0`` for a zero input.

    Raises
    ------
    DomainError
        If ``value`` is negative or the result is outside the nonzero finite
        binary64 range.

    References
    ----------
    .. [1] PYTHON-3.12-MATH. Python 3.12, ``math``: ``isqrt`` and ``ldexp``
           primitives for binary64 statistical exports.
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
    """Sorted exact ``(value, weight)`` rows with normalized moments.

    With ``W = sum(w_i)``, ``mean = sum(w_i x_i) / W`` [1]_ and
    ``variance = sum(w_i (x_i - mean)**2) / W``.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, *Weighted Mean*, equation 2-20,
           normalized weighted mean.
    """

    rows: tuple
    total: Fraction
    mean: Fraction
    variance: Fraction

    def quantile(self, probability):
        """Return the smallest value whose cumulative weight reaches ``p W``.

        This is the lower inverse of the weighted empirical distribution
        function; no interpolation between support points is applied.

        Parameters
        ----------
        probability : real or Fraction
            Probability ``p`` in ``[0, 1]``.

        Returns
        -------
        float
            Support value rounded once to binary64.

        Raises
        ------
        DomainError
            If ``p`` lies outside ``[0, 1]``.
        """
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
    """Build an exact weighted ledger from ``(value, weight)`` pairs.

    Zero-weight rows are dropped; the remaining rows are sorted by value.

    Parameters
    ----------
    values : iterable of (real, real)
        Values and nonnegative probability weights.

    Returns
    -------
    _WeightedLedger
        Rows, total weight ``W``, normalized mean [1]_ and normalized variance.

    Raises
    ------
    DomainError
        If a weight is negative or no row carries positive weight.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, *Weighted Mean*, equation 2-20,
           normalized weighted mean.
    """
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
    """Return the weighted 5th, 50th and 95th percentiles of ``ledger``."""
    return tuple(ledger.quantile(p) for p in (Fraction(1, 20), Fraction(1, 2), Fraction(19, 20)))


def _sample_moments(values):
    """Return the sample mean and ``n - 1`` sample standard deviation.

    ``m = sum(x_i) / n`` and ``s = sqrt(sum((x_i - m)**2) / (n - 1))`` [1]_;
    a single value has ``s = 0``.

    Parameters
    ----------
    values : iterable of real
        At least one finite value.

    Returns
    -------
    tuple of float
        ``(m, s)``, each rounded once to binary64.

    Raises
    ------
    DomainError
        If no value is supplied.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the
           Uncertainty of NIST Measurement Results", 1994 (Type A evaluation).
           https://www.nist.gov/pml/nist-technical-note-1297
    """
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
    """Return ``|left - right| / max(|left|, |right|)``, or zero if both are zero."""
    left, right = _fraction(left), _fraction(right)
    scale = max(abs(left), abs(right))
    ratio = abs(left - right) / scale if scale else Fraction(0)
    return _finite_rational(ratio, "Relative sample difference")


def _joint_weight(scenario, conditional):
    """Return ``P(scenario) * P(draw | scenario)`` rounded once to binary64.

    Raises
    ------
    DomainError
        If the product is outside the nonzero finite binary64 range.
    """
    return _finite_rational(
        _fraction(scenario) * _fraction(conditional), "Joint realization probability"
    )
