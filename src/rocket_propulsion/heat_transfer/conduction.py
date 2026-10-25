"""Heat conduction in walls, bodies and fins: analytic, series and finite-volume models.

The module covers the conduction problems needed around propulsion hardware:

* steady one-dimensional plane, cylindrical and spherical multilayer walls
  with interfacial (contact) conductances, uniform volumetric heat generation,
  temperature-dependent conductivity through the Kirchhoff transform, and
  temperature, flux or convective/radiative surfaces;
* the lumped-capacity transient with its Biot-number check;
* the semi-infinite solid with a step in surface temperature, a step in
  surface flux, a convective surface, and a harmonic surface temperature;
* series (Fourier-Bessel) solutions for convectively cooled slabs, cylinders
  and spheres, including the one-term coefficients;
* a transient finite-volume solver (theta method, Crank-Nicolson included)
  for multilayer walls with temperature-dependent ``k`` and ``rho c``;
* straight, pin and annular fins (efficiency, effectiveness, resistance), with
  the modified Bessel functions implemented here from their DLMF forms.

Conduction theory, examples and tables follow Lienhard and Lienhard [1]_
(``AHTT-V6``; printed page = PDF page - 14). The Kirchhoff transform is the
form stated in [2]_. Bessel and error-function forms are from the NIST DLMF
[3]_ [4]_ [5]_ [6]_, the annular-fin solution from Yovanovich et al. [7]_ and
the fin-tip length correction from NACA Report 158 [8]_. The Stefan-Boltzmann
constant is the CODATA 2022 value [9]_.

All quantities are SI. Temperatures are absolute (K). Models are textbook
one-dimensional idealizations; no result is a qualified hardware prediction.

References
----------
.. [1] AHTT-V6. J. H. Lienhard IV and J. H. Lienhard V, "A Heat Transfer
       Textbook", 6th ed., version 6.00 (19 April 2024), Phlogiston Press,
       Chs. 1, 2, 4 and 5, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
.. [2] ENCIT-2012-KIRCHHOFF. W. F. da Silva, E. D. Correa, R. M. Saldanha da
       Gama and M. L. Martins-Costa, ENCIT 2012 paper ENCIT2012-352, Sec. 3,
       Eq. (5), https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
.. [3] DLMF-10-BESSEL. NIST DLMF Ch. 10, Eqs. 10.2.2, 10.9.2, 10.17.1,
       10.25.2, 10.28.2, 10.29.3, 10.31.1-2, 10.32.9, 10.40.1-2,
       https://dlmf.nist.gov/10
.. [4] DLMF-7-ERFC. NIST DLMF Secs. 7.12 and 7.18, https://dlmf.nist.gov/7.18
.. [5] DLMF-5-GAMMA. NIST DLMF Eqs. 5.2.3 and 5.4.14, https://dlmf.nist.gov/5.4
.. [6] DLMF-3.5-TRAPEZOID. NIST DLMF Sec. 3.5(i), Eqs. 3.5.4-3.5.5,
       https://dlmf.nist.gov/3.5
.. [7] YOVANOVICH-1988-ANNULAR-FIN. M. M. Yovanovich, J. R. Culham and
       T. F. Lemczyk, J. Thermophysics 2(2), 152-157, 1988, Eqs. (26)-(41),
       https://uwaterloo.ca/microelectronics-heat-transfer-laboratory/sites/default/files/2024-03/mhtl88-1.pdf
.. [8] NACA-TR-158. D. R. Harper 3d and W. B. Brown, NACA Report 158, 1923,
       printed pp. 21 and 26,
       https://ntrs.nasa.gov/api/citations/19930091223/downloads/19930091223.pdf
.. [9] NIST-CODATA-2022-STEFAN-BOLTZMANN. NIST CODATA 2022, Stefan-Boltzmann
       constant, https://physics.nist.gov/cgi-bin/cuu/Value?sigma
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import lru_cache

from rocket_propulsion.core.constants import STEFAN_BOLTZMANN_W_M2_K4
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

# Euler's constant, DLMF Eq. 5.2.3: 0.57721 56649 01532 86060...
_EULER_GAMMA = 0.5772156649015329

CONDUCTION_MODEL = "one_dimensional_conduction_v1"
_SHAPES = ("plane", "cylinder", "sphere")
_SERIES_SHAPES = ("slab", "cylinder", "sphere")


# ---------------------------------------------------------------------------
# validation helpers
# ---------------------------------------------------------------------------

def _finite(value: object, label: str) -> float:
    return real(value, label)


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _nonnegative(value: object, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must not be negative.", field=label)
    return number


def _temperature(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be an absolute temperature above 0 K.", field=label)
    return number


def _emissivity(value: object, label: str) -> float:
    number = real(value, label)
    if not 0.0 <= number <= 1.0:
        raise DomainError(f"{label} must lie in [0, 1].", field=label)
    return number


# ---------------------------------------------------------------------------
# special functions
# ---------------------------------------------------------------------------

def _asymptotic_coefficients(order: int, count: int) -> list[float]:
    """Return a_k(order) for k < count from DLMF Eq. 10.17.1 (product form)."""
    values = [1.0]
    four_nu2 = 4.0 * order * order
    for k in range(1, count):
        values.append(values[-1] * (four_nu2 - (2 * k - 1) ** 2) / (8.0 * k))
    return values


def _i_series(order: int, x: float) -> float:
    """Unscaled I_order(x) from the power series, DLMF Eq. 10.25.2."""
    quarter = 0.25 * x * x
    term = (0.5 * x) ** order / math.factorial(order)
    total = term
    k = 0
    while True:
        k += 1
        term *= quarter / (k * (k + order))
        total += term
        if term <= 1.0e-17 * total:
            return total


def _scaled_i(order: int, x: float) -> float:
    """``exp(-x) I_order(x)`` for x >= 0 (series, or DLMF 10.40.1 for x > 50)."""
    if x <= 50.0:
        return _i_series(order, x) * math.exp(-x)
    coefficients = _asymptotic_coefficients(order, 60)
    total = 0.0
    previous = math.inf
    for k, a_k in enumerate(coefficients):
        term = (-1) ** k * a_k / x**k
        if abs(term) >= previous:
            break
        total += term
        previous = abs(term)
        if abs(term) <= 1.0e-17 * abs(total):
            break
    return total / math.sqrt(2.0 * math.pi * x)


def _k_series(order: int, x: float) -> float:
    """Unscaled K_0 or K_1 from DLMF Eqs. 10.31.1-10.31.2 (used for x <= 2).

    ``psi(k + 1) = H_k - gamma`` with ``H_k`` the harmonic number, DLMF 5.4.14.
    """
    quarter = 0.25 * x * x
    log_half = math.log(0.5 * x)
    if order == 0:
        total = -(log_half + _EULER_GAMMA) * _i_series(0, x)
        term = 1.0
        harmonic = 0.0
        k = 0
        while True:
            k += 1
            term *= quarter / (k * k)
            harmonic += 1.0 / k
            total += harmonic * term
            if term * harmonic <= 1.0e-17 * abs(total):
                return total
    # n = 1: 1/x + ln(x/2) I_1(x) - (x/4) sum (psi(k+1) + psi(k+2)) (x^2/4)^k / (k! (k+1)!)
    series = 0.0
    term = 1.0
    harmonic = 0.0  # H_k
    k = 0
    while True:
        psi_sum = (harmonic - _EULER_GAMMA) + (harmonic + 1.0 / (k + 1) - _EULER_GAMMA)
        series += psi_sum * term
        k += 1
        harmonic += 1.0 / k
        term *= quarter / (k * (k + 1))
        if abs(term) * (abs(harmonic) + 3.0) <= 1.0e-17 * max(abs(series), 1e-300):
            break
    return 1.0 / x + log_half * _i_series(1, x) - 0.25 * x * series


def _scaled_k_integral(order: int, x: float) -> float:
    """``exp(x) K_order(x)`` from DLMF Eq. 10.32.9 by the trapezoidal rule (DLMF 3.5.5).

    The integrand ``exp(-x (cosh t - 1)) cosh(order t)`` is entire and even, so
    the trapezoidal sum over the real line converges exponentially in 1/h; the
    discretization error behaves like ``exp(x - pi^2/h)`` (strip half-width
    pi/2), negligible for the step used with x <= 50.
    """
    step = 0.05
    total = 0.5  # half weight at t = 0 where the integrand is 1
    t = 0.0
    while True:
        t += step
        exponent = -2.0 * x * math.sinh(0.5 * t) ** 2 + order * t
        if exponent < -745.0:
            break
        total += math.exp(-2.0 * x * math.sinh(0.5 * t) ** 2) * math.cosh(order * t)
    return step * total


def _scaled_k(order: int, x: float) -> float:
    """``exp(x) K_order(x)`` for x > 0."""
    if x <= 2.0:
        return _k_series(order, x) * math.exp(x)
    if x <= 50.0:
        return _scaled_k_integral(order, x)
    coefficients = _asymptotic_coefficients(order, 80)
    total = 0.0
    previous = math.inf
    for k, a_k in enumerate(coefficients):
        term = a_k / x**k
        if abs(term) >= previous:
            break
        total += term
        previous = abs(term)
        if abs(term) <= 1.0e-17 * abs(total):
            break
    return math.sqrt(math.pi / (2.0 * x)) * total


def _bessel_argument(x: object, label: str, *, positive: bool) -> float:
    value = real(x, label)
    if positive and value <= 0.0:
        raise DomainError(f"{label} must be positive for K functions.", field=label)
    if not positive and value < 0.0:
        raise DomainError(f"{label} must not be negative.", field=label)
    return value


def bessel_i0e(x: float) -> float:
    """Exponentially scaled modified Bessel function ``exp(-x) I_0(x)``, x >= 0.

    Parameters
    ----------
    x : float
        Non-negative real argument.

    Returns
    -------
    float
        ``exp(-x) I_0(x)``.

    Raises
    ------
    DomainError
        If ``x`` is negative or not finite.

    Notes
    -----
    Power series DLMF 10.25.2 for x <= 50 (all terms positive, no
    cancellation) and the Hankel expansion DLMF 10.40.1, truncated at its
    smallest term, beyond [1]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.25.2 and 10.40.1,
           https://dlmf.nist.gov/10.25
    """
    return _scaled_i(0, _bessel_argument(x, "x", positive=False))


def bessel_i1e(x: float) -> float:
    """Exponentially scaled modified Bessel function ``exp(-x) I_1(x)``, x >= 0.

    Parameters
    ----------
    x : float
        Non-negative real argument.

    Returns
    -------
    float
        ``exp(-x) I_1(x)``.

    Raises
    ------
    DomainError
        If ``x`` is negative or not finite.

    Notes
    -----
    Same method as :func:`bessel_i0e` [1]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.25.2 and 10.40.1,
           https://dlmf.nist.gov/10.25
    """
    return _scaled_i(1, _bessel_argument(x, "x", positive=False))


def bessel_k0e(x: float) -> float:
    """Exponentially scaled modified Bessel function ``exp(x) K_0(x)``, x > 0.

    Parameters
    ----------
    x : float
        Positive real argument.

    Returns
    -------
    float
        ``exp(x) K_0(x)``.

    Raises
    ------
    DomainError
        If ``x`` is not positive and finite.

    Notes
    -----
    Series DLMF 10.31.2 for x <= 2, the integral DLMF 10.32.9 evaluated by
    the trapezoidal rule (exponentially convergent, DLMF 3.5.5) for
    2 < x <= 50, and DLMF 10.40.2 truncated at its smallest term for x > 50,
    where DLMF bounds the remainder by the first omitted term [1]_ [2]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.31.2, 10.32.9 and 10.40.2,
           https://dlmf.nist.gov/10.31
    .. [2] DLMF-3.5-TRAPEZOID. NIST DLMF Sec. 3.5(i), Eq. 3.5.5,
           https://dlmf.nist.gov/3.5
    """
    return _scaled_k(0, _bessel_argument(x, "x", positive=True))


def bessel_k1e(x: float) -> float:
    """Exponentially scaled modified Bessel function ``exp(x) K_1(x)``, x > 0.

    Parameters
    ----------
    x : float
        Positive real argument.

    Returns
    -------
    float
        ``exp(x) K_1(x)``.

    Raises
    ------
    DomainError
        If ``x`` is not positive and finite.

    Notes
    -----
    Series DLMF 10.31.1 (n = 1, with ``psi(k+1) = H_k - gamma`` from DLMF
    5.4.14) for x <= 2; otherwise as :func:`bessel_k0e` [1]_ [2]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.31.1, 10.32.9 and 10.40.2,
           https://dlmf.nist.gov/10.31
    .. [2] DLMF-5-GAMMA. NIST DLMF Eqs. 5.2.3 and 5.4.14,
           https://dlmf.nist.gov/5.4
    """
    return _scaled_k(1, _bessel_argument(x, "x", positive=True))


def _unscale(scaled: float, exponent: float, name: str) -> float:
    if exponent > 709.0:
        raise DomainError(f"{name} overflows binary64 at this argument; use the scaled form.",
                          field="x")
    return scaled * math.exp(exponent)


def bessel_i0(x: float) -> float:
    """Modified Bessel function of the first kind and order zero, ``I_0(x)``, x >= 0.

    Parameters
    ----------
    x : float
        Non-negative real argument (at most about 709 before overflow).

    Returns
    -------
    float
        ``I_0(x)``.

    Raises
    ------
    DomainError
        If ``x`` is negative or the value overflows binary64.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eq. 10.25.2, https://dlmf.nist.gov/10.25
    """
    value = _bessel_argument(x, "x", positive=False)
    return _unscale(_scaled_i(0, value), value, "I_0")


def bessel_i1(x: float) -> float:
    """Modified Bessel function of the first kind and order one, ``I_1(x)``, x >= 0.

    Parameters
    ----------
    x : float
        Non-negative real argument (at most about 709 before overflow).

    Returns
    -------
    float
        ``I_1(x)``.

    Raises
    ------
    DomainError
        If ``x`` is negative or the value overflows binary64.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eq. 10.25.2, https://dlmf.nist.gov/10.25
    """
    value = _bessel_argument(x, "x", positive=False)
    return _unscale(_scaled_i(1, value), value, "I_1")


def bessel_k0(x: float) -> float:
    """Modified Bessel function of the second kind and order zero, ``K_0(x)``, x > 0.

    Parameters
    ----------
    x : float
        Positive real argument.

    Returns
    -------
    float
        ``K_0(x)`` (underflows to 0 for very large x).

    Raises
    ------
    DomainError
        If ``x`` is not positive and finite.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.31.2, 10.32.9 and 10.40.2,
           https://dlmf.nist.gov/10.31
    """
    value = _bessel_argument(x, "x", positive=True)
    return _scaled_k(0, value) * math.exp(-value)


def bessel_k1(x: float) -> float:
    """Modified Bessel function of the second kind and order one, ``K_1(x)``, x > 0.

    Parameters
    ----------
    x : float
        Positive real argument.

    Returns
    -------
    float
        ``K_1(x)`` (underflows to 0 for very large x).

    Raises
    ------
    DomainError
        If ``x`` is not positive and finite.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eqs. 10.31.1, 10.32.9 and 10.40.2,
           https://dlmf.nist.gov/10.31
    """
    value = _bessel_argument(x, "x", positive=True)
    return _scaled_k(1, value) * math.exp(-value)


def _bessel_j(order: int, x: float) -> float:
    """J_order(x) from Bessel's integral DLMF 10.9.2, trapezoidal rule over one period.

    ``J_n(x) = (1/2 pi) int_0^{2 pi} cos(x sin t - n t) dt``; the M-point
    trapezoidal sum of this periodic entire integrand differs from J_n by
    aliased terms J_{n +- pM}(x), which are below 1e-17 once M exceeds x by a
    margin of order x^(1/3) (here 12 x^(1/3) + 32).
    """
    magnitude = abs(x)
    points = int(magnitude + 12.0 * magnitude ** (1.0 / 3.0) + 32.0)
    total = 0.0
    for j in range(points):
        angle = 2.0 * math.pi * j / points
        total += math.cos(x * math.sin(angle) - order * angle)
    return total / points


def bessel_j0(x: float) -> float:
    """Bessel function of the first kind and order zero, ``J_0(x)``.

    Parameters
    ----------
    x : float
        Real argument (|x| up to about 1e4 for the stated accuracy).

    Returns
    -------
    float
        ``J_0(x)``; absolute error of order 1e-15.

    Raises
    ------
    DomainError
        If ``x`` is not finite or |x| > 1e4.

    Notes
    -----
    Bessel's integral DLMF 10.9.2 evaluated with the trapezoidal rule over a
    full period, which converges exponentially for periodic analytic
    integrands (DLMF 3.5.4) [1]_ [2]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eq. 10.9.2, https://dlmf.nist.gov/10.9
    .. [2] DLMF-3.5-TRAPEZOID. NIST DLMF Sec. 3.5(i), Eq. 3.5.4,
           https://dlmf.nist.gov/3.5
    """
    value = real(x, "x")
    if abs(value) > 1.0e4:
        raise DomainError("|x| above 1e4 is outside the supported range of J_0.", field="x")
    return _bessel_j(0, value)


def bessel_j1(x: float) -> float:
    """Bessel function of the first kind and order one, ``J_1(x)``.

    Parameters
    ----------
    x : float
        Real argument (|x| up to about 1e4 for the stated accuracy).

    Returns
    -------
    float
        ``J_1(x)``; absolute error of order 1e-15.

    Raises
    ------
    DomainError
        If ``x`` is not finite or |x| > 1e4.

    Notes
    -----
    As :func:`bessel_j0`, with n = 1 in DLMF 10.9.2 [1]_ [2]_.

    References
    ----------
    .. [1] DLMF-10-BESSEL. NIST DLMF Eq. 10.9.2, https://dlmf.nist.gov/10.9
    .. [2] DLMF-3.5-TRAPEZOID. NIST DLMF Sec. 3.5(i), Eq. 3.5.4,
           https://dlmf.nist.gov/3.5
    """
    value = real(x, "x")
    if abs(value) > 1.0e4:
        raise DomainError("|x| above 1e4 is outside the supported range of J_1.", field="x")
    return _bessel_j(1, value)


def _erfcx(z: float) -> float:
    """``exp(z^2) erfc(z)`` for real z; DLMF 7.12.1 for z > 25 (remainder < first omitted term)."""
    if z <= 25.0:
        return math.exp(z * z) * math.erfc(z)
    inverse = 1.0 / (z * z)
    total = 0.0
    term = 1.0
    m = 0
    while True:
        total += term
        m += 1
        new_term = -term * (m - 0.5) * inverse
        if abs(new_term) >= abs(term) or abs(new_term) <= 1.0e-17 * abs(total):
            break
        term = new_term
    return total / (z * math.sqrt(math.pi))


def integrated_erfc(z: float) -> float:
    """First repeated integral of erfc, ``i erfc(z) = exp(-z^2)/sqrt(pi) - z erfc(z)``.

    Parameters
    ----------
    z : float
        Real argument.

    Returns
    -------
    float
        ``i^1 erfc(z) = int_z^inf erfc(t) dt``.

    Raises
    ------
    DomainError
        If ``z`` is not finite.

    Notes
    -----
    Definition DLMF 7.18.2 with n = 1; the closed form follows from the
    recurrence DLMF 7.18.7 with the starting values of 7.18.1 [1]_. For
    z > 25 the difference is formed from the tail (terms m >= 1) of the
    asymptotic series DLMF 7.12.1 to avoid cancellation. ``erfc`` is
    :func:`math.erfc` [2]_.

    References
    ----------
    .. [1] DLMF-7-ERFC. NIST DLMF Eqs. 7.18.1, 7.18.2, 7.18.7 and 7.12.1,
           https://dlmf.nist.gov/7.18
    .. [2] PYTHON-3.12-MATH-ERF. Python 3.12, math.erf and math.erfc,
           https://docs.python.org/3.12/library/math.html#math.erfc
    """
    value = real(z, "z")
    if value > 25.0:
        # 1/sqrt(pi) - z erfcx(z) = (1/sqrt(pi)) sum_{m>=1} (-1)^(m+1) (1/2)_m / z^(2m)
        inverse = 1.0 / (value * value)
        term = 0.5 * inverse
        total = 0.0
        m = 1
        while True:
            total += term
            m += 1
            new_term = -term * (m - 0.5) * inverse
            if abs(new_term) >= abs(term) or abs(new_term) <= 1.0e-17 * abs(total):
                break
            term = new_term
        return math.exp(-value * value) * total / math.sqrt(math.pi)
    return math.exp(-value * value) / math.sqrt(math.pi) - value * math.erfc(value)


# ---------------------------------------------------------------------------
# temperature-dependent properties and the Kirchhoff transform
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class TabulatedProperty:
    """Piecewise-linear positive property of temperature, valid only inside its table.

    Used for thermal conductivity ``k(T)`` (W/(m K)) or volumetric heat
    capacity ``rho c(T)`` (J/(m^3 K)). Values are interpolated linearly
    between nodes and are never extrapolated.

    Attributes
    ----------
    temperatures_k : tuple of float
        Strictly increasing absolute temperatures, at least two.
    values : tuple of float
        Positive property values at those temperatures.
    unit : str
        Unit label of ``values``.

    Raises
    ------
    InputError
        If the two tuples differ in length or have fewer than two entries.
    DomainError
        If a temperature is not positive or not increasing, or a value is not
        positive.

    Notes
    -----
    The integral ``int k dT`` is the Kirchhoff variable of [1]_ (Sec. 3,
    Eq. (5)); for a piecewise-linear table it is exact (trapezoidal per
    segment) and its inverse is a stable quadratic root per segment.
    Lienhard and Lienhard note that k almost always varies with temperature
    and give ``q L = int k dT`` for a plane wall [2]_.

    References
    ----------
    .. [1] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352,
           Sec. 3, Eq. (5),
           https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
    .. [2] AHTT-V6. Lienhard and Lienhard, "A Heat Transfer Textbook", 6th ed.,
           Sec. 2.1, printed p. 51, and Problem 2.2, printed p. 86,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    temperatures_k: tuple[float, ...]
    values: tuple[float, ...]
    unit: str = ""

    def __post_init__(self) -> None:
        if isinstance(self.temperatures_k, str | bytes) or isinstance(self.values, str | bytes):
            raise InputError("Tabulated property tables must be sequences of numbers.",
                             field="values")
        temperatures = tuple(_temperature(t, "temperatures_k") for t in self.temperatures_k)
        values = tuple(_positive(v, "values") for v in self.values)
        if len(temperatures) != len(values) or len(temperatures) < 2:
            raise InputError("A tabulated property needs matching temperature and value "
                             "tuples with at least two entries.", field="values")
        if any(b <= a for a, b in itertools.pairwise(temperatures)):
            raise DomainError("Tabulated temperatures must increase strictly.",
                              field="temperatures_k")
        if not isinstance(self.unit, str):
            raise InputError("unit must be a string.", field="unit")
        object.__setattr__(self, "temperatures_k", temperatures)
        object.__setattr__(self, "values", values)

    @property
    def minimum_temperature_k(self) -> float:
        """Lowest tabulated temperature, K.

        References
        ----------
        .. [1] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352,
               Sec. 3, https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
        """
        return self.temperatures_k[0]

    @property
    def maximum_temperature_k(self) -> float:
        """Highest tabulated temperature, K.

        References
        ----------
        .. [1] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352,
               Sec. 3, https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
        """
        return self.temperatures_k[-1]

    def _cumulative(self) -> tuple[float, ...]:
        return _cumulative_integral(self.temperatures_k, self.values)

    def _segment(self, temperature: float) -> int:
        nodes = self.temperatures_k
        if not nodes[0] <= temperature <= nodes[-1]:
            raise DomainError(
                f"Temperature {temperature:.6g} K is outside the property table "
                f"[{nodes[0]:.6g}, {nodes[-1]:.6g}] K.", code="property_out_of_range",
                field="temperature_k")
        low, high = 0, len(nodes) - 2
        while low < high:
            mid = (low + high + 1) // 2
            if nodes[mid] <= temperature:
                low = mid
            else:
                high = mid - 1
        return low

    def value_at(self, temperature_k: float) -> float:
        """Interpolated property value at ``temperature_k``.

        Parameters
        ----------
        temperature_k : float
            Absolute temperature inside the table.

        Returns
        -------
        float
            Linear interpolation between the bracketing nodes.

        Raises
        ------
        DomainError
            If the temperature is outside the table (code
            ``"property_out_of_range"``).

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 2.1, printed p. 51,
               https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        temperature = float(temperature_k)
        i = self._segment(temperature)
        t0, t1 = self.temperatures_k[i], self.temperatures_k[i + 1]
        v0, v1 = self.values[i], self.values[i + 1]
        return v0 + (v1 - v0) * (temperature - t0) / (t1 - t0)

    def _antiderivative(self, temperature: float) -> float:
        i = self._segment(temperature)
        t0, t1 = self.temperatures_k[i], self.temperatures_k[i + 1]
        v0, v1 = self.values[i], self.values[i + 1]
        u = temperature - t0
        slope = (v1 - v0) / (t1 - t0)
        return self._cumulative()[i] + v0 * u + 0.5 * slope * u * u

    def integral(self, lower_k: float, upper_k: float) -> float:
        """Return ``int_lower^upper p(T) dT`` exactly for the piecewise-linear table.

        Parameters
        ----------
        lower_k, upper_k : float
            Integration limits inside the table, K.

        Returns
        -------
        float
            The integral (negative when ``upper_k < lower_k``).

        Raises
        ------
        DomainError
            If a limit is outside the table.

        References
        ----------
        .. [1] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352,
               Sec. 3, Eq. (5),
               https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
        """
        return self._antiderivative(float(upper_k)) - self._antiderivative(float(lower_k))

    def solve_integral(self, start_k: float, target: float) -> float:
        """Return T with ``int_start^T p dT = target``; +-inf if T leaves the table.

        Parameters
        ----------
        start_k : float
            Lower limit inside the table, K.
        target : float
            Required value of the integral.

        Returns
        -------
        float
            The temperature, or ``math.inf`` / ``-math.inf`` when the target
            lies above / below what the table can reach.

        Raises
        ------
        DomainError
            If ``start_k`` is outside the table.

        References
        ----------
        .. [1] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352,
               Sec. 3, Eqs. (5)-(6),
               https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
        """
        cumulative = self._cumulative()
        wanted = self._antiderivative(float(start_k)) + float(target)
        if wanted < 0.0:
            return -math.inf
        if wanted > cumulative[-1]:
            return math.inf
        i = 0
        while i < len(cumulative) - 2 and cumulative[i + 1] < wanted:
            i += 1
        t0, t1 = self.temperatures_k[i], self.temperatures_k[i + 1]
        v0, v1 = self.values[i], self.values[i + 1]
        slope = (v1 - v0) / (t1 - t0)
        remaining = wanted - cumulative[i]
        # (slope/2) u^2 + v0 u - remaining = 0, numerically stable positive root
        u = 2.0 * remaining / (v0 + math.sqrt(max(v0 * v0 + 2.0 * slope * remaining, 0.0)))
        return min(max(t0 + u, t0), t1)


@lru_cache(maxsize=256)
def _cumulative_integral(temperatures: tuple[float, ...],
                         values: tuple[float, ...]) -> tuple[float, ...]:
    cumulative = [0.0]
    for i in range(len(temperatures) - 1):
        width = temperatures[i + 1] - temperatures[i]
        cumulative.append(cumulative[-1] + 0.5 * width * (values[i] + values[i + 1]))
    return tuple(cumulative)


class _Constant:
    """Constant property with the same interface as :class:`TabulatedProperty`."""

    __slots__ = ("value",)
    minimum_temperature_k = 0.0
    maximum_temperature_k = math.inf

    def __init__(self, value: float) -> None:
        self.value = value

    def value_at(self, temperature_k: float) -> float:
        """Return the constant value."""
        return self.value

    def integral(self, lower_k: float, upper_k: float) -> float:
        """Return ``value * (upper - lower)``."""
        return self.value * (upper_k - lower_k)

    def solve_integral(self, start_k: float, target: float) -> float:
        """Return ``start + target / value``."""
        return start_k + target / self.value


def _property(value: object, label: str):
    if isinstance(value, TabulatedProperty):
        return value
    return _Constant(_positive(value, label))


# ---------------------------------------------------------------------------
# geometry, layers and boundary conditions
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class WallGeometry:
    """One-dimensional wall geometry: plane, cylindrical shell or spherical shell.

    Attributes
    ----------
    kind : {"plane", "cylinder", "sphere"}
        Coordinate system.
    inner_radius_m : float
        Radius of the inner surface (cylinder, sphere); 0 for a solid body
        whose centre is a symmetry point. Must be 0 for a plane wall, whose
        coordinate starts at x = 0 on the inner surface.
    area_m2 : float
        Plane-wall area normal to the heat flow (plane only), m^2.
    length_m : float
        Cylinder length (cylinder only), m; heat rates are per this length.

    Raises
    ------
    InputError
        If ``kind`` is unknown.
    DomainError
        If a dimension is negative, zero where it must be positive, or set for
        a geometry that does not use it.

    Notes
    -----
    The steady and transient models use the plane, cylindrical and spherical
    forms of the conduction equation (2.10) with the Laplacians (2.12),
    (2.13) and (2.14b) of [1]_; heat flows through the areas ``A``,
    ``2 pi r l`` and ``4 pi r^2``.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (2.10)-(2.14), printed
           pp. 55-56, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    kind: str
    inner_radius_m: float = 0.0
    area_m2: float = 1.0
    length_m: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in _SHAPES:
            raise InputError(f"Geometry kind must be one of {_SHAPES}.", field="kind")
        inner = _nonnegative(self.inner_radius_m, "inner_radius_m")
        area = _positive(self.area_m2, "area_m2")
        length = _positive(self.length_m, "length_m")
        if self.kind == "plane" and inner != 0.0:
            raise DomainError("A plane wall has no inner radius.", field="inner_radius_m")
        if self.kind != "plane" and area != 1.0:
            raise DomainError("area_m2 applies to plane walls only.", field="area_m2")
        if self.kind != "cylinder" and length != 1.0:
            raise DomainError("length_m applies to cylinders only.", field="length_m")
        object.__setattr__(self, "inner_radius_m", inner)
        object.__setattr__(self, "area_m2", area)
        object.__setattr__(self, "length_m", length)

    def area(self, position_m: float) -> float:
        """Heat-flow area at ``position_m`` (x or r), m^2.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (2.16), (2.21), printed
               pp. 62-66, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        if self.kind == "plane":
            return self.area_m2
        if self.kind == "cylinder":
            return 2.0 * math.pi * position_m * self.length_m
        return 4.0 * math.pi * position_m * position_m

    def volume(self, start_m: float, end_m: float) -> float:
        """Volume between positions ``start_m`` < ``end_m``, m^3.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.10) control volume,
               printed p. 55, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        if self.kind == "plane":
            return self.area_m2 * (end_m - start_m)
        if self.kind == "cylinder":
            return math.pi * self.length_m * (end_m - start_m) * (end_m + start_m)
        return 4.0 / 3.0 * math.pi * (end_m**3 - start_m**3)

    def shape_integral(self, start_m: float, end_m: float) -> float:
        """``int dr / A(r)`` between two positions, 1/m (conduction resistance times k).

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (2.16) and (2.22), and
               Problem 2.21, printed pp. 62, 66, 90,
               https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        if self.kind == "plane":
            return (end_m - start_m) / self.area_m2
        if start_m == 0.0:
            return math.inf
        if self.kind == "cylinder":
            return math.log(end_m / start_m) / (2.0 * math.pi * self.length_m)
        return (1.0 / start_m - 1.0 / end_m) / (4.0 * math.pi)


@dataclass(frozen=True, slots=True)
class WallLayer:
    """One homogeneous layer of a multilayer wall.

    Attributes
    ----------
    thickness_m : float
        Layer thickness, m.
    thermal_conductivity : float or TabulatedProperty
        ``k`` in W/(m K), constant or tabulated in temperature.
    volumetric_heat_capacity : float or TabulatedProperty or None
        ``rho c`` in J/(m^3 K); required by the transient solver only.
    heat_generation_w_m3 : float
        Uniform volumetric heat release ``q'''`` in W/m^3 (negative for a
        uniform sink).
    name : str
        Optional label.

    Raises
    ------
    DomainError
        If the thickness or a property is not positive and finite.
    InputError
        If ``name`` is not a string.

    Notes
    -----
    The volumetric source enters the conduction equation (2.10) of [1]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.10), printed p. 55,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    thickness_m: float
    thermal_conductivity: float | TabulatedProperty
    volumetric_heat_capacity: float | TabulatedProperty | None = None
    heat_generation_w_m3: float = 0.0
    name: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "thickness_m", _positive(self.thickness_m, "thickness_m"))
        if not isinstance(self.thermal_conductivity, TabulatedProperty):
            object.__setattr__(self, "thermal_conductivity",
                               _positive(self.thermal_conductivity, "thermal_conductivity"))
        capacity = self.volumetric_heat_capacity
        if capacity is not None and not isinstance(capacity, TabulatedProperty):
            object.__setattr__(self, "volumetric_heat_capacity",
                               _positive(capacity, "volumetric_heat_capacity"))
        object.__setattr__(self, "heat_generation_w_m3",
                           _finite(self.heat_generation_w_m3, "heat_generation_w_m3"))
        if not isinstance(self.name, str):
            raise InputError("Layer name must be a string.", field="name")


@dataclass(frozen=True, slots=True)
class SurfaceTemperature:
    """Prescribed surface temperature (boundary condition of the first kind).

    Attributes
    ----------
    temperature_k : float
        Surface temperature, K.

    Raises
    ------
    DomainError
        If the temperature is not above 0 K.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Example 2.4, printed pp. 64-66,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    temperature_k: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "temperature_k", _temperature(self.temperature_k,
                                                               "temperature_k"))


@dataclass(frozen=True, slots=True)
class SurfaceHeatFlux:
    """Prescribed heat flux into the wall through this surface (second kind).

    Attributes
    ----------
    heat_flux_w_m2 : float
        Flux entering the wall, W/m^2; zero for an adiabatic or symmetry
        surface, negative for heat leaving.

    Raises
    ------
    DomainError
        If the flux is not finite.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 5.6, second-kind b.c. before
           Eq. (5.55), printed p. 229,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    heat_flux_w_m2: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "heat_flux_w_m2", _finite(self.heat_flux_w_m2,
                                                           "heat_flux_w_m2"))


@dataclass(frozen=True, slots=True)
class SurfaceConvection:
    """Convective (third-kind) surface, optionally with gray radiation to large surroundings.

    The heat flux into the wall is
    ``h (T_fluid - T_s) + eps sigma (T_surr^4 - T_s^4)``.

    Attributes
    ----------
    heat_transfer_coefficient_w_m2_k : float
        Film coefficient ``h`` (may be 0 when radiation alone acts).
    fluid_temperature_k : float
        Fluid temperature ``T_inf``, K.
    emissivity : float
        Surface emittance in [0, 1]; 0 disables radiation.
    surroundings_temperature_k : float or None
        Radiating-surroundings temperature; defaults to the fluid temperature.

    Raises
    ------
    DomainError
        If ``h`` is negative, the emittance is outside [0, 1], a temperature
        is not above 0 K, or both ``h`` and the emittance are zero.

    Notes
    -----
    Convection follows the energy balance of Example 2.5 of [1]_; radiation
    uses eqn. (1.34) with the small-object transfer factor ``F = eps`` of
    eqn. (1.35), which assumes the surface sees a much larger isothermal
    environment. That assumption does not hold for the concave inner surface
    of a closed tube or sphere. ``sigma`` is the CODATA value [2]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Example 2.5 and Eqs. (1.34)-(1.35),
           printed pp. 33 and 67,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] NIST-CODATA-2022-STEFAN-BOLTZMANN. NIST, Stefan-Boltzmann constant,
           https://physics.nist.gov/cgi-bin/cuu/Value?sigma
    """

    heat_transfer_coefficient_w_m2_k: float
    fluid_temperature_k: float
    emissivity: float = 0.0
    surroundings_temperature_k: float | None = None

    def __post_init__(self) -> None:
        h = _nonnegative(self.heat_transfer_coefficient_w_m2_k,
                         "heat_transfer_coefficient_w_m2_k")
        emissivity = _emissivity(self.emissivity, "emissivity")
        if h == 0.0 and emissivity == 0.0:
            raise DomainError("A convective surface needs h > 0 or a nonzero emittance.",
                              field="heat_transfer_coefficient_w_m2_k")
        fluid = _temperature(self.fluid_temperature_k, "fluid_temperature_k")
        surroundings = (fluid if self.surroundings_temperature_k is None
                        else _temperature(self.surroundings_temperature_k,
                                          "surroundings_temperature_k"))
        object.__setattr__(self, "heat_transfer_coefficient_w_m2_k", h)
        object.__setattr__(self, "emissivity", emissivity)
        object.__setattr__(self, "fluid_temperature_k", fluid)
        object.__setattr__(self, "surroundings_temperature_k", surroundings)

    def flux_into_wall(self, surface_temperature_k: float) -> float:
        """Heat flux entering the wall at surface temperature ``T_s``, W/m^2.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.34)-(1.35), printed
               p. 33, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        t = surface_temperature_k
        flux = self.heat_transfer_coefficient_w_m2_k * (self.fluid_temperature_k - t)
        if self.emissivity:
            flux += self.emissivity * STEFAN_BOLTZMANN_W_M2_K4 * (
                self.surroundings_temperature_k**4 - t**4)
        return flux

    def flux_derivative(self, surface_temperature_k: float) -> float:
        """Derivative of :meth:`flux_into_wall` with respect to ``T_s`` (negative).

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.34)-(1.35), printed
               p. 33, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        return -(self.heat_transfer_coefficient_w_m2_k
                 + 4.0 * self.emissivity * STEFAN_BOLTZMANN_W_M2_K4 * surface_temperature_k**3)


_BOUNDARY_TYPES = (SurfaceTemperature, SurfaceHeatFlux, SurfaceConvection)


def _check_boundary(value: object, label: str):
    if not isinstance(value, _BOUNDARY_TYPES):
        raise InputError(f"{label} must be a SurfaceTemperature, SurfaceHeatFlux or "
                         "SurfaceConvection.", field=label)
    return value


def _check_layers(layers: object) -> tuple[WallLayer, ...]:
    if isinstance(layers, WallLayer):
        layers = (layers,)
    if not isinstance(layers, Sequence) or isinstance(layers, str) or not layers:
        raise InputError("At least one WallLayer is required.", field="layers")
    owned = tuple(layers)
    if not all(isinstance(layer, WallLayer) for layer in owned):
        raise InputError("Every layer must be a WallLayer.", field="layers")
    return owned


def _check_contacts(contacts: object, count: int) -> tuple[float, ...]:
    """Return contact conductances (inf = perfect contact) for the count-1 interfaces."""
    if contacts is None:
        return (math.inf,) * (count - 1)
    if not isinstance(contacts, Sequence) or isinstance(contacts, str) or len(
            contacts) != count - 1:
        raise InputError("Give one contact conductance (or None) per layer interface.",
                         field="contact_conductances_w_m2_k")
    owned = []
    for value in contacts:
        if value is None:
            owned.append(math.inf)
        else:
            owned.append(_positive(value, "contact_conductances_w_m2_k"))
    return tuple(owned)


def _layer_faces(geometry: WallGeometry, layers: tuple[WallLayer, ...]) -> tuple[float, ...]:
    faces = [geometry.inner_radius_m]
    for layer in layers:
        faces.append(faces[-1] + layer.thickness_m)
    return tuple(faces)


# ---------------------------------------------------------------------------
# thermal resistances
# ---------------------------------------------------------------------------

def plane_wall_resistance(thickness_m: float, conductivity_w_m_k: float,
                          area_m2: float) -> float:
    """Conduction resistance of a plane slab, ``R = L / (k A)``.

    Parameters
    ----------
    thickness_m, conductivity_w_m_k, area_m2 : float
        Slab thickness, conductivity and area (all positive).

    Returns
    -------
    float
        Resistance, K/W.

    Raises
    ------
    DomainError
        If an input is not positive and finite.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.16), printed p. 62,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    return (_positive(thickness_m, "thickness_m")
            / (_positive(conductivity_w_m_k, "conductivity_w_m_k")
               * _positive(area_m2, "area_m2")))


def cylindrical_wall_resistance(inner_radius_m: float, outer_radius_m: float,
                                conductivity_w_m_k: float, length_m: float = 1.0) -> float:
    """Conduction resistance of a cylindrical shell, ``ln(r_o/r_i) / (2 pi l k)``.

    Parameters
    ----------
    inner_radius_m, outer_radius_m : float
        Radii with ``0 < r_i < r_o``, m.
    conductivity_w_m_k : float
        Conductivity, W/(m K).
    length_m : float, optional
        Cylinder length, m.

    Returns
    -------
    float
        Resistance, K/W.

    Raises
    ------
    DomainError
        If an input is not positive or ``r_o <= r_i``.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.22), printed p. 66,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    inner = _positive(inner_radius_m, "inner_radius_m")
    outer = _positive(outer_radius_m, "outer_radius_m")
    if outer <= inner:
        raise DomainError("Outer radius must exceed inner radius.", field="outer_radius_m")
    return math.log(outer / inner) / (2.0 * math.pi * _positive(length_m, "length_m")
                                      * _positive(conductivity_w_m_k, "conductivity_w_m_k"))


def spherical_wall_resistance(inner_radius_m: float, outer_radius_m: float,
                              conductivity_w_m_k: float) -> float:
    """Conduction resistance of a spherical shell, ``(1/r_i - 1/r_o) / (4 pi k)``.

    Parameters
    ----------
    inner_radius_m : float
        Inner radius, m (positive).
    outer_radius_m : float
        Outer radius, m; ``math.inf`` gives a sphere in an infinite medium.
    conductivity_w_m_k : float
        Conductivity, W/(m K).

    Returns
    -------
    float
        Resistance, K/W.

    Raises
    ------
    DomainError
        If an input is not positive or ``r_o <= r_i``.

    Notes
    -----
    Obtained by integrating the spherical conduction equation (2.14b) of [1]_
    with the steps of its Example 2.4 (the book poses the derivation as
    Problem 2.21). It reproduces the printed answer Q = 16.9 W of Problem 2.15
    (an isothermal sphere in an infinite medium).

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.14b), printed p. 56;
           Problems 2.15 and 2.21, printed pp. 89-90,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    inner = _positive(inner_radius_m, "inner_radius_m")
    if isinstance(outer_radius_m, float) and outer_radius_m == math.inf:
        outer = math.inf
    else:
        outer = _positive(outer_radius_m, "outer_radius_m")
    if outer <= inner:
        raise DomainError("Outer radius must exceed inner radius.", field="outer_radius_m")
    return (1.0 / inner - 1.0 / outer) / (
        4.0 * math.pi * _positive(conductivity_w_m_k, "conductivity_w_m_k"))


def convection_resistance(heat_transfer_coefficient_w_m2_k: float, area_m2: float) -> float:
    """Convective film resistance, ``R = 1 / (h A)``.

    Parameters
    ----------
    heat_transfer_coefficient_w_m2_k, area_m2 : float
        Film coefficient and wetted area (positive).

    Returns
    -------
    float
        Resistance, K/W.

    Raises
    ------
    DomainError
        If an input is not positive and finite.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.26), printed p. 69,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    return 1.0 / (_positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
                  * _positive(area_m2, "area_m2"))


def contact_resistance(interfacial_conductance_w_m2_k: float, area_m2: float) -> float:
    """Interfacial (contact) resistance, ``R = 1 / (h_c A)``.

    Parameters
    ----------
    interfacial_conductance_w_m2_k, area_m2 : float
        Interfacial conductance ``h_c`` and contact area (positive).

    Returns
    -------
    float
        Resistance, K/W.

    Raises
    ------
    DomainError
        If an input is not positive and finite.

    Notes
    -----
    Lienhard and Lienhard tabulate typical ``h_c`` (Table 2.1) for moderate
    contact pressures; the value depends on finish, materials, pressure and
    interstitial medium and must be supplied by the caller [1]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 2.3 "Contact resistance",
           Table 2.1, printed p. 76,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    return 1.0 / (_positive(interfacial_conductance_w_m2_k, "interfacial_conductance_w_m2_k")
                  * _positive(area_m2, "area_m2"))


def radiation_heat_transfer_coefficient(mean_temperature_k: float, emissivity: float) -> float:
    """Linearized radiation coefficient of a small gray body, ``h_rad = 4 eps sigma T_m^3``.

    Parameters
    ----------
    mean_temperature_k : float
        Mean of the body and surroundings temperatures, K.
    emissivity : float
        Emittance of the small body, in [0, 1].

    Returns
    -------
    float
        ``h_rad``, W/(m^2 K); valid for ``(dT/T_m)^2 / 4 << 1``.

    Raises
    ------
    DomainError
        If the temperature is not above 0 K or the emittance is outside [0, 1].

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (2.29) and (2.31), printed
           p. 72, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] NIST-CODATA-2022-STEFAN-BOLTZMANN. NIST, Stefan-Boltzmann constant,
           https://physics.nist.gov/cgi-bin/cuu/Value?sigma
    """
    temperature = _temperature(mean_temperature_k, "mean_temperature_k")
    return 4.0 * STEFAN_BOLTZMANN_W_M2_K4 * temperature**3 * _emissivity(emissivity, "emissivity")


def critical_insulation_radius(conductivity_w_m_k: float,
                               heat_transfer_coefficient_w_m2_k: float) -> float:
    """Critical radius of insulation on a cylinder, ``r_crit = k / h`` (Bi = 1).

    Parameters
    ----------
    conductivity_w_m_k : float
        Insulation conductivity.
    heat_transfer_coefficient_w_m2_k : float
        Outside film coefficient.

    Returns
    -------
    float
        Radius at which the insulated-cylinder heat loss is maximum, m.

    Raises
    ------
    DomainError
        If an input is not positive and finite.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Example 2.6, Eq. (2.27), printed
           p. 70, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    return (_positive(conductivity_w_m_k, "conductivity_w_m_k")
            / _positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k"))


# ---------------------------------------------------------------------------
# steady multilayer wall
# ---------------------------------------------------------------------------

_STEADY_ASSUMPTIONS = (
    "steady one-dimensional conduction in plane, cylindrical or spherical coordinates",
    "homogeneous layers with uniform volumetric heat generation per layer",
    "temperature-dependent conductivity handled exactly through the Kirchhoff integral",
    "interfaces with a prescribed uniform contact conductance or perfect contact",
    "radiation with transfer factor equal to the surface emittance (large surroundings)",
)


def _kirchhoff_drop(geometry: WallGeometry, start: float, end: float, heat_rate: float,
                    generation: float) -> float:
    """Return ``phi(start) - phi(end) = int Q(r)/A(r) dr`` across a layer slice.

    ``Q(r) = Q(start) + q''' V(start, r)`` is the heat rate in +r; ``phi`` is
    the Kirchhoff variable ``int k dT``.
    """
    if geometry.kind == "plane":
        length = end - start
        return heat_rate * length / geometry.area_m2 + 0.5 * generation * length * length
    if geometry.kind == "cylinder":
        if start == 0.0:
            return 0.25 * generation * end * end
        log_ratio = math.log1p((end - start) / start)
        return (heat_rate * log_ratio / (2.0 * math.pi * geometry.length_m)
                + 0.5 * generation * (0.5 * (end - start) * (end + start)
                                      - start * start * log_ratio))
    if start == 0.0:
        return generation * end * end / 6.0
    return (heat_rate * (1.0 / start - 1.0 / end) / (4.0 * math.pi)
            + generation / 3.0 * (0.5 * (end - start) * (end + start)
                                  - start * start * (end - start) / end))


def _stationary_position(geometry: WallGeometry, start: float, heat_rate: float,
                         generation: float) -> float | None:
    """Position where Q(r) = 0 inside a generating layer (local temperature extremum)."""
    if generation == 0.0:
        return None
    if geometry.kind == "plane":
        target = start - heat_rate / (generation * geometry.area_m2)
        return target
    if geometry.kind == "cylinder":
        square = start * start - heat_rate / (math.pi * geometry.length_m * generation)
        return math.sqrt(square) if square >= 0.0 else None
    cube = start**3 - 3.0 * heat_rate / (4.0 * math.pi * generation)
    return cube ** (1.0 / 3.0) if cube >= 0.0 else None


def _march(geometry, faces, conductivities, generations, contacts, inner_temperature,
           inner_heat_rate, *, samples: int = 0):
    """March temperature and heat rate from the inner to the outer surface.

    Returns ``(T_out, Q_out, layer_in, layer_out, profile)``; ``T_out`` is
    ``+-inf`` when a temperature leaves a conductivity table or falls to 0 K.
    """
    temperature = inner_temperature
    heat_rate = inner_heat_rate
    layer_in, layer_out, profile = [], [], []
    count = len(conductivities)
    for i in range(count):
        start, end = faces[i], faces[i + 1]
        prop, generation = conductivities[i], generations[i]
        if temperature < prop.minimum_temperature_k:
            return -math.inf, math.nan, layer_in, layer_out, profile
        if temperature > prop.maximum_temperature_k:
            return math.inf, math.nan, layer_in, layer_out, profile
        layer_in.append(temperature)
        if samples:
            points = [start + (end - start) * j / samples for j in range(samples + 1)]
            stationary = _stationary_position(geometry, start, heat_rate, generation)
            if stationary is not None and start < stationary < end:
                points.append(stationary)
                points.sort()
            for position in points:
                drop = _kirchhoff_drop(geometry, start, position, heat_rate, generation)
                profile.append((position, prop.solve_integral(temperature, -drop)))
        drop = _kirchhoff_drop(geometry, start, end, heat_rate, generation)
        temperature = prop.solve_integral(temperature, -drop)
        if not math.isfinite(temperature) or temperature <= 0.0:
            return (math.inf if temperature > 0.0 else -math.inf), math.nan, layer_in, \
                layer_out, profile
        heat_rate = heat_rate + generation * geometry.volume(start, end)
        layer_out.append(temperature)
        if i < count - 1 and contacts[i] != math.inf:
            temperature -= heat_rate / (contacts[i] * geometry.area(end))
            if temperature <= 0.0:
                return -math.inf, math.nan, layer_in, layer_out, profile
    return temperature, heat_rate, layer_in, layer_out, profile


def _bisect_sign(function: Callable[[float], float], low: float, high: float,
                 *, iterations: int = 400) -> tuple[float, int]:
    """Bisection on a sign change of ``function`` over ``[low, high]`` (infinities allowed)."""
    f_low = function(low)
    if f_low == 0.0:
        return low, 0
    count = 0
    for count in range(1, iterations + 1):
        middle = 0.5 * (low + high)
        if middle in (low, high):
            return middle, count
        f_middle = function(middle)
        if f_middle == 0.0:
            return middle, count
        if (f_middle > 0.0) == (f_low > 0.0):
            low, f_low = middle, f_middle
        else:
            high = middle
        if high - low <= 4.0e-16 * max(abs(low), abs(high), 1.0e-300):
            return 0.5 * (low + high), count
    raise ConvergenceError("Steady-wall bisection did not converge.",
                           details={"low": low, "high": high})


def _bracket(function: Callable[[float], float], low: float, high: float, *,
             expand_low: bool, expand_high: bool, label: str) -> tuple[float, float]:
    """Expand ``[low, high]`` geometrically until ``function`` changes sign."""
    f_low, f_high = function(low), function(high)
    for _ in range(400):
        if (f_low > 0.0) != (f_high > 0.0) or f_low == 0.0 or f_high == 0.0:
            return low, high
        if expand_high:
            high = high * 4.0 if high > 0.0 else high + 4.0 * (high - low)
            f_high = function(high)
        if expand_low:
            low = low * 4.0 if low < 0.0 else low - 4.0 * (high - low)
            f_low = function(low)
        if not (math.isfinite(low) and math.isfinite(high)) or max(abs(low), abs(high)) > 1e200:
            break
    raise DomainError(f"No steady solution was bracketed for the {label}; the boundary "
                      "conditions may be inconsistent or leave a property table.",
                      code="steady_solution_not_bracketed", field=label)


def _outer_residual(boundary, area: float, temperature: float, heat_rate: float) -> float:
    """Outer-surface energy residual; increasing in T_out, sentinel-aware."""
    if isinstance(boundary, SurfaceTemperature):
        return temperature - boundary.temperature_k
    if isinstance(boundary, SurfaceHeatFlux):
        return heat_rate + boundary.heat_flux_w_m2 * area
    if temperature == math.inf:
        return -math.inf
    if temperature == -math.inf:
        return math.inf
    return heat_rate + area * boundary.flux_into_wall(temperature)


def _solve_surface_for_flux(boundary: SurfaceConvection, flux_into: float) -> float:
    """Surface temperature at which a convective surface passes ``flux_into`` (W/m^2)."""
    def residual(temperature: float) -> float:
        return boundary.flux_into_wall(temperature) - flux_into

    low, high = 1.0e-9, max(boundary.fluid_temperature_k,
                            boundary.surroundings_temperature_k) * 2.0
    low, high = _bracket(residual, low, high, expand_low=False, expand_high=True,
                         label="surface temperature")
    return _bisect_sign(residual, low, high)[0]


@dataclass(frozen=True, slots=True)
class SteadyWallResult:
    """Steady temperature field and heat rates of a multilayer wall.

    Attributes
    ----------
    geometry : WallGeometry
        Geometry used.
    layer_boundaries_m : tuple of float
        Positions (x or r) of the inner surface, the interfaces and the outer
        surface.
    layer_inner_temperatures_k, layer_outer_temperatures_k : tuple of float
        Temperature at each layer's inner and outer face (differences across
        an interface are contact-resistance jumps).
    inner_heat_rate_w : float
        Heat rate entering the wall at the inner surface, flowing toward the
        outer surface (W; per ``length_m`` for a cylinder, per ``area_m2``
        for a plane wall).
    outer_heat_rate_w : float
        Heat rate leaving the wall at the outer surface, W.
    inner_heat_flux_w_m2 : float or None
        Inner-surface flux; None for a solid body (r = 0).
    outer_heat_flux_w_m2 : float
        Outer-surface flux leaving the wall.
    generated_heat_rate_w : float
        Total volumetric generation, W.
    maximum_temperature_k, maximum_temperature_position_m : float
        Largest temperature, from the profile samples and the analytic
        stationary points of generating layers.
    profile_positions_m, profile_temperatures_k : tuple of float
        Sampled temperature profile (exact at the sample points).
    energy_balance_residual_w : float
        ``Q_out - Q_in - generation`` (zero up to rounding).
    boundary_residual_w : float
        Outer boundary-condition residual at the accepted solution, W.
    iterations : int
        Bisection iterations.
    warnings : tuple of WarningMessage
        Model warnings.
    assumptions : tuple of str
        Declared modelling assumptions.
    physical_model_qualified : bool
        Always False.
    model : str
        Model identifier.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Secs. 2.2-2.4, Eqs. (2.10)-(2.26),
           printed pp. 55-80, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352, Sec. 3,
           Eq. (5), https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
    """

    geometry: WallGeometry
    layer_boundaries_m: tuple[float, ...]
    layer_inner_temperatures_k: tuple[float, ...]
    layer_outer_temperatures_k: tuple[float, ...]
    inner_heat_rate_w: float
    outer_heat_rate_w: float
    inner_heat_flux_w_m2: float | None
    outer_heat_flux_w_m2: float
    generated_heat_rate_w: float
    maximum_temperature_k: float
    maximum_temperature_position_m: float
    profile_positions_m: tuple[float, ...]
    profile_temperatures_k: tuple[float, ...]
    energy_balance_residual_w: float
    boundary_residual_w: float
    iterations: int
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...] = _STEADY_ASSUMPTIONS
    physical_model_qualified: bool = False
    model: str = CONDUCTION_MODEL

    @property
    def inner_surface_temperature_k(self) -> float:
        """Temperature of the inner surface (or centre of a solid body), K.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Secs. 2.2-2.3, printed pp. 58-77,
               https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        return self.layer_inner_temperatures_k[0]

    @property
    def outer_surface_temperature_k(self) -> float:
        """Temperature of the outer surface, K.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Secs. 2.2-2.3, printed pp. 58-77,
               https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        return self.layer_outer_temperatures_k[-1]


def solve_steady_wall(
    layers: Sequence[WallLayer] | WallLayer,
    geometry: WallGeometry,
    outer_boundary: SurfaceTemperature | SurfaceHeatFlux | SurfaceConvection,
    inner_boundary: SurfaceTemperature | SurfaceHeatFlux | SurfaceConvection | None = None,
    *,
    contact_conductances_w_m2_k: Sequence[float | None] | None = None,
    profile_points_per_layer: int = 20,
) -> SteadyWallResult:
    """Solve steady 1-D conduction through a plane, cylindrical or spherical multilayer wall.

    Parameters
    ----------
    layers : sequence of WallLayer or WallLayer
        Layers from the inner surface outward.
    geometry : WallGeometry
        Coordinate system and inner radius.
    outer_boundary : SurfaceTemperature, SurfaceHeatFlux or SurfaceConvection
        Condition at the outer surface.
    inner_boundary : SurfaceTemperature, SurfaceHeatFlux, SurfaceConvection or None
        Condition at the inner surface; None means adiabatic (and is required
        for a solid cylinder or sphere, ``inner_radius_m = 0``).
    contact_conductances_w_m2_k : sequence of float or None, optional
        Interfacial conductance ``h_c`` for each of the ``n - 1`` interfaces;
        None (or omitted) means perfect contact.
    profile_points_per_layer : int, optional
        Number of profile intervals per layer.

    Returns
    -------
    SteadyWallResult
        Interface temperatures, heat rates, profile and diagnostics.

    Raises
    ------
    InputError
        If the boundary conditions are both fluxes (no unique solution), a
        solid body has a non-adiabatic centre, or an argument is malformed.
    DomainError
        If no steady solution is bracketed, for example when it would leave a
        conductivity table or require a temperature at or below 0 K.
    ConvergenceError
        If the bisection does not converge.

    Notes
    -----
    In each layer the Kirchhoff variable ``phi = int k dT`` of [2]_ obeys
    ``div grad phi + q''' = 0``, so with the heat rate ``Q(r) = Q(r_a) +
    q''' V(r_a, r)`` (energy balance, eqn. (2.10) of [1]_) and Fourier's law
    ``Q = -A(r) d phi/dr`` the drop across a layer is
    ``phi(r_a) - phi(r_b) = int Q/A dr``, integrated in closed form for
    ``A = A``, ``2 pi r l`` and ``4 pi r^2``. For constant ``k`` and no
    generation this reduces to the resistances (2.16), (2.22) and the shell
    resistance of Problem 2.21; with generation it reproduces eqn. (2.15)
    and the printed answer of Problem 2.20. For a plane wall it is the
    ``q L = int k dT`` of Problem 2.2. A contact conductance adds the jump
    ``Q / (h_c A)`` (Sec. 2.3). The unknown surface temperature or heat rate is
    found by bisection on the outer boundary residual, which is monotonic in
    the unknown because ``k > 0``.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (2.10), (2.14b), (2.15),
           (2.16), (2.22), (2.25), Sec. 2.3 and Problems 2.2, 2.20, 2.21,
           printed pp. 55-90,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352, Sec. 3,
           Eq. (5), https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
    """
    owned = _check_layers(layers)
    if not isinstance(geometry, WallGeometry):
        raise InputError("geometry must be a WallGeometry.", field="geometry")
    outer = _check_boundary(outer_boundary, "outer_boundary")
    inner = SurfaceHeatFlux(0.0) if inner_boundary is None else _check_boundary(
        inner_boundary, "inner_boundary")
    contacts = _check_contacts(contact_conductances_w_m2_k, len(owned))
    if isinstance(profile_points_per_layer, bool) or not isinstance(
            profile_points_per_layer, int) or profile_points_per_layer < 1:
        raise InputError("profile_points_per_layer must be a positive integer.",
                         field="profile_points_per_layer")
    faces = _layer_faces(geometry, owned)
    solid = geometry.kind != "plane" and geometry.inner_radius_m == 0.0
    if solid and not (isinstance(inner, SurfaceHeatFlux) and inner.heat_flux_w_m2 == 0.0):
        raise InputError("A solid cylinder or sphere has a symmetric centre; leave "
                         "inner_boundary as None.", field="inner_boundary")
    if isinstance(inner, SurfaceHeatFlux) and isinstance(outer, SurfaceHeatFlux):
        raise InputError("Two flux boundaries do not fix the temperature level.",
                         code="steady_wall_underdetermined", field="outer_boundary")
    conductivities = tuple(_property(layer.thermal_conductivity, "thermal_conductivity")
                           for layer in owned)
    generations = tuple(layer.heat_generation_w_m3 for layer in owned)
    generated = sum(g * geometry.volume(faces[i], faces[i + 1])
                    for i, g in enumerate(generations))
    inner_area = geometry.area(faces[0])
    outer_area = geometry.area(faces[-1])

    temperatures = [b.temperature_k for b in (inner, outer)
                    if isinstance(b, SurfaceTemperature)]
    temperatures += [t for b in (inner, outer) if isinstance(b, SurfaceConvection)
                     for t in (b.fluid_temperature_k, b.surroundings_temperature_k)]
    reference = max(temperatures) if temperatures else 300.0
    lower_limit = max(1.0e-9, getattr(conductivities[0], "minimum_temperature_k", 0.0))
    upper_cap = getattr(conductivities[0], "maximum_temperature_k", math.inf)

    def state(unknown: float):
        if isinstance(inner, SurfaceTemperature):
            return inner.temperature_k, unknown
        if isinstance(inner, SurfaceHeatFlux):
            return unknown, inner.heat_flux_w_m2 * inner_area
        return unknown, inner.flux_into_wall(unknown) * inner_area

    def residual(unknown: float) -> float:
        t_in, q_in = state(unknown)
        t_out, q_out, *_ = _march(geometry, faces, conductivities, generations, contacts,
                                  t_in, q_in)
        value = _outer_residual(outer, outer_area, t_out, q_out)
        return -value if isinstance(inner, SurfaceTemperature) else value

    iterations = 0
    if isinstance(outer, SurfaceHeatFlux):
        required = -outer.heat_flux_w_m2 * outer_area - generated  # Q_in
        if isinstance(inner, SurfaceTemperature):
            unknown = required
        else:
            unknown = _solve_surface_for_flux(inner, required / inner_area)
    elif isinstance(inner, SurfaceTemperature):
        low, high = _bracket(residual, -1.0, 1.0, expand_low=True, expand_high=True,
                             label="inner heat rate")
        unknown, iterations = _bisect_sign(residual, low, high)
    else:
        high = upper_cap if math.isfinite(upper_cap) else max(2.0 * reference,
                                                              2.0 * lower_limit)
        low, high = _bracket(residual, lower_limit, high, expand_low=False,
                             expand_high=upper_cap == math.inf, label="inner temperature")
        unknown, iterations = _bisect_sign(residual, low, high)

    t_in, q_in = state(unknown)
    t_out, q_out, layer_in, layer_out, profile = _march(
        geometry, faces, conductivities, generations, contacts, t_in, q_in,
        samples=profile_points_per_layer)
    if not math.isfinite(t_out):
        raise DomainError("The steady solution leaves a conductivity table or falls to 0 K.",
                          code="steady_solution_out_of_range", field="layers")
    hottest = max(profile, key=lambda item: item[1])
    warnings = [WarningMessage(
        "conduction_model_unqualified",
        "One-dimensional steady conduction with uniform layer properties and generation; "
        "multidimensional effects are not modelled.", "info")]
    if (geometry.kind != "plane" and isinstance(inner, SurfaceConvection)
            and inner.emissivity > 0.0):
        warnings.append(WarningMessage(
            "concave_surface_radiation",
            "Radiation on a concave inner surface uses F = emittance, which assumes large "
            "isothermal surroundings (AHTT eqn. 1.35).", "warning", "inner_boundary"))
    return SteadyWallResult(
        geometry=geometry,
        layer_boundaries_m=faces,
        layer_inner_temperatures_k=tuple(layer_in),
        layer_outer_temperatures_k=tuple(layer_out),
        inner_heat_rate_w=q_in,
        outer_heat_rate_w=q_out,
        inner_heat_flux_w_m2=None if inner_area == 0.0 else q_in / inner_area,
        outer_heat_flux_w_m2=q_out / outer_area,
        generated_heat_rate_w=generated,
        maximum_temperature_k=hottest[1],
        maximum_temperature_position_m=hottest[0],
        profile_positions_m=tuple(p for p, _ in profile),
        profile_temperatures_k=tuple(t for _, t in profile),
        energy_balance_residual_w=q_out - q_in - generated,
        boundary_residual_w=_outer_residual(outer, outer_area, t_out, q_out),
        iterations=iterations,
        warnings=tuple(warnings),
    )


# ---------------------------------------------------------------------------
# lumped capacitance
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class LumpedCapacitanceResult:
    """Lumped-capacity transient of a body cooled or heated through a film coefficient.

    Attributes
    ----------
    times_s : tuple of float
        Requested times, s.
    temperatures_k : tuple of float
        Body temperature at those times, K.
    time_constant_s : float
        ``T = rho c V / (h A)``, s.
    characteristic_length_m : float
        ``V / A``, m.
    biot_number : float or None
        ``h (V/A) / k`` when a conductivity was given.
    heat_released_j : tuple of float
        Energy given up by the body since t = 0, ``rho c V (T_i - T(t))``, J
        (negative while the body is heated).
    lumped_assumption_satisfied : bool or None
        ``Bi <= biot_limit`` when the Biot number is known.
    warnings : tuple of WarningMessage
        Biot-number and model warnings.
    model : str
        Model identifier.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.19)-(1.23) and (5.7),
           printed pp. 21-24 and 196,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    times_s: tuple[float, ...]
    temperatures_k: tuple[float, ...]
    time_constant_s: float
    characteristic_length_m: float
    biot_number: float | None
    heat_released_j: tuple[float, ...]
    lumped_assumption_satisfied: bool | None
    warnings: tuple[WarningMessage, ...]
    model: str = "lumped_capacitance_v1"


def lumped_capacitance_response(
    *, initial_temperature_k: float, ambient_temperature_k: float,
    heat_transfer_coefficient_w_m2_k: float, surface_area_m2: float, volume_m3: float,
    density_kg_m3: float, specific_heat_j_kg_k: float, times_s: Sequence[float],
    thermal_conductivity_w_m_k: float | None = None, biot_limit: float = 0.1,
) -> LumpedCapacitanceResult:
    """Temperature history of a lumped body, ``(T - T_inf)/(T_i - T_inf) = exp(-t/T)``.

    Parameters
    ----------
    initial_temperature_k, ambient_temperature_k : float
        ``T_i`` and ``T_inf``, K.
    heat_transfer_coefficient_w_m2_k : float
        Film coefficient ``h`` (an effective value may include radiation).
    surface_area_m2, volume_m3 : float
        Cooled surface area ``A`` and volume ``V``.
    density_kg_m3, specific_heat_j_kg_k : float
        ``rho`` and ``c``.
    times_s : sequence of float
        Non-negative times.
    thermal_conductivity_w_m_k : float, optional
        Body conductivity ``k_b``; enables the Biot-number check.
    biot_limit : float, optional
        Largest Biot number accepted as lumped (default 0.1).

    Returns
    -------
    LumpedCapacitanceResult
        Temperatures, time constant, Biot number and warnings.

    Raises
    ------
    DomainError
        If an input is not finite and positive where required.

    Notes
    -----
    Solution (1.22) of the lumped energy balance (1.19)-(1.20) with the time
    constant (1.23) [1]_. The text states the lumped solution is normally
    accurate to about 3 % for ``Bi <~ 0.1`` (printed p. 24); the Biot
    number here uses the volume-to-area length ``V/A`` of eqn. (5.7). The
    book's Example 1.4 instead takes ``L = D/2`` for its bead, which gives a
    Biot number three times larger than ``V/A`` for a sphere.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.19)-(1.23), Example 1.4,
           Eq. (5.7), printed pp. 21-26 and 196,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    t_i = _temperature(initial_temperature_k, "initial_temperature_k")
    t_inf = _temperature(ambient_temperature_k, "ambient_temperature_k")
    h = _positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
    area = _positive(surface_area_m2, "surface_area_m2")
    volume = _positive(volume_m3, "volume_m3")
    capacity = (_positive(density_kg_m3, "density_kg_m3")
                * _positive(specific_heat_j_kg_k, "specific_heat_j_kg_k") * volume)
    limit = _positive(biot_limit, "biot_limit")
    if isinstance(times_s, str) or not isinstance(times_s, Sequence) or not times_s:
        raise InputError("times_s must be a non-empty sequence.", field="times_s")
    times = tuple(_nonnegative(t, "times_s") for t in times_s)
    tau = capacity / (h * area)
    length = volume / area
    temperatures = tuple(t_inf + (t_i - t_inf) * math.exp(-t / tau) for t in times)
    released = tuple(capacity * (t_i - t) for t in temperatures)
    warnings = []
    biot = None
    satisfied = None
    if thermal_conductivity_w_m_k is not None:
        biot = h * length / _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
        satisfied = biot <= limit
        if not satisfied:
            warnings.append(WarningMessage(
                "lumped_biot_too_large",
                f"Bi = h(V/A)/k = {biot:.4g} exceeds {limit:.3g}; internal gradients are not "
                "negligible and the lumped solution is inaccurate.", "warning",
                "thermal_conductivity_w_m_k"))
    else:
        warnings.append(WarningMessage(
            "lumped_biot_unchecked",
            "No body conductivity was given, so the Bi << 1 requirement was not checked.",
            "warning", "thermal_conductivity_w_m_k"))
    return LumpedCapacitanceResult(times, temperatures, tau, length, biot, released, satisfied,
                                   tuple(warnings))


# ---------------------------------------------------------------------------
# semi-infinite solid
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class SemiInfiniteResult:
    """Temperature and heat flux at one depth and time in a semi-infinite solid.

    Attributes
    ----------
    temperature_k : float
        Temperature at the requested depth, K.
    heat_flux_w_m2 : float or None
        Conductive flux at that depth in the +x direction (into the solid),
        W/m^2; None where the source gives no closed form (convective case
        below the surface).
    surface_temperature_k : float
        Surface temperature at the requested time, K.
    surface_heat_flux_w_m2 : float
        Heat flux entering the solid at the surface, W/m^2.
    similarity_variable : float
        ``zeta/2 = x / (2 sqrt(alpha t))``.
    penetration_depth_99_m : float
        ``delta_99 = 3.64 sqrt(alpha t)`` of eqn. (5.51) (step in surface
        temperature; a scale for the other cases).

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 5.6, Eqs. (5.44)-(5.56),
           printed pp. 222-229,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    temperature_k: float
    heat_flux_w_m2: float | None
    surface_temperature_k: float
    surface_heat_flux_w_m2: float
    similarity_variable: float
    penetration_depth_99_m: float


def _semi_infinite_inputs(depth_m, time_s, diffusivity_m2_s):
    depth = _nonnegative(depth_m, "depth_m")
    time = _positive(time_s, "time_s")
    alpha = _positive(diffusivity_m2_s, "thermal_diffusivity_m2_s")
    root = math.sqrt(alpha * time)
    return depth, time, alpha, root, depth / (2.0 * root)


def semi_infinite_surface_temperature_step(
    *, depth_m: float, time_s: float, thermal_diffusivity_m2_s: float,
    thermal_conductivity_w_m_k: float, initial_temperature_k: float,
    surface_temperature_k: float,
) -> SemiInfiniteResult:
    """Semi-infinite solid after a sudden change of surface temperature.

    Parameters
    ----------
    depth_m : float
        Depth below the surface, m (>= 0).
    time_s : float
        Time since the step, s (> 0).
    thermal_diffusivity_m2_s, thermal_conductivity_w_m_k : float
        ``alpha`` and ``k`` of the solid.
    initial_temperature_k, surface_temperature_k : float
        ``T_i`` and the new surface temperature ``T_s``.

    Returns
    -------
    SemiInfiniteResult
        Temperature, local and surface flux.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    ``(T - T_s)/(T_i - T_s) = erf(x / 2 sqrt(alpha t))`` (5.50); the flux
    follows from Fourier's law with ``dTheta/dzeta = exp(-zeta^2/4)/sqrt(pi)``
    (5.48), which at the surface is (5.54) [1]_. ``erf`` is :func:`math.erf`
    [2]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (5.48), (5.50), (5.51),
           (5.54), printed pp. 223-228,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] PYTHON-3.12-MATH-ERF. Python 3.12, math.erf and math.erfc,
           https://docs.python.org/3.12/library/math.html#math.erfc
    """
    _depth, _, _, root, z = _semi_infinite_inputs(depth_m, time_s, thermal_diffusivity_m2_s)
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    t_i = _temperature(initial_temperature_k, "initial_temperature_k")
    t_s = _temperature(surface_temperature_k, "surface_temperature_k")
    surface_flux = k * (t_s - t_i) / (math.sqrt(math.pi) * root)
    return SemiInfiniteResult(
        temperature_k=t_s + (t_i - t_s) * math.erf(z),
        heat_flux_w_m2=surface_flux * math.exp(-z * z),
        surface_temperature_k=t_s,
        surface_heat_flux_w_m2=surface_flux,
        similarity_variable=z,
        penetration_depth_99_m=3.64 * root,
    )


def semi_infinite_surface_flux_step(
    *, depth_m: float, time_s: float, thermal_diffusivity_m2_s: float,
    thermal_conductivity_w_m_k: float, initial_temperature_k: float,
    surface_heat_flux_w_m2: float,
) -> SemiInfiniteResult:
    """Semi-infinite solid after a sudden constant heat flux ``q_w`` is applied at the surface.

    Parameters
    ----------
    depth_m, time_s : float
        Depth (>= 0, m) and time since the flux was applied (> 0, s).
    thermal_diffusivity_m2_s, thermal_conductivity_w_m_k : float
        ``alpha`` and ``k``.
    initial_temperature_k : float
        ``T_i``, K.
    surface_heat_flux_w_m2 : float
        Flux entering the solid, W/m^2 (negative for cooling).

    Returns
    -------
    SemiInfiniteResult
        Temperature, local and surface flux.

    Raises
    ------
    DomainError
        If an input is outside its domain or the solution falls to 0 K.

    Notes
    -----
    The flux obeys ``(q_w - q)/q_w = erf(x / 2 sqrt(alpha t))`` (5.55).
    Integrating Fourier's law from x to infinity, as the source does at the
    surface to obtain (5.56), gives
    ``T - T_i = (2 q_w / k) sqrt(alpha t) ierfc(x / 2 sqrt(alpha t))`` with
    ``ierfc`` the repeated integral of DLMF 7.18.2; at x = 0,
    ``ierfc(0) = 1/sqrt(pi)`` returns (5.56) [1]_ [2]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (5.55)-(5.56), printed p. 229,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] DLMF-7-ERFC. NIST DLMF Eqs. 7.18.1-7.18.2, https://dlmf.nist.gov/7.18
    """
    _, _, _, root, z = _semi_infinite_inputs(depth_m, time_s, thermal_diffusivity_m2_s)
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    t_i = _temperature(initial_temperature_k, "initial_temperature_k")
    q_w = _finite(surface_heat_flux_w_m2, "surface_heat_flux_w_m2")
    temperature = t_i + 2.0 * q_w / k * root * integrated_erfc(z)
    surface = t_i + 2.0 * q_w / k * root / math.sqrt(math.pi)
    if temperature <= 0.0 or surface <= 0.0:
        raise DomainError("The imposed cooling flux drives the solution below 0 K.",
                          field="surface_heat_flux_w_m2")
    return SemiInfiniteResult(
        temperature_k=temperature,
        heat_flux_w_m2=q_w * math.erfc(z),
        surface_temperature_k=surface,
        surface_heat_flux_w_m2=q_w,
        similarity_variable=z,
        penetration_depth_99_m=3.64 * root,
    )


def semi_infinite_convective_surface(
    *, depth_m: float, time_s: float, thermal_diffusivity_m2_s: float,
    thermal_conductivity_w_m_k: float, heat_transfer_coefficient_w_m2_k: float,
    initial_temperature_k: float, fluid_temperature_k: float,
) -> SemiInfiniteResult:
    """Semi-infinite solid suddenly exposed to a fluid through a film coefficient ``h``.

    Parameters
    ----------
    depth_m, time_s : float
        Depth (>= 0, m) and exposure time (> 0, s).
    thermal_diffusivity_m2_s, thermal_conductivity_w_m_k : float
        ``alpha`` and ``k``.
    heat_transfer_coefficient_w_m2_k : float
        ``h``.
    initial_temperature_k, fluid_temperature_k : float
        ``T_i`` and ``T_inf``.

    Returns
    -------
    SemiInfiniteResult
        Temperature at the depth, surface temperature and surface flux
        ``h (T_inf - T_s)``; the local flux below the surface is not given.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    ``Theta = erf(zeta/2) + exp(beta zeta + beta^2) erfc(zeta/2 + beta)``
    with ``beta = h sqrt(alpha t)/k`` (5.52)-(5.53) [1]_. The second term is
    evaluated as ``exp(-(zeta/2)^2) erfcx(zeta/2 + beta)`` (an exact
    rearrangement) to avoid overflow, with ``erfcx`` from DLMF 7.12.1 for
    large arguments [2]_. The book notes the form also serves curved bodies at
    small Fourier number, with under 6 % (cylinder) and 12 % (sphere) error for
    Fo <~ 0.05.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (5.52)-(5.53), Example 5.5,
           printed pp. 225-228,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] DLMF-7-ERFC. NIST DLMF Eq. 7.12.1, https://dlmf.nist.gov/7.12
    """
    _, _, _, root, z = _semi_infinite_inputs(depth_m, time_s, thermal_diffusivity_m2_s)
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    h = _positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
    t_i = _temperature(initial_temperature_k, "initial_temperature_k")
    t_inf = _temperature(fluid_temperature_k, "fluid_temperature_k")
    beta = h * root / k
    theta = math.erf(z) + math.exp(-z * z) * _erfcx(z + beta)
    theta_surface = _erfcx(beta)
    surface = t_inf + (t_i - t_inf) * theta_surface
    return SemiInfiniteResult(
        temperature_k=t_inf + (t_i - t_inf) * theta,
        heat_flux_w_m2=None,
        surface_temperature_k=surface,
        surface_heat_flux_w_m2=h * (t_inf - surface),
        similarity_variable=z,
        penetration_depth_99_m=3.64 * root,
    )


def contact_temperature(*, temperature_1_k: float, effusivity_1: float,
                        temperature_2_k: float, effusivity_2: float) -> float:
    """Interface temperature of two semi-infinite bodies brought suddenly into contact.

    Parameters
    ----------
    temperature_1_k, temperature_2_k : float
        Initial uniform temperatures, K.
    effusivity_1, effusivity_2 : float
        ``sqrt(k rho c_p)`` of each body, W s^0.5/(m^2 K).

    Returns
    -------
    float
        Constant interface temperature ``T_s``, K.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (5.60), printed p. 232,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    t1 = _temperature(temperature_1_k, "temperature_1_k")
    t2 = _temperature(temperature_2_k, "temperature_2_k")
    e1 = _positive(effusivity_1, "effusivity_1")
    e2 = _positive(effusivity_2, "effusivity_2")
    return t2 + (t1 - t2) * e1 / (e1 + e2)


def semi_infinite_periodic_surface(
    *, depth_m: float, time_s: float, thermal_diffusivity_m2_s: float,
    mean_temperature_k: float, amplitude_k: float, angular_frequency_rad_s: float,
) -> float:
    """Steady-periodic temperature below a surface held at ``T_mean + dT cos(omega t)``.

    Parameters
    ----------
    depth_m, time_s : float
        Depth (>= 0, m) and time (s, any real value).
    thermal_diffusivity_m2_s : float
        ``alpha``.
    mean_temperature_k, amplitude_k : float
        Mean surface temperature and oscillation amplitude (``dT < T_mean``).
    angular_frequency_rad_s : float
        ``omega`` (> 0).

    Returns
    -------
    float
        ``T_mean + dT exp(-xi) cos(omega t - xi)``, ``xi = x sqrt(omega/2 alpha)``.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (5.61)-(5.62), Example 5.7,
           printed pp. 233-235,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    depth = _nonnegative(depth_m, "depth_m")
    time = _finite(time_s, "time_s")
    alpha = _positive(thermal_diffusivity_m2_s, "thermal_diffusivity_m2_s")
    mean = _temperature(mean_temperature_k, "mean_temperature_k")
    amplitude = _nonnegative(amplitude_k, "amplitude_k")
    if amplitude >= mean:
        raise DomainError("The amplitude must be below the mean absolute temperature.",
                          field="amplitude_k")
    omega = _positive(angular_frequency_rad_s, "angular_frequency_rad_s")
    xi = depth * math.sqrt(omega / (2.0 * alpha))
    return mean + amplitude * math.exp(-xi) * math.cos(omega * time - xi)


# ---------------------------------------------------------------------------
# series solutions for slabs, cylinders and spheres with surface convection
# ---------------------------------------------------------------------------

def _check_series_shape(shape: object) -> str:
    if shape not in _SERIES_SHAPES:
        raise InputError(f"shape must be one of {_SERIES_SHAPES}.", field="shape")
    return shape


def _check_biot(biot_number: object) -> float:
    if isinstance(biot_number, float) and biot_number == math.inf:
        return math.inf
    return _positive(biot_number, "biot_number")


def _eigen_function(shape: str, biot: float) -> Callable[[float], float]:
    inverse = 0.0 if biot == math.inf else 1.0 / biot
    if shape == "slab":
        return lambda lam: math.cos(lam) - lam * inverse * math.sin(lam)
    if shape == "cylinder":
        return lambda lam: _bessel_j(0, lam) - lam * inverse * _bessel_j(1, lam)

    def sphere(lam: float) -> float:
        sinc = 1.0 if lam == 0.0 else math.sin(lam) / lam
        return math.cos(lam) * inverse + (1.0 - inverse) * sinc
    return sphere


def _refine_root(function: Callable[[float], float], low: float, high: float) -> float:
    f_low = function(low)
    for _ in range(200):
        middle = 0.5 * (low + high)
        if high - low <= 2.0e-16 * high or middle in (low, high):
            return middle
        f_middle = function(middle)
        if f_middle == 0.0:
            return middle
        if (f_middle > 0.0) == (f_low > 0.0):
            low, f_low = middle, f_middle
        else:
            high = middle
    return 0.5 * (low + high)


@lru_cache(maxsize=512)
def _eigenvalues(shape: str, biot: float, count: int) -> tuple[float, ...]:
    if biot == math.inf and shape == "slab":
        return tuple((n - 0.5) * math.pi for n in range(1, count + 1))
    if biot == math.inf and shape == "sphere":
        return tuple(n * math.pi for n in range(1, count + 1))
    function = _eigen_function(shape, biot)
    roots = []
    if shape in ("slab", "sphere"):
        for n in range(1, count + 1):
            if shape == "slab":
                low, high = (n - 1) * math.pi, (n - 0.5) * math.pi
            else:
                low, high = (n - 1) * math.pi, n * math.pi
            roots.append(_refine_root(function, low, high))
        return tuple(roots)
    step = 0.25
    low, f_low = 0.0, function(0.0)
    while len(roots) < count:
        high = low + step
        f_high = function(high)
        if (f_high > 0.0) != (f_low > 0.0) or f_high == 0.0:
            roots.append(_refine_root(function, low, high))
        low, f_low = high, f_high
    return tuple(roots)


def transient_eigenvalues(shape: str, biot_number: float, count: int) -> tuple[float, ...]:
    """First ``count`` positive eigenvalues of the convectively cooled slab, cylinder or sphere.

    Parameters
    ----------
    shape : {"slab", "cylinder", "sphere"}
        Body shape; the slab has half-thickness L (or thickness L with one
        insulated face), the others outer radius r_o.
    biot_number : float
        ``Bi = h L / k`` or ``h r_o / k``; ``math.inf`` gives the
        fixed-surface-temperature limit.
    count : int
        Number of eigenvalues (1 to 2000).

    Returns
    -------
    tuple of float
        ``lambda_1 < lambda_2 < ...``.

    Raises
    ------
    InputError
        If the shape or count is invalid.
    DomainError
        If the Biot number is not positive.

    Notes
    -----
    Roots of the eigenvalue equations of Table 5.1 of [1]_:
    ``cot lambda = lambda / Bi`` (slab), ``lambda J_1(lambda) = Bi J_0(lambda)``
    (cylinder) and ``lambda cot lambda = 1 - Bi`` (sphere), written in forms
    without poles. The slab and sphere roots lie in ``((n-1) pi, (n-1/2) pi]``
    and ``((n-1) pi, n pi)``; the cylinder roots are bracketed by a scan with
    step 0.25 (their spacing exceeds 2). Each root is refined by bisection to
    a relative width of 2e-16. ``J_0`` and ``J_1`` come from DLMF 10.9.2 [2]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Table 5.1 and Eqs. (5.34)-(5.36),
           printed pp. 206-207,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] DLMF-10-BESSEL. NIST DLMF Eq. 10.9.2, https://dlmf.nist.gov/10.9
    """
    name = _check_series_shape(shape)
    biot = _check_biot(biot_number)
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 2000:
        raise InputError("count must be an integer from 1 to 2000.", field="count")
    return _eigenvalues(name, biot, count)


def _series_terms(shape: str, lam: float) -> tuple[float, float]:
    """Return ``(A_n, D_n / A_n)`` for one eigenvalue (Table 5.1 and volume averages)."""
    if shape == "slab":
        s, c = math.sin(lam), math.cos(lam)
        return 2.0 * s / (lam + s * c), s / lam
    if shape == "cylinder":
        j0, j1 = _bessel_j(0, lam), _bessel_j(1, lam)
        return 2.0 * j1 / (lam * (j0 * j0 + j1 * j1)), 2.0 * j1 / lam
    s, c = math.sin(lam), math.cos(lam)
    return 2.0 * (s - lam * c) / (lam - s * c), 3.0 * (s - lam * c) / lam**3


def _mode_shape(shape: str, lam: float, position: float) -> float:
    argument = lam * position
    if shape == "slab":
        return math.cos(argument)
    if shape == "cylinder":
        return _bessel_j(0, argument)
    return 1.0 if argument == 0.0 else math.sin(argument) / argument


@dataclass(frozen=True, slots=True)
class OneTermCoefficients:
    """First-term eigenvalue and coefficients of the transient series.

    Attributes
    ----------
    shape : str
        ``"slab"``, ``"cylinder"`` or ``"sphere"``.
    biot_number : float
        Biot number.
    eigenvalue : float
        ``lambda_1``.
    amplitude : float
        ``A_1`` of eqn. (5.42).
    energy_coefficient : float
        ``D_1`` of eqn. (5.43).

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Table 5.2, Eqs. (5.42)-(5.43),
           printed pp. 216 and 219,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    shape: str
    biot_number: float
    eigenvalue: float
    amplitude: float
    energy_coefficient: float


def one_term_coefficients(shape: str, biot_number: float) -> OneTermCoefficients:
    """Return ``lambda_1``, ``A_1`` and ``D_1`` (the quantities of Table 5.2).

    Parameters
    ----------
    shape : {"slab", "cylinder", "sphere"}
        Body shape.
    biot_number : float
        ``Bi`` (``math.inf`` allowed).

    Returns
    -------
    OneTermCoefficients
        Eigenvalue and coefficients.

    Raises
    ------
    InputError
        If the shape is invalid.
    DomainError
        If the Biot number is not positive.

    Notes
    -----
    ``A_n`` is from Table 5.1 [1]_. ``D_n`` is the volume average of
    ``A_n f_n`` (the mean temperature of eqn. (5.40)): ``A_n sin(lambda)/lambda``
    for the slab, as printed after eqn. (5.38), and, by the same averaging,
    ``2 A_n J_1(lambda)/lambda`` for the cylinder (using
    ``d(z J_1)/dz = z J_0``, DLMF 10.6.6 [2]_) and
    ``3 A_n (sin lambda - lambda cos lambda)/lambda^3`` for the sphere. The
    tests check all three against Table 5.2.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Tables 5.1-5.2, Eqs. (5.36)-(5.43),
           printed pp. 206-219,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] DLMF-10-BESSEL. NIST DLMF Eq. 10.6.6, https://dlmf.nist.gov/10.6
    """
    name = _check_series_shape(shape)
    biot = _check_biot(biot_number)
    lam = _eigenvalues(name, biot, 1)[0]
    amplitude, ratio = _series_terms(name, lam)
    return OneTermCoefficients(name, biot, lam, amplitude, amplitude * ratio)


@dataclass(frozen=True, slots=True)
class SeriesTransientResult:
    """Multi-term series solution for a convectively cooled slab, cylinder or sphere.

    Attributes
    ----------
    shape : str
        Body shape.
    biot_number, fourier_number, position_ratio : float
        ``Bi``, ``Fo = alpha t / L^2`` and ``x/L`` or ``r/r_o``.
    theta : float
        ``(T - T_inf)/(T_i - T_inf)`` at the position.
    mean_theta : float
        Volume-mean ``Theta`` (eqn. 5.40).
    energy_fraction : float
        ``Phi = 1 - mean_theta``, fraction of the initial excess energy removed
        (eqn. 5.38).
    terms : int
        Number of series terms used.
    truncation_estimate : float
        Magnitude bound of the first omitted term of the local series.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (5.36)-(5.40), printed
           pp. 207-212, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    shape: str
    biot_number: float
    fourier_number: float
    position_ratio: float
    theta: float
    mean_theta: float
    energy_fraction: float
    terms: int
    truncation_estimate: float


def series_transient_solution(shape: str, biot_number: float, fourier_number: float,
                              position_ratio: float = 0.0, *, tolerance: float = 1.0e-13,
                              max_terms: int = 600) -> SeriesTransientResult:
    """Evaluate the series (5.36) and (5.38) for a slab, cylinder or sphere cooled by convection.

    Parameters
    ----------
    shape : {"slab", "cylinder", "sphere"}
        Body shape.
    biot_number : float
        ``Bi`` (``math.inf`` for a fixed surface temperature).
    fourier_number : float
        ``Fo = alpha t / L^2`` (> 0).
    position_ratio : float, optional
        ``x/L`` or ``r/r_o`` in [0, 1].
    tolerance : float, optional
        Stop when ``|A_n| exp(-lambda_n^2 Fo)`` falls below this.
    max_terms : int, optional
        Term limit.

    Returns
    -------
    SeriesTransientResult
        Local and mean ``Theta`` and ``Phi``.

    Raises
    ------
    InputError
        If an argument is invalid.
    DomainError
        If a number is outside its domain.
    ConvergenceError
        If more than ``max_terms`` terms would be needed (very small Fo; use
        the semi-infinite solutions there, as the source advises).

    Notes
    -----
    ``Theta = sum A_n exp(-lambda_n^2 Fo) f_n`` with ``A_n``, ``f_n`` and the
    eigenvalue equations of Table 5.1 [1]_; ``Phi = 1 - sum D_n
    exp(-lambda_n^2 Fo)`` (5.38) with ``D_n`` as in
    :func:`one_term_coefficients`. ``|f_n| <= 1`` for all three shapes, so the
    first omitted term bounds the next contribution; the tail is not summed.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Table 5.1, Eqs. (5.33)-(5.40),
           printed pp. 206-212,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    name = _check_series_shape(shape)
    biot = _check_biot(biot_number)
    fourier = _positive(fourier_number, "fourier_number")
    position = _nonnegative(position_ratio, "position_ratio")
    if position > 1.0:
        raise DomainError("position_ratio must lie in [0, 1].", field="position_ratio")
    tol = _positive(tolerance, "tolerance")
    if isinstance(max_terms, bool) or not isinstance(max_terms, int) or not 1 <= max_terms <= 2000:
        raise InputError("max_terms must be an integer from 1 to 2000.", field="max_terms")
    count = 8
    while True:
        lams = _eigenvalues(name, biot, count)
        bound = abs(_series_terms(name, lams[-1])[0]) * math.exp(-lams[-1] ** 2 * fourier)
        if bound < tol or count >= max_terms:
            break
        count = min(2 * count, max_terms)
    if bound >= tol:
        raise ConvergenceError(
            "The series needs more terms than allowed at this Fourier number; use the "
            "semi-infinite solution for very short times.", field="fourier_number",
            details={"fourier_number": fourier, "terms": count, "bound": bound})
    theta = 0.0
    mean = 0.0
    used = 0
    estimate = 0.0
    for lam in lams:
        amplitude, ratio = _series_terms(name, lam)
        decay = math.exp(-lam * lam * fourier)
        if abs(amplitude) * decay < tol * 1e-3:
            estimate = abs(amplitude) * decay
            break
        theta += amplitude * decay * _mode_shape(name, lam, position)
        mean += amplitude * ratio * decay
        used += 1
    return SeriesTransientResult(name, biot, fourier, position, theta, mean, 1.0 - mean, used,
                                 estimate)


# ---------------------------------------------------------------------------
# transient finite-volume wall solver
# ---------------------------------------------------------------------------

_TRANSIENT_ASSUMPTIONS = (
    "one-dimensional transient conduction, cell-centred finite volumes, uniform cells per layer",
    "theta time integration (theta = 0.5 Crank-Nicolson, theta = 1 backward Euler)",
    "temperature-dependent k and rho c updated by Picard iteration within each step",
    "internal face fluxes from the Kirchhoff (secant) conductivity inside a layer",
    "constant boundary conditions over the run; chain runs for time-varying conditions",
)


def _thomas(lower: list[float], diagonal: list[float], upper: list[float],
            rhs: list[float]) -> list[float]:
    """Solve a tridiagonal system (lower[0] and upper[-1] unused)."""
    n = len(diagonal)
    c_prime = [0.0] * n
    d_prime = [0.0] * n
    c_prime[0] = upper[0] / diagonal[0] if n > 1 else 0.0
    d_prime[0] = rhs[0] / diagonal[0]
    for i in range(1, n):
        denominator = diagonal[i] - lower[i] * c_prime[i - 1]
        c_prime[i] = upper[i] / denominator if i < n - 1 else 0.0
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / denominator
    solution = [0.0] * n
    solution[-1] = d_prime[-1]
    for i in range(n - 2, -1, -1):
        solution[i] = d_prime[i] - c_prime[i] * solution[i + 1]
    return solution


def _secant_conductivity(prop, t_a: float, t_b: float) -> float:
    """``(phi(t_a) - phi(t_b)) / (t_a - t_b)``, the exact Kirchhoff mean of k."""
    if isinstance(prop, _Constant):
        return prop.value
    if abs(t_a - t_b) <= 1.0e-9 * max(t_a, t_b):
        return prop.value_at(0.5 * (t_a + t_b))
    return prop.integral(t_b, t_a) / (t_a - t_b)


def _check_in_table(prop, temperature: float, label: str) -> None:
    if not prop.minimum_temperature_k <= temperature <= prop.maximum_temperature_k:
        raise DomainError(f"{label} reached {temperature:.6g} K, outside the property table "
                          f"[{prop.minimum_temperature_k:.6g}, {prop.maximum_temperature_k:.6g}] K.",
                          code="property_out_of_range", field="layers")


class _BoundaryFlux:
    """Heat rate into a boundary cell through its half cell and the surface condition."""

    __slots__ = ("area", "boundary", "prop", "shape")

    def __init__(self, boundary, prop, shape_integral: float, area: float) -> None:
        self.boundary = boundary
        self.prop = prop
        self.shape = shape_integral
        self.area = area

    def evaluate(self, cell_temperature: float) -> tuple[float, float, float]:
        """Return ``(Q_in, dQ_in/dT_cell, T_surface)``.

        References
        ----------
        .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.34)-(1.35) and
               Example 2.5, printed pp. 33 and 67,
               https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
        """
        boundary, prop, shape = self.boundary, self.prop, self.shape
        k_cell = prop.value_at(cell_temperature)
        if isinstance(boundary, SurfaceTemperature):
            surface = boundary.temperature_k
            _check_in_table(prop, surface, "Surface temperature")
            rate = _secant_conductivity(prop, surface, cell_temperature) * (
                surface - cell_temperature) / shape
            return rate, -k_cell / shape, surface
        if isinstance(boundary, SurfaceHeatFlux):
            rate = boundary.heat_flux_w_m2 * self.area
            surface = prop.solve_integral(cell_temperature, rate * shape)
            if not math.isfinite(surface) or surface <= 0.0:
                raise DomainError("The surface temperature leaves the property table or "
                                  "falls to 0 K.", code="property_out_of_range", field="layers")
            return rate, 0.0, surface
        # convective / radiative surface: solve (phi(Ts) - phi(Tc))/S = A g(Ts)
        candidates = (cell_temperature, boundary.fluid_temperature_k,
                      boundary.surroundings_temperature_k)
        low = max(min(candidates), prop.minimum_temperature_k)
        high = min(max(candidates), prop.maximum_temperature_k)

        def residual(surface: float) -> float:
            return (prop.integral(cell_temperature, surface)
                    - shape * self.area * boundary.flux_into_wall(surface))

        if residual(low) > 0.0 or residual(high) < 0.0:
            raise DomainError("The surface temperature leaves the property table.",
                              code="property_out_of_range", field="layers")
        surface = _bisect_sign(residual, low, high)[0] if high > low else low
        derivative = boundary.flux_derivative(surface)
        rate = self.area * boundary.flux_into_wall(surface)
        slope = self.area * derivative * k_cell / (prop.value_at(surface)
                                                   - shape * self.area * derivative)
        return rate, slope, surface


@dataclass(frozen=True, slots=True)
class TransientWallResult:
    """Finite-volume transient solution of a multilayer wall with an energy audit.

    Attributes
    ----------
    geometry : WallGeometry
        Geometry used.
    cell_centers_m : tuple of float
        Cell-centre positions (x or r), m.
    cell_layers : tuple of int
        Layer index of each cell.
    times_s : tuple of float
        Time levels, from 0 to ``step_count * time_step_s``.
    inner_surface_temperatures_k : tuple of float or None
        Inner-surface temperature at each time level; None for a solid centre.
    outer_surface_temperatures_k : tuple of float
        Outer-surface temperature at each time level.
    inner_heat_rates_w, outer_heat_rates_w : tuple of float
        Heat rate entering the wall through the inner surface and leaving
        through the outer surface at each time level, W.
    snapshot_times_s : tuple of float
        Times of the stored temperature fields.
    snapshot_temperatures_k : tuple of tuple of float
        Cell temperatures at the snapshot times.
    final_temperatures_k : tuple of float
        Cell temperatures at the last time level.
    stored_energy_change_j : float
        ``sum V int_{T_0}^{T_end} rho c dT`` over the cells, J.
    boundary_energy_in_j : float
        Net energy entering through both surfaces, theta-weighted as in the
        scheme, J.
    generated_energy_j : float
        Volumetric generation over the run, J.
    energy_residual_j : float
        ``stored - boundary - generated``; zero up to the Picard tolerance.
    max_picard_iterations : int
        Largest number of iterations used in one step.
    time_step_s, theta : float
        Step size and theta.
    step_count : int
        Number of steps.
    warnings : tuple of WarningMessage
        Model warnings.
    assumptions : tuple of str
        Declared assumptions.
    physical_model_qualified : bool
        Always False.
    model : str
        Model identifier.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eq. (2.10) and Sec. 5.4, printed
           pp. 55 and 207, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    geometry: WallGeometry
    cell_centers_m: tuple[float, ...]
    cell_layers: tuple[int, ...]
    times_s: tuple[float, ...]
    inner_surface_temperatures_k: tuple[float | None, ...]
    outer_surface_temperatures_k: tuple[float, ...]
    inner_heat_rates_w: tuple[float, ...]
    outer_heat_rates_w: tuple[float, ...]
    snapshot_times_s: tuple[float, ...]
    snapshot_temperatures_k: tuple[tuple[float, ...], ...]
    final_temperatures_k: tuple[float, ...]
    stored_energy_change_j: float
    boundary_energy_in_j: float
    generated_energy_j: float
    energy_residual_j: float
    max_picard_iterations: int
    time_step_s: float
    theta: float
    step_count: int
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...] = _TRANSIENT_ASSUMPTIONS
    physical_model_qualified: bool = False
    model: str = "finite_volume_theta_wall_v1"


def _positive_int(value: object, label: str, *, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise InputError(f"{label} must be an integer from 1 to {maximum}.", field=label)
    return value


def solve_transient_wall(
    layers: Sequence[WallLayer] | WallLayer,
    geometry: WallGeometry,
    outer_boundary: SurfaceTemperature | SurfaceHeatFlux | SurfaceConvection,
    inner_boundary: SurfaceTemperature | SurfaceHeatFlux | SurfaceConvection | None = None,
    *,
    initial_temperature_k: float | Sequence[float],
    time_step_s: float,
    step_count: int,
    cells_per_layer: int | Sequence[int] = 20,
    theta: float = 0.5,
    contact_conductances_w_m2_k: Sequence[float | None] | None = None,
    snapshot_every: int = 0,
    iteration_tolerance_k: float = 1.0e-10,
    max_iterations: int = 50,
) -> TransientWallResult:
    """Integrate 1-D transient conduction in a multilayer wall with the theta method.

    Parameters
    ----------
    layers : sequence of WallLayer or WallLayer
        Layers from the inner surface outward; each needs
        ``volumetric_heat_capacity``.
    geometry : WallGeometry
        Plane, cylinder or sphere (``inner_radius_m = 0`` for a solid body).
    outer_boundary, inner_boundary : boundary condition
        Surface conditions; ``inner_boundary=None`` is adiabatic (required at
        a solid centre).
    initial_temperature_k : float or sequence of float
        Uniform initial temperature, or one value per cell (for chaining runs
        with changing boundary conditions).
    time_step_s : float
        Step size, s.
    step_count : int
        Number of steps (1 to 1e6).
    cells_per_layer : int or sequence of int, optional
        Uniform cells in each layer.
    theta : float, optional
        Implicitness in [0.5, 1]; 0.5 is Crank-Nicolson (second order in
        time), 1 is backward Euler (first order). Both are unconditionally
        stable for the linear problem; Crank-Nicolson can ring after a
        sudden boundary change when ``alpha dt / dx^2`` is large.
    contact_conductances_w_m2_k : sequence of float or None, optional
        Interfacial conductance per interface; None is perfect contact.
    snapshot_every : int, optional
        Store the field every this many steps (0: initial and final only).
    iteration_tolerance_k : float, optional
        Picard convergence tolerance on the cell temperatures, K.
    max_iterations : int, optional
        Picard iteration limit per step.

    Returns
    -------
    TransientWallResult
        Surface histories, snapshots, final field and the energy audit.

    Raises
    ------
    InputError
        If an argument is malformed or a layer lacks a heat capacity.
    DomainError
        If a temperature leaves a property table or falls to 0 K.
    ConvergenceError
        If the Picard iteration does not converge in a step.

    Notes
    -----
    Each cell obeys the integral energy balance behind eqn. (2.10) of [1]_,
    ``V (e(T^{n+1}) - e(T^n)) / dt = theta F^{n+1} + (1 - theta) F^n +
    q''' V`` with ``e = int rho c dT`` and ``F`` the net face heat rate. The
    face flux between two cells of one layer is the steady Kirchhoff flux
    ``(phi_j - phi_{j+1}) / int dr/A`` [2]_, so the steady profile of a
    single layer without generation is reproduced to rounding for any
    ``k(T)`` at any grid size; across an interface
    the half-cell resistances (with ``k`` at the cell temperatures) and the
    contact resistance add in series (Sec. 2.3 of [1]_). Boundary cells use
    the same half-cell Kirchhoff relation, solved together with the surface
    balance for convection and radiation (eqns. (1.34)-(1.35)). Within a step
    ``k``, ``rho c`` and the boundary rates are re-evaluated (Picard, with
    Newton linearization of storage and boundary terms) until the
    temperature change is below ``iteration_tolerance_k``. Because every
    internal face flux enters its two cells with opposite signs, the summed
    balance is conservative; the result reports the residual of that
    balance.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (1.34)-(1.35), (2.10),
           Sec. 2.3, printed pp. 33, 55, 75-77,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    .. [2] ENCIT-2012-KIRCHHOFF. W. F. da Silva et al., ENCIT2012-352, Sec. 3,
           Eq. (5), https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf
    """
    owned = _check_layers(layers)
    if not isinstance(geometry, WallGeometry):
        raise InputError("geometry must be a WallGeometry.", field="geometry")
    outer = _check_boundary(outer_boundary, "outer_boundary")
    inner = SurfaceHeatFlux(0.0) if inner_boundary is None else _check_boundary(
        inner_boundary, "inner_boundary")
    contacts = _check_contacts(contact_conductances_w_m2_k, len(owned))
    solid = geometry.kind != "plane" and geometry.inner_radius_m == 0.0
    if solid and not (isinstance(inner, SurfaceHeatFlux) and inner.heat_flux_w_m2 == 0.0):
        raise InputError("A solid cylinder or sphere has a symmetric centre; leave "
                         "inner_boundary as None.", field="inner_boundary")
    dt = _positive(time_step_s, "time_step_s")
    steps = _positive_int(step_count, "step_count", maximum=1_000_000)
    theta_value = _finite(theta, "theta")
    if not 0.5 <= theta_value <= 1.0:
        raise DomainError("theta must lie in [0.5, 1] for unconditional stability.",
                          field="theta")
    tolerance = _positive(iteration_tolerance_k, "iteration_tolerance_k")
    iteration_limit = _positive_int(max_iterations, "max_iterations", maximum=1000)
    if isinstance(snapshot_every, bool) or not isinstance(snapshot_every, int) or (
            snapshot_every < 0):
        raise InputError("snapshot_every must be a non-negative integer.",
                         field="snapshot_every")
    if isinstance(cells_per_layer, int) and not isinstance(cells_per_layer, bool):
        counts = (cells_per_layer,) * len(owned)
    elif isinstance(cells_per_layer, Sequence) and len(cells_per_layer) == len(owned):
        counts = tuple(cells_per_layer)
    else:
        raise InputError("cells_per_layer must be an integer or one integer per layer.",
                         field="cells_per_layer")
    counts = tuple(_positive_int(c, "cells_per_layer", maximum=100_000) for c in counts)
    for layer in owned:
        if layer.volumetric_heat_capacity is None:
            raise InputError("Every layer needs volumetric_heat_capacity for a transient run.",
                             field="volumetric_heat_capacity")

    conductivities = [_property(layer.thermal_conductivity, "thermal_conductivity")
                      for layer in owned]
    capacities = [_property(layer.volumetric_heat_capacity, "volumetric_heat_capacity")
                  for layer in owned]
    faces = _layer_faces(geometry, owned)
    centers, volumes, cell_layers, left_faces, right_faces = [], [], [], [], []
    for index, layer in enumerate(owned):
        width = layer.thickness_m / counts[index]
        for j in range(counts[index]):
            left = faces[index] + j * width
            right = faces[index + 1] if j == counts[index] - 1 else left + width
            left_faces.append(left)
            right_faces.append(right)
            centers.append(0.5 * (left + right))
            volumes.append(geometry.volume(left, right))
            cell_layers.append(index)
    n = len(centers)
    sources = [owned[cell_layers[j]].heat_generation_w_m3 * volumes[j] for j in range(n)]

    face_info = []  # (same_layer, shape_total or (shape_left, contact_res, shape_right))
    for j in range(n - 1):
        face = right_faces[j]
        if cell_layers[j] == cell_layers[j + 1]:
            face_info.append((True, geometry.shape_integral(centers[j], centers[j + 1]), 0.0,
                              0.0))
        else:
            index = cell_layers[j]
            contact = contacts[index]
            resistance = 0.0 if contact == math.inf else 1.0 / (contact * geometry.area(face))
            face_info.append((False, geometry.shape_integral(centers[j], face), resistance,
                              geometry.shape_integral(face, centers[j + 1])))

    inner_flux = None if solid else _BoundaryFlux(
        inner, conductivities[0], geometry.shape_integral(faces[0], centers[0]),
        geometry.area(faces[0]))
    outer_flux = _BoundaryFlux(outer, conductivities[-1],
                               geometry.shape_integral(centers[-1], faces[-1]),
                               geometry.area(faces[-1]))

    if isinstance(initial_temperature_k, Sequence) and not isinstance(initial_temperature_k, str):
        if len(initial_temperature_k) != n:
            raise InputError(f"initial_temperature_k needs {n} cell values.",
                             field="initial_temperature_k")
        temperatures = [_temperature(t, "initial_temperature_k") for t in initial_temperature_k]
    else:
        temperatures = [_temperature(initial_temperature_k, "initial_temperature_k")] * n
    for j, t in enumerate(temperatures):
        _check_in_table(conductivities[cell_layers[j]], t, "Initial temperature")
        _check_in_table(capacities[cell_layers[j]], t, "Initial temperature")

    linear = (all(isinstance(p, _Constant) for p in conductivities + capacities)
              and not any(isinstance(b, SurfaceConvection) and b.emissivity > 0.0
                          for b in (inner, outer)))

    def conductances(field: list[float]) -> list[float]:
        values = []
        for j, (same, shape_a, resistance, shape_b) in enumerate(face_info):
            if same:
                prop = conductivities[cell_layers[j]]
                values.append(_secant_conductivity(prop, field[j], field[j + 1]) / shape_a)
            else:
                k_left = conductivities[cell_layers[j]].value_at(field[j])
                k_right = conductivities[cell_layers[j + 1]].value_at(field[j + 1])
                values.append(1.0 / (shape_a / k_left + resistance + shape_b / k_right))
        return values

    def boundary_terms(field: list[float]):
        inner_terms = (0.0, 0.0, None) if inner_flux is None else inner_flux.evaluate(field[0])
        outer_terms = outer_flux.evaluate(field[-1])
        return inner_terms, outer_terms

    def net_rates(field, g_values, inner_terms, outer_terms) -> list[float]:
        rates = [0.0] * n
        for j in range(n - 1):
            flow = g_values[j] * (field[j] - field[j + 1])
            rates[j] -= flow
            rates[j + 1] += flow
        rates[0] += inner_terms[0]
        rates[-1] += outer_terms[0]
        return rates

    g_old = conductances(temperatures)
    inner_old, outer_old = boundary_terms(temperatures)
    initial_field = list(temperatures)
    times = [0.0]
    inner_surface = [inner_old[2]]
    outer_surface = [outer_old[2]]
    inner_rates = [inner_old[0]]
    outer_rates = [-outer_old[0]]
    snapshot_times = [0.0]
    snapshots = [tuple(temperatures)]
    boundary_energy = 0.0
    generated_energy = 0.0
    max_used = 0
    weight_new, weight_old = theta_value, 1.0 - theta_value
    for step in range(1, steps + 1):
        rates_old = net_rates(temperatures, g_old, inner_old, outer_old)
        iterate = list(temperatures)
        converged = False
        for iteration in range(1, iteration_limit + 1):
            g_now = g_old if iteration == 1 and linear else conductances(iterate)
            inner_now, outer_now = boundary_terms(iterate)
            lower, diagonal, upper, rhs = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
            for j in range(n):
                capacity = capacities[cell_layers[j]]
                storage = volumes[j] * capacity.value_at(iterate[j]) / dt
                diagonal[j] = storage
                rhs[j] = (storage * iterate[j]
                          - volumes[j] * capacity.integral(temperatures[j], iterate[j]) / dt
                          + weight_old * rates_old[j] + sources[j])
                if j > 0:
                    diagonal[j] += weight_new * g_now[j - 1]
                    lower[j] = -weight_new * g_now[j - 1]
                if j < n - 1:
                    diagonal[j] += weight_new * g_now[j]
                    upper[j] = -weight_new * g_now[j]
            diagonal[0] -= weight_new * inner_now[1]
            rhs[0] += weight_new * (inner_now[0] - inner_now[1] * iterate[0])
            diagonal[-1] -= weight_new * outer_now[1]
            rhs[-1] += weight_new * (outer_now[0] - outer_now[1] * iterate[-1])
            solution = _thomas(lower, diagonal, upper, rhs)
            change = max(abs(a - b) for a, b in zip(solution, iterate, strict=True))
            for j, value in enumerate(solution):
                if not value > 0.0:
                    raise DomainError("A cell temperature fell to or below 0 K.",
                                      field="layers")
                _check_in_table(conductivities[cell_layers[j]], value, "Cell temperature")
                _check_in_table(capacities[cell_layers[j]], value, "Cell temperature")
            iterate = solution
            if linear or change <= tolerance:
                converged = True
                break
        if not converged:
            raise ConvergenceError("Picard iteration did not converge within the step.",
                                   field="max_iterations",
                                   details={"step": step, "change_k": change})
        max_used = max(max_used, iteration)
        g_new = conductances(iterate)
        inner_new, outer_new = boundary_terms(iterate)
        boundary_energy += dt * (weight_new * (inner_new[0] + outer_new[0])
                                 + weight_old * (inner_old[0] + outer_old[0]))
        generated_energy += dt * sum(sources)
        temperatures, g_old, inner_old, outer_old = iterate, g_new, inner_new, outer_new
        times.append(step * dt)
        inner_surface.append(inner_new[2])
        outer_surface.append(outer_new[2])
        inner_rates.append(inner_new[0])
        outer_rates.append(-outer_new[0])
        if ((snapshot_every and step % snapshot_every == 0) or step == steps) and snapshot_times[-1] != step * dt:
            snapshot_times.append(step * dt)
            snapshots.append(tuple(temperatures))
    stored = sum(volumes[j] * capacities[cell_layers[j]].integral(initial_field[j],
                                                                  temperatures[j])
                 for j in range(n))
    warnings = [WarningMessage(
        "conduction_model_unqualified",
        "One-dimensional finite-volume conduction with constant boundary conditions; "
        "spatial and temporal errors depend on the grid and step chosen.", "info")]
    if any(isinstance(b, SurfaceConvection) and b.emissivity > 0.0 for b in (inner,)) and (
            geometry.kind != "plane"):
        warnings.append(WarningMessage(
            "concave_surface_radiation",
            "Radiation on a concave inner surface uses F = emittance (AHTT eqn. 1.35).",
            "warning", "inner_boundary"))
    return TransientWallResult(
        geometry=geometry,
        cell_centers_m=tuple(centers),
        cell_layers=tuple(cell_layers),
        times_s=tuple(times),
        inner_surface_temperatures_k=tuple(inner_surface),
        outer_surface_temperatures_k=tuple(outer_surface),
        inner_heat_rates_w=tuple(inner_rates),
        outer_heat_rates_w=tuple(outer_rates),
        snapshot_times_s=tuple(snapshot_times),
        snapshot_temperatures_k=tuple(snapshots),
        final_temperatures_k=tuple(temperatures),
        stored_energy_change_j=stored,
        boundary_energy_in_j=boundary_energy,
        generated_energy_j=generated_energy,
        energy_residual_j=stored - boundary_energy - generated_energy,
        max_picard_iterations=max_used,
        time_step_s=dt,
        theta=theta_value,
        step_count=steps,
        warnings=tuple(warnings),
    )


# ---------------------------------------------------------------------------
# extended surfaces
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class FinResult:
    """Performance of a one-dimensional fin.

    Attributes
    ----------
    heat_rate_w : float
        Heat rate through the fin root, positive from the base into the fin
        when the base is hotter than the fluid, W.
    efficiency : float
        ``eta_f`` (4.52): actual heat rate over the rate if the whole
        convecting surface were at the base temperature.
    effectiveness : float
        ``eps_f`` (4.54): heat rate with the fin over the rate from the base
        area it covers, with the fin's ``h``.
    thermal_resistance_k_w : float
        ``(T_0 - T_inf) / Q`` (4.56)-(4.58), K/W.
    fin_parameter_per_m : float
        ``m``, 1/m.
    tip_temperature_k : float
        Tip (or outer-edge) temperature, K.
    transverse_biot_number : float
        ``h (A/P) / k`` (4.27) for a straight fin, ``h (t/2) / k`` for an
        annular fin; the one-dimensional model needs it to be << 1.
    convecting_area_m2 : float
        Surface used in the efficiency, m^2.
    base_area_m2 : float
        Root cross-section, m^2.
    warnings : tuple of WarningMessage
        Validity warnings.
    model : str
        Model identifier.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 4.5, Eqs. (4.27), (4.52)-(4.58),
           printed pp. 163-175,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """

    heat_rate_w: float
    efficiency: float
    effectiveness: float
    thermal_resistance_k_w: float
    fin_parameter_per_m: float
    tip_temperature_k: float
    transverse_biot_number: float
    convecting_area_m2: float
    base_area_m2: float
    warnings: tuple[WarningMessage, ...]
    model: str


def _fin_biot_warning(biot: float) -> list[WarningMessage]:
    if biot <= 0.1:
        return []
    return [WarningMessage(
        "fin_transverse_biot_large",
        f"Transverse fin Biot number {biot:.3g} is not small; the source requires Bi << 1 "
        "for one-dimensional fin conduction (the 0.1 flag level is this library's choice).",
        "warning", "heat_transfer_coefficient_w_m2_k")]


def corrected_fin_length(length_m: float, thickness_m: float) -> float:
    """Plate-fin length with the tip correction ``L_c = L + t/2``.

    Parameters
    ----------
    length_m : float
        Fin length (width from root to tip), m.
    thickness_m : float
        Fin thickness, m.

    Returns
    -------
    float
        ``L + t/2``, to be used with an adiabatic tip.

    Raises
    ------
    DomainError
        If an input is not positive and finite.

    Notes
    -----
    Harper and Brown account for heat lost from the exposed edge "by adding to
    the real fin width one-half of the fin thickness" and then treat the edge
    as insulated (assumption (3), printed p. 26; ``w' = w + t/2`` with
    ``f' = tanh(a w')/(a w')``, eqn. (5)) [1]_.

    References
    ----------
    .. [1] NACA-TR-158. D. R. Harper 3d and W. B. Brown, NACA Report 158,
           1923, printed pp. 21 and 26,
           https://ntrs.nasa.gov/api/citations/19930091223/downloads/19930091223.pdf
    """
    return _positive(length_m, "length_m") + 0.5 * _positive(thickness_m, "thickness_m")


def straight_fin(
    *, thermal_conductivity_w_m_k: float, heat_transfer_coefficient_w_m2_k: float,
    perimeter_m: float, cross_section_area_m2: float, length_m: float,
    base_temperature_k: float, fluid_temperature_k: float,
    tip_heat_transfer_coefficient_w_m2_k: float | None = None,
) -> FinResult:
    """Straight fin of uniform cross-section (rectangular plate, pin or any shape).

    Parameters
    ----------
    thermal_conductivity_w_m_k : float
        Fin conductivity ``k``.
    heat_transfer_coefficient_w_m2_k : float
        Uniform side film coefficient ``h``.
    perimeter_m, cross_section_area_m2 : float
        ``P`` and ``A`` of the cross-section.
    length_m : float
        Fin length ``L`` (use :func:`corrected_fin_length` for a corrected
        adiabatic-tip estimate).
    base_temperature_k, fluid_temperature_k : float
        ``T_0`` and ``T_inf`` (must differ).
    tip_heat_transfer_coefficient_w_m2_k : float, optional
        Tip coefficient ``h_L``; None means an adiabatic tip.

    Returns
    -------
    FinResult
        Heat rate, efficiency, effectiveness, resistance, tip temperature.

    Raises
    ------
    DomainError
        If an input is outside its domain or ``T_0 = T_inf``.

    Notes
    -----
    ``m = sqrt(hP/kA)`` (4.33a). Adiabatic tip: ``Q_0 = sqrt(kAhP) theta_0
    tanh mL`` (4.44), ``Theta_tip = 1/cosh mL`` (4.45), ``eta_f = tanh(mL)/mL``
    (4.53). Convective tip, ``Bi_ax = h_L L / k``: (4.48)-(4.49), with the
    ideal rate including the tip area. Effectiveness and resistance follow
    (4.54)-(4.58) [1]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Eqs. (4.27)-(4.58), Examples
           4.8-4.9, printed pp. 163-176,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    h = _positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
    perimeter = _positive(perimeter_m, "perimeter_m")
    area = _positive(cross_section_area_m2, "cross_section_area_m2")
    length = _positive(length_m, "length_m")
    t_0 = _temperature(base_temperature_k, "base_temperature_k")
    t_inf = _temperature(fluid_temperature_k, "fluid_temperature_k")
    excess = t_0 - t_inf
    if excess == 0.0:
        raise DomainError("Base and fluid temperatures must differ.", field="base_temperature_k")
    m = math.sqrt(h * perimeter / (k * area))
    ml = m * length
    scale = math.sqrt(k * area * h * perimeter)
    tanh_ml = math.tanh(ml)
    side_area = perimeter * length
    if tip_heat_transfer_coefficient_w_m2_k is None:
        ratio = tanh_ml
        tip_theta = 1.0 / math.cosh(ml) if ml < 700.0 else 0.0
        ideal = h * side_area
        model = "straight_fin_adiabatic_tip_v1"
    else:
        h_tip = _nonnegative(tip_heat_transfer_coefficient_w_m2_k,
                             "tip_heat_transfer_coefficient_w_m2_k")
        b = h_tip * length / k / ml
        ratio = (b + tanh_ml) / (1.0 + b * tanh_ml)
        tip_theta = (1.0 / (math.cosh(ml) + b * math.sinh(ml))) if ml < 700.0 else 0.0
        ideal = h * side_area + h_tip * area
        model = "straight_fin_convective_tip_v1"
    rate = scale * excess * ratio
    biot = h * (area / perimeter) / k
    return FinResult(
        heat_rate_w=rate,
        efficiency=rate / (ideal * excess),
        effectiveness=rate / (h * area * excess),
        thermal_resistance_k_w=excess / rate,
        fin_parameter_per_m=m,
        tip_temperature_k=t_inf + tip_theta * excess,
        transverse_biot_number=biot,
        convecting_area_m2=side_area if tip_heat_transfer_coefficient_w_m2_k is None
        else side_area + area,
        base_area_m2=area,
        warnings=tuple(_fin_biot_warning(biot)),
        model=model,
    )


def pin_fin(
    *, thermal_conductivity_w_m_k: float, heat_transfer_coefficient_w_m2_k: float,
    diameter_m: float, length_m: float, base_temperature_k: float,
    fluid_temperature_k: float, tip_heat_transfer_coefficient_w_m2_k: float | None = None,
) -> FinResult:
    """Cylindrical pin fin: :func:`straight_fin` with ``P = pi D`` and ``A = pi D^2/4``.

    Parameters
    ----------
    thermal_conductivity_w_m_k, heat_transfer_coefficient_w_m2_k : float
        ``k`` and ``h``.
    diameter_m, length_m : float
        Pin diameter and length.
    base_temperature_k, fluid_temperature_k : float
        ``T_0`` and ``T_inf``.
    tip_heat_transfer_coefficient_w_m2_k : float, optional
        Tip coefficient; None for an adiabatic tip.

    Returns
    -------
    FinResult
        As :func:`straight_fin`.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    For a circular rod ``A/P = R/2`` [1]_.

    References
    ----------
    .. [1] AHTT-V6. Lienhard and Lienhard, Sec. 4.5, Eq. (4.27) and Example
           4.8, printed pp. 163-172,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    diameter = _positive(diameter_m, "diameter_m")
    return straight_fin(
        thermal_conductivity_w_m_k=thermal_conductivity_w_m_k,
        heat_transfer_coefficient_w_m2_k=heat_transfer_coefficient_w_m2_k,
        perimeter_m=math.pi * diameter, cross_section_area_m2=0.25 * math.pi * diameter**2,
        length_m=length_m, base_temperature_k=base_temperature_k,
        fluid_temperature_k=fluid_temperature_k,
        tip_heat_transfer_coefficient_w_m2_k=tip_heat_transfer_coefficient_w_m2_k)


def annular_fin(
    *, thermal_conductivity_w_m_k: float, heat_transfer_coefficient_w_m2_k: float,
    thickness_m: float, inner_radius_m: float, outer_radius_m: float,
    base_temperature_k: float, fluid_temperature_k: float,
    edge_heat_transfer_coefficient_w_m2_k: float = 0.0,
) -> FinResult:
    """Circular annular fin of uniform thickness on a tube (one-dimensional, radial).

    Parameters
    ----------
    thermal_conductivity_w_m_k : float
        Fin conductivity ``k``.
    heat_transfer_coefficient_w_m2_k : float
        Film coefficient ``h`` on both faces.
    thickness_m : float
        Fin thickness ``2t``.
    inner_radius_m, outer_radius_m : float
        Root radius ``a`` and outer radius ``b``.
    base_temperature_k, fluid_temperature_k : float
        Root temperature (perfect root contact) and fluid temperature.
    edge_heat_transfer_coefficient_w_m2_k : float, optional
        Outer-edge coefficient ``h_e``; 0 is an insulated edge.

    Returns
    -------
    FinResult
        Heat rate, efficiency (4.52), effectiveness, resistance and the
        outer-edge temperature.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    Yovanovich et al. solve ``(1/rho) d/drho (rho dpsi/drho) - m^2 psi = 0``
    (their eqn. (26)) with ``psi = C_1 I_0(m rho) + C_2 K_0(m rho)`` (29), the
    edge condition (28) giving (31), and ``Q = -4 pi k a theta_b dpsi/drho``
    (37); efficiency is ``Q / Q_i`` with the ideal rate (39), which includes
    the edge [1]_. Here the same solution is written in dimensional radius
    with ``m = sqrt(h / (k t))``. With perfect root contact ``psi(a) = 1``;
    the printed eqns. (32) and (35) show a factor ``m`` multiplying
    ``Omega K_0(m alpha)`` that is inconsistent with (29) at ``rho = alpha``,
    so the root condition is applied to (29) directly. The constants are
    obtained from a 2x2 system in exponentially scaled Bessel functions (DLMF
    10.25.2, 10.31, 10.32.9, 10.40) so large ``m b`` does not overflow; the
    derivative uses ``I_0' = I_1``, ``K_0' = -K_1`` (DLMF 10.29.3) [2]_.

    References
    ----------
    .. [1] YOVANOVICH-1988-ANNULAR-FIN. M. M. Yovanovich, J. R. Culham and
           T. F. Lemczyk, "Simplified Solutions to Circular Annular Fins with
           Contact Resistance and End Cooling", J. Thermophysics 2(2),
           152-157, 1988, Eqs. (26)-(41),
           https://uwaterloo.ca/microelectronics-heat-transfer-laboratory/sites/default/files/2024-03/mhtl88-1.pdf
    .. [2] DLMF-10-BESSEL. NIST DLMF Eqs. 10.25.2, 10.29.3, 10.31.1-2,
           10.32.9, 10.40.1-2, https://dlmf.nist.gov/10
    .. [3] AHTT-V6. Lienhard and Lienhard, Eqs. (4.52)-(4.58), Example 4.10,
           printed pp. 174-180,
           https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf
    """
    k = _positive(thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k")
    h = _positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
    thickness = _positive(thickness_m, "thickness_m")
    a = _positive(inner_radius_m, "inner_radius_m")
    b = _positive(outer_radius_m, "outer_radius_m")
    if b <= a:
        raise DomainError("Outer radius must exceed the root radius.", field="outer_radius_m")
    h_edge = _nonnegative(edge_heat_transfer_coefficient_w_m2_k,
                          "edge_heat_transfer_coefficient_w_m2_k")
    t_0 = _temperature(base_temperature_k, "base_temperature_k")
    t_inf = _temperature(fluid_temperature_k, "fluid_temperature_k")
    excess = t_0 - t_inf
    if excess == 0.0:
        raise DomainError("Base and fluid temperatures must differ.", field="base_temperature_k")
    half = 0.5 * thickness
    m = math.sqrt(h / (k * half))
    edge = h_edge / k
    x_a, x_b = m * a, m * b
    decay = math.exp(-(x_b - x_a))
    ie0a, ie1a = _scaled_i(0, x_a), _scaled_i(1, x_a)
    ke0a, ke1a = _scaled_k(0, x_a), _scaled_k(1, x_a)
    ie0b, ie1b = _scaled_i(0, x_b), _scaled_i(1, x_b)
    ke0b, ke1b = _scaled_k(0, x_b), _scaled_k(1, x_b)
    # P = C1 exp(m b), R = C2 exp(-m a):
    #   [decay*ie0a, ke0a] [P]   [1]
    #   [alpha2, decay*beta2] [R] = [0]
    alpha2 = m * ie1b + edge * ie0b
    beta2 = -m * ke1b + edge * ke0b
    determinant = decay * decay * ie0a * beta2 - ke0a * alpha2
    p_value = decay * beta2 / determinant
    r_value = -alpha2 / determinant
    gradient = m * (p_value * decay * ie1a - r_value * ke1a)  # dpsi/dr at r = a
    rate = -2.0 * math.pi * a * thickness * k * excess * gradient
    edge_theta = p_value * ie0b + r_value * decay * ke0b
    ideal = h * 2.0 * math.pi * (b * b - a * a) + h_edge * 2.0 * math.pi * b * thickness
    base_area = 2.0 * math.pi * a * thickness
    biot = h * half / k
    return FinResult(
        heat_rate_w=rate,
        efficiency=rate / (ideal * excess),
        effectiveness=rate / (h * base_area * excess),
        thermal_resistance_k_w=excess / rate,
        fin_parameter_per_m=m,
        tip_temperature_k=t_inf + edge_theta * excess,
        transverse_biot_number=biot,
        convecting_area_m2=2.0 * math.pi * (b * b - a * a)
        + (2.0 * math.pi * b * thickness if h_edge else 0.0),
        base_area_m2=base_area,
        warnings=tuple(_fin_biot_warning(biot)),
        model="annular_fin_bessel_v1",
    )
