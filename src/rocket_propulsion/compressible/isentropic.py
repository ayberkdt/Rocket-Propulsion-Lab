"""Isentropic relations for a calorically perfect gas.

The module evaluates the one-dimensional isentropic ratios of NACA Report 1135
[1]_ and inverts them for Mach number.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, degrees, isfinite, radians, sin, sqrt

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import real, text

from .isentropic_inverse import ratio_inverse
from .isentropic_numerics import precise_area, precise_ratios, requires_precise_ratios
from .sonic_numerics import sonic_area, sonic_inverse, sonic_pm

#: Supported numerical inverse bracket; this is an implementation limit,
#: not a physical validity boundary for a perfect gas.
_MINIMUM_MACH = 1.0e-8
_MAXIMUM_MACH = 1000.0


@dataclass(frozen=True, slots=True)
class IsentropicResult:
    """Dimensionless isentropic properties associated with a Mach number [1]_.

    Attributes
    ----------
    mach : float
        Local Mach number.
    mach_angle_deg : float or None
        Mach angle ``asin(1/M)`` in degrees; ``None`` for subsonic flow.
    prandtl_meyer_deg : float or None
        Prandtl-Meyer function in degrees; ``None`` for subsonic flow.
    pressure_total_ratio, density_total_ratio, temperature_total_ratio : float
        Static-to-stagnation ratios ``p/p0``, ``rho/rho0`` and ``T/T0``.
    pressure_critical_ratio, density_critical_ratio, temperature_critical_ratio : float
        Ratios ``p/p*``, ``rho/rho*`` and ``T/T*`` to the sonic state at the
        same stagnation condition.
    area_critical_ratio : float
        Area ratio ``A/A*`` for the same mass flow.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    mach: float
    mach_angle_deg: float | None
    prandtl_meyer_deg: float | None
    pressure_total_ratio: float
    density_total_ratio: float
    temperature_total_ratio: float
    pressure_critical_ratio: float
    density_critical_ratio: float
    temperature_critical_ratio: float
    area_critical_ratio: float


def _validate_gamma(gamma: float) -> None:
    """Raise :class:`DomainError` unless ``gamma`` is finite and greater than one."""
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.")


def _validate_mach(mach: float) -> None:
    """Raise :class:`DomainError` unless ``mach`` is finite and positive."""
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and greater than 0.")


def calculate_isentropic(mach: float, gamma: float = 1.4) -> IsentropicResult:
    """Calculate perfect-gas isentropic flow relations at one Mach number.

    Parameters
    ----------
    mach : float
        Local Mach number, finite and strictly positive.
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    IsentropicResult
        Immutable set of flow ratios and, for ``M >= 1``, angles in degrees.

    Raises
    ------
    DomainError
        If ``mach`` or ``gamma`` is outside its domain, or a positive output
        ratio is not representable as a finite nonzero float. The latter
        carries ``numeric_range`` and the affected computed field.

    Notes
    -----
    With ``X = 1 + (gamma - 1) M^2 / 2`` the stagnation ratios are
    ``T/T0 = 1/X``, ``p/p0 = X^(-gamma/(gamma-1))`` and
    ``rho/rho0 = X^(-1/(gamma-1))`` [1]_. Starred properties are referenced to
    the sonic state at the same stagnation condition, and
    ``A/A* = (1/M) (2 X/(gamma+1))^((gamma+1)/(2(gamma-1)))``. The Mach angle
    ``mu = asin(1/M)`` and the Prandtl-Meyer function
    ``nu = sqrt((gamma+1)/(gamma-1)) atan(sqrt((gamma-1)/(gamma+1) (M^2-1)))
    - atan(sqrt(M^2-1))`` [1]_ are defined only for ``M >= 1`` and are
    returned as ``None`` below sonic speed.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    gamma = real(gamma, 'gamma')
    mach = real(mach, 'mach')
    _validate_gamma(gamma)
    _validate_mach(mach)

    if requires_precise_ratios(mach, gamma):
        ratios = precise_ratios(mach, gamma)
    else:
        energy_factor = 1.0 + 0.5 * (gamma - 1.0) * mach * mach
        critical_factor = (gamma + 1.0) / (2.0 * energy_factor)
        ratios = {
            'temperature_total_ratio': energy_factor**-1.0,
            'pressure_total_ratio': energy_factor ** (-gamma / (gamma - 1.0)),
            'density_total_ratio': energy_factor ** (-1.0 / (gamma - 1.0)),
            'temperature_critical_ratio': critical_factor,
            'pressure_critical_ratio': critical_factor ** (gamma / (gamma - 1.0)),
            'density_critical_ratio': critical_factor ** (1.0 / (gamma - 1.0)),
            'area_critical_ratio': _area_ratio(mach, gamma),
        }
        if any(not isfinite(value) or value <= 0 for value in ratios.values()):
            ratios = precise_ratios(mach, gamma)
    if mach >= 1.0:
        mach_angle = degrees(asin(1.0 / mach))
        prandtl_meyer = _prandtl_meyer_deg(mach, gamma)
    else:
        mach_angle = None
        prandtl_meyer = None

    return IsentropicResult(
        mach=mach,
        mach_angle_deg=mach_angle,
        prandtl_meyer_deg=prandtl_meyer,
        **ratios,
    )


def mach_from_isentropic(
    input_kind: str,
    value: float,
    gamma: float = 1.4,
) -> float:
    """Recover Mach number from one isentropic property.

    Parameters
    ----------
    input_kind : str
        Name of the supplied property: ``mach``, ``mach_angle_deg``,
        ``prandtl_meyer_deg``, ``pressure_total_ratio``,
        ``density_total_ratio``, ``temperature_total_ratio``,
        ``pressure_critical_ratio``, ``density_critical_ratio``,
        ``temperature_critical_ratio``, ``area_ratio_subsonic`` or
        ``area_ratio_supersonic``.
    value : float
        Finite value of that property (angles in degrees).
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    float
        Mach inferred from the supplied rounded float, with the conditioning
        limits described below. A rounded forward value need not recover the
        Mach that generated it.

    Raises
    ------
    DomainError
        If ``input_kind`` is unsupported, ``value`` is not finite, or
        ``value`` lies outside the range of the selected property: Mach angle
        in ``(0, 90]`` degrees, Prandtl-Meyer angle below its ``M -> infinity``
        limit, positive flow ratios, or ``A/A* >= 1``.

    Notes
    -----
    The implementation uses the following numerical inverses of [1]_. Rounded
    forward plateaus and conditioning can limit recoverable Mach precision:

    * Mach angle: ``M = 1/sin(mu)``;
    * stagnation ratios, with ``X = 1 + (gamma-1) M^2/2``:
      ``M^2 = 2 (X - 1)/(gamma - 1)`` and ``X - 1`` from ``T/T0 = 1/X``,
      ``p/p0 = X^(-gamma/(gamma-1))`` or ``rho/rho0 = X^(-1/(gamma-1))``
      in closed form using private complete logarithmic expressions;
    * critical ratios: retain the sonic factor ``(gamma+1)/2`` before the
      logarithm, without first rounding the sonic temperature ratio;
    * area ratio and Prandtl-Meyer angle: Newton iteration with the analytic
      derivatives ``d(A/A*)/dM = (A/A*)(M^2-1)/(M X)`` and
      ``d nu/dM = sqrt(M^2-1)/(M X)``, kept inside a shrinking bracket
      (bisection step whenever Newton would leave it) and stopped when the
      bracket or the step falls below 4 ulp of ``M``.

    Derived results lie in ``[1e-8, 1000]``; the six scalar flow-ratio inverses
    compare their complete logarithmic root with these bounds before float
    export. Direct ``input_kind='mach'`` validates and returns any positive
    finite Mach without imposing this inverse bracket. The area ratio is
    double-valued, so the subsonic or supersonic branch must be named explicitly;
    ``A/A* = 1`` returns
    ``M = 1`` on either branch.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    gamma = real(gamma, 'gamma')
    value = real(value, 'value')
    input_kind = text(input_kind, 'input_kind')
    _validate_gamma(gamma)
    if input_kind == "mach":
        _validate_mach(value)
        return value
    if input_kind == "mach_angle_deg":
        if not 0.0 < value <= 90.0:
            raise DomainError("Mach angle must be in the interval (0, 90].")
        return _bounded(1.0 / sin(radians(value)), input_kind, value)
    if input_kind == "prandtl_meyer_deg":
        maximum = 180.0/(gamma-1)/(sqrt((gamma+1)/(gamma-1))+1)
        if not 0.0 <= value < maximum:
            raise DomainError(f"Prandtl-Meyer angle must be in [0, {maximum:g}).")
        if value == 0.0:
            return 1.0
        if value <= 1.0e-4:
            return sonic_inverse(value, gamma, area=False,
                                  minimum_mach=_MINIMUM_MACH, maximum_mach=_MAXIMUM_MACH)
        return _newton_bracketed(lambda m: _prandtl_meyer_deg(m, gamma), _prandtl_meyer_slope(gamma),
                                 value, 1.0, _MAXIMUM_MACH, input_kind)

    if input_kind in {
        'pressure_total_ratio', 'pressure_critical_ratio',
        'density_total_ratio', 'density_critical_ratio',
        'temperature_total_ratio', 'temperature_critical_ratio',
    }:
        if value <= 0.0:
            raise DomainError('Flow ratios must be greater than zero.', field='value')
        return ratio_inverse(input_kind, value, gamma, _MINIMUM_MACH, _MAXIMUM_MACH)

    if input_kind in {"area_ratio_subsonic", "area_ratio_supersonic"}:
        if value < 1.0:
            raise DomainError("A/A* must be greater than or equal to 1.")
        if value == 1.0:
            return 1.0
        if value <= 1.001:
            return sonic_inverse(value, gamma, area=True,
                                  subsonic=input_kind.endswith('subsonic'),
                                  minimum_mach=_MINIMUM_MACH, maximum_mach=_MAXIMUM_MACH)
        lower, upper = (_MINIMUM_MACH, 1.0) if input_kind.endswith("subsonic") else (1.0, _MAXIMUM_MACH)
        return _newton_bracketed(lambda m: _area_ratio(m, gamma), _area_ratio_slope(gamma),
                                 value, lower, upper, input_kind)

    raise DomainError(f"Unsupported isentropic input kind: {input_kind!r}.")


def _bounded(mach: float, kind: str, value: float) -> float:
    """Return ``mach`` if it lies in the inverse bracket, else refuse."""
    if not (_MINIMUM_MACH <= mach <= _MAXIMUM_MACH):
        raise DomainError(f"{kind} = {value:g} corresponds to a Mach number outside "
                          f"[{_MINIMUM_MACH:g}, {_MAXIMUM_MACH:g}].", code="outside_inverse_range",
                          field="value")
    return mach


def _area_ratio(mach: float, gamma: float) -> float:
    if abs(mach-1.0) <= .01:
        return sonic_area(mach, gamma)
    if requires_precise_ratios(mach, gamma):
        return precise_area(mach, gamma)
    x = 1.0 + 0.5 * (gamma - 1.0) * mach * mach
    result = (2.0 * x / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))) / mach
    return result if isfinite(result) and result > 0 else precise_area(mach, gamma)


def _area_ratio_slope(gamma: float):
    def slope(mach: float, area: float) -> float:
        return area * (mach * mach - 1.0) / (mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach))
    return slope


def _prandtl_meyer_deg(mach: float, gamma: float) -> float:
    if 1.0 <= mach <= 1.01 or gamma >= 1e6:
        return sonic_pm(mach, gamma)
    root = sqrt(max(0.0, (mach-1.0)*(mach+1.0)))
    return degrees(sqrt((gamma + 1.0) / (gamma - 1.0)) * atan(sqrt((gamma - 1.0) / (gamma + 1.0)) * root)
                   - atan(root))


def _prandtl_meyer_slope(gamma: float):
    def slope(mach: float, _value: float) -> float:
        return degrees(sqrt(max(0.0, mach * mach - 1.0)) / (mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach)))
    return slope


def _newton_bracketed(function, slope, target: float, lower: float, upper: float, kind: str) -> float:
    """Solve ``function(M) = target`` for monotonic ``function`` on ``[lower, upper]``.

    Newton steps with an analytic slope, replaced by a bisection step whenever
    they would leave the current bracket; the bracket is updated from the sign
    of every residual. Newton convergence is local and depends on a smooth,
    well-resolved callback. Stops at a sampled zero or when
    the bracket or step is below 4 ulp of ``M``. Rounded callback plateaus are
    not a mathematical enclosure; exhausting the iteration budget is refused.
    """
    f_low, f_high = function(lower) - target, function(upper) - target
    if f_low == 0.0:
        return lower
    if f_high == 0.0:
        return upper
    if (f_low > 0.0) == (f_high > 0.0):
        raise DomainError(f"{kind} = {target:g} corresponds to a Mach number outside "
                          f"[{lower:g}, {upper:g}].", code="outside_inverse_range", field="value")
    increasing = f_high > 0.0
    mach = 0.5 * (lower + upper) if upper < 10.0 else sqrt(lower * upper)
    for _ in range(200):
        value = function(mach)
        residual = value - target
        if residual == 0.0:
            return mach
        if (residual > 0.0) == increasing:
            upper = mach
        else:
            lower = mach
        derivative = slope(mach, value)
        step_ok = derivative != 0.0 and isfinite(derivative)
        candidate = mach - residual / derivative if step_ok else 0.5 * (lower + upper)
        if not (lower < candidate < upper):
            candidate = 0.5 * (lower + upper)
        if abs(candidate - mach) <= 4.0 * 2.220446049250313e-16 * mach or upper - lower <= 4.0 * 2.220446049250313e-16 * upper:
            return candidate
        mach = candidate
    raise ConvergenceError('Isentropic inverse exhausted its iteration budget.',
                           code='isentropic_inverse_not_converged', field='value')
