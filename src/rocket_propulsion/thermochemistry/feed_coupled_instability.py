"""Feed-coupled low-frequency (chugging) stability: injection admittances and the Nyquist test.

The single- and double-time-lag chamber models of NASA SP-194 Sect. 5.2 are
coupled to the dynamic response of each propellant's feed system, written as
an injection admittance

    ``G = -(p_bar / m_dot)(m_i' / p')``

(normalised by the total flow ``m_dot``) and built component by component
from the propellant tank to the injector orifices with the admittance ratios
of Sect. 5.4.1 [1]_:

* constant-area line, Eq. (5.4.1-2), with characteristic impedance
  ``Z_l = a_L m_dot/(A_l p_c)`` and time ``t_l = l/a_L`` (Eqs. 5.4.1-3/4) and
  the wall-compliant sound speed of Eq. (5.4.1-5);
* orifice or valve with inertance, Eq. (5.4.1-9), with
  ``Z_o = (2 Dp/m_L)(m_dot/p_c)`` and ``t_o = l_o m_L/(2 Dp A_o)``;
* manifold capacitance, Eq. (5.4.1-16), ``C_m = rho_L V_m p_c/(K_L m_dot)``;
* propellant tank as a capacitance, starting from ``G_0 = 0`` at the top of
  the tank (Sect. 5.4.2.1), or as a constant-pressure tank (``C_t -> inf``).

Stability of ``F(s) = 1 + f(s)`` with

    ``f(s) = exp(-s tau_m)/(1 + s theta_c) {[1 + ((1+r)/c*) dc*/dr] G_ox e^(-s tau_vox)
    + [1 - (r(1+r)/c*) dc*/dr] G_F e^(-s tau_vF)}``

(Eq. 5.2.3-3; ``tau_vox = tau_vF = 0`` gives the single-lag Eq. 5.2.2-6) is
decided by the Nyquist criterion SP-194 describes in Sect. 5.5.1: with no
open-loop poles in the right half-plane, the number of right-half-plane zeros
of ``F`` equals the number of clockwise encirclements of the origin by
``F(i omega)``. The encirclements are counted from the unwrapped phase of
``F(i omega)`` on a contour resolved with amplitude and derivative bounds.
An unresolved contour raises an error and has no stability verdict.

Limits: linear (small-amplitude) theory; one combustion time lag per
propellant, constant in time; low Mach number in the feed lines (Eq. 5.4.1-2,
not the Mach-corrected 5.4.1-2a); pumps, branch lines and lines with area
change are not modelled; passive feed components are assumed to add no
right-half-plane poles to ``f``.

References
----------
.. [1] NASA-SP-194-1972. D. T. Harrje and F. H. Reardon (eds.), Liquid
       Propellant Rocket Combustion Instability, NASA SP-194, 1972, Sects.
       5.2, 5.4 and 5.5.1, pp. 234-253,
       https://ntrs.nasa.gov/citations/19720026079
"""

from __future__ import annotations

import cmath
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import cos, hypot, pi, sin, sqrt
from typing import TypedDict

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import numeric_range_guard, real
from rocket_propulsion.thermochemistry.nyquist_numerics import resolve_contour

__all__ = [
    "AdmittanceBounds",
    "FeedCoupledStability",
    "FeedSystem",
    "FeedSystemResponse",
    "feed_coupled_stability",
    "line_admittance_ratio",
    "manifold_admittance_ratio",
    "orifice_admittance_ratio",
    "pressurized_tank_injection_admittance",
    "wall_compliant_sound_speed",
]


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _nonnegative(value: object, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be nonnegative.", field=label)
    return number


# ---------------------------------------------------------------------------
# Component admittance ratios (SP-194 Sect. 5.4.1)
# ---------------------------------------------------------------------------


def wall_compliant_sound_speed(*, liquid_sound_speed_m_s: float, liquid_density_kg_m3: float,
                               line_diameter_m: float, wall_thickness_m: float,
                               wall_elastic_modulus_pa: float) -> float:
    """Effective sound speed in a liquid line with elastic walls, Eq. (5.4.1-5).

    ``a_L = a_L,inf / sqrt(1 + rho_L a_L,inf^2 D_l / (E delta_w))``.

    Parameters
    ----------
    liquid_sound_speed_m_s : float
        Sound speed in the unconfined liquid, ``a_L,inf``.
    liquid_density_kg_m3 : float
        ``rho_L``.
    line_diameter_m, wall_thickness_m : float
        Mean pipe diameter ``D_l`` and wall thickness ``delta_w``.
    wall_elastic_modulus_pa : float
        Elastic modulus ``E`` of the pipe material.

    Returns
    -------
    float
        ``a_L`` in m/s.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (5.4.1-5), p. 242,
           https://ntrs.nasa.gov/citations/19720026079
    """
    a_inf = _positive(liquid_sound_speed_m_s, "liquid_sound_speed_m_s")
    rho = _positive(liquid_density_kg_m3, "liquid_density_kg_m3")
    d = _positive(line_diameter_m, "line_diameter_m")
    t = _positive(wall_thickness_m, "wall_thickness_m")
    e = _positive(wall_elastic_modulus_pa, "wall_elastic_modulus_pa")
    return a_inf / sqrt(1.0 + rho * a_inf * a_inf * d / (e * t))


def line_admittance_ratio(upstream_admittance: complex | None, s: complex, *,
                          line_impedance: float, line_time_s: float) -> complex:
    """Admittance at the downstream end of a constant-area lossless line, Eq. (5.4.1-2).

    ``G_2 = (G_1 + tanh(s t_l)/Z_l) / (1 + G_1 Z_l tanh(s t_l))``, the printed
    ratio ``G_2/G_1`` multiplied out so that ``G_1 = 0`` and ``G_1 -> inf``
    (``upstream_admittance=None``, a constant-pressure end) are both exact.

    Parameters
    ----------
    upstream_admittance : complex or None
        ``G_1`` at the upstream end; None for an infinite admittance.
    s : complex
        Laplace variable ``lambda + i omega`` in 1/s.
    line_impedance : float
        ``Z_l = a_L m_dot / (A_l p_c)``.
    line_time_s : float
        ``t_l = l / a_L``.

    Returns
    -------
    complex
        ``G_2``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.1-2)-(5.4.1-4),
           p. 242, https://ntrs.nasa.gov/citations/19720026079
    """
    z = _positive(line_impedance, "line_impedance")
    t = _nonnegative(line_time_s, "line_time_s")
    th = cmath.tanh(s * t)
    if upstream_admittance is None:
        if th == 0:
            return complex("inf")
        return 1.0 / (z * th)
    g1 = complex(upstream_admittance)
    return (g1 + th / z) / (1.0 + g1 * z * th)


def orifice_admittance_ratio(upstream_admittance: complex | None, s: complex, *,
                             orifice_impedance: float, orifice_time_s: float = 0.0) -> complex:
    """Admittance downstream of an orifice or valve, Eq. (5.4.1-9).

    ``G_2 = G_1 / (1 + Z_o G_1 (1 + s t_o))``; ``t_o = 0`` is the short orifice
    of Eq. (5.4.1-12), and ``G_1 -> inf`` gives ``1/(Z_o (1 + s t_o))``.

    Parameters
    ----------
    upstream_admittance : complex or None
        ``G_1``; None for an infinite admittance.
    s : complex
        Laplace variable in 1/s.
    orifice_impedance : float
        ``Z_o = (2 Dp / m_L)(m_dot / p_c)``.
    orifice_time_s : float, optional
        ``t_o = l_o m_L / (2 Dp A_o)``.

    Returns
    -------
    complex
        ``G_2``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.1-9)-(5.4.1-12),
           p. 243, https://ntrs.nasa.gov/citations/19720026079
    """
    z = _positive(orifice_impedance, "orifice_impedance")
    t = _nonnegative(orifice_time_s, "orifice_time_s")
    if upstream_admittance is None:
        return 1.0 / (z * (1.0 + s * t))
    g1 = complex(upstream_admittance)
    return g1 / (1.0 + z * g1 * (1.0 + s * t))


def manifold_admittance_ratio(upstream_admittance: complex | None, s: complex, *,
                              manifold_capacitance_s: float) -> complex | None:
    """Admittance downstream of a manifold capacitance, Eq. (5.4.1-16): ``G_2 = G_1 + s C_m``.

    Parameters
    ----------
    upstream_admittance : complex or None
        ``G_1``; None for an infinite admittance (returned unchanged).
    s : complex
        Laplace variable in 1/s.
    manifold_capacitance_s : float
        ``C_m = rho_L V_m p_c / (K_L m_dot)``.

    Returns
    -------
    complex or None
        ``G_2``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.1-16)-(5.4.1-17),
           p. 243, https://ntrs.nasa.gov/citations/19720026079
    """
    c = _nonnegative(manifold_capacitance_s, "manifold_capacitance_s")
    if upstream_admittance is None:
        return None
    return complex(upstream_admittance) + s * c


def pressurized_tank_injection_admittance(s: complex, *, injector_impedance: float,
                                          line_impedance: float, line_time_s: float,
                                          manifold_capacitance_s: float,
                                          tank_capacitance_s: float | None = None) -> complex:
    """Closed-form injection admittance of a tank-line-manifold-orifice feed, Eq. (5.4.2-1).

    ``G_i = (1/Z_i) [s Z_i (C_t + C_m) + ((Z_i/Z_l) + s^2 C_m C_t Z_l Z_i) tanh s t_l]
    / [1 + s Z_i (C_t + C_m) + ((Z_i/Z_l) + s C_t Z_l + s^2 C_m C_t Z_l Z_i) tanh s t_l]``

    This is the composition of the component ratios from ``G_0 = 0`` (tank
    capacitance, line, manifold, short orifice). The scan's denominator shows
    the last term as ``s^2 C_m Z_l Z_i``; without ``C_t`` there the form does
    not reduce to the printed Eq. (5.4.2-2) as ``C_t -> inf``, so the composed
    form is used. For ``C_t -> inf`` (``tank_capacitance_s=None``) it is Eq.
    (5.4.2-2), ``G_i = (1/Z_i)(1 + s C_m Z_l tanh s t_l) /
    (1 + [(Z_l/Z_i) + s C_m Z_l] tanh s t_l)``.

    Parameters
    ----------
    s : complex
        Laplace variable in 1/s.
    injector_impedance : float
        ``Z_i = 2 Dp_i m_dot / (m_i p_c)``.
    line_impedance, line_time_s : float
        ``Z_l`` and ``t_l``.
    manifold_capacitance_s : float
        ``C_m``.
    tank_capacitance_s : float, optional
        ``C_t``; None for a constant-pressure tank.

    Returns
    -------
    complex
        ``G_i``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.2-1)-(5.4.2-2),
           p. 244, and Fig. 5.4.2a, p. 245,
           https://ntrs.nasa.gov/citations/19720026079
    """
    zi = _positive(injector_impedance, "injector_impedance")
    zl = _positive(line_impedance, "line_impedance")
    tl = _nonnegative(line_time_s, "line_time_s")
    cm = _nonnegative(manifold_capacitance_s, "manifold_capacitance_s")
    th = cmath.tanh(s * tl)
    if tank_capacitance_s is None:
        return (1.0 / zi) * (1.0 + s * cm * zl * th) / (1.0 + ((zl / zi) + s * cm * zl) * th)
    ct = _positive(tank_capacitance_s, "tank_capacitance_s")
    num = s * zi * (ct + cm) + ((zi / zl) + s * s * cm * ct * zl * zi) * th
    den = 1.0 + s * zi * (ct + cm) + ((zi / zl) + s * ct * zl + s * s * cm * ct * zl * zi) * th
    return (1.0 / zi) * num / den


# ---------------------------------------------------------------------------
# Physical feed description
# ---------------------------------------------------------------------------


class NormalizedFeed(TypedDict):
    """Normalised feed parameters returned by :meth:`FeedSystem.normalized`.

    The keys are the keyword arguments of
    :func:`pressurized_tank_injection_admittance` [1]_.

    Attributes
    ----------
    injector_impedance, line_impedance : float
        ``Z_i`` and ``Z_l``.
    line_time_s, manifold_capacitance_s : float
        ``t_l`` and ``C_m``.
    tank_capacitance_s : float or None
        ``C_t``; None for a constant-pressure tank.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.2-1)-(5.4.2-2),
           p. 244, https://ntrs.nasa.gov/citations/19720026079
    """

    injector_impedance: float
    line_impedance: float
    line_time_s: float
    manifold_capacitance_s: float
    tank_capacitance_s: float | None


@dataclass(frozen=True, slots=True)
class FeedSystem:
    """One propellant's pressurized-tank feed (SP-194 Fig. 5.4.2a), in physical units.

    Attributes
    ----------
    propellant_mass_flow_kg_s : float
        ``m_i`` (or ``m_L``), this propellant's flow.
    injector_pressure_drop_pa : float
        ``Dp_i`` across the injector orifices.
    line_length_m, line_area_m2 : float
        Feed line length ``l`` and flow area ``A_l``.
    liquid_sound_speed_m_s : float
        Effective sound speed ``a_L`` in the line (see
        :func:`wall_compliant_sound_speed`).
    liquid_density_kg_m3 : float
        ``rho_L``.
    liquid_bulk_modulus_pa : float
        ``K_L``.
    manifold_volume_m3 : float
        ``V_m``.
    tank_volume_m3 : float or None
        ``V_t``; None for a constant-pressure tank (``C_t -> inf``).
    vaporization_lag_s : float
        ``tau_v`` of the double-time-lag model (0 for a single lag or a gas).

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Fig. 5.4.2a, p. 245, and
           Eqs. (5.4.1-3), (5.4.1-17), (5.2.3-2),
           https://ntrs.nasa.gov/citations/19720026079
    """

    propellant_mass_flow_kg_s: float
    injector_pressure_drop_pa: float
    line_length_m: float
    line_area_m2: float
    liquid_sound_speed_m_s: float
    liquid_density_kg_m3: float
    liquid_bulk_modulus_pa: float
    manifold_volume_m3: float
    tank_volume_m3: float | None = None
    vaporization_lag_s: float = 0.0

    def normalized(self, *, total_mass_flow_kg_s: float, chamber_pressure_pa: float) -> NormalizedFeed:
        """Return ``Z_i, Z_l, t_l, C_m, C_t`` normalised by the total flow and chamber pressure."""
        m = _positive(total_mass_flow_kg_s, "total_mass_flow_kg_s")
        pc = _positive(chamber_pressure_pa, "chamber_pressure_pa")
        mi = _positive(self.propellant_mass_flow_kg_s, "propellant_mass_flow_kg_s")
        if mi > m:
            raise DomainError("A propellant flow cannot exceed the total flow.",
                              field="propellant_mass_flow_kg_s")
        dp = _positive(self.injector_pressure_drop_pa, "injector_pressure_drop_pa")
        a = _positive(self.liquid_sound_speed_m_s, "liquid_sound_speed_m_s")
        rho = _positive(self.liquid_density_kg_m3, "liquid_density_kg_m3")
        k = _positive(self.liquid_bulk_modulus_pa, "liquid_bulk_modulus_pa")
        tank = None if self.tank_volume_m3 is None else rho * _positive(
            self.tank_volume_m3, "tank_volume_m3") * pc / (k * m)
        return {
            "injector_impedance": 2.0 * dp * m / (mi * pc),
            "line_impedance": a * m / (_positive(self.line_area_m2, "line_area_m2") * pc),
            "line_time_s": _nonnegative(self.line_length_m, "line_length_m") / a,
            "manifold_capacitance_s": rho * _nonnegative(self.manifold_volume_m3, "manifold_volume_m3") * pc / (k * m),
            "tank_capacitance_s": tank,
        }

    def admittance(self, s: complex, *, total_mass_flow_kg_s: float, chamber_pressure_pa: float) -> complex:
        """Injection admittance ``G_i(s)`` from Eq. (5.4.2-1) or (5.4.2-2)."""
        return pressurized_tank_injection_admittance(
            s, **self.normalized(total_mass_flow_kg_s=total_mass_flow_kg_s,
                                 chamber_pressure_pa=chamber_pressure_pa))


@dataclass(frozen=True, slots=True)
class FeedSystemResponse:
    """Injection admittance of one feed at one frequency.

    Attributes
    ----------
    frequency_hz : float
    admittance : complex
        ``G_i(i omega)``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Fig. 5.4.2b, p. 245,
           https://ntrs.nasa.gov/citations/19720026079
    """

    frequency_hz: float
    admittance: complex


@dataclass(frozen=True, slots=True)
class FeedCoupledStability:
    """Nyquist verdict for the feed-coupled chamber.

    Attributes
    ----------
    stable : bool
        True when ``F(s) = 1 + f(s)`` has no zero in the right half-plane.
    right_half_plane_zeros : int
        Count from the clockwise encirclements of the origin (complex
        roots come in conjugate pairs, so an unstable mode counts 2).
    gain_crossover_hz : tuple of float
        Refined gain crossovers bracketed by the resolved contour samples.
        Tangencies and unbracketed crossings are not enumerated.
    phase_margins_deg : tuple of float
        ``arg f(i omega) - (-180 deg)`` at each crossover, wrapped to
        (-180, 180]; positive is stable at a single crossover (SP-194 p. 251).
    minimum_distance_to_minus_one : float
        Conservative lower bound on ``|1 + f(i omega)|`` over the contour
        and its entire omitted tail, conditional on the admittance bounds.
    frequencies_hz : tuple of float
        Evaluated frequencies, including zero. The final frequency starts
        a tail whose global envelope excludes any further encirclement.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Sect. 5.5.1, pp. 250-253,
           https://ntrs.nasa.gov/citations/19720026079
    """

    stable: bool
    right_half_plane_zeros: int
    gain_crossover_hz: tuple[float, ...]
    phase_margins_deg: tuple[float, ...]
    minimum_distance_to_minus_one: float
    frequencies_hz: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class AdmittanceBounds:
    """Caller-supplied bounds for a real, stable custom injection admittance.

    ``maximum_magnitude`` bounds |G(i omega)| for every omega >= 0.
    ``derivative_bound`` is a nonnegative constant, or a function (a,b)
    bounding |dG(i omega)/domega| throughout [a,b], in seconds. The callback
    must also have no RHP poles, be finite on the imaginary axis and obey
    conjugate symmetry. These are caller obligations, not facts inferred
    from a finite set of samples. Built-in physical feeds derive their
    bounds automatically. A constant rigid admittance uses a zero slope.
    """

    maximum_magnitude: float
    derivative_bound: float | Callable[[float, float], float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "maximum_magnitude", _nonnegative(self.maximum_magnitude, "maximum_magnitude"))
        if not callable(self.derivative_bound):
            object.__setattr__(self, "derivative_bound", _nonnegative(self.derivative_bound, "derivative_bound"))

    def slope(self, lower: float, upper: float) -> float:
        value = self.derivative_bound(lower, upper) if callable(self.derivative_bound) else self.derivative_bound
        return _nonnegative(value, "derivative_bound")


def _feed_admittance_bounds(feed: NormalizedFeed) -> AdmittanceBounds:
    """Bounds for the passive lossless line/capacitance plus resistive injector.

    On the imaginary axis G = A/(Zi*A +/- i*B), with real A,B, hence
    |G| <= 1/Zi and |dG/domega| = |A'B-AB'|/|D|^2. The line transfer
    matrix's determinant bounds |D| away from zero; a local midpoint
    Lipschitz bound sharpens it as intervals shrink. No resonance is
    excluded merely because it lies between frequency samples.
    """
    zi, zl = feed["injector_impedance"], feed["line_impedance"]
    tl, cm, ct = feed["line_time_s"], feed["manifold_capacitance_s"], feed["tank_capacitance_s"]

    def derivative(lower: float, upper: float) -> float:
        mid = lower + (upper - lower) / 2
        c, s = cos(mid * tl), sin(mid * tl)
        if ct is None:
            a, b = c - mid * cm * zl * s, zl * s
            frobenius = hypot(zi, zi * upper * cm * zl, zl)
            determinant = zi * zl
            a_prime = tl + cm * zl + upper * cm * zl * tl
            b_prime = zl * tl
            numerator = zl * (tl + cm * zl)
        else:
            coefficient = max(abs(1 / zl - lower**2 * cm * ct * zl),
                              abs(1 / zl - upper**2 * cm * ct * zl))
            a = mid * (ct + cm) * c + (1 / zl - mid**2 * cm * ct * zl) * s
            b = c - mid * ct * zl * s
            frobenius = hypot(zi * upper * (ct + cm), zi * coefficient, 1, upper * ct * zl)
            determinant = zi * (1 / zl + lower**2 * ct**2 * zl)
            a_prime = ct + cm + upper * (ct + cm) * tl + 2 * upper * cm * ct * zl + coefficient * tl
            b_prime = tl + ct * zl + upper * ct * zl * tl
            # Exact Wronskian: tl*(1/Zl + omega^2*Ct^2*Zl) + Ct + Cm*B^2.
            numerator = tl * (1 / zl + upper**2 * ct**2 * zl) + ct + cm * (1 + upper * ct * zl)**2
        denominator = max(determinant / frobenius,
                          hypot(zi * a, b) - hypot(zi * a_prime, b_prime) * ((upper - lower) / 2))
        return numerator / denominator**2 * (1 + 1e-12)

    return AdmittanceBounds(1 / zi, derivative)


def _open_loop(s: complex, *, k_ox: float, k_f: float, g_ox, g_f, tau_m: float,
               tau_vox: float, tau_vf: float, theta: float) -> complex:
    return cmath.exp(-s * tau_m) / (1.0 + s * theta) * (
        k_ox * g_ox(s) * cmath.exp(-s * tau_vox) + k_f * g_f(s) * cmath.exp(-s * tau_vf))


@numeric_range_guard("feed_coupled_stability")
def feed_coupled_stability(*, mixture_ratio: float, characteristic_velocity_m_s: float,
                           cstar_mixture_ratio_derivative_m_s: float, chamber_pressure_pa: float,
                           gas_residence_time_s: float, combustion_time_lag_s: float,
                           oxidizer_feed: FeedSystem | None = None, fuel_feed: FeedSystem | None = None,
                           oxidizer_admittance=None, fuel_admittance=None,
                           oxidizer_admittance_bounds: AdmittanceBounds | None = None,
                           fuel_admittance_bounds: AdmittanceBounds | None = None,
                           maximum_frequency_hz: float | None = None,
                           points: int = 40000) -> FeedCoupledStability:
    """Nyquist stability of the single- or double-time-lag chamber with dynamic feeds.

    Each propellant's admittance comes from a :class:`FeedSystem` (Eq.
    5.4.2-1/2) or from a callable ``G(s)`` (any chain of the component
    ratios, or a measured admittance). ``combustion_time_lag_s`` is ``tau_m``
    (or the single lag ``tau_T`` when both vaporization lags are 0).

    Parameters
    ----------
    mixture_ratio, characteristic_velocity_m_s, cstar_mixture_ratio_derivative_m_s : float
        ``r``, ``c*`` and ``dc*/dr``.
    chamber_pressure_pa : float
        ``p_c``.
    gas_residence_time_s : float
        ``theta_c``.
    combustion_time_lag_s : float
        ``tau_m`` (common mixing-reaction lag).
    oxidizer_feed, fuel_feed : FeedSystem, optional
        Physical feeds; the total flow is the sum of the two propellant flows.
    oxidizer_admittance, fuel_admittance : callable, optional
        ``G(s)`` used instead of a feed (both needed when no feeds are given).
    oxidizer_admittance_bounds, fuel_admittance_bounds : AdmittanceBounds, optional
        Required for custom admittances. Global amplitude and interval
        derivative bounds cannot be established from sampling a callback.
    maximum_frequency_hz : float, optional
        Maximum frequency available to close the contour. Must place its
        entire remaining tail in |f| < 1 by the global amplitude envelope.
        Defaults to at least 50 / theta_c; the solve stops earlier when
        that envelope is <= 1/2.
    points : int, optional
        Maximum response evaluations, 1000 to 1000000 (default 40000).
        Exhaustion raises ``nyquist_resolution`` rather than a verdict.

    Returns
    -------
    FeedCoupledStability

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.2.2-6), (5.2.3-3),
           (5.4.2-1)-(5.4.2-2) and Sect. 5.5.1, pp. 236-253,
           https://ntrs.nasa.gov/citations/19720026079
    """
    r = _positive(mixture_ratio, "mixture_ratio")
    cstar = _positive(characteristic_velocity_m_s, "characteristic_velocity_m_s")
    slope = real(cstar_mixture_ratio_derivative_m_s, "cstar_mixture_ratio_derivative_m_s")
    pc = _positive(chamber_pressure_pa, "chamber_pressure_pa")
    theta = _positive(gas_residence_time_s, "gas_residence_time_s")
    tau_m = _nonnegative(combustion_time_lag_s, "combustion_time_lag_s")
    k_ox = 1.0 + (1.0 + r) / cstar * slope
    k_f = 1.0 - r * (1.0 + r) / cstar * slope
    tau_vox = tau_vf = 0.0
    if oxidizer_feed is not None or fuel_feed is not None:
        if any(value is not None for value in (oxidizer_admittance, fuel_admittance,
                                               oxidizer_admittance_bounds, fuel_admittance_bounds)):
            raise DomainError("Give physical feeds or bounded callbacks, not both.", field="oxidizer_feed")
        if oxidizer_feed is None or fuel_feed is None:
            raise DomainError("Give both feeds or neither.", field="fuel_feed")
        total = oxidizer_feed.propellant_mass_flow_kg_s + fuel_feed.propellant_mass_flow_kg_s
        ratio = oxidizer_feed.propellant_mass_flow_kg_s / fuel_feed.propellant_mass_flow_kg_s
        if abs(ratio / r - 1.0) > 1e-6:
            raise DomainError("The feed flows do not match the mixture ratio.", field="mixture_ratio",
                              details={"feed_ratio": ratio, "mixture_ratio": r})
        tau_vox = _nonnegative(oxidizer_feed.vaporization_lag_s, "vaporization_lag_s")
        tau_vf = _nonnegative(fuel_feed.vaporization_lag_s, "vaporization_lag_s")
        ox_norm = oxidizer_feed.normalized(total_mass_flow_kg_s=total, chamber_pressure_pa=pc)
        fu_norm = fuel_feed.normalized(total_mass_flow_kg_s=total, chamber_pressure_pa=pc)
        ox_bounds, fu_bounds = _feed_admittance_bounds(ox_norm), _feed_admittance_bounds(fu_norm)

        def g_ox(s):
            return pressurized_tank_injection_admittance(s, **ox_norm)

        def g_f(s):
            return pressurized_tank_injection_admittance(s, **fu_norm)
    else:
        if oxidizer_admittance is None or fuel_admittance is None:
            raise DomainError("Give feeds or both admittance functions.", field="oxidizer_admittance")
        if not (callable(oxidizer_admittance) and callable(fuel_admittance)):
            raise DomainError("Admittances must be callable.", field="oxidizer_admittance")
        if not (isinstance(oxidizer_admittance_bounds, AdmittanceBounds)
                and isinstance(fuel_admittance_bounds, AdmittanceBounds)):
            raise DomainError("Custom admittances require global amplitude and derivative bounds.",
                              code="admittance_bounds_required", field="oxidizer_admittance_bounds")
        ox_bounds, fu_bounds = oxidizer_admittance_bounds, fuel_admittance_bounds
        g_ox, g_f = oxidizer_admittance, fuel_admittance
    f_max = None if maximum_frequency_hz is None else _positive(maximum_frequency_hz, "maximum_frequency_hz")

    def open_loop(s: complex) -> complex:
        return _open_loop(s, k_ox=k_ox, k_f=k_f, g_ox=g_ox, g_f=g_f, tau_m=tau_m,
                          tau_vox=tau_vox, tau_vf=tau_vf, theta=theta)

    amplitude = _nonnegative(abs(k_ox) * ox_bounds.maximum_magnitude + abs(k_f) * fu_bounds.maximum_magnitude,
                             "admittance_amplitude_bound")

    def slope_bound(lower: float, upper: float) -> float:
        denominator = hypot(1, lower * theta)
        return sum(abs(k) * ((bounds.slope(lower, upper) + (tau_m + lag) * bounds.maximum_magnitude) / denominator
                             + theta * bounds.maximum_magnitude / denominator**2)
                   for k, bounds, lag in ((k_ox, ox_bounds, tau_vox), (k_f, fu_bounds, tau_vf)))

    zeros, crossings, margins, distance, freqs = resolve_contour(
        response=open_loop, slope_bound=slope_bound, amplitude_bound=amplitude,
        residence_time_s=theta, maximum_frequency_hz=f_max, maximum_evaluations=points)
    return FeedCoupledStability(zeros == 0, zeros, crossings, margins, distance, freqs)


def feed_admittance_sweep(feed: FeedSystem, *, total_mass_flow_kg_s: float,
                          chamber_pressure_pa: float,
                          frequencies_hz: Sequence[float]) -> tuple[FeedSystemResponse, ...]:
    """``G_i(i omega)`` of a feed at the given frequencies (the curves of SP-194 Fig. 5.4.2b).

    Parameters
    ----------
    feed : FeedSystem
    total_mass_flow_kg_s, chamber_pressure_pa : float
        Normalisation of Eq. (5.4.2-1).
    frequencies_hz : sequence of float

    Returns
    -------
    tuple of FeedSystemResponse

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.4.2-1)-(5.4.2-3) and
           Fig. 5.4.2b, pp. 244-245, https://ntrs.nasa.gov/citations/19720026079
    """
    out = []
    for f in frequencies_hz:
        freq = _positive(f, "frequencies_hz")
        out.append(FeedSystemResponse(freq, feed.admittance(
            2j * pi * freq, total_mass_flow_kg_s=total_mass_flow_kg_s,
            chamber_pressure_pa=chamber_pressure_pa)))
    return tuple(out)


__all__.append("feed_admittance_sweep")
