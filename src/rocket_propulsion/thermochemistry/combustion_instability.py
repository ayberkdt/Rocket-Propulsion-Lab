"""Combustion instability: chamber acoustics, time-lag stability limits, damping.

The relations are those of NASA SP-194, *Liquid Propellant Rocket Combustion
Instability* (Harrje and Reardon, eds., 1972) [1]_, read from the scanned
report:

* **Chamber acoustics.** Transverse eigenvalues ``s_nu,eta``, the roots of
  ``J_nu'(s) = 0`` (Table 3.5.1), and the frequency of a combined
  tangential, radial and longitudinal mode of a closed cylindrical chamber
  without mean flow, ``omega = sqrt(s^2 + (j pi r_c/L)^2)`` (Eq. 3.5.1-19,
  Eq. 8.2.1-11), with the first-longitudinal mean-flow factor
  ``(1 - M^2)`` of Sect. 4.2.2.3.
* **Low-frequency (chugging) stability.** The single-time-lag model with
  purely resistive injectors and decoupled feed lines: gain ``K`` and the
  neutral-stability time lag of Eq. (5.5.1-6), Example 5.5.1.1.
* **Intrinsic instability.** Crocco and Cheng's critical sensitive time lag
  for ``G_ox = G_F = 0``, Eqs. (5.3.1-18)-(5.3.1-19).
* **High-frequency n-tau neutral curve.** ``n (1 - cos omega tau) = G_R``,
  Eq. (4.2.2-23), with minimum ``n = G_R/2`` (Fig. 4.2.2a).
* **Damping devices and measures.** Helmholtz-resonator frequency
  (Eq. 8.3.1-2), effective orifice length (Eqs. 8.3.2-1 and 8.3.2-2), the
  liner-array form (Eq. 8.3.3-1), the response factor (Eqs. 4.4.1-5 and the
  integral form on p. 210), the quality factor and the decay rate of a
  ringing mode (Sect. 9.6.6).

Conventions: SI units. Frequencies in Hz unless named ``angular``. The
transverse index ``eta`` follows Table 3.5.1, where the first root of
``J_0'`` (``s = 0``) is ``eta = 1``, so the first radial mode is
``(nu, eta) = (0, 2)``.

Limits: the chamber is a rigid-walled circular cylinder; the nozzle end is
treated as a closed (velocity-node) end, the lowest-order approximation in
SP-194 Sect. 3.5.1. Mean-flow corrections other than the 1L factor need the
nozzle admittance and the axial profiles of SP-194 Eqs. (3.5.1-20)-(3.5.1-22)
and are not modelled; transverse frequencies with a nonzero Mach number carry
a warning. Baffles, liners, annular chambers and nonlinear (triggered)
instability are not modelled. SP-194 misprints the 1T2R root as 8.5263; the
computed root, 8.5363, is used.

References
----------
.. [1] NASA-SP-194-1972. D. T. Harrje and F. H. Reardon (eds.), Liquid
       Propellant Rocket Combustion Instability, NASA SP-194, 1972,
       https://ntrs.nasa.gov/citations/19720026079
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from math import acos, atan, cos, isfinite, log, pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import numeric_range_guard, real
from rocket_propulsion.heat_transfer.conduction import bessel_jn

__all__ = [
    "AcousticMode",
    "ChuggingStability",
    "acoustic_mode_frequency",
    "chamber_acoustic_modes",
    "chugging_stability",
    "damping_rate_from_decay",
    "helmholtz_resonator_frequency",
    "intrinsic_critical_time_lag",
    "liner_resonant_frequency",
    "neutral_interaction_index",
    "orifice_effective_length",
    "quality_factor",
    "response_factor",
    "transverse_eigenvalue",
]

_SP194 = "NASA SP-194 (1972)"

#: Printed Table 3.5.1 of SP-194 (p. 143): ``(nu, eta) -> s_nu,eta``. Kept as
#: printed for the tests; the 1T2R entry is a misprint of 8.5363.
SP194_TABLE_3_5_1: dict[tuple[int, int], float] = {
    (1, 1): 1.8413, (2, 1): 3.0543, (0, 2): 3.8317, (3, 1): 4.2012,
    (0, 3): 7.0156, (1, 2): 5.3313, (1, 3): 8.5263, (2, 2): 6.7060,
}


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


def _index(value: object, label: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise DomainError(f"{label} must be an integer >= {minimum}.", field=label)
    return value


# ---------------------------------------------------------------------------
# Chamber acoustics
# ---------------------------------------------------------------------------


def _bessel_jn_derivative(order: int, x: float) -> float:
    """``J_n'(x)`` from DLMF 10.6.1 (``J_0' = -J_1``)."""
    if order == 0:
        return -bessel_jn(1, x)
    return 0.5 * (bessel_jn(order - 1, x) - bessel_jn(order + 1, x))


def transverse_eigenvalue(tangential_order: int, radial_index: int) -> float:
    """Transverse eigenvalue ``s_nu,eta``: the ``eta``-th root of ``J_nu'(s) = 0``.

    Parameters
    ----------
    tangential_order : int
        ``nu >= 0``, the number of nodal diameters.
    radial_index : int
        ``eta >= 1``; for ``nu = 0`` the root ``s = 0`` counts as ``eta = 1``
        (Table 3.5.1 numbering), so the first radial mode is ``(0, 2)``.

    Returns
    -------
    float
        ``s_nu,eta``; the roots are bracketed on a 0.05 grid and bisected to
        an interval of 1e-13. Transverse frequency
        ``f = s a / (2 pi r_c)``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Table 3.5.1, p. 143,
           https://ntrs.nasa.gov/citations/19720026079
    .. [2] DLMF-10-BESSEL. NIST DLMF Eqs. 10.6.1 and 10.9.2,
           https://dlmf.nist.gov/10.6
    """
    nu = _index(tangential_order, "tangential_order", 0)
    eta = _index(radial_index, "radial_index", 1)
    if nu > 30 or eta > 30:
        raise DomainError("Orders above 30 are outside the supported range.",
                          field="tangential_order")
    return _eigenvalue(nu, eta)


@lru_cache(maxsize=1024)
def _eigenvalue(nu: int, eta: int) -> float:
    """Cached root search behind :func:`transverse_eigenvalue` (validated arguments)."""
    wanted = eta - 1 if nu == 0 else eta
    if wanted == 0:
        return 0.0
    step = 0.05
    found = 0
    x0 = step
    f0 = _bessel_jn_derivative(nu, x0)
    while True:
        x1 = x0 + step
        f1 = _bessel_jn_derivative(nu, x1)
        if f0 == 0.0 or f0 * f1 < 0.0:
            found += 1
            if found == wanted:
                lo, hi, flo = x0, x1, f0
                if flo == 0.0:
                    return lo
                while hi - lo > 1e-13:
                    mid = 0.5 * (lo + hi)
                    fm = _bessel_jn_derivative(nu, mid)
                    if (fm < 0.0) == (flo < 0.0):
                        lo, flo = mid, fm
                    else:
                        hi = mid
                return 0.5 * (lo + hi)
        x0, f0 = x1, f1


@dataclass(frozen=True, slots=True)
class AcousticMode:
    """One acoustic mode of a closed cylindrical chamber.

    Attributes
    ----------
    name : str
        Mode label such as ``"1T"``, ``"1R"``, ``"2L"`` or ``"1T1L"``.
    tangential_order, radial_order, longitudinal_order : int
        ``nu``, the radial order ``eta - 1`` and ``j``.
    eigenvalue : float
        ``s_nu,eta``.
    frequency_hz : float
        Mode frequency.
    warnings : tuple of WarningMessage
        Mean-flow notes.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (3.5.1-19), p. 144, and
           Table 3.5.1, p. 143, https://ntrs.nasa.gov/citations/19720026079
    """

    name: str
    tangential_order: int
    radial_order: int
    longitudinal_order: int
    eigenvalue: float
    frequency_hz: float
    warnings: tuple[WarningMessage, ...] = ()


def _mode_name(nu: int, radial: int, j: int) -> str:
    parts = []
    if nu:
        parts.append(f"{nu}T")
    if radial:
        parts.append(f"{radial}R")
    if j:
        parts.append(f"{j}L")
    return "".join(parts) or "0"


@numeric_range_guard("chamber_length_m")
def acoustic_mode_frequency(*, sound_speed_m_s: float, chamber_radius_m: float,
                            chamber_length_m: float, tangential_order: int = 0,
                            radial_index: int = 1, longitudinal_order: int = 0,
                            mean_mach_number: float = 0.0) -> AcousticMode:
    """Frequency of an acoustic mode of a closed cylindrical chamber.

    ``f = (a / 2 pi) sqrt((s_nu,eta / r_c)^2 + (j pi / L)^2)``, the
    dimensional form of SP-194 Eq. (3.5.1-19) (lengths there normalised by
    ``L``; Eq. 8.2.1-11 by ``r_c``). For a pure longitudinal mode with a
    mean flow the frequency is multiplied by ``(1 - M^2)`` (SP-194 p. 182,
    there for the 1L mode with ``L`` the effective length including part of
    the subsonic nozzle).

    Parameters
    ----------
    sound_speed_m_s : float
        Mean speed of sound of the chamber gas.
    chamber_radius_m : float
        Chamber radius ``r_c``.
    chamber_length_m : float
        Acoustic (effective) chamber length ``L``.
    tangential_order, radial_index : int, optional
        ``nu`` and ``eta`` as in :func:`transverse_eigenvalue`.
    longitudinal_order : int, optional
        ``j >= 0``.
    mean_mach_number : float, optional
        Mean chamber Mach number, in [0, 1).

    Returns
    -------
    AcousticMode

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (3.5.1-19), p. 144;
           Table 3.5.1, p. 143; Sect. 4.2.2.3, p. 182,
           https://ntrs.nasa.gov/citations/19720026079
    """
    a = _positive(sound_speed_m_s, "sound_speed_m_s")
    radius = _positive(chamber_radius_m, "chamber_radius_m")
    length = _positive(chamber_length_m, "chamber_length_m")
    j = _index(longitudinal_order, "longitudinal_order", 0)
    mach = _nonnegative(mean_mach_number, "mean_mach_number")
    if mach >= 1.0:
        raise DomainError("mean_mach_number must be below 1.", field="mean_mach_number")
    s = transverse_eigenvalue(tangential_order, radial_index)
    nu = tangential_order
    radial = radial_index - 1
    if s == 0.0 and j == 0:
        raise DomainError("The (0, 1, 0) mode is the trivial mean state.",
                          field="longitudinal_order")
    frequency = a / (2.0 * pi) * sqrt((s / radius) ** 2 + (j * pi / length) ** 2)
    warnings: list[WarningMessage] = []
    if mach > 0.0:
        if s == 0.0:
            frequency *= 1.0 - mach * mach
        else:
            warnings.append(WarningMessage(
                "mach_correction_not_modelled",
                "SP-194 gives the (1 - M^2) factor only for longitudinal modes; the "
                "transverse and mixed frequency is the no-flow value.", "warning",
                "mean_mach_number"))
    return AcousticMode(_mode_name(nu, radial, j), nu, radial, j, s, frequency, tuple(warnings))


def chamber_acoustic_modes(*, sound_speed_m_s: float, chamber_radius_m: float,
                           chamber_length_m: float, maximum_frequency_hz: float,
                           mean_mach_number: float = 0.0,
                           max_order: int = 6) -> tuple[AcousticMode, ...]:
    """All modes up to ``maximum_frequency_hz``, sorted by frequency.

    Parameters
    ----------
    sound_speed_m_s, chamber_radius_m, chamber_length_m, mean_mach_number : float
        As in :func:`acoustic_mode_frequency`.
    maximum_frequency_hz : float
        Upper frequency bound.
    max_order : int, optional
        Largest ``nu``, radial order and ``j`` searched.

    Returns
    -------
    tuple of AcousticMode

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (3.5.1-19), p. 144,
           https://ntrs.nasa.gov/citations/19720026079
    """
    limit = _positive(maximum_frequency_hz, "maximum_frequency_hz")
    top = _index(max_order, "max_order", 1)
    modes = []
    for nu in range(top + 1):
        for eta in range(1, top + 2):
            for j in range(top + 1):
                if nu == 0 and eta == 1 and j == 0:
                    continue
                mode = acoustic_mode_frequency(
                    sound_speed_m_s=sound_speed_m_s, chamber_radius_m=chamber_radius_m,
                    chamber_length_m=chamber_length_m, tangential_order=nu,
                    radial_index=eta, longitudinal_order=j, mean_mach_number=mean_mach_number)
                if mode.frequency_hz <= limit:
                    modes.append(mode)
    return tuple(sorted(modes, key=lambda m: m.frequency_hz))


# ---------------------------------------------------------------------------
# Time-lag stability
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ChuggingStability:
    """Single-time-lag low-frequency stability with resistive injectors.

    Attributes
    ----------
    gain : float
        ``K = K_ox + K_F`` of SP-194 Example 5.5.1.1.
    oxidizer_gain, fuel_gain : float
        The oxidizer and fuel terms of ``K``.
    unconditionally_stable : bool
        True when ``K <= 1``.
    critical_time_lag_s : float or None
        Largest stable total time lag ``tau_T`` from Eq. (5.5.1-6); None when
        unconditionally stable.
    neutral_frequency_hz : float or None
        ``sqrt(K^2 - 1)/(2 pi theta_c)``, where ``|f(i omega)| = 1``; derived
        here from ``f(s) = K exp(-s tau)/(1 + s theta_c)``.
    stable : bool or None
        Verdict for the supplied ``time_lag_s``; None when no lag was given.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Sect. 5.5.1, Eq. (5.5.1-6)
           and Example 5.5.1.1, p. 253,
           https://ntrs.nasa.gov/citations/19720026079
    """

    gain: float
    oxidizer_gain: float
    fuel_gain: float
    unconditionally_stable: bool
    critical_time_lag_s: float | None
    neutral_frequency_hz: float | None
    stable: bool | None


def chugging_stability(*, mixture_ratio: float, characteristic_velocity_m_s: float,
                       cstar_mixture_ratio_derivative_m_s: float,
                       chamber_pressure_pa: float, oxidizer_injector_drop_pa: float,
                       fuel_injector_drop_pa: float, gas_residence_time_s: float,
                       time_lag_s: float | None = None) -> ChuggingStability:
    """Chugging limit of the single-time-lag model with rigid, resistive feeds.

    ``f(s) = K exp(-s tau_T)/(1 + s theta_c)`` with

        ``K = [1 + ((1+r)/c*) dc*/dr] (r/(1+r)) p_c/(2 Dp_ox)
        + [1 - (r(1+r)/c*) dc*/dr] (1/(1+r)) p_c/(2 Dp_F)``.

    ``K < 1`` is unconditionally stable; otherwise the engine is stable for
    ``tau_T/theta_c < (pi - atan sqrt(K^2-1))/sqrt(K^2-1)`` (Eq. 5.5.1-6).

    Parameters
    ----------
    mixture_ratio : float
        ``r = O/F``.
    characteristic_velocity_m_s : float
        ``c*``.
    cstar_mixture_ratio_derivative_m_s : float
        ``dc*/dr`` at the operating point.
    chamber_pressure_pa : float
        ``p_c``.
    oxidizer_injector_drop_pa, fuel_injector_drop_pa : float
        Injector pressure drops.
    gas_residence_time_s : float
        ``theta_c = rho V / m_dot``.
    time_lag_s : float, optional
        Total combustion time lag ``tau_T`` to classify.

    Returns
    -------
    ChuggingStability

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.2.2-6) and
           (5.5.1-6), Example 5.5.1.1, pp. 236-237 and 253,
           https://ntrs.nasa.gov/citations/19720026079
    """
    r = _positive(mixture_ratio, "mixture_ratio")
    cstar = _positive(characteristic_velocity_m_s, "characteristic_velocity_m_s")
    slope = real(cstar_mixture_ratio_derivative_m_s, "cstar_mixture_ratio_derivative_m_s")
    pc = _positive(chamber_pressure_pa, "chamber_pressure_pa")
    dp_ox = _positive(oxidizer_injector_drop_pa, "oxidizer_injector_drop_pa")
    dp_f = _positive(fuel_injector_drop_pa, "fuel_injector_drop_pa")
    theta = _positive(gas_residence_time_s, "gas_residence_time_s")
    k_ox = (1.0 + (1.0 + r) / cstar * slope) * r / (1.0 + r) * pc / (2.0 * dp_ox)
    k_f = (1.0 - r * (1.0 + r) / cstar * slope) / (1.0 + r) * pc / (2.0 * dp_f)
    k = k_ox + k_f
    if k <= 0.0:
        raise DomainError("The injection gain K is not positive; the single-lag model "
                          "does not apply.", field="cstar_mixture_ratio_derivative_m_s")
    lag = None if time_lag_s is None else _nonnegative(time_lag_s, "time_lag_s")
    if k <= 1.0:
        return ChuggingStability(k, k_ox, k_f, True, None, None, None if lag is None else True)
    root = sqrt(k * k - 1.0)
    critical = theta * (pi - atan(root)) / root
    return ChuggingStability(k, k_ox, k_f, False, critical, root / (2.0 * pi * theta),
                             None if lag is None else lag < critical)


def intrinsic_critical_time_lag(interaction_index: float, gas_residence_time_s: float) -> float | None:
    """Crocco-Cheng critical sensitive time lag for intrinsic low-frequency instability.

    With rigid feeds (``G_ox = G_F = 0``) the characteristic equation is
    ``n (1 - exp(-s tau)) = 1 + s theta_c`` (Eq. 5.3.1-18). Instability is
    possible only for ``n > 1/2``, and then for ``tau > tau_cr`` with

        ``tau_cr = theta_c [pi - arccos((1 - n)/n)] / sqrt(2n - 1)``

    (Eq. 5.3.1-19; the scan prints the argument as "1-n/n"). The neutral
    angular frequency is ``sqrt(2n - 1)/theta_c``.

    Parameters
    ----------
    interaction_index : float
        Crocco's ``n``.
    gas_residence_time_s : float
        ``theta_c``.

    Returns
    -------
    float or None
        ``tau_cr`` in s; None when ``n <= 1/2`` (no intrinsic instability).

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (5.3.1-17)-(5.3.1-19),
           p. 239, https://ntrs.nasa.gov/citations/19720026079
    """
    n = _positive(interaction_index, "interaction_index")
    theta = _positive(gas_residence_time_s, "gas_residence_time_s")
    if n <= 0.5:
        return None
    return theta * (pi - acos((1.0 - n) / n)) / sqrt(2.0 * n - 1.0)


def neutral_interaction_index(phase_lag_rad: float, stability_factor: float) -> float:
    """Neutral interaction index of the linear n-tau model, ``n = G_R / (1 - cos omega tau)``.

    SP-194 Eq. (4.2.2-23): ``n (1 - cos omega_0 tau) = G_R``; the minimum,
    ``n = G_R/2``, lies at ``omega_0 tau = pi`` (``tau/t_w = 0.5``,
    Fig. 4.2.2a) and repeats every ``2 pi``.

    Parameters
    ----------
    phase_lag_rad : float
        ``omega_0 tau``; must not be a multiple of ``2 pi``.
    stability_factor : float
        ``G_R`` (Eq. 4.2.2-23a for longitudinal and Eq. 4.2.2-29 for
        transverse modes), which collects the nozzle, mean-flow and drag
        damping and must be supplied.

    Returns
    -------
    float
        ``n`` on the stability limit; the mode is unstable for larger ``n``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (4.2.2-23) and Fig.
           4.2.2a, pp. 180-182, https://ntrs.nasa.gov/citations/19720026079
    """
    phase = real(phase_lag_rad, "phase_lag_rad")
    g_r = _positive(stability_factor, "stability_factor")
    denominator = 1.0 - cos(phase)
    if denominator < 1e-12:
        raise DomainError("omega tau is a multiple of 2 pi; no finite n reaches the limit.",
                          field="phase_lag_rad")
    return g_r / denominator


# ---------------------------------------------------------------------------
# Resonators and damping measures
# ---------------------------------------------------------------------------


def orifice_effective_length(*, orifice_length_m: float, orifice_diameter_m: float,
                             open_area_ratio: float,
                             regime: str = "linear") -> float:
    """Effective orifice length ``l_eff = l_o + delta`` of a resonator aperture.

    ``"linear"``: ``delta = 0.85 d_o (1 - 0.7 sqrt(a))`` (Eq. 8.3.2-1).
    ``"high_amplitude"``: ``l_o + 0.375 delta`` for high mean flow or high
    wave amplitude (Eq. 8.3.2-2).

    Parameters
    ----------
    orifice_length_m : float
        Aperture thickness ``l_o``.
    orifice_diameter_m : float
        ``d_o``.
    open_area_ratio : float
        ``a``, aperture area over resonator cross-section area, as a fraction
        in (0, 1].
    regime : {"linear", "high_amplitude"}, optional

    Returns
    -------
    float
        ``l_eff`` in m.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (8.3.2-1) and
           (8.3.2-2), pp. 410-411, https://ntrs.nasa.gov/citations/19720026079
    """
    length = _nonnegative(orifice_length_m, "orifice_length_m")
    diameter = _positive(orifice_diameter_m, "orifice_diameter_m")
    ratio = _positive(open_area_ratio, "open_area_ratio")
    if ratio > 1.0:
        raise DomainError("open_area_ratio is a fraction and must not exceed 1.",
                          field="open_area_ratio")
    if regime not in ("linear", "high_amplitude"):
        raise DomainError("regime must be 'linear' or 'high_amplitude'.", field="regime")
    delta = 0.85 * diameter * (1.0 - 0.7 * sqrt(ratio))
    return length + (delta if regime == "linear" else 0.375 * delta)


def helmholtz_resonator_frequency(*, sound_speed_m_s: float, orifice_area_m2: float,
                                  effective_length_m: float, cavity_volume_m3: float) -> float:
    """Helmholtz resonator frequency, ``omega_r = a sqrt(A_o/(l_eff V_cav))``.

    Valid when the largest cavity dimension is below a quarter wavelength
    (lumped-element assumption, SP-194 p. 401).

    Parameters
    ----------
    sound_speed_m_s : float
        Speed of sound of the gas in the resonator.
    orifice_area_m2 : float
        ``A_o``.
    effective_length_m : float
        ``l_eff``, see :func:`orifice_effective_length`.
    cavity_volume_m3 : float
        ``V_cav``.

    Returns
    -------
    float
        ``f_r = omega_r / (2 pi)`` in Hz.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (8.3.1-2), p. 401,
           https://ntrs.nasa.gov/citations/19720026079
    """
    a = _positive(sound_speed_m_s, "sound_speed_m_s")
    area = _positive(orifice_area_m2, "orifice_area_m2")
    length = _positive(effective_length_m, "effective_length_m")
    volume = _positive(cavity_volume_m3, "cavity_volume_m3")
    return a * sqrt(area / (length * volume)) / (2.0 * pi)


def liner_resonant_frequency(*, sound_speed_m_s: float, open_area_ratio: float,
                             cavity_depth_m: float, effective_length_m: float) -> float:
    """Resonant frequency of a large flat liner array, ``omega_r = a sqrt(a/(L_cav l_eff))``.

    Parameters
    ----------
    sound_speed_m_s : float
        Speed of sound of the cavity gas, ``a_3``.
    open_area_ratio : float
        Open-area fraction ``a`` in (0, 1].
    cavity_depth_m : float
        Backing-cavity depth ``L_cav``.
    effective_length_m : float
        ``l_eff``.

    Returns
    -------
    float
        ``f_r`` in Hz.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eq. (8.3.3-1), p. 415,
           https://ntrs.nasa.gov/citations/19720026079
    """
    a = _positive(sound_speed_m_s, "sound_speed_m_s")
    ratio = _positive(open_area_ratio, "open_area_ratio")
    if ratio > 1.0:
        raise DomainError("open_area_ratio must not exceed 1.", field="open_area_ratio")
    depth = _positive(cavity_depth_m, "cavity_depth_m")
    length = _positive(effective_length_m, "effective_length_m")
    return a * sqrt(ratio / (depth * length)) / (2.0 * pi)


def response_factor(pressure_fluctuation: Sequence[float],
                    source_fluctuation: Sequence[float]) -> float:
    """Response factor ``N = sum(m' p') / sum(p'^2)`` of sampled fluctuations.

    The discrete form of SP-194's ``N = int int m' p_c' dt dV /
    int int (p_c')^2 dt dV`` (p. 210) for equally spaced samples over whole
    cycles at one location; for sinusoids it equals
    ``(m°'/p_c°') cos phi`` (Eq. 4.4.1-5). ``N > 0`` means the source is in
    phase with the pressure and drives the oscillation (Rayleigh's
    criterion, p. 208).

    Parameters
    ----------
    pressure_fluctuation : sequence of float
        ``p'`` samples, normalised as the caller chooses (e.g. ``p'/p_bar``).
    source_fluctuation : sequence of float
        ``m'`` (or heat-release) samples on the same normalisation.

    Returns
    -------
    float
        ``N``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Eqs. (4.4.1-4)-(4.4.1-5) and
           the integral form, pp. 208-210,
           https://ntrs.nasa.gov/citations/19720026079
    """
    p = [real(v, "pressure_fluctuation") for v in pressure_fluctuation]
    m = [real(v, "source_fluctuation") for v in source_fluctuation]
    if len(p) != len(m) or len(p) < 2:
        raise DomainError("Both series need the same length of at least 2.",
                          field="source_fluctuation")
    energy = sum(v * v for v in p)
    if energy == 0.0:
        raise DomainError("The pressure fluctuation is identically zero.",
                          field="pressure_fluctuation")
    return sum(a * b for a, b in zip(p, m)) / energy


def quality_factor(resonant_frequency_hz: float, lower_half_power_hz: float,
                   upper_half_power_hz: float) -> float:
    """Quality factor ``Q = f_0 / (f_2 - f_1)`` from the half-power bandwidth.

    ``f_1`` and ``f_2`` are where the response amplitude is ``0.707`` of the
    peak (-3 dB).

    Parameters
    ----------
    resonant_frequency_hz, lower_half_power_hz, upper_half_power_hz : float
        ``f_0``, ``f_1`` and ``f_2`` with ``f_1 < f_0 < f_2``.

    Returns
    -------
    float
        ``Q``.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Sect. 9.6.6.2, p. 523,
           https://ntrs.nasa.gov/citations/19720026079
    """
    f0 = _positive(resonant_frequency_hz, "resonant_frequency_hz")
    f1 = _positive(lower_half_power_hz, "lower_half_power_hz")
    f2 = _positive(upper_half_power_hz, "upper_half_power_hz")
    if not f1 < f0 < f2:
        raise DomainError("Need f_1 < f_0 < f_2.", field="lower_half_power_hz")
    return f0 / (f2 - f1)


def damping_rate_from_decay(times_s: Sequence[float], amplitudes: Sequence[float]) -> float:
    """Damping rate ``alpha`` from a ringing-mode envelope ``p_t = p_max e^(-alpha t)``.

    Least-squares slope of ``ln p_t`` against time (SP-194 fits the line
    graphically); ``alpha / f_0`` is the damping rate per cycle.

    Parameters
    ----------
    times_s : sequence of float
        Sample times of the envelope peaks.
    amplitudes : sequence of float
        Positive envelope amplitudes at those times.

    Returns
    -------
    float
        ``alpha`` in 1/s (positive for a decaying mode).

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Sect. 9.6.6.2, p. 524,
           https://ntrs.nasa.gov/citations/19720026079
    """
    t = [real(v, "times_s") for v in times_s]
    y = [log(_positive(v, "amplitudes")) for v in amplitudes]
    if len(t) != len(y) or len(t) < 2:
        raise DomainError("Need at least two (time, amplitude) pairs.", field="amplitudes")
    mean_t = sum(t) / len(t)
    mean_y = sum(y) / len(y)
    sxx = sum((v - mean_t) ** 2 for v in t)
    if sxx == 0.0:
        raise DomainError("Sample times must differ.", field="times_s")
    slope = sum((a - mean_t) * (b - mean_y) for a, b in zip(t, y)) / sxx
    if not isfinite(slope):
        raise DomainError("Decay fit is not finite.", code="non_finite_result")
    return -slope
