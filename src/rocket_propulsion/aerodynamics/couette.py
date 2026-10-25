"""Plane Couette flow: incompressible closed form and compressible quadrature.

The lower plate (``y = 0``) is at rest at ``Tw``; the upper plate (``y = D``)
moves at ``ue`` and is held at ``Te``; the pressure is uniform [1]_.

* Incompressible (constant ``rho``, ``mu``, ``k``): linear velocity and the
  parabolic enthalpy profile of Anderson Eqs. (16.6)-(16.25), recovery factor
  ``r = Pr`` (16.50) and Reynolds analogy ``CH/cf = Pr^-1/2`` (16.59).
* Compressible, constant ``Pr`` and ``cp``, any viscosity law ``mu(T)``: the
  energy integral (16.80) gives ``h(u)`` in closed form, so the constant shear
  follows from ``D = (1/tau) int_0^ue mu(T(u)) du`` (Eq. (16.65) integrated) by
  composite Simpson quadrature, and the wall heat flux from Eq. (16.95). This
  replaces the textbook's nested shooting iteration with an exact reduction
  that is valid under the same constant-``Pr`` assumption (Sec. 16.4.4).

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson, Fundamentals of Aerodynamics, 6th
       ed., McGraw-Hill, 2017, Ch. 16, printed pp. 964-995.
"""

from __future__ import annotations

from dataclasses import dataclass

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

from .boundary_layer import ANDERSON_AIR_SUTHERLAND, SutherlandViscosity

__all__ = [
    "CompressibleCouetteResult",
    "IncompressibleCouetteResult",
    "compressible_couette",
    "incompressible_couette",
]


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class IncompressibleCouetteResult:
    """Constant-property Couette flow solution.

    Attributes
    ----------
    shear_stress_pa : float
        ``tau = mu ue / D``, Eq. (16.9), constant across the gap.
    reynolds : float
        ``Re = rho ue D / mu`` (Sec. 16.3.5).
    skin_friction : float
        ``cf = 2/Re``, Eq. (16.54).
    lower_wall_heat_flux_w_m2, upper_wall_heat_flux_w_m2 : float
        ``(mu/Pr)(he - hw + Pr ue^2/2)/D`` and ``(mu/Pr)(he - hw - Pr ue^2/2)/D``
        (Eqs. (16.24), (16.25)), the enthalpy-gradient values ``(mu/Pr) dh/dy``
        at each wall. At the lower wall a positive value means heat flows into
        that wall.
    adiabatic_wall_temperature_k : float
        ``Taw = Te + Pr ue^2/(2 cp)``, Eq. (16.40).
    recovery_factor : float
        ``r = Pr``, Eq. (16.50).
    stanton : float
        ``CH = 1/(Re Pr)``, Eq. (16.58).
    maximum_temperature_k : float
        Largest temperature on ``[0, D]`` of Eq. (16.16).
    inputs : tuple of float
        ``(ue, D, mu, rho, cp, Pr, Te, Tw)`` replayed by :meth:`profile`.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.6)-(16.59), printed pp. 964-978.
    """

    shear_stress_pa: float
    reynolds: float
    skin_friction: float
    lower_wall_heat_flux_w_m2: float
    upper_wall_heat_flux_w_m2: float
    adiabatic_wall_temperature_k: float
    recovery_factor: float
    stanton: float
    maximum_temperature_k: float
    inputs: tuple[float, ...]

    def profile(self, y_over_gap: float) -> tuple[float, float]:
        """Return ``(u, T)`` at ``y/D`` from Eqs. (16.6) and (16.16).

        Parameters
        ----------
        y_over_gap : float
            Fractional height in ``[0, 1]``.

        Returns
        -------
        tuple of float
            Velocity in m/s and temperature in K.

        References
        ----------
        .. [1] ANDERSON-FOA-2017. Eqs. (16.6) and (16.16), printed pp. 964-968.
        """

        eta = real(y_over_gap, "y_over_gap")
        if not 0.0 <= eta <= 1.0:
            raise DomainError("y_over_gap must lie in [0, 1].", field="y_over_gap")
        ue, _gap, _mu, _rho, cp, pr, te, tw = self.inputs
        temperature = tw + ((te - tw) + 0.5 * pr * ue * ue / cp) * eta - 0.5 * pr * ue * ue / cp * eta * eta
        return ue * eta, temperature


def incompressible_couette(
    *,
    upper_plate_velocity_m_s: float,
    gap_m: float,
    viscosity_pa_s: float,
    density_kg_m3: float,
    specific_heat_j_kg_k: float,
    prandtl: float,
    upper_wall_temperature_k: float,
    lower_wall_temperature_k: float,
) -> IncompressibleCouetteResult:
    """Return the constant-property Couette solution of Anderson Sec. 16.3.

    Parameters
    ----------
    upper_plate_velocity_m_s, gap_m : float
        ``ue`` and ``D``, positive.
    viscosity_pa_s, density_kg_m3, specific_heat_j_kg_k, prandtl : float
        Constant ``mu``, ``rho``, ``cp`` and ``Pr``.
    upper_wall_temperature_k, lower_wall_temperature_k : float
        ``Te`` and ``Tw``.

    Returns
    -------
    IncompressibleCouetteResult
        Shear, heat fluxes, recovery and analogy quantities.

    Notes
    -----
    Valid only while viscous heating leaves ``rho``, ``mu`` and ``k`` nearly
    constant (Anderson's Example 16.1 at Mach 0.18).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.6)-(16.59), Example 16.1, printed
           pp. 964-980.
    """

    ue = _positive(upper_plate_velocity_m_s, "upper_plate_velocity_m_s")
    gap = _positive(gap_m, "gap_m")
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    rho = _positive(density_kg_m3, "density_kg_m3")
    cp = _positive(specific_heat_j_kg_k, "specific_heat_j_kg_k")
    pr = _positive(prandtl, "prandtl")
    te = _positive(upper_wall_temperature_k, "upper_wall_temperature_k")
    tw = _positive(lower_wall_temperature_k, "lower_wall_temperature_k")
    dh = cp * (te - tw)
    dissipation = 0.5 * pr * ue * ue
    reynolds = rho * ue * gap / mu
    # Eq. (16.16) in T: T = Tw + (a) eta - (b) eta^2, a = (Te-Tw) + b, b = Pr ue^2/(2 cp)
    b = dissipation / cp
    a = (te - tw) + b
    eta_peak = min(1.0, max(0.0, a / (2.0 * b)))
    peak = max(tw, te, tw + a * eta_peak - b * eta_peak * eta_peak)
    return IncompressibleCouetteResult(
        shear_stress_pa=mu * ue / gap,
        reynolds=reynolds,
        skin_friction=2.0 / reynolds,
        lower_wall_heat_flux_w_m2=mu / pr * (dh + dissipation) / gap,
        upper_wall_heat_flux_w_m2=mu / pr * (dh - dissipation) / gap,
        adiabatic_wall_temperature_k=te + dissipation / cp,
        recovery_factor=pr,
        stanton=1.0 / (reynolds * pr),
        maximum_temperature_k=peak,
        inputs=(ue, gap, mu, rho, cp, pr, te, tw),
    )


@dataclass(frozen=True, slots=True)
class CompressibleCouetteResult:
    """Compressible constant-``Pr`` Couette solution.

    Attributes
    ----------
    shear_stress_pa : float
        Constant shear ``tau = (1/D) int_0^ue mu(T(u)) du``.
    lower_wall_heat_flux_w_m2 : float
        ``qw = tau/(ue Pr) (he - hw + Pr ue^2/2)``, Eq. (16.95); zero for an
        adiabatic lower wall. Positive means heat flows into the lower wall.
    lower_wall_temperature_k : float
        Prescribed ``Tw`` or, for an adiabatic wall, ``Taw``.
    adiabatic_wall_temperature_k : float
        ``Taw = Te + Pr ue^2/(2 cp)``, Eq. (16.82) (``r = Pr``, Eq. (16.83)).
    skin_friction : float
        ``cf = tau/(rho_e ue^2/2)`` with ``rho_e = p/(R Te)``, Eq. (16.51).
    stanton : float
        ``CH = qw/(rho_e ue (haw - hw))``, Eq. (16.55); ``nan`` if adiabatic.
    maximum_temperature_k : float
        Peak of ``T(u)`` on ``[0, ue]``.
    mach : float
        Upper-plate Mach number ``ue/sqrt(gamma R Te)``.
    quadrature_intervals : int
        Number of Simpson intervals used.
    profile : tuple of (float, float, float)
        ``(y/D, u/ue, T)`` samples at the quadrature nodes.
    warnings : tuple of WarningMessage
        Model-applicability notes.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.51)-(16.98), printed pp. 977-995.
    """

    shear_stress_pa: float
    lower_wall_heat_flux_w_m2: float
    lower_wall_temperature_k: float
    adiabatic_wall_temperature_k: float
    skin_friction: float
    stanton: float
    maximum_temperature_k: float
    mach: float
    quadrature_intervals: int
    profile: tuple[tuple[float, float, float], ...]
    warnings: tuple[WarningMessage, ...]


def compressible_couette(
    *,
    upper_plate_velocity_m_s: float,
    gap_m: float,
    pressure_pa: float,
    upper_wall_temperature_k: float,
    lower_wall_temperature_k: float | None,
    gas_constant_j_kg_k: float,
    gamma: float = 1.4,
    prandtl: float = 0.71,
    viscosity: SutherlandViscosity = ANDERSON_AIR_SUTHERLAND,
    intervals: int = 400,
) -> CompressibleCouetteResult:
    """Solve compressible Couette flow for constant ``Pr`` and ``cp``.

    Parameters
    ----------
    upper_plate_velocity_m_s, gap_m, pressure_pa : float
        ``ue``, ``D`` and the uniform pressure.
    upper_wall_temperature_k : float
        ``Te``.
    lower_wall_temperature_k : float or None
        ``Tw``; ``None`` makes the lower wall adiabatic (``hw = haw``).
    gas_constant_j_kg_k, gamma, prandtl : float
        Calorically perfect gas with ``cp = gamma R/(gamma - 1)``.
    viscosity : SutherlandViscosity, optional
        Viscosity law ``mu(T)`` (Anderson uses Sutherland's law, Eq. (16.63)).
    intervals : int, optional
        Even number of Simpson intervals in ``u`` (at least 20).

    Returns
    -------
    CompressibleCouetteResult
        Shear, wall heat flux, coefficients and a ``(y/D, u/ue, T)`` profile.

    Notes
    -----
    From Eq. (16.80), ``h(u) = hw + (u/ue)(he - hw) + (Pr/2) u (ue - u)``,
    which holds for any ``mu(T)`` because only ``tau = const`` and constant
    ``Pr`` were used. Since ``du/dy = tau/mu`` (Eq. (16.65)),
    ``y(u) = (1/tau) int_0^u mu du'`` and ``y(ue) = D`` fixes ``tau``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (16.61), (16.65), (16.80)-(16.98),
           Example 16.2, printed pp. 981-995.
    """

    ue = _positive(upper_plate_velocity_m_s, "upper_plate_velocity_m_s")
    gap = _positive(gap_m, "gap_m")
    pressure = _positive(pressure_pa, "pressure_pa")
    te = _positive(upper_wall_temperature_k, "upper_wall_temperature_k")
    gas = _positive(gas_constant_j_kg_k, "gas_constant_j_kg_k")
    g = real(gamma, "gamma")
    if g <= 1.0:
        raise DomainError("gamma must exceed 1.", field="gamma")
    pr = _positive(prandtl, "prandtl")
    if isinstance(intervals, bool) or not isinstance(intervals, int) or intervals < 20 or intervals % 2:
        raise DomainError("intervals must be an even integer of at least 20.", field="intervals")
    cp = g * gas / (g - 1.0)
    he = cp * te
    taw = te + pr * ue * ue / (2.0 * cp)
    adiabatic = lower_wall_temperature_k is None
    tw = taw if adiabatic else _positive(lower_wall_temperature_k, "lower_wall_temperature_k")
    hw = cp * tw

    def temperature(u: float) -> float:
        return (hw + (u / ue) * (he - hw) + 0.5 * pr * u * (ue - u)) / cp

    step = ue / intervals
    nodes = [index * step for index in range(intervals + 1)]
    temperatures = [temperature(u) for u in nodes]
    if min(temperatures) <= 0.0:
        raise DomainError("The energy integral gives a non-positive temperature.", field="lower_wall_temperature_k")
    mus = [viscosity.viscosity(t) for t in temperatures]
    # Cumulative Simpson integral of mu du at even nodes; odd nodes use the
    # quadratic through the same three points (h/12 (5 f0 + 8 f1 - f2)).
    cumulative = [0.0] * (intervals + 1)
    for index in range(2, intervals + 1, 2):
        cumulative[index] = cumulative[index - 2] + step / 3.0 * (mus[index - 2] + 4.0 * mus[index - 1] + mus[index])
        cumulative[index - 1] = cumulative[index - 2] + step / 12.0 * (
            5.0 * mus[index - 2] + 8.0 * mus[index - 1] - mus[index]
        )
    integral = cumulative[-1]
    tau = integral / gap
    qw = tau / (ue * pr) * (he - hw + 0.5 * pr * ue * ue)
    rho_e = pressure / (gas * te)
    cf = tau / (0.5 * rho_e * ue * ue)
    stanton = float("nan") if adiabatic else qw / (rho_e * ue * (cp * taw - hw))
    profile = tuple((cumulative[i] / integral, nodes[i] / ue, temperatures[i]) for i in range(intervals + 1))
    peak = max(temperatures)
    warnings: list[WarningMessage] = []
    if peak > 1000.0:
        warnings.append(WarningMessage(
            "calorically_perfect_limit",
            "Peak temperature above about 1000 K; Anderson justifies constant cp and Pr only below roughly 1000 K.",
        ))
    return CompressibleCouetteResult(
        shear_stress_pa=tau,
        lower_wall_heat_flux_w_m2=0.0 if adiabatic else qw,
        lower_wall_temperature_k=tw,
        adiabatic_wall_temperature_k=taw,
        skin_friction=cf,
        stanton=stanton,
        maximum_temperature_k=peak,
        mach=ue / (g * gas * te) ** 0.5,
        quadrature_intervals=intervals,
        profile=profile,
        warnings=tuple(warnings),
    )
