"""Ballistic entry: Allen-Eggers closed form and a planar point-mass integrator.

Allen-Eggers [1]_ assume a constant drag coefficient, an isothermal
exponential atmosphere ``rho = rho0 exp(-beta y)``, negligible gravity and a
straight flight path at the entry angle ``theta_E`` (measured below the
horizontal). Their closed forms give velocity, deceleration and the altitudes
and values of peak deceleration (Eqs. (13)-(17), (20)-(21)), peak average
heating rate (Eqs. (36)-(39), (41)) and peak stagnation-point heating rate
for ``qs = C sqrt(rho/sigma) V^3`` (Eqs. (45)-(48), (50)); the total heat input
follows Eq. (33) (equivalently Anderson Eq. (14.60) [2]_).

The numerical model integrates the two polar-coordinate equations of motion
of a point mass over a spherical, non-rotating planet in a spherically
symmetric atmosphere, Chapman Eqs. (3)-(7) [3]_, with lift-to-drag ratio,
inverse-square gravity ``g = g0 (r0/r)^2`` (U.S. Standard Atmosphere 1976
Eq. (17), ``r0 = 6 356 766`` m [4]_) and either an exponential atmosphere or
the 1976 standard atmosphere below 86 km geometric altitude (geopotential
Eq. (18), density Eq. (42); the pressure and layer table are reused from
:mod:`rocket_propulsion.compressible.off_design_nozzle`). Planet rotation,
oblateness, winds, variable ``C_D`` and ablation mass loss are not modelled.

References
----------
.. [1] ALLEN-EGGERS-NACA-R1381. H. J. Allen and A. J. Eggers Jr., NACA
       Report 1381 (supersedes NACA TN 4047, 1957),
       https://ntrs.nasa.gov/citations/19930091020
.. [2] ANDERSON-FOA-2017. J. D. Anderson, Fundamentals of Aerodynamics, 6th
       ed., McGraw-Hill, 2017, Sec. 14.8, printed pp. 902-905.
.. [3] CHAPMAN-1958-TN4276. D. R. Chapman, NACA TN 4276, May 1958, printed
       pp. 10-11, https://ntrs.nasa.gov/citations/19930085059
.. [4] USSA-1976. U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562,
       printed pp. 3-15, https://ntrs.nasa.gov/citations/19770009539
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from math import atan2, cos, degrees, exp, log, radians, sin, sqrt

from rocket_propulsion.compressible.off_design_nozzle import (
    _ATMOSPHERE_LAYERS,
    standard_atmosphere_pressure_pa,
)
from rocket_propulsion.core.constants import (
    STANDARD_GRAVITY_M_S2,
    USSA1976_GAS_CONSTANT_J_MOL_K,
    USSA1976_MOLAR_MASS_KG_MOL,
)
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

__all__ = [
    "USSA1976_EARTH_RADIUS_M",
    "USSA1976_MAX_GEOMETRIC_ALTITUDE_M",
    "AllenEggersEntry",
    "EntrySample",
    "EntryTrajectory",
    "ExponentialAtmosphere",
    "allen_eggers_entry",
    "integrate_planar_entry",
    "total_heat_input_reynolds_analogy",
    "ussa1976_density",
    "ussa1976_geopotential_altitude",
]

USSA1976_EARTH_RADIUS_M = 6.356766e6
"""Effective earth radius ``r0`` of the 1976 standard (Table 2, p. 2)."""

_GEOPOTENTIAL_CEILING_M = 84_852.0
USSA1976_MAX_GEOMETRIC_ALTITUDE_M = (
    USSA1976_EARTH_RADIUS_M * _GEOPOTENTIAL_CEILING_M / (USSA1976_EARTH_RADIUS_M - _GEOPOTENTIAL_CEILING_M)
)
"""Geometric altitude of the 84.852 km' layer top (about 86 km), Eq. (19)."""


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


# ---------------------------------------------------------------------------
# Atmospheres
# ---------------------------------------------------------------------------


def ussa1976_geopotential_altitude(geometric_altitude_m: float) -> float:
    """Return geopotential altitude ``H = r0 Z/(r0 + Z)`` (Gamma = 1 m'/m).

    Parameters
    ----------
    geometric_altitude_m : float
        Geometric altitude ``Z`` in m, non-negative.

    Returns
    -------
    float
        Geopotential altitude in m'.

    References
    ----------
    .. [1] USSA-1976. Eq. (18) and ``r0`` of Table 2, printed pp. 2, 8.
    """

    z = _nonnegative(geometric_altitude_m, "geometric_altitude_m")
    return USSA1976_EARTH_RADIUS_M * z / (USSA1976_EARTH_RADIUS_M + z)


def ussa1976_density(geometric_altitude_m: float) -> float:
    """Return U.S. Standard Atmosphere 1976 density below 86 km geometric altitude.

    Parameters
    ----------
    geometric_altitude_m : float
        Geometric altitude in m, from 0 to
        :data:`USSA1976_MAX_GEOMETRIC_ALTITUDE_M` (about 86 km).

    Returns
    -------
    float
        Mass density in kg/m^3.

    Raises
    ------
    DomainError
        Above 86 km, where the 1976 model is defined by species number
        densities that are not implemented here.

    Notes
    -----
    ``rho = P M0/(R* T_M)``, Eq. (42), with ``T_M = T_M,b + L_M,b (H - H_b)``,
    Eq. (23), from the repository's layer table and pressure function
    (Eqs. (33a)/(33b)). Molecular-scale temperature is used, so no
    ``M/M0`` correction is needed between 80 and 86 km.

    References
    ----------
    .. [1] USSA-1976. Eqs. (18), (23), (33a), (33b), (42), Table 4, printed
           pp. 3, 8, 10, 12, 15.
    """

    z = _nonnegative(geometric_altitude_m, "geometric_altitude_m")
    if z > USSA1976_MAX_GEOMETRIC_ALTITUDE_M:
        raise DomainError("The 1976 standard atmosphere is implemented only up to 86 km geometric altitude.",
                          field="geometric_altitude_m")
    h = min(ussa1976_geopotential_altitude(z), _GEOPOTENTIAL_CEILING_M)
    base = _ATMOSPHERE_LAYERS[0]
    for layer in _ATMOSPHERE_LAYERS:
        if layer[0] <= h:
            base = layer
        else:
            break
    molecular_temperature = base[1] + base[3] * (h - base[0])
    pressure = standard_atmosphere_pressure_pa(h)
    return pressure * USSA1976_MOLAR_MASS_KG_MOL / (USSA1976_GAS_CONSTANT_J_MOL_K * molecular_temperature)


@dataclass(frozen=True, slots=True)
class ExponentialAtmosphere:
    """Isothermal exponential atmosphere ``rho = rho0 exp(-beta y)``.

    Attributes
    ----------
    surface_density_kg_m3 : float
        ``rho0``.
    inverse_scale_height_1_m : float
        ``beta`` in 1/m.

    Notes
    -----
    Allen and Eggers fit ``rho0 = 0.0034`` slug/ft^3 and ``beta = 1/22000``
    1/ft to the then-standard atmosphere between 20 000 and 180 000 ft
    (their Eq. (7)) [1]_.

    References
    ----------
    .. [1] ALLEN-EGGERS-NACA-R1381. Eqs. (2) and (7), printed p. 2.
    """

    surface_density_kg_m3: float
    inverse_scale_height_1_m: float

    def __post_init__(self) -> None:
        _positive(self.surface_density_kg_m3, "surface_density_kg_m3")
        _positive(self.inverse_scale_height_1_m, "inverse_scale_height_1_m")

    def density(self, altitude_m: float) -> float:
        """Return density at ``altitude_m``.

        Parameters
        ----------
        altitude_m : float
            Altitude in m (negative values allowed for the analytic model).

        Returns
        -------
        float
            Density in kg/m^3.

        References
        ----------
        .. [1] ALLEN-EGGERS-NACA-R1381. Eq. (2), printed p. 2.
        """

        return self.surface_density_kg_m3 * exp(-self.inverse_scale_height_1_m * real(altitude_m, "altitude_m"))


# ---------------------------------------------------------------------------
# Allen-Eggers closed form
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AllenEggersEntry:
    """Allen-Eggers ballistic-entry closed-form results.

    Attributes
    ----------
    entry_velocity_m_s, entry_angle_deg : float
        ``V_E`` and ``theta_E`` (below horizontal).
    ballistic_parameter : float
        ``B = C_D rho0 A/(beta m sin theta_E)``, Eq. (51).
    peak_deceleration_altitude_m : float
        ``y1 = ln(B)/beta``, Eq. (15); negative means the maximum in flight
        occurs at sea level.
    peak_deceleration_velocity_m_s : float
        ``V_E e^-1/2`` (Eq. (16)) if ``y1 >= 0``, else the sea-level speed (20).
    peak_deceleration_g : float
        ``beta V_E^2 sin theta_E/(2 g e)``, Eq. (17), or Eq. (21) at sea level.
    peak_average_heating_altitude_m, peak_average_heating_velocity_m_s : float
        ``y2 = ln(3B/2)/beta`` and ``V_E e^-1/3``, Eqs. (38)-(39).
    peak_stagnation_heating_altitude_m, peak_stagnation_heating_velocity_m_s : float
        ``y3 = ln(3B)/beta`` and ``V_E e^-1/6``, Eqs. (47)-(48).
    peak_stagnation_heating_w_m2 : float or None
        Eq. (46) (or (50) when ``y3 < 0``) for the supplied heating constant.
    impact_velocity_m_s : float
        ``V_0 = V_E exp(-B/2)``, the velocity at ``y = 0``.
    standard_gravity_m_s2 : float
        ``g`` used to express decelerations in g.

    References
    ----------
    .. [1] ALLEN-EGGERS-NACA-R1381. Eqs. (13)-(21), (33)-(51), printed
           pp. 2-8.
    """

    entry_velocity_m_s: float
    entry_angle_deg: float
    surface_density_kg_m3: float
    inverse_scale_height_1_m: float
    ballistic_parameter: float
    peak_deceleration_altitude_m: float
    peak_deceleration_velocity_m_s: float
    peak_deceleration_g: float
    peak_average_heating_altitude_m: float
    peak_average_heating_velocity_m_s: float
    peak_stagnation_heating_altitude_m: float
    peak_stagnation_heating_velocity_m_s: float
    peak_stagnation_heating_w_m2: float | None
    impact_velocity_m_s: float
    standard_gravity_m_s2: float

    def velocity(self, altitude_m: float) -> float:
        """Return ``V = V_E exp(-(B/2) e^(-beta y))``, Eq. (13).

        Parameters
        ----------
        altitude_m : float
            Altitude ``y`` in m.

        Returns
        -------
        float
            Velocity in m/s.

        References
        ----------
        .. [1] ALLEN-EGGERS-NACA-R1381. Eq. (13), printed p. 4.
        """

        y = real(altitude_m, "altitude_m")
        return self.entry_velocity_m_s * exp(-0.5 * self.ballistic_parameter * exp(-self.inverse_scale_height_1_m * y))

    def deceleration_g(self, altitude_m: float) -> float:
        """Return ``-dV/dt / g`` from Eq. (14).

        Parameters
        ----------
        altitude_m : float
            Altitude ``y`` in m.

        Returns
        -------
        float
            Deceleration in units of ``standard_gravity_m_s2``.

        References
        ----------
        .. [1] ALLEN-EGGERS-NACA-R1381. Eq. (14), printed p. 4.
        """

        y = real(altitude_m, "altitude_m")
        beta = self.inverse_scale_height_1_m
        sin_theta = sin(radians(self.entry_angle_deg))
        # C_D rho0 A / (2 m) = B beta sin(theta)/2
        k = 0.5 * self.ballistic_parameter * beta * sin_theta
        return k * self.entry_velocity_m_s**2 * exp(-beta * y) * exp(
            -self.ballistic_parameter * exp(-beta * y)
        ) / self.standard_gravity_m_s2


def allen_eggers_entry(
    *,
    entry_velocity_m_s: float,
    entry_angle_deg: float,
    mass_kg: float,
    drag_coefficient: float,
    reference_area_m2: float,
    atmosphere: ExponentialAtmosphere,
    stagnation_heating_constant: float | None = None,
    nose_radius_m: float | None = None,
    standard_gravity_m_s2: float = STANDARD_GRAVITY_M_S2,
) -> AllenEggersEntry:
    """Return the Allen-Eggers closed-form ballistic entry.

    Parameters
    ----------
    entry_velocity_m_s : float
        ``V_E`` at the top of the sensible atmosphere.
    entry_angle_deg : float
        ``theta_E`` below the horizontal, in ``(0, 90]`` degrees.
    mass_kg, drag_coefficient, reference_area_m2 : float
        ``m``, constant ``C_D`` and ``A``.
    atmosphere : ExponentialAtmosphere
        ``rho0`` and ``beta``.
    stagnation_heating_constant : float, optional
        ``C`` in ``qs = C sqrt(rho/sigma) V^3`` (SI). Allen-Eggers derive
        ``C = 6.8e-6`` in English units (Eq. (44)); Tauber's Eq. (40) gives
        ``C = 1.83e-4`` in SI for a cold wall. Required with ``nose_radius_m``
        to evaluate the peak stagnation heating rate.
    nose_radius_m : float, optional
        ``sigma``.
    standard_gravity_m_s2 : float, optional
        ``g`` used to express decelerations in g.

    Returns
    -------
    AllenEggersEntry
        Closed-form peak conditions and the velocity/deceleration functions.

    References
    ----------
    .. [1] ALLEN-EGGERS-NACA-R1381. Eqs. (13)-(21), (36)-(51), printed
           pp. 2-8.
    """

    ve = _positive(entry_velocity_m_s, "entry_velocity_m_s")
    theta_deg = real(entry_angle_deg, "entry_angle_deg")
    if not 0.0 < theta_deg <= 90.0:
        raise DomainError("entry_angle_deg must lie in (0, 90].", field="entry_angle_deg")
    mass = _positive(mass_kg, "mass_kg")
    cd = _positive(drag_coefficient, "drag_coefficient")
    area = _positive(reference_area_m2, "reference_area_m2")
    if not isinstance(atmosphere, ExponentialAtmosphere):
        raise InputError("atmosphere must be an ExponentialAtmosphere.", field="atmosphere")
    g = _positive(standard_gravity_m_s2, "standard_gravity_m_s2")
    rho0 = atmosphere.surface_density_kg_m3
    beta = atmosphere.inverse_scale_height_1_m
    sin_theta = sin(radians(theta_deg))
    b = cd * rho0 * area / (beta * mass * sin_theta)
    y1 = log(b) / beta
    impact_velocity = ve * exp(-0.5 * b)
    if y1 >= 0.0:
        v1 = ve * exp(-0.5)
        peak_g = beta * ve * ve * sin_theta / (2.0 * g * exp(1.0))
    else:
        v1 = impact_velocity
        peak_g = cd * rho0 * area * ve * ve / (2.0 * mass * g) * exp(-b)
    y2 = log(1.5 * b) / beta
    y3 = log(3.0 * b) / beta
    peak_q = None
    if (stagnation_heating_constant is None) != (nose_radius_m is None):
        raise InputError("Give both stagnation_heating_constant and nose_radius_m, or neither.",
                         field="stagnation_heating_constant")
    if stagnation_heating_constant is not None:
        c = _positive(stagnation_heating_constant, "stagnation_heating_constant")
        sigma = _positive(nose_radius_m, "nose_radius_m")
        if y3 >= 0.0:
            peak_q = c * sqrt(beta * mass * sin_theta / (3.0 * exp(1.0) * cd * sigma * area)) * ve**3
        else:
            peak_q = c * sqrt(rho0 / sigma) * ve**3 * exp(-1.5 * b)
    return AllenEggersEntry(
        entry_velocity_m_s=ve,
        entry_angle_deg=theta_deg,
        surface_density_kg_m3=rho0,
        inverse_scale_height_1_m=beta,
        ballistic_parameter=b,
        peak_deceleration_altitude_m=y1,
        peak_deceleration_velocity_m_s=v1,
        peak_deceleration_g=peak_g,
        peak_average_heating_altitude_m=y2,
        peak_average_heating_velocity_m_s=ve * exp(-1.0 / 3.0),
        peak_stagnation_heating_altitude_m=y3,
        peak_stagnation_heating_velocity_m_s=ve * exp(-1.0 / 6.0),
        peak_stagnation_heating_w_m2=peak_q,
        impact_velocity_m_s=impact_velocity,
        standard_gravity_m_s2=g,
    )


def total_heat_input_reynolds_analogy(
    *, mass_kg: float, entry_velocity_m_s: float, final_velocity_m_s: float,
    friction_to_drag_ratio: float,
) -> float:
    """Return total convective heat input ``Q = (1/4)(Cf/CD) m (VE^2 - V0^2)``.

    Parameters
    ----------
    mass_kg : float
        Vehicle mass.
    entry_velocity_m_s, final_velocity_m_s : float
        ``V_E`` and the final (impact) velocity ``V_0 < V_E``.
    friction_to_drag_ratio : float
        ``Cf/CD`` (Allen-Eggers' ``C_f' S/(C_D A)``), in ``(0, 1]``.

    Returns
    -------
    float
        Heat input in J.

    Notes
    -----
    Allen-Eggers Eq. (34) (Reynolds analogy, ``Pr = 1``); with ``V_0 = 0`` it
    is Anderson Eq. (14.60), ``Q = (1/2)(Cf/CD)(m V_E^2/2)``.

    References
    ----------
    .. [1] ALLEN-EGGERS-NACA-R1381. Eqs. (33)-(34), printed p. 6.
    .. [2] ANDERSON-FOA-2017. Eq. (14.60), printed p. 904.
    """

    mass = _positive(mass_kg, "mass_kg")
    ve = _positive(entry_velocity_m_s, "entry_velocity_m_s")
    v0 = _nonnegative(final_velocity_m_s, "final_velocity_m_s")
    ratio = _positive(friction_to_drag_ratio, "friction_to_drag_ratio")
    if ratio > 1.0:
        raise DomainError("friction_to_drag_ratio must not exceed 1.", field="friction_to_drag_ratio")
    if v0 >= ve:
        raise DomainError("final_velocity_m_s must be below entry_velocity_m_s.", field="final_velocity_m_s")
    return 0.25 * ratio * mass * (ve * ve - v0 * v0)


# ---------------------------------------------------------------------------
# Planar point-mass integrator
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class EntrySample:
    """One trajectory sample.

    Attributes
    ----------
    time_s, altitude_m, downrange_m, velocity_m_s, flight_path_angle_deg : float
        Time, altitude above ``planet_radius_m``, surface arc length, speed and
        flight-path angle (negative descending, Chapman's ``phi``).
    density_kg_m3 : float
        Atmospheric density.
    deceleration_g : float
        Aerodynamic specific force ``sqrt(D^2 + L^2)/(m g0)``.
    stagnation_heating_w_m2 : float or None
        ``C sqrt(rho/r_n) V^3`` if a heating constant was supplied.

    References
    ----------
    .. [1] CHAPMAN-1958-TN4276. Eqs. (3)-(7), printed pp. 10-11.
    """

    time_s: float
    altitude_m: float
    downrange_m: float
    velocity_m_s: float
    flight_path_angle_deg: float
    density_kg_m3: float
    deceleration_g: float
    stagnation_heating_w_m2: float | None


@dataclass(frozen=True, slots=True)
class EntryTrajectory:
    """Integrated planar entry trajectory and its peaks.

    Attributes
    ----------
    samples : tuple of EntrySample
        Samples at every accepted step, ending at the termination event.
    termination : str
        ``"altitude"`` (reached ``final_altitude_m``), ``"velocity"``
        (speed fell below ``final_velocity_m_s``) or ``"time"``.
    peak_deceleration : EntrySample
        Sample with the largest deceleration.
    peak_heating : EntrySample or None
        Sample with the largest stagnation heating, if computed.
    integrated_heat_load_j_m2 : float or None
        Trapezoidal time integral of the stagnation heating.
    time_step_s : float
        RK4 step.
    warnings : tuple of WarningMessage
        Modelling notes.

    References
    ----------
    .. [1] CHAPMAN-1958-TN4276. Eqs. (3)-(7), printed pp. 10-11.
    """

    samples: tuple[EntrySample, ...]
    termination: str
    peak_deceleration: EntrySample
    peak_heating: EntrySample | None
    integrated_heat_load_j_m2: float | None
    time_step_s: float
    warnings: tuple[WarningMessage, ...]


_EntryState = tuple[float, float, float, float]


def _entry_state(values: Iterable[float]) -> _EntryState:
    """Pack ``(r, theta, v, u)`` from an RK stage expression into a 4-tuple."""
    r, theta, v, u = values
    return (r, theta, v, u)


def _stagnation_heating(item: EntrySample) -> float:
    """Return the stagnation heating of a sample integrated with heating on."""
    if item.stagnation_heating_w_m2 is None:  # pragma: no cover - set for every sample
        raise DomainError("Entry sample has no stagnation heating rate.",
                          field="stagnation_heating_w_m2")
    return item.stagnation_heating_w_m2


def integrate_planar_entry(
    *,
    entry_altitude_m: float,
    entry_velocity_m_s: float,
    entry_flight_path_angle_deg: float,
    mass_kg: float,
    drag_coefficient: float,
    reference_area_m2: float,
    atmosphere: str | ExponentialAtmosphere | Callable[[float], float] = "ussa1976",
    lift_to_drag: float = 0.0,
    gravity: str = "inverse_square",
    curvature: bool = True,
    planet_radius_m: float = USSA1976_EARTH_RADIUS_M,
    surface_gravity_m_s2: float = STANDARD_GRAVITY_M_S2,
    nose_radius_m: float | None = None,
    stagnation_heating_constant: float | None = None,
    time_step_s: float = 0.05,
    final_altitude_m: float = 0.0,
    final_velocity_m_s: float = 1.0,
    max_time_s: float = 5000.0,
) -> EntryTrajectory:
    """Integrate a planar point-mass entry with fixed-step RK4.

    Parameters
    ----------
    entry_altitude_m, entry_velocity_m_s : float
        Initial altitude above ``planet_radius_m`` and speed.
    entry_flight_path_angle_deg : float
        Initial flight-path angle in degrees, negative for descent, in
        ``[-90, 0)``.
    mass_kg, drag_coefficient, reference_area_m2 : float
        Constant mass, ``C_D`` and ``A``.
    atmosphere : {"ussa1976"} or ExponentialAtmosphere or callable
        Density model. ``"ussa1976"`` requires all altitudes to stay below
        about 86 km (:data:`USSA1976_MAX_GEOMETRIC_ALTITUDE_M`). A callable
        maps altitude in m to density in kg/m^3.
    lift_to_drag : float, optional
        Constant ``L/D`` with lift in the trajectory plane (positive up).
    gravity : {"inverse_square", "none"}
        ``g = g0 (r0/r)^2`` or no gravity (the Allen-Eggers assumption).
    curvature : bool, optional
        ``False`` drops the ``u^2/r`` and ``u v/r`` terms and integrates
        downrange as ``dx/dt = u``: the flat-earth equations of Allen-Eggers
        Eq. (1). Use with ``gravity="none"`` to reproduce their closed form.
    planet_radius_m, surface_gravity_m_s2 : float, optional
        ``r0`` and ``g0``; defaults are the 1976 standard's values.
    nose_radius_m, stagnation_heating_constant : float, optional
        Enable ``qs = C sqrt(rho/r_n) V^3`` along the path (e.g. Tauber's
        cold-wall ``C = 1.83e-4``).
    time_step_s : float, optional
        Fixed RK4 step.
    final_altitude_m, final_velocity_m_s, max_time_s : float, optional
        Termination events; the last sample is linearly interpolated onto the
        altitude event.

    Returns
    -------
    EntryTrajectory
        Samples, termination reason, peak deceleration and heating.

    Notes
    -----
    State ``(r, theta, v, u)`` with radial and tangential velocity
    components; with ``tan(phi) = v/u`` Chapman's Eqs. (6)-(7) give
    ``dv/dt = -g + u^2/r + (L cos(phi) - D sin(phi))/m`` and
    ``du/dt = -u v/r - (D cos(phi) + L sin(phi))/m`` [1]_. Chapman then drops
    ``u v/r`` for his analytic theory; this integrator keeps every term.

    References
    ----------
    .. [1] CHAPMAN-1958-TN4276. Eqs. (3)-(7), printed pp. 10-11.
    .. [2] USSA-1976. Eq. (17) and Table 2, printed pp. 2, 8.
    .. [3] ALLEN-EGGERS-NACA-R1381. Eq. (1), printed p. 2.
    """

    h0 = _nonnegative(entry_altitude_m, "entry_altitude_m")
    v_entry = _positive(entry_velocity_m_s, "entry_velocity_m_s")
    gamma_deg = real(entry_flight_path_angle_deg, "entry_flight_path_angle_deg")
    if not -90.0 <= gamma_deg < 0.0:
        raise DomainError("entry_flight_path_angle_deg must lie in [-90, 0).", field="entry_flight_path_angle_deg")
    mass = _positive(mass_kg, "mass_kg")
    cd = _positive(drag_coefficient, "drag_coefficient")
    area = _positive(reference_area_m2, "reference_area_m2")
    lod = real(lift_to_drag, "lift_to_drag")
    radius0 = _positive(planet_radius_m, "planet_radius_m")
    g0 = _positive(surface_gravity_m_s2, "surface_gravity_m_s2")
    dt = _positive(time_step_s, "time_step_s")
    h_final = real(final_altitude_m, "final_altitude_m")
    v_final = _nonnegative(final_velocity_m_s, "final_velocity_m_s")
    t_max = _positive(max_time_s, "max_time_s")
    if h_final >= h0:
        raise DomainError("final_altitude_m must be below entry_altitude_m.", field="final_altitude_m")
    if gravity not in {"inverse_square", "none"}:
        raise InputError("gravity must be 'inverse_square' or 'none'.", field="gravity")
    if not isinstance(curvature, bool):
        raise InputError("curvature must be a bool.", field="curvature")
    if isinstance(atmosphere, ExponentialAtmosphere):
        density: Callable[[float], float] = atmosphere.density
    elif atmosphere == "ussa1976":
        if h0 > USSA1976_MAX_GEOMETRIC_ALTITUDE_M:
            raise DomainError("Entry altitude is above the implemented 1976 atmosphere (86 km).",
                              field="entry_altitude_m")
        if h_final < 0.0:
            raise DomainError("final_altitude_m must be >= 0 with the 1976 atmosphere.", field="final_altitude_m")
        density = ussa1976_density
    elif callable(atmosphere):
        density = atmosphere
    else:
        raise InputError("atmosphere must be 'ussa1976', an ExponentialAtmosphere or a callable.",
                         field="atmosphere")
    heating = None
    if (nose_radius_m is None) != (stagnation_heating_constant is None):
        raise InputError("Give both nose_radius_m and stagnation_heating_constant, or neither.",
                         field="nose_radius_m")
    if nose_radius_m is not None:
        heating = (_positive(stagnation_heating_constant, "stagnation_heating_constant"),
                   _positive(nose_radius_m, "nose_radius_m"))
    k_drag = 0.5 * cd * area / mass

    def derivatives(state: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
        r, _theta, v, u = state
        # An RK stage of the terminal step may dip below the surface; it then sees the
        # surface density. The step itself is truncated at the altitude crossing below.
        rho = density(max(r - radius0, 0.0))
        speed = sqrt(u * u + v * v)
        drag = k_drag * rho * speed * speed
        lift = lod * drag
        phi = atan2(v, u)
        g = g0 * (radius0 / r) ** 2 if gravity == "inverse_square" else 0.0
        centripetal = u * u / r if curvature else 0.0
        coriolis = u * v / r if curvature else 0.0
        dv = -g + centripetal + lift * cos(phi) - drag * sin(phi)
        du = -coriolis - (drag * cos(phi) + lift * sin(phi))
        return (v, u / (r if curvature else radius0), dv, du)

    def sample(time: float, state: tuple[float, float, float, float]) -> EntrySample:
        r, theta, v, u = state
        rho = density(r - radius0)
        speed = sqrt(u * u + v * v)
        drag = k_drag * rho * speed * speed
        q = None if heating is None else heating[0] * sqrt(rho / heating[1]) * speed**3
        return EntrySample(
            time_s=time,
            altitude_m=r - radius0,
            downrange_m=theta * radius0,
            velocity_m_s=speed,
            flight_path_angle_deg=degrees(atan2(v, u)),
            density_kg_m3=rho,
            deceleration_g=drag * sqrt(1.0 + lod * lod) / g0,
            stagnation_heating_w_m2=q,
        )

    gamma = radians(gamma_deg)
    state = (radius0 + h0, 0.0, v_entry * sin(gamma), v_entry * cos(gamma))
    time = 0.0
    samples = [sample(time, state)]
    termination = "time"
    while time < t_max:
        k1 = derivatives(state)
        k2 = derivatives(_entry_state(s + 0.5 * dt * k for s, k in zip(state, k1)))
        k3 = derivatives(_entry_state(s + 0.5 * dt * k for s, k in zip(state, k2)))
        k4 = derivatives(_entry_state(s + dt * k for s, k in zip(state, k3)))
        new_state = _entry_state(s + dt / 6.0 * (a + 2.0 * b + 2.0 * c + d)
                                 for s, a, b, c, d in zip(state, k1, k2, k3, k4))
        new_altitude = new_state[0] - radius0
        if new_altitude <= h_final:
            old_altitude = state[0] - radius0
            fraction = (old_altitude - h_final) / (old_altitude - new_altitude)
            crossing = _entry_state(s + fraction * (n - s) for s, n in zip(state, new_state))
            crossing = (radius0 + h_final, crossing[1], crossing[2], crossing[3])
            time += fraction * dt
            samples.append(sample(time, crossing))
            termination = "altitude"
            break
        state = new_state
        time += dt
        samples.append(sample(time, state))
        if samples[-1].velocity_m_s <= v_final:
            termination = "velocity"
            break
    peak_decel = max(samples, key=lambda item: item.deceleration_g)
    peak_heat = None
    load = None
    if heating is not None:
        rates = [_stagnation_heating(item) for item in samples]
        peak_heat = samples[max(range(len(samples)), key=rates.__getitem__)]
        load = sum(0.5 * (qa + qb) * (b.time_s - a.time_s)
                   for (a, qa), (b, qb) in itertools.pairwise(zip(samples, rates, strict=True)))
    warnings = [WarningMessage(
        "point_mass_planar_model",
        "Non-rotating spherical planet, constant C_D and mass, no winds, no ablation.",
        severity="info",
    )]
    if termination == "time":
        warnings.append(WarningMessage("max_time_reached", "Integration stopped at max_time_s before an event."))
    if heating is not None:
        warnings.append(WarningMessage(
            "stagnation_heating_cold_wall",
            "Stagnation heating uses qs = C sqrt(rho/rn) V^3 without a wall-enthalpy factor (cold wall).",
            severity="info",
        ))
    return EntryTrajectory(
        samples=tuple(samples),
        termination=termination,
        peak_deceleration=peak_decel,
        peak_heating=peak_heat,
        integrated_heat_load_j_m2=load,
        time_step_s=dt,
        warnings=tuple(warnings),
    )
