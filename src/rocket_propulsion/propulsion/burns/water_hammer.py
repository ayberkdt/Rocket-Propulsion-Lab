"""Distributed single-pipe water hammer by the method of characteristics.

A uniform liquid line runs from a constant-pressure reservoir to a valve whose
relative flow coefficient follows a prescribed opening or closing law.  The
continuity and momentum equations of one-dimensional hydraulic transients are
transformed into compatibility equations along the characteristic lines
``dx/dt = +a`` and ``dx/dt = -a`` and solved on a Courant-number-one grid,
``dt = dx / a`` [1]_.  This is the pressure/volume-flow form of the classical
method of characteristics; the NESC report used here as the primary source
gives it in its Appendix B [1]_.

The module reports pressure and flow histories at chosen stations, the peak
pressure, the Joukowsky surge ``rho a dV`` [1]_ for comparison, and mass and
acoustic-energy balance diagnostics.  The liquid density and wave speed are
constant; the wave speed may be computed from the pipe-elasticity formula and
boundary factors of the same report [1]_.

Priming of an evacuated or gas-filled line is deliberately not modelled.  The
NASA sources read for this module describe priming qualitatively and bound the
dead-end surge by the Joukowsky relation for a known impact velocity [1]_,
[2]_, but give no closed-form model of the advancing liquid front; one is not
invented here.

References
----------
.. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
       Transient Pressure Events in Space Flight Pressurized Systems",
       NASA/TM-20210022275, 2021. https://ntrs.nasa.gov/citations/20210022275
.. [2] NASA-20190030443-PRIMING. A. LeClair, A. Boehm and A. Majumdar, "Fluid
       Transient Analysis of Propellant Feedlines during a Priming Event",
       TFAWS 2019. https://ntrs.nasa.gov/citations/20190030443
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from math import isfinite, pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, integer, real

WATER_HAMMER_SCHEMA = "rocket_propulsion_water_hammer_moc_v1"
WATER_HAMMER_REFERENCE_IDS = (
    "NASA-TM-20210022275-TRANSIENTS",
    "NASA-19980236009-PRESSURIZATION",
    "NASA-20190030443-PRIMING",
)

#: Diameter-to-thickness ratio at and below which the thick-walled boundary
#: factors apply (NASA/TM-20210022275, printed p. 19, below eq. 2).
THICK_WALL_DIAMETER_RATIO = 25.0
#: Upper bound on the number of reaches; keeps a run bounded in memory and time.
MAXIMUM_REACH_COUNT = 20_000
#: Upper bound on the number of characteristic time steps of one run.
MAXIMUM_TIME_STEPS = 2_000_000


def _positive(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and positive."""
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.", field=label)
    return number


def _nonnegative(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and non-negative."""
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.", field=label)
    return number


class PipeAnchoring(str, Enum):
    """Longitudinal restraint of the pipe, selecting the boundary factor ``C``.

    Members follow the four boundary cases of the NESC report [1]_: anchored
    at the upstream end only (case I), anchored against longitudinal movement
    (case II), with expansion joints throughout (case III), and rigid
    (case IV, ``C = 0``).

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed p. 20, Figures 4-5.
           https://ntrs.nasa.gov/citations/20210022275
    """

    UPSTREAM_ANCHORED = "upstream_anchored"
    NO_LONGITUDINAL_MOVEMENT = "no_longitudinal_movement"
    EXPANSION_JOINTS = "expansion_joints"
    RIGID = "rigid"


def pipe_boundary_factor(
    anchoring: PipeAnchoring,
    *,
    poisson_ratio: float,
    inner_diameter_m: float,
    wall_thickness_m: float,
) -> float:
    """Return the pipe boundary factor ``C`` of the wave-speed formula.

    Parameters
    ----------
    anchoring : PipeAnchoring
        Longitudinal restraint case.
    poisson_ratio : float
        Poisson's ratio ``mu`` of the pipe material, in ``[0, 0.5)``.
    inner_diameter_m : float
        Positive inner diameter ``D``, m.
    wall_thickness_m : float
        Positive wall thickness ``t``, m.

    Returns
    -------
    float
        Dimensionless boundary factor ``C``.

    Raises
    ------
    DomainError
        If an input is outside its domain or ``anchoring`` is not a
        `PipeAnchoring`.

    Notes
    -----
    For thin walls (``D / t > 25``) the factors are ``1 - mu/2``,
    ``1 - mu**2``, ``1`` and ``0`` for cases I to IV.  For ``D / t <= 25``
    the thick-walled forms ``(2t/D)(1 + mu) + D/(D + t) * k`` are used, with
    ``k`` equal to the corresponding thin-walled factor; the rigid case stays
    zero [1]_.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed pp. 19-20, Figure 5.
           https://ntrs.nasa.gov/citations/20210022275
    """

    if not isinstance(anchoring, PipeAnchoring):
        raise DomainError("Pipe anchoring must be a PipeAnchoring member.", field="anchoring")
    mu = _nonnegative(poisson_ratio, "Poisson ratio")
    if mu >= 0.5:
        raise DomainError("Poisson ratio must be below 0.5.", field="Poisson ratio")
    diameter = _positive(inner_diameter_m, "Pipe inner diameter")
    thickness = _positive(wall_thickness_m, "Pipe wall thickness")
    if anchoring is PipeAnchoring.RIGID:
        return 0.0
    thin = {
        PipeAnchoring.UPSTREAM_ANCHORED: 1.0 - 0.5 * mu,
        PipeAnchoring.NO_LONGITUDINAL_MOVEMENT: 1.0 - mu * mu,
        PipeAnchoring.EXPANSION_JOINTS: 1.0,
    }[anchoring]
    if diameter / thickness > THICK_WALL_DIAMETER_RATIO:
        return thin
    return 2.0 * thickness / diameter * (1.0 + mu) + diameter / (diameter + thickness) * thin


def pipe_wave_speed_m_s(
    *,
    bulk_modulus_pa: float,
    density_kg_m3: float,
    inner_diameter_m: float,
    wall_thickness_m: float,
    youngs_modulus_pa: float,
    poisson_ratio: float,
    anchoring: PipeAnchoring,
) -> float:
    """Return the pressure-wave speed in a liquid-filled elastic pipe.

    Parameters
    ----------
    bulk_modulus_pa : float
        Positive liquid bulk modulus ``beta``, Pa.
    density_kg_m3 : float
        Positive liquid density ``rho``, kg/m^3.
    inner_diameter_m : float
        Positive inner diameter ``D``, m.
    wall_thickness_m : float
        Positive wall thickness ``t``, m.
    youngs_modulus_pa : float
        Positive Young's modulus ``E`` of the pipe, Pa.
    poisson_ratio : float
        Poisson's ratio of the pipe material, in ``[0, 0.5)``.
    anchoring : PipeAnchoring
        Longitudinal restraint case selecting ``C``.

    Returns
    -------
    float
        Wave speed ``a``, m/s.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    The wave speed is [1]_

    ``a = sqrt(beta / rho) / sqrt(1 + beta D C / (E t))``

    with ``C`` from `pipe_boundary_factor`.  A rigid pipe (``C = 0``) gives
    the liquid sound speed ``sqrt(beta / rho)``.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed p. 19, eq. 2.
           https://ntrs.nasa.gov/citations/20210022275
    """

    beta = _positive(bulk_modulus_pa, "Liquid bulk modulus")
    rho = _positive(density_kg_m3, "Liquid density")
    diameter = _positive(inner_diameter_m, "Pipe inner diameter")
    thickness = _positive(wall_thickness_m, "Pipe wall thickness")
    modulus = _positive(youngs_modulus_pa, "Pipe Young's modulus")
    factor = pipe_boundary_factor(
        anchoring,
        poisson_ratio=poisson_ratio,
        inner_diameter_m=diameter,
        wall_thickness_m=thickness,
    )
    speed = sqrt(beta / rho) / sqrt(1.0 + (beta / modulus) * (diameter / thickness) * factor)
    return _positive(speed, "Pipe wave speed")


def joukowsky_surge_pa(
    *, density_kg_m3: float, wave_speed_m_s: float, velocity_change_m_s: float
) -> float:
    """Return the Joukowsky pressure change for an instantaneous velocity change.

    Parameters
    ----------
    density_kg_m3 : float
        Positive liquid density ``rho``, kg/m^3.
    wave_speed_m_s : float
        Positive wave speed ``a``, m/s.
    velocity_change_m_s : float
        Finite change ``dV = V_final - V_initial`` of the liquid velocity, m/s.

    Returns
    -------
    float
        Pressure change ``dP = -rho a dV``, Pa; positive for a deceleration
        such as a valve closure.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    The relation holds for changes faster than the wave travel time and is
    conservative for slower valve motion [1]_.  The same relation bounds the
    dead-end surge of a priming liquid front of known impact velocity [2]_.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed p. 19, eq. 1, and p. 54,
           eqs. A-1 to A-6. https://ntrs.nasa.gov/citations/20210022275
    .. [2] NASA-20190030443-PRIMING. A. LeClair, A. Boehm and A. Majumdar,
           "Fluid Transient Analysis of Propellant Feedlines during a Priming
           Event", TFAWS 2019, slide 3.
           https://ntrs.nasa.gov/citations/20190030443
    """

    rho = _positive(density_kg_m3, "Liquid density")
    speed = _positive(wave_speed_m_s, "Wave speed")
    change = real(velocity_change_m_s, "Velocity change")
    return real(-rho * speed * change, "Joukowsky surge")


@dataclass(frozen=True, slots=True)
class WaterHammerPipe:
    """Uniform liquid-filled pipe for the characteristic solution.

    Parameters
    ----------
    length_m : float
        Positive pipe length ``L`` from reservoir to valve, m.
    inner_diameter_m : float
        Positive inner diameter ``D``, m.
    wave_speed_m_s : float
        Positive pressure-wave speed ``a``, m/s; see `pipe_wave_speed_m_s`.
    liquid_density_kg_m3 : float
        Positive constant liquid density ``rho``, kg/m^3.
    darcy_friction_factor : float
        Non-negative steady Darcy friction factor ``f``.
    axial_acceleration_m_s2 : float, default 0.0
        Finite body-force acceleration component ``alpha sin(theta)`` along
        the flow direction, m/s^2, as it appears in the momentum equation.

    Raises
    ------
    DomainError
        If a value is outside its domain.

    Notes
    -----
    The governing equations are [1]_

    ``dp/dt + (rho a**2 / A) dQ/dx = 0``

    ``(1/A) dQ/dt + (1/rho) dp/dx + alpha sin(theta) + f Q abs(Q) / (2 D A**2) = 0``

    with volume flow ``Q`` and flow area ``A = pi D**2 / 4``.  Friction is the
    steady Darcy law evaluated with the flow at the foot of each
    characteristic.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed p. 55, eqs. B-1 and B-2.
           https://ntrs.nasa.gov/citations/20210022275
    """

    length_m: float
    inner_diameter_m: float
    wave_speed_m_s: float
    liquid_density_kg_m3: float
    darcy_friction_factor: float
    axial_acceleration_m_s2: float = 0.0

    def __post_init__(self) -> None:
        coerce_fields(self, "Water-hammer pipe")
        _positive(self.length_m, "Pipe length")
        _positive(self.inner_diameter_m, "Pipe inner diameter")
        _positive(self.wave_speed_m_s, "Wave speed")
        _positive(self.liquid_density_kg_m3, "Liquid density")
        _nonnegative(self.darcy_friction_factor, "Darcy friction factor")
        _positive(self.flow_area_m2, "Pipe flow area")

    @property
    def flow_area_m2(self) -> float:
        """Circular flow area ``A = pi D**2 / 4``, m^2 [1]_.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 55.
               https://ntrs.nasa.gov/citations/20210022275
        """

        return 0.25 * pi * self.inner_diameter_m * self.inner_diameter_m


@dataclass(frozen=True, slots=True)
class LinearValveLaw:
    """Piecewise-linear relative valve flow coefficient ``tau(t)``.

    Parameters
    ----------
    start_time_s : float
        Non-negative time at which the valve starts to move, s.
    stroke_time_s : float
        Non-negative stroke duration, s; zero gives a step change just after
        ``start_time_s``.
    initial_opening : float, default 1.0
        Opening before the stroke, in ``[0, 1]``.
    final_opening : float, default 0.0
        Opening after the stroke, in ``[0, 1]``.

    Raises
    ------
    DomainError
        If a value is outside its domain.

    Notes
    -----
    ``tau`` is the valve flow at a fixed pressure drop relative to the fully
    open valve, so the valve resistance is ``K_open / tau**2``.  Linear area
    laws are among the opening profiles compared with priming experiments in
    [1]_; closure faster than the wave travel time produces the Joukowsky
    surge and slower closure a smaller one [2]_.

    References
    ----------
    .. [1] NASA-20190030443-PRIMING. A. LeClair, A. Boehm and A. Majumdar,
           "Fluid Transient Analysis of Propellant Feedlines during a Priming
           Event", TFAWS 2019, slide 17.
           https://ntrs.nasa.gov/citations/20190030443
    .. [2] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed pp. 13 and 24.
           https://ntrs.nasa.gov/citations/20210022275
    """

    start_time_s: float
    stroke_time_s: float
    initial_opening: float = 1.0
    final_opening: float = 0.0

    def __post_init__(self) -> None:
        coerce_fields(self, "Valve law")
        _nonnegative(self.start_time_s, "Valve stroke start time")
        _nonnegative(self.stroke_time_s, "Valve stroke time")
        for value, label in (
            (self.initial_opening, "Initial valve opening"),
            (self.final_opening, "Final valve opening"),
        ):
            if not 0.0 <= value <= 1.0:
                raise DomainError(f"{label} must be in [0, 1].", field=label)

    def __call__(self, time_s: float) -> float:
        """Return the opening ``tau`` at ``time_s`` [1]_.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 13.
               https://ntrs.nasa.gov/citations/20210022275
        """

        if time_s <= self.start_time_s:
            return self.initial_opening
        if self.stroke_time_s == 0.0 or time_s >= self.start_time_s + self.stroke_time_s:
            return self.final_opening
        fraction = (time_s - self.start_time_s) / self.stroke_time_s
        return self.initial_opening + (self.final_opening - self.initial_opening) * fraction


@dataclass(frozen=True, slots=True)
class WaterHammerConfiguration:
    """Reservoir-pipe-valve system and run control.

    Parameters
    ----------
    pipe : WaterHammerPipe
        Uniform pipe and liquid.
    reservoir_pressure_pa : float
        Positive constant pressure at the pipe inlet, Pa.
    downstream_pressure_pa : float
        Positive constant pressure downstream of the valve, Pa.
    valve_full_open_resistance_pa_s2_m6 : float
        Positive quadratic valve resistance ``K_open`` at ``tau = 1``, so
        that ``p_valve - p_down = (K_open / tau**2) Q abs(Q)``, Pa s^2/m^6.
    valve_opening : Callable[[float], float]
        Relative valve flow coefficient ``tau(t)`` in ``[0, 1]``, for example
        a `LinearValveLaw`.
    reach_count : int
        Number of equal reaches ``N`` (1 to `MAXIMUM_REACH_COUNT`).
    duration_s : float
        Positive simulated time, s.
    station_indices : tuple of int, optional
        Node indices ``0..N`` whose histories are recorded; default is
        reservoir, mid-pipe and valve.
    vapor_pressure_pa : float, default 0.0
        Non-negative liquid vapour pressure, Pa.  A computed pressure at or
        below it means column separation, which is not modelled.
    schema : str, default WATER_HAMMER_SCHEMA
        Configuration schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If a value is outside its domain, the opening law is not callable,
        a station index is out of range or repeated, or the run would exceed
        `MAXIMUM_TIME_STEPS`.

    Notes
    -----
    The quadratic loss law follows the network-code branch resistance
    ``K_f mdot abs(mdot)`` [2]_, written here for volume flow.  The grid is
    ``dx = L / N`` and ``dt = dx / a``, so the Courant number is exactly one
    and the characteristics pass through grid nodes [1]_.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed pp. 55-57, Figure B-1.
           https://ntrs.nasa.gov/citations/20210022275
    .. [2] NASA-19980236009-PRESSURIZATION. A. Majumdar and T. Steadman,
           "Numerical Modeling of Pressurization of a Propellant Tank",
           AIAA 99-0879 draft, p. 3, eq. 2.
           https://ntrs.nasa.gov/citations/19980236009
    """

    pipe: WaterHammerPipe
    reservoir_pressure_pa: float
    downstream_pressure_pa: float
    valve_full_open_resistance_pa_s2_m6: float
    valve_opening: Callable[[float], float]
    reach_count: int
    duration_s: float
    station_indices: tuple[int, ...] | None = None
    vapor_pressure_pa: float = 0.0
    schema: str = WATER_HAMMER_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Water-hammer configuration")
        if self.schema != WATER_HAMMER_SCHEMA:
            raise DomainError("Unsupported water-hammer schema.", field="schema")
        if not isinstance(self.pipe, WaterHammerPipe):
            raise DomainError("Water-hammer pipe must be a WaterHammerPipe.", field="pipe")
        _positive(self.reservoir_pressure_pa, "Reservoir pressure")
        _positive(self.downstream_pressure_pa, "Downstream pressure")
        _positive(self.valve_full_open_resistance_pa_s2_m6, "Valve full-open resistance")
        if not callable(self.valve_opening):
            raise DomainError("Valve opening law must be callable.", field="valve_opening")
        if not 1 <= self.reach_count <= MAXIMUM_REACH_COUNT:
            raise DomainError(
                f"Reach count must be between 1 and {MAXIMUM_REACH_COUNT}.", field="reach_count")
        _positive(self.duration_s, "Water-hammer duration")
        _nonnegative(self.vapor_pressure_pa, "Vapour pressure")
        if self.step_count > MAXIMUM_TIME_STEPS:
            raise DomainError("Water-hammer run exceeds the time-step limit.", field="duration_s")
        if self.station_indices is None:
            stations = (0, self.reach_count // 2, self.reach_count)
            stations = tuple(dict.fromkeys(stations))
        else:
            if not isinstance(self.station_indices, (list, tuple)) or not self.station_indices:
                raise DomainError("Station indices must be a non-empty sequence.",
                                  field="station_indices")
            stations = tuple(integer(item, "Station index") for item in self.station_indices)
            if len(set(stations)) != len(stations):
                raise DomainError("Station indices must be unique.", field="station_indices")
            if any(not 0 <= item <= self.reach_count for item in stations):
                raise DomainError("Station index outside 0..reach_count.", field="station_indices")
        object.__setattr__(self, "station_indices", stations)

    @property
    def reach_length_m(self) -> float:
        """Reach length ``dx = L / N``, m [1]_.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 56.
               https://ntrs.nasa.gov/citations/20210022275
        """

        return self.pipe.length_m / self.reach_count

    @property
    def time_step_s(self) -> float:
        """Courant-number-one step ``dt = dx / a``, s [1]_.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 56, eq. B-8.
               https://ntrs.nasa.gov/citations/20210022275
        """

        return self.reach_length_m / self.pipe.wave_speed_m_s

    @property
    def step_count(self) -> int:
        """Number of characteristic steps covering ``duration_s`` (rounded up) [1]_.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 56.
               https://ntrs.nasa.gov/citations/20210022275
        """

        ratio = self.duration_s / self.time_step_s
        count = int(ratio)
        return count if count >= ratio else count + 1


@dataclass(frozen=True, slots=True)
class WaterHammerResult:
    """Histories, peak values and balance diagnostics of one MOC run.

    Attributes
    ----------
    time_s : tuple of float
        Sample times ``n dt``, s, starting at zero.
    station_positions_m : tuple of float
        Distance of each recorded station from the reservoir, m.
    station_pressure_pa : tuple of tuple of float
        Pressure history per station, Pa.
    station_volume_flow_m3_s : tuple of tuple of float
        Volume-flow history per station, m^3/s.
    initial_volume_flow_m3_s : float
        Steady initial flow ``Q0``, m^3/s.
    initial_velocity_m_s : float
        ``V0 = Q0 / A``, m/s.
    joukowsky_surge_pa : float
        ``rho a V0``: the surge of an instantaneous full closure [1]_, Pa.
    valve_initial_pressure_pa : float
        Steady valve-node pressure before the transient, Pa.
    valve_peak_surge_pa : float
        Maximum valve-node pressure minus its initial value, Pa.
    peak_pressure_pa : float
        Maximum pressure over all nodes and times, Pa.
    peak_pressure_time_s : float
        Time of the peak, s.
    peak_pressure_position_m : float
        Position of the peak from the reservoir, m.
    minimum_pressure_pa : float
        Minimum pressure over all nodes and times, Pa.
    time_step_s : float
        Characteristic time step, s.
    reach_length_m : float
        Reach length, m.
    maximum_mass_balance_residual_kg : float
        Largest absolute residual of the stored-liquid mass balance, kg.
    final_mass_balance_residual_kg : float
        Residual at the final time, kg.
    initial_acoustic_energy_j : float
        Kinetic plus elastic energy at ``t = 0``, J.
    final_acoustic_energy_j : float
        Same at the final time, J.
    friction_dissipation_j : float
        Integrated Darcy dissipation, J.
    boundary_work_j : float
        Integrated acoustic power entering through both pipe ends, J.
    body_force_work_j : float
        Integrated work against the axial body force, J.
    energy_balance_residual_j : float
        ``E_final - E_initial - boundary_work + dissipation + body_work``, J.

    Notes
    -----
    Stored volume follows from the continuity equation as
    ``(A / (rho a**2)) * integral(p dx)``; its change must equal the net
    inflow volume.  Multiplying the momentum equation by ``rho Q`` and the
    continuity equation by ``A p' / (rho a**2)``, with ``p' = p - p_res``,
    gives the acoustic energy balance whose residual is reported [1]_.
    Spatial integrals use the trapezoidal rule over nodes and time integrals
    the trapezoidal rule over steps, so residuals are discretization
    diagnostics, not exact invariants.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed p. 19 eq. 1 and p. 55
           eqs. B-1, B-2. https://ntrs.nasa.gov/citations/20210022275
    """

    time_s: tuple[float, ...]
    station_positions_m: tuple[float, ...]
    station_pressure_pa: tuple[tuple[float, ...], ...]
    station_volume_flow_m3_s: tuple[tuple[float, ...], ...]
    initial_volume_flow_m3_s: float
    initial_velocity_m_s: float
    joukowsky_surge_pa: float
    valve_initial_pressure_pa: float
    valve_peak_surge_pa: float
    peak_pressure_pa: float
    peak_pressure_time_s: float
    peak_pressure_position_m: float
    minimum_pressure_pa: float
    time_step_s: float
    reach_length_m: float
    maximum_mass_balance_residual_kg: float
    final_mass_balance_residual_kg: float
    initial_acoustic_energy_j: float
    final_acoustic_energy_j: float
    friction_dissipation_j: float
    boundary_work_j: float
    body_force_work_j: float
    energy_balance_residual_j: float
    reference_ids: tuple[str, ...] = field(default=WATER_HAMMER_REFERENCE_IDS)

    @property
    def surge_to_joukowsky_ratio(self) -> float:
        """Valve peak surge divided by ``rho a V0`` [1]_; ``nan`` if ``V0 = 0``.

        References
        ----------
        .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment
               of Transient Pressure Events in Space Flight Pressurized
               Systems", NASA/TM-20210022275, 2021, printed p. 19, eq. 1.
               https://ntrs.nasa.gov/citations/20210022275
        """

        if self.joukowsky_surge_pa == 0.0:
            return float("nan")
        return self.valve_peak_surge_pa / self.joukowsky_surge_pa


def _opening(law: Callable[[float], float], time_s: float) -> float:
    """Evaluate the valve law and require a finite value in ``[0, 1]``."""
    try:
        value = law(time_s)
    except DomainError:
        raise
    except (TypeError, ValueError, ArithmeticError) as error:
        raise DomainError("Valve opening law failed.", field="valve_opening") from error
    value = real(value, "Valve opening")
    if not 0.0 <= value <= 1.0:
        raise DomainError("Valve opening must be in [0, 1].", field="valve_opening")
    return value


def _valve_flow(c_plus: float, impedance: float, downstream: float, resistance: float,
                opening: float) -> float:
    """Solve ``Cp - B Q - p_down = (K/tau**2) Q abs(Q)`` for ``Q`` without cancellation."""
    if opening == 0.0:
        return 0.0
    drive = c_plus - downstream
    k = resistance / (opening * opening)
    return 2.0 * drive / (impedance + sqrt(impedance * impedance + 4.0 * k * abs(drive)))


def _trapezoid(values: list[float], spacing: float) -> float:
    """Trapezoidal integral of equally spaced samples."""
    if len(values) == 1:
        return 0.0
    return spacing * (sum(values) - 0.5 * (values[0] + values[-1]))


def simulate_water_hammer(configuration: WaterHammerConfiguration) -> WaterHammerResult:
    """Integrate the reservoir-pipe-valve transient by the method of characteristics.

    Parameters
    ----------
    configuration : WaterHammerConfiguration
        System and run control.

    Returns
    -------
    WaterHammerResult
        Station histories, peaks, Joukowsky comparison and diagnostics.

    Raises
    ------
    DomainError
        If the configuration is invalid, the valve law leaves ``[0, 1]``, a
        computed value is not finite, or a pressure reaches the vapour
        pressure (column separation is outside the model).

    Notes
    -----
    With impedance ``B = rho a / A``, reach friction
    ``R = rho f dx / (2 A**2 D)`` and body term ``S alpha = rho dx alpha sin(theta)``,
    the interior update is [1]_

    ``Cp = p[i-1] + (B - R abs(Q[i-1])) Q[i-1] - S alpha``

    ``Cm = p[i+1] - (B - R abs(Q[i+1])) Q[i+1] + S alpha``

    ``p = (Cp + Cm) / 2``,  ``Q = (Cp - Cm) / (2 B)``.

    The reservoir node keeps ``p = p_res`` and takes ``Q = (p_res - Cm) / B``.
    The valve node solves ``Cp - B Q - p_down = (K_open / tau**2) Q abs(Q)``;
    a closed valve gives ``Q = 0``, ``p = Cp``.  The initial state is the
    steady solution of the same discrete equations at ``tau(0)``, so an
    unchanging valve leaves the state exactly steady.

    References
    ----------
    .. [1] NASA-TM-20210022275-TRANSIENTS. K. Imtiaz et al., "Treatment of
           Transient Pressure Events in Space Flight Pressurized Systems",
           NASA/TM-20210022275, 2021, printed pp. 55-57, eqs. B-11 to B-17.
           https://ntrs.nasa.gov/citations/20210022275
    """

    if not isinstance(configuration, WaterHammerConfiguration):
        raise DomainError("Configuration must be a WaterHammerConfiguration.")
    pipe = configuration.pipe
    n = configuration.reach_count
    dx = configuration.reach_length_m
    dt = configuration.time_step_s
    rho = pipe.liquid_density_kg_m3
    a = pipe.wave_speed_m_s
    area = pipe.flow_area_m2
    impedance = rho * a / area
    friction = rho * pipe.darcy_friction_factor * dx / (2.0 * area * area * pipe.inner_diameter_m)
    body = rho * dx * pipe.axial_acceleration_m_s2
    p_res = configuration.reservoir_pressure_pa
    p_down = configuration.downstream_pressure_pa
    k_open = configuration.valve_full_open_resistance_pa_s2_m6

    tau0 = _opening(configuration.valve_opening, 0.0)
    drive0 = p_res - n * body - p_down
    if tau0 == 0.0:
        q0 = 0.0
    else:
        total = n * friction + k_open / (tau0 * tau0)
        q0 = (1.0 if drive0 >= 0.0 else -1.0) * sqrt(abs(drive0) / total)
    pressure = [p_res - i * (friction * q0 * abs(q0) + body) for i in range(n + 1)]
    flow = [q0] * (n + 1)
    if min(pressure) <= configuration.vapor_pressure_pa:
        raise DomainError("Initial steady pressure reaches the vapour pressure.")

    stored_factor = area / (rho * a * a)
    elastic_factor = area / (2.0 * rho * a * a)
    kinetic_factor = rho / (2.0 * area)
    dissipation_factor = rho * pipe.darcy_friction_factor / (
        2.0 * pipe.inner_diameter_m * area * area)
    body_factor = rho * pipe.axial_acceleration_m_s2

    def energy(p_nodes: list[float], q_nodes: list[float]) -> float:
        return _trapezoid([kinetic_factor * q * q + elastic_factor * (p - p_res) ** 2
                           for p, q in zip(p_nodes, q_nodes, strict=True)], dx)

    def rates(p_nodes: list[float], q_nodes: list[float]) -> tuple[float, float, float, float]:
        net_in = q_nodes[0] - q_nodes[-1]
        boundary = (p_nodes[0] - p_res) * q_nodes[0] - (p_nodes[-1] - p_res) * q_nodes[-1]
        dissipation = _trapezoid([dissipation_factor * abs(q) ** 3 for q in q_nodes], dx)
        body_work = _trapezoid([body_factor * q for q in q_nodes], dx)
        return net_in, boundary, dissipation, body_work

    stations = configuration.station_indices
    times = [0.0]
    station_p = [[pressure[i]] for i in stations]
    station_q = [[flow[i]] for i in stations]
    stored0 = stored_factor * _trapezoid(pressure, dx)
    energy0 = energy(pressure, flow)
    previous_rates = rates(pressure, flow)
    inflow_volume = boundary_work = dissipation = body_work = 0.0
    max_mass_residual = 0.0
    valve_initial = pressure[n]
    valve_peak = valve_initial
    peak, peak_time, peak_node = max(pressure), 0.0, pressure.index(max(pressure))
    minimum = min(pressure)
    mass_residual = 0.0

    for step in range(1, configuration.step_count + 1):
        time = step * dt
        new_p = [0.0] * (n + 1)
        new_q = [0.0] * (n + 1)
        for i in range(1, n):
            q_left, q_right = flow[i - 1], flow[i + 1]
            c_plus = pressure[i - 1] + (impedance - friction * abs(q_left)) * q_left - body
            c_minus = pressure[i + 1] - (impedance - friction * abs(q_right)) * q_right + body
            new_p[i] = 0.5 * (c_plus + c_minus)
            new_q[i] = (c_plus - c_minus) / (2.0 * impedance)
        q1 = flow[1]
        c_minus = pressure[1] - (impedance - friction * abs(q1)) * q1 + body
        new_p[0] = p_res
        new_q[0] = (p_res - c_minus) / impedance
        q_last = flow[n - 1]
        c_plus = pressure[n - 1] + (impedance - friction * abs(q_last)) * q_last - body
        tau = _opening(configuration.valve_opening, time)
        valve_q = _valve_flow(c_plus, impedance, p_down, k_open, tau)
        new_q[n] = valve_q
        new_p[n] = c_plus - impedance * valve_q
        if not all(isfinite(v) for v in new_p) or not all(isfinite(v) for v in new_q):
            raise DomainError("Water-hammer state became non-finite.")
        pressure, flow = new_p, new_q

        current_rates = rates(pressure, flow)
        inflow_volume += 0.5 * dt * (previous_rates[0] + current_rates[0])
        boundary_work += 0.5 * dt * (previous_rates[1] + current_rates[1])
        dissipation += 0.5 * dt * (previous_rates[2] + current_rates[2])
        body_work += 0.5 * dt * (previous_rates[3] + current_rates[3])
        previous_rates = current_rates
        stored = stored_factor * _trapezoid(pressure, dx)
        mass_residual = rho * (stored - stored0 - inflow_volume)
        max_mass_residual = max(max_mass_residual, abs(mass_residual))

        step_max = max(pressure)
        if step_max > peak:
            peak, peak_time, peak_node = step_max, time, pressure.index(step_max)
        minimum = min(minimum, min(pressure))
        if minimum <= configuration.vapor_pressure_pa:
            raise DomainError(
                "Pressure reached the vapour pressure: column separation is not modelled.")
        valve_peak = max(valve_peak, pressure[n])
        times.append(time)
        for series, index in zip(station_p, stations, strict=True):
            series.append(pressure[index])
        for series, index in zip(station_q, stations, strict=True):
            series.append(flow[index])

    final_energy = energy(pressure, flow)
    velocity0 = q0 / area
    return WaterHammerResult(
        time_s=tuple(times),
        station_positions_m=tuple(index * dx for index in stations),
        station_pressure_pa=tuple(tuple(series) for series in station_p),
        station_volume_flow_m3_s=tuple(tuple(series) for series in station_q),
        initial_volume_flow_m3_s=q0,
        initial_velocity_m_s=velocity0,
        joukowsky_surge_pa=rho * a * velocity0,
        valve_initial_pressure_pa=valve_initial,
        valve_peak_surge_pa=valve_peak - valve_initial,
        peak_pressure_pa=peak,
        peak_pressure_time_s=peak_time,
        peak_pressure_position_m=peak_node * dx,
        minimum_pressure_pa=minimum,
        time_step_s=dt,
        reach_length_m=dx,
        maximum_mass_balance_residual_kg=max_mass_residual,
        final_mass_balance_residual_kg=mass_residual,
        initial_acoustic_energy_j=energy0,
        final_acoustic_energy_j=final_energy,
        friction_dissipation_j=dissipation,
        boundary_work_j=boundary_work,
        body_force_work_j=body_work,
        energy_balance_residual_j=final_energy - energy0 - boundary_work + dissipation + body_work,
    )
