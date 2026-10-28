"""Back-pressure-driven converging-diverging nozzle operating regimes.

A fixed converging-diverging nozzle with perfect-gas, quasi-one-dimensional,
adiabatic and inviscid flow is solved at an imposed back pressure. The regime
follows from comparing the back pressure with three analytical thresholds:
the just-choked subsonic exit pressure, the pressure behind a normal shock
standing at the exit plane, and the design exit pressure [1]_, [2]_. Flow
separation is screened with an empirical exit-to-ambient pressure criterion
[3]_; it is not solved. Altitude sweeps use the U.S. Standard Atmosphere 1976
[4]_.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
.. [3] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976, p. 17.
       https://ntrs.nasa.gov/citations/19710021390
.. [4] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
       NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log, log1p, sqrt

from rocket_propulsion.core.constants import (
    STANDARD_GRAVITY_M_S2,
    USSA1976_GAS_CONSTANT_J_MOL_K,
    USSA1976_MOLAR_MASS_KG_MOL,
    USSA1976_SEA_LEVEL_PRESSURE_PA,
    USSA1976_SEA_LEVEL_TEMPERATURE_K,
)
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .isentropic import calculate_isentropic, mach_from_isentropic
from .normal_shock import calculate_normal_shock
from .nozzle import calculate_ideal_nozzle


@dataclass(frozen=True, slots=True)
class NozzleRegimeThresholds:
    """Back-pressure boundaries for a fixed perfect-gas nozzle [1]_, [2]_.

    Attributes
    ----------
    chamber_pressure_pa : float
        Chamber stagnation pressure in Pa.
    just_choked_back_pressure_pa : float
        Back pressure at which the throat first chokes with a fully subsonic
        diverging section (subsonic isentropic branch at the exit).
    shock_at_exit_back_pressure_pa : float
        Back pressure that places a normal shock exactly at the exit plane.
    design_exit_pressure_pa : float
        Exit pressure on the supersonic isentropic branch (shock-free design
        condition).

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    chamber_pressure_pa: float
    just_choked_back_pressure_pa: float
    shock_at_exit_back_pressure_pa: float
    design_exit_pressure_pa: float


@dataclass(frozen=True, slots=True)
class NozzleFlowStation:
    """One normalized station used by the off-design contour plot.

    Station Mach numbers follow the isentropic area relation [1]_.

    Attributes
    ----------
    x_fraction : float
        Normalized position from throat (0) to exit (1).
    area_ratio : float
        Local area ratio ``A/A_t``, linear in ``x_fraction``.
    mach, pressure_pa, total_pressure_pa : float
        Local Mach number, static pressure and stagnation pressure.
    region : str
        ``"subsonic"``, ``"supersonic"``, ``"pre-shock supersonic"`` or
        ``"post-shock subsonic"``.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    x_fraction: float
    area_ratio: float
    mach: float
    pressure_pa: float
    total_pressure_pa: float
    region: str


@dataclass(frozen=True, slots=True)
class ShockLocation:
    """Solved normal-shock state within the divergent section [1]_.

    Attributes
    ----------
    area_ratio : float
        Area ratio ``A_s/A_t`` at the shock.
    x_fraction : float
        Normalized shock position from throat (0) to exit (1).
    upstream_mach, downstream_mach : float
        Mach numbers immediately before and after the shock.
    stagnation_pressure_ratio : float
        Total-pressure recovery across the shock.
    downstream_exit_mach : float
        Subsonic exit Mach number behind the shock.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    area_ratio: float
    x_fraction: float
    upstream_mach: float
    downstream_mach: float
    stagnation_pressure_ratio: float
    downstream_exit_mach: float


@dataclass(frozen=True, slots=True)
class SeparationAssessment:
    """Empirical separation screening with an explicit validity statement [1]_.

    Attributes
    ----------
    predicted : bool
        Whether either criterion indicates possible separation.
    correlation, criterion : str
        Name and inequality of the Summerfield screening correlation.
    exit_to_ambient_pressure_ratio : float or None
        ``p_exit/p_ambient``; ``None`` for vacuum back pressure.
    valid_for : str
        Stated validity envelope of the Summerfield correlation.
    source, source_url : str
        Source document of the Summerfield correlation.
    schmucker_predicted : bool
        Whether the Schmucker criterion places the separation point inside
        the divergent section.
    schmucker_separation_mach, schmucker_separation_area_ratio : float or None
        One-dimensional wall Mach number and area ratio of the incipient
        separation point ``i``; ``None`` when the flow stays attached or the
        criterion is not applicable.
    schmucker_separation_pressure_pa : float or None
        Wall pressure ``p_i`` at that point.
    schmucker_within_range : bool
        False when the search had to stop at ``M_i = 5``, the end of the data
        behind the criterion; separation beyond it is not assessed.

    Notes
    -----
    Summerfield's ``p_i/p_a = 0.4`` is the oldest criterion, a conservative
    guideline for ``p_c/p_a`` of 15 to 20 [1]_, [2]_. Schmucker fitted the
    averaged separation data of chemical engines with a Mach-number form,
    ``p_i/p_a = (1.88 M_i - 1)^(-0.64)``, which removes the a-priori
    isentropic exponent of pressure-ratio criteria and locates the point of
    incipient separation; it is not to be used above ``M_i = 5`` [2]_.

    References
    ----------
    .. [1] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976, p. 17.
           https://ntrs.nasa.gov/citations/19770009165
    .. [2] NASA-TM-77396-SCHMUCKER. R. H. Schmucker, "Flow Processes in
           Overexpanded Chemical Rocket Nozzles. Part 1: Flow Separation",
           NASA TM-77396, 1984, Eqs. (7) and (11), printed pp. 37-39.
           https://ntrs.nasa.gov/citations/19840011402
    """

    predicted: bool
    correlation: str
    criterion: str
    exit_to_ambient_pressure_ratio: float | None
    valid_for: str
    source: str
    source_url: str
    schmucker_predicted: bool = False
    schmucker_separation_mach: float | None = None
    schmucker_separation_area_ratio: float | None = None
    schmucker_separation_pressure_pa: float | None = None
    schmucker_within_range: bool = True


@dataclass(frozen=True, slots=True)
class OffDesignNozzleResult:
    """Complete one-dimensional nozzle state at one ambient pressure [1]_, [2]_.

    Attributes
    ----------
    regime : str
        One of ``"no-flow"``, ``"unchoked"``, ``"just-choked"``,
        ``"internal-normal-shock"``, ``"shock-at-exit"``, ``"overexpanded"``,
        ``"design"`` or ``"underexpanded"``.
    chamber_pressure_pa, ambient_pressure_pa : float
        Imposed chamber stagnation and ambient pressures in Pa.
    area_ratio, throat_area_m2, exit_area_m2 : float
        Nozzle geometry.
    throat_mach, exit_mach : float
        Mach numbers at throat and exit.
    exit_pressure_pa, exit_temperature_k, exit_velocity_m_s : float
        Exit-plane static state.
    exit_total_pressure_pa : float
        Exit stagnation pressure, reduced by any normal shock.
    mass_flow_kg_s : float
        Mass flow set at the throat (or at the exit when unchoked).
    exit_mass_flow_kg_s : float
        Mass flow recomputed at the exit plane.
    mass_flow_relative_residual : float
        ``|mdot_exit - mdot| / mdot`` as a continuity check; zero for no flow.
    momentum_thrust_n, pressure_thrust_n, thrust_n : float
        ``mdot V_e``, ``(p_e - p_a) A_e`` and their sum, in N.
    specific_impulse_s : float or None
        ``F / (mdot g0)``; ``None`` when there is no flow.
    thresholds : NozzleRegimeThresholds
        Back-pressure boundaries used for the classification.
    shock : ShockLocation or None
        Normal-shock state when a shock stands in or at the exit of the
        divergent section.
    separation : SeparationAssessment
        Empirical separation screening.
    warnings : tuple of WarningMessage
        Model-limit warnings for this state.
    profile : tuple of NozzleFlowStation
        Station profile from throat to exit; empty for no flow.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    regime: str
    chamber_pressure_pa: float
    ambient_pressure_pa: float
    area_ratio: float
    throat_area_m2: float
    exit_area_m2: float
    throat_mach: float
    exit_mach: float
    exit_pressure_pa: float
    exit_temperature_k: float
    exit_velocity_m_s: float
    exit_total_pressure_pa: float
    mass_flow_kg_s: float
    exit_mass_flow_kg_s: float
    mass_flow_relative_residual: float
    momentum_thrust_n: float
    pressure_thrust_n: float
    thrust_n: float
    specific_impulse_s: float | None
    thresholds: NozzleRegimeThresholds
    shock: ShockLocation | None
    separation: SeparationAssessment
    warnings: tuple[WarningMessage, ...]
    profile: tuple[NozzleFlowStation, ...]


@dataclass(frozen=True, slots=True)
class NozzleAltitudePoint:
    """One atmosphere/nozzle result in an altitude sweep.

    Ambient pressure is the U.S. Standard Atmosphere 1976 value [1]_.

    Attributes
    ----------
    altitude_m, ambient_pressure_pa : float
        Geopotential altitude and standard-atmosphere pressure.
    regime : str
        Operating regime from :func:`calculate_off_design_nozzle`.
    thrust_n : float
        Thrust in N.
    specific_impulse_s : float or None
        Specific impulse; ``None`` when there is no flow.
    chamber_to_ambient_pressure_ratio : float or None
        ``p_c/p_a``; ``None`` for zero ambient pressure.
    separation_predicted : bool
        Separation screening flag.
    shock_x_fraction : float or None
        Normalized shock position, if a shock is present.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    """

    altitude_m: float
    ambient_pressure_pa: float
    regime: str
    thrust_n: float
    specific_impulse_s: float | None
    chamber_to_ambient_pressure_ratio: float | None
    separation_predicted: bool
    shock_x_fraction: float | None


@dataclass(frozen=True, slots=True)
class NozzleAltitudeSweep:
    """Off-design performance across the 1976 standard atmosphere [1]_.

    Attributes
    ----------
    points : tuple of NozzleAltitudePoint
        Samples in increasing altitude.
    atmosphere_model : str
        Name of the atmosphere model used.
    warning : str
        Statement of the separation-screening limitation.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    """

    points: tuple[NozzleAltitudePoint, ...]
    atmosphere_model: str
    warning: str


def _validate_inputs(
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    throat_area_m2: float,
    area_ratio: float,
    ambient_pressure_pa: float,
    gamma: float,
    gas_constant_j_kg_k: float,
) -> None:
    """Raise :class:`DomainError` with the offending field for invalid nozzle inputs."""

    positive = (
        (chamber_pressure_pa, "chamber_pressure_pa"),
        (chamber_temperature_k, "chamber_temperature_k"),
        (throat_area_m2, "throat_area_m2"),
        (area_ratio, "area_ratio"),
        (gamma, "gamma"),
        (gas_constant_j_kg_k, "gas_constant_j_kg_k"),
    )
    for value, field in positive:
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if area_ratio < 1.0:
        raise DomainError("Area ratio must be at least one.", field="area_ratio")
    if gamma <= 1.0:
        raise DomainError("Gamma must exceed one.", field="gamma")
    if not isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0.0:
        raise DomainError(
            "Ambient pressure must be finite and non-negative.",
            field="ambient_pressure_pa",
        )


def nozzle_regime_thresholds(
    *, chamber_pressure_pa: float, area_ratio: float, gamma: float
) -> NozzleRegimeThresholds:
    """Return analytical back-pressure boundaries for a fixed nozzle.

    Parameters
    ----------
    chamber_pressure_pa : float
        Chamber stagnation pressure in Pa.
    area_ratio : float
        Exit-to-throat area ratio, at least one.
    gamma : float
        Ratio of specific heats, greater than one.

    Returns
    -------
    NozzleRegimeThresholds
        Just-choked, shock-at-exit and design exit pressures.

    Raises
    ------
    DomainError
        If ``area_ratio`` or ``gamma`` is outside the isentropic domain.

    Notes
    -----
    The just-choked pressure is ``p_c (p/p0)`` at the subsonic solution of
    ``A/A* = area_ratio``; the design pressure is ``p_c (p/p0)`` at the
    supersonic solution; the shock-at-exit pressure multiplies the design
    pressure by the normal-shock static-pressure ratio at the supersonic exit
    Mach number [1]_, [2]_. Inputs are not validated here beyond the checks of
    the isentropic inversion.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    exit_subsonic_mach = mach_from_isentropic("area_ratio_subsonic", area_ratio, gamma)
    exit_supersonic_mach = mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
    subsonic = calculate_isentropic(exit_subsonic_mach, gamma)
    supersonic = calculate_isentropic(exit_supersonic_mach, gamma)
    shock = calculate_normal_shock(exit_supersonic_mach, gamma)
    design_exit_pressure = chamber_pressure_pa * supersonic.pressure_total_ratio
    return NozzleRegimeThresholds(
        chamber_pressure_pa=chamber_pressure_pa,
        just_choked_back_pressure_pa=(
            chamber_pressure_pa * subsonic.pressure_total_ratio
        ),
        shock_at_exit_back_pressure_pa=design_exit_pressure * shock.pressure_ratio,
        design_exit_pressure_pa=design_exit_pressure,
    )


def _choked_mass_flow(
    pressure_pa: float,
    temperature_k: float,
    area_m2: float,
    gamma: float,
    gas_constant: float,
) -> float:
    """Return the choked mass flow through ``area_m2`` from stagnation conditions.

    Notes
    -----
    ``mdot = p0 A / sqrt(T0) sqrt(gamma/R) (2/(gamma+1))^((gamma+1)/(2(gamma-1)))``
    for a perfect gas with a sonic section of area ``A`` [1]_, [2]_.

    References
    ----------
    .. [1] NASA-GLENN-MASS-FLOW. NASA Glenn Research Center, isentropic mass-flow
           conservation and density/velocity context.
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    return (
        pressure_pa
        * area_m2
        / sqrt(temperature_k)
        * sqrt(gamma / gas_constant)
        * (2.0 / (gamma + 1.0))
        ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )


def _mass_flow_at_station(
    *,
    total_pressure_pa: float,
    total_temperature_k: float,
    area_m2: float,
    mach: float,
    gamma: float,
    gas_constant: float,
) -> float:
    """Return ``rho V A`` at a station of known Mach number and stagnation state.

    Notes
    -----
    Static pressure and temperature follow from the isentropic ratios [1]_;
    density from the perfect-gas law and velocity from
    ``V = M sqrt(gamma R T)`` [2]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    .. [2] NASA-GLENN-MASS-FLOW. NASA Glenn Research Center, isentropic mass-flow
           conservation and density/velocity context.
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/
    """

    ratios = calculate_isentropic(mach, gamma)
    pressure = total_pressure_pa * ratios.pressure_total_ratio
    temperature = total_temperature_k * ratios.temperature_total_ratio
    density = pressure / (gas_constant * temperature)
    velocity = mach * sqrt(gamma * gas_constant * temperature)
    return density * velocity * area_m2


def _shock_exit_state(
    shock_area_ratio: float,
    *,
    exit_area_ratio: float,
    chamber_pressure_pa: float,
    gamma: float,
) -> tuple[float, float, float, float, float]:
    """Return the exit state for a normal shock at area ratio ``shock_area_ratio``.

    Returns
    -------
    tuple of float
        ``(exit_pressure_pa, upstream_mach, downstream_mach,
        stagnation_pressure_ratio, exit_mach)``.

    Notes
    -----
    Upstream of the shock the flow is on the supersonic isentropic branch.
    Behind it the sonic reference area grows by ``1/(p02/p01)``, so the exit
    sees ``A_e/A2* = (A_e/A_s) (A_s/A2*)`` and is solved on the subsonic
    branch with the reduced stagnation pressure [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    upstream_mach = mach_from_isentropic(
        "area_ratio_supersonic", shock_area_ratio, gamma
    )
    shock = calculate_normal_shock(upstream_mach, gamma)
    downstream_shock_area_ratio = calculate_isentropic(
        shock.downstream_mach, gamma
    ).area_critical_ratio
    exit_to_downstream_critical = (
        exit_area_ratio / shock_area_ratio * downstream_shock_area_ratio
    )
    exit_mach = mach_from_isentropic(
        "area_ratio_subsonic", exit_to_downstream_critical, gamma
    )
    exit_total_pressure = chamber_pressure_pa * shock.stagnation_pressure_ratio
    exit_pressure = (
        exit_total_pressure
        * calculate_isentropic(exit_mach, gamma).pressure_total_ratio
    )
    return (
        exit_pressure,
        upstream_mach,
        shock.downstream_mach,
        shock.stagnation_pressure_ratio,
        exit_mach,
    )


def _solve_internal_shock(
    *,
    ambient_pressure_pa: float,
    chamber_pressure_pa: float,
    area_ratio: float,
    gamma: float,
) -> ShockLocation:
    """Locate the normal shock whose subsonic exit pressure equals the back pressure.

    Notes
    -----
    The exit pressure of :func:`_shock_exit_state` decreases monotonically as
    the shock moves downstream, so the shock area ratio is found by bracketed
    bisection between just downstream of the throat and the exit plane
    (tolerance 1e-10, at most 200 iterations) [1]_, [2]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    epsilon = max(1.0e-9, (area_ratio - 1.0) * 1.0e-9)
    result = solve_bisect_monotonic(
        lambda shock_area: _shock_exit_state(
            shock_area,
            exit_area_ratio=area_ratio,
            chamber_pressure_pa=chamber_pressure_pa,
            gamma=gamma,
        )[0],
        ambient_pressure_pa,
        1.0 + epsilon,
        area_ratio,
        tolerance=1.0e-10,
        max_iterations=200,
    )
    _, upstream, downstream, total_ratio, exit_mach = _shock_exit_state(
        result.value,
        exit_area_ratio=area_ratio,
        chamber_pressure_pa=chamber_pressure_pa,
        gamma=gamma,
    )
    return ShockLocation(
        area_ratio=result.value,
        x_fraction=(result.value - 1.0) / max(area_ratio - 1.0, 1.0e-12),
        upstream_mach=upstream,
        downstream_mach=downstream,
        stagnation_pressure_ratio=total_ratio,
        downstream_exit_mach=exit_mach,
    )


#: Upper end of the separation data behind Schmucker's criterion (NASA TM-77396, p. 39).
SCHMUCKER_MAXIMUM_MACH = 5.0
#: Mach grid spacing of the sign-change search before bisection (library choice).
_SCHMUCKER_SEARCH_STEP = 0.01


def schmucker_separation_point(
    *, chamber_pressure_pa: float, ambient_pressure_pa: float, exit_mach: float, gamma: float
) -> tuple[float | None, bool]:
    """Return the wall Mach number of incipient separation by Schmucker's criterion.

    Parameters
    ----------
    chamber_pressure_pa : float
        Stagnation pressure of the attached supersonic flow.
    ambient_pressure_pa : float
        Back pressure ``p_a``; zero gives no separation.
    exit_mach : float
        One-dimensional exit Mach number of the full-flowing nozzle (>= 1).
    gamma : float
        Ratio of specific heats.

    Returns
    -------
    tuple of (float or None, bool)
        The Mach number ``M_i`` in ``[1, min(exit_mach, 5)]`` where the
        isentropic wall pressure first falls to ``p_a (1.88 M - 1)^(-0.64)``
        (``None`` if it never does), and whether the search covered the
        whole divergent section (False when it stopped at ``M = 5``).

    Notes
    -----
    The residual ``ln(p/p_a) + 0.64 ln(1.88 M - 1)`` is not monotone just
    downstream of the throat, so its first sign change is found on a
    ``0.01`` Mach grid and refined by bisection to a relative 1e-12. The
    wall pressure is the one-dimensional isentropic value; a contoured
    bell's wall Mach number differs from it, which is not modelled here [1]_.

    References
    ----------
    .. [1] NASA-TM-77396-SCHMUCKER. R. H. Schmucker, NASA TM-77396, 1984,
           Eq. (11), printed pp. 38-39.
           https://ntrs.nasa.gov/citations/19840011402
    """

    if not (isfinite(chamber_pressure_pa) and chamber_pressure_pa > 0.0):
        raise DomainError("Chamber pressure must be positive and finite.", field="chamber_pressure_pa")
    if not (isfinite(exit_mach) and exit_mach >= 1.0):
        raise DomainError("Exit Mach number must be finite and at least 1.", field="exit_mach")
    if not (isfinite(gamma) and gamma > 1.0):
        raise DomainError("gamma must be finite and greater than 1.", field="gamma")
    if not isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0.0:
        raise DomainError("Ambient pressure must be finite and non-negative.", field="ambient_pressure_pa")
    if ambient_pressure_pa == 0.0:
        return None, True
    upper = min(exit_mach, SCHMUCKER_MAXIMUM_MACH)
    log_ratio = log(chamber_pressure_pa / ambient_pressure_pa)
    exponent = gamma / (gamma - 1.0)

    def residual(mach: float) -> float:
        # ln(p/p_a) - ln((1.88 M - 1)^(-0.64)); positive while the flow is attached
        wall = log_ratio - exponent * log1p(0.5 * (gamma - 1.0) * mach * mach)
        return wall + 0.64 * log(1.88 * mach - 1.0)

    lower = 1.0
    if residual(lower) <= 0.0:
        return 1.0, True
    while lower < upper:
        step = min(upper, lower + _SCHMUCKER_SEARCH_STEP)
        if residual(step) <= 0.0:
            low, high = lower, step
            while high - low > 1.0e-12 * high:
                middle = 0.5 * (low + high)
                if residual(middle) > 0.0:
                    low = middle
                else:
                    high = middle
            return high, True
        lower = step
    return None, exit_mach <= SCHMUCKER_MAXIMUM_MACH


def _separation_assessment(
    *, exit_pressure_pa: float, ambient_pressure_pa: float, chamber_pressure_pa: float,
    exit_mach: float | None, gamma: float,
) -> SeparationAssessment:
    """Assess flow separation with the Summerfield and Schmucker criteria.

    Notes
    -----
    Summerfield flags ``p_exit/p_ambient <= 0.4`` [1]_, an onset screen for
    roughly 15 degree conical divergent sections at chamber-to-ambient
    pressure ratios above 16. Schmucker's Mach-number criterion [2]_ is
    evaluated along the attached supersonic branch when the throat is
    choked (``exit_mach`` given) and locates the incipient separation
    point. Neither is a viscous solution or a side-load estimate. Vacuum
    back pressure yields no ratio and no flag.

    References
    ----------
    .. [1] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976, p. 17.
           https://ntrs.nasa.gov/citations/19770009165
    .. [2] NASA-TM-77396-SCHMUCKER. R. H. Schmucker, NASA TM-77396, 1984,
           Eq. (11), printed pp. 38-39. https://ntrs.nasa.gov/citations/19840011402
    """

    ratio = (
        exit_pressure_pa / ambient_pressure_pa if ambient_pressure_pa > 0.0 else None
    )
    summerfield = ratio is not None and ratio <= 0.4
    mach_i: float | None = None
    within = True
    if exit_mach is not None and exit_mach > 1.0:
        mach_i, within = schmucker_separation_point(
            chamber_pressure_pa=chamber_pressure_pa, ambient_pressure_pa=ambient_pressure_pa,
            exit_mach=exit_mach, gamma=gamma)
    area_i: float | None = None
    pressure_i: float | None = None
    if mach_i is not None:
        point = calculate_isentropic(mach_i, gamma)
        area_i = point.area_critical_ratio
        pressure_i = chamber_pressure_pa * point.pressure_total_ratio
    return SeparationAssessment(
        predicted=summerfield or mach_i is not None,
        correlation="Summerfield onset screening",
        criterion="p_exit / p_ambient <= 0.4",
        exit_to_ambient_pressure_ratio=ratio,
        valid_for=(
            "approximately 15 degree conical divergent nozzles and "
            "chamber-to-ambient pressure ratio above 16"
        ),
        source="NASA SP-8120, Liquid Rocket Engine Nozzles (1976), p. 17",
        source_url="https://ntrs.nasa.gov/citations/19770009165",
        schmucker_predicted=mach_i is not None,
        schmucker_separation_mach=mach_i,
        schmucker_separation_area_ratio=area_i,
        schmucker_separation_pressure_pa=pressure_i,
        schmucker_within_range=within,
    )


def _profile(
    *,
    area_ratio: float,
    chamber_pressure_pa: float,
    gamma: float,
    regime: str,
    throat_mach: float,
    shock: ShockLocation | None,
    station_count: int,
) -> tuple[NozzleFlowStation, ...]:
    """Return evenly spaced quasi-one-dimensional stations from throat to exit.

    Notes
    -----
    The local area ratio is linear in the normalized axial coordinate. Mach
    number at each station follows from the isentropic area relation [1]_ on
    the branch implied by the regime: subsonic for unchoked and just-choked
    flow (referenced to the effective sonic area), supersonic for shock-free
    choked flow, and supersonic/subsonic on either side of a normal shock with
    the post-shock stagnation pressure behind it.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    stations: list[NozzleFlowStation] = []
    if regime == "unchoked":
        throat_area_critical = calculate_isentropic(
            throat_mach, gamma
        ).area_critical_ratio
    else:
        throat_area_critical = 1.0
    downstream_critical_ratio = None
    if shock is not None:
        downstream_critical_ratio = (
            shock.area_ratio
            / calculate_isentropic(shock.downstream_mach, gamma).area_critical_ratio
        )
    for index in range(station_count):
        x_fraction = index / (station_count - 1)
        local_area_ratio = 1.0 + x_fraction * (area_ratio - 1.0)
        if shock is not None and local_area_ratio > shock.area_ratio:
            assert downstream_critical_ratio is not None
            local_effective_ratio = local_area_ratio / downstream_critical_ratio
            mach = mach_from_isentropic(
                "area_ratio_subsonic", local_effective_ratio, gamma
            )
            total_pressure = chamber_pressure_pa * shock.stagnation_pressure_ratio
            region = "post-shock subsonic"
        elif shock is not None:
            mach = mach_from_isentropic(
                "area_ratio_supersonic", local_area_ratio, gamma
            )
            total_pressure = chamber_pressure_pa
            region = "pre-shock supersonic"
        elif regime in {"unchoked", "just-choked"}:
            effective_ratio = local_area_ratio * throat_area_critical
            mach = mach_from_isentropic(
                "area_ratio_subsonic", effective_ratio, gamma
            )
            total_pressure = chamber_pressure_pa
            region = "subsonic"
        else:
            mach = mach_from_isentropic(
                "area_ratio_supersonic", local_area_ratio, gamma
            )
            total_pressure = chamber_pressure_pa
            region = "supersonic"
        pressure = (
            total_pressure * calculate_isentropic(mach, gamma).pressure_total_ratio
        )
        stations.append(
            NozzleFlowStation(
                x_fraction=x_fraction,
                area_ratio=local_area_ratio,
                mach=mach,
                pressure_pa=pressure,
                total_pressure_pa=total_pressure,
                region=region,
            )
        )
    return tuple(stations)


def calculate_off_design_nozzle(
    *,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    throat_area_m2: float,
    area_ratio: float,
    ambient_pressure_pa: float,
    gamma: float = 1.22,
    gas_constant_j_kg_k: float = 355.0,
    station_count: int = 81,
) -> OffDesignNozzleResult:
    """Classify and solve a perfect-gas C-D nozzle at imposed back pressure.

    Parameters
    ----------
    chamber_pressure_pa : float
        Chamber stagnation pressure ``p_c`` in Pa, finite and positive.
    chamber_temperature_k : float
        Chamber stagnation temperature in K, finite and positive.
    throat_area_m2 : float
        Throat area in m^2, finite and positive.
    area_ratio : float
        Exit-to-throat area ratio, at least one.
    ambient_pressure_pa : float
        Imposed back pressure ``p_a`` in Pa, finite and non-negative.
    gamma : float, default 1.22
        Frozen ratio of specific heats, greater than one.
    gas_constant_j_kg_k : float, default 355.0
        Specific gas constant in J/(kg K), finite and positive.
    station_count : int, default 81
        Number of profile stations, between 21 and 401.

    Returns
    -------
    OffDesignNozzleResult
        Regime, exit state, thrust, separation screening, warnings and the
        station profile.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    With thresholds from :func:`nozzle_regime_thresholds` and a pressure
    tolerance of ``1e-7 p_c``, the regimes are, in order of decreasing back
    pressure [1]_, [2]_:

    * ``no-flow`` for ``p_a >= p_c``;
    * ``unchoked``: subsonic isentropic flow with ``p_exit = p_a``, mass flow
      evaluated at the exit;
    * ``just-choked``: sonic throat, subsonic diverging section;
    * ``internal-normal-shock``: a normal shock located by bisection so that
      the subsonic exit pressure equals ``p_a``;
    * ``shock-at-exit``: a normal shock standing at the exit plane;
    * ``overexpanded``, ``design`` and ``underexpanded``: shock-free
      supersonic exit flow from
      :func:`~rocket_propulsion.compressible.nozzle.calculate_ideal_nozzle`,
      with external waves not resolved.

    Thrust is ``F = mdot V_e + (p_e - p_a) A_e``. The exit mass flow is
    recomputed independently and reported as a relative continuity residual.
    Separation is screened, not solved [3]_. Viscous, divergence, chemistry
    and discharge losses are not included.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    .. [3] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976, p. 17.
           https://ntrs.nasa.gov/citations/19710021390
    """

    _validate_inputs(
        chamber_pressure_pa,
        chamber_temperature_k,
        throat_area_m2,
        area_ratio,
        ambient_pressure_pa,
        gamma,
        gas_constant_j_kg_k,
    )
    if not 21 <= station_count <= 401:
        raise DomainError("Station count must be between 21 and 401.")
    thresholds = nozzle_regime_thresholds(
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratio=area_ratio,
        gamma=gamma,
    )
    tolerance = 1.0e-7 * chamber_pressure_pa
    shock: ShockLocation | None = None
    exit_area = throat_area_m2 * area_ratio

    if ambient_pressure_pa >= chamber_pressure_pa:
        regime = "no-flow"
        throat_mach = 0.0
        exit_mach = 0.0
        exit_pressure = ambient_pressure_pa
        exit_temperature = chamber_temperature_k
        exit_velocity = 0.0
        exit_total_pressure = chamber_pressure_pa
        mass_flow = 0.0
        exit_mass_flow = 0.0
    elif ambient_pressure_pa > thresholds.just_choked_back_pressure_pa + tolerance:
        regime = "unchoked"
        exit_pressure = ambient_pressure_pa
        exit_pressure_ratio = exit_pressure / chamber_pressure_pa
        exit_mach = mach_from_isentropic(
            "pressure_total_ratio", exit_pressure_ratio, gamma
        )
        exit_area_critical = calculate_isentropic(
            exit_mach, gamma
        ).area_critical_ratio
        throat_area_critical = exit_area_critical / area_ratio
        throat_mach = mach_from_isentropic(
            "area_ratio_subsonic", throat_area_critical, gamma
        )
        ratios = calculate_isentropic(exit_mach, gamma)
        exit_temperature = chamber_temperature_k * ratios.temperature_total_ratio
        exit_velocity = exit_mach * sqrt(
            gamma * gas_constant_j_kg_k * exit_temperature
        )
        exit_total_pressure = chamber_pressure_pa
        mass_flow = _mass_flow_at_station(
            total_pressure_pa=chamber_pressure_pa,
            total_temperature_k=chamber_temperature_k,
            area_m2=exit_area,
            mach=exit_mach,
            gamma=gamma,
            gas_constant=gas_constant_j_kg_k,
        )
        exit_mass_flow = mass_flow
    elif abs(ambient_pressure_pa - thresholds.just_choked_back_pressure_pa) <= tolerance:
        regime = "just-choked"
        throat_mach = 1.0
        exit_mach = mach_from_isentropic("area_ratio_subsonic", area_ratio, gamma)
        ratios = calculate_isentropic(exit_mach, gamma)
        exit_pressure = chamber_pressure_pa * ratios.pressure_total_ratio
        exit_temperature = chamber_temperature_k * ratios.temperature_total_ratio
        exit_velocity = exit_mach * sqrt(
            gamma * gas_constant_j_kg_k * exit_temperature
        )
        exit_total_pressure = chamber_pressure_pa
        mass_flow = _choked_mass_flow(
            chamber_pressure_pa,
            chamber_temperature_k,
            throat_area_m2,
            gamma,
            gas_constant_j_kg_k,
        )
        exit_mass_flow = _mass_flow_at_station(
            total_pressure_pa=exit_total_pressure,
            total_temperature_k=chamber_temperature_k,
            area_m2=exit_area,
            mach=exit_mach,
            gamma=gamma,
            gas_constant=gas_constant_j_kg_k,
        )
    elif ambient_pressure_pa > thresholds.shock_at_exit_back_pressure_pa + tolerance:
        regime = "internal-normal-shock"
        throat_mach = 1.0
        shock = _solve_internal_shock(
            ambient_pressure_pa=ambient_pressure_pa,
            chamber_pressure_pa=chamber_pressure_pa,
            area_ratio=area_ratio,
            gamma=gamma,
        )
        exit_mach = shock.downstream_exit_mach
        exit_total_pressure = chamber_pressure_pa * shock.stagnation_pressure_ratio
        ratios = calculate_isentropic(exit_mach, gamma)
        exit_pressure = exit_total_pressure * ratios.pressure_total_ratio
        exit_temperature = chamber_temperature_k * ratios.temperature_total_ratio
        exit_velocity = exit_mach * sqrt(
            gamma * gas_constant_j_kg_k * exit_temperature
        )
        mass_flow = _choked_mass_flow(
            chamber_pressure_pa,
            chamber_temperature_k,
            throat_area_m2,
            gamma,
            gas_constant_j_kg_k,
        )
        exit_mass_flow = _mass_flow_at_station(
            total_pressure_pa=exit_total_pressure,
            total_temperature_k=chamber_temperature_k,
            area_m2=exit_area,
            mach=exit_mach,
            gamma=gamma,
            gas_constant=gas_constant_j_kg_k,
        )
    elif abs(ambient_pressure_pa - thresholds.shock_at_exit_back_pressure_pa) <= tolerance:
        regime = "shock-at-exit"
        throat_mach = 1.0
        upstream_mach = mach_from_isentropic(
            "area_ratio_supersonic", area_ratio, gamma
        )
        normal_shock = calculate_normal_shock(upstream_mach, gamma)
        shock = ShockLocation(
            area_ratio=area_ratio,
            x_fraction=1.0,
            upstream_mach=upstream_mach,
            downstream_mach=normal_shock.downstream_mach,
            stagnation_pressure_ratio=normal_shock.stagnation_pressure_ratio,
            downstream_exit_mach=normal_shock.downstream_mach,
        )
        exit_mach = normal_shock.downstream_mach
        exit_total_pressure = (
            chamber_pressure_pa * normal_shock.stagnation_pressure_ratio
        )
        ratios = calculate_isentropic(exit_mach, gamma)
        exit_pressure = exit_total_pressure * ratios.pressure_total_ratio
        exit_temperature = chamber_temperature_k * ratios.temperature_total_ratio
        exit_velocity = exit_mach * sqrt(
            gamma * gas_constant_j_kg_k * exit_temperature
        )
        mass_flow = _choked_mass_flow(
            chamber_pressure_pa,
            chamber_temperature_k,
            throat_area_m2,
            gamma,
            gas_constant_j_kg_k,
        )
        exit_mass_flow = _mass_flow_at_station(
            total_pressure_pa=exit_total_pressure,
            total_temperature_k=chamber_temperature_k,
            area_m2=exit_area,
            mach=exit_mach,
            gamma=gamma,
            gas_constant=gas_constant_j_kg_k,
        )
    else:
        ideal = calculate_ideal_nozzle(
            chamber_pressure_pa=chamber_pressure_pa,
            chamber_temperature_k=chamber_temperature_k,
            throat_area_m2=throat_area_m2,
            area_ratio=area_ratio,
            ambient_pressure_pa=ambient_pressure_pa,
            gamma=gamma,
            gas_constant_j_kg_k=gas_constant_j_kg_k,
        )
        throat_mach = 1.0
        exit_mach = ideal.exit_mach
        exit_pressure = ideal.exit_pressure_pa
        exit_temperature = ideal.exit_temperature_k
        exit_velocity = ideal.exit_velocity_m_s
        exit_total_pressure = chamber_pressure_pa
        mass_flow = ideal.mass_flow_kg_s
        exit_mass_flow = _mass_flow_at_station(
            total_pressure_pa=chamber_pressure_pa,
            total_temperature_k=chamber_temperature_k,
            area_m2=exit_area,
            mach=exit_mach,
            gamma=gamma,
            gas_constant=gas_constant_j_kg_k,
        )
        if abs(ambient_pressure_pa - thresholds.design_exit_pressure_pa) <= tolerance:
            regime = "design"
        elif ambient_pressure_pa > thresholds.design_exit_pressure_pa:
            regime = "overexpanded"
        else:
            regime = "underexpanded"

    momentum_thrust = mass_flow * exit_velocity
    pressure_thrust = (exit_pressure - ambient_pressure_pa) * exit_area
    thrust = momentum_thrust + pressure_thrust
    relative_residual = (
        abs(exit_mass_flow - mass_flow) / mass_flow if mass_flow > 0.0 else 0.0
    )
    # Schmucker's criterion needs an attached supersonic wall flow, i.e. a choked throat
    # followed by the supersonic branch (up to the shock, if there is one).
    supersonic_branch = regime in {"internal-normal-shock", "shock-at-exit", "overexpanded",
                                   "design", "underexpanded"}
    separation = _separation_assessment(
        exit_pressure_pa=exit_pressure, ambient_pressure_pa=ambient_pressure_pa,
        chamber_pressure_pa=chamber_pressure_pa,
        exit_mach=(mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
                   if supersonic_branch else None),
        gamma=gamma,
    )
    warnings: list[WarningMessage] = []
    if separation.predicted:
        warnings.append(
            WarningMessage(
                code="possible_flow_separation",
                severity="critical",
                message=(
                    "Empirical separation criteria ("
                    + "; ".join(name for name, hit in (
                        ("Summerfield", separation.exit_to_ambient_pressure_ratio is not None
                         and separation.exit_to_ambient_pressure_ratio <= 0.4),
                        (f"Schmucker at A/A* = {separation.schmucker_separation_area_ratio:.4g}"
                         if separation.schmucker_separation_area_ratio is not None else "", 
                         separation.schmucker_predicted)) if hit)
                    + ") indicate possible boundary-layer separation; "
                    "this is not a viscous flow solution."
                ),
            )
        )
    if regime in {"overexpanded", "shock-at-exit"}:
        warnings.append(
            WarningMessage(
                code="external_compression_structure",
                message="External shocks are not resolved by the one-dimensional model.",
            )
        )
    profile = () if regime == "no-flow" else _profile(
        area_ratio=area_ratio,
        chamber_pressure_pa=chamber_pressure_pa,
        gamma=gamma,
        regime=regime,
        throat_mach=throat_mach,
        shock=shock,
        station_count=station_count,
    )
    return OffDesignNozzleResult(
        regime=regime,
        chamber_pressure_pa=chamber_pressure_pa,
        ambient_pressure_pa=ambient_pressure_pa,
        area_ratio=area_ratio,
        throat_area_m2=throat_area_m2,
        exit_area_m2=exit_area,
        throat_mach=throat_mach,
        exit_mach=exit_mach,
        exit_pressure_pa=exit_pressure,
        exit_temperature_k=exit_temperature,
        exit_velocity_m_s=exit_velocity,
        exit_total_pressure_pa=exit_total_pressure,
        mass_flow_kg_s=mass_flow,
        exit_mass_flow_kg_s=exit_mass_flow,
        mass_flow_relative_residual=relative_residual,
        momentum_thrust_n=momentum_thrust,
        pressure_thrust_n=pressure_thrust,
        thrust_n=thrust,
        specific_impulse_s=(
            thrust / (mass_flow * STANDARD_GRAVITY_M_S2) if mass_flow > 0.0 else None
        ),
        thresholds=thresholds,
        shock=shock,
        separation=separation,
        warnings=tuple(warnings),
        profile=profile,
    )


# U.S. Standard Atmosphere 1976 layers below 84.852 km geopotential altitude:
# (base geopotential altitude H_b [m'], base temperature T_b [K], base
# pressure P_b [Pa], molecular-scale temperature gradient L_b [K/m']).  H_b,
# T_b and L_b are the defining values of the 1976 model; P_b are the base
# pressures that its hydrostatic equations give from the sea-level pressure
# USSA1976_SEA_LEVEL_PRESSURE_PA with R* and M0, as tabulated in the 1976
# report.  Source: https://ntrs.nasa.gov/citations/19770009539
_ATMOSPHERE_LAYERS = (
    (0.0, USSA1976_SEA_LEVEL_TEMPERATURE_K, USSA1976_SEA_LEVEL_PRESSURE_PA, -0.0065),
    (11_000.0, 216.65, 22_632.06, 0.0),
    (20_000.0, 216.65, 5_474.889, 0.0010),
    (32_000.0, 228.65, 868.0187, 0.0028),
    (47_000.0, 270.65, 110.9063, 0.0),
    (51_000.0, 270.65, 66.93887, -0.0028),
    (71_000.0, 214.65, 3.956420, -0.0020),
)


def standard_atmosphere_pressure_pa(altitude_m: float) -> float:
    """Return U.S. Standard Atmosphere 1976 pressure at a geopotential altitude.

    Parameters
    ----------
    altitude_m : float
        Geopotential altitude ``H`` in m', from 0 to 84 852 m' inclusive, the
        range of the model's seven linear temperature layers.

    Returns
    -------
    float
        Static pressure in Pa.

    Raises
    ------
    DomainError
        If ``altitude_m`` is not finite or lies outside the model range.

    Notes
    -----
    Within a layer with base ``H_b``, ``T_b``, ``P_b`` and gradient ``L_b``
    [1]_::

        P = P_b (T_b / (T_b + L_b (H - H_b)))^(g0 M0 / (R* L_b))   (L_b != 0)
        P = P_b exp(-g0 M0 (H - H_b) / (R* T_b))                  (L_b == 0)

    with the model's own defining constants: ``g0`` is
    ``STANDARD_GRAVITY_M_S2``, ``M0`` is ``USSA1976_MOLAR_MASS_KG_MOL`` and
    ``R*`` is ``USSA1976_GAS_CONSTANT_J_MOL_K`` from
    :mod:`rocket_propulsion.core.constants`. Using a different ``R*`` (for
    example the 2019 SI value) makes the pressure jump by about 3e-5 at every
    layer base.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    """

    if not isfinite(altitude_m) or not 0.0 <= altitude_m <= 84_852.0:
        raise DomainError("Altitude must be between 0 and 84,852 m.", field="altitude_m")
    gravity = STANDARD_GRAVITY_M_S2
    molar_mass = USSA1976_MOLAR_MASS_KG_MOL
    universal_gas_constant = USSA1976_GAS_CONSTANT_J_MOL_K
    layer = _ATMOSPHERE_LAYERS[0]
    for candidate in _ATMOSPHERE_LAYERS:
        if candidate[0] <= altitude_m:
            layer = candidate
        else:
            break
    base_altitude, base_temperature, base_pressure, lapse = layer
    delta = altitude_m - base_altitude
    if lapse == 0.0:
        return base_pressure * exp(
            -gravity * molar_mass * delta
            / (universal_gas_constant * base_temperature)
        )
    local_temperature = base_temperature + lapse * delta
    return base_pressure * (
        base_temperature / local_temperature
    ) ** (gravity * molar_mass / (universal_gas_constant * lapse))


def generate_nozzle_altitude_sweep(
    *,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    gas_constant_j_kg_k: float = 355.0,
    maximum_altitude_m: float = 50_000.0,
    point_count: int = 51,
) -> NozzleAltitudeSweep:
    """Evaluate regime, thrust, and Isp along a standard-atmosphere climb.

    Parameters
    ----------
    chamber_pressure_pa, chamber_temperature_k, throat_area_m2, area_ratio : float
        Nozzle and chamber definition, as in
        :func:`calculate_off_design_nozzle`.
    gamma : float, default 1.22
        Frozen ratio of specific heats.
    gas_constant_j_kg_k : float, default 355.0
        Specific gas constant in J/(kg K).
    maximum_altitude_m : float, default 50000.0
        Top of the sweep in geopotential metres, in ``(0, 84852]``.
    point_count : int, default 51
        Number of evenly spaced altitudes from sea level, between 2 and 201.

    Returns
    -------
    NozzleAltitudeSweep
        Per-altitude regime, thrust, Isp and separation flags.

    Raises
    ------
    DomainError
        If the sweep bounds or nozzle inputs are outside their domains.

    Notes
    -----
    Ambient pressure comes from :func:`standard_atmosphere_pressure_pa` [1]_;
    each point is solved with :func:`calculate_off_design_nozzle` using 41
    profile stations. Separation flags are a screening result only [2]_.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    .. [2] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976, p. 17.
           https://ntrs.nasa.gov/citations/19710021390
    """

    if not 2 <= point_count <= 201:
        raise DomainError("Altitude point count must be between 2 and 201.")
    if not 0.0 < maximum_altitude_m <= 84_852.0:
        raise DomainError("Maximum altitude must be in (0, 84,852] m.")
    points: list[NozzleAltitudePoint] = []
    for index in range(point_count):
        altitude = maximum_altitude_m * index / (point_count - 1)
        ambient = standard_atmosphere_pressure_pa(altitude)
        result = calculate_off_design_nozzle(
            chamber_pressure_pa=chamber_pressure_pa,
            chamber_temperature_k=chamber_temperature_k,
            throat_area_m2=throat_area_m2,
            area_ratio=area_ratio,
            ambient_pressure_pa=ambient,
            gamma=gamma,
            gas_constant_j_kg_k=gas_constant_j_kg_k,
            station_count=41,
        )
        points.append(
            NozzleAltitudePoint(
                altitude_m=altitude,
                ambient_pressure_pa=ambient,
                regime=result.regime,
                thrust_n=result.thrust_n,
                specific_impulse_s=result.specific_impulse_s,
                chamber_to_ambient_pressure_ratio=(
                    chamber_pressure_pa / ambient if ambient > 0.0 else None
                ),
                separation_predicted=result.separation.predicted,
                shock_x_fraction=(
                    result.shock.x_fraction if result.shock is not None else None
                ),
            )
        )
    return NozzleAltitudeSweep(
        points=tuple(points),
        atmosphere_model="U.S. Standard Atmosphere 1976, geopotential layers",
        warning=(
            "Separation flags use the Summerfield p_exit/p_ambient <= 0.4 screen; "
            "no viscous or side-load solution is performed."
        ),
    )

