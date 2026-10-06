"""Back-pressure-driven converging-diverging nozzle operating regimes."""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .isentropic import calculate_isentropic, mach_from_isentropic
from .normal_shock import calculate_normal_shock
from .nozzle import STANDARD_GRAVITY_M_S2, calculate_ideal_nozzle


@dataclass(frozen=True, slots=True)
class NozzleRegimeThresholds:
    """Back-pressure boundaries for a fixed perfect-gas nozzle."""

    chamber_pressure_pa: float
    just_choked_back_pressure_pa: float
    shock_at_exit_back_pressure_pa: float
    design_exit_pressure_pa: float


@dataclass(frozen=True, slots=True)
class NozzleFlowStation:
    """One normalized station used by the off-design contour plot."""

    x_fraction: float
    area_ratio: float
    mach: float
    pressure_pa: float
    total_pressure_pa: float
    region: str


@dataclass(frozen=True, slots=True)
class ShockLocation:
    """Solved normal-shock state within the divergent section."""

    area_ratio: float
    x_fraction: float
    upstream_mach: float
    downstream_mach: float
    stagnation_pressure_ratio: float
    downstream_exit_mach: float


@dataclass(frozen=True, slots=True)
class SeparationAssessment:
    """Empirical separation screening with an explicit validity statement."""

    predicted: bool
    correlation: str
    criterion: str
    exit_to_ambient_pressure_ratio: float | None
    valid_for: str
    source: str
    source_url: str


@dataclass(frozen=True, slots=True)
class OffDesignNozzleResult:
    """Complete one-dimensional nozzle state at one ambient pressure."""

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
    """One atmosphere/nozzle result in an altitude sweep."""

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
    """Off-design performance across the 1976 standard atmosphere."""

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
    """Return analytical back-pressure boundaries for a fixed nozzle."""

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


def _separation_assessment(
    *, exit_pressure_pa: float, ambient_pressure_pa: float
) -> SeparationAssessment:
    ratio = (
        exit_pressure_pa / ambient_pressure_pa if ambient_pressure_pa > 0.0 else None
    )
    predicted = ratio is not None and ratio <= 0.4
    return SeparationAssessment(
        predicted=predicted,
        correlation="Summerfield onset screening",
        criterion="p_exit / p_ambient <= 0.4",
        exit_to_ambient_pressure_ratio=ratio,
        valid_for=(
            "approximately 15 degree conical divergent nozzles and "
            "chamber-to-ambient pressure ratio above 16"
        ),
        source="NASA SP-8041, Solid Rocket Motor Nozzles",
        source_url="https://ntrs.nasa.gov/citations/19710021390",
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
    """Classify and solve a perfect-gas C-D nozzle at imposed back pressure."""

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
    separation = _separation_assessment(
        exit_pressure_pa=exit_pressure, ambient_pressure_pa=ambient_pressure_pa
    )
    warnings: list[WarningMessage] = []
    if separation.predicted:
        warnings.append(
            WarningMessage(
                code="possible_flow_separation",
                severity="critical",
                message=(
                    "Summerfield screening indicates possible boundary-layer separation; "
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


_ATMOSPHERE_LAYERS = (
    (0.0, 288.15, 101_325.0, -0.0065),
    (11_000.0, 216.65, 22_632.06, 0.0),
    (20_000.0, 216.65, 5_474.889, 0.0010),
    (32_000.0, 228.65, 868.0187, 0.0028),
    (47_000.0, 270.65, 110.9063, 0.0),
    (51_000.0, 270.65, 66.93887, -0.0028),
    (71_000.0, 214.65, 3.956420, -0.0020),
)


def standard_atmosphere_pressure_pa(altitude_m: float) -> float:
    """Return 1976 standard-atmosphere pressure through 84.852 km."""

    if not isfinite(altitude_m) or not 0.0 <= altitude_m <= 84_852.0:
        raise DomainError("Altitude must be between 0 and 84,852 m.", field="altitude_m")
    gravity = 9.80665
    molar_mass = 0.0289644
    universal_gas_constant = 8.3144598
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
    """Evaluate regime, thrust, and Isp along a standard-atmosphere climb."""

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

