"""Reference-backed pressure-fed blowdown propulsion model.

This preliminary L2 model resolves rigid-tank ullage expansion, a declared
polytropic gas law, pressure-correlated thrust/Isp, and a converged tabulated
force/mass-flow history.  It is propulsion-only and contains no orbit or frame
logic.

References
----------
NASA-SP-8112, *Pressurization Systems for Liquid Rockets*:
https://ntrs.nasa.gov/citations/19760015212
NASA-CR-131400, *Reliability Model of a Monopropellant Auxiliary Propulsion
System*: https://ntrs.nasa.gov/citations/19730012094
NASA-CR-122347, *Monopropellant Hydrazine Resisto Jet*:
https://ntrs.nasa.gov/citations/19720011120
"""

from __future__ import annotations

import hashlib
from bisect import bisect_right
from dataclasses import dataclass
from itertools import pairwise
from math import ceil, isfinite
from typing import Protocol

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .models import EngineOperatingPoint, PropellantFlow
from .serialization import canonical_json_bytes
from .tabulated import TabulatedBurnArtifact, TabulatedBurnSegment

BLOWDOWN_REFERENCE_IDS = (
    "NASA-SP-8112",
    "NASA-CR-131400",
    "NASA-CR-122347",
    "NASA-GRC-THRUST-EQUATION",
)


class PressurePerformanceModel(Protocol):
    """Structural contract for a pressure-to-operating-point model.

    Implementations must declare their pressure validity interval and refuse
    extrapolation. This keeps the blowdown integrator independent of whether
    performance came from a compact correlation or a measured hot-fire map.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    minimum_supply_pressure_pa: float

    @property
    def maximum_pressure_pa(self) -> float:
        """Return the inclusive upper validity limit."""

        ...

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Return a stream-closed point at the requested supply pressure."""

        ...


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _scaled_operating_point(
    rated_point: EngineOperatingPoint,
    *,
    supply_pressure_pa: float,
    delivered_thrust_n: float,
    system_specific_impulse_s: float,
    chamber_pressure_pa: float,
    model_id: str,
    source: str,
    warning: str,
) -> EngineOperatingPoint:
    total_flow = delivered_thrust_n / (
        STANDARD_GRAVITY_M_S2 * system_specific_impulse_s
    )
    flow_ratio = total_flow / rated_point.total_tank_flow_kg_s
    streams = tuple(
        PropellantFlow(
            stream_id=stream.stream_id,
            tank_id=stream.tank_id,
            propellant_key=stream.propellant_key,
            role=stream.role,
            destination=stream.destination,
            mass_flow_kg_s=stream.mass_flow_kg_s * flow_ratio,
            tank_depleting=stream.tank_depleting,
            contributes_to_delivered_thrust=stream.contributes_to_delivered_thrust,
        )
        for stream in rated_point.streams
    )
    chamber_flow = rated_point.chamber_mass_flow_kg_s * flow_ratio
    return EngineOperatingPoint(
        name=f"{rated_point.name} @ {supply_pressure_pa:.9g} Pa supply",
        commanded_throttle=rated_point.commanded_throttle,
        realized_throttle=rated_point.realized_throttle,
        active_engine_count=rated_point.active_engine_count,
        ideal_thrust_n=rated_point.ideal_thrust_n
        * delivered_thrust_n
        / rated_point.delivered_thrust_n,
        delivered_thrust_n=delivered_thrust_n,
        chamber_specific_impulse_s=(
            delivered_thrust_n / (STANDARD_GRAVITY_M_S2 * chamber_flow)
        ),
        system_specific_impulse_s=system_specific_impulse_s,
        chamber_mass_flow_kg_s=chamber_flow,
        total_tank_flow_kg_s=total_flow,
        oxidizer_mass_flow_kg_s=rated_point.oxidizer_mass_flow_kg_s * flow_ratio,
        fuel_mass_flow_kg_s=rated_point.fuel_mass_flow_kg_s * flow_ratio,
        mixture_ratio=rated_point.mixture_ratio,
        chamber_pressure_pa=chamber_pressure_pa,
        ambient_pressure_pa=rated_point.ambient_pressure_pa,
        streams=streams,
        model_id=model_id,
        source=source,
        warnings=rated_point.warnings + (warning,),
    )


@dataclass(frozen=True, slots=True)
class BlowdownTankModel:
    """Rigid tank with incompressible propellant and polytropic ullage gas.

    The pressure law is ``p = p0 * (V0 / V)**n``. ``n=1`` represents the
    isothermal ideal-gas limit; a larger declared exponent approximates reduced
    heat transfer. Propellant density is constant and dissolved-gas, diaphragm,
    vapor-pressure, line-loss, and thermal-soak effects are outside this model.

    References
    ----------
    NASA-SP-8112, blowdown-system section:
    https://ntrs.nasa.gov/citations/19760015212
    NASA-CR-131400, unregulated monopropellant performance history:
    https://ntrs.nasa.gov/citations/19730012094
    """

    initial_pressure_pa: float
    initial_ullage_volume_m3: float
    initial_propellant_mass_kg: float
    propellant_density_kg_m3: float
    polytropic_exponent: float = 1.0

    def __post_init__(self) -> None:
        for value, label in (
            (self.initial_pressure_pa, "Initial ullage pressure"),
            (self.initial_ullage_volume_m3, "Initial ullage volume"),
            (self.initial_propellant_mass_kg, "Initial propellant mass"),
            (self.propellant_density_kg_m3, "Propellant density"),
            (self.polytropic_exponent, "Polytropic exponent"),
        ):
            _positive(value, label)

    @property
    def initial_propellant_volume_m3(self) -> float:
        """Return initial liquid volume from the incompressible density model.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        """

        return self.initial_propellant_mass_kg / self.propellant_density_kg_m3

    @property
    def tank_volume_m3(self) -> float:
        """Return rigid internal tank volume represented by the model.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        """

        return self.initial_ullage_volume_m3 + self.initial_propellant_volume_m3

    def ullage_volume_m3(self, remaining_propellant_mass_kg: float) -> float:
        """Return ullage volume after expelling a declared propellant mass.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        """

        self._validate_remaining_mass(remaining_propellant_mass_kg)
        expelled_mass = self.initial_propellant_mass_kg - remaining_propellant_mass_kg
        return self.initial_ullage_volume_m3 + expelled_mass / self.propellant_density_kg_m3

    def pressure_pa(self, remaining_propellant_mass_kg: float) -> float:
        """Evaluate the declared polytropic ullage-pressure law.

        Notes
        -----
        This is a preliminary thermodynamic closure, not a heat-transfer or
        pressurant-dissolution model. The exponent is an explicit input rather
        than an inferred gas property.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
        """

        ullage = self.ullage_volume_m3(remaining_propellant_mass_kg)
        return self.initial_pressure_pa * (
            self.initial_ullage_volume_m3 / ullage
        ) ** self.polytropic_exponent

    def remaining_mass_at_pressure(self, pressure_pa: float) -> float:
        """Invert the polytropic law for a pressure-triggered cutoff.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        """

        _positive(pressure_pa, "Requested ullage pressure")
        if pressure_pa >= self.initial_pressure_pa:
            return self.initial_propellant_mass_kg
        required_ullage = self.initial_ullage_volume_m3 * (
            self.initial_pressure_pa / pressure_pa
        ) ** (1.0 / self.polytropic_exponent)
        expelled_mass = (
            required_ullage - self.initial_ullage_volume_m3
        ) * self.propellant_density_kg_m3
        return min(
            self.initial_propellant_mass_kg,
            max(0.0, self.initial_propellant_mass_kg - expelled_mass),
        )

    def _validate_remaining_mass(self, remaining_propellant_mass_kg: float) -> None:
        if (
            not isfinite(remaining_propellant_mass_kg)
            or remaining_propellant_mass_kg < 0.0
            or remaining_propellant_mass_kg > self.initial_propellant_mass_kg
        ):
            raise DomainError(
                "Remaining propellant mass must be within the tank's initial inventory."
            )


@dataclass(frozen=True, slots=True)
class PressurePerformanceLaw:
    """Declared pressure correlation around a measured/reference operating point.

    ``F/Fref = (p/pref)**a`` and ``Isp/Isp_ref = (p/pref)**b``. The exponents
    are calibration inputs, not universal rocket laws. NASA-CR-122347 reports
    nearly linear thrust decline with supply pressure for its tested device and
    only small Isp variation over a 200-to-100 psi blowdown; other hardware
    requires its own data fit.

    References
    ----------
    NASA-CR-122347, sections 3.2 and 3.3:
    https://ntrs.nasa.gov/citations/19720011120
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    NASA rocket thrust equation:
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    """

    reference_supply_pressure_pa: float
    minimum_supply_pressure_pa: float
    thrust_pressure_exponent: float = 1.0
    specific_impulse_pressure_exponent: float = 0.0
    chamber_pressure_exponent: float = 1.0
    maximum_supply_pressure_pa: float | None = None

    def __post_init__(self) -> None:
        _positive(self.reference_supply_pressure_pa, "Reference supply pressure")
        _positive(self.minimum_supply_pressure_pa, "Minimum supply pressure")
        maximum = self.maximum_supply_pressure_pa or self.reference_supply_pressure_pa
        _positive(maximum, "Maximum supply pressure")
        if not self.minimum_supply_pressure_pa < maximum:
            raise DomainError("Minimum supply pressure must be below maximum pressure.")
        if not self.minimum_supply_pressure_pa <= self.reference_supply_pressure_pa <= maximum:
            raise DomainError("Reference supply pressure must lie inside the validity range.")
        for value, label in (
            (self.thrust_pressure_exponent, "Thrust pressure exponent"),
            (self.chamber_pressure_exponent, "Chamber-pressure exponent"),
        ):
            _nonnegative(value, label)
        if not isfinite(self.specific_impulse_pressure_exponent):
            raise DomainError("Specific-impulse pressure exponent must be finite.")

    @property
    def maximum_pressure_pa(self) -> float:
        """Return the declared upper validity limit.

        References
        ----------
        NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
        """

        return self.maximum_supply_pressure_pa or self.reference_supply_pressure_pa

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Scale a stream-closed point while preserving system-Isp identity.

        All tank-depleting and chamber streams retain their rated proportions.
        Total flow is recomputed from ``F/(g0*Isp)``; it is never inferred from
        chamber flow alone.

        References
        ----------
        NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
        NASA rocket thrust equation:
        https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
        """

        tolerance = max(1e-9, self.maximum_pressure_pa * 1e-12)
        if not isfinite(supply_pressure_pa) or (
            supply_pressure_pa < self.minimum_supply_pressure_pa - tolerance
            or supply_pressure_pa > self.maximum_pressure_pa + tolerance
        ):
            raise DomainError("Supply pressure lies outside the declared correlation range.")
        supply_pressure_pa = min(
            self.maximum_pressure_pa,
            max(self.minimum_supply_pressure_pa, supply_pressure_pa),
        )
        pressure_ratio = supply_pressure_pa / self.reference_supply_pressure_pa
        thrust_ratio = pressure_ratio**self.thrust_pressure_exponent
        isp_ratio = pressure_ratio**self.specific_impulse_pressure_exponent
        return _scaled_operating_point(
            rated_point,
            supply_pressure_pa=supply_pressure_pa,
            delivered_thrust_n=rated_point.delivered_thrust_n * thrust_ratio,
            system_specific_impulse_s=rated_point.system_specific_impulse_s * isp_ratio,
            chamber_pressure_pa=(
                rated_point.chamber_pressure_pa
                * pressure_ratio**self.chamber_pressure_exponent
            ),
            model_id="pressure_correlated_blowdown_v1",
            source=f"{rated_point.source}; pressure law calibrated by user",
            warning=(
                "Pressure exponents are hardware-specific calibration inputs; "
                "do not extrapolate beyond the declared supply-pressure range."
            ),
        )


@dataclass(frozen=True, slots=True)
class PressurePerformanceSample:
    """One measured or validated pressure-performance map row.

    Values are absolute SI quantities rather than normalized multipliers, so
    the map remains auditable without access to the original rated point.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    supply_pressure_pa: float
    delivered_thrust_n: float
    system_specific_impulse_s: float
    chamber_pressure_pa: float

    def __post_init__(self) -> None:
        for value, label in (
            (self.supply_pressure_pa, "Map supply pressure"),
            (self.delivered_thrust_n, "Map delivered thrust"),
            (self.system_specific_impulse_s, "Map system specific impulse"),
            (self.chamber_pressure_pa, "Map chamber pressure"),
        ):
            _positive(value, label)


@dataclass(frozen=True, slots=True)
class PressurePerformanceTable:
    """Piecewise-linear, non-extrapolating pressure-performance map.

    The caller supplies pressure-ordered test or validated-model samples.
    Interpolation is linear in each measured quantity. The table deliberately
    refuses extrapolation and labels the resulting operating point with the
    declared source identifier.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    samples: tuple[PressurePerformanceSample, ...]
    source_id: str

    def __post_init__(self) -> None:
        if len(self.samples) < 2:
            raise DomainError("A pressure-performance table requires at least two rows.")
        if not self.source_id.strip():
            raise DomainError("Pressure-performance table source must not be empty.")
        pressures = tuple(sample.supply_pressure_pa for sample in self.samples)
        if any(end <= start for start, end in pairwise(pressures)):
            raise DomainError("Pressure-performance samples must be strictly increasing.")

    @property
    def minimum_supply_pressure_pa(self) -> float:
        """Return the inclusive lower validity limit.

        References
        ----------
        NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
        """

        return self.samples[0].supply_pressure_pa

    @property
    def maximum_pressure_pa(self) -> float:
        """Return the inclusive upper validity limit.

        References
        ----------
        NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
        """

        return self.samples[-1].supply_pressure_pa

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Interpolate the measured map and preserve tank-stream closure.

        References
        ----------
        NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
        NASA rocket thrust equation:
        https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
        """

        tolerance = max(1e-9, self.maximum_pressure_pa * 1e-12)
        if not isfinite(supply_pressure_pa) or (
            supply_pressure_pa < self.minimum_supply_pressure_pa - tolerance
            or supply_pressure_pa > self.maximum_pressure_pa + tolerance
        ):
            raise DomainError("Supply pressure lies outside the tabulated map range.")
        supply_pressure_pa = min(
            self.maximum_pressure_pa,
            max(self.minimum_supply_pressure_pa, supply_pressure_pa),
        )
        pressures = tuple(sample.supply_pressure_pa for sample in self.samples)
        upper_index = min(bisect_right(pressures, supply_pressure_pa), len(self.samples) - 1)
        lower_index = max(0, upper_index - 1)
        lower = self.samples[lower_index]
        upper = self.samples[upper_index]
        fraction = (
            0.0
            if lower is upper
            else (supply_pressure_pa - lower.supply_pressure_pa)
            / (upper.supply_pressure_pa - lower.supply_pressure_pa)
        )

        def interpolate(start: float, end: float) -> float:
            return start + fraction * (end - start)

        return _scaled_operating_point(
            rated_point,
            supply_pressure_pa=supply_pressure_pa,
            delivered_thrust_n=interpolate(
                lower.delivered_thrust_n, upper.delivered_thrust_n
            ),
            system_specific_impulse_s=interpolate(
                lower.system_specific_impulse_s, upper.system_specific_impulse_s
            ),
            chamber_pressure_pa=interpolate(
                lower.chamber_pressure_pa, upper.chamber_pressure_pa
            ),
            model_id="tabulated_pressure_performance_v1",
            source=f"{rated_point.source}; pressure map {self.source_id}",
            warning=(
                "Pressure performance was linearly interpolated inside the declared "
                "test/model map; no extrapolation was permitted."
            ),
        )


@dataclass(frozen=True, slots=True)
class BlowdownSample:
    """One state in a pressure-fed blowdown history.

    References
    ----------
    NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
    """

    time_s: float
    remaining_propellant_mass_kg: float
    ullage_volume_m3: float
    supply_pressure_pa: float
    delivered_thrust_n: float
    axial_thrust_n: float
    system_specific_impulse_s: float
    total_mass_flow_kg_s: float


@dataclass(frozen=True, slots=True)
class BlowdownIntegrationEvidence:
    """Mesh-refinement evidence for the preliminary L2 history.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    method: str
    mass_step_kg: float
    refinement_count: int
    relative_tolerance: float
    consumed_propellant_kg: float
    delivered_impulse_n_s: float
    mass_closure_error_kg: float
    impulse_closure_error_n_s: float
    duration_change_s: float
    centroid_change_s: float


@dataclass(frozen=True, slots=True)
class BlowdownSimulationResult:
    """Reference-labelled L2 blowdown profile and exchange artifact.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    """

    samples: tuple[BlowdownSample, ...]
    artifact: TabulatedBurnArtifact
    evidence: BlowdownIntegrationEvidence
    cutoff_reason: str
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]


def _tank_flows(point: EngineOperatingPoint) -> tuple[tuple[str, float], ...]:
    flows: dict[str, float] = {}
    for stream in point.streams:
        if stream.tank_depleting:
            flows[stream.tank_id] = flows.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s
    return tuple(sorted(flows.items()))


def _simpson(start: float, middle: float, end: float, width: float) -> float:
    return width * (start + 4.0 * middle + end) / 6.0


def _build_blowdown_mesh(
    tank: BlowdownTankModel,
    rated_point: EngineOperatingPoint,
    pressure_law: PressurePerformanceModel,
    *,
    final_remaining_mass_kg: float,
    cant_efficiency: float,
    maximum_mass_step_kg: float,
    source_hash: str,
) -> tuple[tuple[BlowdownSample, ...], TabulatedBurnArtifact, float]:
    consumed = tank.initial_propellant_mass_kg - final_remaining_mass_kg
    segment_count = max(1, ceil(consumed / maximum_mass_step_kg))
    mass_step = consumed / segment_count

    def state(expelled_mass_kg: float) -> tuple[float, EngineOperatingPoint]:
        remaining = tank.initial_propellant_mass_kg - expelled_mass_kg
        pressure = tank.pressure_pa(remaining)
        return remaining, pressure_law.scale_operating_point(rated_point, pressure)

    samples: list[BlowdownSample] = []
    segments: list[TabulatedBurnSegment] = []
    elapsed = 0.0
    simpson_impulse = 0.0
    for index in range(segment_count + 1):
        expelled = index * mass_step
        remaining, point = state(expelled)
        pressure = tank.pressure_pa(remaining)
        samples.append(
            BlowdownSample(
                time_s=elapsed,
                remaining_propellant_mass_kg=remaining,
                ullage_volume_m3=tank.ullage_volume_m3(remaining),
                supply_pressure_pa=pressure,
                delivered_thrust_n=point.delivered_thrust_n,
                axial_thrust_n=point.delivered_thrust_n * cant_efficiency,
                system_specific_impulse_s=point.system_specific_impulse_s,
                total_mass_flow_kg_s=point.total_tank_flow_kg_s,
            )
        )
        if index == segment_count:
            break

        _, end_point = state(expelled + mass_step)
        _, middle_point = state(expelled + 0.5 * mass_step)
        duration = _simpson(
            1.0 / point.total_tank_flow_kg_s,
            1.0 / middle_point.total_tank_flow_kg_s,
            1.0 / end_point.total_tank_flow_kg_s,
            mass_step,
        )
        impulse = _simpson(
            STANDARD_GRAVITY_M_S2 * point.system_specific_impulse_s,
            STANDARD_GRAVITY_M_S2 * middle_point.system_specific_impulse_s,
            STANDARD_GRAVITY_M_S2 * end_point.system_specific_impulse_s,
            mass_step,
        )
        simpson_impulse += impulse
        segments.append(
            TabulatedBurnSegment(
                start_time_s=elapsed,
                end_time_s=elapsed + duration,
                start_delivered_thrust_n=point.delivered_thrust_n,
                end_delivered_thrust_n=end_point.delivered_thrust_n,
                start_axial_thrust_n=point.delivered_thrust_n * cant_efficiency,
                end_axial_thrust_n=end_point.delivered_thrust_n * cant_efficiency,
                start_total_mass_flow_kg_s=point.total_tank_flow_kg_s,
                end_total_mass_flow_kg_s=end_point.total_tank_flow_kg_s,
                start_stream_mass_flows_kg_s=_tank_flows(point),
                end_stream_mass_flows_kg_s=_tank_flows(end_point),
                phase="pressure-fed-blowdown",
            )
        )
        elapsed += duration
    artifact = TabulatedBurnArtifact(
        segments=tuple(segments),
        source_id="pressure-fed-blowdown",
        source_sha256=source_hash,
        reference_ids=BLOWDOWN_REFERENCE_IDS,
        warnings=(
            (
                "Preliminary polytropic blowdown model; validate pressure exponents "
                "against hardware-specific hot-fire data."
            ),
        ),
    )
    return tuple(samples), artifact, simpson_impulse


def simulate_pressure_fed_blowdown(
    tank: BlowdownTankModel,
    rated_point: EngineOperatingPoint,
    pressure_law: PressurePerformanceModel,
    *,
    reserve_mass_kg: float = 0.0,
    cant_efficiency: float = 1.0,
    maximum_mass_step_kg: float = 0.1,
    relative_tolerance: float = 1e-7,
    maximum_refinements: int = 12,
) -> BlowdownSimulationResult:
    """Simulate an unregulated pressure-fed burn to reserve or pressure cutoff.

    Integration uses expelled propellant mass as the independent variable:
    ``dt/dm = 1/mdot`` and ``dJ/dm = g0*Isp``. Composite Simpson integration
    supplies the reference duration/impulse, while the exchange artifact is
    independently integrated and mesh-refined until mass, impulse, duration,
    and centroid residuals satisfy the requested tolerance.

    This model is preliminary. The pressure-performance exponents must be fit
    to a thruster map or hot-fire data; the function refuses extrapolation.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    _nonnegative(reserve_mass_kg, "Blowdown reserve mass")
    if reserve_mass_kg >= tank.initial_propellant_mass_kg:
        raise DomainError("Blowdown reserve must be below initial propellant mass.")
    if not isfinite(cant_efficiency) or not 0.0 < cant_efficiency <= 1.0:
        raise DomainError("Blowdown cant efficiency must be finite and in (0, 1].")
    _positive(maximum_mass_step_kg, "Maximum blowdown mass step")
    _positive(relative_tolerance, "Blowdown relative tolerance")
    if maximum_refinements < 0:
        raise DomainError("Maximum blowdown refinements must be non-negative.")
    if tank.initial_pressure_pa > pressure_law.maximum_pressure_pa * (1.0 + 1e-12):
        raise DomainError("Initial tank pressure exceeds the pressure-law validity range.")
    if tank.initial_pressure_pa < pressure_law.minimum_supply_pressure_pa:
        raise DomainError("Initial tank pressure is below the minimum operating pressure.")

    pressure_cutoff_mass = tank.remaining_mass_at_pressure(
        pressure_law.minimum_supply_pressure_pa
    )
    final_remaining = max(reserve_mass_kg, pressure_cutoff_mass)
    cutoff_reason = (
        "minimum_supply_pressure"
        if pressure_cutoff_mass >= reserve_mass_kg
        else "propellant_reserve"
    )
    consumed = tank.initial_propellant_mass_kg - final_remaining
    if consumed <= 0.0:
        raise DomainError("Declared blowdown limits leave no usable propellant.")

    source_hash = hashlib.sha256(
        canonical_json_bytes(
            {
                "tank": tank,
                "rated_point": rated_point,
                "pressure_law": pressure_law,
                "reserve_mass_kg": reserve_mass_kg,
                "cant_efficiency": cant_efficiency,
            }
        )
    ).hexdigest()
    step = min(maximum_mass_step_kg, consumed)
    previous_duration: float | None = None
    previous_centroid: float | None = None
    samples: tuple[BlowdownSample, ...] | None = None
    artifact: TabulatedBurnArtifact | None = None
    simpson_impulse = 0.0
    residuals: tuple[float, float, float, float] | None = None
    refinement = 0
    for refinement in range(maximum_refinements + 1):
        samples, artifact, simpson_impulse = _build_blowdown_mesh(
            tank,
            rated_point,
            pressure_law,
            final_remaining_mass_kg=final_remaining,
            cant_efficiency=cant_efficiency,
            maximum_mass_step_kg=step,
            source_hash=source_hash,
        )
        mass_error = abs(artifact.consumed_propellant_kg - consumed)
        impulse_error = abs(artifact.delivered_total_impulse_n_s - simpson_impulse)
        duration_change = (
            0.0 if previous_duration is None else abs(artifact.duration_s - previous_duration)
        )
        centroid_change = (
            0.0
            if previous_centroid is None
            else abs(artifact.thrust_centroid_time_s - previous_centroid)
        )
        residuals = (mass_error, impulse_error, duration_change, centroid_change)
        scales = (
            consumed,
            simpson_impulse,
            artifact.duration_s,
            artifact.thrust_centroid_time_s,
        )
        refinement_stable = previous_duration is not None
        if refinement_stable and all(
            error <= max(1e-11, abs(scale) * relative_tolerance)
            for error, scale in zip(residuals, scales, strict=True)
        ):
            break
        previous_duration = artifact.duration_s
        previous_centroid = artifact.thrust_centroid_time_s
        step *= 0.5
    else:
        raise DomainError("Blowdown mesh did not meet the requested convergence tolerance.")
    assert samples is not None and artifact is not None and residuals is not None
    evidence = BlowdownIntegrationEvidence(
        method="mass_domain_composite_simpson_with_tabulation_refinement_v1",
        mass_step_kg=step,
        refinement_count=refinement,
        relative_tolerance=relative_tolerance,
        consumed_propellant_kg=consumed,
        delivered_impulse_n_s=simpson_impulse,
        mass_closure_error_kg=residuals[0],
        impulse_closure_error_n_s=residuals[1],
        duration_change_s=residuals[2],
        centroid_change_s=residuals[3],
    )
    warnings = (
        "Preliminary L2 model: ullage follows a user-declared polytropic exponent.",
        "Pressure-performance exponents require hardware-specific calibration data.",
        (
            "No regulator, line thermal state, dissolved gas, vapor pressure, or diaphragm "
            "mechanics are modeled."
        ),
    )
    return BlowdownSimulationResult(
        samples=samples,
        artifact=artifact,
        evidence=evidence,
        cutoff_reason=cutoff_reason,
        reference_ids=BLOWDOWN_REFERENCE_IDS + ("NASA-20040000363",),
        warnings=warnings,
    )
