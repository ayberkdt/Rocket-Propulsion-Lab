"""Pressure-fed blowdown propulsion model.

This preliminary L2 model resolves rigid-tank ullage expansion under a declared
polytropic gas law [1]_, pressure-correlated thrust and specific impulse
calibrated against hardware data [3]_, and a mesh-converged tabulated
force/mass-flow history.  The unregulated monopropellant time histories of [2]_
illustrate the pressure decay this model represents.  The module is
propulsion-only and contains no orbit or reference-frame logic.

References
----------
.. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
       https://ntrs.nasa.gov/citations/19760015212
.. [2] NASA CR-131400. NASA, "Reliability Model of a Monopropellant Auxiliary
       Propulsion System". https://ntrs.nasa.gov/citations/19730012094
.. [3] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
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
    """Structural contract for a supply-pressure-to-operating-point model.

    Implementations declare an inclusive supply-pressure validity interval
    and refuse extrapolation outside it.  The blowdown integrator depends
    only on this contract, so performance may come from a compact calibrated
    correlation or from a measured hot-fire map [1]_, [2]_.

    Attributes
    ----------
    minimum_supply_pressure_pa : float
        Inclusive lower validity limit of the model, Pa.

    References
    ----------
    .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    # Read-only, so frozen dataclass implementations satisfy the protocol.
    @property
    def minimum_supply_pressure_pa(self) -> float:
        """Inclusive lower supply-pressure validity limit, Pa."""

        ...

    @property
    def maximum_pressure_pa(self) -> float:
        """Inclusive upper supply-pressure validity limit, Pa."""

        ...

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Return a stream-closed operating point at a supply pressure.

        Parameters
        ----------
        rated_point : EngineOperatingPoint
            Reference operating point whose stream proportions are retained.
        supply_pressure_pa : float
            Tank supply pressure inside the declared validity interval, Pa.

        Returns
        -------
        EngineOperatingPoint
            Operating point at ``supply_pressure_pa``.

        Raises
        ------
        DomainError
            If ``supply_pressure_pa`` lies outside the validity interval.
        """

        ...


def _positive(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and strictly positive."""
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and non-negative."""
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
    """Rescale a rated point to a new thrust and system specific impulse.

    Total tank flow follows the specific-impulse definition [1]_,
    ``mdot = F / (g0 * Isp)`` with ``g0`` the registered
    ``STANDARD_GRAVITY_M_S2``.  Every stream and the chamber flow are scaled by
    the same ratio ``mdot / mdot_rated``, so stream proportions, mixture ratio
    and tank-flow closure are preserved; ideal thrust scales with delivered
    thrust and chamber specific impulse is recomputed from the scaled chamber
    flow.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
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

    Parameters
    ----------
    initial_pressure_pa : float
        Initial ullage pressure ``p0``, Pa.
    initial_ullage_volume_m3 : float
        Initial ullage volume ``V0``, m^3.
    initial_propellant_mass_kg : float
        Initial liquid propellant mass, kg.
    propellant_density_kg_m3 : float
        Constant liquid density, kg/m^3.
    polytropic_exponent : float, default 1.0
        Declared ullage polytropic exponent ``n``.

    Raises
    ------
    DomainError
        If any input is non-finite or not strictly positive.

    Notes
    -----
    The ullage follows the blowdown pressure law [1]_
    ``p = p0 * (V0 / V)**n``, where ``V`` grows by the expelled liquid volume.
    ``n = 1`` is the isothermal ideal-gas limit; a larger declared exponent
    approximates reduced gas-wall heat transfer.  Propellant density is
    constant.  Dissolved gas, diaphragm mechanics, vapour pressure, line losses
    and thermal soak are outside the model.  The unregulated monopropellant
    histories of [2]_ show the pressure decay this model idealizes.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA CR-131400. NASA, "Reliability Model of a Monopropellant
           Auxiliary Propulsion System".
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
        """Initial liquid volume ``m0 / rho`` of the incompressible propellant, m^3 [1]_.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               1975. https://ntrs.nasa.gov/citations/19760015212
        """

        return self.initial_propellant_mass_kg / self.propellant_density_kg_m3

    @property
    def tank_volume_m3(self) -> float:
        """Rigid internal tank volume, initial ullage plus initial liquid, m^3 [1]_.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               1975. https://ntrs.nasa.gov/citations/19760015212
        """

        return self.initial_ullage_volume_m3 + self.initial_propellant_volume_m3

    def ullage_volume_m3(self, remaining_propellant_mass_kg: float) -> float:
        """Return the ullage volume after part of the propellant is expelled.

        Parameters
        ----------
        remaining_propellant_mass_kg : float
            Liquid mass still in the tank, within ``[0, initial mass]``, kg.

        Returns
        -------
        float
            ``V0 + (m0 - m) / rho``, m^3; a rigid tank gains the expelled
            liquid volume as ullage [1]_.

        Raises
        ------
        DomainError
            If the remaining mass is non-finite or outside the initial inventory.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               1975. https://ntrs.nasa.gov/citations/19760015212
        """

        self._validate_remaining_mass(remaining_propellant_mass_kg)
        expelled_mass = self.initial_propellant_mass_kg - remaining_propellant_mass_kg
        return self.initial_ullage_volume_m3 + expelled_mass / self.propellant_density_kg_m3

    def pressure_pa(self, remaining_propellant_mass_kg: float) -> float:
        """Evaluate the declared polytropic ullage-pressure law.

        Parameters
        ----------
        remaining_propellant_mass_kg : float
            Liquid mass still in the tank, within ``[0, initial mass]``, kg.

        Returns
        -------
        float
            Ullage pressure ``p0 * (V0 / V)**n``, Pa.

        Raises
        ------
        DomainError
            If the remaining mass is non-finite or outside the initial inventory.

        Notes
        -----
        This is a preliminary thermodynamic closure [1]_, not a heat-transfer
        or pressurant-dissolution model.  The exponent ``n`` is an explicit
        input rather than an inferred gas property; [2]_ documents the pressure
        decay of an unregulated system of this kind.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               1975. https://ntrs.nasa.gov/citations/19760015212
        .. [2] NASA CR-131400. NASA, "Reliability Model of a Monopropellant
               Auxiliary Propulsion System".
               https://ntrs.nasa.gov/citations/19730012094
        """

        ullage = self.ullage_volume_m3(remaining_propellant_mass_kg)
        return self.initial_pressure_pa * (
            self.initial_ullage_volume_m3 / ullage
        ) ** self.polytropic_exponent

    def remaining_mass_at_pressure(self, pressure_pa: float) -> float:
        """Invert the polytropic law to locate a pressure-triggered cutoff.

        Parameters
        ----------
        pressure_pa : float
            Ullage pressure at which the cutoff occurs, Pa.

        Returns
        -------
        float
            Remaining liquid mass at that pressure, kg.  The required ullage is
            ``V0 * (p0 / p)**(1/n)`` [1]_; the result is clipped to
            ``[0, initial mass]`` and equals the initial mass for
            ``pressure_pa >= p0``.

        Raises
        ------
        DomainError
            If ``pressure_pa`` is non-finite or not strictly positive.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               1975. https://ntrs.nasa.gov/citations/19760015212
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
        """Raise `DomainError` unless the mass lies in ``[0, initial mass]``."""
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
    """Declared power-law pressure correlation around a reference operating point.

    Parameters
    ----------
    reference_supply_pressure_pa : float
        Supply pressure ``p_ref`` of the rated operating point, Pa.
    minimum_supply_pressure_pa : float
        Inclusive lower validity limit, Pa.
    thrust_pressure_exponent : float, default 1.0
        Non-negative thrust exponent ``a``.
    specific_impulse_pressure_exponent : float, default 0.0
        Finite specific-impulse exponent ``b``.
    chamber_pressure_exponent : float, default 1.0
        Non-negative chamber-pressure exponent ``c``.
    maximum_supply_pressure_pa : float or None, default None
        Inclusive upper validity limit, Pa; ``None`` uses ``p_ref``.

    Raises
    ------
    DomainError
        If a pressure is not finite and positive, the validity interval is
        empty, ``p_ref`` lies outside it, or an exponent violates its sign
        or finiteness requirement.

    Notes
    -----
    With ``r = p / p_ref`` the correlation is ``F / F_ref = r**a``,
    ``Isp / Isp_ref = r**b`` and ``pc / pc_ref = r**c``.  The exponents are
    hardware-specific calibration inputs, not universal rocket laws.  [1]_
    (sections 3.2 and 3.3) reports a nearly linear thrust decline with supply
    pressure for its tested device and only a small Isp variation over a
    200-to-100 psi blowdown; other hardware requires its own data fit.  Thrust
    and specific impulse remain linked through the thrust equation [2]_, [3]_.

    References
    ----------
    .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [2] NASA SP-125. NASA, "Design of Liquid Propellant Rocket Engines".
           https://ntrs.nasa.gov/citations/19710019929
    .. [3] NASA Glenn Research Center, "Rocket Thrust Equation".
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
        """Inclusive upper supply-pressure validity limit, Pa.

        Equals ``maximum_supply_pressure_pa`` when declared, otherwise the
        reference supply pressure; the calibration range follows the tested
        pressure span of the source data [1]_.

        References
        ----------
        .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        """

        return self.maximum_supply_pressure_pa or self.reference_supply_pressure_pa

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Scale a stream-closed point while preserving the system-Isp identity.

        Parameters
        ----------
        rated_point : EngineOperatingPoint
            Operating point measured or modelled at the reference pressure.
        supply_pressure_pa : float
            Tank supply pressure, Pa.  Values within a relative ``1e-12``
            (absolute ``1e-9`` Pa) of the validity limits are clamped to them.

        Returns
        -------
        EngineOperatingPoint
            Operating point with thrust, system Isp and chamber pressure scaled
            by the declared power laws [1]_.

        Raises
        ------
        DomainError
            If ``supply_pressure_pa`` is non-finite or outside the validity
            interval beyond the clamping tolerance.

        Notes
        -----
        All tank-depleting and chamber streams retain their rated proportions.
        Total flow is recomputed as ``F / (g0 * Isp)`` [2]_ with ``g0`` the
        registered ``STANDARD_GRAVITY_M_S2``; it is never inferred from chamber
        flow alone.

        References
        ----------
        .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        .. [2] NASA Glenn Research Center, "Rocket Thrust Equation".
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
    """One measured or validated row of a pressure-performance map.

    Parameters
    ----------
    supply_pressure_pa : float
        Tank supply pressure, Pa.
    delivered_thrust_n : float
        Delivered thrust at that pressure, N.
    system_specific_impulse_s : float
        System specific impulse at that pressure, s.
    chamber_pressure_pa : float
        Chamber pressure at that pressure, Pa.

    Raises
    ------
    DomainError
        If any value is non-finite or not strictly positive.

    Notes
    -----
    Values are absolute SI quantities rather than normalized multipliers, so
    the map is traceable to its test or model source [1]_, [2]_ without
    access to the rated operating point.

    References
    ----------
    .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [2] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
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

    Parameters
    ----------
    samples : tuple of PressurePerformanceSample
        At least two rows with strictly increasing supply pressure.
    source_id : str
        Non-blank identifier of the test or model source of the rows.

    Raises
    ------
    DomainError
        If fewer than two rows are supplied, the source is blank, or the
        pressures are not strictly increasing.

    Notes
    -----
    The rows come from hot-fire tests or a validated engine model [1]_, [2]_.
    Each quantity is interpolated linearly in supply pressure between the
    bracketing rows.  The table refuses extrapolation and labels the resulting
    operating point with ``source_id``.  [3]_ gives the transient-modelling
    context in which such component maps are used.

    References
    ----------
    .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [2] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    .. [3] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
        """Inclusive lower validity limit: the first row's supply pressure, Pa [1]_.

        References
        ----------
        .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        """

        return self.samples[0].supply_pressure_pa

    @property
    def maximum_pressure_pa(self) -> float:
        """Inclusive upper validity limit: the last row's supply pressure, Pa [1]_.

        References
        ----------
        .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        """

        return self.samples[-1].supply_pressure_pa

    def scale_operating_point(
        self, rated_point: EngineOperatingPoint, supply_pressure_pa: float
    ) -> EngineOperatingPoint:
        """Interpolate the map and preserve tank-stream closure.

        Parameters
        ----------
        rated_point : EngineOperatingPoint
            Operating point whose stream proportions are retained.
        supply_pressure_pa : float
            Tank supply pressure, Pa.  Values within a relative ``1e-12``
            (absolute ``1e-9`` Pa) of the table limits are clamped to them.

        Returns
        -------
        EngineOperatingPoint
            Operating point carrying the linearly interpolated thrust, system
            Isp and chamber pressure [1]_.

        Raises
        ------
        DomainError
            If ``supply_pressure_pa`` is non-finite or outside the table range
            beyond the clamping tolerance.

        Notes
        -----
        Total flow is recomputed as ``F / (g0 * Isp)`` [2]_ from the
        interpolated values, with ``g0`` the registered
        ``STANDARD_GRAVITY_M_S2``, and all streams keep their rated proportions.

        References
        ----------
        .. [1] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        .. [2] NASA Glenn Research Center, "Rocket Thrust Equation".
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
    """One state of a pressure-fed blowdown history.

    A sequence of samples records the pressure decay of an unregulated system
    of the kind documented in [1]_.

    Attributes
    ----------
    time_s : float
        Elapsed burn time, s.
    remaining_propellant_mass_kg : float
        Liquid mass still in the tank, kg.
    ullage_volume_m3 : float
        Ullage volume, m^3.
    supply_pressure_pa : float
        Ullage (supply) pressure, Pa.
    delivered_thrust_n : float
        Delivered thrust magnitude, N.
    axial_thrust_n : float
        Delivered thrust multiplied by the cant efficiency, N.
    system_specific_impulse_s : float
        System specific impulse, s.
    total_mass_flow_kg_s : float
        Total tank outflow, kg/s.

    References
    ----------
    .. [1] NASA CR-131400. NASA, "Reliability Model of a Monopropellant
           Auxiliary Propulsion System".
           https://ntrs.nasa.gov/citations/19730012094
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
    """Mesh-refinement and closure evidence for a blowdown history.

    The residuals document the numerical closure of the integration; [1]_
    gives the transient-modelling context.

    Attributes
    ----------
    method : str
        Identifier of the integration and refinement scheme.
    mass_step_kg : float
        Expelled-mass step of the accepted mesh, kg.
    refinement_count : int
        Number of step halvings performed.
    relative_tolerance : float
        Requested relative tolerance.
    consumed_propellant_kg : float
        Propellant expelled between start and cutoff, kg.
    delivered_impulse_n_s : float
        Composite-Simpson delivered impulse, N s.
    mass_closure_error_kg : float
        Absolute difference between tabulated and expected consumed mass, kg.
    impulse_closure_error_n_s : float
        Absolute difference between tabulated and Simpson impulse, N s.
    duration_change_s : float
        Burn-duration change from the previous refinement, s.
    centroid_change_s : float
        Thrust-centroid-time change from the previous refinement, s.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
    """Blowdown profile, exchange artifact and source labels of one simulation.

    The history combines the ullage law of [1]_ with a pressure-performance
    calibration of the kind reported in [3]_; [2]_ documents comparable
    unregulated histories.

    Attributes
    ----------
    samples : tuple of BlowdownSample
        State history on the accepted expelled-mass mesh.
    artifact : TabulatedBurnArtifact
        Piecewise-linear force and mass-flow exchange artifact.
    evidence : BlowdownIntegrationEvidence
        Closure and mesh-refinement evidence.
    cutoff_reason : str
        ``"minimum_supply_pressure"`` or ``"propellant_reserve"``.
    reference_ids : tuple of str
        Identifiers of the sources behind the model.
    warnings : tuple of str
        Model-scope warnings.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA CR-131400. NASA, "Reliability Model of a Monopropellant
           Auxiliary Propulsion System".
           https://ntrs.nasa.gov/citations/19730012094
    .. [3] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    """

    samples: tuple[BlowdownSample, ...]
    artifact: TabulatedBurnArtifact
    evidence: BlowdownIntegrationEvidence
    cutoff_reason: str
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]


def _tank_flows(point: EngineOperatingPoint) -> tuple[tuple[str, float], ...]:
    """Sum tank-depleting stream flows per tank, sorted by tank identifier."""
    flows: dict[str, float] = {}
    for stream in point.streams:
        if stream.tank_depleting:
            flows[stream.tank_id] = flows.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s
    return tuple(sorted(flows.items()))


def _simpson(start: float, middle: float, end: float, width: float) -> float:
    """Return the one-panel Simpson integral ``h * (f0 + 4 * f_mid + f1) / 6``."""
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
    """Build one uniform expelled-mass mesh of the blowdown history.

    The consumed mass is split into ``ceil(consumed / maximum_mass_step_kg)``
    equal steps.  On each step the duration is the Simpson integral of
    ``dt/dm = 1 / mdot`` and the impulse the Simpson integral of
    ``dJ/dm = g0 * Isp`` [2]_, with ``g0`` the registered
    ``STANDARD_GRAVITY_M_S2``; pressure at every node follows the polytropic
    tank law [1]_.  Each step becomes one linear segment of the exchange
    artifact.

    Returns
    -------
    tuple
        ``(samples, artifact, simpson_impulse_n_s)``.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
    consumed =tank.initial_propellant_mass_kg - final_remaining_mass_kg
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

    Parameters
    ----------
    tank : BlowdownTankModel
        Rigid tank with polytropic ullage.
    rated_point : EngineOperatingPoint
        Operating point at the pressure model's reference pressure.
    pressure_law : PressurePerformanceModel
        Supply-pressure-to-performance model with a declared validity range.
    reserve_mass_kg : float, default 0.0
        Non-negative propellant mass left unburned, kg.
    cant_efficiency : float, default 1.0
        Axial fraction of delivered thrust, in ``(0, 1]``.
    maximum_mass_step_kg : float, default 0.1
        Initial upper bound on the expelled-mass step, kg.
    relative_tolerance : float, default 1e-7
        Relative tolerance on the mass, impulse, duration and centroid
        residuals (absolute floor ``1e-11``).
    maximum_refinements : int, default 12
        Maximum number of step halvings.

    Returns
    -------
    BlowdownSimulationResult
        Sampled history, exchange artifact, convergence evidence, cutoff
        reason, source identifiers and model warnings.

    Raises
    ------
    DomainError
        If an argument is outside its domain, the initial tank pressure lies
        outside the pressure model's validity range, the declared limits leave
        no usable propellant, or the mesh does not converge within
        ``maximum_refinements`` halvings.

    Notes
    -----
    Cutoff occurs at the larger of the reserve mass and the remaining mass at
    which the polytropic ullage pressure [1]_ reaches the model's minimum
    supply pressure.  Expelled propellant mass ``m`` is the independent
    variable, ``dt/dm = 1 / mdot`` and ``dJ/dm = g0 * Isp`` [4]_, with ``g0``
    the registered ``STANDARD_GRAVITY_M_S2``.  Composite Simpson integration
    supplies the reference duration and impulse; the exchange artifact is
    integrated independently, and the mass step is halved until the mass,
    impulse, duration and centroid residuals satisfy the requested tolerance.

    The model is preliminary.  Pressure-performance exponents are fitted to a
    thruster map or hot-fire data [3]_; the function refuses extrapolation.
    [2]_ documents the pressure and performance histories of an unregulated
    monopropellant system, and [5]_ the transient-modelling context.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA CR-131400. NASA, "Reliability Model of a Monopropellant
           Auxiliary Propulsion System".
           https://ntrs.nasa.gov/citations/19730012094
    .. [3] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [4] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [5] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
