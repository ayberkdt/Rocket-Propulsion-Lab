"""Optional reference integration of one liquid feed, a motor and translation.

The twelve states are position, velocity and the six regulated-feed states.
Vehicle mass is reconstructed from liquid inventory and constant nonliquid
mass. SciPy is loaded only when propagation is requested. Profile boundaries
split the integration and protective margins stop the commanded thrust.

References
----------
.. [1] SciPy 1.14.1, DOP853 accepted steps and dense output.
       https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.DOP853.html
.. [2] NASA SP-8112, Pressurization Systems for Liquid Rockets.
       https://ntrs.nasa.gov/citations/19760015212
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import astuple, dataclass, replace
from functools import lru_cache
from importlib import import_module
from itertools import pairwise
from math import hypot, isfinite, nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, real_mapping, vector3
from rocket_propulsion.propulsion.burns.feed_system import (
    CoupledFeedPropulsionEvaluation,
    FeedEventKind,
    RegulatedFeedState,
    RegulatedFeedSystem,
    evaluate_coupled_feed_propulsion,
)
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational, _fraction
from rocket_propulsion.propulsion.burns.propagator import PropagatorPropulsionBridge
from rocket_propulsion.propulsion.burns.tabulated import TabulatedBurnSegment


def _first_margin_crossing(
    start: float,
    end: float,
    sample: Callable[[float], tuple[tuple[float, ...], tuple[float, ...]]],
    brentq: Callable,
) -> tuple[float, int] | None:
    """Scan an accepted step in time order, testing bracketed margin minima."""
    # A fixed absolute tolerance in seconds can span the entire bracket.
    # Let Brent's relative tolerance and the binary64 time grid set resolution.
    root_absolute_tolerance = nextafter(0.0, 1.0)
    points = sorted({start, end, *(start + (end - start) * (i / 4) for i in (1, 2, 3))})
    left_values, left_rates = sample(start)
    for left, right in pairwise(points):
        right_values, right_rates = sample(right)
        roots: list[tuple[float, int]] = []
        for index, right_value in enumerate(right_values):
            if left_values[index] <= 0.0:
                roots.append((left, index))
                continue

            def value(t: float, i: int = index) -> float:
                return sample(t)[0][i]

            def rate(t: float, i: int = index) -> float:
                return sample(t)[1][i]

            bracket_end = right
            bracket_value = right_value
            if left_rates[index] < 0.0 <= right_rates[index]:
                trough = float(brentq(rate, left, right, xtol=root_absolute_tolerance))
                trough_value = value(trough)
                if trough_value <= 0.0:
                    bracket_end, bracket_value = trough, trough_value
            if bracket_value <= 0.0:
                root = float(brentq(value, left, bracket_end, xtol=root_absolute_tolerance))
                roots.append((root, index))
        if roots:
            return min(roots)
        left_values, left_rates = right_values, right_rates
    return None


@dataclass(frozen=True, slots=True)
class FeedTrajectoryState:
    """State at an accepted epoch, with mass derived from liquid inventory.

    Attributes
    ----------
    time_s : float
        Time since ignition, s.
    position_m, velocity_m_s : tuple of float
        Cartesian position and velocity in the caller's fixed inertial frame.
    vehicle_mass_kg : float
        Constant nonliquid mass plus the current liquid mass, kg.
    feed : RegulatedFeedState
        Physical liquid, gas and regulator state [1]_.

    References
    ----------
    .. [1] NASA SP-8112. https://ntrs.nasa.gov/citations/19760015212
    """

    time_s: float
    position_m: tuple[float, float, float]
    velocity_m_s: tuple[float, float, float]
    vehicle_mass_kg: float
    feed: RegulatedFeedState


@dataclass(frozen=True, slots=True)
class RegulatedFeedTrajectory:
    """Final state and first commanded stop in the reference integration.

    Attributes
    ----------
    final_state : FeedTrajectoryState
        State at the requested final epoch.
    cutoff_state : FeedTrajectoryState or None
        State when the first protective margin or nominal profile end stopped
        thrust. None when the requested interval ended during the burn.
    cutoff_id : str or None
        First terminal margin, or ``burn:stop`` for the nominal profile end.
        Equal reported root times use detector order [1]_.
    final_feed_margins : tuple of (str, float)
        Feed margin IDs and signed values at the final state. If thrust has
        stopped, these use zero demand and can expose a later coast violation.

    References
    ----------
    .. [1] SciPy DOP853 accepted steps and dense output.
           https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.DOP853.html
    """

    final_state: FeedTrajectoryState
    cutoff_state: FeedTrajectoryState | None
    cutoff_id: str | None
    final_feed_margins: tuple[tuple[str, float], ...]


def propagate_regulated_feed(
    bridge: PropagatorPropulsionBridge,
    feed_system: RegulatedFeedSystem,
    initial_feed_state: RegulatedFeedState,
    *,
    duration_s: float,
    initial_vehicle_mass_kg: float,
    initial_position_m: tuple[float, float, float],
    initial_velocity_m_s: tuple[float, float, float],
    direction: tuple[float, float, float],
    pressure_condition_name: str,
    max_step_s: float,
    mu_m3_s2: float = 0.0,
    other_operating_conditions: Mapping[str, float] | None = None,
    relative_tolerance: float = 1e-10,
    absolute_tolerance: float = 1e-11,
) -> RegulatedFeedTrajectory:
    """Propagate a pressure-coupled single-tank motor through cutoff and coast.

    Parameters
    ----------
    bridge, feed_system, initial_feed_state
        Production motor bridge, regulated liquid-feed model and initial state.
        The bridge must have a total-flow pressure performance surface and no
        independent inventory bindings, nozzle ambient correction or ignition
        bias. The artifact must have a single propellant stream.
    duration_s : float
        Positive final time after ignition, s.
    initial_vehicle_mass_kg : float
        Total mass greater than initial liquid mass, kg. Pressurant remains
        onboard. There are no modeled leaks, stage drops or other mass sinks.
    initial_position_m, initial_velocity_m_s, direction : tuple of float
        Inertial initial conditions and fixed, nonzero thrust direction.
    pressure_condition_name : str
        Performance axis in Pa driven by the physical injector pressure.
    max_step_s : float
        Positive maximum solver step, s. Choose it for the relevant event and
        regulator bandwidth. Multiple extrema inside a pressure scan interval
        may escape detection.
    mu_m3_s2 : float, default 0.0
        Nonnegative central gravity parameter. Zero gives force-free motion.
    other_operating_conditions : Mapping or None
        Constant remaining performance-axis values.
    relative_tolerance, absolute_tolerance : float
        Positive solve_ivp tolerances. The absolute tolerance applies to each
        component in its state units.

    Returns
    -------
    RegulatedFeedTrajectory
        Final state and first thrust-stop state. Coast advances the regulator
        and gas states at zero liquid demand, retaining central gravity.

    Raises
    ------
    DomainError
        On invalid inputs, unsupported bridge configuration, invalid solver
        stages or integration failure. Missing SciPy requires the dev extra.

    Notes
    -----
    DOP853 supplies dense output per accepted step [1]_. Signed margins and
    their rates are checked at quarter-step spacing. A falling-to-rising
    rate brackets a possible minimum even when both pressures are safe.
    Each profile segment ends with its left limit. The next checks the right
    limit before thrust resumes, including demand-induced pressure jumps.
    Pressure queries beyond a map edge use endpoint continuation solely for
    bracketing; the physical pressure margin still determines the stop.
    Trial liquid masses at or below zero use the limiting empty-tank geometry
    with a positive sentinel solely for root location. Accepted depletion
    states contain the exact configured reserve. Shutdown is instantaneous
    commanded thrust removal. Valve closure, tailoff, attitude and slosh need
    additional models and hardware measurements [2]_.

    References
    ----------
    .. [1] SciPy 1.14.1, DOP853.
           https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.integrate.DOP853.html
    .. [2] NASA SP-8112. https://ntrs.nasa.gov/citations/19760015212
    """

    if not isinstance(bridge, PropagatorPropulsionBridge):
        raise DomainError("Trajectory requires a PropagatorPropulsionBridge.")
    if not isinstance(feed_system, RegulatedFeedSystem):
        raise DomainError("Trajectory requires a RegulatedFeedSystem.")
    if not isinstance(initial_feed_state, RegulatedFeedState):
        raise DomainError("Trajectory requires a RegulatedFeedState.")
    parameters = bridge.resolved_parameters()
    if (bridge.stream_bindings or bridge.state_constraints or bridge.nozzle_exit_area_m2 is not None
            or parameters["ignition_time_bias_s"] != 0.0):
        raise DomainError("Reference feed trajectory requires one inventory and ignition at zero, "
                          "without independent state bindings or nozzle ambient correction.")
    stream_ids = {name for segment in bridge.artifact.segments
                  for name, _ in (*segment.start_stream_mass_flows_kg_s,
                                  *segment.end_stream_mass_flows_kg_s)}
    if len(stream_ids) != 1:
        raise DomainError("Reference feed trajectory requires a single propellant stream.")
    if bridge.performance_surface is None:
        raise DomainError("Reference feed trajectory requires a total-flow performance surface.")
    if not isinstance(pressure_condition_name, str):
        raise DomainError("Pressure condition must name a performance axis.")
    axis = next((axis for axis in bridge.performance_surface.axes
                 if axis.name == pressure_condition_name), None)
    if axis is None or axis.unit != "Pa" or axis.values[0] < 0.0:
        raise DomainError("Pressure condition requires a nonnegative pressure axis in Pa.")
    duration_s = real(duration_s, "Trajectory duration")
    max_step_s = real(max_step_s, "Trajectory maximum step")
    initial_vehicle_mass_kg = real(initial_vehicle_mass_kg, "Initial vehicle mass")
    mu_m3_s2 = real(mu_m3_s2, "Central gravity parameter")
    relative_tolerance = real(relative_tolerance, "Trajectory relative tolerance")
    absolute_tolerance = real(absolute_tolerance, "Trajectory absolute tolerance")
    if min(duration_s, max_step_s, relative_tolerance, absolute_tolerance) <= 0.0:
        raise DomainError("Duration, step and solver tolerances must be positive.")
    if initial_vehicle_mass_kg <= initial_feed_state.propellant_mass_kg or mu_m3_s2 < 0.0:
        raise DomainError("Nonliquid mass must be positive and central gravity nonnegative.")
    position = vector3(initial_position_m, "Initial position")
    velocity = vector3(initial_velocity_m_s, "Initial velocity")
    direction = vector3(direction, "Thrust direction")
    if hypot(*direction) == 0.0 or (mu_m3_s2 > 0.0 and hypot(*position) == 0.0):
        raise DomainError("Thrust direction and a gravitating position must be nonzero.")
    conditions = real_mapping(other_operating_conditions, "Operating conditions")
    if pressure_condition_name in conditions:
        raise DomainError("Pressure condition is owned by the feed-system coupling.")
    bridge.performance_surface.evaluate({**conditions, pressure_condition_name: axis.values[0]})
    try:
        integrators = import_module("scipy.integrate")
        solve_ivp = integrators.solve_ivp
        dop853 = integrators.DOP853
        brentq = import_module("scipy.optimize").brentq
    except ImportError as error:
        raise DomainError("Reference feed trajectory requires SciPy; install the dev extra.") from error

    nonliquid_mass = _fraction(initial_vehicle_mass_kg) - _fraction(
        initial_feed_state.propellant_mass_kg)
    config = feed_system.configuration
    surfaces = feed_system.event_surfaces()
    # Reuse the exact pressure-margin arithmetic with the two map thresholds.
    injector = next(surface for surface in surfaces
                    if surface.kind is FeedEventKind.INJECTOR_PRESSURE)
    map_lower = RegulatedFeedSystem(replace(config, minimum_injector_pressure_pa=axis.values[0]))
    map_upper = RegulatedFeedSystem(replace(config, minimum_injector_pressure_pa=axis.values[-1]))
    ids = [surface.event_id for surface in surfaces]
    ids.extend(("performance:pressure-minimum", "performance:pressure-maximum"))
    dry_mass = bridge.protected_dry_mass_kg
    if dry_mass is not None:
        ids.append("burn:dry-mass")

    def mass(liquid: float) -> float:
        return _finite_rational(nonliquid_mass + _fraction(liquid), "Vehicle mass")

    def dry_margin(liquid: float, limit: float) -> float:
        return _finite_rational(nonliquid_mass + _fraction(liquid) - _fraction(limit),
                                "Dry-mass margin")

    def feed_state(values: Sequence[float], *, trial: bool = False) -> RegulatedFeedState:
        liquid = values[6]
        if trial and liquid <= 0.0:
            liquid = nextafter(0.0, 1.0)
        return RegulatedFeedState(liquid, *values[7:12])

    def saved(time: float, values: Sequence[float]) -> FeedTrajectoryState:
        state = feed_state(values)
        return FeedTrajectoryState(time, vector3(values[:3], "Position"),
                                   vector3(values[3:6], "Velocity"),
                                   mass(state.propellant_mass_kg), state)

    def gravity(values: Sequence[float]) -> tuple[float, float, float]:
        if mu_m3_s2 == 0.0:
            return (0.0, 0.0, 0.0)
        x, y, z = values[:3]
        radius = hypot(x, y, z)
        if not isfinite(radius) or radius == 0.0:
            raise DomainError("Central gravity requires a finite, nonzero radius.")
        scale = -_fraction(mu_m3_s2) / _fraction(radius)**3
        return (_finite_rational(scale * _fraction(x), "Gravity x"),
                _finite_rational(scale * _fraction(y), "Gravity y"),
                _finite_rational(scale * _fraction(z), "Gravity z"))

    def coast_rhs(time: float, values: Sequence[float]) -> tuple[float, ...]:
        gas = feed_system.evaluate(feed_state(values), 0.0)
        return (*values[3:6], *gravity(values), *gas.derivative.as_tuple())

    def integrate_burn(
        rhs: Callable,
        sample: Callable,
        start: float,
        end: float,
        values: tuple[float, ...],
    ) -> tuple[float, tuple[float, ...], int | None]:
        solver = dop853(rhs, start, values, end, max_step=max_step_s,
                        rtol=relative_tolerance, atol=absolute_tolerance)
        while solver.status == "running":
            left = float(solver.t)
            message = solver.step()
            if solver.status == "failed":
                raise DomainError(f"Reference feed trajectory integration failed: {message}")
            right = float(solver.t)
            dense = solver.dense_output()

            @lru_cache(maxsize=32)
            def dense_sample(epoch: float, interpolate: Callable = dense):
                return sample(epoch, tuple(float(value) for value in interpolate(epoch)))

            crossing = _first_margin_crossing(left, right, dense_sample, brentq)
            if crossing is not None:
                epoch, index = crossing
                return epoch, tuple(float(value) for value in dense(epoch)), index
        return float(solver.t), tuple(float(value) for value in solver.y), None

    current = (*position, *velocity, *astuple(initial_feed_state))
    time = 0.0
    cutoff_id: str | None = None
    cutoff_state: FeedTrajectoryState | None = None
    flow_scale = parameters["mass_flow_scale"]
    final_demand = 0.0
    for segment in bridge.artifact.segments:
        end = min(segment.end_time_s, duration_s)
        reference_slope = (
            (_fraction(segment.end_total_mass_flow_kg_s) - _fraction(segment.start_total_mass_flow_kg_s))
            / (_fraction(segment.end_time_s) - _fraction(segment.start_time_s))
        )

        @lru_cache(maxsize=1)
        def coupled(epoch: float, values: tuple[float, ...],
                    active_segment: TabulatedBurnSegment = segment) -> CoupledFeedPropulsionEvaluation:
            state = feed_state(values, trial=True)
            side = "left" if epoch == active_segment.end_time_s else "right"
            reference_flow = bridge.artifact.state_at(epoch, side=side)["total_mass_flow_kg_s"]
            return evaluate_coupled_feed_propulsion(
                bridge, feed_system, state, relative_time_s=epoch,
                vehicle_mass_kg=mass(state.propellant_mass_kg), direction=direction,
                pressure_condition_name=pressure_condition_name, other_operating_conditions=conditions,
                apply_inventory_limits=False, relative_tolerance=1e-12,
                side=side, initial_mass_flow_kg_s=reference_flow * flow_scale,
                for_pressure_event_location=True,
                mass_flow_scale_floor_kg_s=0.0,
            )

        @lru_cache(maxsize=1)
        def margins(epoch: float, values: tuple[float, ...]) -> tuple[float, ...]:
            state = feed_state(values, trial=True)
            evaluation = coupled(epoch, tuple(values))
            demand = -evaluation.feed.derivative.propellant_mass_kg_s
            result = [values[6] - config.minimum_propellant_mass_kg
                      if surface.kind is FeedEventKind.PROPELLANT_RESERVE
                      else feed_system.event_value(surface, state, demand) for surface in surfaces]
            result.extend((map_lower.event_value(injector, state, demand),
                           -map_upper.event_value(injector, state, demand)))
            if dry_mass is not None:
                result.append(dry_margin(values[6], dry_mass))
            return tuple(result)

        def margin_sample(epoch: float, values: tuple[float, ...],
                          active_segment: TabulatedBurnSegment = segment,
                          flow_slope=reference_slope) -> tuple[tuple[float, ...], tuple[float, ...]]:
            result = margins(epoch, values)
            evaluation = coupled(epoch, values)
            state = feed_state(values, trial=True)
            derivative = evaluation.feed.derivative
            liquid_rate = _fraction(derivative.propellant_mass_kg_s)
            gas_rate = _fraction(derivative.tank_pressurant_mass_kg_s)
            temperature_rate = _fraction(derivative.tank_gas_temperature_k_s)
            partials = dict(evaluation.feed.pressure_partials.injector_pressure_wrt_state)
            tank_rate = (
                _fraction(partials["propellant_mass_kg"]) * (liquid_rate if values[6] > 0.0 else 0)
                + _fraction(partials["tank_pressurant_mass_kg"]) * gas_rate
                + _fraction(partials["tank_gas_temperature_k"]) * temperature_rate
            )
            bottle_rate = (
                _fraction(config.gas.specific_gas_constant_j_kg_k)
                / _fraction(config.bottle_volume_m3)
                * (_fraction(state.bottle_gas_temperature_k) * _fraction(derivative.bottle_pressurant_mass_kg_s)
                   + _fraction(state.bottle_pressurant_mass_kg)
                   * _fraction(derivative.bottle_gas_temperature_k_s))
            )
            # q = q_profile(t) * scale(p). Hold p fixed for the explicit
            # time term, then close dp/dt through the same algebraic loop.
            line_term = _fraction(0.0)
            if liquid_rate and flow_slope:
                reference = active_segment.state_at(epoch)["total_mass_flow_kg_s"]
                if reference == 0.0:
                    raise DomainError("Nonzero feed flow requires a positive profile flow.")
                line_term = (_fraction(evaluation.feed.pressure_partials.injector_pressure_wrt_propellant_flow)
                             * flow_slope * (-liquid_rate / _fraction(reference)))
            injector_rate = ((tank_rate + line_term)
                             / _fraction(evaluation.partials.algebraic_denominator))
            by_kind = {
                FeedEventKind.PROPELLANT_RESERVE: liquid_rate,
                FeedEventKind.INJECTOR_PRESSURE: injector_rate,
                FeedEventKind.BOTTLE_PRESSURE_MARGIN: bottle_rate - tank_rate,
                FeedEventKind.TANK_OVERPRESSURE: -tank_rate,
            }
            rates = [by_kind[surface.kind] for surface in surfaces]
            rates.extend((injector_rate, -injector_rate))
            if dry_mass is not None:
                rates.append(liquid_rate)
            return result, tuple(_finite_rational(rate, "Feed cutoff margin rate") for rate in rates)

        # Check inventory and gas-only margins before requesting any liquid.
        if current[6] <= config.minimum_propellant_mass_kg:
            cutoff_id = "feed:propellant-reserve"
        elif dry_mass is not None and dry_margin(current[6], dry_mass) <= 0.0:
            cutoff_id = "burn:dry-mass"
        else:
            cutoff_id = next((surface.event_id for surface in surfaces
                              if surface.kind in (FeedEventKind.BOTTLE_PRESSURE_MARGIN,
                                                  FeedEventKind.TANK_OVERPRESSURE)
                              and feed_system.event_value(surface, feed_state(current), 0.0) <= 0.0),
                             None)
        if cutoff_id is None:
            cutoff_id = next((identifier for identifier, value in zip(ids, margins(time, current),
                                                                     strict=True) if value <= 0.0), None)
        if cutoff_id is not None:
            break
        if time == duration_s:
            final_demand = -coupled(time, current).feed.derivative.propellant_mass_kg_s
            break

        def active_rhs(epoch: float, values: Sequence[float]) -> tuple[float, ...]:
            evaluation = coupled(epoch, tuple(values))
            acceleration = tuple(a + g for a, g in zip(evaluation.propulsion.acceleration_m_s2,
                                                       gravity(values), strict=True))
            return (*values[3:6], *acceleration, *evaluation.feed.derivative.as_tuple())

        time, current, triggered = integrate_burn(active_rhs, margin_sample, time, end, current)
        if triggered is not None:
            cutoff_id = ids[triggered]
            if cutoff_id == "feed:propellant-reserve":
                current = (*current[:6], config.minimum_propellant_mass_kg, *current[7:])
            elif cutoff_id == "burn:dry-mass":
                assert dry_mass is not None
                liquid = _finite_rational(_fraction(dry_mass) - nonliquid_mass,
                                          "Dry-mass liquid boundary")
                current = (*current[:6], liquid, *current[7:])
            break
        if any(value <= 0.0 for value in margins(time, current)):
            raise DomainError("Trajectory solver reached a violated burn margin without a root.")
        if time == duration_s and time != segment.end_time_s:
            final_demand = -coupled(time, current).feed.derivative.propellant_mass_kg_s
            break

    if cutoff_id is None and time >= bridge.artifact.duration_s:
        cutoff_id = "burn:stop"
    if cutoff_id is not None:
        cutoff_state = saved(time, current)
    if time < duration_s:
        coast = solve_ivp(coast_rhs, (time, duration_s), current, method="DOP853",
                          max_step=max_step_s, rtol=relative_tolerance, atol=absolute_tolerance)
        if not coast.success:
            raise DomainError("Reference feed trajectory integration failed: " + coast.message)
        current = tuple(float(value) for value in coast.y[:, -1])
    final = saved(duration_s, current)
    final_margins = tuple((surface.event_id, feed_system.event_value(surface, final.feed, final_demand))
                          for surface in surfaces)
    return RegulatedFeedTrajectory(final, cutoff_state, cutoff_id, final_margins)
