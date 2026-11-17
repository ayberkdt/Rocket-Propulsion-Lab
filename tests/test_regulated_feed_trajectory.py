"""End-to-end feed cutoffs against analytic mass, impulse and coast oracles."""

import math
from dataclasses import replace

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.integrations.regulated_feed_trajectory import propagate_regulated_feed
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    evaluate_coupled_feed_propulsion,
)
from rocket_propulsion.propulsion.burns.feed_system import RegulatedFeedSystem
from tests.test_burn_feed_system import example_feed_system
from tests.test_burn_propagator import example_artifact


def fixture(*, minimum_pressure=1.5e6, reserve=0.0, resistance=1000.0,
            pressure_grid=(1e5, 3e6), rates=(1.0,), force=50.0, segment_duration=20.0,
            pressure_dependent=False, dry_mass=None):
    system = RegulatedFeedSystem(replace(
        example_feed_system().configuration,
        gas=replace(example_feed_system().configuration.gas, gamma=1.4),
        regulator_set_pressure_pa=1e5,
        minimum_injector_pressure_pa=minimum_pressure,
        minimum_propellant_mass_kg=reserve,
        feed_line_resistance_pa_s2_kg2=resistance,
    ))
    state = system.initial_state(
        propellant_mass_kg=10.0, tank_pressure_pa=2e6, tank_gas_temperature_k=300.0,
        bottle_pressure_pa=20e6, bottle_gas_temperature_k=300.0, regulator_opening=0.0,
    )
    segments = tuple(TabulatedBurnSegment(
        start_time_s=i * segment_duration, end_time_s=(i + 1) * segment_duration,
        start_delivered_thrust_n=force, end_delivered_thrust_n=force,
        start_axial_thrust_n=force, end_axial_thrust_n=force,
        start_total_mass_flow_kg_s=rate, end_total_mass_flow_kg_s=rate,
        start_stream_mass_flows_kg_s=(("propellant", rate),),
        end_stream_mass_flows_kg_s=(("propellant", rate),), phase="firing",
    ) for i, rate in enumerate(rates))
    surface = RectilinearPerformanceSurface(
        (PerformanceAxis("injector_pa", "Pa", pressure_grid),),
        tuple(PerformanceScalePoint((pressure,), pressure / 2e6 if pressure_dependent else 1.0,
                                    pressure / 2e6 if pressure_dependent else 1.0)
              for pressure in pressure_grid),
        source_id="analytic-test-field", source_sha256="a" * 64,
    )
    bridge = PropagatorPropulsionBridge(
        TabulatedBurnArtifact(segments, source_id="analytic-test-command", source_sha256="b" * 64),
        performance_surface=surface, require_complete_stream_binding=False,
        protected_dry_mass_kg=dry_mass,
    )
    return bridge, system, state


def run_case(bridge, system, state, duration=5.0, **options):
    arguments = {
        "duration_s": duration, "initial_vehicle_mass_kg": 100.0,
        "initial_position_m": (0.0, 0.0, 0.0), "initial_velocity_m_s": (0.0, 0.0, 0.0),
        "direction": (1.0, 0.0, 0.0), "pressure_condition_name": "injector_pa", "max_step_s": 0.17,
    }
    arguments.update(options)
    try:
        return propagate_regulated_feed(bridge, system, state, **arguments)
    except DomainError as error:
        raise AssertionError("A physical burn, event and coast must propagate.") from error


def constant_flow_oracle(cutoff, *, initial_mass=100.0, rate=1.0, force=50.0, final=5.0):
    remaining = initial_mass - rate * cutoff
    speed = force / rate * math.log(initial_mass / remaining)
    position = force / rate * (cutoff + remaining / rate * math.log(remaining / initial_mass))
    return remaining, speed, position + speed * (final - cutoff)


@pytest.mark.parametrize("max_step", [0.07, 0.83])
def test_pressure_cutoff_and_coast_match_sealed_adiabatic_tank_and_rocket_equation(max_step):
    bridge, system, state = fixture()
    stop = 10.0 * ((2e6 / (1.5e6 + 1000.0)) ** (1 / 1.4) - 1)
    result = run_case(bridge, system, state, max_step_s=max_step)
    mass, speed, position = constant_flow_oracle(stop)
    assert result.cutoff_id == "feed:injector-pressure"
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=2e-9)
    assert result.cutoff_state.feed.propellant_mass_kg == pytest.approx(10.0 - stop, abs=2e-9)
    assert result.final_state.vehicle_mass_kg == pytest.approx(mass, abs=2e-9)
    assert result.final_state.feed == result.cutoff_state.feed
    assert result.final_state.velocity_m_s == pytest.approx((speed, 0.0, 0.0), abs=2e-9)
    assert result.final_state.position_m == pytest.approx((position, 0.0, 0.0), abs=2e-8)


def test_pressure_dependent_motor_flow_changes_cutoff_time_and_delivered_impulse():
    bridge, system, state = fixture(resistance=0.0, pressure_dependent=True)
    # q = q0 (V0/V)^gamma. Integrating dV/dt gives V^(gamma+1) affine in time.
    stop = 10.0 / 2.4 * ((2e6 / 1.5e6) ** (2.4 / 1.4) - 1)
    liquid_at_stop = 20.0 - 10.0 * (2e6 / 1.5e6) ** (1 / 1.4)
    result = run_case(bridge, system, state)
    assert result.cutoff_id == "feed:injector-pressure"
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=2e-9)
    assert result.final_state.feed.propellant_mass_kg == pytest.approx(liquid_at_stop, abs=2e-9)
    expected_speed = 50.0 * math.log(100.0 / (90.0 + liquid_at_stop))
    assert result.final_state.velocity_m_s[0] == pytest.approx(expected_speed, abs=2e-9)


@pytest.mark.parametrize(("max_step", "depth_pa", "competing_stop_s", "rtol", "time_error", "time_scale"), [
    pytest.param(0.01, 10.0, None, 1e-10, 2e-9, 1.0, id="fine"),
    pytest.param(1.0, 10.0, None, 1e-10, 1e-8, 1.0, id="coarse"),
    pytest.param(1.0, 0.01, None, 1e-10, 5e-7, 1.0, id="grazing"),
    pytest.param(1.0, 0.01, None, 1e-13, 1e-8, 1.0, id="grazing-tight-tolerance"),
    pytest.param(1.0, 0.01, 0.39, 1e-10, 5e-7, 1.0, id="grazing-before-dry-mass"),
    pytest.param(1.0, -0.01, None, 1e-10, 2e-9, 1.0, id="near-miss"),
    pytest.param(1.0, 10.0, None, 1e-10, 1e-8, 1e-15, id="rescaled-time"),
])
def test_pressure_cutoff_during_underpressure_excursion_with_safe_segment_endpoints(
    max_step, depth_pa, competing_stop_s, rtol, time_error, time_scale,
):
    from scipy.integrate import quad
    from scipy.optimize import brentq

    def ullage_as_displaced_liquid_kg(t):
        return 10.0 + t - t**2 / 2.0

    def injector_pressure(t):
        # Sealed adiabatic tank, q = (1 - t) kg/s, rho*V0 = 10 kg.
        gas_pressure = 2e6 * (10.0 / ullage_as_displaced_liquid_kg(t))**1.4
        return gas_pressure - 130_000.0 * (1.0 - t)**2

    minimum = brentq(
        lambda t: 260_000.0 - 1.4 * 2e6 * (10.0 / ullage_as_displaced_liquid_kg(t))**1.4
        / ullage_as_displaced_liquid_kg(t),
        0.0, 0.9,
    )
    assert minimum == pytest.approx(0.389426534333, abs=1e-12)
    limit = injector_pressure(minimum) + depth_pa

    def margin(t):
        return injector_pressure(t) - limit

    assert margin(0.0) > 0.0 and margin(1.0) > 0.0
    if depth_pa > 0.0:
        stop = brentq(margin, 0.0, minimum)
        recovery = brentq(margin, minimum, 1.0)
        assert recovery - stop < 0.06
    else:
        stop = 1.0
    dry_mass = None
    if competing_stop_s is not None:
        assert stop < competing_stop_s < recovery
        dry_mass = 110.0 - ullage_as_displaced_liquid_kg(competing_stop_s)

    # Changing the time scale preserves liquid loss, pressure history in
    # normalized time, and impulse. The shortest case probes root numerics.
    bridge, system, state = fixture(minimum_pressure=limit, resistance=130_000.0 * time_scale**2,
                                    segment_duration=time_scale, dry_mass=dry_mass,
                                    rates=(1.0 / time_scale,), force=50.0 / time_scale)
    segment = replace(bridge.artifact.segments[0], end_total_mass_flow_kg_s=0.0,
                      end_stream_mass_flows_kg_s=(("propellant", 0.0),))
    bridge = replace(bridge, artifact=replace(bridge.artifact, segments=(segment,)))
    result = run_case(bridge, system, state, duration=time_scale, max_step_s=max_step * time_scale,
                      relative_tolerance=rtol, absolute_tolerance=rtol / 10.0)
    assert result.cutoff_id == ("feed:injector-pressure" if depth_pa > 0.0 else "burn:stop")
    # A grazing root amplifies pressure interpolation errors by 1/abs(dp/dt).
    assert result.cutoff_state.time_s / time_scale == pytest.approx(stop, abs=time_error, rel=0.0)
    expected_liquid = 20.0 - ullage_as_displaced_liquid_kg(stop)
    # Per normalized second, flow is at most 1 kg and speed changes by at
    # most 0.503 m/s. These bound the propagated error due to root timing.
    assert result.final_state.feed.propellant_mass_kg == pytest.approx(expected_liquid, abs=time_error, rel=0.0)

    def acceleration(t):
        return 50.0 / (110.0 - ullage_as_displaced_liquid_kg(t))

    speed = quad(acceleration, 0.0, stop)[0]
    position = quad(lambda t: (1.0 - t) * acceleration(t), 0.0, stop)[0]
    assert result.final_state.velocity_m_s == pytest.approx((speed, 0.0, 0.0), abs=time_error, rel=0.0)
    normalized_position = tuple(value / time_scale for value in result.final_state.position_m)
    assert normalized_position == pytest.approx((position, 0.0, 0.0), abs=time_error, rel=0.0)


def test_bottle_pressure_margin_cuts_at_first_of_two_thermal_crossings():
    from scipy.optimize import brentq

    # Closed, fixed-volume ideal gas in each vessel. The bottle cools faster
    # than the tank and the pressure difference has an interior minimum.
    minimum_time = math.log(8.0) / 1.5

    def head(t):
        return 17e6 + 2e6 * math.exp(-2.0 * t) - 1e6 * math.exp(-0.5 * t)

    limit = head(minimum_time) + 0.01
    assert head(0.0) > limit and head(4.0) > limit
    stop = brentq(lambda t: head(t) - limit, 0.0, minimum_time)
    bridge, system, state = fixture(rates=(0.0,), force=0.0, minimum_pressure=1e5)
    bottle_volume = system.configuration.bottle_volume_m3
    config = replace(system.configuration, minimum_bottle_pressure_margin_pa=limit,
                     tank_wall_temperature_k=150.0, bottle_wall_temperature_k=270.0,
                     tank_heat_transfer_w_k=0.5 * 2e6 * 0.01 / (0.4 * 300.0),
                     bottle_heat_transfer_w_k=2.0 * 20e6 * bottle_volume / (0.4 * 300.0))
    result = run_case(bridge, RegulatedFeedSystem(config), state, duration=4.0,
                      max_step_s=4.0, relative_tolerance=1e-13, absolute_tolerance=1e-14)
    assert result.cutoff_id == "feed:bottle-pressure-margin"
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=1e-7, rel=0.0)
    assert result.final_state.feed.tank_gas_temperature_k == pytest.approx(
        150.0 + 150.0 * math.exp(-2.0), abs=1e-9, rel=0.0)
    assert result.final_state.feed.bottle_gas_temperature_k == pytest.approx(
        270.0 + 30.0 * math.exp(-8.0), abs=1e-9, rel=0.0)


def test_pressure_rate_turns_with_ramped_command_and_pressure_dependent_motor():
    from scipy.integrate import quad, solve_ivp
    from scipy.optimize import brentq, minimize_scalar

    resistance = 140_000.0

    def pressure_and_flow(t, used):
        tank = 2e6 * (10.0 / (10.0 + used))**1.4
        a = (1.0 - t) / 2e6
        # q = a p and p = p_tank - R q². Use the positive quadratic root.
        flow = 2.0 * a * tank / (1.0 + math.sqrt(1.0 + 4.0 * a * a * resistance * tank))
        return tank - resistance * flow**2, flow

    oracle = solve_ivp(lambda t, y: [pressure_and_flow(t, y[0])[1]], (0.0, 1.0), [0.0],
                       method="RK45", max_step=0.002, rtol=3e-14, atol=1e-15,
                       dense_output=True)
    assert oracle.success

    def pressure(t):
        return pressure_and_flow(t, oracle.sol(t)[0])[0]

    minimum = minimize_scalar(pressure, bounds=(0.0, 1.0), method="bounded")
    assert minimum.success
    limit = minimum.fun + 0.05
    assert pressure(0.0) > limit and pressure(1.0) > limit
    stop = brentq(lambda t: pressure(t) - limit, 0.0, minimum.x)
    bridge, system, state = fixture(minimum_pressure=limit, resistance=resistance,
                                    segment_duration=1.0, pressure_dependent=True)
    segment = replace(bridge.artifact.segments[0], end_total_mass_flow_kg_s=0.0,
                      end_stream_mass_flows_kg_s=(("propellant", 0.0),))
    bridge = replace(bridge, artifact=replace(bridge.artifact, segments=(segment,)))
    result = run_case(bridge, system, state, duration=1.0, max_step_s=1.0,
                      relative_tolerance=1e-12, absolute_tolerance=1e-13)
    assert result.cutoff_id == "feed:injector-pressure"
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=1e-8, rel=0.0)
    assert result.final_state.vehicle_mass_kg == pytest.approx(
        100.0 - oracle.sol(stop)[0], abs=1e-8, rel=0.0)
    speed = quad(lambda t: 50.0 * pressure(t) / (2e6 * (100.0 - oracle.sol(t)[0])), 0.0, stop)[0]
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=1e-8, rel=0.0)


@pytest.mark.parametrize("duration", [1.0, 4.0])
def test_demand_jump_trips_pressure_cutoff_at_exact_profile_boundary(duration, monkeypatch):
    bridge, system, state = fixture(rates=(0.5, 2.0), segment_duration=1.0, resistance=2e5)

    def guarded_evaluation(*args, **kwargs):
        assert kwargs["relative_time_s"] <= 1.0, "Unsafe command reached a later motor stage."
        return evaluate_coupled_feed_propulsion(*args, **kwargs)

    monkeypatch.setattr(
        "rocket_propulsion.integrations.regulated_feed_trajectory.evaluate_coupled_feed_propulsion",
        guarded_evaluation,
    )
    result = run_case(bridge, system, state, duration=duration)
    assert result.cutoff_id == "feed:injector-pressure"
    assert result.cutoff_state.time_s == 1.0
    assert result.final_state.feed.propellant_mass_kg == pytest.approx(9.5, abs=2e-10)
    _, speed, position = constant_flow_oracle(1.0, rate=0.5, final=duration)
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=2e-9)
    assert result.final_state.position_m[0] == pytest.approx(position, abs=2e-9)


def test_safe_horizon_at_profile_boundary_reports_right_side_pressure_margin():
    bridge, system, state = fixture(rates=(0.5, 1.0), segment_duration=1.0, resistance=1e5)
    result = run_case(bridge, system, state, duration=1.0)
    assert result.cutoff_id is None and result.cutoff_state is None
    expected_tank = 2e6 * (10.0 / 10.5) ** 1.4
    assert dict(result.final_feed_margins)["feed:injector-pressure"] == pytest.approx(
        expected_tank - 1e5 - 1.5e6, abs=2e-6)


@pytest.mark.parametrize(("reserve", "expected_time"), [(0.0, 10.0), (2.0, 8.0)])
def test_reserve_root_with_solver_overshoot_cannot_publish_negative_or_consumed_reserve(reserve, expected_time):
    bridge, system, state = fixture(minimum_pressure=1e4, reserve=reserve)
    result = run_case(bridge, system, state, duration=12.0, max_step_s=0.93)
    assert result.cutoff_id == "feed:propellant-reserve"
    assert result.cutoff_state.time_s == pytest.approx(expected_time, abs=3e-10)
    assert result.cutoff_state.feed.propellant_mass_kg == reserve
    assert result.final_state.feed.propellant_mass_kg == reserve
    assert result.final_state.vehicle_mass_kg == 90.0 + reserve
    _, speed, position = constant_flow_oracle(expected_time, final=12.0)
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=3e-9)
    assert result.final_state.position_m[0] == pytest.approx(position, abs=2e-8)


def test_pressure_map_lower_edge_stops_before_injector_protection_and_allows_bracketing():
    bridge, system, state = fixture(minimum_pressure=1e5, pressure_grid=(1.6e6, 3e6))
    stop = 10.0 * ((2e6 / (1.6e6 + 1000.0)) ** (1 / 1.4) - 1)
    result = run_case(bridge, system, state)
    assert result.cutoff_id == "performance:pressure-minimum"
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=2e-9)
    _, speed, position = constant_flow_oracle(stop)
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=2e-9)
    assert result.final_state.position_m[0] == pytest.approx(position, abs=2e-8)


def test_dry_mass_stops_before_feed_reserve_and_coast_retains_exact_inventory():
    bridge, system, state = fixture(minimum_pressure=1e4, reserve=1.0, dry_mass=97.0)
    result = run_case(bridge, system, state)
    assert result.cutoff_id == "burn:dry-mass"
    assert result.cutoff_state.time_s == pytest.approx(3.0, abs=2e-10)
    assert result.final_state.vehicle_mass_kg == 97.0
    assert result.final_state.feed.propellant_mass_kg == 7.0
    _, speed, position = constant_flow_oracle(3.0)
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=2e-9)
    assert result.final_state.position_m[0] == pytest.approx(position, abs=2e-8)


def test_nominal_stop_uses_left_limit_and_burn_cannot_restart_during_coast():
    bridge, system, state = fixture(segment_duration=0.5)
    result = run_case(bridge, system, state, duration=2.0)
    assert result.cutoff_id == "burn:stop"
    assert result.cutoff_state.time_s == 0.5
    _, speed, position = constant_flow_oracle(0.5, final=2.0)
    assert result.final_state.feed.propellant_mass_kg == pytest.approx(9.5, abs=2e-10)
    assert result.final_state.velocity_m_s[0] == pytest.approx(speed, abs=2e-9)
    assert result.final_state.position_m[0] == pytest.approx(position, abs=2e-9)
    short = run_case(bridge, system, state, duration=0.25)
    assert short.cutoff_id is None and short.cutoff_state is None
    assert short.final_state.feed.propellant_mass_kg == pytest.approx(9.75, abs=2e-10)


@pytest.mark.parametrize("reason", ["reserve", "empty", "bottle", "tank", "dry"])
def test_initial_unsafe_inventory_and_gas_cannot_request_engine_flow(reason, monkeypatch):
    bridge, system, state = fixture()
    if reason == "reserve":
        system = RegulatedFeedSystem(replace(system.configuration, minimum_propellant_mass_kg=10.0))
        expected = "feed:propellant-reserve"
    elif reason == "empty":
        state = replace(state, propellant_mass_kg=0.0)
        expected = "feed:propellant-reserve"
    elif reason == "bottle":
        state = replace(state, bottle_pressurant_mass_kg=state.bottle_pressurant_mass_kg / 20)
        expected = "feed:bottle-pressure-margin"
    elif reason == "tank":
        state = replace(state, tank_gas_temperature_k=400.0)
        expected = "feed:tank-overpressure"
    else:
        bridge = replace(bridge, protected_dry_mass_kg=100.0)
        expected = "burn:dry-mass"

    def forbidden(*args, **kwargs):
        raise AssertionError("An initial protective stop must precede any engine request.")

    monkeypatch.setattr("rocket_propulsion.integrations.regulated_feed_trajectory."
                        "evaluate_coupled_feed_propulsion", forbidden)
    result = run_case(bridge, system, state, duration=0.1)
    assert result.cutoff_id == expected
    assert result.cutoff_state.time_s == 0.0
    assert result.final_state.velocity_m_s == (0.0, 0.0, 0.0)
    assert result.final_state.feed.propellant_mass_kg == state.propellant_mass_kg


def test_initial_injector_underpressure_and_outside_map_cannot_burn():
    bridge, system, state = fixture(minimum_pressure=1_999_500.0)
    result = run_case(bridge, system, state, duration=0.1)
    assert result.cutoff_id == "feed:injector-pressure" and result.cutoff_state.time_s == 0.0
    assert result.final_state.feed.propellant_mass_kg == 10.0
    bridge, system, state = fixture(pressure_grid=(1e5, 1.8e6))
    result = run_case(bridge, system, state, duration=0.1)
    assert result.cutoff_id == "performance:pressure-maximum" and result.cutoff_state.time_s == 0.0
    assert result.final_state.velocity_m_s == (0.0, 0.0, 0.0)


def test_small_nonzero_flow_cannot_accelerate_without_draining_the_tank():
    bridge, system, state = fixture(minimum_pressure=1e4, rates=(1e-13,), force=1e-11,
                                    segment_duration=25.0)
    state = replace(state, propellant_mass_kg=1e-12)
    surface = replace(bridge.performance_surface, points=tuple(
        PerformanceScalePoint(point.coordinates, 0.5, 0.5)
        for point in bridge.performance_surface.points))
    bridge = replace(bridge, performance_surface=surface)
    result = run_case(bridge, system, state, duration=22.0, initial_vehicle_mass_kg=1e-9,
                      max_step_s=1.0)
    assert result.cutoff_id == "feed:propellant-reserve"
    assert result.cutoff_state.time_s == pytest.approx(20.0, abs=2e-9)
    assert result.final_state.feed.propellant_mass_kg == 0.0
    expected = 100.0 * math.log(1e-9 / (1e-9 - 1e-12))
    assert result.final_state.velocity_m_s[0] == pytest.approx(expected, abs=2e-9)


def test_regulator_keeps_evolving_after_motor_stop_with_closed_pressurant_inventory():
    bridge, system, state = fixture(rates=(0.01,), segment_duration=0.2)
    system = RegulatedFeedSystem(replace(system.configuration, regulator_set_pressure_pa=2.1e6))
    state = replace(state, regulator_opening=0.5)
    result = run_case(bridge, system, state, duration=0.5, max_step_s=0.01)
    assert result.cutoff_id == "burn:stop"
    final = result.final_state.feed
    stopped = result.cutoff_state.feed
    assert final.propellant_mass_kg == pytest.approx(9.998, abs=1e-11)
    assert final.propellant_mass_kg == stopped.propellant_mass_kg
    assert final.tank_pressurant_mass_kg > stopped.tank_pressurant_mass_kg
    assert final.bottle_pressurant_mass_kg < stopped.bottle_pressurant_mass_kg
    assert final.tank_pressurant_mass_kg + final.bottle_pressurant_mass_kg == pytest.approx(
        state.tank_pressurant_mass_kg + state.bottle_pressurant_mass_kg, abs=2e-14)
    expected_speed = 50 / 0.01 * math.log(100 / 99.998)
    assert result.final_state.velocity_m_s[0] == pytest.approx(expected_speed, abs=1e-10)


def test_coast_retains_central_gravity_against_exact_circular_orbit():
    bridge, system, state = fixture(rates=(0.0,), force=0.0, segment_duration=0.1)
    mu, radius, elapsed = 3.986004418e14, 7e6, 600.0
    speed = math.sqrt(mu / radius)
    result = run_case(bridge, system, state, duration=elapsed, mu_m3_s2=mu,
                      initial_position_m=(radius, 0.0, 0.0),
                      initial_velocity_m_s=(0.0, speed, 0.0), max_step_s=20.0)
    angle = math.sqrt(mu / radius**3) * elapsed
    assert result.cutoff_id == "burn:stop"
    assert result.final_state.position_m == pytest.approx(
        (radius * math.cos(angle), radius * math.sin(angle), 0.0), abs=2e-5)
    assert result.final_state.velocity_m_s == pytest.approx(
        (-speed * math.sin(angle), speed * math.cos(angle), 0.0), abs=2e-8)


def test_trial_pressure_extension_is_flagged_has_constant_map_derivatives_and_preserves_default_refusal():
    bridge, system, state = fixture(pressure_grid=(1e5, 1.8e6), pressure_dependent=True)
    args = {"relative_time_s": 0.0, "vehicle_mass_kg": 100.0, "direction": (1.0, 0.0, 0.0),
            "pressure_condition_name": "injector_pa"}
    with pytest.raises(DomainError, match="extrapolation refused"):
        evaluate_coupled_feed_propulsion(bridge, system, state, **args)
    evaluation = evaluate_coupled_feed_propulsion(bridge, system, state, **args,
                                                for_pressure_event_location=True)
    assert evaluation.evidence.pressure_boundary_extension
    assert evaluation.feed.injector_pressure_pa > 1.8e6
    assert evaluation.propulsion.axial_thrust_n == pytest.approx(45.0)
    assert all(value == 0.0 for _, value in evaluation.partials.propellant_flow_wrt_feed_state)
    assert all(value == (0.0, 0.0, 0.0)
               for _, value in evaluation.partials.acceleration_wrt_feed_state)
    for flag in (1, "true", None):
        with pytest.raises(DomainError, match="boolean"):
            evaluate_coupled_feed_propulsion(bridge, system, state, **args,
                                            for_pressure_event_location=flag)


@pytest.mark.parametrize(("expected_id", "threshold"), [
    ("feed:tank-overpressure", 2.3e6),
    ("feed:bottle-pressure-margin", 2.5e6),
    ("performance:pressure-maximum", 2.2e6),
])
def test_heated_sealed_tank_crosses_the_earliest_upper_pressure_bound(expected_id, threshold):
    pressure_grid = (1e5, threshold if expected_id == "performance:pressure-maximum" else 3e6)
    bridge, system, state = fixture(rates=(0.0,), force=0.0, pressure_grid=pressure_grid)
    config = replace(system.configuration, tank_heat_transfer_w_k=10.0,
                     tank_wall_temperature_k=600.0, minimum_bottle_pressure_margin_pa=0.5e6,
                     maximum_tank_pressure_pa=2.3e6 if expected_id == "feed:tank-overpressure"
                     else 3.5e6)
    system = RegulatedFeedSystem(config)
    if expected_id == "feed:bottle-pressure-margin":
        state = replace(state, bottle_pressurant_mass_kg=0.15 * state.bottle_pressurant_mass_kg)
    # Sealed fixed-volume tank. Newton cooling/heating makes p approach 4 MPa
    # with time constant m*cv/UA. Gas mass cancels via pV=mRT.
    tau = 2e6 * 0.01 / (300.0 * 0.4 * 10.0)
    stop = -tau * math.log((4e6 - threshold) / 2e6)
    result = run_case(bridge, system, state, duration=8.0, max_step_s=0.5)
    assert result.cutoff_id == expected_id
    assert result.cutoff_state.time_s == pytest.approx(stop, abs=2e-9)
    assert result.cutoff_state.feed.tank_gas_temperature_k == pytest.approx(
        300.0 * threshold / 2e6, abs=2e-9)
    assert result.final_state.feed.tank_gas_temperature_k == pytest.approx(
        600.0 - 300.0 * math.exp(-8.0 / tau), abs=2e-8)
    if expected_id.startswith("feed:"):
        assert dict(result.final_feed_margins)[expected_id] < 0.0


def test_unbound_multiple_propellant_streams_cannot_be_merged_into_one_liquid():
    bridge, system, state = fixture()
    bridge = replace(bridge, artifact=example_artifact())
    with pytest.raises(DomainError, match="single propellant stream"):
        propagate_regulated_feed(
            bridge, system, state, duration_s=1.0, initial_vehicle_mass_kg=100.0,
            initial_position_m=(0.0, 0.0, 0.0), initial_velocity_m_s=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0), pressure_condition_name="injector_pa", max_step_s=0.1,
        )


def test_pressure_axis_with_incompatible_units_cannot_silently_trip_or_fire():
    bridge, system, state = fixture()
    surface = replace(bridge.performance_surface,
                      axes=(replace(bridge.performance_surface.axes[0], unit="bar"),))
    bridge = replace(bridge, performance_surface=surface)
    with pytest.raises(DomainError, match="pressure axis in Pa"):
        propagate_regulated_feed(
            bridge, system, state, duration_s=1.0, initial_vehicle_mass_kg=100.0,
            initial_position_m=(0.0, 0.0, 0.0), initial_velocity_m_s=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0), pressure_condition_name="injector_pa", max_step_s=0.1,
        )
