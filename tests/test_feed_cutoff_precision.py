"""Adversarial tests for public feed cutoff surfaces and root location."""

import math
from dataclasses import astuple, replace
from decimal import Context, Decimal, localcontext

import pytest
from scipy.integrate import solve_ivp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    DynamicFeedLineEventKind,
    DynamicFeedLineEventSurface,
    DynamicFeedLineState,
    FeedEventKind,
    FeedSystemEventSurface,
    IdealPressurantGas,
    RegulatedFeedState,
    RegulatedFeedSystem,
)
from tests.test_burn_feed_line import example_line
from tests.test_burn_feed_system import example_feed_system


def binary_fixture():
    """Use an exactly unit-sized ullage and gas constant to isolate rounding."""
    system = RegulatedFeedSystem(replace(
        example_feed_system().configuration,
        gas=IdealPressurantGas(1.0, 1.4, "synthetic"),
        tank_internal_volume_m3=2.0,
        propellant_density_kg_m3=1.0,
        bottle_volume_m3=1.0,
        regulator_set_pressure_pa=1.0,
        feed_line_resistance_pa_s2_kg2=math.ldexp(1.0, -54),
        minimum_propellant_mass_kg=0.5,
        minimum_injector_pressure_pa=2.0,
        maximum_tank_pressure_pa=2.0,
        minimum_bottle_pressure_margin_pa=math.ulp(2.0),
    ))
    return system, RegulatedFeedState(1.0, 1.0, 2.0, 1.0, math.nextafter(2.0, math.inf), 0.0)


def surface_for(system, kind):
    return next(surface for surface in system.event_surfaces() if surface.kind is kind)


def test_zero_reserve_admits_empty_tank_and_zero_demand_but_no_withdrawal():
    nominal, _ = binary_fixture()
    system = RegulatedFeedSystem(replace(nominal.configuration, minimum_propellant_mass_kg=0.0))
    try:
        empty = system.initial_state(
            propellant_mass_kg=0.0, tank_pressure_pa=2.0, tank_gas_temperature_k=2.0,
            bottle_pressure_pa=2.0, bottle_gas_temperature_k=2.0,
        )
    except DomainError as error:
        raise AssertionError("The configured reserve boundary must be constructible.") from error
    reserve = surface_for(system, FeedEventKind.PROPELLANT_RESERVE)
    assert system.event_value(reserve, empty, 0.0) == 0.0
    assert system.event_value(reserve, empty, 1.0) == 0.0
    last_positive = replace(empty, propellant_mass_kg=math.ulp(0.0))
    assert system.event_value(reserve, last_positive, 0.0) == math.ulp(0.0)
    result = system.evaluate(empty, 0.0)
    assert result.tank_ullage_volume_m3 == system.configuration.tank_internal_volume_m3
    assert result.tank_pressure_pa == 2.0
    assert result.derivative.propellant_mass_kg_s == 0.0
    try:
        system.evaluate(empty, math.ulp(0.0))
    except DomainError as error:
        assert error.field == "Propellant mass flow"
    else:
        raise AssertionError("Positive withdrawal from an empty tank was admitted.")
    with pytest.raises(DomainError, match="empty tank"):
        system.evaluate(empty, -math.ulp(0.0), for_event_location=True)
    with pytest.raises(DomainError):
        replace(empty, propellant_mass_kg=-math.ulp(0.0))


@pytest.mark.parametrize(
    "kind",
    [FeedEventKind.INJECTOR_PRESSURE, FeedEventKind.BOTTLE_PRESSURE_MARGIN,
     FeedEventKind.TANK_OVERPRESSURE],
)
def test_sub_display_pressure_margin_retains_both_signs_and_exact_zero(kind):
    system, boundary = binary_fixture()
    flow = 1.0 if kind is FeedEventKind.INJECTOR_PRESSURE else 0.0
    config = system.configuration
    if kind is FeedEventKind.INJECTOR_PRESSURE:
        states = [boundary, replace(boundary, tank_pressurant_mass_kg=math.nextafter(1.0, math.inf))]
    else:
        unsafe = replace(
            boundary, tank_pressurant_mass_kg=math.nextafter(1.0, math.inf),
            tank_gas_temperature_k=math.nextafter(2.0, 0.0),
        )
        states = [unsafe, replace(boundary, tank_pressurant_mass_kg=math.nextafter(1.0, 0.0))]
    surface = surface_for(system, kind)
    assert system.event_value(surface, boundary, 0.0) == 0.0

    with localcontext(Context(prec=1200)):
        for index, state in enumerate(states):
            ullage = Decimal(config.tank_internal_volume_m3) - (
                Decimal(state.propellant_mass_kg) / Decimal(config.propellant_density_kg_m3)
            )
            tank = Decimal(state.tank_pressurant_mass_kg) * Decimal(
                state.tank_gas_temperature_k
            ) / ullage
            if kind is FeedEventKind.INJECTOR_PRESSURE:
                expected = (tank - Decimal(config.feed_line_resistance_pa_s2_kg2)
                            * Decimal(flow)**2 - Decimal(config.minimum_injector_pressure_pa))
            elif kind is FeedEventKind.BOTTLE_PRESSURE_MARGIN:
                expected = (Decimal(state.bottle_pressurant_mass_kg)
                            * Decimal(state.bottle_gas_temperature_k)
                            / Decimal(config.bottle_volume_m3) - tank
                            - Decimal(config.minimum_bottle_pressure_margin_pa))
            else:
                expected = Decimal(config.maximum_tank_pressure_pa) - tank
            assert expected < 0 if index == 0 else expected > 0
            assert system.event_value(surface, state, flow) == float(expected)


def test_reserve_and_overpressure_queries_ignore_unrelated_static_loss_overflow():
    nominal, state = binary_fixture()
    system = RegulatedFeedSystem(replace(
        nominal.configuration, feed_line_resistance_pa_s2_kg2=1e308,
        maximum_tank_pressure_pa=3.0,
    ))
    with pytest.raises(DomainError):
        system.evaluate(state, 2.0)
    reserve = surface_for(system, FeedEventKind.PROPELLANT_RESERVE)
    overpressure = surface_for(system, FeedEventKind.TANK_OVERPRESSURE)
    injector = surface_for(system, FeedEventKind.INJECTOR_PRESSURE)
    try:
        reserve_margin = system.event_value(reserve, state, 2.0)
        overpressure_margin = system.event_value(overpressure, state, 2.0)
    except DomainError as error:
        raise AssertionError("An unrelated RHS failure hid a protective margin.") from error
    assert reserve_margin == 0.5
    assert overpressure_margin == 1.0
    try:
        system.event_value(injector, state, 2.0)
    except DomainError as error:
        assert error.field == "Feed cutoff margin"
    else:
        pytest.fail("An unrepresentable signed injector margin must be refused.")


def test_unconfigured_tank_cutoff_returns_a_domain_failure():
    nominal, state = binary_fixture()
    system = RegulatedFeedSystem(replace(nominal.configuration, maximum_tank_pressure_pa=None))
    overpressure = FeedSystemEventSurface("test:overpressure", FeedEventKind.TANK_OVERPRESSURE)
    with pytest.raises(DomainError, match="requires a configured maximum"):
        system.event_value(overpressure, state, 0.0)


@pytest.mark.parametrize(
    ("record", "kind"),
    [(FeedSystemEventSurface, FeedEventKind.PROPELLANT_RESERVE),
     (DynamicFeedLineEventSurface, DynamicFeedLineEventKind.REVERSE_FLOW)],
)
def test_cutoff_metadata_cannot_forge_a_stop_flag_or_an_event_kind(record, kind):
    assert record("test:cutoff", kind.value, False).kind is kind
    for args in (("", kind, True), (" ", kind, True), ("event", "typo", True),
                 ("event", kind, "false"), ("event", kind, 1)):
        try:
            record(*args)
        except DomainError:
            continue
        raise AssertionError(f"Invalid cutoff metadata accepted: {args!r}")


def test_cutoff_entry_points_reject_wrong_state_and_surface_records():
    system, feed_state = binary_fixture()
    feed_event = surface_for(system, FeedEventKind.PROPELLANT_RESERVE)
    line = example_line()
    line_event = surface_for(line, DynamicFeedLineEventKind.REVERSE_FLOW)
    line_state = DynamicFeedLineState(0.25, 1.8e6)
    for evaluator in (
        lambda: system.event_value(line_event, feed_state, 0.0),
        lambda: system.event_value(feed_event, line_state, 0.0),
        lambda: line.event_value(feed_event, line_state),
        lambda: line.event_value(line_event, feed_state),
        lambda: system.event_value(feed_event, feed_state, -0.1),
    ):
        with pytest.raises(DomainError):
            evaluator()


def test_adiabatic_tank_reaches_injector_cutoff_at_closed_form_time():
    system = RegulatedFeedSystem(replace(
        example_feed_system().configuration, regulator_set_pressure_pa=1e6,
    ))
    p0 = 2e6
    q = 0.05
    config = system.configuration
    initial = system.initial_state(
        propellant_mass_kg=10.0, tank_pressure_pa=p0, tank_gas_temperature_k=300.0,
        bottle_pressure_pa=20e6, bottle_gas_temperature_k=300.0, regulator_opening=0.0,
    )
    cutoff = surface_for(system, FeedEventKind.INJECTOR_PRESSURE)
    ullage0 = config.tank_internal_volume_m3 - 10.0 / config.propellant_density_kg_m3
    cutoff_tank = config.minimum_injector_pressure_pa + config.feed_line_resistance_pa_s2_kg2*q*q
    crossing_time = (
        ullage0 * config.propellant_density_kg_m3 / q
        * ((p0/cutoff_tank)**(1/config.gas.gamma) - 1)
    )

    def rhs(_time, values):
        return system.evaluate(RegulatedFeedState(*map(float, values)), q).derivative.as_tuple()

    def reaches_cutoff(_time, values):
        return system.event_value(cutoff, RegulatedFeedState(*map(float, values)), q)

    reaches_cutoff.terminal = True
    reaches_cutoff.direction = -1
    trace = solve_ivp(
        rhs, (0.0, 1.3 * crossing_time), astuple(initial), events=reaches_cutoff,
        rtol=1e-10, atol=1e-12, max_step=0.3,
    )
    assert trace.status == 1
    assert trace.t_events[0] == pytest.approx([crossing_time], abs=2e-8)
    assert trace.y_events[0][0][0] == pytest.approx(10.0 - q*crossing_time, abs=1e-9)
