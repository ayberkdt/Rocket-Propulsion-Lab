"""Adversarial controls for the public lumped liquid-line model."""

import math
from dataclasses import replace
from decimal import Context, Decimal, localcontext

import pytest
from scipy.integrate import solve_ivp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    DynamicFeedLine,
    DynamicFeedLineEventKind,
    DynamicFeedLineState,
)
from tests.test_burn_feed_line import example_line


@pytest.mark.parametrize(
    ("flow", "upstream", "manifold"),
    [(0.1, 2.0, 1.0), (-0.1, 1.0, 2.0)],
)
def test_signed_acceleration_survives_head_loss_cancellation(flow, upstream, manifold):
    line = example_line(resistance=100.0)
    with localcontext(Context(prec=1200)):
        q = Decimal(flow)
        expected = (
            Decimal(upstream) - Decimal(manifold) - Decimal(100) * q * abs(q)
        ) / Decimal(line.configuration.inertance_pa_s2_kg)
    actual = line.evaluate(
        DynamicFeedLineState(flow, manifold),
        upstream_pressure_pa=upstream,
        engine_mass_flow_kg_s=0.0,
    )
    assert expected != 0
    assert actual.derivative.mass_flow_kg_s2 == float(expected)


def test_finite_jacobian_when_doubling_resistance_would_overflow():
    line = DynamicFeedLine(replace(
        example_line(resistance=1e308).configuration,
        inertance_pa_s2_kg=1e308,
        compliance_kg_pa=1.0,
    ))
    actual = line.evaluate(
        DynamicFeedLineState(0.5, 2.0), upstream_pressure_pa=3.0, engine_mass_flow_kg_s=0.0
    )
    assert actual.friction_pressure_drop_pa == 2.5e307
    assert actual.derivative.mass_flow_kg_s2 == -0.25
    assert dict(actual.partials.derivative_wrt_state)["mass_flow_kg_s"][0] == -1.0
    assert actual.damping_ratio == pytest.approx(5e153, rel=3e-16)


@pytest.mark.parametrize(
    ("inertance", "compliance", "resistance", "flow", "field"),
    [
        (1e308, 1.0, 100.0, 0.1, "Feed-line flow derivative"),
        (1.0, 2.0, 0.0, math.ulp(0.0), "Feed-line pressure derivative"),
        (1.0, 1.0, 1.0, math.ulp(0.0), "Feed-line friction pressure drop"),
    ],
)
def test_nonzero_physical_output_cannot_be_exported_as_zero(
    inertance, compliance, resistance, flow, field
):
    line = DynamicFeedLine(replace(
        example_line(resistance=resistance).configuration,
        inertance_pa_s2_kg=inertance,
        compliance_kg_pa=compliance,
    ))
    with pytest.raises(DomainError) as error:
        line.evaluate(
            DynamicFeedLineState(flow, 1.0),
            upstream_pressure_pa=2.0,
            engine_mass_flow_kg_s=0.0,
        )
    assert error.value.field == field


def test_mass_residual_audits_the_exported_pressure_rate():
    line = DynamicFeedLine(replace(
        example_line(resistance=0.0).configuration, compliance_kg_pa=0.3
    ))
    actual = line.evaluate(
        DynamicFeedLineState(0.1, 1.0), upstream_pressure_pa=1.0, engine_mass_flow_kg_s=0.0
    )
    with localcontext(Context(prec=1200)):
        residual = abs(
            Decimal.from_float(0.3) * Decimal(actual.derivative.manifold_pressure_pa_s)
            - Decimal.from_float(0.1)
        )
    assert residual > 0
    assert actual.external_mass_closure_error_kg_s == float(residual)


def test_demand_step_reaches_pressure_cutoff_at_closed_form_time():
    line = example_line(resistance=0.0)
    p_up = 2e6
    demand = 0.2
    cutoff = line.configuration.minimum_manifold_pressure_pa
    inertance = line.configuration.inertance_pa_s2_kg
    compliance = line.configuration.compliance_kg_pa
    omega = 1 / math.sqrt(inertance * compliance)
    amplitude = demand / (compliance * omega)
    angle = math.asin((p_up - cutoff) / amplitude)
    pressure_event = next(
        event for event in line.event_surfaces()
        if event.kind is DynamicFeedLineEventKind.MINIMUM_MANIFOLD_PRESSURE
    )

    def rhs(_time, values):
        return line.evaluate(
            DynamicFeedLineState(*map(float, values)),
            upstream_pressure_pa=p_up,
            engine_mass_flow_kg_s=demand,
        ).derivative.as_tuple()

    def reaches_cutoff(_time, values):
        return line.event_value(pressure_event, DynamicFeedLineState(*map(float, values)))

    reaches_cutoff.terminal = True
    reaches_cutoff.direction = -1
    # A missed event must finish inside the positive-pressure numerical domain.
    trace = solve_ivp(
        rhs, (0, 1.5 * angle / omega), (0, p_up), events=reaches_cutoff,
        method="DOP853", rtol=1e-11, atol=(1e-13, 1e-6), max_step=0.005,
    )
    assert trace.status == 1
    assert trace.t_events[0] == pytest.approx([angle / omega], abs=1e-11)
    flow, pressure = trace.y_events[0][0]
    assert flow == pytest.approx(demand * (1 - math.cos(angle)), abs=1e-11)
    assert pressure == pytest.approx(cutoff, abs=1e-5)
