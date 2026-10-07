"""Dynamic liquid feed-line physics, events, and propagator coupling tests."""

from __future__ import annotations

import unittest
from dataclasses import replace

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    DynamicFeedLine,
    DynamicFeedLineConfiguration,
    DynamicFeedLineEventKind,
    DynamicFeedLineState,
    PropagatorPropulsionBridge,
    RegulatedFeedSystem,
    evaluate_dynamic_feed_propulsion,
)
from tests.test_burn_feed_system import example_feed_system
from tests.test_burn_performance_surface import example_surface
from tests.test_burn_propagator import example_artifact


def example_line(*, resistance: float = 4_000_000.0) -> DynamicFeedLine:
    """Return a deterministic two-state line used by the analytic checks."""

    return DynamicFeedLine(
        DynamicFeedLineConfiguration(
            inertance_pa_s2_kg=1_000_000.0,
            compliance_kg_pa=1.0e-8,
            resistance_pa_s2_kg2=resistance,
            minimum_manifold_pressure_pa=1_000_000.0,
            maximum_manifold_pressure_pa=3_000_000.0,
            maximum_absolute_flow_kg_s=5.0,
        )
    )


class DynamicFeedLineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.line = example_line()
        self.state = DynamicFeedLineState(0.2, 1_840_000.0)

    def evaluate(self, state: DynamicFeedLineState):
        return self.line.evaluate(
            state,
            upstream_pressure_pa=2_000_000.0,
            engine_mass_flow_kg_s=0.2,
        )

    def test_steady_state_modal_diagnostics_and_mass_closure(self) -> None:
        value = self.evaluate(self.state)
        self.assertAlmostEqual(value.friction_pressure_drop_pa, 160_000.0)
        self.assertAlmostEqual(value.derivative.mass_flow_kg_s2, 0.0)
        self.assertAlmostEqual(value.derivative.manifold_pressure_pa_s, 0.0)
        self.assertAlmostEqual(value.stored_liquid_mass_rate_kg_s, 0.0)
        self.assertEqual(value.external_mass_closure_error_kg_s, 0.0)
        self.assertAlmostEqual(value.natural_angular_frequency_rad_s, 10.0)
        self.assertAlmostEqual(value.damping_ratio, 0.08)

    def test_local_partials_match_central_finite_differences(self) -> None:
        nominal = self.evaluate(self.state)
        state_partials = dict(nominal.partials.derivative_wrt_state)
        cases = (
            ("mass_flow_kg_s", 1.0e-6),
            ("manifold_pressure_pa", 1.0),
        )
        for name, step in cases:
            plus = self.evaluate(replace(self.state, **{name: getattr(self.state, name) + step}))
            minus = self.evaluate(replace(self.state, **{name: getattr(self.state, name) - step}))
            finite = tuple(
                (right - left) / (2.0 * step)
                for right, left in zip(
                    plus.derivative.as_tuple(), minus.derivative.as_tuple(), strict=True
                )
            )
            for analytic, numeric in zip(state_partials[name], finite, strict=True):
                self.assertAlmostEqual(analytic, numeric, delta=max(1e-9, abs(numeric) * 1e-8))

        pressure_step = 1.0
        plus_pressure = self.line.evaluate(
            self.state,
            upstream_pressure_pa=2_000_000.0 + pressure_step,
            engine_mass_flow_kg_s=0.2,
        )
        minus_pressure = self.line.evaluate(
            self.state,
            upstream_pressure_pa=2_000_000.0 - pressure_step,
            engine_mass_flow_kg_s=0.2,
        )
        finite_pressure = tuple(
            (right - left) / (2.0 * pressure_step)
            for right, left in zip(
                plus_pressure.derivative.as_tuple(),
                minus_pressure.derivative.as_tuple(),
                strict=True,
            )
        )
        self.assertEqual(
            nominal.partials.derivative_wrt_upstream_pressure,
            finite_pressure,
        )

        flow_step = 1.0e-7
        plus_demand = self.line.evaluate(
            self.state,
            upstream_pressure_pa=2_000_000.0,
            engine_mass_flow_kg_s=0.2 + flow_step,
        )
        minus_demand = self.line.evaluate(
            self.state,
            upstream_pressure_pa=2_000_000.0,
            engine_mass_flow_kg_s=0.2 - flow_step,
        )
        finite_demand = tuple(
            (right - left) / (2.0 * flow_step)
            for right, left in zip(
                plus_demand.derivative.as_tuple(),
                minus_demand.derivative.as_tuple(),
                strict=True,
            )
        )
        for analytic, numeric in zip(
            nominal.partials.derivative_wrt_engine_flow,
            finite_demand,
            strict=True,
        ):
            self.assertAlmostEqual(analytic, numeric, delta=max(1e-8, abs(numeric) * 1e-9))

    def test_cutoff_surfaces_are_positive_safe_and_reach_zero(self) -> None:
        surfaces = self.line.event_surfaces()
        self.assertEqual(
            tuple(surface.kind for surface in surfaces),
            (
                DynamicFeedLineEventKind.REVERSE_FLOW,
                DynamicFeedLineEventKind.MINIMUM_MANIFOLD_PRESSURE,
                DynamicFeedLineEventKind.MAXIMUM_MANIFOLD_PRESSURE,
                DynamicFeedLineEventKind.MAXIMUM_ABSOLUTE_FLOW,
            ),
        )
        self.assertTrue(
            all(self.line.event_value(surface, self.state) > 0.0 for surface in surfaces)
        )
        boundary_states = (
            DynamicFeedLineState(0.0, self.state.manifold_pressure_pa),
            DynamicFeedLineState(self.state.mass_flow_kg_s, 1_000_000.0),
            DynamicFeedLineState(self.state.mass_flow_kg_s, 3_000_000.0),
            DynamicFeedLineState(5.0, self.state.manifold_pressure_pa),
        )
        for surface, boundary in zip(surfaces, boundary_states, strict=True):
            self.assertEqual(self.line.event_value(surface, boundary), 0.0)

    def test_manifest_binds_configuration_and_references(self) -> None:
        manifest = self.line.manifest()
        self.assertEqual(manifest["model_hash"], self.line.model_hash)
        self.assertEqual(manifest["schema"], "rocket_propulsion_dynamic_feed_line_v1")
        self.assertIn("NASA-19740028545", manifest["reference_ids"])
        changed = example_line(resistance=4_000_001.0)
        self.assertNotEqual(changed.model_hash, self.line.model_hash)


class DynamicFeedPropulsionCouplingTests(unittest.TestCase):
    def setUp(self) -> None:
        static_feed = example_feed_system()
        self.feed = RegulatedFeedSystem(
            replace(
                static_feed.configuration,
                feed_line_resistance_pa_s2_kg2=0.0,
            )
        )
        self.feed_state = self.feed.initial_state(
            propellant_mass_kg=10.0,
            tank_pressure_pa=2_000_000.0,
            tank_gas_temperature_k=300.0,
            bottle_pressure_pa=20_000_000.0,
            bottle_gas_temperature_k=300.0,
            regulator_opening=0.5,
        )
        self.bridge = PropagatorPropulsionBridge(
            example_artifact(),
            performance_surface=example_surface(),
            require_complete_stream_binding=False,
        )
        self.engine_flow = 3.255
        resistance = (2_000_000.0 - 1_800_000.0) / self.engine_flow**2
        self.line = example_line(resistance=resistance)
        self.line_state = DynamicFeedLineState(self.engine_flow, 1_800_000.0)

    def evaluate(
        self,
        line_state: DynamicFeedLineState | None = None,
        *,
        feed_state=None,
        relative_time_s: float = 0.5,
    ):
        return evaluate_dynamic_feed_propulsion(
            self.bridge,
            self.feed,
            self.feed_state if feed_state is None else feed_state,
            self.line,
            self.line_state if line_state is None else line_state,
            relative_time_s=relative_time_s,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            pressure_condition_name="supply_pressure_pa",
            other_operating_conditions={"ambient_pressure_pa": 50_000.0},
        )

    def test_dynamic_line_closes_engine_demand_and_vehicle_mass(self) -> None:
        value = self.evaluate()
        self.assertTrue(value.propulsion.active)
        self.assertAlmostEqual(-value.propulsion.mass_derivative_kg_s, self.engine_flow)
        self.assertAlmostEqual(
            value.regulated_feed.derivative.propellant_mass_kg_s,
            -self.line_state.mass_flow_kg_s,
        )
        self.assertAlmostEqual(value.feed_line.stored_liquid_mass_rate_kg_s, 0.0)
        self.assertAlmostEqual(value.feed_line.derivative.mass_flow_kg_s2, 0.0)
        self.assertAlmostEqual(value.feed_line.derivative.manifold_pressure_pa_s, 0.0)
        self.assertEqual(value.feed_line.external_mass_closure_error_kg_s, 0.0)
        self.assertEqual(
            dict(value.propulsion.operating_conditions)["supply_pressure_pa"],
            self.line_state.manifold_pressure_pa,
        )

    def test_composed_line_and_acceleration_jacobians_match_finite_differences(self) -> None:
        nominal = self.evaluate()
        jacobian = nominal.partials.line_state_jacobian
        acceleration_partials = dict(nominal.partials.acceleration_wrt_line_state)
        for column, (name, step) in enumerate(
            (("mass_flow_kg_s", 1.0e-5), ("manifold_pressure_pa", 1.0))
        ):
            plus = self.evaluate(
                replace(self.line_state, **{name: getattr(self.line_state, name) + step})
            )
            minus = self.evaluate(
                replace(self.line_state, **{name: getattr(self.line_state, name) - step})
            )
            finite_rhs = tuple(
                (right - left) / (2.0 * step)
                for right, left in zip(
                    plus.feed_line.derivative.as_tuple(),
                    minus.feed_line.derivative.as_tuple(),
                    strict=True,
                )
            )
            self.assertAlmostEqual(jacobian[0][column], finite_rhs[0], delta=1e-9)
            self.assertAlmostEqual(jacobian[1][column], finite_rhs[1], delta=1e-3)
            finite_acceleration = tuple(
                (right - left) / (2.0 * step)
                for right, left in zip(
                    plus.propulsion.acceleration_m_s2,
                    minus.propulsion.acceleration_m_s2,
                    strict=True,
                )
            )
            for analytic, numeric in zip(
                acceleration_partials[name], finite_acceleration, strict=True
            ):
                self.assertAlmostEqual(analytic, numeric, delta=1e-11)

    def test_feed_state_to_line_partial_matches_resolved_finite_difference(self) -> None:
        nominal = self.evaluate()
        partials = dict(nominal.partials.line_derivative_wrt_feed_state)
        step = 1.0e-8
        plus = self.evaluate(
            feed_state=replace(
                self.feed_state,
                tank_pressurant_mass_kg=self.feed_state.tank_pressurant_mass_kg + step,
            )
        )
        minus = self.evaluate(
            feed_state=replace(
                self.feed_state,
                tank_pressurant_mass_kg=self.feed_state.tank_pressurant_mass_kg - step,
            )
        )
        finite = tuple(
            (right - left) / (2.0 * step)
            for right, left in zip(
                plus.feed_line.derivative.as_tuple(),
                minus.feed_line.derivative.as_tuple(),
                strict=True,
            )
        )
        for analytic, numeric in zip(partials["tank_pressurant_mass_kg"], finite, strict=True):
            self.assertAlmostEqual(analytic, numeric, delta=max(1e-7, abs(numeric) * 2e-8))

    def test_inactive_engine_has_zero_pressure_gain_while_line_relaxes(self) -> None:
        value = self.evaluate(relative_time_s=3.0)
        self.assertFalse(value.propulsion.active)
        self.assertEqual(value.propulsion.mass_derivative_kg_s, 0.0)
        self.assertEqual(value.partials.engine_flow_wrt_manifold_pressure, 0.0)
        self.assertEqual(
            dict(value.partials.acceleration_wrt_line_state)["manifold_pressure_pa"],
            (0.0, 0.0, 0.0),
        )
        self.assertGreater(value.feed_line.derivative.manifold_pressure_pa_s, 0.0)

    def test_static_and_dynamic_line_resistance_cannot_be_double_counted(self) -> None:
        with self.assertRaisesRegex(DomainError, "static line resistance"):
            evaluate_dynamic_feed_propulsion(
                self.bridge,
                example_feed_system(),
                self.feed_state,
                self.line,
                self.line_state,
                relative_time_s=0.5,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
                pressure_condition_name="supply_pressure_pa",
                other_operating_conditions={"ambient_pressure_pa": 50_000.0},
            )


if __name__ == "__main__":
    unittest.main()
