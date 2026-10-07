"""Regulated feed dynamics, conservation, events, and bridge coupling tests."""

from __future__ import annotations

import unittest
from dataclasses import replace
from math import sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    FeedEventKind,
    IdealPressurantGas,
    PropagatorPropulsionBridge,
    RegulatedFeedConfiguration,
    RegulatedFeedState,
    RegulatedFeedSystem,
    evaluate_coupled_feed_propulsion,
)
from tests.test_burn_performance_surface import example_surface
from tests.test_burn_propagator import example_artifact


def example_feed_system() -> RegulatedFeedSystem:
    return RegulatedFeedSystem(
        RegulatedFeedConfiguration(
            gas=IdealPressurantGas(2_077.1, 1.667, "helium"),
            tank_internal_volume_m3=0.02,
            propellant_density_kg_m3=1_000.0,
            bottle_volume_m3=0.01,
            regulator_set_pressure_pa=2_100_000.0,
            regulator_full_open_error_pa=200_000.0,
            regulator_maximum_area_m2=1.0e-7,
            regulator_discharge_coefficient=0.8,
            valve_time_constant_s=0.05,
            feed_line_resistance_pa_s2_kg2=1_000.0,
            tank_heat_transfer_w_k=0.0,
            tank_wall_temperature_k=300.0,
            bottle_heat_transfer_w_k=0.0,
            bottle_wall_temperature_k=300.0,
            minimum_propellant_mass_kg=1.0,
            minimum_injector_pressure_pa=1_500_000.0,
            minimum_bottle_pressure_margin_pa=100_000.0,
            maximum_tank_pressure_pa=2_300_000.0,
        )
    )


class RegulatedFeedSystemTests(unittest.TestCase):
    def setUp(self) -> None:
        self.system = example_feed_system()
        self.state = self.system.initial_state(
            propellant_mass_kg=10.0,
            tank_pressure_pa=2_000_000.0,
            tank_gas_temperature_k=300.0,
            bottle_pressure_pa=20_000_000.0,
            bottle_gas_temperature_k=300.0,
            regulator_opening=0.5,
        )

    def test_initial_state_reproduces_pressures_and_regulator_command(self) -> None:
        result = self.system.evaluate(self.state, 0.05)
        self.assertAlmostEqual(result.tank_ullage_volume_m3, 0.01)
        self.assertAlmostEqual(result.tank_pressure_pa, 2_000_000.0)
        self.assertAlmostEqual(result.bottle_pressure_pa, 20_000_000.0)
        self.assertAlmostEqual(result.regulator_command, 0.5)
        self.assertAlmostEqual(result.derivative.regulator_opening_s, 0.0)
        self.assertTrue(result.regulator_choked)
        self.assertGreater(result.regulator_mass_flow_kg_s, 0.0)
        self.assertAlmostEqual(result.feed_line_pressure_drop_pa, 2.5)
        self.assertAlmostEqual(result.injector_pressure_pa, 1_999_997.5)

    def test_mass_and_open_system_energy_balances_close(self) -> None:
        result = self.system.evaluate(self.state, 0.05)
        derivative = result.derivative
        self.assertAlmostEqual(
            derivative.tank_pressurant_mass_kg_s
            + derivative.bottle_pressurant_mass_kg_s,
            0.0,
        )
        self.assertAlmostEqual(derivative.propellant_mass_kg_s, -0.05)
        gas = self.system.configuration.gas
        cv = gas.specific_heat_cv_j_kg_k
        cp = gas.specific_heat_cp_j_kg_k
        ullage_rate = 0.05 / self.system.configuration.propellant_density_kg_m3
        tank_energy_rate = cv * (
            self.state.tank_pressurant_mass_kg
            * derivative.tank_gas_temperature_k_s
            + self.state.tank_gas_temperature_k
            * derivative.tank_pressurant_mass_kg_s
        )
        tank_rhs = (
            result.regulator_mass_flow_kg_s
            * cp
            * self.state.bottle_gas_temperature_k
            - result.tank_pressure_pa * ullage_rate
        )
        self.assertAlmostEqual(tank_energy_rate, tank_rhs)
        bottle_energy_rate = cv * (
            self.state.bottle_pressurant_mass_kg
            * derivative.bottle_gas_temperature_k_s
            + self.state.bottle_gas_temperature_k
            * derivative.bottle_pressurant_mass_kg_s
        )
        bottle_rhs = (
            -result.regulator_mass_flow_kg_s
            * cp
            * self.state.bottle_gas_temperature_k
        )
        self.assertAlmostEqual(bottle_energy_rate, bottle_rhs)

    def test_injector_pressure_partials_match_central_differences(self) -> None:
        nominal = self.system.evaluate(self.state, 0.05)
        partials = dict(nominal.pressure_partials.injector_pressure_wrt_state)
        fields = (
            "propellant_mass_kg",
            "tank_pressurant_mass_kg",
            "tank_gas_temperature_k",
        )
        steps = (1.0e-5, 1.0e-8, 1.0e-3)
        values = self.state.__dict__ if hasattr(self.state, "__dict__") else {
            name: getattr(self.state, name)
            for name in self.state.__dataclass_fields__
        }
        for field, step in zip(fields, steps, strict=True):
            plus_values = {**values, field: values[field] + step}
            minus_values = {**values, field: values[field] - step}
            plus = self.system.evaluate(RegulatedFeedState(**plus_values), 0.05)
            minus = self.system.evaluate(RegulatedFeedState(**minus_values), 0.05)
            finite = (plus.injector_pressure_pa - minus.injector_pressure_pa) / (2.0 * step)
            self.assertAlmostEqual(partials[field], finite, delta=abs(finite) * 2e-7)
        flow_step = 1.0e-6
        plus = self.system.evaluate(self.state, 0.05 + flow_step)
        minus = self.system.evaluate(self.state, 0.05 - flow_step)
        finite_flow = (plus.injector_pressure_pa - minus.injector_pressure_pa) / (
            2.0 * flow_step
        )
        self.assertAlmostEqual(
            nominal.pressure_partials.injector_pressure_wrt_propellant_flow,
            finite_flow,
            delta=1e-3,
        )

    def test_cutoff_surfaces_are_positive_safe_and_reach_zero(self) -> None:
        surfaces = self.system.event_surfaces()
        self.assertEqual(
            [surface.kind for surface in surfaces],
            [
                FeedEventKind.PROPELLANT_RESERVE,
                FeedEventKind.INJECTOR_PRESSURE,
                FeedEventKind.BOTTLE_PRESSURE_MARGIN,
                FeedEventKind.TANK_OVERPRESSURE,
            ],
        )
        self.assertTrue(
            all(self.system.event_value(surface, self.state, 0.05) > 0.0 for surface in surfaces)
        )
        reserve_state = RegulatedFeedState(
            propellant_mass_kg=1.0,
            tank_pressurant_mass_kg=self.state.tank_pressurant_mass_kg,
            tank_gas_temperature_k=self.state.tank_gas_temperature_k,
            bottle_pressurant_mass_kg=self.state.bottle_pressurant_mass_kg,
            bottle_gas_temperature_k=self.state.bottle_gas_temperature_k,
            regulator_opening=self.state.regulator_opening,
        )
        self.assertEqual(self.system.event_value(surfaces[0], reserve_state, 0.0), 0.0)
        cutoff_flow = sqrt(
            (2_000_000.0 - 1_500_000.0)
            / self.system.configuration.feed_line_resistance_pa_s2_kg2
        )
        self.assertAlmostEqual(
            self.system.event_value(surfaces[1], self.state, cutoff_flow), 0.0
        )

    def test_model_manifest_binds_configuration(self) -> None:
        manifest = self.system.manifest()
        self.assertEqual(manifest["model_hash"], self.system.model_hash)
        self.assertEqual(manifest["schema"], "rocket_propulsion_regulated_feed_system_v1")
        self.assertIn("NASA-20240003493", manifest["reference_ids"])


class CoupledFeedPerformanceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.feed_system = example_feed_system()
        self.feed_state = self.feed_system.initial_state(
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

    def test_feed_pressure_and_engine_flow_close_iteratively(self) -> None:
        result = evaluate_coupled_feed_propulsion(
            self.bridge,
            self.feed_system,
            self.feed_state,
            relative_time_s=0.5,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            pressure_condition_name="supply_pressure_pa",
            other_operating_conditions={"ambient_pressure_pa": 50_000.0},
            initial_mass_flow_kg_s=3.0,
            relaxation=0.8,
        )
        self.assertTrue(result.evidence.converged)
        self.assertLess(result.evidence.mass_flow_residual_kg_s, 1e-9)
        self.assertAlmostEqual(
            -result.propulsion.mass_derivative_kg_s,
            -result.feed.derivative.propellant_mass_kg_s,
            places=9,
        )
        conditions = dict(result.propulsion.operating_conditions)
        self.assertAlmostEqual(
            conditions["supply_pressure_pa"], result.feed.injector_pressure_pa
        )
        self.assertGreater(result.partials.algebraic_denominator, 1.0)

    def test_implicit_coupled_state_partials_match_resolved_finite_difference(self) -> None:
        nominal = evaluate_coupled_feed_propulsion(
            self.bridge,
            self.feed_system,
            self.feed_state,
            relative_time_s=0.5,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            pressure_condition_name="supply_pressure_pa",
            other_operating_conditions={"ambient_pressure_pa": 50_000.0},
            initial_mass_flow_kg_s=3.0,
            relative_tolerance=1e-12,
            relaxation=0.8,
        )
        step = 1.0e-8
        plus_state = replace(
            self.feed_state,
            tank_pressurant_mass_kg=self.feed_state.tank_pressurant_mass_kg + step,
        )
        minus_state = replace(
            self.feed_state,
            tank_pressurant_mass_kg=self.feed_state.tank_pressurant_mass_kg - step,
        )

        def solve(state: RegulatedFeedState):
            return evaluate_coupled_feed_propulsion(
                self.bridge,
                self.feed_system,
                state,
                relative_time_s=0.5,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
                pressure_condition_name="supply_pressure_pa",
                other_operating_conditions={"ambient_pressure_pa": 50_000.0},
                initial_mass_flow_kg_s=-nominal.propulsion.mass_derivative_kg_s,
                relative_tolerance=1e-12,
                relaxation=0.8,
            )

        plus = solve(plus_state)
        minus = solve(minus_state)
        finite_pressure = (
            plus.feed.injector_pressure_pa - minus.feed.injector_pressure_pa
        ) / (2.0 * step)
        finite_flow = (
            -plus.propulsion.mass_derivative_kg_s
            + minus.propulsion.mass_derivative_kg_s
        ) / (2.0 * step)
        finite_acceleration = (
            plus.propulsion.acceleration_m_s2[0]
            - minus.propulsion.acceleration_m_s2[0]
        ) / (2.0 * step)
        pressure_partials = dict(nominal.partials.injector_pressure_wrt_feed_state)
        flow_partials = dict(nominal.partials.propellant_flow_wrt_feed_state)
        acceleration_partials = dict(nominal.partials.acceleration_wrt_feed_state)
        name = "tank_pressurant_mass_kg"
        self.assertAlmostEqual(
            pressure_partials[name], finite_pressure, delta=abs(finite_pressure) * 2e-6
        )
        self.assertAlmostEqual(
            flow_partials[name], finite_flow, delta=abs(finite_flow) * 2e-6
        )
        self.assertAlmostEqual(
            acceleration_partials[name][0],
            finite_acceleration,
            delta=abs(finite_acceleration) * 2e-6,
        )

    def test_coupling_configuration_errors_and_nonconvergence_fail_closed(self) -> None:
        bridge_without_surface = PropagatorPropulsionBridge(
            example_artifact(), require_complete_stream_binding=False
        )
        with self.assertRaisesRegex(DomainError, "requires a bridge performance surface"):
            evaluate_coupled_feed_propulsion(
                bridge_without_surface,
                self.feed_system,
                self.feed_state,
                relative_time_s=0.5,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
                pressure_condition_name="supply_pressure_pa",
            )
        with self.assertRaisesRegex(DomainError, "did not converge"):
            evaluate_coupled_feed_propulsion(
                self.bridge,
                self.feed_system,
                self.feed_state,
                relative_time_s=0.5,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
                pressure_condition_name="supply_pressure_pa",
                other_operating_conditions={"ambient_pressure_pa": 50_000.0},
                maximum_iterations=1,
            )


if __name__ == "__main__":
    unittest.main()
