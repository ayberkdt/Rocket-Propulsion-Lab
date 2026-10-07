"""State-dependent performance surface and propagator-coupling tests."""

from __future__ import annotations

import json
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
)
from tests.test_burn_propagator import example_artifact


def example_surface() -> RectilinearPerformanceSurface:
    return RectilinearPerformanceSurface(
        axes=(
            PerformanceAxis("supply_pressure_pa", "Pa", (1_000_000.0, 2_000_000.0)),
            PerformanceAxis("ambient_pressure_pa", "Pa", (0.0, 100_000.0)),
        ),
        points=(
            PerformanceScalePoint((1_000_000.0, 0.0), 0.8, 0.9),
            PerformanceScalePoint((1_000_000.0, 100_000.0), 0.7, 0.95),
            PerformanceScalePoint((2_000_000.0, 0.0), 1.0, 1.1),
            PerformanceScalePoint((2_000_000.0, 100_000.0), 0.9, 1.15),
        ),
        source_id="qualification-grid-A",
        source_sha256="b" * 64,
    )


class PerformanceSurfaceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = example_surface()
        self.conditions = {
            "supply_pressure_pa": 1_500_000.0,
            "ambient_pressure_pa": 50_000.0,
        }

    def test_multilinear_value_and_analytic_gradient(self) -> None:
        result = self.surface.evaluate(self.conditions)
        self.assertAlmostEqual(result.force_scale, 0.85)
        self.assertAlmostEqual(result.mass_flow_scale, 1.025)
        force_gradient = dict(result.force_scale_gradient)
        flow_gradient = dict(result.mass_flow_scale_gradient)
        self.assertAlmostEqual(force_gradient["supply_pressure_pa"], 2.0e-7)
        self.assertAlmostEqual(force_gradient["ambient_pressure_pa"], -1.0e-6)
        self.assertAlmostEqual(flow_gradient["supply_pressure_pa"], 2.0e-7)
        self.assertAlmostEqual(flow_gradient["ambient_pressure_pa"], 5.0e-7)

    def test_three_dimensional_affine_field_is_reproduced_exactly(self) -> None:
        axes = (
            PerformanceAxis("x", "1", (0.0, 1.0)),
            PerformanceAxis("y", "1", (0.0, 2.0)),
            PerformanceAxis("z", "1", (-1.0, 1.0)),
        )
        points = tuple(
            PerformanceScalePoint(
                (x, y, z),
                2.0 + 0.2 * x + 0.3 * y - 0.1 * z,
                1.0 + 0.1 * x + 0.05 * y + 0.2 * z,
            )
            for x in axes[0].values
            for y in axes[1].values
            for z in axes[2].values
        )
        surface = RectilinearPerformanceSurface(
            axes,
            points,
            source_id="three-dimensional-affine-check",
            source_sha256="d" * 64,
        )
        result = surface.evaluate({"x": 0.25, "y": 1.5, "z": 0.4})
        self.assertAlmostEqual(result.force_scale, 2.46)
        self.assertAlmostEqual(result.mass_flow_scale, 1.18)
        for name, expected in {"x": 0.2, "y": 0.3, "z": -0.1}.items():
            self.assertAlmostEqual(dict(result.force_scale_gradient)[name], expected)
        for name, expected in {"x": 0.1, "y": 0.05, "z": 0.2}.items():
            self.assertAlmostEqual(dict(result.mass_flow_scale_gradient)[name], expected)

    def test_surface_round_trip_and_tamper_refusal(self) -> None:
        reopened = RectilinearPerformanceSurface.from_json(self.surface.to_json())
        self.assertEqual(reopened, self.surface)
        self.assertEqual(reopened.surface_hash, self.surface.surface_hash)
        payload = self.surface.to_dict()
        payload["points"][0]["force_scale"] = 99.0  # type: ignore[index]
        with self.assertRaises(InputError):
            RectilinearPerformanceSurface.from_dict(payload)
        with self.assertRaises(InputError):
            RectilinearPerformanceSurface.from_json(json.dumps([]))

    def test_complete_grid_and_validation_domain_fail_closed(self) -> None:
        with self.assertRaisesRegex(DomainError, "complete Cartesian grid"):
            RectilinearPerformanceSurface(
                axes=self.surface.axes,
                points=self.surface.points[:-1],
                source_id="incomplete",
                source_sha256="c" * 64,
            )
        with self.assertRaisesRegex(DomainError, "extrapolation refused"):
            self.surface.evaluate(
                {**self.conditions, "supply_pressure_pa": 2_000_001.0}
            )
        with self.assertRaisesRegex(DomainError, "match surface axes"):
            self.surface.evaluate({"supply_pressure_pa": 1_500_000.0})
        with self.assertRaisesRegex(DomainError, "match surface axes"):
            self.surface.evaluate({**self.conditions, "ignored": 1.0})

    def test_validity_surfaces_are_positive_inside_and_zero_on_boundary(self) -> None:
        surfaces = self.surface.validity_surfaces()
        self.assertEqual(len(surfaces), 4)
        self.assertEqual([surface.bound for surface in surfaces], ["lower", "upper"] * 2)
        for surface in surfaces:
            self.assertGreater(surface.value(self.conditions), 0.0)
            boundary = dict(self.conditions)
            boundary[surface.condition_name] = surface.boundary_value
            self.assertEqual(surface.value(boundary), 0.0)


class PerformanceSurfaceBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = example_surface()
        self.bridge = PropagatorPropulsionBridge(
            example_artifact(),
            performance_surface=self.surface,
            require_complete_stream_binding=False,
        )
        self.conditions = {
            "supply_pressure_pa": 1_500_000.0,
            "ambient_pressure_pa": 50_000.0,
        }

    def evaluate(self, conditions: dict[str, float]):
        return self.bridge.evaluate(
            relative_time_s=0.5,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            operating_conditions=conditions,
        )

    def test_surface_scales_force_and_flow_without_collapsing_isp_change(self) -> None:
        value = self.evaluate(self.conditions)
        self.assertAlmostEqual(value.axial_thrust_n, 150.0 * 0.85)
        self.assertAlmostEqual(value.acceleration_m_s2[0], 1.275)
        self.assertAlmostEqual(value.mass_derivative_kg_s, -3.0 * 1.025)
        self.assertEqual(value.operating_conditions, tuple(self.conditions.items()))
        self.assertAlmostEqual(value.performance_force_scale, 0.85)
        self.assertAlmostEqual(value.performance_mass_flow_scale, 1.025)

    def test_condition_partials_match_central_finite_differences(self) -> None:
        nominal = self.evaluate(self.conditions)
        acceleration_partials = dict(nominal.partials.acceleration_wrt_conditions)
        mass_partials = dict(nominal.partials.mass_rate_wrt_conditions)
        steps = {"supply_pressure_pa": 1.0, "ambient_pressure_pa": 0.1}
        for name, step in steps.items():
            plus_conditions = {**self.conditions, name: self.conditions[name] + step}
            minus_conditions = {**self.conditions, name: self.conditions[name] - step}
            plus = self.evaluate(plus_conditions)
            minus = self.evaluate(minus_conditions)
            finite_acceleration = (
                plus.acceleration_m_s2[0] - minus.acceleration_m_s2[0]
            ) / (2.0 * step)
            finite_mass = (
                plus.mass_derivative_kg_s - minus.mass_derivative_kg_s
            ) / (2.0 * step)
            self.assertAlmostEqual(
                acceleration_partials[name][0], finite_acceleration, places=11
            )
            self.assertAlmostEqual(mass_partials[name], finite_mass, places=11)

    def test_surface_domain_is_exposed_as_propagator_events(self) -> None:
        surfaces = self.bridge.operating_condition_surfaces()
        self.assertEqual(surfaces, self.surface.validity_surfaces())
        self.assertTrue(all(surface.stops_burn for surface in surfaces))
        self.assertEqual(
            self.bridge.manifest()["performance_surface_hash"],
            self.surface.surface_hash,
        )

    def test_conditions_are_required_only_while_surface_burn_is_active(self) -> None:
        with self.assertRaisesRegex(DomainError, "match surface axes"):
            self.bridge.evaluate(
                relative_time_s=0.5,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
            )
        off = self.bridge.evaluate(
            relative_time_s=-1.0,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
        )
        self.assertFalse(off.active)
        self.assertEqual(off.mass_derivative_kg_s, 0.0)


if __name__ == "__main__":
    unittest.main()
