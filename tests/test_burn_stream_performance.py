"""Independent propellant-stream performance maps and bridge tests."""

from __future__ import annotations

import json
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PropagatorPropulsionBridge,
    RectilinearStreamPerformanceSurface,
    StreamPerformanceScalePoint,
    StreamStateBinding,
)
from tests.test_burn_performance_surface import example_surface
from tests.test_burn_propagator import example_artifact


def example_stream_surface(*, fuel_stream_id: str = "fuel") -> RectilinearStreamPerformanceSurface:
    axes = (
        PerformanceAxis("oxidizer_pressure_pa", "Pa", (1.0e6, 2.0e6)),
        PerformanceAxis("fuel_pressure_pa", "Pa", (1.0e6, 2.0e6)),
    )
    return RectilinearStreamPerformanceSurface(
        axes=axes,
        points=tuple(
            StreamPerformanceScalePoint(
                coordinates=(oxidizer_pressure, fuel_pressure),
                force_scale=(
                    0.8 + 1.0e-7 * (oxidizer_pressure - 1.0e6) + 5.0e-8 * (fuel_pressure - 1.0e6)
                ),
                stream_flow_scales=(
                    (
                        "oxidizer",
                        0.9 + 2.0e-7 * (oxidizer_pressure - 1.0e6),
                    ),
                    (
                        fuel_stream_id,
                        1.1 + 1.0e-7 * (fuel_pressure - 1.0e6),
                    ),
                ),
            )
            for oxidizer_pressure in axes[0].values
            for fuel_pressure in axes[1].values
        ),
        source_id="bipropellant-qualification-grid-A",
        source_sha256="e" * 64,
    )


class StreamPerformanceSurfaceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = example_stream_surface()
        self.conditions = {
            "oxidizer_pressure_pa": 1.8e6,
            "fuel_pressure_pa": 1.6e6,
        }

    def test_independent_stream_values_and_gradients_are_interpolated(self) -> None:
        result = self.surface.evaluate(self.conditions)
        self.assertAlmostEqual(result.force_scale, 0.91)
        scales = dict(result.stream_flow_scales)
        self.assertAlmostEqual(scales["fuel"], 1.16)
        self.assertAlmostEqual(scales["oxidizer"], 1.06)
        force_gradient = dict(result.force_scale_gradient)
        self.assertAlmostEqual(force_gradient["oxidizer_pressure_pa"], 1.0e-7)
        self.assertAlmostEqual(force_gradient["fuel_pressure_pa"], 5.0e-8)
        gradients = {
            stream_id: dict(values) for stream_id, values in result.stream_flow_scale_gradients
        }
        self.assertAlmostEqual(gradients["oxidizer"]["oxidizer_pressure_pa"], 2.0e-7)
        self.assertAlmostEqual(gradients["oxidizer"]["fuel_pressure_pa"], 0.0)
        self.assertAlmostEqual(gradients["fuel"]["oxidizer_pressure_pa"], 0.0)
        self.assertAlmostEqual(gradients["fuel"]["fuel_pressure_pa"], 1.0e-7)

    def test_round_trip_domain_events_and_tamper_refusal(self) -> None:
        reopened = RectilinearStreamPerformanceSurface.from_json(self.surface.to_json())
        self.assertEqual(reopened, self.surface)
        self.assertEqual(reopened.surface_hash, self.surface.surface_hash)
        self.assertEqual(len(self.surface.validity_surfaces()), 4)
        self.assertTrue(
            all(
                surface.value(self.conditions) > 0.0 for surface in self.surface.validity_surfaces()
            )
        )
        payload = self.surface.to_dict()
        payload["points"][0]["stream_flow_scales"][0][1] = 99.0  # type: ignore[index]
        with self.assertRaises(InputError):
            RectilinearStreamPerformanceSurface.from_dict(payload)
        with self.assertRaises(InputError):
            RectilinearStreamPerformanceSurface.from_json(json.dumps([]))

    def test_grid_and_channel_identity_fail_closed(self) -> None:
        with self.assertRaisesRegex(DomainError, "complete Cartesian grid"):
            RectilinearStreamPerformanceSurface(
                self.surface.axes,
                self.surface.points[:-1],
                source_id="incomplete",
                source_sha256="f" * 64,
            )
        changed = list(self.surface.points)
        changed[-1] = StreamPerformanceScalePoint(
            changed[-1].coordinates,
            changed[-1].force_scale,
            (("fuel", 1.0), ("unmatched", 1.0)),
        )
        with self.assertRaisesRegex(DomainError, "identical streams"):
            RectilinearStreamPerformanceSurface(
                self.surface.axes,
                tuple(changed),
                source_id="channel-mismatch",
                source_sha256="f" * 64,
            )


class StreamPerformanceBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = example_stream_surface()
        self.bridge = PropagatorPropulsionBridge(
            example_artifact(),
            stream_performance_surface=self.surface,
            stream_bindings=(
                StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
                StreamStateBinding("fuel", "fuel_mass_kg"),
            ),
        )
        self.conditions = {
            "oxidizer_pressure_pa": 1.8e6,
            "fuel_pressure_pa": 1.6e6,
        }
        self.states = {"oxidizer_mass_kg": 20.0, "fuel_mass_kg": 10.0}

    def evaluate(self, conditions: dict[str, float]):
        return self.bridge.evaluate(
            relative_time_s=0.5,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            additional_states=self.states,
            operating_conditions=conditions,
        )

    def test_bridge_preserves_independent_stream_scaling_and_total_mass(self) -> None:
        value = self.evaluate(self.conditions)
        drains = dict(value.additional_state_derivatives)
        self.assertAlmostEqual(drains["oxidizer_mass_kg"], -2.25 * 1.06)
        self.assertAlmostEqual(drains["fuel_mass_kg"], -0.75 * 1.16)
        self.assertAlmostEqual(value.mass_derivative_kg_s, sum(drains.values()))
        self.assertAlmostEqual(value.performance_mass_flow_scale, 1.085)
        self.assertAlmostEqual(value.performance_force_scale, 0.91)
        self.assertAlmostEqual(value.acceleration_m_s2[0], 1.365)
        self.assertEqual(
            self.bridge.manifest()["stream_performance_surface_hash"],
            self.surface.surface_hash,
        )

    def test_total_component_and_acceleration_partials_match_finite_difference(self) -> None:
        nominal = self.evaluate(self.conditions)
        mass_partials = dict(nominal.partials.mass_rate_wrt_conditions)
        acceleration_partials = dict(nominal.partials.acceleration_wrt_conditions)
        state_partials = {
            condition: dict(values)
            for condition, values in nominal.partials.additional_state_rates_wrt_conditions
        }
        for condition in self.conditions:
            step = 1.0
            plus = self.evaluate({**self.conditions, condition: self.conditions[condition] + step})
            minus = self.evaluate({**self.conditions, condition: self.conditions[condition] - step})
            finite_mass = (plus.mass_derivative_kg_s - minus.mass_derivative_kg_s) / (2.0 * step)
            finite_acceleration = (plus.acceleration_m_s2[0] - minus.acceleration_m_s2[0]) / (
                2.0 * step
            )
            self.assertAlmostEqual(mass_partials[condition], finite_mass, places=11)
            self.assertAlmostEqual(
                acceleration_partials[condition][0], finite_acceleration, places=11
            )
            for state_name in self.states:
                finite_state = (
                    dict(plus.additional_state_derivatives)[state_name]
                    - dict(minus.additional_state_derivatives)[state_name]
                ) / (2.0 * step)
                self.assertAlmostEqual(
                    state_partials[condition][state_name], finite_state, places=11
                )

    def test_incompatible_or_ambiguous_surfaces_are_refused(self) -> None:
        with self.assertRaisesRegex(DomainError, "either a total-flow or stream-resolved"):
            PropagatorPropulsionBridge(
                example_artifact(),
                performance_surface=example_surface(),
                stream_performance_surface=self.surface,
            )
        with self.assertRaisesRegex(DomainError, "outputs must match"):
            PropagatorPropulsionBridge(
                example_artifact(),
                stream_performance_surface=example_stream_surface(fuel_stream_id="unmatched"),
                stream_bindings=(
                    StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
                    StreamStateBinding("fuel", "fuel_mass_kg"),
                ),
            )


if __name__ == "__main__":
    unittest.main()
