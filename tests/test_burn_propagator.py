"""Propagation contract tests: events, conservation, cutoffs, and Jacobians."""

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    BurnEventKind,
    PropagatorPropulsionBridge,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def example_artifact() -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                0.0,
                1.0,
                120.0,
                240.0,
                100.0,
                200.0,
                2.0,
                4.0,
                (("oxidizer", 1.5), ("fuel", 0.5)),
                (("oxidizer", 3.0), ("fuel", 1.0)),
                "ramp",
            ),
            TabulatedBurnSegment(
                1.0,
                2.0,
                180.0,
                180.0,
                150.0,
                150.0,
                3.0,
                3.0,
                (("oxidizer", 2.25), ("fuel", 0.75)),
                (("oxidizer", 2.25), ("fuel", 0.75)),
                "steady",
            ),
        ),
        source_id="propagator-test",
        source_sha256="a" * 64,
        reference_ids=("NASA-20040000363",),
    )


def example_bridge() -> PropagatorPropulsionBridge:
    return PropagatorPropulsionBridge(
        example_artifact(),
        stream_bindings=(
            StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
            StreamStateBinding("fuel", "fuel_mass_kg"),
        ),
        state_constraints=(
            AdditionalStateConstraint("oxidizer_mass_kg", 2.0),
            AdditionalStateConstraint("fuel_mass_kg", 1.0),
        ),
        protected_dry_mass_kg=50.0,
    )


class PropagatorBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.bridge = example_bridge()
        self.states = {"oxidizer_mass_kg": 10.0, "fuel_mass_kg": 5.0}

    def evaluate(self, time_s: float, **kwargs: object):
        arguments = {
            "relative_time_s": time_s,
            "vehicle_mass_kg": 100.0,
            "direction": (3.0, 4.0, 0.0),
            "additional_states": self.states,
        }
        arguments.update(kwargs)
        return self.bridge.evaluate(**arguments)  # type: ignore[arg-type]

    def test_force_total_mass_and_bound_tank_derivatives_close(self) -> None:
        value = self.evaluate(0.5)
        self.assertTrue(value.active)
        self.assertEqual(value.phase, "ramp")
        self.assertAlmostEqual(value.axial_thrust_n, 150.0)
        self.assertEqual(value.force_vector_n, (90.0, 120.0, 0.0))
        self.assertEqual(value.acceleration_m_s2, (0.9, 1.2, 0.0))
        self.assertAlmostEqual(value.mass_derivative_kg_s, -3.0)
        drains = dict(value.additional_state_derivatives)
        self.assertAlmostEqual(drains["oxidizer_mass_kg"], -2.25)
        self.assertAlmostEqual(drains["fuel_mass_kg"], -0.75)
        self.assertAlmostEqual(sum(drains.values()), value.mass_derivative_kg_s)
        for actual, expected in zip(
            value.partials.acceleration_wrt_mass, (-0.009, -0.012, -0.0), strict=True
        ):
            self.assertAlmostEqual(actual, expected)

    def test_discontinuities_have_explicit_left_and_right_limits(self) -> None:
        left = self.evaluate(1.0, side="left")
        right = self.evaluate(1.0, side="right")
        self.assertAlmostEqual(left.axial_thrust_n, 200.0)
        self.assertAlmostEqual(right.axial_thrust_n, 150.0)
        self.assertEqual(left.partials.nonsmooth_parameters, ("ignition_time_bias_s",))
        self.assertEqual(right.partials.nonsmooth_parameters, ("ignition_time_bias_s",))
        self.assertFalse(self.evaluate(0.0, side="left").active)
        self.assertTrue(self.evaluate(0.0, side="right").active)
        self.assertTrue(self.evaluate(2.0, side="left").active)
        self.assertFalse(self.evaluate(2.0, side="right").active)

    def test_event_surfaces_cover_time_dry_mass_and_each_reserve(self) -> None:
        surfaces = self.bridge.event_surfaces()
        self.assertEqual(
            [surface.kind for surface in surfaces],
            [
                BurnEventKind.START,
                BurnEventKind.PROFILE_BOUNDARY,
                BurnEventKind.STOP,
                BurnEventKind.DRY_MASS,
                BurnEventKind.TANK_RESERVE,
                BurnEventKind.TANK_RESERVE,
            ],
        )
        self.assertAlmostEqual(
            surfaces[1].value(
                relative_time_s=1.2,
                vehicle_mass_kg=100.0,
                additional_states=self.states,
                ignition_time_bias_s=0.2,
            ),
            0.0,
        )
        self.assertEqual(
            surfaces[3].value(
                relative_time_s=0.0,
                vehicle_mass_kg=50.0,
                additional_states=self.states,
            ),
            0.0,
        )
        self.assertEqual(
            surfaces[4].value(
                relative_time_s=0.0,
                vehicle_mass_kg=100.0,
                additional_states={**self.states, "oxidizer_mass_kg": 2.0},
            ),
            0.0,
        )

    def test_state_constraints_inhibit_burn_without_stopping_propagation(self) -> None:
        depleted = self.evaluate(
            0.5,
            additional_states={"oxidizer_mass_kg": 2.0, "fuel_mass_kg": 5.0},
        )
        self.assertFalse(depleted.active)
        self.assertEqual(depleted.inhibited_by, ("oxidizer_mass_kg",))
        self.assertEqual(depleted.force_vector_n, (0.0, 0.0, 0.0))
        self.assertEqual(depleted.mass_derivative_kg_s, 0.0)
        dry = self.evaluate(0.5, vehicle_mass_kg=50.0)
        self.assertFalse(dry.active)
        self.assertEqual(dry.inhibited_by, ("protected_dry_mass_kg",))

    def test_parameter_jacobians_match_central_finite_differences(self) -> None:
        nominal = self.evaluate(0.5)
        acceleration_partials = dict(nominal.partials.acceleration_wrt_parameters)
        mass_partials = dict(nominal.partials.mass_rate_wrt_parameters)
        step = 1e-6
        for parameter in ("thrust_scale", "mass_flow_scale", "ignition_time_bias_s"):
            reference = dict(nominal.parameters)[parameter]
            plus = self.evaluate(0.5, parameter_overrides={parameter: reference + step})
            minus = self.evaluate(0.5, parameter_overrides={parameter: reference - step})
            finite_acceleration = tuple(
                (upper - lower) / (2.0 * step)
                for upper, lower in zip(
                    plus.acceleration_m_s2, minus.acceleration_m_s2, strict=True
                )
            )
            for actual, expected in zip(
                acceleration_partials[parameter], finite_acceleration, strict=True
            ):
                self.assertAlmostEqual(actual, expected, places=8)
            finite_mass = (
                plus.mass_derivative_kg_s - minus.mass_derivative_kg_s
            ) / (2.0 * step)
            self.assertAlmostEqual(mass_partials[parameter], finite_mass, places=8)

    def test_time_bias_shifts_profile_and_all_time_events_together(self) -> None:
        nominal = self.evaluate(0.5)
        shifted = self.evaluate(0.7, parameter_overrides={"ignition_time_bias_s": 0.2})
        self.assertEqual(nominal.axial_thrust_n, shifted.axial_thrust_n)
        self.assertEqual(nominal.mass_derivative_kg_s, shifted.mass_derivative_kg_s)
        for surface in self.bridge.event_surfaces()[:3]:
            assert surface.relative_time_s is not None
            self.assertEqual(
                surface.value(
                    relative_time_s=surface.relative_time_s + 0.2,
                    vehicle_mass_kg=100.0,
                    additional_states=self.states,
                    ignition_time_bias_s=0.2,
                ),
                0.0,
            )

    def test_manifest_binds_artifact_semantics_and_configuration(self) -> None:
        manifest = self.bridge.manifest()
        self.assertEqual(manifest["model_hash"], self.bridge.model_hash)
        self.assertEqual(manifest["artifact_hash"], self.bridge.artifact.artifact_hash)
        self.assertIn("OREKIT-13.1.5-PROPULSION-MODEL", manifest["reference_ids"])
        changed = PropagatorPropulsionBridge(
            example_artifact(),
            stream_bindings=self.bridge.stream_bindings,
            state_constraints=self.bridge.state_constraints,
            protected_dry_mass_kg=51.0,
        )
        self.assertNotEqual(changed.model_hash, self.bridge.model_hash)

    def test_invalid_or_incomplete_contracts_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            PropagatorPropulsionBridge(
                example_artifact(),
                stream_bindings=(StreamStateBinding("oxidizer", "oxidizer_mass_kg"),),
            )
        with self.assertRaises(DomainError):
            self.evaluate(0.5, additional_states={"oxidizer_mass_kg": 10.0})
        with self.assertRaises(DomainError):
            self.evaluate(
                0.5,
                additional_states={
                    "oxidizer_mass_kg": "invalid",
                    "fuel_mass_kg": 5.0,
                },
            )
        with self.assertRaises(DomainError):
            self.evaluate(0.5, parameter_overrides={"unknown": 1.0})
        with self.assertRaises(DomainError):
            self.evaluate(0.5, parameter_overrides={"thrust_scale": 3.0})
        with self.assertRaises(DomainError):
            self.evaluate(0.5, parameter_overrides={"thrust_scale": "invalid"})
        with self.assertRaises(DomainError):
            self.evaluate(0.5, direction=(0.0, 0.0, 0.0))
        with self.assertRaises(DomainError):
            self.bridge.event_surfaces()[0].value(
                relative_time_s=0.0,
                vehicle_mass_kg=100.0,
                additional_states=self.states,
                ignition_time_bias_s=float("nan"),
            )


if __name__ == "__main__":
    unittest.main()
