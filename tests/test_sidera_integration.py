"""Boundary tests for exact constant-burn Sidera hand-off."""

import unittest
from math import sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.integrations import (
    export_sidera_artifact,
    inspect_sidera_capabilities,
    to_native_finite_burn,
    to_native_maneuver_plan,
)
from rocket_propulsion.propulsion import STANDARD_GRAVITY_M_S2
from rocket_propulsion.propulsion.burns import BurnTarget, BurnTargetKind, simulate_constant_burn
from tests.test_burns import example_definition, example_operating_point


class SideraAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
            example_operating_point(),
        )

    def test_artifact_preserves_force_mass_law_and_source_hash(self) -> None:
        artifact = export_sidera_artifact(
            self.result,
            t_start_s=50.0,
            direction=(2.0, 0.0, 0.0),
            frame="ric",
        )
        burn = artifact["maneuvers"][0]
        self.assertEqual(artifact["schema"], "sidera_maneuver_plan_v1")
        self.assertEqual(artifact["manifest"]["source_result_hash"], self.result.result_hash)
        self.assertEqual(burn["direction"], [1.0, 0.0, 0.0])
        self.assertEqual(burn["throttle"], 1.0)
        sidera_mass_flow = burn["thrust_n"] / (STANDARD_GRAVITY_M_S2 * burn["isp_s"])
        self.assertAlmostEqual(
            sidera_mass_flow, self.result.operating_point.total_tank_flow_kg_s
        )

    def test_direction_is_normalized_and_invalid_frame_refused(self) -> None:
        artifact = export_sidera_artifact(
            self.result,
            t_start_s=0.0,
            direction=(1.0, 2.0, 2.0),
        )
        direction = artifact["maneuvers"][0]["direction"]
        self.assertAlmostEqual(sqrt(sum(value * value for value in direction)), 1.0)
        with self.assertRaises(DomainError):
            export_sidera_artifact(
                self.result,
                t_start_s=0.0,
                direction=(1.0, 0.0, 0.0),
                frame="body",  # type: ignore[arg-type]
            )

    def test_canted_burn_exports_net_force_without_losing_tank_drain(self) -> None:
        result = simulate_constant_burn(
            example_definition(
                BurnTarget(BurnTargetKind.DURATION, 1.0),
                cant_efficiency=0.8,
            ),
            example_operating_point(),
        )
        artifact = export_sidera_artifact(
            result,
            t_start_s=0.0,
            direction=(1.0, 0.0, 0.0),
        )
        burn = artifact["maneuvers"][0]
        self.assertAlmostEqual(
            burn["thrust_n"],
            result.summary.axial_total_impulse_n_s / result.summary.duration_s,
        )
        self.assertAlmostEqual(
            burn["thrust_n"] / (STANDARD_GRAVITY_M_S2 * burn["isp_s"]),
            result.operating_point.total_tank_flow_kg_s,
        )
        self.assertEqual(
            artifact["manifest"]["force_semantics"],
            "net_axial_thrust_after_hardware_cant",
        )

    def test_capability_probe_is_typed_even_without_sidera(self) -> None:
        capabilities = inspect_sidera_capabilities()
        self.assertIn(capabilities.status, {"unavailable", "incompatible", "supported"})
        self.assertEqual(
            capabilities.constant_burn_exact,
            capabilities.constructor_compatible and capabilities.maneuver_plan_available,
        )

    def test_native_public_contract_when_sidera_is_on_python_path(self) -> None:
        capabilities = inspect_sidera_capabilities()
        if not capabilities.constant_burn_exact:
            self.skipTest("Compatible Sidera installation is not on the explicit test path.")
        burn = to_native_finite_burn(
            self.result,
            t_start_s=12.0,
            direction=(0.0, 1.0, 0.0),
            frame="vnb",
        )
        plan = to_native_maneuver_plan(
            self.result,
            t_start_s=12.0,
            direction=(0.0, 1.0, 0.0),
            frame="vnb",
        )
        self.assertAlmostEqual(
            burn.mass_flow_kg_s, self.result.operating_point.total_tank_flow_kg_s
        )
        self.assertEqual(burn.throttle, 1.0)
        self.assertTrue(plan.has_finite_burns)


if __name__ == "__main__":
    unittest.main()
