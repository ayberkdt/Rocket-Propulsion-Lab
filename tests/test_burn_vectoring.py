"""Force, moment, integrity, and derivative tests for mounted engine clusters."""

from __future__ import annotations

import copy
import math
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns.vectoring import (
    BodyFrameEngineMount,
    EngineVectorCommand,
    MountedEngineCluster,
)


def _component(vector: tuple[float, float, float], index: int) -> float:
    return float(vector[index])


class MountedEngineClusterTests(unittest.TestCase):
    """Verify exact wrench geometry and fail-closed public contracts."""

    def setUp(self) -> None:
        self.mount = BodyFrameEngineMount(
            engine_id="main",
            application_point_body_m=(0.0, 1.0, 0.0),
            nominal_direction_body=(1.0, 0.0, 0.0),
            gimbal_axis_1_body=(0.0, 1.0, 0.0),
            gimbal_axis_1_limit_rad=0.2,
            gimbal_axis_2_limit_rad=0.2,
        )
        self.cluster = MountedEngineCluster((self.mount,))

    def test_single_engine_force_and_moment_about_moving_centre_of_mass(self) -> None:
        origin = self.cluster.evaluate((EngineVectorCommand("main", 100.0),))
        self.assertEqual(origin.force_body_n, (100.0, 0.0, 0.0))
        self.assertEqual(origin.moment_body_n_m, (0.0, 0.0, -100.0))
        self.assertEqual(origin.resultant_direction_body, (1.0, 0.0, 0.0))
        self.assertAlmostEqual(origin.vector_cancellation_loss_n, 0.0)

        at_engine = self.cluster.evaluate(
            (EngineVectorCommand("main", 100.0),),
            reference_point_body_m=(0.0, 1.0, 0.0),
        )
        self.assertEqual(at_engine.moment_body_n_m, (0.0, 0.0, 0.0))

    def test_symmetric_cluster_cancels_moment_and_reports_vector_loss(self) -> None:
        mounts = (
            self.mount,
            BodyFrameEngineMount(
                engine_id="opposite-mount",
                application_point_body_m=(0.0, -1.0, 0.0),
                nominal_direction_body=(1.0, 0.0, 0.0),
                gimbal_axis_1_body=(0.0, 1.0, 0.0),
                gimbal_axis_1_limit_rad=0.2,
                gimbal_axis_2_limit_rad=0.2,
            ),
        )
        cluster = MountedEngineCluster(mounts)
        aligned = cluster.evaluate(
            (
                EngineVectorCommand("main", 100.0),
                EngineVectorCommand("opposite-mount", 100.0),
            )
        )
        self.assertEqual(aligned.force_body_n, (200.0, 0.0, 0.0))
        self.assertEqual(aligned.moment_body_n_m, (0.0, 0.0, 0.0))

        canted = cluster.evaluate(
            (
                EngineVectorCommand("main", 100.0, gimbal_axis_2_rad=0.1),
                EngineVectorCommand("opposite-mount", 100.0, gimbal_axis_2_rad=-0.1),
            )
        )
        self.assertLess(canted.resultant_thrust_n, canted.scalar_thrust_sum_n)
        self.assertGreater(canted.vector_cancellation_loss_n, 0.0)

    def test_opposed_engines_can_have_zero_resultant_without_fake_direction(self) -> None:
        fixed = MountedEngineCluster(
            (
                BodyFrameEngineMount("plus", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                BodyFrameEngineMount("minus", (0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            )
        )
        result = fixed.evaluate(
            (EngineVectorCommand("plus", 50.0), EngineVectorCommand("minus", 50.0))
        )
        self.assertEqual(result.force_body_n, (0.0, 0.0, 0.0))
        self.assertEqual(result.resultant_thrust_n, 0.0)
        self.assertIsNone(result.resultant_direction_body)
        self.assertEqual(result.vector_cancellation_loss_n, 100.0)

    def test_gimbal_force_and_moment_partials_match_finite_difference(self) -> None:
        nominal = EngineVectorCommand(
            "main",
            123.0,
            gimbal_axis_1_rad=0.07,
            gimbal_axis_2_rad=-0.04,
        )
        result = self.cluster.evaluate((nominal,), reference_point_body_m=(0.1, -0.2, 0.3))
        engine = result.engines[0]

        cases = (
            (
                "thrust_n",
                1e-4,
                engine.partials.force_wrt_thrust,
                engine.partials.moment_wrt_thrust,
            ),
            (
                "gimbal_axis_1_rad",
                1e-6,
                engine.partials.force_wrt_gimbal_axis_1_rad,
                engine.partials.moment_wrt_gimbal_axis_1_rad,
            ),
            (
                "gimbal_axis_2_rad",
                1e-6,
                engine.partials.force_wrt_gimbal_axis_2_rad,
                engine.partials.moment_wrt_gimbal_axis_2_rad,
            ),
        )
        values = {
            "thrust_n": nominal.thrust_n,
            "gimbal_axis_1_rad": nominal.gimbal_axis_1_rad,
            "gimbal_axis_2_rad": nominal.gimbal_axis_2_rad,
        }
        for name, step, analytic_force, analytic_moment in cases:
            plus_values = dict(values)
            minus_values = dict(values)
            plus_values[name] += step
            minus_values[name] -= step
            plus = self.cluster.evaluate(
                (EngineVectorCommand("main", **plus_values),),
                reference_point_body_m=(0.1, -0.2, 0.3),
            )
            minus = self.cluster.evaluate(
                (EngineVectorCommand("main", **minus_values),),
                reference_point_body_m=(0.1, -0.2, 0.3),
            )
            for index in range(3):
                force_fd = (
                    _component(plus.force_body_n, index) - _component(minus.force_body_n, index)
                ) / (2.0 * step)
                moment_fd = (
                    _component(plus.moment_body_n_m, index)
                    - _component(minus.moment_body_n_m, index)
                ) / (2.0 * step)
                self.assertAlmostEqual(force_fd, analytic_force[index], places=7)
                self.assertAlmostEqual(moment_fd, analytic_moment[index], places=7)

    def test_disabled_engine_has_zero_force_and_command_partials(self) -> None:
        result = self.cluster.evaluate(
            (
                EngineVectorCommand(
                    "main",
                    100.0,
                    gimbal_axis_1_rad=0.1,
                    enabled=False,
                ),
            )
        )
        engine = result.engines[0]
        self.assertFalse(engine.active)
        self.assertEqual(engine.force_body_n, (0.0, 0.0, 0.0))
        self.assertEqual(engine.partials.force_wrt_thrust, (0.0, 0.0, 0.0))
        self.assertEqual(engine.partials.force_wrt_gimbal_axis_1_rad, (0.0, 0.0, 0.0))

    def test_round_trip_manifest_and_semantic_tamper_refusal(self) -> None:
        reopened = MountedEngineCluster.from_json(self.cluster.to_json())
        self.assertEqual(reopened, self.cluster)
        self.assertEqual(reopened.model_sha256, self.cluster.model_sha256)
        manifest = self.cluster.model_manifest()
        self.assertEqual(manifest["model_sha256"], self.cluster.model_sha256)
        self.assertEqual(
            manifest["reference_ids"],
            ["NASA-SP-8059", "NASA-TM-20240012525"],
        )
        self.assertFalse(manifest["owns_attitude_propagation"])

        tampered = copy.deepcopy(self.cluster.to_dict())
        tampered["mounts"][0]["application_point_body_m"][1] = 2.0
        with self.assertRaises(InputError):
            MountedEngineCluster.from_dict(tampered)

        wrong_convention = copy.deepcopy(self.cluster.to_dict())
        wrong_convention["rotation_convention"] = "implicit-unspecified"
        with self.assertRaises(InputError):
            MountedEngineCluster.from_dict(wrong_convention)

    def test_invalid_geometry_limits_and_command_sets_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            BodyFrameEngineMount("bad", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0))
        with self.assertRaises(DomainError):
            BodyFrameEngineMount(
                "bad-limit",
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                gimbal_axis_1_limit_rad=math.pi / 2.0,
            )
        with self.assertRaises(DomainError):
            MountedEngineCluster((self.mount, self.mount))
        with self.assertRaises(DomainError):
            self.cluster.evaluate(())
        with self.assertRaises(DomainError):
            self.cluster.evaluate((EngineVectorCommand("unknown", 1.0),))
        with self.assertRaises(DomainError):
            self.cluster.evaluate((EngineVectorCommand("main", 1.0, gimbal_axis_1_rad=0.21),))
        with self.assertRaises(DomainError):
            self.cluster.evaluate(
                (EngineVectorCommand("main", 1.0),),
                reference_point_body_m=(float("nan"), 0.0, 0.0),
            )


if __name__ == "__main__":
    unittest.main()
