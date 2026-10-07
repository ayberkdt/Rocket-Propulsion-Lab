"""Primary Sidera provider hand-off, state-layout, and body-wrench tests."""

from __future__ import annotations

import copy
import unittest
from dataclasses import replace

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.integrations import (
    SideraCapabilities,
    SideraPropulsionProviderHandoff,
    missing_sidera_propulsion_capabilities,
    require_sidera_propulsion_capabilities,
)
from rocket_propulsion.propulsion.burns import (
    BodyFrameEngineMount,
    EngineVectorCommand,
    MountedEngineCluster,
    PropagatorPropulsionBridge,
    StreamStateBinding,
)
from tests.test_burn_propagator import example_artifact
from tests.test_burn_stream_performance import example_stream_surface


def example_handoff(*, with_cluster: bool = False) -> SideraPropulsionProviderHandoff:
    bridge = PropagatorPropulsionBridge(
        example_artifact(),
        stream_performance_surface=example_stream_surface(),
        stream_bindings=(
            StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
            StreamStateBinding("fuel", "fuel_mass_kg"),
        ),
    )
    cluster = None
    if with_cluster:
        cluster = MountedEngineCluster(
            (
                BodyFrameEngineMount(
                    "main",
                    (0.0, 1.0, 0.0),
                    (1.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    0.2,
                    0.2,
                ),
            )
        )
    return SideraPropulsionProviderHandoff(
        bridge,
        ignition_epoch_s=1_000.0,
        additional_state_order=("oxidizer_mass_kg", "fuel_mass_kg"),
        engine_cluster=cluster,
    )


class SideraPropulsionProviderHandoffTests(unittest.TestCase):
    def setUp(self) -> None:
        self.conditions = {
            "oxidizer_pressure_pa": 1.8e6,
            "fuel_pressure_pa": 1.6e6,
        }
        self.states = {"oxidizer_mass_kg": 20.0, "fuel_mass_kg": 10.0}

    def test_absolute_epoch_state_order_and_nested_round_trip(self) -> None:
        handoff = example_handoff()
        value = handoff.evaluate(
            epoch_s=1_000.5,
            vehicle_mass_kg=100.0,
            direction=(2.0, 0.0, 0.0),
            additional_states=self.states,
            operating_conditions=self.conditions,
        )
        self.assertAlmostEqual(value.effective_profile_time_s, 0.5)
        self.assertEqual(
            handoff.ordered_additional_state_derivatives(value),
            (
                dict(value.additional_state_derivatives)["oxidizer_mass_kg"],
                dict(value.additional_state_derivatives)["fuel_mass_kg"],
            ),
        )
        reopened = SideraPropulsionProviderHandoff.from_json(handoff.to_json())
        self.assertEqual(reopened.handoff_hash, handoff.handoff_hash)
        self.assertEqual(reopened.bridge.model_hash, handoff.bridge.model_hash)
        self.assertEqual(
            reopened.required_operating_conditions,
            ("oxidizer_pressure_pa", "fuel_pressure_pa"),
        )

    def test_manifest_exposes_events_parameters_and_host_ownership(self) -> None:
        handoff = example_handoff(with_cluster=True)
        manifest = handoff.manifest()
        self.assertEqual(manifest["force_frame"], "vehicle_body")
        self.assertEqual(manifest["attitude_owned_by"], "sidera")
        self.assertIn("additional_states", manifest["required_sidera_capabilities"])
        self.assertIn("body_frame_force_moment", manifest["required_sidera_capabilities"])
        self.assertTrue(manifest["event_ids"])
        self.assertEqual(
            manifest["estimable_parameters"],
            ["thrust_scale", "mass_flow_scale", "ignition_time_bias_s"],
        )

    def test_body_wrench_closes_scalar_thrust_and_moving_centre_of_mass(self) -> None:
        handoff = example_handoff(with_cluster=True)
        scalar = handoff.evaluate(
            epoch_s=1_000.5,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
            additional_states=self.states,
            operating_conditions=self.conditions,
        )
        command = EngineVectorCommand("main", scalar.axial_thrust_n, gimbal_axis_2_rad=0.1)
        value = handoff.evaluate_body_frame(
            epoch_s=1_000.5,
            vehicle_mass_kg=100.0,
            commands=(command,),
            reference_point_body_m=(0.0, 0.0, 0.0),
            additional_states=self.states,
            operating_conditions=self.conditions,
        )
        self.assertAlmostEqual(value.force_closure_residual_n, 0.0)
        self.assertAlmostEqual(value.wrench.scalar_thrust_sum_n, scalar.axial_thrust_n)
        self.assertNotEqual(value.wrench.moment_body_n_m, (0.0, 0.0, 0.0))

        with self.assertRaisesRegex(DomainError, "does not close"):
            handoff.evaluate_body_frame(
                epoch_s=1_000.5,
                vehicle_mass_kg=100.0,
                commands=(EngineVectorCommand("main", scalar.axial_thrust_n + 1.0),),
                reference_point_body_m=(0.0, 0.0, 0.0),
                additional_states=self.states,
                operating_conditions=self.conditions,
            )

    def test_layout_and_hash_tampering_fail_closed(self) -> None:
        handoff = example_handoff(with_cluster=True)
        with self.assertRaisesRegex(DomainError, "cover every bound"):
            SideraPropulsionProviderHandoff(
                handoff.bridge,
                ignition_epoch_s=0.0,
                additional_state_order=("fuel_mass_kg",),
            )
        payload = copy.deepcopy(handoff.to_dict())
        payload["ignition_epoch_s"] = 2_000.0
        with self.assertRaises(InputError):
            SideraPropulsionProviderHandoff.from_dict(payload)

    def test_capability_gate_reports_exact_host_gaps(self) -> None:
        handoff = example_handoff(with_cluster=True)
        unavailable = SideraCapabilities(
            status="incompatible",
            version="test",
            finite_burn_available=True,
            maneuver_plan_available=True,
            constructor_compatible=True,
            supported_frames=("inertial",),
            constant_burn_exact=True,
            tabulated_burn_available=False,
            independent_mass_flow_available=False,
        )
        missing = missing_sidera_propulsion_capabilities(unavailable, handoff)
        self.assertEqual(
            missing,
            (
                "propulsion_provider_available",
                "independent_mass_flow_available",
                "additional_states_available",
                "event_detectors_available",
                "analytic_partials_available",
                "body_frame_force_moment_available",
            ),
        )
        with self.assertRaises(FeatureUnavailableError):
            require_sidera_propulsion_capabilities(unavailable, handoff)

        supported = replace(
            unavailable,
            independent_mass_flow_available=True,
            propulsion_provider_available=True,
            additional_states_available=True,
            event_detectors_available=True,
            analytic_partials_available=True,
            body_frame_force_moment_available=True,
        )
        self.assertEqual(missing_sidera_propulsion_capabilities(supported, handoff), ())
        require_sidera_propulsion_capabilities(supported, handoff)


if __name__ == "__main__":
    unittest.main()
