"""Integrity and physical-channel tests for future Sidera ensemble hand-off."""

from __future__ import annotations

import copy
import unittest

from rocket_propulsion.core.errors import FeatureUnavailableError, InputError
from rocket_propulsion.integrations import (
    SideraCapabilities,
    SideraEnsembleHandoff,
    export_sidera_ensemble_handoff,
    missing_sidera_ensemble_capabilities,
    require_sidera_ensemble_capabilities,
)
from rocket_propulsion.propulsion.burns import (
    ConstantPerformanceProvider,
    DistributionKind,
    EngineClusterMember,
    EngineTransientCommand,
    PropulsionScenario,
    PropulsionScenarioSet,
    ScenarioConditionalUncertainty,
    ScenarioEvent,
    ScenarioEventType,
    ThrottleSchedule,
    ThrottleSegment,
    UncertainParameter,
    UncertaintyClass,
    build_hierarchical_propulsion_ensemble,
    build_propulsion_scenario_ensemble,
    hierarchical_ensemble_hash,
)


class SideraEnsembleHandoffTests(unittest.TestCase):
    def setUp(self) -> None:
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=200.0,
            chamber_pressure_pa=1.0e6,
            source="Sidera ensemble fixture",
        ).operating_point()
        members = (
            EngineClusterMember("a", point, cant_efficiency=0.8),
            EngineClusterMember("b", point, cant_efficiency=0.8),
        )
        schedule = ThrottleSchedule((ThrottleSegment(5.0, 1.0, 1.0),))
        commands = tuple(
            EngineTransientCommand(member.engine_id, schedule) for member in members
        )
        scenario_set = PropulsionScenarioSet(
            scenarios=(
                PropulsionScenario(
                    "nominal",
                    0.9,
                    (),
                    "handoff event-tree fixture",
                    "both engines operate",
                ),
                PropulsionScenario(
                    "common-cause-out",
                    0.1,
                    (
                        ScenarioEvent(
                            ScenarioEventType.ENGINE_DISABLED,
                            ("a", "b"),
                            None,
                            "handoff failure fixture",
                            "shared initiation unavailable",
                        ),
                    ),
                    "handoff event-tree fixture",
                    "both engines unavailable",
                ),
            ),
            objective="verify future Sidera weighted hand-off",
        )
        discrete = build_propulsion_scenario_ensemble(
            members,
            commands,
            scenario_set,
        )
        thrust_scale = UncertainParameter(
            name="thrust_scale",
            nominal=1.0,
            lower=0.98,
            upper=1.02,
            distribution=DistributionKind.UNIFORM,
            uncertainty_class=UncertaintyClass.EPISTEMIC,
            source="handoff qualification fixture",
            rationale="exercise conditional performance histories",
        )
        definitions = (
            ScenarioConditionalUncertainty(
                "nominal",
                (thrust_scale,),
                4,
                "nominal handoff fixture",
                "nominal performance variation",
            ),
            ScenarioConditionalUncertainty(
                "common-cause-out",
                (),
                1,
                "zero-thrust handoff fixture",
                "disabled engines remain deterministic",
            ),
        )
        self.ensemble = build_hierarchical_propulsion_ensemble(
            discrete,
            definitions,
            seed=83,
            replicate_tolerance=0.5,
        )

    def test_export_preserves_weights_zero_thrust_and_physical_channels(self) -> None:
        handoff = export_sidera_ensemble_handoff(
            self.ensemble,
            t_start_s=120.0,
            direction=(0.0, 3.0, 4.0),
            frame="vnb",
        )
        payload = handoff.to_dict()
        self.assertEqual(handoff.source_ensemble_hash, hierarchical_ensemble_hash(self.ensemble))
        self.assertEqual(len(handoff.jobs), 5)
        self.assertAlmostEqual(handoff.joint_probability_sum, 1.0)
        self.assertEqual(payload["force_channel"], "axial_thrust_n")
        self.assertEqual(
            payload["mass_flow_channel"],
            "total_tank_mass_flow_kg_s",
        )
        self.assertEqual(handoff.jobs[0].direction, (0.0, 0.6, 0.8))
        self.assertTrue(
            any(
                item.artifact.delivered_total_impulse_n_s == 0.0
                and item.joint_probability == 0.1
                for item in handoff.jobs
            )
        )
        nominal = next(
            item for item in handoff.jobs if item.scenario_id == "nominal"
        )
        self.assertLess(
            nominal.artifact.axial_total_impulse_n_s,
            nominal.artifact.delivered_total_impulse_n_s,
        )

    def test_round_trip_and_nested_tamper_refusal(self) -> None:
        handoff = export_sidera_ensemble_handoff(
            self.ensemble,
            t_start_s=10.0,
            direction=(1.0, 0.0, 0.0),
        )
        reopened = SideraEnsembleHandoff.from_json(handoff.to_json())
        self.assertEqual(reopened, handoff)
        self.assertEqual(reopened.handoff_hash, handoff.handoff_hash)

        tampered = copy.deepcopy(handoff.to_dict())
        tampered["jobs"][0]["artifact"]["warnings"].append(
            "unhashed tamper"
        )
        with self.assertRaisesRegex(InputError, "artifact hash"):
            SideraEnsembleHandoff.from_dict(tampered)

        changed_semantics = copy.deepcopy(handoff.to_dict())
        changed_semantics["force_channel"] = "delivered_thrust_n"
        with self.assertRaisesRegex(InputError, "physical semantics"):
            SideraEnsembleHandoff.from_dict(changed_semantics)

    def test_capability_gate_names_every_missing_native_feature(self) -> None:
        capabilities = SideraCapabilities(
            status="supported",
            version="fixture",
            finite_burn_available=True,
            maneuver_plan_available=True,
            constructor_compatible=True,
            supported_frames=("inertial", "ric", "vnb"),
            constant_burn_exact=True,
            tabulated_burn_available=False,
            independent_mass_flow_available=False,
            reason=None,
            weighted_maneuver_ensemble_available=False,
            exact_weighted_tabulated_ensemble=False,
        )
        self.assertEqual(
            missing_sidera_ensemble_capabilities(capabilities),
            (
                "tabulated_finite_burn",
                "independent_mass_flow",
                "weighted_maneuver_ensemble",
            ),
        )
        with self.assertRaises(FeatureUnavailableError) as context:
            require_sidera_ensemble_capabilities(capabilities)
        self.assertEqual(
            context.exception.details["missing_capabilities"],
            [
                "tabulated_finite_burn",
                "independent_mass_flow",
                "weighted_maneuver_ensemble",
            ],
        )


if __name__ == "__main__":
    unittest.main()

