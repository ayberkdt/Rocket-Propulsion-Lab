"""Tests for scenario-conditioned continuous propulsion ensembles."""

from __future__ import annotations

import copy
import hashlib
import json
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
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
    hierarchical_ensemble_from_dict,
    hierarchical_ensemble_from_json,
    hierarchical_ensemble_hash,
    hierarchical_ensemble_to_dict,
    hierarchical_ensemble_to_json,
)


def uncertain_scale(
    name: str,
    lower: float,
    upper: float,
    distribution: DistributionKind = DistributionKind.UNIFORM,
) -> UncertainParameter:
    return UncertainParameter(
        name=name,
        nominal=1.0,
        lower=lower,
        upper=upper,
        distribution=distribution,
        uncertainty_class=UncertaintyClass.EPISTEMIC,
        source="conditional qualification fixture",
        rationale="verify conditional propulsion uncertainty propagation",
    )


class HierarchicalPropulsionEnsembleTests(unittest.TestCase):
    def setUp(self) -> None:
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=200.0,
            chamber_pressure_pa=1_000_000.0,
            tank_id="shared",
            source="hierarchical test point",
        ).operating_point()
        members = (
            EngineClusterMember("left", point),
            EngineClusterMember("right", point),
        )
        schedule = ThrottleSchedule((ThrottleSegment(10.0, 1.0, 1.0),))
        commands = (
            EngineTransientCommand("left", schedule),
            EngineTransientCommand("right", schedule),
        )
        scenarios = PropulsionScenarioSet(
            scenarios=(
                PropulsionScenario(
                    "nominal",
                    0.90,
                    (),
                    "event-tree fixture",
                    "both engines complete the burn",
                ),
                PropulsionScenario(
                    "left-out",
                    0.08,
                    (
                        ScenarioEvent(
                            ScenarioEventType.ENGINE_DISABLED,
                            ("left",),
                            None,
                            "reliability fixture",
                            "left engine does not ignite",
                        ),
                    ),
                    "event-tree fixture",
                    "single-engine outcome",
                ),
                PropulsionScenario(
                    "common-cause-out",
                    0.02,
                    (
                        ScenarioEvent(
                            ScenarioEventType.ENGINE_DISABLED,
                            ("left", "right"),
                            None,
                            "common-cause fixture",
                            "shared initiation failure",
                        ),
                    ),
                    "event-tree fixture",
                    "both engines unavailable",
                ),
            ),
            objective="joint propulsion consequence verification",
        )
        self.discrete = build_propulsion_scenario_ensemble(
            members,
            commands,
            scenarios,
        )
        shared_parameters = (
            uncertain_scale("thrust_scale", 0.98, 1.02),
            uncertain_scale(
                "specific_impulse_scale",
                0.99,
                1.01,
                DistributionKind.TRIANGULAR,
            ),
            uncertain_scale("time_scale", 0.99, 1.01),
        )
        self.definitions = (
            ScenarioConditionalUncertainty(
                "nominal",
                shared_parameters,
                8,
                "nominal conditional fixture",
                "nominal performance scatter",
            ),
            ScenarioConditionalUncertainty(
                "left-out",
                shared_parameters,
                8,
                "engine-out conditional fixture",
                "remaining engine performance scatter",
            ),
            ScenarioConditionalUncertainty(
                "common-cause-out",
                (),
                1,
                "deterministic zero-thrust fixture",
                "continuous scales cannot restore disabled engines",
            ),
        )

    def test_joint_weights_close_once_and_artifacts_are_reproducible(self) -> None:
        first = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=71,
            replicate_tolerance=0.20,
        )
        second = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=71,
            replicate_tolerance=0.20,
        )
        self.assertEqual(first.realization_count, 17)
        self.assertAlmostEqual(first.joint_probability_sum, 1.0)
        self.assertLessEqual(first.joint_probability_closure_error, 1e-12)
        self.assertAlmostEqual(first.zero_thrust_probability, 0.02)
        for scenario_id, expected in (
            ("nominal", 0.90),
            ("left-out", 0.08),
            ("common-cause-out", 0.02),
        ):
            actual = sum(
                item.joint_probability
                for item in first.realizations
                if item.scenario_id == scenario_id
            )
            self.assertAlmostEqual(actual, expected)
        self.assertEqual(
            tuple(item.artifact.artifact_hash for item in first.realizations),
            tuple(item.artifact.artifact_hash for item in second.realizations),
        )
        self.assertTrue(
            all(
                item.primary_seed != item.audit_seed
                for item in first.sampling_evidence
            )
        )

        revised_definitions = (
            ScenarioConditionalUncertainty(
                "nominal",
                self.definitions[0].parameters,
                8,
                "revised traceability source",
                "same mathematics with a changed evidence basis",
            ),
            self.definitions[1],
            self.definitions[2],
        )
        revised = build_hierarchical_propulsion_ensemble(
            self.discrete,
            revised_definitions,
            seed=71,
            replicate_tolerance=0.20,
        )
        first_nominal = next(
            item for item in first.realizations if item.scenario_id == "nominal"
        )
        revised_nominal = next(
            item for item in revised.realizations if item.scenario_id == "nominal"
        )
        self.assertNotEqual(
            first_nominal.artifact.artifact_hash,
            revised_nominal.artifact.artifact_hash,
        )

    def test_physical_scaling_and_scenario_contributions_close(self) -> None:
        result = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=19,
            replicate_tolerance=0.20,
        )
        nominal_base = next(
            item.artifact for item in self.discrete.realizations if item.scenario_id == "nominal"
        )
        realization = next(
            item
            for item in result.realizations
            if item.scenario_id == "nominal" and item.conditional_index == 0
        )
        inputs = dict(realization.inputs)
        self.assertAlmostEqual(
            realization.artifact.duration_s,
            nominal_base.duration_s * inputs["time_scale"],
        )
        self.assertAlmostEqual(
            realization.artifact.system_equivalent_specific_impulse_s,
            nominal_base.system_equivalent_specific_impulse_s
            * inputs["specific_impulse_scale"],
        )
        impulse = next(
            item for item in result.statistics if item.metric == "delivered_impulse_n_s"
        )
        contribution_sum = sum(
            item.joint_expected_value_contribution
            for item in result.contributions
            if item.metric == "delivered_impulse_n_s"
        )
        self.assertAlmostEqual(contribution_sum, impulse.expected_value)
        self.assertIn("NASA-STD-7009B", realization.artifact.reference_ids)

    def test_definition_coverage_and_axial_physics_fail_closed(self) -> None:
        with self.assertRaisesRegex(DomainError, "exactly one"):
            build_hierarchical_propulsion_ensemble(
                self.discrete,
                self.definitions[:-1],
            )
        invalid_axial = ScenarioConditionalUncertainty(
            "nominal",
            (uncertain_scale("axial_efficiency_scale", 0.99, 1.01),),
            8,
            "invalid axial fixture",
            "prove axial thrust cannot exceed delivered thrust",
        )
        definitions = (invalid_axial, self.definitions[1], self.definitions[2])
        with self.assertRaisesRegex(DomainError, "axial thrust exceed"):
            build_hierarchical_propulsion_ensemble(
                self.discrete,
                definitions,
            )

    def test_total_variance_closes_between_and_within_scenarios(self) -> None:
        result = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=29,
            replicate_tolerance=0.20,
        )
        decomposition = next(
            item
            for item in result.variance_decomposition
            if item.metric == "delivered_impulse_n_s"
        )
        self.assertAlmostEqual(
            decomposition.total_variance,
            decomposition.within_scenario_variance
            + decomposition.between_scenario_variance,
        )
        self.assertLessEqual(decomposition.closure_error, 1e-9)
        self.assertAlmostEqual(
            decomposition.within_fraction + decomposition.between_fraction,
            1.0,
        )
        self.assertGreater(
            decomposition.between_scenario_variance,
            decomposition.within_scenario_variance,
        )

    def test_ensemble_round_trip_and_semantic_tamper_refusal(self) -> None:
        result = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=41,
            replicate_tolerance=0.20,
        )
        text = hierarchical_ensemble_to_json(result)
        reopened = hierarchical_ensemble_from_json(text)
        self.assertEqual(reopened, result)
        self.assertEqual(
            hierarchical_ensemble_hash(reopened),
            hierarchical_ensemble_to_dict(result)["ensemble_hash"],
        )

        nested_tamper = copy.deepcopy(hierarchical_ensemble_to_dict(result))
        nested_tamper["realizations"][0]["artifact"]["segments"][0][
            "start_delivered_thrust_n"
        ] += 1.0
        with self.assertRaisesRegex(InputError, "artifact hash"):
            hierarchical_ensemble_from_dict(nested_tamper)

        semantic_tamper = copy.deepcopy(hierarchical_ensemble_to_dict(result))
        semantic_tamper["statistics"][0]["expected_value"] += 1.0
        semantic_tamper.pop("ensemble_hash")
        semantic_tamper["ensemble_hash"] = hashlib.sha256(
            json.dumps(
                semantic_tamper,
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        with self.assertRaisesRegex(InputError, "statistic"):
            hierarchical_ensemble_from_dict(semantic_tamper)


if __name__ == "__main__":
    unittest.main()
