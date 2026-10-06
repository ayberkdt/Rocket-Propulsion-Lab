"""Checks for probability-labelled discrete propulsion scenarios."""

from __future__ import annotations

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    ConstantPerformanceProvider,
    EngineClusterMember,
    EngineTransientCommand,
    PropulsionScenario,
    PropulsionScenarioSet,
    ScenarioEvent,
    ScenarioEventType,
    ThrottleSchedule,
    ThrottleSegment,
    build_propulsion_scenario_ensemble,
)


def event(
    event_type: ScenarioEventType,
    engine_ids: tuple[str, ...],
    value: float | None = None,
) -> ScenarioEvent:
    return ScenarioEvent(
        event_type=event_type,
        engine_ids=engine_ids,
        value=value,
        source="reviewed reliability fixture",
        rationale="exercise an explicit discrete propulsion outcome",
    )


def scenario(
    scenario_id: str,
    probability: float,
    events: tuple[ScenarioEvent, ...] = (),
) -> PropulsionScenario:
    return PropulsionScenario(
        scenario_id=scenario_id,
        probability=probability,
        events=events,
        source="reviewed event-tree fixture",
        rationale="mutually exclusive software verification scenario",
    )


class PropulsionScenarioTests(unittest.TestCase):
    def setUp(self) -> None:
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=200.0,
            chamber_pressure_pa=1_000_000.0,
            tank_id="shared-propellant",
            source="scenario test rated point",
        ).operating_point()
        self.members = (
            EngineClusterMember("left", point),
            EngineClusterMember("right", point),
        )
        schedule = ThrottleSchedule(
            (ThrottleSegment(10.0, 1.0, 1.0, "steady"),),
            name="ten-second firing",
        )
        self.commands = (
            EngineTransientCommand("left", schedule),
            EngineTransientCommand("right", schedule),
        )

    def test_engine_out_and_common_cause_scenarios_keep_exact_probabilities(self) -> None:
        scenario_set = PropulsionScenarioSet(
            scenarios=(
                scenario("nominal", 0.90),
                scenario(
                    "left-engine-out",
                    0.08,
                    (event(ScenarioEventType.ENGINE_DISABLED, ("left",)),),
                ),
                scenario(
                    "common-cause-no-ignition",
                    0.02,
                    (
                        event(
                            ScenarioEventType.ENGINE_DISABLED,
                            ("left", "right"),
                        ),
                    ),
                ),
            ),
            objective="quantify delivered-impulse consequences",
        )
        result = build_propulsion_scenario_ensemble(
            self.members,
            self.commands,
            scenario_set,
        )
        impulse = next(
            item for item in result.statistics if item.metric == "delivered_impulse_n_s"
        )
        self.assertAlmostEqual(impulse.expected_value, 1880.0)
        self.assertAlmostEqual(impulse.p05, 1000.0)
        self.assertAlmostEqual(impulse.p50, 2000.0)
        self.assertAlmostEqual(result.zero_thrust_probability, 0.02)
        self.assertEqual(result.scenarios, scenario_set.scenarios)
        zero = next(
            item for item in result.realizations if item.scenario_id == "common-cause-no-ignition"
        )
        self.assertIsNone(zero.evidence)
        self.assertEqual(zero.artifact.delivered_total_impulse_n_s, 0.0)
        self.assertIn("NASA-20160007073", zero.artifact.reference_ids)

    def test_cutoff_delay_and_performance_degradation_modify_physical_history(self) -> None:
        degraded_events = (
            event(ScenarioEventType.EARLY_CUTOFF, ("left",), 4.0),
            event(ScenarioEventType.IGNITION_DELAY, ("right",), 2.0),
            event(ScenarioEventType.THRUST_SCALE, ("right",), 0.8),
            event(ScenarioEventType.SPECIFIC_IMPULSE_SCALE, ("right",), 0.95),
            event(ScenarioEventType.RESPONSE_TIME_SCALE, ("right",), 2.0),
        )
        result = build_propulsion_scenario_ensemble(
            self.members,
            self.commands,
            PropulsionScenarioSet(
                scenarios=(
                    scenario("nominal", 0.5),
                    scenario("degraded", 0.5, degraded_events),
                ),
                objective="compare execution degradations",
            ),
        )
        nominal = next(item for item in result.realizations if item.scenario_id == "nominal")
        degraded = next(item for item in result.realizations if item.scenario_id == "degraded")
        self.assertLess(
            degraded.artifact.delivered_total_impulse_n_s,
            nominal.artifact.delivered_total_impulse_n_s,
        )
        self.assertGreater(degraded.artifact.duration_s, nominal.artifact.duration_s)
        self.assertLess(
            degraded.artifact.system_equivalent_specific_impulse_s,
            nominal.artifact.system_equivalent_specific_impulse_s,
        )

    def test_incomplete_probability_and_conflicting_events_fail_closed(self) -> None:
        with self.assertRaisesRegex(DomainError, "sum to one"):
            PropulsionScenarioSet(
                scenarios=(scenario("nominal", 0.9),),
                objective="invalid incomplete event tree",
            )
        conflicting = PropulsionScenarioSet(
            scenarios=(
                scenario(
                    "conflict",
                    1.0,
                    (
                        event(ScenarioEventType.ENGINE_DISABLED, ("left",)),
                        event(ScenarioEventType.THRUST_SCALE, ("left",), 0.5),
                    ),
                ),
            ),
            objective="prove ambiguous outcomes are refused",
        )
        with self.assertRaisesRegex(DomainError, "disabled engine"):
            build_propulsion_scenario_ensemble(
                self.members,
                self.commands,
                conflicting,
            )


if __name__ == "__main__":
    unittest.main()
