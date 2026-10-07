"""Fail-closed input handling and randomized invariants of the propagation API.

1. Fuzzing: malformed values (None, text, bool, containers, NaN, infinity)
   at every public propagation-facing boundary must raise DomainError or
   InputError, never a bare TypeError/ValueError/ZeroDivisionError.
2. Canonical numbers: artifacts built with integers must reopen from JSON
   (they used to fail their own hash check).
3. Property tests: seeded random profiles, tables and geometries must satisfy
   conservation, continuity and consistency invariants.
"""

import dataclasses
import json
import random
import unittest
from itertools import pairwise
from math import isclose

from test_burn_propagator_integration import STATES, EventAwareRk4

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    DiscreteThrottlePath,
    DynamicFeedLine,
    DynamicFeedLineState,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    FirstOrderGaussMarkov,
    GatesExecutionErrorModel,
    MassJettison,
    MonotoneCubicCurve,
    PerformanceAxis,
    PerformanceScalePoint,
    PowerLimitedPropulsionModel,
    PowerSwitchSurface,
    PropagatorPropulsionBridge,
    PropellantMassModel,
    PropellantStation,
    PropulsionParameterPrior,
    PropulsionTimeline,
    RectilinearPerformanceSurface,
    RegulatedFeedSystem,
    ScheduledBurn,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    ThrottleableEngineModel,
    ThrusterLayout,
    ThrusterMount,
    delta_v_log_sensitivities,
    free_space_burn_metrics,
)
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256

ACCEPTED = (DomainError, InputError)
BAD_VALUES = (
    None, "1.0", "x", float("nan"), float("inf"), -float("inf"), True, False,
    [], {}, (1.0, 2.0), object(),
)


def simple_artifact() -> TabulatedBurnArtifact:
    flows = (("ox", 0.8), ("fu", 0.2))
    return TabulatedBurnArtifact(
        segments=(TabulatedBurnSegment(0.0, 2.0, 100.0, 100.0, 90.0, 90.0, 1.0, 1.0,
                                       flows, flows, "s"),),
        source_id="fuzz",
        source_sha256="a" * 64,
        reference_ids=("NASA-20040000363",),
    )


def simple_bridge() -> PropagatorPropulsionBridge:
    return PropagatorPropulsionBridge(
        simple_artifact(),
        stream_bindings=(StreamStateBinding("ox", "o"), StreamStateBinding("fu", "f")),
        state_constraints=(AdditionalStateConstraint("o", 0.1),),
        protected_dry_mass_kg=10.0,
        nozzle_exit_area_m2=0.001,
    )


def electric_model() -> PowerLimitedPropulsionModel:
    table = ElectricThrottleTable(
        levels=(
            ElectricThrottleLevel("a", 1000.0, 0.04, 1.5e-6),
            ElectricThrottleLevel("b", 2000.0, 0.08, 2.6e-6),
        ),
        source_id="t",
        source_sha256="b" * 64,
    )
    return PowerLimitedPropulsionModel(
        DiscreteThrottlePath.from_table(table), propellant_reserve_kg=1.0
    )


class FuzzTests(unittest.TestCase):
    def assert_fails_closed(self, label: str, call) -> None:
        try:
            call()
        except ACCEPTED:
            pass
        except Exception as error:  # noqa: BLE001 - any other exception is the failure under test
            self.fail(f"{label} raised {type(error).__name__}: {error}")

    def test_bridge_evaluation_and_events(self) -> None:
        bridge = simple_bridge()
        base = {
            "relative_time_s": 1.0,
            "vehicle_mass_kg": 50.0,
            "direction": (1.0, 0.0, 0.0),
            "additional_states": {"o": 5.0, "f": 2.0},
            "ambient_pressure_pa": 1000.0,
        }
        keys = (*base, "parameter_overrides", "operating_conditions", "side",
                "apply_inventory_limits")
        extra = ((1.0, "a", 0.0), {"o": "x", "f": 1.0}, {"o": None, "f": 1.0}, {1: 2.0})
        for key in keys:
            for value in (*BAD_VALUES, *extra):
                arguments = dict(base, **{key: value})
                self.assert_fails_closed(
                    f"evaluate {key}={value!r}", lambda a=arguments: bridge.evaluate(**a)
                )
        for value in BAD_VALUES:
            self.assert_fails_closed(
                f"mass_closure {value!r}",
                lambda v=value: bridge.mass_closure(
                    vehicle_mass_kg=v, additional_states={"o": 1.0, "f": 1.0},
                    non_propellant_mass_kg=1.0,
                ),
            )
            self.assert_fails_closed(
                f"segment_index_at {value!r}",
                lambda v=value: bridge.artifact.segment_index_at(v),
            )
            self.assert_fails_closed(
                f"state_at {value!r}", lambda v=value: bridge.artifact.state_at(v)
            )
            for event in bridge.event_surfaces():
                self.assert_fails_closed(
                    f"{event.event_id} {value!r}",
                    lambda v=value, e=event: e.value(
                        relative_time_s=v, vehicle_mass_kg=1.0, additional_states={"o": v}
                    ),
                )

    def test_constructors(self) -> None:
        artifact = simple_artifact()
        bindings = (StreamStateBinding("ox", "o"), StreamStateBinding("fu", "f"))
        for value in BAD_VALUES:
            cases = {
                "segment time": lambda v=value: TabulatedBurnSegment(v, 1.0, 1, 1, 1, 1, 1, 1),
                "segment thrust": lambda v=value: TabulatedBurnSegment(0.0, 1.0, v, 1, 1, 1, 1, 1),
                "segment streams": lambda v=value: TabulatedBurnSegment(
                    0.0, 1.0, 1, 1, 1, 1, 1, 1, v, ()
                ),
                "artifact segments": lambda v=value: TabulatedBurnArtifact(
                    segments=v, source_id="s", source_sha256="a" * 64
                ),
                "bridge area": lambda v=value: PropagatorPropulsionBridge(
                    artifact, stream_bindings=bindings, nozzle_exit_area_m2=v
                ),
                "bridge dry mass": lambda v=value: PropagatorPropulsionBridge(
                    artifact, stream_bindings=bindings, protected_dry_mass_kg=v
                ),
                "bridge bindings": lambda v=value: PropagatorPropulsionBridge(
                    artifact, stream_bindings=v
                ),
                "axis": lambda v=value: PerformanceAxis("p", "Pa", (0.0, v)),
                "point": lambda v=value: PerformanceScalePoint((0.0,), v, 1.0),
                "level": lambda v=value: ElectricThrottleLevel("a", v, 0.04, 1.5e-6),
                "gates": lambda v=value: GatesExecutionErrorModel(v, 0.0, 0.0, 0.0),
                "gauss-markov": lambda v=value: FirstOrderGaussMarkov("x", 0.1, v),
                "mount": lambda v=value: ThrusterMount("a", (v, 0.0, 0.0), (0.0, 0.0, 1.0)),
                "prior": lambda v=value: PropulsionParameterPrior.from_standard_deviations(
                    {"thrust_scale": v}, source="s", rationale="r"
                ),
                "mass model": lambda v=value: PropellantMassModel(
                    v, (0.0, 0.0, 0.0), ((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)), ()
                ),
            }
            for label, call in cases.items():
                self.assert_fails_closed(f"{label} {value!r}", call)

    def test_electric_and_navigation_helpers(self) -> None:
        model = electric_model()
        base = {
            "level_index": 0,
            "available_power_w": 1500.0,
            "vehicle_mass_kg": 100.0,
            "direction": (0.0, 0.0, 1.0),
            "additional_states": {"xenon_mass_kg": 5.0},
        }
        for key in (*base, "parameter_overrides", "apply_inventory_limits"):
            for value in (*BAD_VALUES, 1.5, -2, 5):
                arguments = dict(base, **{key: value})
                self.assert_fails_closed(
                    f"electric {key}={value!r}", lambda a=arguments: model.evaluate(**a)
                )
        gates = GatesExecutionErrorModel(0.1, 0.1, 0.1, 0.1)
        process = FirstOrderGaussMarkov("x", 0.1, 10.0)
        layout = ThrusterLayout((ThrusterMount("a", (0, 0, 0), (0, 0, 1)),))
        masses = PropellantMassModel(
            1.0, (0, 0, 0), ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            (PropellantStation("o", (0, 0, 1)),),
        )
        timeline = PropulsionTimeline(burns=(ScheduledBurn("b", simple_bridge(), 0.0),))
        throttle_engine = ThrottleableEngineModel(100.0, 0.1, 0.4, 1.0, (("o", 1.0),))
        for value in BAD_VALUES:
            calls = {
                "select_level": lambda v=value: model.select_level(v),
                "select_level current": lambda v=value: model.select_level(1500.0, v),
                "switching": lambda v=value: model.switching_surfaces(v),
                "gates covariance": lambda v=value: gates.covariance((v, 0.0, 1.0)),
                "gates sample": lambda v=value: gates.sample((1.0, 0.0, 0.0), v),
                "gm propagate": lambda v=value: process.propagate(v, 0.1, 1.0),
                "gm transition": lambda v=value: process.transition(v),
                "gm path": lambda v=value: process.sample_path((0.0, v), random.Random(1)),
                "dv sensitivity": lambda v=value: delta_v_log_sensitivities(v),
                "resultant": lambda v=value: layout.resultant({"a": v}, (0, 0, 0)),
                "resultant thrusts": lambda v=value: layout.resultant(v, (0, 0, 0)),
                "mass model": lambda v=value: masses.evaluate({"o": v}),
                "table json": lambda v=value: ElectricThrottleTable.from_json(v),
                "surface json": lambda v=value: RectilinearPerformanceSurface.from_json(v),
                "artifact json": lambda v=value: TabulatedBurnArtifact.from_json(v),
                "prior dict": lambda v=value: PropulsionParameterPrior.from_dict(v),
                "allocate force": lambda v=value: layout.allocate(
                    force_n=v, torque_n_m=(0, 0, 0), center_of_mass_m=(0, 0, 0)
                ),
                "allocate max": lambda v=value: layout.allocate(
                    force_n=(0, 0, 1), torque_n_m=(0, 0, 0), center_of_mass_m=(0, 0, 0),
                    maximum_thrust_n=v,
                ),
                "bridge json": lambda v=value: PropagatorPropulsionBridge.from_json(v),
                "electric json": lambda v=value: PowerLimitedPropulsionModel.from_json(v),
                "timeline json": lambda v=value: PropulsionTimeline.from_json(v),
                "timeline evaluate": lambda v=value: timeline.evaluate(
                    time_s=v, vehicle_mass_kg=10.0, directions={"b": (1, 0, 0)},
                    additional_states={"o": 1.0, "f": 1.0},
                ),
                "timeline directions": lambda v=value: timeline.evaluate(
                    time_s=1.0, vehicle_mass_kg=10.0, directions=v,
                    additional_states={"o": 1.0, "f": 1.0},
                ),
                "timeline stopped": lambda v=value: timeline.evaluate(
                    time_s=1.0, vehicle_mass_kg=10.0, directions={"b": (1, 0, 0)},
                    additional_states={"o": 1.0, "f": 1.0}, stopped_burns=v,
                ),
                "jettison": lambda v=value: MassJettison("j", 1.0, v),
                "throttle": lambda v=value: throttle_engine.evaluate(
                    throttle=v, vehicle_mass_kg=100.0, direction=(1, 0, 0),
                    additional_states={"o": 1.0},
                ),
                "throttle states": lambda v=value: throttle_engine.evaluate(
                    throttle=0.7, vehicle_mass_kg=100.0, direction=(1, 0, 0), additional_states=v,
                ),
                "throttle fractions": lambda v=value: ThrottleableEngineModel(
                    100.0, 0.1, 0.4, 1.0, v
                ),
                "throttle json": lambda v=value: ThrottleableEngineModel.from_json(v),
                "curve": lambda v=value: MonotoneCubicCurve((0.0, 1.0), (0.0, v)),
                "curve evaluate": lambda v=value: MonotoneCubicCurve((0.0, 1.0), (0.0, 1.0))
                .evaluate(v),
                "metrics mass": lambda v=value: free_space_burn_metrics(simple_artifact(), v),
                "metrics points": lambda v=value: free_space_burn_metrics(
                    simple_artifact(), 100.0, quadrature_points=v
                ),
                "event jump": lambda v=value: simple_bridge().time_event_jump(
                    relative_time_s=v, vehicle_mass_kg=10.0, direction=(1, 0, 0),
                    additional_states={"o": 1.0, "f": 1.0}, ambient_pressure_pa=0.0,
                ),
                "scheduled": lambda v=value: ScheduledBurn("b", simple_bridge(), v),
                "table dict": lambda v=value: ElectricThrottleTable.from_dict(v),
                "artifact dict": lambda v=value: TabulatedBurnArtifact.from_dict(v),
            }
            for label, call in calls.items():
                self.assert_fails_closed(f"{label} {value!r}", call)


class FeedAndAdapterFuzzTests(unittest.TestCase):
    """Feed-system, feed-line and Lunaris-adapter boundaries fail closed too."""

    def assert_fails_closed(self, label: str, call) -> None:
        try:
            call()
        except ACCEPTED:
            pass
        except Exception as error:  # noqa: BLE001 - any other exception is the failure under test
            self.fail(f"{label} raised {type(error).__name__}: {error}")

    def test_feed_system_and_line(self) -> None:
        from tests.test_burn_feed_line import example_line
        from tests.test_burn_feed_system import example_feed_system

        system = example_feed_system()
        state = system.initial_state(
            propellant_mass_kg=10.0, tank_pressure_pa=2e6, tank_gas_temperature_k=300.0,
            bottle_pressure_pa=2e7, bottle_gas_temperature_k=300.0, regulator_opening=0.5,
        )
        line = example_line()
        line_state = DynamicFeedLineState(0.2, 1.84e6)
        initial = {
            "propellant_mass_kg": 10.0, "tank_pressure_pa": 2e6, "tank_gas_temperature_k": 300.0,
            "bottle_pressure_pa": 2e7, "bottle_gas_temperature_k": 300.0,
            "regulator_opening": 0.5,
        }
        for value in BAD_VALUES:
            self.assert_fails_closed("feed flow", lambda v=value: system.evaluate(state, v))
            self.assert_fails_closed("feed state", lambda v=value: system.evaluate(v, 0.1))
            for field in dataclasses.fields(state):
                self.assert_fails_closed(
                    f"feed state {field.name}",
                    lambda v=value, f=field: system.evaluate(
                        dataclasses.replace(state, **{f.name: v}), 0.1
                    ),
                )
            for name in initial:
                self.assert_fails_closed(
                    f"initial {name}",
                    lambda v=value, n=name: system.initial_state(**dict(initial, **{n: v})),
                )
            for field in dataclasses.fields(system.configuration):
                self.assert_fails_closed(
                    f"feed config {field.name}",
                    lambda v=value, f=field: RegulatedFeedSystem(
                        dataclasses.replace(system.configuration, **{f.name: v})
                    ),
                )
            self.assert_fails_closed(
                "line upstream",
                lambda v=value: line.evaluate(
                    line_state, upstream_pressure_pa=v, engine_mass_flow_kg_s=0.2
                ),
            )
            self.assert_fails_closed(
                "line flow",
                lambda v=value: line.evaluate(
                    line_state, upstream_pressure_pa=2e6, engine_mass_flow_kg_s=v
                ),
            )
            self.assert_fails_closed(
                "line state",
                lambda v=value: line.evaluate(
                    v, upstream_pressure_pa=2e6, engine_mass_flow_kg_s=0.2
                ),
            )
            for field in dataclasses.fields(line.configuration):
                self.assert_fails_closed(
                    f"line config {field.name}",
                    lambda v=value, f=field: DynamicFeedLine(
                        dataclasses.replace(line.configuration, **{f.name: v})
                    ),
                )

    def test_sidera_adapters(self) -> None:
        from rocket_propulsion.integrations import (
            SideraEnsembleBurnJob,
            export_sidera_artifact,
        )
        from rocket_propulsion.propulsion.burns import (
            BurnTarget,
            BurnTargetKind,
            simulate_constant_burn,
        )
        from tests.test_burns import example_definition, example_operating_point

        result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
            example_operating_point(),
        )
        base = {"t_start_s": 0.0, "direction": (1.0, 0.0, 0.0), "frame": "inertial"}
        for value in BAD_VALUES:
            self.assert_fails_closed(
                "sidera result", lambda v=value: export_sidera_artifact(v, **base)
            )
            for key in (*base, "mode"):
                self.assert_fails_closed(
                    f"sidera {key}",
                    lambda v=value, k=key: export_sidera_artifact(result, **dict(base, **{k: v})),
                )
            job = {
                "scenario_id": "s", "conditional_index": 0, "scenario_probability": 1.0,
                "conditional_probability": 1.0, "joint_probability": 1.0, "t_start_s": 0.0,
                "direction": (1.0, 0.0, 0.0), "frame": "ric", "artifact": simple_artifact(),
            }
            for key in job:
                self.assert_fails_closed(
                    f"sidera job {key}",
                    lambda v=value, k=key, j=job: SideraEnsembleBurnJob(**dict(j, **{k: v})),
                )

    def test_lunaris_adapter(self) -> None:
        from rocket_propulsion.integrations import (
            LunarisFiniteBurnDynamics,
            lunaris_impulsive_equivalents,
        )

        bridge = PropagatorPropulsionBridge(simple_artifact(), require_complete_stream_binding=False)
        timeline = PropulsionTimeline(burns=(ScheduledBurn("b", bridge, 0.0),))
        for value in BAD_VALUES:
            self.assert_fails_closed(
                "impulsive mass",
                lambda v=value: lunaris_impulsive_equivalents(
                    timeline, initial_mass_kg=v, directions={"b": ((1, 0, 0), "inertial")}
                ),
            )
            self.assert_fails_closed(
                "impulsive directions",
                lambda v=value: lunaris_impulsive_equivalents(
                    timeline, initial_mass_kg=10.0, directions=v
                ),
            )
            self.assert_fails_closed(
                "impulsive entry",
                lambda v=value: lunaris_impulsive_equivalents(
                    timeline, initial_mass_kg=10.0, directions={"b": v}
                ),
            )
            self.assert_fails_closed(
                "impulsive frame",
                lambda v=value: lunaris_impulsive_equivalents(
                    timeline, initial_mass_kg=10.0, directions={"b": ((1, 0, 0), v)}
                ),
            )
            self.assert_fails_closed(
                "finite timeline", lambda v=value: LunarisFiniteBurnDynamics(object(), v)
            )


class CanonicalNumberTests(unittest.TestCase):
    def test_integer_built_artifacts_reopen_from_json(self) -> None:
        segment = TabulatedBurnSegment(0, 2, 100, 100, 90, 90, 1, 1, (("ox", 1),), (("ox", 1),))
        artifact = TabulatedBurnArtifact(segments=(segment,), source_id="i", source_sha256="a" * 64)
        self.assertEqual(TabulatedBurnArtifact.from_json(artifact.to_json()), artifact)
        self.assertIsInstance(segment.start_time_s, float)
        table = ElectricThrottleTable(
            levels=(ElectricThrottleLevel("a", 1000, 0.04, 0.0000015),),
            source_id="t",
            source_sha256="b" * 64,
        )
        self.assertEqual(ElectricThrottleTable.from_json(table.to_json()), table)
        surface = RectilinearPerformanceSurface(
            axes=(PerformanceAxis("p", "Pa", (0, 1)),),
            points=(PerformanceScalePoint((0,), 1, 1), PerformanceScalePoint((1,), 1, 1)),
            source_id="s",
            source_sha256="c" * 64,
        )
        self.assertEqual(RectilinearPerformanceSurface.from_json(surface.to_json()), surface)
        prior = PropulsionParameterPrior(("thrust_scale",), (1,), ((0,),), "s", "r")
        reopened = PropulsionParameterPrior.from_dict(json.loads(json.dumps(prior.to_dict())))
        self.assertEqual(reopened, prior)

    def test_integer_and_float_configurations_share_one_identity(self) -> None:
        as_int = PropagatorPropulsionBridge(simple_artifact(), stream_bindings=(
            StreamStateBinding("ox", "o"), StreamStateBinding("fu", "f")),
            protected_dry_mass_kg=10, nozzle_exit_area_m2=1)
        as_float = PropagatorPropulsionBridge(simple_artifact(), stream_bindings=(
            StreamStateBinding("ox", "o"), StreamStateBinding("fu", "f")),
            protected_dry_mass_kg=10.0, nozzle_exit_area_m2=1.0)
        self.assertEqual(as_int.model_hash, as_float.model_hash)

    def test_booleans_are_not_numbers(self) -> None:
        with self.assertRaises(DomainError):
            TabulatedBurnSegment(0.0, True, 1, 1, 1, 1, 1, 1)
        with self.assertRaises(DomainError):
            electric_model().evaluate(
                level_index=True, available_power_w=1500.0, vehicle_mass_kg=10.0,
                direction=(0, 0, 1), additional_states={"xenon_mass_kg": 5.0},
            )


class ReviewRegressionTests(unittest.TestCase):
    """One test per finding of the code review of the propagation interface."""

    def test_non_string_identifiers_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            PropagatorPropulsionBridge(
                simple_artifact(),
                stream_bindings=(StreamStateBinding("ox", "o"), StreamStateBinding("fu", "f")),
                model_id=None,  # type: ignore[arg-type]
            )
        with self.assertRaises(DomainError):
            PowerSwitchSurface(None, 1.0, "increasing", 0, 1)  # type: ignore[arg-type]
        with self.assertRaises(DomainError):
            PowerSwitchSurface("a", 1.0, ["x"], 0, 1)  # type: ignore[arg-type]
        with self.assertRaises(DomainError):
            RectilinearPerformanceSurface(
                axes=(PerformanceAxis("p", "Pa", (0.0, 1.0)),),
                points=(PerformanceScalePoint((0.0,), 1, 1), PerformanceScalePoint((1.0,), 1, 1)),
                source_id="s",
                source_sha256="c" * 64,
                reference_ids=5,  # type: ignore[arg-type]
            )

    def test_nested_tamper_diagnosis_survives_bridge_reopen(self) -> None:
        document = simple_bridge().to_dict()
        document["artifact"]["segments"][0]["end_axial_thrust_n"] = 91.0
        with self.assertRaises(InputError) as caught:
            PropagatorPropulsionBridge.from_dict(document)
        self.assertIn("hash does not match", str(caught.exception))

    def test_readers_do_not_coerce_text_to_numbers(self) -> None:
        prior = PropulsionParameterPrior.from_standard_deviations(
            {"thrust_scale": 0.01}, source="s", rationale="r"
        )
        payload = prior.to_dict(include_hash=False)
        payload["mean"] = ["1.0"]
        payload["prior_hash"] = canonical_sha256(payload)
        with self.assertRaises(InputError):
            PropulsionParameterPrior.from_dict(payload)
        artifact = simple_artifact().to_dict(include_hash=False)
        artifact["segments"][0]["start_time_s"] = "0.0"
        artifact["artifact_hash"] = canonical_sha256(artifact)
        with self.assertRaises(InputError):
            TabulatedBurnArtifact.from_dict(artifact)

    def test_timeline_event_namespaces_cannot_collide(self) -> None:
        with self.assertRaises(DomainError):
            ScheduledBurn("jettison", simple_bridge(), 0.0)
        with self.assertRaises(DomainError):
            MassJettison("a/b", 1.0, 1.0)

    def test_timeline_validates_timing_driver_even_when_burn_looks_inactive(self) -> None:
        timeline = PropulsionTimeline(burns=(ScheduledBurn("b", simple_bridge(), 0.0),))
        with self.assertRaises(DomainError):
            timeline.evaluate(
                time_s=100.0, vehicle_mass_kg=10.0, directions={},
                additional_states={"o": 1.0, "f": 1.0},
                parameter_overrides={"b": {"ignition_time_bias_s": 1e9}},
            )

    def test_burn_metrics_error_estimate_is_never_a_fake_zero(self) -> None:
        with self.assertRaises(DomainError):
            free_space_burn_metrics(simple_artifact(), 100.0, quadrature_points=101)
        metrics = free_space_burn_metrics(simple_artifact(), 100.0, quadrature_points=100)
        self.assertEqual(metrics.quadrature_points, 200)

    def test_net_thrust_margin_is_none_only_when_not_firing(self) -> None:
        bridge = simple_bridge()
        common = {"vehicle_mass_kg": 50.0, "direction": (1.0, 0.0, 0.0),
                  "additional_states": {"o": 5.0, "f": 2.0}, "ambient_pressure_pa": 1000.0}
        self.assertIsNone(bridge.evaluate(relative_time_s=5.0, **common).net_thrust_margin_n)
        self.assertAlmostEqual(
            bridge.evaluate(relative_time_s=1.0, **common).net_thrust_margin_n, 90.0 - 1.0
        )


def random_artifact(generator: random.Random) -> TabulatedBurnArtifact:
    segment_count = generator.randint(1, 12)
    boundaries = [0.0]
    for _ in range(segment_count):
        boundaries.append(boundaries[-1] + generator.uniform(0.05, 3.0))
    exhaust_velocity = 50.0
    segments = []
    previous_end = None
    for index, (start, end) in enumerate(pairwise(boundaries)):
        # Jumps at some boundaries, continuity at others.
        start_thrust = (
            previous_end if previous_end is not None and generator.random() < 0.5
            else generator.uniform(0.0, 300.0)
        )
        end_thrust = generator.uniform(0.0, 300.0)
        split = generator.uniform(0.6, 0.9)

        def streams(thrust: float, split: float = split) -> tuple[tuple[str, float], ...]:
            flow = thrust / exhaust_velocity
            return (("ox", flow * split), ("fu", flow - flow * split))

        segments.append(
            TabulatedBurnSegment(
                start, end, start_thrust, end_thrust, start_thrust, end_thrust,
                start_thrust / exhaust_velocity, end_thrust / exhaust_velocity,
                streams(start_thrust), streams(end_thrust), f"p{index}",
            )
        )
        previous_end = end_thrust
    return TabulatedBurnArtifact(
        segments=tuple(segments), source_id="random", source_sha256="d" * 64,
        reference_ids=("NASA-20040000363",),
    )


class RandomizedBridgeInvariantTests(unittest.TestCase):
    def test_conservation_continuity_and_round_trip(self) -> None:
        generator = random.Random(20261007)
        for _ in range(25):
            artifact = random_artifact(generator)
            self.assertEqual(TabulatedBurnArtifact.from_json(artifact.to_json()), artifact)
            bridge = PropagatorPropulsionBridge(
                artifact,
                stream_bindings=(
                    StreamStateBinding("ox", STATES[0]), StreamStateBinding("fu", STATES[1])
                ),
            )
            times = [event.relative_time_s for event in bridge.event_surfaces()]
            self.assertEqual(times, sorted(times))
            self.assertEqual(len(times), len(set(times)))
            for segment in artifact.segments:
                for time_s, side, expected in (
                    (segment.start_time_s, "right", segment.start_axial_thrust_n),
                    (segment.end_time_s, "left", segment.end_axial_thrust_n),
                ):
                    value = bridge.evaluate(
                        relative_time_s=time_s, vehicle_mass_kg=500.0,
                        direction=(0.0, 1.0, 0.0),
                        additional_states={STATES[0]: 100.0, STATES[1]: 50.0}, side=side,
                    )
                    self.assertAlmostEqual(value.axial_thrust_n, expected, places=9)
                    drains = dict(value.additional_state_derivatives)
                    self.assertAlmostEqual(
                        sum(drains.values()), value.mass_derivative_kg_s, places=12
                    )
                    scaled = bridge.evaluate(
                        relative_time_s=time_s, vehicle_mass_kg=500.0,
                        direction=(0.0, 1.0, 0.0),
                        additional_states={STATES[0]: 100.0, STATES[1]: 50.0}, side=side,
                        parameter_overrides={"thrust_scale": 1.25},
                    )
                    self.assertAlmostEqual(scaled.axial_thrust_n, 1.25 * expected, places=9)

            # Event-aware RK4 consumes exactly the artifact's propellant.
            propagator = EventAwareRk4(bridge, (0.0, 1.0, 0.0))
            initial = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2000.0, 900.0, 300.0]
            final = propagator.propagate(-0.1, artifact.duration_s + 0.5, initial, 0.05)
            self.assertTrue(
                isclose(2000.0 - final[6], artifact.consumed_propellant_kg,
                        rel_tol=1e-9, abs_tol=1e-9)
            )
            closure = bridge.mass_closure(
                vehicle_mass_kg=final[6],
                additional_states=dict(zip(STATES, final[7:9], strict=True)),
                non_propellant_mass_kg=800.0,
            )
            self.assertLess(abs(closure.residual_kg), 1e-9)


class RandomizedElectricInvariantTests(unittest.TestCase):
    def test_paths_and_level_selection(self) -> None:
        generator = random.Random(99)
        for _ in range(40):
            levels = []
            for index in range(generator.randint(1, 15)):
                power = generator.uniform(300.0, 7000.0)
                isp = generator.uniform(1200.0, 4500.0)
                efficiency = generator.uniform(0.3, 0.7)
                exhaust = 9.80665 * isp
                thrust = 2.0 * efficiency * power / exhaust
                levels.append(ElectricThrottleLevel(f"L{index}", power, thrust, thrust / exhaust))
            table = ElectricThrottleTable(
                levels=tuple(levels), source_id="r", source_sha256="e" * 64
            )
            for objective in ("thrust", "specific_impulse"):
                path = DiscreteThrottlePath.from_table(table, objective=objective)
                path_levels = path.levels
                powers = [level.input_power_w for level in path_levels]
                scores = [
                    level.thrust_n if objective == "thrust" else level.specific_impulse_s
                    for level in path_levels
                ]
                self.assertEqual(powers, sorted(set(powers)))
                self.assertEqual(scores, sorted(set(scores)))
                model = PowerLimitedPropulsionModel(path)
                previous = -1
                for power in sorted(generator.uniform(0.0, 8000.0) for _ in range(30)):
                    selected = model.select_level(power)
                    self.assertGreaterEqual(selected, previous)  # monotonic in power
                    previous = selected
                    if selected >= 0:
                        self.assertLessEqual(path_levels[selected].input_power_w, power + 1e-9)
                    if selected + 1 < len(path_levels):
                        self.assertGreater(path_levels[selected + 1].input_power_w, power)
                    for current in range(-1, len(path_levels)):
                        # Without hysteresis the mode is memoryless.
                        self.assertEqual(model.select_level(power, current), selected)


class RandomizedGeometryInvariantTests(unittest.TestCase):
    def test_torque_is_invariant_for_zero_net_force_and_mass_moments_close(self) -> None:
        generator = random.Random(5)
        for _ in range(30):
            position = tuple(generator.uniform(-1, 1) for _ in range(3))
            axis = tuple(generator.uniform(-1, 1) for _ in range(3))
            layout = ThrusterLayout(
                (
                    ThrusterMount("a", position, axis),
                    ThrusterMount(
                        "b", tuple(-value for value in position), tuple(-v for v in axis)
                    ),
                )
            )
            thrusts = {"a": 3.0, "b": 3.0}
            first = layout.resultant(thrusts, (0.0, 0.0, 0.0))
            second = layout.resultant(thrusts, tuple(generator.uniform(-2, 2) for _ in range(3)))
            for a, b in zip(first.torque_n_m, second.torque_n_m, strict=True):
                self.assertAlmostEqual(a, b, places=10)  # pure couple
            masses = {"o": generator.uniform(0, 50), "f": generator.uniform(0, 50)}
            model = PropellantMassModel(
                generator.uniform(10, 100),
                tuple(generator.uniform(-1, 1) for _ in range(3)),
                ((5.0, 0.0, 0.0), (0.0, 6.0, 0.0), (0.0, 0.0, 7.0)),
                (
                    PropellantStation("o", tuple(generator.uniform(-1, 1) for _ in range(3))),
                    PropellantStation("f", tuple(generator.uniform(-1, 1) for _ in range(3))),
                ),
            )
            value = model.evaluate(masses)
            moment = [model.dry_mass_kg * c for c in model.dry_center_of_mass_m]
            for station in model.stations:
                for index in range(3):
                    moment[index] += masses[station.state_name] * station.position_m[index]
            for index in range(3):
                self.assertAlmostEqual(
                    value.center_of_mass_m[index] * value.mass_kg, moment[index], places=9
                )
            inertia = value.inertia_kg_m2
            for i in range(3):
                self.assertAlmostEqual(inertia[i][(i + 1) % 3], inertia[(i + 1) % 3][i])
                self.assertGreaterEqual(
                    inertia[i][i] + 1e-9, 0.0
                )


if __name__ == "__main__":
    unittest.main()
