"""Reopenable, hash-verified documents for configured propagation models."""

import copy
import json
import unittest

from rocket_propulsion.core.errors import InputError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    DiscreteThrottlePath,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    EstimableParameter,
    FixedEfficiencyThrottle,
    PerformanceAxis,
    PerformanceScalePoint,
    PowerLimitedPropulsionModel,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def configured_bridge() -> PropagatorPropulsionBridge:
    flows = (("ox", 0.8), ("fu", 0.2))
    artifact = TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(0.0, 1.0, 0.0, 100.0, 0.0, 95.0, 0.0, 1.0,
                                 (("ox", 0.0), ("fu", 0.0)), flows, "ramp"),
            TabulatedBurnSegment(1.0, 9.0, 100.0, 100.0, 95.0, 95.0, 1.0, 1.0,
                                 flows, flows, "steady"),
        ),
        source_id="persistence",
        source_sha256="1" * 64,
        reference_ids=("NASA-20040000363",),
    )
    surface = RectilinearPerformanceSurface(
        axes=(PerformanceAxis("supply_pressure_pa", "Pa", (1.0e6, 2.0e6, 3.0e6)),),
        points=(
            PerformanceScalePoint((1.0e6,), 0.8, 0.85),
            PerformanceScalePoint((2.0e6,), 1.0, 1.0),
            PerformanceScalePoint((3.0e6,), 1.1, 1.08),
        ),
        source_id="map",
        source_sha256="2" * 64,
    )
    return PropagatorPropulsionBridge(
        artifact,
        performance_surface=surface,
        stream_bindings=(StreamStateBinding("ox", "m_ox"), StreamStateBinding("fu", "m_fu")),
        state_constraints=(AdditionalStateConstraint("m_ox", 0.5, "ox reserve"),),
        protected_dry_mass_kg=40.0,
        parameters=(
            EstimableParameter("thrust_scale", 1.0, 1e-3, 0.5, 1.5, "1"),
            EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.5, 1.5, "1"),
            EstimableParameter("ignition_time_bias_s", 0.0, 1e-3, -2.0, 2.0, "s"),
        ),
        nozzle_exit_area_m2=0.002,
        model_id="stage-2-main",
    )


def evaluate(bridge: PropagatorPropulsionBridge):
    return bridge.evaluate(
        relative_time_s=3.0,
        vehicle_mass_kg=80.0,
        direction=(0.0, 0.6, 0.8),
        additional_states={"m_ox": 10.0, "m_fu": 3.0},
        operating_conditions={"supply_pressure_pa": 2.4e6},
        ambient_pressure_pa=5_000.0,
        parameter_overrides={"thrust_scale": 1.02},
    )


class BridgePersistenceTests(unittest.TestCase):
    def test_round_trip_preserves_identity_and_behaviour(self) -> None:
        bridge = configured_bridge()
        reopened = PropagatorPropulsionBridge.from_json(bridge.to_json())
        self.assertEqual(reopened.model_hash, bridge.model_hash)
        self.assertEqual(reopened, bridge)
        self.assertEqual(evaluate(reopened), evaluate(bridge))
        self.assertEqual(reopened.event_surfaces(), bridge.event_surfaces())
        vacuum = PropagatorPropulsionBridge(
            configured_bridge().artifact,
            stream_bindings=configured_bridge().stream_bindings,
        )
        self.assertEqual(PropagatorPropulsionBridge.from_dict(vacuum.to_dict()), vacuum)

    def test_tampering_is_detected_at_every_level(self) -> None:
        document = configured_bridge().to_dict()
        mutations = {
            "parameter bound": lambda d: d["parameters"][0].__setitem__("maximum", 1.6),
            "binding": lambda d: d["stream_bindings"][0].__setitem__("state_name", "m_x"),
            "exit area": lambda d: d.__setitem__("nozzle_exit_area_m2", 0.003),
            "artifact thrust": lambda d: d["artifact"]["segments"][1].__setitem__(
                "end_axial_thrust_n", 96.0
            ),
            "surface scale": lambda d: d["performance_surface"]["points"][0].__setitem__(
                "force_scale", 0.81
            ),
            "missing field": lambda d: d.pop("parameters"),
            "model hash": lambda d: d.__setitem__("model_hash", "0" * 64),
        }
        for label, mutate in mutations.items():
            tampered = copy.deepcopy(document)
            mutate(tampered)
            with self.subTest(label), self.assertRaises(InputError):
                PropagatorPropulsionBridge.from_dict(tampered)
        for bad in (None, "[]", 3):
            with self.assertRaises(InputError):
                PropagatorPropulsionBridge.from_json(bad)  # type: ignore[arg-type]
        with self.assertRaises(InputError):
            PropagatorPropulsionBridge.from_dict([])  # type: ignore[arg-type]


def discrete_model() -> PowerLimitedPropulsionModel:
    table = ElectricThrottleTable(
        levels=(
            ElectricThrottleLevel("L1", 1000.0, 0.040, 1.5e-6),
            ElectricThrottleLevel("L2", 2000.0, 0.080, 2.6e-6),
        ),
        source_id="table",
        source_sha256="3" * 64,
    )
    return PowerLimitedPropulsionModel(
        DiscreteThrottlePath.from_table(table, up_switch_hysteresis_w=25.0),
        propellant_reserve_kg=1.0,
        protected_dry_mass_kg=300.0,
        power_margin_w=40.0,
        duty_cycle=0.9,
        model_id="hall-string-a",
    )


class ElectricPersistenceTests(unittest.TestCase):
    def test_round_trip_for_both_laws(self) -> None:
        for model in (
            discrete_model(),
            PowerLimitedPropulsionModel(FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 4000.0, 20.0)),
        ):
            reopened = PowerLimitedPropulsionModel.from_json(model.to_json())
            self.assertEqual(reopened.model_hash, model.model_hash)
            self.assertEqual(reopened, model)
            arguments = {
                "level_index": 0,
                "available_power_w": 1800.0,
                "vehicle_mass_kg": 400.0,
                "direction": (1.0, 0.0, 0.0),
                "additional_states": {"xenon_mass_kg": 20.0},
            }
            self.assertEqual(reopened.evaluate(**arguments), model.evaluate(**arguments))
            self.assertEqual(reopened.switching_surfaces(0), model.switching_surfaces(0))

    def test_tampering_is_detected(self) -> None:
        document = discrete_model().to_dict()
        mutations = {
            "table level": lambda d: d["throttle_table"]["levels"][0].__setitem__(
                "thrust_n", 0.041
            ),
            "path hysteresis": lambda d: d["throttle_law"].__setitem__(
                "up_switch_hysteresis_w", 30.0
            ),
            "margin": lambda d: d.__setitem__("power_margin_w", 0.0),
            "unknown law": lambda d: d["throttle_law"].__setitem__("law", "magic"),
            "table swap": lambda d: d["throttle_law"].__setitem__("table_hash", "0" * 64),
        }
        for label, mutate in mutations.items():
            tampered = json.loads(json.dumps(document))
            mutate(tampered)
            with self.subTest(label), self.assertRaises(InputError):
                PowerLimitedPropulsionModel.from_dict(tampered)


if __name__ == "__main__":
    unittest.main()
