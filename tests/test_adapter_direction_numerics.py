"""Public adapter scale/sign, representability and native-boundary attacks."""

import importlib.util
import json
import random
import unittest
from decimal import Decimal, localcontext
from itertools import permutations
from math import hypot, isfinite, ldexp, ulp
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import (
    SideraEnsembleBurnJob,
    SideraEnsembleHandoff,
    export_sidera_artifact,
    export_sidera_ensemble_handoff,
    lunaris_impulsive_equivalents,
    to_native_finite_burn,
    to_native_maneuver_plan,
)
from rocket_propulsion.propulsion.burns import (
    PropulsionTimeline,
    ScheduledBurn,
    simulate_constant_burn,
)
from tests.test_burn_constant_numerics import mono
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256
from tests.test_lunaris_integration import mass_only_bridge
from tests import test_sidera_ensemble_handoff as ensemble_fixtures


def decimal_vector(direction, magnitude, precision):
    with localcontext() as context:
        context.prec = precision
        coordinates = tuple(Decimal.from_float(float(value)) for value in direction)
        length = sum((value * value for value in coordinates), Decimal(0)).sqrt()
        return tuple(
            float(Decimal.from_float(float(magnitude)) * value / length) for value in coordinates
        )


class AdapterDirectionAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        definition, point = mono(cant=0.8)
        cls.burn = simulate_constant_burn(definition, point)
        cls.bridge = mass_only_bridge(100.0, 0.1, 1.0)
        cls.timeline = PropulsionTimeline((ScheduledBurn("burn", cls.bridge, 0.0),))
        fixture = ensemble_fixtures.SideraEnsembleHandoffTests()
        fixture.setUp()
        cls.ensemble = fixture.ensemble

    def unit_outputs(self, direction):
        exported = export_sidera_artifact(self.burn, t_start_s=0.0, direction=direction)
        job = SideraEnsembleBurnJob(
            "synthetic", 0, 1.0, 1.0, 1.0, 0.0, direction, "inertial", self.bridge.artifact
        )
        return {"constant": tuple(exported["maneuvers"][0]["direction"]), "ensemble": job.direction}

    def impulse(self, direction, bridge=None):
        timeline = (
            self.timeline
            if bridge is None
            else PropulsionTimeline((ScheduledBurn("burn", bridge, 0.0),))
        )
        return lunaris_impulsive_equivalents(
            timeline, initial_mass_kg=1000.0, directions={"burn": (direction, "inertial")}
        )[0]

    def close_vector(self, actual, expected):
        for value, reference in zip(actual, expected, strict=True):
            self.assertTrue(isfinite(value))
            self.assertLessEqual(abs(value - reference), 8 * ulp(reference), (actual, expected))

    def reference(self, direction, magnitude=1.0):
        value = decimal_vector(direction, magnitude, 300)
        self.assertEqual(value, decimal_vector(direction, magnitude, 200))
        return value

    def test_huge_and_subnormal_axes_preserve_signed_units(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        for scale in (float_info.max, 1e-300, tiny):
            for axis in range(3):
                for sign in (-1.0, 1.0):
                    direction = tuple(sign * scale if index == axis else 0.0 for index in range(3))
                    expected = tuple(sign if index == axis else 0.0 for index in range(3))
                    for name, value in self.unit_outputs(direction).items():
                        with self.subTest(scale=scale, axis=axis, sign=sign, adapter=name):
                            self.assertEqual(value, expected)

    def test_unrepresentable_raw_norm_keeps_finite_signed_components(self):
        direction = (float_info.max, -float_info.max, float_info.max)
        for name, value in self.unit_outputs(direction).items():
            with self.subTest(adapter=name):
                self.close_vector(value, self.reference(direction))
                self.assertLessEqual(abs(hypot(*value) - 1.0), 2 * ulp(1.0))

    def test_power_of_two_scaling_and_permutation_preserve_direction(self):
        for base in permutations((1.0, 2.0, -2.0)):
            expected = self.reference(base)
            for exponent in (-1000, 0, 1000):
                direction = tuple(ldexp(value, exponent) for value in base)
                for name, value in self.unit_outputs(direction).items():
                    with self.subTest(base=base, exponent=exponent, adapter=name):
                        self.close_vector(value, expected)

    def test_seeded_anisotropic_vectors_match_independent_decimal_paths(self):
        random_source = random.Random(83017)
        for _ in range(32):
            common = random_source.randint(-700, 700)
            direction = tuple(
                ldexp(
                    random_source.choice((-1.0, 1.0)) * random_source.uniform(0.5, 1.0),
                    common - random_source.randint(0, 300),
                )
                for _ in range(3)
            )
            for name, value in self.unit_outputs(direction).items():
                with self.subTest(direction=direction, adapter=name):
                    self.close_vector(value, self.reference(direction))

    def test_unit_component_range_loss_is_explicit_in_artifact_adapters(self):
        direction = (float_info.max, float.fromhex("0x0.0000000000001p-1022"), 0.0)
        for call in (
            lambda: export_sidera_artifact(self.burn, t_start_s=0.0, direction=direction),
            lambda: SideraEnsembleBurnJob(
                "s", 0, 1.0, 1.0, 1.0, 0.0, direction, "inertial", self.bridge.artifact
            ),
        ):
            with self.assertRaises(DomainError) as caught:
                call()
            self.assertIn("direction", (caught.exception.field or "").lower())

    def test_zero_malformed_and_boolean_directions_are_rejected_at_public_boundaries(self):
        for direction in (
            (0.0, 0.0, 0.0),
            None,
            (True, 0.0, 1.0),
            (1.0, 2.0),
            (1.0, float("nan"), 0.0),
            (1.0, float("inf"), 0.0),
            "123",
        ):
            for name, call in (
                (
                    "constant",
                    lambda: export_sidera_artifact(self.burn, t_start_s=0.0, direction=direction),
                ),
                (
                    "ensemble",
                    lambda: SideraEnsembleBurnJob(
                        "s", 0, 1.0, 1.0, 1.0, 0.0, direction, "inertial", self.bridge.artifact
                    ),
                ),
                ("lunaris", lambda: self.impulse(direction)),
            ):
                with (
                    self.subTest(adapter=name, direction=direction),
                    self.assertRaises(DomainError),
                ):
                    call()

    def test_public_ensemble_export_reopens_without_direction_or_hash_drift(self):
        for direction in ((-1.3, 4.2, 2.1), (1e308, -1e308, 1e308), (1e-300, 2e-300, -2e-300)):
            with self.subTest(direction=direction):
                handoff = export_sidera_ensemble_handoff(
                    self.ensemble, t_start_s=12.0, direction=direction
                )
                self.close_vector(handoff.jobs[0].direction, self.reference(direction))
                encoded = handoff.to_json()
                for _ in range(4):
                    reopened = SideraEnsembleHandoff.from_json(encoded)
                    self.assertEqual(reopened.to_json(), encoded)
                    self.assertEqual(reopened.handoff_hash, handoff.handoff_hash)
                    encoded = reopened.to_json()

    def test_constant_artifact_json_preserves_source_force_and_mass_channels(self):
        artifact = export_sidera_artifact(
            self.burn, t_start_s=50.0, direction=(1e308, -1e308, 0.0), frame="ric"
        )
        reopened = json.loads(json.dumps(artifact, allow_nan=False))
        burn = reopened["maneuvers"][0]
        self.close_vector(burn["direction"], self.reference((1.0, -1.0, 0.0)))
        self.assertEqual(reopened["manifest"]["source_result_hash"], self.burn.result_hash)
        self.assertEqual(burn["thrust_n"], 8.0)
        self.assertEqual(burn["throttle"], 1.0)
        self.assertEqual(burn["frame"], "ric")

    def test_direction_policy_is_bound_and_reader_does_not_repair_invalid_stored_units(self):
        handoff = export_sidera_ensemble_handoff(
            self.ensemble, t_start_s=0.0, direction=(1.0, 0.0, 0.0)
        )
        payload = handoff.to_dict()
        constant = export_sidera_artifact(self.burn, t_start_s=0.0, direction=(1.0, 0.0, 0.0))
        self.assertIsNotNone(payload.get("direction_numerics"))
        self.assertEqual(
            payload["direction_numerics"], constant["manifest"].get("direction_numerics")
        )
        for mode in ("omit", "replace", "nonunit", "boolean"):
            modified = json.loads(json.dumps(payload))
            if mode == "omit":
                del modified["direction_numerics"]
            elif mode == "replace":
                modified["direction_numerics"] = "unsupported"
            else:
                modified["jobs"][0]["direction"] = (
                    [2.0, 0.0, 0.0] if mode == "nonunit" else [True, 0.0, 0.0]
                )
            modified["handoff_hash"] = canonical_sha256(
                {key: value for key, value in modified.items() if key != "handoff_hash"}
            )
            with self.subTest(mode=mode), self.assertRaises(InputError):
                SideraEnsembleHandoff.from_dict(modified)

    def test_lunaris_true_component_range_loss_is_a_typed_field_error(self):
        direction = (float_info.max, float.fromhex("0x0.0000000000001p-1022"), 0.0)
        with self.assertRaises(DomainError) as caught:
            self.impulse(direction)
        self.assertIn("Lunaris impulse delta-v", caught.exception.field or "")

    def test_lunaris_delta_v_projection_is_scale_and_sign_invariant(self):
        for direction in ((1e308, -1e308, 1e308), (1e-300, 2e-300, -2e-300), (-3.0, 4.0, 0.0)):
            with self.subTest(direction=direction):
                result = self.impulse(direction)
                self.close_vector(result.dv_mps, self.reference(direction, result.delta_v_m_s))
                self.assertLessEqual(
                    abs(hypot(*result.dv_mps) - result.delta_v_m_s), 8 * ulp(result.delta_v_m_s)
                )

    def test_lunaris_projection_retains_small_component_rescued_by_delta_v(self):
        bridge = mass_only_bridge(1e300, 0.1, 1.0)
        direction = (1e300, 1e-200, 0.0)
        result = self.impulse(direction, bridge)
        expected = self.reference(direction, result.delta_v_m_s)
        self.assertGreater(expected[1], 0.0)
        self.close_vector(result.dv_mps, expected)

    @unittest.skipUnless(
        importlib.util.find_spec("sidera") is not None,
        "Explicit pinned Sidera source path required",
    )
    def test_native_sidera_acceleration_and_plan_keep_direction_in_supported_frames(self):
        import numpy as np
        from sidera.core.propagation.plans import thrust_acceleration_and_mass_rate

        state = np.array((7e6, 0.0, 0.0, 0.0, 1.0, 0.0, 100.0))
        direction = (1e308, -1e308, 1e308)
        unit = self.reference(direction)
        for frame, expected in (
            ("inertial", unit),
            ("ric", unit),
            ("vnb", (unit[2], unit[0], unit[1])),
        ):
            with self.subTest(frame=frame):
                native = to_native_finite_burn(
                    self.burn, t_start_s=0.0, direction=direction, frame=frame
                )
                acceleration, rate = thrust_acceleration_and_mass_rate(state, native)
                self.close_vector(acceleration, tuple(0.08 * value for value in expected))
                self.assertAlmostEqual(rate, -1.0, places=14)
                plan = to_native_maneuver_plan(
                    self.burn, t_start_s=0.0, direction=direction, frame=frame
                )
                self.assertTrue(plan.has_finite_burns)

    @unittest.skipUnless(
        importlib.util.find_spec("lunaris") is not None,
        "Explicit pinned Lunaris source path required",
    )
    def test_native_lunaris_impulse_preserves_delta_v_and_mass_law(self):
        import numpy as np
        from lunaris.core.propagation.plans import apply_impulsive_maneuver

        result = self.impulse((1e308, -1e308, 1e308))
        native = result.to_lunaris()
        state = np.array((7e6, 0.0, 0.0, 0.0, 1.0, 0.0, 1000.0))
        final, record = apply_impulsive_maneuver(state, native)
        expected = self.reference((1.0, -1.0, 1.0), result.delta_v_m_s)
        self.close_vector(
            final[3:6], tuple(old + delta for old, delta in zip(state[3:6], expected, strict=True))
        )
        self.assertLessEqual(abs(final[6] - result.final_mass_kg), 8 * ulp(result.final_mass_kg))
        self.assertLessEqual(
            abs(record["dv_norm_mps"] - result.delta_v_m_s), 8 * ulp(result.delta_v_m_s)
        )
