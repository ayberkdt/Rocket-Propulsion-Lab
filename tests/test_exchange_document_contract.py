"""Adversarial raw schema and declared-metadata contracts at public readers."""

import copy
import json
import unittest

from rocket_propulsion.core.errors import InputError
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    RectilinearStreamPerformanceSurface,
    StreamPerformanceScalePoint,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
)
from tests.test_burn_model_persistence import configured_bridge, evaluate
from tests.test_burn_timeline import two_burn_timeline


def stream_surface():
    return RectilinearStreamPerformanceSurface(
        axes=(PerformanceAxis("pressure", "Pa", (1.0, 2.0)),),
        points=(
            StreamPerformanceScalePoint((1.0,), 1.0, (("ox", 1.0), ("fu", 1.0))),
            StreamPerformanceScalePoint((2.0,), 1.2, (("ox", 0.9), ("fu", 1.1))),
        ),
        source_id="synthetic-map",
        source_sha256="a" * 64,
    )


def documents():
    bridge = configured_bridge()
    return (bridge, two_burn_timeline(), bridge.performance_surface, stream_surface())


def at(document, path):
    for key in path:
        document = document[key]
    return document


class ExchangeDocumentAttacks(unittest.TestCase):
    def reject(self, model, mutate):
        raw = copy.deepcopy(model.to_dict())
        mutate(raw)
        for encoding in ("mapping", "json"):
            with self.subTest(model=type(model).__name__, encoding=encoding):
                with self.assertRaises(InputError):
                    if encoding == "mapping":
                        type(model).from_dict(raw)
                    else:
                        type(model).from_json(json.dumps(raw))

    def test_unknown_root_fields_cannot_disappear(self):
        for model in documents():
            with self.subTest(model=type(model).__name__):
                self.reject(model, lambda d: d.update(unmodeled={"units": "lbf"}))

    def test_bridge_unknown_nested_fields_cannot_disappear(self):
        model = configured_bridge()
        for path in (
            ("stream_bindings", 0),
            ("state_constraints", 0),
            ("parameters", 0),
            ("operating_condition_axes", 0),
            ("performance_surface", "axes", 0),
            ("performance_surface", "points", 0),
            ("artifact", "segments", 0),
        ):
            with self.subTest(path=path):
                self.reject(model, lambda d: at(d, path).update(unmodeled=42))

    def test_surface_unknown_axis_and_point_fields_cannot_disappear(self):
        for model in documents()[2:]:
            for path in (("axes", 0), ("points", 0)):
                with self.subTest(model=type(model).__name__, path=path):
                    self.reject(model, lambda d: at(d, path).update(unmodeled=42))

    def test_timeline_unknown_burn_jettison_and_bridge_fields_cannot_disappear(self):
        model = two_burn_timeline()
        burn = model.burns[0].burn_id
        for path in (
            ("burns", 0),
            ("jettisons", 0),
            ("bridges", burn),
            ("bridges", burn, "parameters", 0),
        ):
            with self.subTest(path=path):
                self.reject(model, lambda d: at(d, path).update(unmodeled=42))

    def test_bridge_declared_hash_links_must_match_nested_models(self):
        for key in ("artifact_hash", "performance_surface_hash", "stream_performance_surface_hash"):
            with self.subTest(key=key):
                self.reject(configured_bridge(), lambda d: d.__setitem__(key, "0" * 64))

    def test_bridge_declared_semantics_and_reference_ids_are_checked(self):
        for key in ("force_semantics", "mass_semantics", "time_semantics", "reference_ids"):
            with self.subTest(key=key):
                self.reject(
                    configured_bridge(),
                    lambda d: d.__setitem__(
                        key,
                        ["unreviewed-source"]
                        if key == "reference_ids"
                        else "unsupported-semantics",
                    ),
                )

    def test_bridge_operating_axes_cannot_claim_different_units_or_limits(self):
        for key, value in (("name", "wrong"), ("unit", "psi"), ("minimum", 0.0), ("maximum", 1.0)):
            with self.subTest(key=key):
                self.reject(
                    configured_bridge(),
                    lambda d: d["operating_condition_axes"][0].__setitem__(key, value),
                )

    def test_timeline_burn_bridge_hash_link_is_checked(self):
        self.reject(
            two_burn_timeline(), lambda d: d["burns"][0].__setitem__("bridge_model_hash", "0" * 64)
        )

    def test_timeline_bridge_registry_must_exactly_match_scheduled_burns(self):
        self.reject(
            two_burn_timeline(),
            lambda d: d["bridges"].__setitem__("orphan", {"arbitrary": "unused"}),
        )
        self.reject(
            two_burn_timeline(),
            lambda d: d["bridges"].__setitem__(
                "unused-valid-bridge", copy.deepcopy(next(iter(d["bridges"].values())))
            ),
        )
        self.reject(two_burn_timeline(), lambda d: d["bridges"].pop(d["burns"][0]["burn_id"]))

    def test_required_writer_metadata_may_not_be_removed(self):
        cases = (
            (
                configured_bridge(),
                (
                    "artifact_hash",
                    "performance_surface_hash",
                    "operating_condition_axes",
                    "force_semantics",
                    "mass_semantics",
                    "time_semantics",
                    "reference_ids",
                ),
            ),
            (two_burn_timeline(), ("timeline_hash",)),
        )
        for model, keys in cases:
            for key in keys:
                with self.subTest(key=key):
                    self.reject(model, lambda d: d.pop(key))
        self.reject(two_burn_timeline(), lambda d: d["burns"][0].pop("bridge_model_hash"))

    def test_json_array_fields_are_not_programmatic_tuple_or_mapping_shortcuts(self):
        for model, key in (
            (configured_bridge(), "stream_bindings"),
            (two_burn_timeline(), "burns"),
            (configured_bridge().performance_surface, "axes"),
            (stream_surface(), "points"),
        ):
            for replacement in ((), {}, ""):
                raw = copy.deepcopy(model.to_dict())
                raw[key] = replacement
                with self.subTest(
                    model=type(model).__name__, replacement=type(replacement).__name__
                ):
                    with self.assertRaises(InputError):
                        type(model).from_dict(raw)
        model = configured_bridge()
        raw = model.to_dict()
        raw["stream_bindings"] = tuple(raw["stream_bindings"])
        with self.assertRaises(InputError):
            type(model).from_dict(raw)

    def test_ignored_metadata_cannot_hide_invalid_json_scalar_types(self):
        for key, value in (
            ("artifact_hash", False),
            ("operating_condition_axes", ()),
            ("reference_ids", [1]),
            ("mass_semantics", float("nan")),
            ("force_semantics", "\ud800"),
        ):
            raw = configured_bridge().to_dict()
            raw[key] = value
            with self.subTest(key=key):
                with self.assertRaises(InputError):
                    PropagatorPropulsionBridge.from_dict(raw)

    def test_all_canonical_documents_roundtrip_without_identity_changes(self):
        for model in documents():
            with self.subTest(model=type(model).__name__):
                self.assertEqual(type(model).from_dict(model.to_dict()), model)
                self.assertEqual(type(model).from_json(model.to_json()), model)
        model = configured_bridge()
        self.assertEqual(
            evaluate(PropagatorPropulsionBridge.from_json(model.to_json())), evaluate(model)
        )

    def test_vacuum_and_stream_surface_documents_remain_supported(self):
        base = configured_bridge()
        for model in (
            PropagatorPropulsionBridge(base.artifact, stream_bindings=base.stream_bindings),
            PropagatorPropulsionBridge(
                base.artifact,
                stream_bindings=base.stream_bindings,
                stream_performance_surface=stream_surface(),
            ),
        ):
            with self.subTest(surface=model.stream_performance_surface is not None):
                self.assertEqual(type(model).from_json(model.to_json()), model)

    def test_legacy_absent_null_stream_fields_and_explicit_multilinear_remain_supported(self):
        model = configured_bridge()
        raw = model.to_dict()
        raw.pop("stream_performance_surface")
        raw.pop("stream_performance_surface_hash")
        self.assertEqual(type(model).from_dict(raw), model)
        surface = model.performance_surface
        raw = surface.to_dict()
        raw["interpolation"] = "multilinear"
        self.assertEqual(type(surface).from_dict(raw), surface)

    def test_valid_changed_model_with_its_own_hash_is_not_authentication_failure(self):
        # Content identity is not an origin signature or proof of physical validity.
        model = two_burn_timeline()
        raw = model.to_dict()
        other = PropulsionTimeline(
            model.burns, model.jettisons, model.allow_overlapping_burns, "another-valid-id"
        )
        self.assertNotEqual(model.timeline_hash, other.timeline_hash)
        self.assertEqual(type(other).from_dict(other.to_dict()), other)
        self.assertEqual(raw, model.to_dict())

    def test_active_stream_surface_hash_link_is_checked_and_required(self):
        base = configured_bridge()
        model = PropagatorPropulsionBridge(
            base.artifact,
            stream_bindings=base.stream_bindings,
            stream_performance_surface=stream_surface(),
        )
        self.reject(model, lambda d: d.pop("stream_performance_surface_hash"))
        self.reject(model, lambda d: d.__setitem__("stream_performance_surface_hash", "0" * 64))

    def test_missing_nested_metadata_fields_cannot_be_invented(self):
        for key in ("name", "unit", "minimum", "maximum"):
            with self.subTest(key=key):
                self.reject(
                    configured_bridge(), lambda d: d["operating_condition_axes"][0].pop(key)
                )

    def test_direct_reader_metadata_does_not_conflate_boolean_with_one(self):
        base = configured_bridge()
        model = PropagatorPropulsionBridge(
            base.artifact,
            stream_bindings=base.stream_bindings,
            stream_performance_surface=stream_surface(),
        )
        self.reject(model, lambda d: d["operating_condition_axes"][0].__setitem__("minimum", True))

    def test_unordered_object_keys_and_surface_point_order_do_not_change_identity(self):
        for model in documents():
            raw = dict(reversed(list(model.to_dict().items())))
            if "points" in raw:
                raw["points"].reverse()
            reopened = type(model).from_dict(raw)
            self.assertEqual(reopened.to_json(), model.to_json())
            if "points" in raw:
                conditions = {axis.name: (axis.minimum + axis.maximum) / 2 for axis in model.axes}
                self.assertEqual(reopened.evaluate(conditions), model.evaluate(conditions))

    def test_nested_artifact_input_errors_keep_public_type_and_field_context(self):
        for model in (configured_bridge(), two_burn_timeline()):
            for key, value in (("source_id", "\ud800"), ("summary", {"invalid": 1})):
                raw = model.to_dict()
                bridge = (
                    raw
                    if isinstance(model, PropagatorPropulsionBridge)
                    else next(iter(raw["bridges"].values()))
                )
                bridge["artifact"][key] = value
                with self.subTest(model=type(model).__name__, key=key):
                    with self.assertRaises(InputError) as caught:
                        type(model).from_dict(raw)
                    self.assertTrue(caught.exception.field)
