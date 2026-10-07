"""Attack every persisted-document reader with mutated and hostile JSON.

Rule under test: a reader of a hash-protected document either refuses it
with ``DomainError`` / ``InputError`` or returns a model whose
re-serialization equals the original.  Anything else is a defect: a crash
with another exception type, or (worse) a changed model that slipped past the
content hash.  The study workspace is a user-editable file without a hash: a
changed document may give a changed workspace, but it must never crash and
what it accepts must re-serialize stably.

Mutations are generated from each valid document: every leaf replaced by
values of the wrong type or range, every key deleted, extra keys added, lists
emptied or truncated, NaN/Infinity tokens and over-long integer literals.
Raw hostile texts cover deep nesting, BOM, truncation and non-text input.

Sources for the threat model:

* RFC 8259, The JavaScript Object Notation (JSON) Data Interchange Format
  (no NaN/Infinity; duplicate names make behaviour unpredictable):
  https://www.rfc-editor.org/rfc/rfc8259
* Python ``json`` and the 4300-digit integer-conversion limit
  (CVE-2020-10735):
  https://docs.python.org/3/library/stdtypes.html#integer-string-conversion-length-limitation
* Content hashes of canonical JSON follow the repository's
  ``serialization.canonical_sha256``.
"""

import json
import random
import unittest
from copy import deepcopy

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import SideraEnsembleHandoff
from rocket_propulsion.propulsion.burns import (
    BurnTarget,
    BurnTargetKind,
    ElectricThrottleTable,
    HotFireMultivariateCalibration,
    PowerLimitedPropulsionModel,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    PulseModeThruster,
    RectilinearPerformanceSurface,
    ScaleDiscrepancyCalibration,
    TabulatedBurnArtifact,
    ThrottleableEngineModel,
    build_hierarchical_propulsion_ensemble,
    burn_result_from_json,
    burn_result_to_json,
    calibrate_hot_fire_multivariate,
    hierarchical_ensemble_from_json,
    hierarchical_ensemble_to_json,
    simulate_constant_burn,
)
from rocket_propulsion.studies import create_workspace, workspace_from_json, workspace_to_json

ACCEPTED = (DomainError, InputError)
unhashed: set[str] = set()
MUTANTS_PER_READER = 160


def fixtures():
    """(name, valid JSON text, reader, writer) for every persisted model."""

    from rocket_propulsion.integrations import export_sidera_ensemble_handoff
    from tests import (
        test_burn_discrepancy,
        test_burn_electric,
        test_burn_hierarchical_ensemble,
        test_burn_model_persistence,
        test_burn_multivariate_calibration,
        test_burn_performance_surface_cubic,
        test_burn_pulse_mode,
        test_burn_throttleable,
        test_burn_timeline,
        test_burns,
        test_sidera_ensemble_handoff,
        test_studies_systems,
    )

    items = []

    def add(name, model, reader, writer=None):
        writer = writer or (lambda value: value.to_json())
        items.append((name, writer(model), reader, writer))

    unhashed.add("workspace")

    bridge = test_burn_model_persistence.configured_bridge()
    add("tabulated artifact", bridge.artifact, TabulatedBurnArtifact.from_json)
    add("propagator bridge", bridge, PropagatorPropulsionBridge.from_json)
    cubic = test_burn_performance_surface_cubic
    add("performance surface", cubic.surface_from(cubic.smooth, cubic.GRID),
        RectilinearPerformanceSurface.from_json)
    add("electric table", test_burn_electric.synthetic_table(), ElectricThrottleTable.from_json)
    add("electric model", test_burn_electric.thrust_model(),
        PowerLimitedPropulsionModel.from_json)
    add("timeline", test_burn_timeline.two_burn_timeline(), PropulsionTimeline.from_json)
    add("throttleable engine", test_burn_throttleable.engine(), ThrottleableEngineModel.from_json)
    add("pulse thruster", test_burn_pulse_mode.synthetic_thruster(), PulseModeThruster.from_json)
    result = simulate_constant_burn(
        test_burns.example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
        test_burns.example_operating_point(),
    )
    add("burn result", result, burn_result_from_json, burn_result_to_json)
    add("multivariate calibration",
        calibrate_hot_fire_multivariate(
            test_burn_multivariate_calibration.campaign(), source="adversarial",
            rationale="mutation fixture"),
        HotFireMultivariateCalibration.from_json)
    discrepancy_case = test_burn_discrepancy.DiscrepancyCalibrationTests()
    discrepancy_case.setUp()
    add("discrepancy calibration", discrepancy_case.model_form,
        ScaleDiscrepancyCalibration.from_json)
    ensemble_case = test_burn_hierarchical_ensemble.HierarchicalPropulsionEnsembleTests()
    ensemble_case.setUp()
    ensemble = build_hierarchical_propulsion_ensemble(
        ensemble_case.discrete, ensemble_case.definitions, seed=41, replicate_tolerance=0.20)
    add("hierarchical ensemble", ensemble, hierarchical_ensemble_from_json,
        hierarchical_ensemble_to_json)
    handoff_case = test_sidera_ensemble_handoff.SideraEnsembleHandoffTests()
    handoff_case.setUp()
    handoff = export_sidera_ensemble_handoff(handoff_case.ensemble, t_start_s=10.0,
                                             direction=(1.0, 0.0, 0.0))
    add("sidera handoff", handoff, SideraEnsembleHandoff.from_json)
    workspace_case = test_studies_systems.WorkspaceTests()
    workspace_case.setUp()
    add("workspace", create_workspace("Demo", (workspace_case.case,)), workspace_from_json,
        workspace_to_json)
    return items


def paths(node, prefix=()):
    yield prefix, node
    if isinstance(node, dict):
        for key, value in node.items():
            yield from paths(value, (*prefix, key))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from paths(value, (*prefix, index))


def replaced(document, path, value):
    copy = deepcopy(document)
    target = copy
    for step in path[:-1]:
        target = target[step]
    target[path[-1]] = value
    return copy


def deleted(document, path):
    copy = deepcopy(document)
    target = copy
    for step in path[:-1]:
        target = target[step]
    del target[path[-1]]
    return copy


HOSTILE_LEAVES = (None, "x", "", True, False, 0, -1, 1e308, -1e308, [], {}, [[[]]],
                  float("nan"), float("inf"), 10 ** 400)


def mutants(document, generator):
    """Yield (label, JSON text) mutants of a decoded document."""

    every = [(path, node) for path, node in paths(document) if path]
    generator.shuffle(every)
    produced = 0
    for path, node in every:
        if produced >= MUTANTS_PER_READER:
            return
        value = generator.choice(HOSTILE_LEAVES)
        if value == node and type(value) is type(node):
            continue
        # allow_nan=True emits NaN / Infinity tokens, which RFC 8259 forbids.
        yield f"{path} -> {value!r:.40}", json.dumps(replaced(document, path, value),
                                                    allow_nan=True)
        produced += 1
        if isinstance(path[-1], str):
            yield f"delete {path}", json.dumps(deleted(document, path))
            produced += 1
        if isinstance(node, list) and node:
            yield f"truncate {path}", json.dumps(replaced(document, path, node[:-1]))
            yield f"duplicate {path}", json.dumps(replaced(document, path, [*node, node[0]]))
            produced += 2
        if isinstance(node, dict):
            yield f"extra key {path}", json.dumps(replaced(document, path,
                                                           {**node, "zz_unknown": 1}))
            produced += 1


RAW_HOSTILE = (
    "", " ", "null", "[]", "{}", "0", "\"text\"", "9" * 5000, "[" * 50_000 + "]" * 50_000,
    "{\"a\": " * 20_000 + "1" + "}" * 20_000, "NaN", "{\"schema\": NaN}", "\x00",
)


class ReaderAttackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixtures = fixtures()

    def test_valid_documents_round_trip_exactly(self) -> None:
        for name, text, reader, writer in self.fixtures:
            self.assertEqual(writer(reader(text)), text, name)

    def test_mutated_documents_are_refused_or_unchanged(self) -> None:
        generator = random.Random(20261007)
        total = 0
        for name, text, reader, writer in self.fixtures:
            document = json.loads(text)
            for label, mutant in mutants(document, generator):
                total += 1
                try:
                    model = reader(mutant)
                except ACCEPTED:
                    continue
                except Exception as error:  # noqa: BLE001 - any other type is the defect
                    self.fail(f"{name}: {label}: {type(error).__name__}: {error}")
                if name in unhashed:
                    once = writer(model)
                    self.assertEqual(writer(reader(once)), once, f"{name}: {label} unstable")
                else:
                    self.assertEqual(writer(model), text, f"{name}: {label} changed the model")
        self.assertGreater(total, 1500)

    def test_hostile_raw_inputs_are_refused(self) -> None:
        for name, text, reader, _ in self.fixtures:
            for raw in (*RAW_HOSTILE, "﻿" + text, text[: len(text) // 2], text + "x"):
                with self.subTest(reader=name, raw=raw[:20]):
                    try:
                        reader(raw)
                    except ACCEPTED:
                        pass
                    except Exception as error:  # noqa: BLE001
                        self.fail(f"{name}: {raw[:20]!r}: {type(error).__name__}: {error}")
                    else:
                        self.fail(f"{name} accepted {raw[:20]!r}")
            for value in (None, 5, b"{}", text.encode(), ["{}"]):
                with self.subTest(reader=name, value=type(value).__name__), \
                        self.assertRaises(ACCEPTED):
                    reader(value)

    def test_duplicate_names_cannot_smuggle_a_value_past_the_hash(self) -> None:
        for name, text, reader, writer in self.fixtures:
            document = json.loads(text)
            for key, value in document.items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    smuggled = text.rstrip()[:-1] + f', "{key}": {json.dumps(value * 1.5 + 1)}}}'
                    with self.subTest(reader=name, key=key):
                        try:
                            model = reader(smuggled)
                        except ACCEPTED:
                            continue
                        if name not in unhashed:
                            self.assertEqual(writer(model), text)
                    break


if __name__ == "__main__":
    unittest.main()
