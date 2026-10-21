"""Public counterexamples for ensemble schema, identity and redundant evidence."""

from copy import deepcopy
from dataclasses import replace
import hashlib
import json
from math import nextafter
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import (
    SideraEnsembleHandoff,
    export_sidera_ensemble_handoff,
)
from rocket_propulsion.propulsion.burns import (
    CorrelationModel,
    build_hierarchical_propulsion_ensemble,
    hierarchical_ensemble_from_dict,
    hierarchical_ensemble_from_json,
    hierarchical_ensemble_to_dict,
    hierarchical_ensemble_to_json,
)
from tests import test_burn_hierarchical_ensemble as ensemble_fixture


def at(document, path):
    for part in path:
        document = document[part]
    return document


def changed(document, path, value):
    result = deepcopy(document)
    at(result, path[:-1])[path[-1]] = value
    return result


def rehash(document, field="ensemble_hash"):
    result = deepcopy(document)
    result.pop(field)
    result[field] = hashlib.sha256(
        json.dumps(
            result,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode()
    ).hexdigest()
    return result


class EnsembleDocumentAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = ensemble_fixture.HierarchicalPropulsionEnsembleTests()
        fixture.setUp()
        correlation = CorrelationModel(
            ("thrust_scale", "specific_impulse_scale", "time_scale"),
            ((1.0, 0.2, 0.0), (0.2, 1.0, 0.1), (0.0, 0.1, 1.0)),
            "synthetic correlation",
            "exercise the full persisted structure",
        )
        cls.definitions = tuple(
            replace(d, sample_count=4, correlation=correlation) if d.parameters else d
            for d in fixture.definitions
        )
        cls.discrete = fixture.discrete
        cls.ensemble = build_hierarchical_propulsion_ensemble(
            cls.discrete,
            cls.definitions,
            seed=71,
            replicate_tolerance=0.2,
        )
        cls.raw = hierarchical_ensemble_to_dict(cls.ensemble)
        cls.handoff = export_sidera_ensemble_handoff(
            cls.ensemble,
            t_start_s=10.0,
            direction=(1.0, 0.0, 0.0),
        )
        cls.hand = cls.handoff.to_dict()

    def refuse(self, document, reader=None):
        if reader is None:
            reader = hierarchical_ensemble_from_dict
        with self.assertRaises(InputError):
            reader(document)

    def test_unknown_fields_cannot_disappear_at_any_owned_object_depth(self):
        paths = [
            (),
            ("definitions", 0),
            ("definitions", 0, "parameters", 0),
            ("definitions", 0, "correlation"),
            ("realizations", 0),
            ("statistics", 0),
            ("variance_decomposition", 0),
            ("contributions", 0),
            ("sampling_evidence", 0),
        ]
        for path in paths:
            raw = deepcopy(self.raw)
            at(raw, path)["unmodelled_acceptance_claim"] = True
            with self.subTest(path=path):
                self.refuse(raw)

    def test_count_and_seed_fields_refuse_fractional_or_text_coercion(self):
        paths = [(key,) for key in ("base_seed", "scenario_count", "realization_count")]
        paths += [("definitions", 0, "sample_count"), ("realizations", 0, "conditional_index")]
        paths += [
            ("sampling_evidence", 0, key) for key in ("sample_count", "primary_seed", "audit_seed")
        ]
        for path in paths:
            value = at(self.raw, path)
            bad_values = [str(value), float(value)]
            if abs(value) < 2**52:
                bad_values.append(value + 0.25)
            for bad in bad_values:
                with self.subTest(path=path, bad=bad):
                    self.refuse(changed(self.raw, path, bad))

    def test_real_fields_refuse_numeric_text_and_booleans(self):
        paths = [("definitions", 0, "parameters", 0, key) for key in ("nominal", "lower", "upper")]
        paths += [
            ("definitions", 0, "correlation", "matrix", 0, 0),
            ("realizations", 0, "scenario_probability"),
            ("realizations", 0, "conditional_probability"),
            ("realizations", 0, "joint_probability"),
            ("realizations", 0, "inputs", 0, 1),
            ("statistics", 0, "expected_value"),
            ("variance_decomposition", 0, "total_variance"),
            ("contributions", 0, "conditional_expected_value"),
            ("sampling_evidence", 0, "requested_tolerance"),
            ("joint_probability_sum",),
        ]
        for path in paths:
            for bad in (str(at(self.raw, path)), True):
                with self.subTest(path=path, bad=bad):
                    self.refuse(changed(self.raw, path, bad))

    def test_arrays_and_matrix_rows_cannot_be_programmatic_tuples(self):
        paths = [
            (key,)
            for key in (
                "definitions",
                "realizations",
                "statistics",
                "variance_decomposition",
                "contributions",
                "sampling_evidence",
                "reference_ids",
                "warnings",
            )
        ]
        paths += [
            ("definitions", 0, "parameters"),
            ("definitions", 0, "correlation", "parameter_names"),
            ("definitions", 0, "correlation", "matrix"),
            ("definitions", 0, "correlation", "matrix", 0),
            ("realizations", 0, "inputs"),
        ]
        for path in paths:
            with self.subTest(path=path):
                self.refuse(changed(self.raw, path, tuple(at(self.raw, path))))

    def test_input_pairs_cannot_hide_third_components_or_mapping_values(self):
        pair = self.raw["realizations"][0]["inputs"][0]
        for bad in ((*pair,), [*pair, "discarded"], {pair[0]: "ignored", str(pair[1]): 999}):
            with self.subTest(bad=bad):
                self.refuse(changed(self.raw, ("realizations", 0, "inputs", 0), bad))

    def test_raw_literal_subclasses_cannot_supply_schema_equality_callbacks(self):
        class Forged(str):
            def __eq__(self, other):
                return True

        self.refuse(changed(self.raw, ("schema",), Forged("unsupported")))
        self.refuse(
            changed(self.hand, ("schema",), Forged("unsupported")), SideraEnsembleHandoff.from_dict
        )

    def test_duplicate_metric_contribution_and_evidence_identities_are_rejected(self):
        for key in ("statistics", "variance_decomposition", "contributions", "sampling_evidence"):
            raw = deepcopy(self.raw)
            raw[key].append(deepcopy(raw[key][0]))
            with self.subTest(key=key):
                self.refuse(rehash(raw))

    def test_rehashed_draws_require_unique_complete_declared_inputs(self):
        pairs = self.raw["realizations"][0]["inputs"]
        variants = [pairs[:-1], [*pairs, pairs[0]], [["unknown_scale", pairs[0][1]], *pairs[1:]]]
        for bad in variants:
            with self.subTest(bad=bad):
                self.refuse(rehash(changed(self.raw, ("realizations", 0, "inputs"), bad)))

    def test_rehashed_draw_value_cannot_escape_declared_bounds(self):
        upper = self.raw["definitions"][0]["parameters"][0]["upper"]
        raw = changed(self.raw, ("realizations", 0, "inputs", 0, 1), nextafter(upper, float("inf")))
        self.refuse(rehash(raw))

    def test_rehashed_primary_audit_and_base_seeds_must_remain_linked(self):
        for path in (
            ("base_seed",),
            ("sampling_evidence", 0, "primary_seed"),
            ("sampling_evidence", 0, "audit_seed"),
        ):
            with self.subTest(path=path):
                self.refuse(rehash(changed(self.raw, path, at(self.raw, path) + 1)))

    def test_equal_weight_method_cannot_hide_unequal_zero_thrust_samples(self):
        # Identical zero-thrust outputs make recomputed moments alone unable
        # to expose a false conditional-weight claim.
        definitions = (
            *self.definitions[:-1],
            replace(
                self.definitions[-1],
                parameters=(self.definitions[0].parameters[0],),
                sample_count=4,
            ),
        )
        ensemble = build_hierarchical_propulsion_ensemble(
            self.discrete,
            definitions,
            seed=71,
            replicate_tolerance=0.2,
        )
        raw = hierarchical_ensemble_to_dict(ensemble)
        self.assertEqual(hierarchical_ensemble_to_dict(hierarchical_ensemble_from_dict(raw)), raw)
        zero = [item for item in raw["realizations"] if item["scenario_id"] == "common-cause-out"]
        self.assertEqual(len(zero), 4)
        for item, probability in zip(zero, (0.5, 0.125, 0.125, 0.25)):
            item["conditional_probability"] = probability
            item["joint_probability"] = item["scenario_probability"] * probability
        self.refuse(rehash(raw))

    def test_variance_and_probability_domains_do_not_use_closeness_as_permission(self):
        deterministic = tuple(
            replace(d, parameters=(), sample_count=1, correlation=None) for d in self.definitions
        )
        raw = hierarchical_ensemble_to_dict(
            build_hierarchical_propulsion_ensemble(
                self.discrete,
                deterministic,
                seed=71,
                replicate_tolerance=0.2,
            )
        )
        self.assertEqual(hierarchical_ensemble_to_dict(hierarchical_ensemble_from_dict(raw)), raw)
        zero = next(i for i, row in enumerate(raw["statistics"]) if row["metric"] == "duration_s")
        self.assertEqual(raw["statistics"][zero]["standard_deviation"], 0.0)
        raw = changed(raw, ("statistics", zero, "standard_deviation"), -1e-15)
        self.refuse(rehash(raw))
        raw = changed(self.raw, ("joint_probability_closure_error",), -1e-15)
        self.refuse(rehash(raw))

    def test_required_fields_are_not_optional_while_explicit_nullable_legacy_fields_are(self):
        for path in (
            ("method",),
            ("definitions", 0, "source"),
            ("sampling_evidence", 0, "audit_seed"),
        ):
            raw = deepcopy(self.raw)
            del at(raw, path[:-1])[path[-1]]
            with self.subTest(path=path):
                self.refuse(raw)
        raw = deepcopy(self.raw)
        for definition in raw["definitions"]:
            for parameter in definition["parameters"]:
                del parameter["standard_deviation"]
            if definition["correlation"] is None:
                del definition["correlation"]
        self.assertEqual(
            hierarchical_ensemble_to_dict(hierarchical_ensemble_from_dict(raw)), self.raw
        )

    def test_canonical_correlated_json_reopens_and_valid_new_seed_is_not_frozen(self):
        self.assertEqual(
            hierarchical_ensemble_from_json(hierarchical_ensemble_to_json(self.ensemble)),
            self.ensemble,
        )
        other = build_hierarchical_propulsion_ensemble(
            self.discrete,
            self.definitions,
            seed=2**64 + 17,
            method="monte-carlo",
            replicate_tolerance=0.2,
        )
        self.assertEqual(
            hierarchical_ensemble_from_dict(hierarchical_ensemble_to_dict(other)), other
        )

    def test_single_scenario_writer_normalizes_zero_probability_without_hash_drift(self):
        discrete = replace(
            self.discrete, realizations=(replace(self.discrete.realizations[0], probability=1.0),)
        )
        ensemble = build_hierarchical_propulsion_ensemble(
            discrete,
            (self.definitions[0],),
            seed=0,
            replicate_tolerance=0.2,
        )
        self.assertEqual(ensemble.zero_thrust_probability, 0.0)
        self.assertIs(type(ensemble.zero_thrust_probability), float)
        try:
            reopened = hierarchical_ensemble_from_dict(hierarchical_ensemble_to_dict(ensemble))
        except InputError as error:
            self.fail(f"A public writer must reopen its own zero-probability field: {error}")
        self.assertEqual(reopened, ensemble)

    def test_programmatic_ensemble_records_refuse_boolean_numeric_fields(self):
        cases = [
            (self.ensemble, "base_seed"),
            (self.definitions[0], "sample_count"),
            (self.ensemble.realizations[0], "conditional_index"),
            (self.ensemble.statistics[0], "standard_deviation"),
            (self.ensemble.variance_decomposition[0], "total_variance"),
            (self.ensemble.contributions[0], "scenario_probability"),
            (self.ensemble.sampling_evidence[0], "primary_seed"),
        ]
        for record, field in cases:
            with (
                self.subTest(record=type(record).__name__, field=field),
                self.assertRaises(DomainError),
            ):
                replace(record, **{field: True})

    def test_programmatic_collection_snapshots_do_not_follow_later_caller_edits(self):
        definitions = list(self.ensemble.definitions)
        realizations = list(self.ensemble.realizations)
        ensemble = replace(self.ensemble, definitions=definitions, realizations=realizations)
        definition_inputs = list(self.definitions[0].parameters)
        definition = replace(self.definitions[0], parameters=definition_inputs)
        pairs = [list(pair) for pair in self.ensemble.realizations[0].inputs]
        realization = replace(self.ensemble.realizations[0], inputs=pairs)
        definitions.clear()
        realizations.clear()
        definition_inputs.clear()
        pairs[0][1] = 999.0
        self.assertEqual(ensemble.definitions, self.ensemble.definitions)
        self.assertEqual(ensemble.realizations, self.ensemble.realizations)
        self.assertEqual(definition.parameters, self.definitions[0].parameters)
        self.assertEqual(realization.inputs, self.ensemble.realizations[0].inputs)

    def test_builder_refuses_malformed_controls_before_sampling(self):
        for kwargs in (
            {"seed": True},
            {"seed": 0.25},
            {"seed": "71"},
            {"method": None},
            {"method": True},
            {"replicate_tolerance": True},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(DomainError):
                build_hierarchical_propulsion_ensemble(self.discrete, self.definitions, **kwargs)
        for scenarios, definitions in ((None, self.definitions), (self.discrete, (None,))):
            with self.subTest(scenarios=scenarios), self.assertRaises(DomainError):
                build_hierarchical_propulsion_ensemble(scenarios, definitions)

    def test_reordered_records_and_input_pairs_keep_their_valid_meaning(self):
        raw = deepcopy(self.raw)
        for key in (
            "definitions",
            "realizations",
            "statistics",
            "variance_decomposition",
            "contributions",
            "sampling_evidence",
        ):
            raw[key].reverse()
        for item in raw["realizations"]:
            item["inputs"].reverse()
        raw = rehash(raw)
        self.assertEqual(hierarchical_ensemble_to_dict(hierarchical_ensemble_from_dict(raw)), raw)

    def test_handoff_root_and_job_extensions_cannot_be_discarded(self):
        for path in ((), ("jobs", 0)):
            raw = deepcopy(self.hand)
            at(raw, path)["unmodelled_acceptance_claim"] = True
            with self.subTest(path=path):
                self.refuse(raw, SideraEnsembleHandoff.from_dict)

    def test_handoff_job_integer_numeric_and_collection_kinds_are_preserved(self):
        variants = [
            (("jobs", 0, "conditional_index"), 0.25),
            (
                ("jobs", 0, "scenario_probability"),
                str(self.hand["jobs"][0]["scenario_probability"]),
            ),
            (("jobs", 0, "direction"), tuple(self.hand["jobs"][0]["direction"])),
            (("jobs", 0, "direction", 0), True),
            (("jobs",), tuple(self.hand["jobs"])),
            (("required_capabilities",), tuple(self.hand["required_capabilities"])),
        ]
        for path, bad in variants:
            with self.subTest(path=path):
                self.refuse(changed(self.hand, path, bad), SideraEnsembleHandoff.from_dict)

    def test_handoff_derived_end_time_rejects_even_one_display_ulp(self):
        raw = changed(
            self.hand,
            ("jobs", 0, "t_end_s"),
            nextafter(self.hand["jobs"][0]["t_end_s"], float("inf")),
        )
        self.refuse(raw, SideraEnsembleHandoff.from_dict)

    def test_handoff_indices_cannot_have_a_rehashed_gap(self):
        raw = changed(self.hand, ("jobs", 0, "conditional_index"), 100)
        self.refuse(rehash(raw, "handoff_hash"), SideraEnsembleHandoff.from_dict)

    def test_handoff_roundtrip_reordering_and_matching_integer_real_values_remain_valid(self):
        self.assertEqual(SideraEnsembleHandoff.from_json(self.handoff.to_json()), self.handoff)
        raw = replace(self.handoff, jobs=tuple(reversed(self.handoff.jobs))).to_dict()
        raw["jobs"][0]["direction"] = [1, 0, 0]
        self.assertEqual(SideraEnsembleHandoff.from_dict(raw).to_dict(), raw)

    def test_handoff_owned_literal_identity_and_programmatic_records_are_typed(self):
        class Label(str):
            pass

        job = replace(self.handoff.jobs[0], scenario_id=Label(self.handoff.jobs[0].scenario_id))
        self.assertIs(type(job.scenario_id), str)
        hand = replace(self.handoff, schema=Label(self.handoff.schema))
        self.assertIs(type(hand.schema), str)
        for changes in ({"jobs": (None,)}, {"warnings": (True,)}, {"joint_probability_sum": True}):
            with self.subTest(changes=changes), self.assertRaises(DomainError):
                replace(self.handoff, **changes)
