"""Independent Decimal oracles for weighted moments and public rare-event weights."""

from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import ldexp, nextafter
import sys
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    ConstantPerformanceProvider,
    EngineClusterMember,
    EngineTransientCommand,
    ThrottleSchedule,
    ThrottleSegment,
    PropulsionScenario,
    PropulsionScenarioSet,
    ScenarioEvent,
    ScenarioEventType,
    ScenarioConditionalUncertainty,
    UncertainParameter,
    UncertaintyClass,
    DistributionKind,
    build_propulsion_scenario_ensemble,
    build_hierarchical_propulsion_ensemble,
    hierarchical_ensemble_from_dict,
    hierarchical_ensemble_to_dict,
)
from rocket_propulsion.propulsion.burns import ensemble as hierarchy
from rocket_propulsion.propulsion.burns import scenarios as discrete

MAX = sys.float_info.max
TINY = ldexp(1.0, -1074)


def decimal_statistics(values):
    with localcontext() as context:
        context.prec = 2200
        rows = [(Decimal.from_float(x), Decimal.from_float(w)) for x, w in values if w > 0]
        total = sum(w for x, w in rows)
        mean = sum(x * w for x, w in rows) / total
        variance = sum(w * (x - mean) ** 2 for x, w in rows) / total
        quantiles = []
        for probability in (Decimal(".05"), Decimal(".5"), Decimal(".95")):
            cumulative = Decimal(0)
            for x, w in sorted(rows):
                cumulative += w
                if cumulative >= probability * total:
                    quantiles.append(float(x))
                    break
        return (
            float(mean),
            float(variance.sqrt()),
            min(float(x) for x, w in rows),
            quantiles,
            max(float(x) for x, w in rows),
        )


def public_case(probability=0.02, samples=4, *, thrust=100.0):
    point = ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=thrust,
        system_specific_impulse_s=200.0,
        chamber_pressure_pa=1e6,
        source="synthetic numerical probability fixture",
    ).operating_point()
    members = (EngineClusterMember("a", point),)
    commands = (EngineTransientCommand("a", ThrottleSchedule((ThrottleSegment(5.0, 1.0, 1.0),))),)
    source = PropulsionScenarioSet(
        (
            PropulsionScenario("nominal", 1.0 - probability, (), "synthetic", "nominal branch"),
            PropulsionScenario(
                "rare-off",
                probability,
                (
                    ScenarioEvent(
                        ScenarioEventType.ENGINE_DISABLED, ("a",), None, "synthetic", "rare branch"
                    ),
                ),
                "synthetic",
                "positive declared branch",
            ),
        ),
        "probability-range contract",
    )
    scenarios = build_propulsion_scenario_ensemble(members, commands, source)
    parameter = UncertainParameter(
        "thrust_scale",
        1.0,
        0.98,
        1.02,
        DistributionKind.UNIFORM,
        UncertaintyClass.EPISTEMIC,
        "synthetic",
        "conditional samples",
    )
    definitions = (
        ScenarioConditionalUncertainty("nominal", (), 1, "synthetic", "deterministic"),
        ScenarioConditionalUncertainty(
            "rare-off", (parameter,) if samples > 1 else (), samples, "synthetic", "conditional"
        ),
    )
    return scenarios, definitions


class WeightedKernelAttacks(unittest.TestCase):
    def compare(self, values):
        expected = decimal_statistics(values)
        try:
            result = hierarchy._weighted_statistics("independent-kernel", values)
            scenario = discrete._metric_statistics("independent-kernel", values)
        except (DomainError, ArithmeticError) as error:
            self.fail(f"Representable moments must recover before output: {error}")
        for actual in (result, scenario):
            self.assertEqual(actual.expected_value, expected[0])
            self.assertEqual(actual.minimum, expected[2])
            self.assertEqual((actual.p05, actual.p50, actual.p95), tuple(expected[3]))
            self.assertEqual(actual.maximum, expected[4])
        self.assertEqual(result.standard_deviation, expected[1])

    def test_weighted_products_that_underflow_individually_recover_in_the_sum(self):
        self.compare(((TINY, 0.5), (TINY, 0.5)))

    def test_tiny_and_huge_dispersion_do_not_square_in_binary_float(self):
        for values in (
            ((-TINY, 0.5), (TINY, 0.5)),
            ((1e-300, 0.5), (2e-300, 0.5)),
            ((MAX, 0.5), (nextafter(MAX, 0.0), 0.5)),
            ((-MAX, 0.5), (MAX, 0.5)),
        ):
            with self.subTest(values=values):
                self.compare(values)

    def test_all_weight_scales_produce_their_normalized_decimal_distribution(self):
        for scale in (1e-300, 1.0, 1e300):
            with self.subTest(scale=scale):
                self.compare(((2.0, 0.25 * scale), (4.0, 0.75 * scale)))

    def test_constant_value_is_not_biased_by_accepted_probability_closure(self):
        self.compare(((1e100, 0.5), (1e100, 0.5 - ldexp(1.0, -42))))

    def test_percentile_uses_exact_mass_without_absolute_probability_floor(self):
        weight = nextafter(0.05, 0.0)
        values = ((0.0, weight), (1.0, 1.0 - weight))
        self.assertEqual(decimal_statistics(values)[3][0], 1.0)
        self.compare(values)
        self.compare(((0.0, 1.0), (1.0, 19.0)))

    def test_zero_weight_outliers_do_not_enter_distribution_support(self):
        self.compare(((-MAX, 0.0), (2.0, 0.5), (4.0, 0.5), (MAX, 0.0)))

    def test_permutations_do_not_change_exact_weighted_outputs(self):
        values = ((1e100, 0.5), (-1e100, 0.5), (3.0, TINY))
        # The mean is representable even though the tiny weighted term is
        # lost by ordinary summation after cancellation of the large terms.
        self.compare(values)
        self.compare(tuple(reversed(values)))

    def test_unbiased_sample_deviation_recovers_after_extreme_squaring(self):
        for values in ((1e-300, 2e-300), (MAX, nextafter(MAX, 0.0))):
            with localcontext() as context:
                context.prec = 2200
                rows = [Decimal.from_float(x) for x in values]
                mean = sum(rows) / len(rows)
                deviation = (sum((x - mean) ** 2 for x in rows) / (len(rows) - 1)).sqrt()
            try:
                actual = hierarchy._mean_and_deviation(values)
            except (DomainError, ArithmeticError) as error:
                self.fail(f"Representable unbiased deviation was lost: {error}")
            self.assertEqual(actual, (float(mean), float(deviation)))

    def test_relative_diagnostics_are_scale_invariant_and_do_not_overflow_difference(self):
        for x, y, expected in (
            (1e-300, 2e-300, 0.5),
            (1.0, 2.0, 0.5),
            (1e300, 2e300, 0.5),
            (MAX, -MAX, 2.0),
            (0.0, 0.0, 0.0),
        ):
            with self.subTest(x=x, y=y):
                self.assertEqual(hierarchy._relative_difference(x, y), expected)

    def test_empty_negative_and_zero_mass_distributions_fail_typed(self):
        for values in ((), ((1.0, 0.0),), ((1.0, -0.5), (2.0, 1.5))):
            with self.subTest(values=values), self.assertRaises(DomainError):
                hierarchy._weighted_statistics("invalid", values)

    def test_genuinely_unrepresentable_mandatory_mean_is_not_erased(self):
        with self.assertRaises(DomainError):
            hierarchy._weighted_statistics("unrepresentable", ((TINY, TINY), (0.0, 1.0)))

    def test_audit_impulse_and_qualified_scale_components_recover_before_export(self):
        draw = hierarchy.UncertaintyDraw(0, (("thrust_scale", 16.0), ("time_scale", 1.0 / 16)))
        self.assertEqual(hierarchy._impulses_for_draws(MAX, (draw,)), (MAX,))
        components = tuple(
            (f"thrust_scale::epistemic:{i}", value)
            for i, value in enumerate(
                (ldexp(1.0, 800), ldexp(1.0, 800), ldexp(1.0, -800), ldexp(1.0, -800))
            )
        )
        draw = hierarchy.UncertaintyDraw(0, components)
        self.assertEqual(hierarchy._scale_values(draw)["thrust_scale"], 1.0)
        self.assertEqual(hierarchy._impulses_for_draws(3.0, (draw,)), (3.0,))
        draw = hierarchy.UncertaintyDraw(0, (("thrust_scale", 16.0),))
        with self.assertRaises(DomainError):
            hierarchy._impulses_for_draws(MAX, (draw,))

    def test_square_root_midpoint_ties_use_destination_even_significand(self):
        first = 1.0
        second = nextafter(first, float("inf"))
        third = nextafter(second, float("inf"))
        cases = [
            ((Fraction(first) + Fraction(second)) / 2, first),
            ((Fraction(second) + Fraction(third)) / 2, third),
            (Fraction(3) * Fraction(TINY) / 2, 2 * TINY),
            (
                (Fraction(nextafter(sys.float_info.min, 0.0)) + Fraction(sys.float_info.min)) / 2,
                sys.float_info.min,
            ),
        ]
        for root, expected in cases:
            with self.subTest(root=root):
                actual = hierarchy._weighted_statistics("root-tie", ((-root, 0.5), (root, 0.5)))
                self.assertEqual(actual.standard_deviation, expected)
        for root in (Fraction(TINY) / 2, Fraction(2) ** 1024 - Fraction(2) ** 970):
            with self.subTest(root=root), self.assertRaises(DomainError):
                hierarchy._weighted_statistics("root-range", ((-root, 0.5), (root, 0.5)))


class PublicProbabilityAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scenarios, cls.definitions = public_case()
        cls.ensemble = build_hierarchical_propulsion_ensemble(
            cls.scenarios,
            cls.definitions,
            seed=71,
            replicate_tolerance=0.2,
        )

    def test_unrepresentable_joint_weight_is_a_typed_boundary_failure(self):
        scenarios, definitions = public_case(TINY)
        with self.assertRaises(DomainError) as caught:
            build_hierarchical_propulsion_ensemble(scenarios, definitions, seed=71)
        self.assertEqual(caught.exception.field, "Joint realization probability")

    def test_positive_representable_rare_event_is_retained_without_a_floor(self):
        scenarios, definitions = public_case(1e-20)
        result = build_hierarchical_propulsion_ensemble(scenarios, definitions, seed=71)
        rows = [x for x in result.realizations if x.scenario_id == "rare-off"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(x.joint_probability > 0 for x in rows))
        self.assertEqual(tuple(x.joint_probability for x in rows), (2.5e-21,) * 4)
        self.assertEqual(
            hierarchical_ensemble_from_dict(hierarchical_ensemble_to_dict(result)), result
        )

    def test_minimum_positive_single_draw_survives_when_metric_outputs_are_representable(self):
        # A larger mass/impulse difference keeps all mandatory variances
        # representable as well as the smallest positive joint probability.
        scenarios, definitions = public_case(TINY, samples=1, thrust=1000.0)
        result = build_hierarchical_propulsion_ensemble(scenarios, definitions, seed=71)
        rare = [row for row in result.realizations if row.scenario_id == "rare-off"]
        self.assertEqual(len(rare), 1)
        self.assertEqual(rare[0].joint_probability, TINY)
        self.assertEqual(result.zero_thrust_probability, TINY)
        self.assertEqual(result.joint_probability_sum, 1.0)
        self.assertEqual(result.joint_probability_closure_error, TINY)
        self.assertEqual(
            hierarchical_ensemble_from_dict(hierarchical_ensemble_to_dict(result)), result
        )

    def test_realization_record_cannot_claim_zero_for_positive_product(self):
        row = self.ensemble.realizations[0]
        with self.assertRaises(DomainError):
            replace(row, scenario_probability=1e-20, joint_probability=0.0)

    def test_total_variance_and_contributions_match_normalized_decimal_group_ledger(self):
        metric = "delivered_impulse_n_s"
        rows = [
            (x.artifact.delivered_total_impulse_n_s, x.joint_probability)
            for x in self.ensemble.realizations
        ]
        expected = decimal_statistics(rows)
        statistics = next(x for x in self.ensemble.statistics if x.metric == metric)
        self.assertEqual(statistics.expected_value, expected[0])
        with localcontext() as context:
            context.prec = 2200
            total = sum(Decimal.from_float(w) for x, w in rows)
            mean = sum(Decimal.from_float(x) * Decimal.from_float(w) for x, w in rows) / total
            variance = (
                sum(Decimal.from_float(w) * (Decimal.from_float(x) - mean) ** 2 for x, w in rows)
                / total
            )
        decomposition = next(x for x in self.ensemble.variance_decomposition if x.metric == metric)
        self.assertEqual(decomposition.total_variance, float(variance))
        displayed = [
            decomposition.total_variance,
            decomposition.within_scenario_variance,
            decomposition.between_scenario_variance,
        ]
        with localcontext() as context:
            context.prec = 2200
            residual = abs(
                Decimal.from_float(displayed[0])
                - Decimal.from_float(displayed[1])
                - Decimal.from_float(displayed[2])
            )
        self.assertEqual(decomposition.closure_error, float(residual))

    def test_statistics_policy_is_required_and_hash_bound(self):
        raw = hierarchical_ensemble_to_dict(self.ensemble)
        self.assertIn("statistics_numerics", raw)
        for value in (None, "unsupported"):
            changed = dict(raw)
            if value is None:
                del changed["statistics_numerics"]
            else:
                changed["statistics_numerics"] = value
            with self.assertRaises(InputError) as caught:
                hierarchical_ensemble_from_dict(changed)
            self.assertEqual(caught.exception.field, "statistics_numerics")

    def test_group_variance_matches_independent_conditional_decimal_sums(self):
        for scale in (1e-150, 1.0, 1e150):
            rows = []
            metric_values = {}
            for index, (value, scenario, probability) in enumerate(
                ((1.0, "a", 0.25), (1.5, "a", 0.25), (1.1, "b", 0.75), (1.9, "b", 0.75))
            ):
                artifact = replace(
                    self.ensemble.realizations[0].artifact, source_id=f"kernel:{index}"
                )
                row = replace(
                    self.ensemble.realizations[0],
                    scenario_id=scenario,
                    scenario_probability=probability,
                    conditional_probability=0.5,
                    joint_probability=probability * 0.5,
                    artifact=artifact,
                )
                rows.append(row)
                metric_values[artifact.source_id] = value * scale
            getter = lambda artifact: metric_values[artifact.source_id]
            with localcontext() as context:
                context.prec = 2200
                data = [
                    (
                        Decimal.from_float(getter(x.artifact)),
                        Decimal.from_float(x.joint_probability),
                        x.scenario_id,
                    )
                    for x in rows
                ]
                total = sum(w for x, w, g in data)
                mean = sum(x * w for x, w, g in data) / total
                variance = sum(w * (x - mean) ** 2 for x, w, g in data) / total
                within = Decimal(0)
                between = Decimal(0)
                for group in ("a", "b"):
                    subset = [(x, w) for x, w, g in data if g == group]
                    mass = sum(w for x, w in subset)
                    group_mean = sum(x * w for x, w in subset) / mass
                    within += sum(w * (x - group_mean) ** 2 for x, w in subset) / total
                    between += mass * (group_mean - mean) ** 2 / total
                expected = tuple(
                    float(x)
                    for x in (variance, within, between, within / variance, between / variance)
                )
            try:
                result = hierarchy._variance_decomposition(
                    "conditional-kernel", tuple(rows), getter
                )
            except (DomainError, ArithmeticError) as error:
                self.fail(f"Representable grouped variance must recover: {error}")
            self.assertEqual(
                (
                    result.total_variance,
                    result.within_scenario_variance,
                    result.between_scenario_variance,
                    result.within_fraction,
                    result.between_fraction,
                ),
                expected,
            )
            with localcontext() as context:
                context.prec = 2200
                residual = abs(
                    Decimal.from_float(result.total_variance)
                    - Decimal.from_float(result.within_scenario_variance)
                    - Decimal.from_float(result.between_scenario_variance)
                )
            self.assertEqual(result.closure_error, float(residual))

    def test_decomposition_does_not_erase_unrepresentable_mandatory_variance(self):
        rows = self.ensemble.realizations
        getter = lambda artifact: 1e-200 if artifact.delivered_total_impulse_n_s else 2e-200
        with self.assertRaises(DomainError):
            hierarchy._variance_decomposition("variance-range", rows, getter)

    def test_rehashed_derived_statistics_cannot_drift_by_one_output_ulp(self):
        from tests.test_ensemble_document_contract import rehash

        raw = hierarchical_ensemble_to_dict(self.ensemble)
        for family, field in (
            ("statistics", "expected_value"),
            ("variance_decomposition", "total_variance"),
            ("contributions", "joint_expected_value_contribution"),
        ):
            import copy

            changed = copy.deepcopy(raw)
            changed[family][0][field] = nextafter(changed[family][0][field], float("inf"))
            with self.subTest(family=family), self.assertRaises(InputError):
                hierarchical_ensemble_from_dict(rehash(changed))
