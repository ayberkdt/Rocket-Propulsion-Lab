"""Attack affine decision ties and integrity of hypothetical mission scaling.

References: NASA-SE-6 decision/input traceability; NASA-7009B numerical
verification only. Preference arithmetic is not hardware qualification.
"""

import copy
import math
import unittest
from decimal import Context, Decimal, localcontext

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import (
    Criterion, DesignCase, design_sensitivity, project_from_dict, weighted_trade,
)
from tests.test_sass_studies import example_request, q


class DecisionBoundaryAttacks(unittest.TestCase):
    def test_equal_affine_scores_cannot_gain_a_false_order_from_term_rounding(self):
        criteria = tuple(Criterion(str(k), q(1), q(0), q(1 / 3)) for k in range(3))
        utilities = {"a": (0., .5, 1.), "b": (0., .625, .875)}
        rows = weighted_trade({name: {str(k): q(value) for k, value in enumerate(values)}
                               for name, values in utilities.items()}, criteria, ())
        with localcontext(Context(prec=115)):
            expected = float(Decimal(1 / 3) * sum(map(Decimal, utilities["a"])))
        self.assertEqual(rows[0].score.value, expected)
        self.assertEqual(rows[1].score.value, expected)
        self.assertEqual(rows[0].score.value, rows[1].score.value)

    def test_small_weighted_terms_must_not_be_erased_before_the_score_sum(self):
        smallest = math.nextafter(0., 1.)
        criteria = (Criterion("a", q(1), q(0), q(.5)), Criterion("b", q(1), q(0), q(.5)))
        row = weighted_trade({"x": {"a": q(smallest), "b": q(smallest)}}, criteria, ())[0]
        self.assertEqual(row.score.value, smallest)

    def test_positive_scaled_delta_v_cannot_silently_become_zero(self):
        request = example_request()
        for item in request["items"]:
            item["nominal_delta_v"]["value"] = 1e-300
            item["margin_fraction"]["value"] = 0
            item["count"] = 1
        self.assertGreater(project_from_dict(request).budget.get("allocated_delta_v").value, 0)
        with self.assertRaises(DomainError):
            design_sensitivity(request, (DesignCase("unrepresentable positive change", q(1e-30), q(1)),))

    def test_explicit_zero_scale_is_distinct_from_lost_positive_requirement(self):
        request = example_request()
        original = copy.deepcopy(request)
        result = design_sensitivity(request, (DesignCase("zero mission allocation", q(0), q(1)),))[0]
        self.assertEqual(result.budget.get("allocated_delta_v").value, 0)
        self.assertEqual(request, original)


if __name__ == "__main__":
    unittest.main()
