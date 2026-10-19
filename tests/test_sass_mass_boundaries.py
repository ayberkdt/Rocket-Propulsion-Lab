"""Mission closure attacks at numerical feasibility boundaries.

Independent oracle: forward rocket-equation propagation plus a bracketed root
search, not the production closed-form propellant expression. These synthetic
cases establish numerical behavior only, not a feasible hardware design.
References: NASA-ROCKET; NASA-7009B numerical verification and domain records.
"""

import math
import unittest
from decimal import Context, Decimal, ROUND_DOWN, localcontext

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2 as G0
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import size_mission
from tests.test_sass_studies import maneuver, q


def forward_root(dry, reserve, fraction, burns):
    """NASA-ROCKET: solve final mass = dry + residual by forward propagation."""
    with localcontext(Context(prec=115)):
        survival = tuple((-Decimal(b.delta_v.value) / (Decimal(G0)
                          * Decimal(b.system_isp.value) * Decimal(b.axial_efficiency.value))).exp()
                         for b in burns)
        total_survival = math.prod(survival)
        if total_survival <= Decimal(fraction):
            return None
        def residual(loaded):
            mass = Decimal(dry) + loaded
            for ratio in survival:
                mass *= ratio
            return mass - Decimal(dry) - Decimal(reserve) - Decimal(fraction) * loaded
        lower, upper = Decimal(0), max(Decimal(1), Decimal(dry) + Decimal(reserve))
        for _ in range(2200):
            if residual(upper) >= 0:
                break
            upper *= 2
        else:
            raise AssertionError("Independent oracle failed to bracket the mission.")
        for _ in range(360):
            middle = (lower + upper) / 2
            if residual(middle) < 0:
                lower = middle
            else:
                upper = middle
        return float((lower + upper) / 2)


class MissionBoundaryAttacks(unittest.TestCase):
    def test_near_unit_residual_cannot_hide_fifty_percent_propellant_error(self):
        fraction = math.nextafter(1., 0.)
        burns = (maneuver("boundary", (1 - fraction) * .5 * G0 * 200),)
        expected = forward_root(10., 0., fraction, burns)
        result = size_mission(q(10, "kg"), burns, residual_fraction=q(fraction),
                              fixed_reserve=q(0, "kg"))
        self.assertTrue(math.isclose(result.summary.get("loaded_propellant").value,
                                     expected, rel_tol=3e-14, abs_tol=0))

    def test_adjacent_feasibility_sides_agree_with_independent_forward_root(self):
        for fraction in (.01, .2, .8):
            with localcontext(Context(prec=115)):
                boundary = float(-Decimal(fraction).ln() * Decimal(G0) * Decimal(200))
            for direction in (0., math.inf):
                dv = boundary
                for _ in range(4):
                    dv = math.nextafter(dv, direction)
                burns = (maneuver("boundary", dv),)
                expected = forward_root(10., .3, fraction, burns)
                with self.subTest(fraction=fraction, direction=direction):
                    if expected is None:
                        with self.assertRaises(DomainError):
                            size_mission(q(10, "kg"), burns, residual_fraction=q(fraction),
                                         fixed_reserve=q(.3, "kg"))
                    else:
                        result = size_mission(q(10, "kg"), burns, residual_fraction=q(fraction),
                                              fixed_reserve=q(.3, "kg"))
                        self.assertTrue(math.isclose(result.summary.get("loaded_propellant").value,
                                                     expected, rel_tol=3e-14, abs_tol=0))

    def test_large_mass_ratio_is_allowed_when_reported_quantities_are_representable(self):
        burns = (maneuver("large ratio", 1000 * G0 * 200),)
        expected = forward_root(1e-300, 0., 0., burns)
        result = size_mission(q(1e-300, "kg"), burns, residual_fraction=q(0),
                              fixed_reserve=q(0, "kg"))
        self.assertTrue(math.isclose(result.summary.get("loaded_propellant").value,
                                     expected, rel_tol=3e-14, abs_tol=0))
        self.assertTrue(math.isclose(result.summary.get("final_mass").value,
                                     1e-300, rel_tol=3e-14, abs_tol=0))

    def test_nonzero_burn_cannot_be_reported_as_zero_consumption_from_underflow(self):
        with self.assertRaises(DomainError):
            size_mission(q(1, "kg"), (maneuver("underflow", math.nextafter(0., 1.)),),
                         residual_fraction=q(0), fixed_reserve=q(0, "kg"))

    def test_tiny_exponent_with_large_mass_preserves_a_representable_propellant_load(self):
        result = size_mission(q(1e300, "kg"), (maneuver("small exponent", 1e-200),),
                              residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        # Taylor bound: expm1(s)/s differs from one by <2*s here, far below
        # binary rounding. This oracle does not subtract Decimal exp(s)-1.
        expected = 1e300 * (1e-200 / (G0 * 200))
        self.assertTrue(math.isclose(result.summary.get("loaded_propellant").value,
                                     expected, rel_tol=3e-14, abs_tol=0))

    def test_caller_decimal_context_cannot_change_the_mission(self):
        burns = (maneuver("a", 100),)
        expected = size_mission(q(100, "kg"), burns, residual_fraction=q(.02),
                                fixed_reserve=q(1, "kg"))
        with localcontext(Context(prec=6, rounding=ROUND_DOWN, Emax=9, Emin=-9)):
            actual = size_mission(q(100, "kg"), burns, residual_fraction=q(.02),
                                  fixed_reserve=q(1, "kg"))
        self.assertEqual(actual.to_dict(), expected.to_dict())

    def test_closure_evidence_distinguishes_reported_rounding_from_working_residuals(self):
        result = size_mission(q(10, "kg"), (maneuver("rounding", 100),),
                              residual_fraction=q(.02), fixed_reserve=q(.5, "kg"))
        values = {name: value.value for name, value in result.summary.values}
        self.assertEqual(values["mass_closure_residual"], math.fsum((
            values["final_mass"], -values["dry_mass"], -values["retained_propellant"])))
        self.assertEqual(values["propellant_closure_residual"], math.fsum((
            values["loaded_propellant"], -values["consumed_propellant"],
            -values["retained_propellant"])))
        evidence = result.to_dict()["numerics"]
        self.assertGreaterEqual(len(evidence["checked_precision_digits"].split(",")), 2)
        self.assertLess(abs(Decimal(evidence["mass_closure_residual_decimal_kg"])), Decimal("1e-70"))
        self.assertGreater(Decimal(evidence["denominator_subtraction_condition_decimal"]), 1)
        self.assertIn("physical uncertainty is not quantified", evidence["interpretation"])


if __name__ == "__main__":
    unittest.main()
