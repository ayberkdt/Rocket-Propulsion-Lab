"""Hostile numerical and evidence checks; see docs/PROPULSION_ASSURANCE.md.

Decimal oracles use 75 digits and binary inputs converted exactly, independently
of the production implementation. Synthetic cases are code verification only.
"""

from __future__ import annotations

import json
import math
import random
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2 as G0
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass import (
    Candidate,
    CitedScalar,
    assumption,
    bipropellant_split,
    circular_spiral_delta_v,
    hohmann_transfer,
    ideal_cold_gas_storage,
    list_thrusters,
    nominal_operating_point,
    perigee_lowering_delta_v,
    project_from_json,
    size_mission,
)
from rocket_propulsion.studies.sass.evidence import AnalysisResult
from tests.test_sass_studies import example_request, maneuver


def q(value: float, unit: str = "1") -> CitedScalar:
    """Synthetic input only; NASA-SE-6 assumption recording."""
    return assumption(value, unit, "Synthetic adversarial verification input.")


class HostileRegressionTests(unittest.TestCase):
    def test_neighbouring_orbits_preserve_both_small_signed_impulses(self):
        for r0 in (7e6, 4.2e7):
            for r1 in (math.nextafter(r0, math.inf), math.nextafter(r0, 0)):
                with localcontext() as ctx:
                    ctx.prec = 75
                    a, b, mu = map(Decimal, (r0, r1, 3.986004418e14))
                    v0, v1 = (mu / a).sqrt(), (mu / b).sqrt()
                    expected = ((2 * mu * b / (a * (a + b))).sqrt() - v0,
                                v1 - (2 * mu * a / (b * (a + b))).sqrt())
                    spiral = abs(v0 - v1)
                result = hohmann_transfer(q(r0, "m"), q(r1, "m"), q(float(mu), "m^3/s^2"))
                for channel, target in zip(("departure_delta_v", "arrival_delta_v"), expected):
                    with self.subTest(r0=r0, r1=r1, channel=channel):
                        self.assertTrue(math.isclose(result.get(channel).value, float(target),
                                                     rel_tol=3e-14, abs_tol=0))
                actual = circular_spiral_delta_v(q(r0, "m"), q(r1, "m"),
                                                 q(float(mu), "m^3/s^2"))
                self.assertTrue(math.isclose(actual.get("total_delta_v").value, float(spiral),
                                             rel_tol=3e-14, abs_tol=0))

    def test_large_mass_ratio_keeps_positive_dry_mass(self):
        result = size_mission(q(10, "kg"), (maneuver("severe", 100 * G0 * 200),),
                              residual_fraction=q(0), fixed_reserve=q(0, "kg"))
        self.assertAlmostEqual(result.summary.get("final_mass").value, 10, delta=2e-12)
        self.assertGreater(result.burns[0][1].get("final_mass").value, 0)

    def test_small_oxidizer_inventory_does_not_disappear(self):
        result = bipropellant_split(q(1, "kg"), q(1e-20))
        self.assertTrue(math.isclose(result.get("oxidizer_mass").value, 1e-20,
                                     rel_tol=1e-15, abs_tol=0))

    def test_small_gas_residual_does_not_disappear(self):
        result = ideal_cold_gas_storage(q(1, "kg"), q(1e20, "Pa"), q(1, "Pa"),
                                       q(300, "K"), q(300, "J/(kg*K)"))
        self.assertTrue(math.isclose(result.get("remaining_gas").value, 1e-20,
                                     rel_tol=1e-15, abs_tol=0))

    def test_result_does_not_retain_mutable_caller_containers(self):
        values, limits = [["mass", q(1, "kg")]], ["Synthetic model."]
        result = AnalysisResult("synthetic", values, limits)
        values[0][1] = q(2, "kg")
        limits[0] = "Changed."
        self.assertEqual(result.get("mass").value, 1)
        self.assertEqual(result.limitations, ("Synthetic model.",))

    def test_assumption_cannot_be_relabelled_as_published(self):
        with self.assertRaises(DomainError):
            replace(q(99, "N"), status="published", note="")

    def test_duplicate_json_fields_cannot_override_evidence(self):
        payload = json.dumps(example_request())
        payload = payload.replace('"schema":', '"schema": "tampered", "schema":', 1)
        with self.assertRaises((DomainError, InputError)):
            project_from_json(payload)

    def test_unrepresentable_counts_fail_with_domain_error(self):
        with self.assertRaises(DomainError):
            Candidate("bad", "moog-058e151", q(1, "kg"), 10**400, 1, q(1), q(1))

    def test_no_catalog_entry_can_claim_validated_joint_performance(self):
        for item in list_thrusters():
            with self.subTest(item=item.key), self.assertRaises(DomainError):
                nominal_operating_point(item.key, require_validated=True)

    def test_single_impulse_does_not_require_a_representable_orbit_period(self):
        with localcontext() as ctx:
            ctx.prec = 75
            ra, rp = Decimal.from_float(2e210), Decimal.from_float(1e210)
            expected = (1 / ra).sqrt() - (2 * rp / (ra * (ra + rp))).sqrt()
        result = perigee_lowering_delta_v(q(float(ra), "m"), q(float(rp), "m"),
                                          q(1, "m^3/s^2"))
        self.assertTrue(math.isclose(result.get("total_delta_v").value, float(expected),
                                     rel_tol=3e-14, abs_tol=0))


class IndependentMassOracleTests(unittest.TestCase):
    def test_seeded_ordered_burns_against_backward_decimal_propagation(self):
        rng = random.Random(792681)
        for case in range(250):
            dry, reserve, fraction = 10**rng.uniform(-2, 4), rng.random(), rng.random() * .08
            burns = tuple(maneuver(str(i), rng.uniform(0, 400), rng.uniform(60, 400),
                                   10**rng.uniform(-3, 2), rng.uniform(.3, 1),
                                   rng.uniform(.05, 1)) for i in range(rng.randrange(1, 7)))
            with localcontext() as ctx:
                ctx.prec = 75
                # Solve for loaded propellant by bisection, propagating masses
                # backward from dry + retained, not the production closure formula.
                lo, hi = Decimal(0), Decimal(dry + reserve + 1)
                def needed(loaded, *, dry=dry, reserve=reserve, fraction=fraction, burns=burns):
                    mass = Decimal(dry) + Decimal(reserve) + Decimal(fraction) * loaded
                    for burn in reversed(burns):
                        exponent = Decimal(burn.delta_v.value) / (
                            Decimal(G0) * Decimal(burn.system_isp.value)
                            * Decimal(burn.axial_efficiency.value))
                        mass *= exponent.exp()
                    return mass - Decimal(dry) - loaded
                for _ in range(100):
                    if needed(hi) <= 0:
                        break
                    hi *= 2
                else:
                    with self.assertRaises(DomainError):
                        size_mission(q(dry, "kg"), burns, residual_fraction=q(fraction),
                                     fixed_reserve=q(reserve, "kg"))
                    continue
                for _ in range(260):
                    mid = (lo + hi) / 2
                    if needed(mid) > 0:
                        lo = mid
                    else:
                        hi = mid
                expected = float((lo + hi) / 2)
            result = size_mission(q(dry, "kg"), burns, residual_fraction=q(fraction),
                                  fixed_reserve=q(reserve, "kg"))
            with self.subTest(case=case):
                self.assertTrue(math.isclose(result.summary.get("loaded_propellant").value,
                                             expected, rel_tol=2e-13, abs_tol=1e-13))
                self.assertTrue(math.isclose(result.summary.get("final_mass").value,
                                             dry + fraction * expected + reserve,
                                             rel_tol=2e-13, abs_tol=1e-13))


if __name__ == "__main__":
    unittest.main()
