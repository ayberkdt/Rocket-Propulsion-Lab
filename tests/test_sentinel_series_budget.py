"""Independent numerical and bookkeeping controls for the assignment companion."""
import copy
import importlib.util
import json
import math
import subprocess
import sys
import unittest
from decimal import Decimal, localcontext
from itertools import pairwise
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/sentinel_series_budget.py"
DATA = ROOT / "docs/missions/sentinel/SERI_VERILERI.json"
spec = importlib.util.spec_from_file_location("sentinel_budget", SCRIPT)
assert spec is not None and spec.loader is not None
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


class SentinelSeriesBudgetTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(DATA.read_text(encoding="utf-8"))

    def test_rocket_equation_against_high_precision_decimal(self):
        with localcontext() as ctx:
            ctx.prec = 60
            for mass, fuel, isp in [(2186, 154, 205), (1143, 133, 220), (1000, 1e-10, 200)]:
                m, p, i = (Decimal(str(v)) for v in (mass, fuel, isp))
                expected = Decimal("9.80665") * i * (m / (m - p)).ln()
                self.assertTrue(math.isclose(model.capacity(mass, fuel, isp), float(expected),
                                             rel_tol=3e-15))

    def test_tiny_burn_does_not_cancel_to_zero(self):
        with localcontext() as ctx:
            ctx.prec = 60
            expected = Decimal(1000) * (1 - (-Decimal("1e-12")
                                            / (Decimal("9.80665") * 220)).exp())
        self.assertTrue(math.isclose(model.propellant(1000, 1e-12, 220), float(expected),
                                     rel_tol=3e-15))

    def test_fuel_mass_inverse_and_conservation(self):
        for mass, dv, isp in [(2186, 100, 205), (1143, 7.97, 220), (900, 73, 200)]:
            used = model.propellant(mass, dv, isp)
            self.assertAlmostEqual(model.capacity(mass, used, isp), dv, places=12)

    def test_phases_agree_with_one_total_exponential(self):
        rows = model.ledger(1143, 133, self.data["s2c_events"], 220)
        dv = sum(math.hypot(*p["dv_components_m_s"]) for p in self.data["s2c_events"])
        self.assertAlmostEqual(sum(r["propellant_ideal_kg"] for r in rows),
                               model.propellant(1143, dv, 220), places=12)
        for previous, following in pairwise(rows):
            self.assertEqual(previous["mass_after_kg"], following["mass_before_kg"])
        self.assertAlmostEqual(rows[-1]["mass_after_kg"]
                               - rows[-1]["fuel_remaining_ideal_kg"], 1010, places=10)

    def test_opposite_cam_burns_do_not_cancel_fuel_cost(self):
        events = [{"id": "out", "dv_components_m_s": [-0.02]},
                  {"id": "back", "dv_components_m_s": [0.02]}]
        rows = model.ledger(1000, 10, events, 220)
        self.assertAlmostEqual(sum(r["dv_m_s"] for r in rows), 0.04)
        self.assertGreater(sum(r["propellant_ideal_kg"] for r in rows), 0.018)

    def test_combined_oop_ip_is_vector_norm(self):
        rows = model.ledger(1143, 133, self.data["s2c_events"], 220)
        self.assertAlmostEqual(rows[4]["dv_m_s"], math.sqrt(.9344**2 + .0537**2), places=14)
        self.assertNotAlmostEqual(rows[4]["dv_m_s"], .9344 + .0537, places=3)

    def test_october31_not_double_counted(self):
        ids = [p["id"] for p in self.data["s2c_events"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(sum("2024-10-31" in key for key in ids), 1)
        answer = model.calculate(self.data)
        # Independent scalar total plus the single orthogonal two-axis event.
        expected = 7.0325 + math.sqrt(.9344**2 + .0537**2)
        self.assertAlmostEqual(answer["s2c_total_dv_m_s"], expected, places=12)

    def test_commanded_not_substituted_for_calibrated(self):
        record = next(p for p in self.data["s2c_events"] if "2024-11-06" in p["id"])
        self.assertEqual(record["dv_components_m_s"], [.9344, .0537])
        self.assertEqual(record["duration_s"], 314.1)

    def test_duplicate_phase_is_rejected(self):
        events = self.data["s2c_events"]
        with self.assertRaisesRegex(ValueError, "double-counted"):
            model.ledger(1143, 133, events + events[:1], 220)

    def test_redundant_vehicle_fuel_mixing_is_rejected(self):
        self.data["capacity_cases"]["S1C"]["fuel_vehicle"] = "S1A"
        with self.assertRaisesRegex(ValueError, "different vehicles"):
            model.calculate(self.data)

    def test_schema_prevents_silent_unit_convention_change(self):
        self.data["schema"] = "km-and-grams"
        with self.assertRaisesRegex(ValueError, "schema"):
            model.calculate(self.data)

    def test_cant_penalty_increases_scalar_propellant(self):
        aligned = model.propellant(1143, 8, 220, 1)
        canted = model.propellant(1143, 8, 220, .8)
        self.assertGreater(canted, aligned)
        self.assertAlmostEqual(model.capacity(1143, canted, 220) * .8, 8)

    def test_finite_burn_inferred_thrust_uses_post_acquisition_mass(self):
        rows = model.ledger(1143, 133, self.data["s2c_events"], 220)
        manoeuvre = rows[4]
        self.assertLess(manoeuvre["mass_before_kg"], 1143)
        self.assertTrue(3.38 < manoeuvre["inferred_average_axial_thrust_N"] < 3.41)

    def test_hohmann_against_energy_and_angular_momentum(self):
        # Independent transfer energy, then h = sqrt(mu * semilatus rectum).
        for h1, h2 in [(786000, 653000), (693000, 530000)]:
            r1, r2 = model.EARTH_RADIUS + h1, model.EARTH_RADIUS + h2
            momentum = math.sqrt(2 * model.MU * r1 * r2 / (r1 + r2))
            expected = abs(momentum / r1 - math.sqrt(model.MU / r1)) + abs(
                momentum / r2 - math.sqrt(model.MU / r2))
            self.assertAlmostEqual(model.lowering_dv(h1, h2, h2), expected, places=10)

    def test_final_ellipse_not_old_circular_disposal_goal(self):
        final = model.lowering_dv(693000, 580000, 610000)
        old = model.lowering_dv(693000, 530000, 530000)
        self.assertTrue(50 < final < 55)
        self.assertTrue(88 < old < 89)
        self.assertEqual(model.lowering_dv(693000, 693000, 693000), 0)

    def test_orbit_control_linearisation_against_published_rounded_examples(self):
        # #617 p4: 1 cm/s -> about 160 m and .342 s/day. Model differs <2%.
        result = model.drift(.01, 86400, 786000)
        self.assertLess(abs(result["ground_track_cross_track_m"] / 160 - 1), .02)
        self.assertLess(abs(result["along_track_time_s"] / .342 - 1), .02)
        result = model.drift(.0385, 86400, 786000)
        self.assertLess(abs(result["ground_track_cross_track_m"] / 614 - 1), .02)

    def test_decimal_gigabytes_and_per_carrier_rates(self):
        result = model.calculate(self.data)
        self.assertEqual(result["downlink_600s_520Mbps_GB"], 39)
        self.assertEqual(result["downlink_600s_560Mbps_GB"], 42)
        self.assertAlmostEqual(result["s2_300GB_emptying_at_560Mbps_s"], 4285.714285714286)

    def test_catalogue_nominals_are_not_one_measured_operating_point(self):
        result = model.calculate(self.data)
        self.assertAlmostEqual(result["catalogue_0_44g_s_220s_implied_N"], .94928372)
        self.assertNotAlmostEqual(result["one_N_220s_mass_flow_kg_s"], .00044, places=6)

    def test_period_eclipse_and_heater_dimensions(self):
        result = model.calculate(self.data)
        self.assertAlmostEqual(result["orbit_s2_786km"]["period_s"] / 60, 100.575, delta=.01)
        self.assertAlmostEqual(result["orbit_s1_693km"]["beta_zero_eclipse_s"] / 60,
                               35.30, delta=.01)
        self.assertAlmostEqual(result["rafael_heater_28V_257ohm_W_per_element"], 3.05058,
                               delta=.00001)

    def test_radiator_scaling(self):
        a = model.thermal_radiator(300, .8, 300)
        self.assertAlmostEqual(model.thermal_radiator(600, .8, 300), 2 * a)
        self.assertAlmostEqual(model.thermal_radiator(300, .8, 600), a / 16)

    def test_invalid_domains_do_not_silently_produce_budgets(self):
        calls = [lambda v=v: model.propellant(1000, v, 220)
                 for v in [True, "1", float("nan"), float("inf"), -1]]
        calls += [lambda: model.propellant(1000, 1, 0),
                  lambda: model.propellant(1000, 1, 220, 1.1),
                  lambda: model.capacity(1000, 1000, 220),
                  lambda: model.capacity(1000, -1, 220),
                  lambda: model.lowering_dv(500000, 600000, 700000),
                  lambda: model.thermal_radiator(300, 1.1, 300),
                  lambda: model.ledger(1000, .0001, self.data["s2c_events"], 220)]
        for call in calls:
            with self.subTest(call=call), self.assertRaises((TypeError, ValueError)):
                call()

    def test_empty_and_nonfinite_vectors_rejected(self):
        for vector in [[], [float("nan")], [True], ["0.02"]]:
            with self.subTest(vector=vector), self.assertRaises(ValueError):
                model.ledger(1143, 133, [{"id": "bad", "dv_components_m_s": vector}], 220)

    def test_insufficient_fuel_rejected_even_when_net_dv_zero(self):
        events = [{"id": "one", "dv_components_m_s": [1]},
                  {"id": "two", "dv_components_m_s": [-1]}]
        with self.assertRaisesRegex(ValueError, "insufficient"):
            model.ledger(1000, .5, events, 220)

    def test_published_data_are_not_mutated(self):
        original = copy.deepcopy(self.data)
        model.calculate(self.data)
        self.assertEqual(self.data, original)

    def test_cli_runs_without_site_packages_and_reproduces_saved_results(self):
        completed = subprocess.run([sys.executable, "-I", "-S", str(SCRIPT), "--data", str(DATA)],
                                   check=True, capture_output=True, text=True, encoding="utf-8")
        answer = json.loads(completed.stdout)
        expected = json.loads((DATA.parent / "SERI_SONUCLARI.json").read_text(encoding="utf-8"))
        self.assertEqual(answer, expected)


if __name__ == "__main__":
    unittest.main()
