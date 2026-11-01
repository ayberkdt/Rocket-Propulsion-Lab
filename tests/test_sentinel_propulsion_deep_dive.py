"""Independent physical balances and failure cases for the academic V5 model."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from sentinel_propulsion_deep_dive import burn, cutoff_residual, finite_burn


class DeepDiveTests(unittest.TestCase):
    def test_feed_cutoff_closes_pv_and_mass(self):
        r = cutoff_residual(.204, 154, 1004, 22, 5.5, .5, 1.2)
        self.assertAlmostEqual(22*(r["initial_gas_L"])**1.2,
                               6*r["cutoff_gas_L"]**1.2, places=9)
        self.assertAlmostEqual(r["accessible_kg"]+r["residual_kg"],154)
        self.assertAlmostEqual(r["residual_kg"],1004*(.204-r["cutoff_gas_L"]/1000))

    def test_zero_residual_at_reachable_empty_tank(self):
        self.assertEqual(cutoff_residual(.2, 100, 1000, 22, 5.5, 0, 1)["residual_kg"], 0)

    def test_adiabatic_and_loss_reduce_accessible_mass(self):
        baseline=cutoff_residual(.204,154,1004,22,5.5,0,1)["accessible_kg"]
        cooled=cutoff_residual(.204,154,1004,22,5.5,.5,5/3)["accessible_kg"]
        self.assertLess(cooled,baseline)

    def test_cant_changes_propellant_but_duty_changes_time(self):
        a=burn(1143,1,220,4,1,30)
        b=burn(1143,1,220,4,1,30,.5)
        c=burn(1143,1,220,4,1,0)
        self.assertEqual(a["propellant_kg"],b["propellant_kg"])
        self.assertAlmostEqual(b["duration_s"],2*a["duration_s"])
        self.assertGreater(a["propellant_kg"],c["propellant_kg"])
        self.assertAlmostEqual(a["mass_flow_kg_s"]*a["duration_s"],a["propellant_kg"])
        self.assertAlmostEqual(a["axial_thrust_N"]*a["duration_s"],
                               a["propellant_kg"]*220*9.80665*math.cos(math.pi/6))

    def test_force_changes_time_not_scalar_rocket_budget(self):
        a=burn(2186,1,205,1,1,0); b=burn(2186,1,205,1,.5,0)
        self.assertEqual(a["propellant_kg"],b["propellant_kg"])
        self.assertEqual(b["duration_s"],2*a["duration_s"])

    def test_finite_burn_conserves_work_momentum_and_mass(self):
        a=finite_burn(693000,2100,.7,2.18,205)
        self.assertLess(abs(a["energy_work_residual_J_kg"]),1e-5)
        self.assertLess(abs(a["angular_momentum_residual_m2_s"]),.001)
        self.assertLess(abs(a["mass_residual_kg"]),1e-7)
        self.assertLess(a["perigee_km"],a["apogee_km"])
        self.assertLess(a["apogee_km"],693)

    def test_step_halving_preserves_osculating_elements(self):
        a=finite_burn(693000,2100,.7,2.18,205)
        b=finite_burn(693000,2100,.7,2.18,205,.5)
        for k in ("perigee_km","apogee_km"):
            self.assertAlmostEqual(a[k],b[k],places=7)

    def test_shorter_burn_approaches_impulsive_eccentricity(self):
        a=finite_burn(693000,2100,.7,1,205)
        b=finite_burn(693000,2100,.7,20,205,.1)
        self.assertLess(abs(b["perigee_km"]-b["impulsive_perigee_km"]),
                        abs(a["perigee_km"]-a["impulsive_perigee_km"]))

    def test_invalid_domains_refused(self):
        for p in (float("nan"),True,-1):
            with self.subTest(p=p),self.assertRaises((TypeError,ValueError)):
                cutoff_residual(.204,154,1004,p,5.5,0,1)
        with self.assertRaises(ValueError): cutoff_residual(.1,154,1004,22,5.5,0,1)
        with self.assertRaises(ValueError): cutoff_residual(.204,154,1004,5,5.5,0,1)
        with self.assertRaises(ValueError): burn(1143,1,220,4,1,90)
        with self.assertRaises(ValueError): burn(1143,1,220,True,1,0)
        with self.assertRaises(ValueError): finite_burn(693000,2100,1000,1,205)


if __name__ == "__main__": unittest.main()
