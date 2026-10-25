"""HTTP routes for heat transfer, regenerative cooling and combustion products/losses."""

from __future__ import annotations

import json
import unittest

from rocket_propulsion.api.routes import ROUTE_METADATA, dispatch_calculation
from rocket_propulsion.core.errors import DomainError, InputError

ROUTES = ("heat-conduction", "heat-convection", "heat-radiation", "regen-chamber",
          "combustion-products", "combustion-losses")


class ThermalRouteDefaultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.responses = {name: dispatch_calculation(f"/api/v1/{name}", {}) for name in ROUTES}

    def test_every_route_answers_with_strict_json_and_metadata(self) -> None:
        for name, response in self.responses.items():
            with self.subTest(route=name):
                json.dumps(response, allow_nan=False)
                self.assertIn(f"/api/v1/{name}", ROUTE_METADATA)

    def test_conduction_energy_audit_closes(self) -> None:
        response = self.responses["heat-conduction"]
        self.assertLess(abs(response["energy_residual_j_m2"]),
                        1.0e-9 * abs(response["stored_energy_change_j_m2"]))
        self.assertGreater(response["final_hot_surface_k"], response["final_cold_surface_k"])

    def test_radiation_wall_cools_as_emissivity_rises(self) -> None:
        temperatures = [row["wall_temperature_k"] for row in self.responses["heat-radiation"]["sweep"]]
        self.assertEqual(temperatures, sorted(temperatures, reverse=True))

    def test_regen_route_conserves_energy(self) -> None:
        response = self.responses["regen-chamber"]
        self.assertLess(abs(response["energy_balance_residual_w"]), 1.0e-9 * response["total_heat_load_w"])

    def test_products_show_recombination_and_losses_chain_multiplies(self) -> None:
        radicals = [row["radical_mole_fraction"] for row in self.responses["combustion-products"]["stations"]]
        self.assertEqual(radicals, sorted(radicals, reverse=True))
        chain = self.responses["combustion-losses"]["chain"]
        self.assertAlmostEqual(chain["overall_efficiency"],
                               chain["delivered_specific_impulse_s"] / chain["ode_specific_impulse_s"], delta=1.0e-12)
        self.assertLess(chain["kinetic_efficiency"], 1.0)


class ThermalRouteInputTests(unittest.TestCase):
    def test_invalid_inputs_are_refused(self) -> None:
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/heat-conduction", {"duration_s": 0})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/combustion-losses", {"propellant_pair": "LOX_CH4", "compute_kinetics": True})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/combustion-products", {"propellant_pair": "UNKNOWN"})


if __name__ == "__main__":
    unittest.main()
