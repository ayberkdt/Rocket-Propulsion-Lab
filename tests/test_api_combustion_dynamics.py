"""HTTP routes for transpiration cooling, combustion stability and propellant vaporization."""

from __future__ import annotations

import json
import unittest
from math import pi

from rocket_propulsion.api.routes import ROUTE_METADATA, dispatch_calculation
from rocket_propulsion.core.errors import DomainError, InputError

ROUTES = ("transpiration-cooling", "combustion-stability", "propellant-vaporization")


class CombustionDynamicsRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.responses = {name: dispatch_calculation(f"/api/v1/{name}", {}) for name in ROUTES}

    def test_strict_json_and_metadata(self) -> None:
        for name, response in self.responses.items():
            with self.subTest(route=name):
                json.dumps(response, allow_nan=False)
                self.assertIn(f"/api/v1/{name}", ROUTE_METADATA)

    def test_transpiration_balance_and_monotone_sweep(self) -> None:
        response = self.responses["transpiration-cooling"]
        required = response["required"]
        self.assertLess(abs(required["absorbed_heat_flux_w_m2"] / required["convective_heat_flux_w_m2"] - 1), 1e-10)
        self.assertAlmostEqual(required["wall_temperature_k"], 900.0, places=9)
        walls = [row["wall_temperature_k"] for row in response["sweep"]]
        self.assertEqual(walls, sorted(walls, reverse=True))
        self.assertAlmostEqual(walls[0], 3300.0, places=6)  # no coolant: recovery temperature

    def test_stability_route_is_self_consistent(self) -> None:
        response = self.responses["combustion-stability"]
        chamber, chug, nyquist = response["chamber"], response["chugging"], response["nyquist"]
        # Equal injector drops make K = p_c / (2 Dp) exactly.
        self.assertAlmostEqual(chug["gain"], 1.0 / (2 * 0.2), places=12)
        self.assertEqual(chug["stable"], nyquist["stable"])
        self.assertLess(abs(nyquist["gain_crossover_hz"][0] / chug["neutral_frequency_hz"] - 1), 2e-3)
        names = [m["name"] for m in response["modes"]]
        self.assertEqual(names[0], "1L")
        self.assertIn("1T", names)
        one_t = next(m for m in response["modes"] if m["name"] == "1T")
        expected = one_t["eigenvalue"] * chamber["sound_speed_m_s"] / (2 * pi * 0.10)
        self.assertLess(abs(one_t["frequency_hz"] / expected - 1), 1e-12)

    def test_longer_lag_than_critical_is_unstable(self) -> None:
        critical = self.responses["combustion-stability"]["chugging"]["critical_time_lag_s"]
        response = dispatch_calculation("/api/v1/combustion-stability",
                                        {"combustion_time_lag_s": 1.05 * critical})
        self.assertFalse(response["chugging"]["stable"])
        self.assertEqual(response["nyquist"]["right_half_plane_zeros"], 2)

    def test_vaporization_route(self) -> None:
        response = self.responses["propellant-vaporization"]
        result = response["result"]
        for key in ("oxidizer_percent_vaporized", "fuel_percent_vaporized"):
            self.assertGreater(result[key], 0.0)
            self.assertLessEqual(result[key], 100.0)
        self.assertLess(result["cstar_efficiency"], 1.0)
        self.assertAlmostEqual(response["properties"]["oxidizer"]["critical_temperature_k"], 154.581, places=3)
        sweep = [row["oxidizer_percent_vaporized"] for row in response["sweep"]]
        self.assertEqual(sweep, sorted(sweep))
        hydrogen = dispatch_calculation("/api/v1/propellant-vaporization", {"propellant_pair": "LOX_LH2"})
        self.assertIsNone(hydrogen["properties"]["fuel"])
        self.assertEqual(hydrogen["result"]["fuel_percent_vaporized"], 99.99)

    def test_invalid_inputs_are_refused(self) -> None:
        with self.assertRaises((InputError, DomainError)):
            dispatch_calculation("/api/v1/propellant-vaporization", {"propellant_pair": "LOX_RP1"})
        with self.assertRaises((InputError, DomainError)):
            dispatch_calculation("/api/v1/combustion-stability", {"contraction_ratio": 1.0})
        with self.assertRaises((InputError, DomainError)):
            dispatch_calculation("/api/v1/transpiration-cooling", {"target_wall_temperature_k": 5000})


if __name__ == "__main__":
    unittest.main()
