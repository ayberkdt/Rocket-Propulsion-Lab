"""HTTP routes for CEA-equivalent rocket performance and chemistry-driven nozzles."""

from __future__ import annotations

import json
import unittest
from unittest.mock import patch

from rocket_propulsion.api.routes import ROUTE_METADATA, dispatch_calculation
from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.thermochemistry import cea_library_info


class RocketPerformanceRouteTests(unittest.TestCase):
    def test_default_request_returns_sorted_exits_and_cea_status(self) -> None:
        response = dispatch_calculation("/api/v1/rocket-performance", {})
        result = response["result"]
        self.assertEqual(response["pair"]["key"], "LOX_LH2")
        self.assertEqual([item["area_ratio"] for item in result["exits"]], [5.0, 10.0, 25.0, 40.0, 80.0])
        isp = [item["vacuum_specific_impulse_s"] for item in result["exits"]]
        self.assertEqual(isp, sorted(isp))
        self.assertIn("available", response["cea"])
        self.assertIn("/api/v1/rocket-performance", ROUTE_METADATA)

    def test_response_is_strict_json(self) -> None:
        response = dispatch_calculation("/api/v1/rocket-performance", {"area_ratios": [10]})
        json.dumps(response, allow_nan=False)
        self.assertIsNone(response["result"]["chamber"]["area_ratio"])

    def test_frozen_chemistry_gives_lower_isp_than_shifting(self) -> None:
        payload = {"propellant_pair": "LOX_LH2", "area_ratios": [40]}
        shifting = dispatch_calculation("/api/v1/rocket-performance", payload)["result"]
        frozen = dispatch_calculation("/api/v1/rocket-performance", {**payload, "chemistry": "frozen"})["result"]
        self.assertLess(frozen["exits"][0]["vacuum_specific_impulse_s"],
                        shifting["exits"][0]["vacuum_specific_impulse_s"])

    def test_invalid_requests_are_refused(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/rocket-performance", {"provider": "other"})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/rocket-performance", {"area_ratios": list(range(2, 12))})
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/rocket-performance", {"area_ratios": [0.5]})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/rocket-performance", {"chamber_pressure_pa": True})

    def test_cea_provider_without_package_is_a_feature_error(self) -> None:
        cea_library_info.cache_clear()
        try:
            with patch('rocket_propulsion.thermochemistry.cea_library._import_cea',
                       side_effect=ModuleNotFoundError('deliberately unavailable CEA')):
                self.assertFalse(cea_library_info().available)
                with self.assertRaises(FeatureUnavailableError):
                    dispatch_calculation("/api/v1/rocket-performance", {"provider": "compare"})
        finally:
            cea_library_info.cache_clear()


class ChemistryNozzleRouteTests(unittest.TestCase):
    def test_default_request_is_a_sea_level_rp1_booster(self) -> None:
        response = dispatch_calculation("/api/v1/nozzle-chemistry", {"point_count": 6})
        self.assertEqual(response["propellant_pair"], "LOX_RP1")
        self.assertEqual(response["off_design"]["regime"], "overexpanded")
        self.assertEqual(len(response["altitude_sweep"]["points"]), 6)
        self.assertLess(abs(response["equivalent_gas"]["vacuum_isp_relative_residual"]), 0.02)

    def test_fractional_point_count_is_refused(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/nozzle-chemistry", {"point_count": 6.5})


if __name__ == "__main__":
    unittest.main()
