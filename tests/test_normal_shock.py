"""Reference tests for normal-shock flow."""

import unittest

from rocket_propulsion.compressible.normal_shock import calculate_normal_shock
from rocket_propulsion.core.errors import DomainError


class NormalShockTests(unittest.TestCase):
    def test_mach_two_reference_case(self) -> None:
        result = calculate_normal_shock(2.0, 1.4)

        self.assertAlmostEqual(result.downstream_mach, 0.5773502692, places=9)
        self.assertAlmostEqual(result.pressure_ratio, 4.5, places=10)
        self.assertAlmostEqual(result.density_ratio, 8.0 / 3.0, places=10)
        self.assertAlmostEqual(result.temperature_ratio, 1.6875, places=10)
        self.assertAlmostEqual(result.stagnation_pressure_ratio, 0.7208738615, places=9)

    def test_requires_supersonic_upstream_flow(self) -> None:
        with self.assertRaises(DomainError):
            calculate_normal_shock(1.0)


if __name__ == "__main__":
    unittest.main()

