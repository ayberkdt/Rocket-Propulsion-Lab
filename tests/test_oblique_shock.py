"""Reference tests for attached oblique shocks."""

import unittest

from rocket_propulsion.compressible.oblique_shock import calculate_oblique_shock
from rocket_propulsion.core.errors import DomainError


class ObliqueShockTests(unittest.TestCase):
    def test_mach_two_ten_degree_weak_solution(self) -> None:
        result = calculate_oblique_shock(
            upstream_mach=2.0,
            turning_angle_deg=10.0,
            gamma=1.4,
            branch="weak",
        )
        self.assertAlmostEqual(result.wave_angle_deg, 39.3139, places=3)
        self.assertAlmostEqual(result.downstream_mach, 1.6405, places=3)
        self.assertAlmostEqual(result.pressure_ratio, 1.7066, places=3)
        self.assertGreater(result.maximum_turning_angle_deg, 22.0)

    def test_strong_branch_is_subsonic_downstream(self) -> None:
        result = calculate_oblique_shock(
            upstream_mach=2.0,
            turning_angle_deg=10.0,
            branch="strong",
        )
        self.assertGreater(result.wave_angle_deg, 80.0)
        self.assertLess(result.downstream_mach, 1.0)

    def test_detached_request_is_rejected_with_limit(self) -> None:
        with self.assertRaisesRegex(DomainError, "Detached shock"):
            calculate_oblique_shock(upstream_mach=2.0, turning_angle_deg=30.0)


if __name__ == "__main__":
    unittest.main()

