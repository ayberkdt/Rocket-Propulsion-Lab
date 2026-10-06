"""Reference and inverse tests for isentropic flow."""

import unittest

from rocket_propulsion.compressible.isentropic import (
    calculate_isentropic,
    mach_from_isentropic,
)
from rocket_propulsion.core.errors import DomainError


class IsentropicTests(unittest.TestCase):
    def test_mach_two_reference_case(self) -> None:
        result = calculate_isentropic(2.0, 1.4)

        self.assertAlmostEqual(result.mach_angle_deg or 0.0, 30.0, places=10)
        self.assertAlmostEqual(result.prandtl_meyer_deg or 0.0, 26.3797608134, places=8)
        self.assertAlmostEqual(result.pressure_total_ratio, 0.1278045255, places=9)
        self.assertAlmostEqual(result.density_total_ratio, 0.2300481458, places=9)
        self.assertAlmostEqual(result.temperature_total_ratio, 5.0 / 9.0, places=10)
        self.assertAlmostEqual(result.area_critical_ratio, 1.6875, places=10)

    def test_subsonic_angles_are_not_defined(self) -> None:
        result = calculate_isentropic(0.5)
        self.assertIsNone(result.mach_angle_deg)
        self.assertIsNone(result.prandtl_meyer_deg)

    def test_inverse_properties_recover_mach(self) -> None:
        result = calculate_isentropic(2.5)
        cases = {
            "mach_angle_deg": result.mach_angle_deg,
            "prandtl_meyer_deg": result.prandtl_meyer_deg,
            "pressure_total_ratio": result.pressure_total_ratio,
            "density_total_ratio": result.density_total_ratio,
            "temperature_total_ratio": result.temperature_total_ratio,
            "pressure_critical_ratio": result.pressure_critical_ratio,
            "density_critical_ratio": result.density_critical_ratio,
            "temperature_critical_ratio": result.temperature_critical_ratio,
            "area_ratio_supersonic": result.area_critical_ratio,
        }
        for kind, value in cases.items():
            with self.subTest(kind=kind):
                assert value is not None
                self.assertAlmostEqual(mach_from_isentropic(kind, value), 2.5, places=7)

    def test_area_ratio_branch_is_explicit(self) -> None:
        area_ratio = 1.6875
        self.assertAlmostEqual(
            mach_from_isentropic("area_ratio_supersonic", area_ratio), 2.0, places=8
        )
        self.assertLess(mach_from_isentropic("area_ratio_subsonic", area_ratio), 1.0)

    def test_rejects_nonphysical_gamma_and_mach(self) -> None:
        with self.assertRaises(DomainError):
            calculate_isentropic(0.0)
        with self.assertRaises(DomainError):
            calculate_isentropic(1.0, gamma=1.0)


if __name__ == "__main__":
    unittest.main()

