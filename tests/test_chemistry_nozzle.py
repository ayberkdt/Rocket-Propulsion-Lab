"""Off-design nozzle regimes calibrated to shifting-equilibrium chemistry."""

from __future__ import annotations

import unittest
from functools import cache

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.thermochemistry.chemistry_nozzle import solve_chemistry_coupled_nozzle


@cache
def _solve(ambient: float, sweep: bool = False):
    return solve_chemistry_coupled_nozzle(
        propellant_pair="LOX_LH2", chamber_pressure_pa=7.0e6, area_ratio=40.0,
        ambient_pressure_pa=ambient, throat_area_m2=0.01, include_altitude_sweep=sweep,
        point_count=11,
    )


class ChemistryCoupledNozzleTests(unittest.TestCase):
    def test_mass_flow_and_design_exit_pressure_match_the_equilibrium_solution(self) -> None:
        result = _solve(0.0)
        equilibrium = result.equilibrium
        expected_flow = 7.0e6 * 0.01 / equilibrium.characteristic_velocity_m_s
        self.assertAlmostEqual(result.off_design.mass_flow_kg_s / expected_flow, 1.0, delta=1.0e-9)
        self.assertAlmostEqual(
            result.off_design.exit_pressure_pa / equilibrium.exits[0].station.pressure_pa, 1.0, delta=1.0e-9
        )
        self.assertEqual(result.off_design.regime, "underexpanded")

    def test_calibration_residual_is_reported_and_small(self) -> None:
        gas = _solve(0.0).equivalent_gas
        self.assertGreater(gas.gamma, 1.05)
        self.assertLess(gas.gamma, 1.67)
        self.assertLess(abs(gas.vacuum_isp_relative_residual), 0.02)
        self.assertAlmostEqual(
            gas.equivalent_vacuum_specific_impulse_s / gas.equilibrium_vacuum_specific_impulse_s - 1.0,
            gas.vacuum_isp_relative_residual, delta=1.0e-12,
        )

    def test_sea_level_back_pressure_overexpands_and_screens_separation(self) -> None:
        result = _solve(101_325.0)
        self.assertEqual(result.off_design.regime, "overexpanded")
        self.assertTrue(result.off_design.separation.predicted)
        self.assertTrue(any(warning.code == "equivalent_perfect_gas" for warning in result.warnings))

    def test_altitude_sweep_uses_the_same_calibrated_gas(self) -> None:
        result = _solve(101_325.0, True)
        sweep = result.altitude_sweep
        self.assertIsNotNone(sweep)
        thrusts = [point.thrust_n for point in sweep.points]
        self.assertEqual(thrusts, sorted(thrusts))

    def test_invalid_inputs_are_rejected(self) -> None:
        with self.assertRaises(InputError):
            solve_chemistry_coupled_nozzle(propellant_pair="LOX_XYZ", chamber_pressure_pa=7.0e6,
                                           area_ratio=40.0, ambient_pressure_pa=0.0)
        with self.assertRaises(DomainError):
            solve_chemistry_coupled_nozzle(propellant_pair="LOX_LH2", chamber_pressure_pa=7.0e6,
                                           area_ratio=0.5, ambient_pressure_pa=0.0)


if __name__ == "__main__":
    unittest.main()
