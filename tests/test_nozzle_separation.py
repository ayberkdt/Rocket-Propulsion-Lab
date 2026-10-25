"""Schmucker's separation criterion (NASA TM-77396, 1984, Eq. (11), printed p. 38).

``p_i/p_a = (1.88 M_i - 1)^(-0.64)`` with the isentropic wall pressure
``p_i = p_c (1 + (gamma-1) M_i^2/2)^(-gamma/(gamma-1))``. The tests re-solve
the criterion with SciPy ``brentq``, check that the reported point is the
*first* one where the wall pressure meets it, and check the trends the
physics requires: a higher back pressure or a lower chamber pressure moves
separation upstream; vacuum never separates; the data range ends at M = 5.
"""

from __future__ import annotations

import math
import unittest

from scipy import optimize

from rocket_propulsion.compressible.isentropic import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.compressible.off_design_nozzle import (
    SCHMUCKER_MAXIMUM_MACH,
    calculate_off_design_nozzle,
    schmucker_separation_point,
)
from rocket_propulsion.core.errors import DomainError


def wall_to_ambient(mach: float, pc: float, pa: float, gamma: float) -> float:
    return pc / pa * (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (-gamma / (gamma - 1.0))


def criterion(mach: float) -> float:
    return (1.88 * mach - 1.0) ** -0.64


class CriterionAgainstBrentq(unittest.TestCase):
    def test_point_satisfies_the_criterion_and_is_the_first_crossing(self) -> None:
        for gamma in (1.15, 1.22, 1.3, 1.4):
            for pc, pa in ((7.0e6, 101325.0), (2.0e7, 101325.0), (3.0e6, 101325.0), (1.0e6, 2.0e4)):
                mach, within = schmucker_separation_point(chamber_pressure_pa=pc, ambient_pressure_pa=pa,
                                                          exit_mach=SCHMUCKER_MAXIMUM_MACH, gamma=gamma)
                with self.subTest(gamma=gamma, pc=pc, pa=pa):
                    self.assertTrue(within)
                    if mach is None:
                        # attached all the way to M = 5
                        self.assertGreater(wall_to_ambient(5.0, pc, pa, gamma), criterion(5.0))
                        continue
                    def gap(m: float, pc=pc, pa=pa, gamma=gamma) -> float:
                        return math.log(wall_to_ambient(m, pc, pa, gamma)) - math.log(criterion(m))
                    self.assertLess(abs(gap(mach)), 1e-10)
                    reference = optimize.brentq(gap, max(1.0, mach - 0.02), mach + 0.02, xtol=1e-14)
                    self.assertLess(abs(mach - reference), 1e-10 * mach)
                    # attached everywhere upstream
                    for k in range(1, 400):
                        m = 1.0 + (mach - 1.0) * k / 400
                        self.assertGreater(gap(m), 0.0)

    def test_back_pressure_and_chamber_pressure_trends(self) -> None:
        previous = math.inf
        for pa in (2.0e4, 4.0e4, 6.0e4, 8.0e4, 1.0e5):
            mach, _ = schmucker_separation_point(chamber_pressure_pa=5.0e6, ambient_pressure_pa=pa,
                                                 exit_mach=5.0, gamma=1.2)
            self.assertIsNotNone(mach)
            assert mach is not None
            self.assertLess(mach, previous)
            previous = mach
        previous = 0.0
        for pc in (2.0e6, 4.0e6, 8.0e6, 1.6e7):
            mach, _ = schmucker_separation_point(chamber_pressure_pa=pc, ambient_pressure_pa=101325.0,
                                                 exit_mach=5.0, gamma=1.2)
            assert mach is not None
            self.assertGreater(mach, previous)
            previous = mach

    def test_vacuum_range_end_and_low_chamber_pressure(self) -> None:
        self.assertEqual(schmucker_separation_point(chamber_pressure_pa=7e6, ambient_pressure_pa=0.0,
                                                    exit_mach=8.0, gamma=1.2), (None, True))
        # nearly vacuum, exit beyond the data: not assessed past M = 5
        self.assertEqual(schmucker_separation_point(chamber_pressure_pa=7e6, ambient_pressure_pa=1.0,
                                                    exit_mach=8.0, gamma=1.2), (None, False))
        # back pressure above the throat pressure criterion: separation at the throat
        mach, within = schmucker_separation_point(chamber_pressure_pa=1.0e5, ambient_pressure_pa=9.0e4,
                                                  exit_mach=3.0, gamma=1.2)
        self.assertEqual((mach, within), (1.0, True))

    def test_hostile_inputs_are_refused(self) -> None:
        good = {"chamber_pressure_pa": 7e6, "ambient_pressure_pa": 1e5, "exit_mach": 4.0, "gamma": 1.2}
        for field, bad in (("chamber_pressure_pa", 0.0), ("chamber_pressure_pa", math.nan),
                           ("ambient_pressure_pa", -1.0), ("ambient_pressure_pa", math.inf),
                           ("exit_mach", 0.9), ("exit_mach", math.nan), ("gamma", 1.0), ("gamma", math.inf)):
            with self.subTest(field=field, value=bad), self.assertRaises(DomainError):
                schmucker_separation_point(**{**good, field: bad})


class OffDesignIntegration(unittest.TestCase):
    def test_separation_point_lies_on_the_reported_wall_profile(self) -> None:
        gamma, pc, pa = 1.22, 7.0e6, 101325.0
        result = calculate_off_design_nozzle(chamber_pressure_pa=pc, chamber_temperature_k=3500.0,
                                             throat_area_m2=0.01, area_ratio=40.0, ambient_pressure_pa=pa,
                                             gamma=gamma)
        sep = result.separation
        self.assertTrue(sep.schmucker_predicted and sep.predicted)
        assert sep.schmucker_separation_mach is not None and sep.schmucker_separation_area_ratio is not None
        assert sep.schmucker_separation_pressure_pa is not None
        self.assertLess(1.0, sep.schmucker_separation_area_ratio)
        self.assertLess(sep.schmucker_separation_area_ratio, 40.0)
        back = mach_from_isentropic("area_ratio_supersonic", sep.schmucker_separation_area_ratio, gamma)
        self.assertLess(abs(back - sep.schmucker_separation_mach), 1e-10)
        expected = pc * calculate_isentropic(sep.schmucker_separation_mach, gamma).pressure_total_ratio
        self.assertLess(abs(sep.schmucker_separation_pressure_pa - expected), 1e-9 * expected)
        self.assertLess(abs(sep.schmucker_separation_pressure_pa / pa - criterion(sep.schmucker_separation_mach)),
                        1e-9)
        self.assertTrue(any(w.code == "possible_flow_separation" and "Schmucker" in w.message
                            for w in result.warnings))

    def test_full_flowing_and_unchoked_cases(self) -> None:
        attached = calculate_off_design_nozzle(chamber_pressure_pa=2.0e7, chamber_temperature_k=3500.0,
                                               throat_area_m2=0.01, area_ratio=8.0, ambient_pressure_pa=101325.0)
        self.assertFalse(attached.separation.schmucker_predicted)
        self.assertIsNone(attached.separation.schmucker_separation_area_ratio)
        subsonic = calculate_off_design_nozzle(chamber_pressure_pa=1.05e5, chamber_temperature_k=300.0,
                                               throat_area_m2=0.01, area_ratio=4.0, ambient_pressure_pa=1.04e5)
        self.assertIn(subsonic.regime, {"unchoked", "just-choked"})
        self.assertFalse(subsonic.separation.schmucker_predicted)


if __name__ == "__main__":
    unittest.main()
