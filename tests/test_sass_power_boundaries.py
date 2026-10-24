"""Independent phase-energy checks at duty and binary scale boundaries.

References: BIPM-SI phase energy; NASA-7009B numerical verification only.
"""

import math
import unittest
from decimal import Context, Decimal, localcontext

from rocket_propulsion.studies.sass import power_resources
from tests.test_sass_studies import q


class PowerBoundaryAttacks(unittest.TestCase):
    def test_near_continuous_duty_preserves_small_standby_energy(self):
        for time in (1., 1.3, 100., 1e-310):
            duty = math.nextafter(1., 0.)
            with localcontext(Context(prec=115)):
                expected = float(Decimal.from_float(1e300) * Decimal(time)
                                 * (1 - Decimal(duty)) / Decimal(duty))
            result = power_resources(q(0, "W"), q(time, "s"), q(duty),
                                     q(1e300, "W"), q(0, "W"), q(0, "s"))
            with self.subTest(time=time):
                self.assertTrue(math.isclose(result.get("energy").value, expected,
                                             rel_tol=3e-14, abs_tol=0))
                self.assertEqual(result.get("peak_bus_power").value, 1e300)

    def test_phase_energy_survives_underflow_of_individual_float_products(self):
        # Each energy term rounds to zero separately; their sum is representable.
        smallest = math.nextafter(0., 1.)
        result = power_resources(q(smallest, "W"), q(.4, "s"), q(1),
                                 q(0, "W"), q(smallest, "W"), q(.4, "s"))
        self.assertEqual(result.get("energy").value, smallest)
        self.assertEqual(result.get("average_bus_power").value, smallest)


if __name__ == "__main__":
    unittest.main()
