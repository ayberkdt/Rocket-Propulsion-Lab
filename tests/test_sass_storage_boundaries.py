"""Storage attacks with independently solved ideal-gas endpoint inventories.

References: NASA-GAS ideal equation of state; NASA-SP8112 inventory boundaries;
NASA-7009B numerical verification only. Extreme scales are synthetic, not gas
property extrapolation or pressure-vessel design recommendations.
"""

import math
import unittest
from decimal import Context, Decimal, localcontext

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import (
    ideal_cold_gas_storage, ideal_regulated_pressurant, liquid_tank_volume,
)
from tests.test_sass_studies import q


def bottle_root(target, p0, pf, temperature, gas_r):
    """Solve endpoint gas inventory difference for V, not a closed volume formula."""
    with localcontext(Context(prec=115)):
        m, high, low, t, r = map(Decimal, (target, p0, pf, temperature, gas_r))
        def delivered(volume):
            return high * volume / (r * t) - low * volume / (r * t)
        upper = Decimal(1)
        for _ in range(4096):
            if delivered(upper) >= m:
                break
            upper *= 2
        else:
            raise AssertionError("Independent inventory root could not bracket upward.")
        for _ in range(4096):
            if delivered(upper / 2) < m:
                break
            upper /= 2
        else:
            raise AssertionError("Independent inventory root could not bracket downward.")
        lower = Decimal(0)
        for _ in range(360):
            middle = (lower + upper) / 2
            if delivered(middle) < m:
                lower = middle
            else:
                upper = middle
        volume = (lower + upper) / 2
        return tuple(map(float, (volume, high * volume / (r * t), low * volume / (r * t))))


class StorageBoundaryAttacks(unittest.TestCase):
    def test_cold_gas_volume_survives_overflow_and_underflow_of_intermediate_product(self):
        for mass, high, low, temperature, gas_r in (
                (1e200, 2e200, 1e200, 1., 1e200),
                (1e-200, 2e-200, 1e-200, 1., 1e-200),
                (1e-200, 2e200, 1e200, 1e200, 1e200)):
            expected = bottle_root(mass, high, low, temperature, gas_r)
            with self.subTest(inputs=(mass, high, low, temperature, gas_r)):
                result = ideal_cold_gas_storage(q(mass, "kg"), q(high, "Pa"), q(low, "Pa"),
                                               q(temperature, "K"), q(gas_r, "J/(kg*K)"))
                for name, reference in zip(("tank_volume", "loaded_gas", "remaining_gas"), expected):
                    self.assertTrue(math.isclose(result.get(name).value, reference,
                                                 rel_tol=3e-14, abs_tol=0), name)

    def test_regulated_transfer_survives_underflow_of_density_before_displacement(self):
        displaced, pt, high, low, temperature, gas_r = 1e200, 1e-200, 4e-200, 2e-200, 1., 1e200
        with localcontext(Context(prec=115)):
            target = Decimal(pt) * Decimal(displaced) / (Decimal(temperature) * Decimal(gas_r))
        volume, loaded, remaining = bottle_root(target, high, low, temperature, gas_r)
        result = ideal_regulated_pressurant(q(displaced, "m^3"), q(pt, "Pa"), q(high, "Pa"),
                                           q(low, "Pa"), q(temperature, "K"), q(gas_r, "J/(kg*K)"))
        for name, expected in (("bottle_volume", volume), ("loaded_pressurant", loaded),
                               ("remaining_pressurant", remaining), ("delivered_pressurant", float(target))):
            self.assertTrue(math.isclose(result.get(name).value, expected, rel_tol=3e-14, abs_tol=0), name)

    def test_regulated_transfer_survives_overflow_of_density_before_displacement(self):
        displaced, pt, high, low, temperature, gas_r = 1e-200, 1e200, 4e200, 2e200, 1., 1e-200
        with localcontext(Context(prec=115)):
            target = Decimal(pt) * Decimal(displaced) / (Decimal(temperature) * Decimal(gas_r))
        volume, loaded, remaining = bottle_root(target, high, low, temperature, gas_r)
        result = ideal_regulated_pressurant(q(displaced, "m^3"), q(pt, "Pa"), q(high, "Pa"),
                                           q(low, "Pa"), q(temperature, "K"), q(gas_r, "J/(kg*K)"))
        for name, expected in (("bottle_volume", volume), ("loaded_pressurant", loaded),
                               ("remaining_pressurant", remaining), ("delivered_pressurant", float(target))):
            self.assertTrue(math.isclose(result.get(name).value, expected, rel_tol=3e-14, abs_tol=0), name)

    def test_positive_unrepresentable_ullage_cannot_be_reported_as_empty(self):
        with self.assertRaises(DomainError):
            liquid_tank_volume(q(math.nextafter(0., 1.), "kg"), q(1, "kg/m^3"), q(.1))

    def test_zero_inventory_stays_zero_and_unrepresentable_final_volumes_are_rejected(self):
        result = ideal_cold_gas_storage(q(0, "kg"), q(2, "Pa"), q(1, "Pa"),
                                       q(1e200, "K"), q(1e200, "J/(kg*K)"))
        self.assertTrue(all(value.value == 0 for _, value in result.values))
        with self.assertRaises(DomainError):
            ideal_cold_gas_storage(q(1e200, "kg"), q(2, "Pa"), q(1, "Pa"),
                                   q(1, "K"), q(1e200, "J/(kg*K)"))


if __name__ == "__main__":
    unittest.main()
