"""Independent n=1 blowdown curves plus bounded repeated-evaluation controls.

For F=F0*p/p0 and constant Isp, dt/dx=(g0*Isp/F0)*(1+x/(rho*V0)).
The analytic primitive checks the entire sampled time history. These are
synthetic equation checks (NASA SP-8112), not hardware-validation evidence.
"""

from __future__ import annotations

import math
import unittest
from dataclasses import dataclass
from unittest.mock import patch

from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    simulate_pressure_fed_blowdown,
)


@dataclass(frozen=True)
class ForwardingLaw:
    base: PressurePerformanceLaw

    @property
    def minimum_supply_pressure_pa(self):
        return self.base.minimum_supply_pressure_pa

    @property
    def maximum_pressure_pa(self):
        return self.base.maximum_pressure_pa

    def scale_operating_point(self, point, pressure):
        return self.base.scale_operating_point(point, pressure)


class DerivedLaw(PressurePerformanceLaw):
    """A user subclass must retain custom-provider evaluation semantics."""


class BlowdownEvaluationReuseTests(unittest.TestCase):
    def simulate(self, thrust, exponent, law=None):
        tank = BlowdownTankModel(2e6, .01, 10., 1000., 1.)
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=thrust, system_specific_impulse_s=220.,
            chamber_pressure_pa=1.5e6, source='analytic evaluation-reuse fixture').operating_point()
        if law is None:
            law = PressurePerformanceLaw(2e6, 1e6, thrust_pressure_exponent=exponent,
                                         chamber_pressure_exponent=1.)
        return simulate_pressure_fed_blowdown(tank, point, law, reserve_mass_kg=2.,
            maximum_mass_step_kg=2., relative_tolerance=1e-5)

    def assert_analytic_curve(self, result, thrust, exponent):
        ce = 9.80665 * 220.
        for sample in result.samples:
            expelled = 10. - sample.remaining_propellant_mass_kg
            expected_pressure = 2e6 * .01 / (.01 + expelled/1000.)
            expected_thrust = thrust * (expected_pressure/2e6)**exponent
            primitive = expelled + (expelled**2/20. if exponent == 1. else 0.)
            self.assertTrue(math.isclose(sample.time_s,ce/thrust*primitive,
                                         rel_tol=3e-12,abs_tol=1e-12))
            self.assertTrue(math.isclose(sample.supply_pressure_pa,expected_pressure,rel_tol=3e-14))
            self.assertTrue(math.isclose(sample.delivered_thrust_n,expected_thrust,rel_tol=3e-14))
        self.assertTrue(math.isclose(result.evidence.delivered_impulse_n_s,ce*8.,rel_tol=3e-13))
        self.assertTrue(math.isclose(result.artifact.consumed_propellant_kg,8.,rel_tol=1e-5))

    def test_builtin_law_avoids_duplicate_solver_work_and_keeps_analytic_history(self):
        calls = []
        original = PressurePerformanceLaw.scale_operating_point
        def measured(model, point, pressure):
            calls.append(pressure)
            return original(model,point,pressure)
        with patch.object(PressurePerformanceLaw,'scale_operating_point',measured):
            result = self.simulate(100.,1.)
        self.assert_analytic_curve(result,100.,1.)
        self.assertGreater(len(calls), 3)
        self.assertEqual(len(calls),len(set(calls)), 'Repeated built-in pressure queries redo validation.')

    def test_custom_provider_keeps_its_evaluation_behavior_and_analytic_history(self):
        for law in (ForwardingLaw(PressurePerformanceLaw(2e6,1e6)), DerivedLaw(2e6,1e6)):
            with self.subTest(provider=type(law).__name__):
                calls = []
                original = type(law).scale_operating_point
                def measured(model, point, pressure, observed=calls, delegate=original):
                    observed.append(pressure)
                    return delegate(model,point,pressure)
                with patch.object(type(law),'scale_operating_point',measured):
                    result = self.simulate(100.,1.,law)
                self.assert_analytic_curve(result,100.,1.)
                self.assertGreater(len(calls),len(set(calls)), 'Custom providers were implicitly memoized.')

    def test_reuse_cannot_leak_across_rated_points_or_pressure_laws(self):
        for thrust,exponent in [(100.,1.),(125.,1.),(100.,0.),(125.,0.)]:
            with self.subTest(thrust=thrust,exponent=exponent):
                self.assert_analytic_curve(self.simulate(thrust,exponent),thrust,exponent)


if __name__ == '__main__':
    unittest.main()
