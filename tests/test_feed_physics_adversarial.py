"""Independent feed-physics oracles and ill-conditioned boundary attacks.

Synthetic inputs verify equations and numerical implementation only. Oracles:
NASA Glenn mass-flow/choking and isentropic relations; NASA SP-8112 open-system
mass/energy conservation. No experimental or hardware-validation claim.
https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
https://ntrs.nasa.gov/citations/19760015212
"""

from __future__ import annotations

import math
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    DynamicFeedLine,
    DynamicFeedLineState,
    IdealPressurantGas,
    RegulatedFeedState,
    RegulatedFeedSystem,
)
from rocket_propulsion.propulsion.burns.regulator_hysteresis import (
    _is_choked,
    compressible_valve_mass_flow_kg_s,
)
from tests.test_burn_feed_line import example_line
from tests.test_burn_feed_system import example_feed_system


def decimal_orifice(config, upstream, downstream, temperature, opening):
    """75-digit pressure-form isentropic flow, without production helper calls."""
    with localcontext() as ctx:
        ctx.prec = 75
        gamma, gas_r = map(Decimal, (config.gas.gamma, config.gas.specific_gas_constant_j_kg_k))
        pu, pd, t, a, cd, fraction = map(Decimal, (
            upstream, downstream, temperature, config.regulator_maximum_area_m2,
            config.regulator_discharge_coefficient, opening))
        if pd >= pu or fraction == 0:
            return 0., False
        ratio = pd / pu
        critical = (Decimal(2) / (gamma + 1)) ** (gamma / (gamma - 1))
        if ratio <= critical:
            speed = gamma.sqrt() * (Decimal(2) / (gamma + 1)) ** (
                (gamma + 1) / (2 * (gamma - 1)))
        else:
            speed = (2 * gamma / (gamma - 1) * (ratio ** (2 / gamma)
                     - ratio ** ((gamma + 1) / gamma))).sqrt()
        return float(cd * a * fraction * pu * speed / (gas_r * t).sqrt()), ratio <= critical


def step_rk4(rhs, values, dt):
    """Independent classical RK4 integration of the public six-state RHS."""
    first = rhs(values)
    second = rhs(tuple(x + dt * .5 * k for x, k in zip(values, first)))
    third = rhs(tuple(x + dt * .5 * k for x, k in zip(values, second)))
    fourth = rhs(tuple(x + dt * k for x, k in zip(values, third)))
    return tuple(x + dt / 6 * (a + 2*b + 2*c + d)
                 for x, a, b, c, d in zip(values, first, second, third, fourth))


class FeedBoundaryAttacks(unittest.TestCase):
    def test_adjacent_sonic_states_match_independent_pressure_ratio(self):
        with localcontext() as context:
            context.prec = 110
            for gamma in (math.nextafter(1.0, math.inf), 1.01, 1.1, 1.4, 1.5, 1.667, 2.0, 5.0):
                gas = IdealPressurantGas(1.0, gamma, "synthetic")
                config = replace(example_feed_system().configuration, gas=gas)
                system = RegulatedFeedSystem(config)
                exact_gamma = Decimal(gamma)
                critical = (Decimal(2) / (exact_gamma + 1))**(exact_gamma / (exact_gamma - 1))
                for upstream in (1.0, 2.0, 2e6, 1e-100, 1e100):
                    center = float(Decimal(upstream) * critical)
                    for downstream in (
                        math.nextafter(center, 0.0), center, math.nextafter(center, math.inf),
                    ):
                        with self.subTest(gamma=gamma, upstream=upstream, downstream=downstream):
                            expected = Decimal(downstream) / Decimal(upstream) <= critical
                            _, dynamic_flag = system._regulator_flow(
                                upstream_pressure_pa=upstream, downstream_pressure_pa=downstream,
                                upstream_temperature_k=300.0, opening=0.5,
                            )
                            actual = (dynamic_flag, _is_choked(gas, upstream, downstream))
                            self.assertEqual(actual, (expected, expected))

    def test_exact_sonic_ratio_is_choked_and_has_distinct_adjacent_sides(self):
        # Exact fractions are (2/3)^2, (4/5)^3 and (8/9)^5.
        for gamma, upstream, center in (
            (2.0, 9.0, 4.0), (1.5, 125.0, 64.0), (1.25, 59049.0, 32768.0),
        ):
            gas = IdealPressurantGas(1.0, gamma, "synthetic")
            system = RegulatedFeedSystem(replace(example_feed_system().configuration, gas=gas))
            for downstream, expected in (
                (math.nextafter(center, 0.0), True), (center, True),
                (math.nextafter(center, math.inf), False),
            ):
                with self.subTest(gamma=gamma, downstream=downstream):
                    flow, choked = system._regulator_flow(
                        upstream_pressure_pa=upstream, downstream_pressure_pa=downstream,
                        upstream_temperature_k=300.0, opening=0.5,
                    )
                    self.assertGreater(flow, 0.0)
                    self.assertEqual((choked, _is_choked(gas, upstream, downstream)),
                                     (expected, expected))

    def test_direct_valve_flow_does_not_vanish_at_one_ulp_forward_head(self):
        config = replace(
            example_feed_system().configuration,
            gas=IdealPressurantGas(300.0, math.nextafter(1.0, math.inf), "synthetic"),
            regulator_discharge_coefficient=1.0,
        )
        upstream = 2e6
        downstream = math.nextafter(upstream, 0.0)
        expected, _ = decimal_orifice(config, upstream, downstream, 300.0, 1.0)
        actual = compressible_valve_mass_flow_kg_s(
            gas=config.gas, effective_area_m2=config.regulator_maximum_area_m2,
            upstream_pressure_pa=upstream, downstream_pressure_pa=downstream,
            upstream_temperature_k=300.0,
        )
        self.assertGreater(actual, 0.0)
        self.assertTrue(math.isclose(actual, expected, rel_tol=8e-14))

    def test_sonic_regime_resolves_pressure_ratios_closer_than_initial_precision(self):
        gamma = 1.4
        gas = IdealPressurantGas(1.0, gamma, "synthetic")
        system = RegulatedFeedSystem(replace(example_feed_system().configuration, gas=gas))
        with localcontext() as context:
            context.prec = 380
            g = Decimal(gamma)
            critical = (Decimal(2) / (1 + g))**(g / (g - 1))
            scale = 10**160
            floor = int(critical * scale)
            self.assertGreater(critical * scale - floor, 0)
        for offset, expected in ((0, True), (1, False)):
            ratio = Fraction(floor + offset, scale)
            with self.subTest(offset=offset):
                flow, choked = system._regulator_flow(
                    upstream_pressure_pa=1.0, downstream_pressure_pa=ratio,
                    upstream_temperature_k=300.0, opening=0.5,
                )
                self.assertGreater(flow, 0.0)
                self.assertEqual(choked, expected)

    def test_exact_integer_power_sonic_equality_after_first_bound(self):
        gamma = 1.015625
        gas = IdealPressurantGas(1.0, gamma, "synthetic")
        system = RegulatedFeedSystem(replace(example_feed_system().configuration, gas=gas))
        exact_limit = Fraction(128, 129)**65
        try:
            flow, choked = system._regulator_flow(
                upstream_pressure_pa=1.0, downstream_pressure_pa=exact_limit,
                upstream_temperature_k=300.0, opening=0.5,
            )
        except DomainError as error:
            raise AssertionError("Exact equality was lost in rounded interval bounds.") from error
        self.assertGreater(flow, 0.0)
        self.assertTrue(choked)

    def test_unresolved_sonic_side_is_refused_instead_of_guessed(self):
        gamma = 1.4
        gas = IdealPressurantGas(1.0, gamma, "synthetic")
        system = RegulatedFeedSystem(replace(example_feed_system().configuration, gas=gas))
        with localcontext() as context:
            context.prec = 1600
            g = Decimal(gamma)
            near_sonic = Fraction((Decimal(2) / (1 + g))**(g / (g - 1)))
        try:
            system._regulator_flow(
                upstream_pressure_pa=1.0, downstream_pressure_pa=near_sonic,
                upstream_temperature_k=300.0, opening=0.5,
            )
        except DomainError as error:
            self.assertEqual(error.field, "Sonic pressure ratio")
        else:
            raise AssertionError("An unresolved sonic-side decision was published.")

    def test_equal_display_pressures_still_resolve_forward_head(self):
        config = replace(
            example_feed_system().configuration,
            gas=IdealPressurantGas(1.0, 1.4, "synthetic"),
            tank_internal_volume_m3=2.0, propellant_density_kg_m3=1.0,
            bottle_volume_m3=1.0,
        )
        system = RegulatedFeedSystem(config)
        gas_mass = math.nextafter(1.0, math.inf)
        gas_temperature = math.nextafter(2.0, 0.0)
        state = RegulatedFeedState(1.0, 1.0, 2.0, gas_mass, gas_temperature, 0.5)
        with localcontext() as context:
            context.prec = 100
            upstream = Decimal(gas_mass) * Decimal(gas_temperature)
            self.assertGreater(upstream, 2)
            expected, _ = decimal_orifice(config, upstream, 2.0, gas_temperature, 0.5)
        result = system.evaluate(state, 0.0)
        self.assertEqual(result.tank_pressure_pa, result.bottle_pressure_pa)
        self.assertGreater(result.regulator_mass_flow_kg_s, 0.0)
        self.assertTrue(math.isclose(result.regulator_mass_flow_kg_s, expected, rel_tol=8e-14))
        self.assertFalse(result.regulator_choked)
        self.assertEqual(result.derivative.tank_pressurant_mass_kg_s, result.regulator_mass_flow_kg_s)
        self.assertEqual(result.derivative.bottle_pressurant_mass_kg_s, -result.regulator_mass_flow_kg_s)
        reverse = replace(state, tank_pressurant_mass_kg=gas_mass, tank_gas_temperature_k=gas_temperature,
                          bottle_pressurant_mass_kg=1.0, bottle_gas_temperature_k=2.0)
        equal = replace(reverse, tank_pressurant_mass_kg=1.0, tank_gas_temperature_k=2.0)
        self.assertEqual(system.evaluate(reverse, 0.0).regulator_mass_flow_kg_s, 0.0)
        self.assertEqual(system.evaluate(equal, 0.0).regulator_mass_flow_kg_s, 0.0)

    def test_valve_motion_uses_unrounded_pressure_error_and_command(self):
        config = replace(
            example_feed_system().configuration,
            gas=IdealPressurantGas(1.0, 1.4, "synthetic"),
            tank_internal_volume_m3=2.0, propellant_density_kg_m3=1.0,
            regulator_set_pressure_pa=2.5, regulator_full_open_error_pa=1.0,
            valve_time_constant_s=0.5,
        )
        mass = math.nextafter(1.0, math.inf)
        temperature = math.nextafter(math.nextafter(2.0, 0.0), 0.0)
        state = RegulatedFeedState(1.0, mass, temperature, 1.0, 2.0, 0.5)
        result = RegulatedFeedSystem(config).evaluate(state, 0.0)
        with localcontext() as context:
            context.prec = 100
            pressure = Decimal(mass) * Decimal(temperature)
            expected = (Decimal("2.5") - pressure - Decimal("0.5")) / Decimal("0.5")
        self.assertGreater(expected, 0.0)
        self.assertEqual(result.tank_pressure_pa, 2.0)
        self.assertEqual(result.regulator_command, 0.5)
        self.assertEqual(result.derivative.regulator_opening_s, float(expected))

    def test_regulator_flow_survives_compensated_area_and_pressure_scales(self):
        for area, opening, pressure, gas_r, temperature in (
            (1e300, 1.0, 1e300, 1e307, 1e307),
            (1e-300, 1e-100, 1e200, 300.0, 300.0),
        ):
            with self.subTest(area=area, opening=opening):
                config = replace(
                    example_feed_system().configuration,
                    gas=IdealPressurantGas(gas_r, 1.4, "synthetic"),
                    regulator_maximum_area_m2=area,
                )
                expected, choked = decimal_orifice(config, pressure, pressure / 2, temperature, opening)
                self.assertGreater(expected, 0.0)
                try:
                    flow, actual_choked = RegulatedFeedSystem(config)._regulator_flow(
                        upstream_pressure_pa=pressure, downstream_pressure_pa=pressure / 2,
                        upstream_temperature_k=temperature, opening=opening,
                    )
                except DomainError as error:
                    raise AssertionError("A representable flow was lost in an intermediate.") from error
                self.assertEqual(actual_choked, choked)
                self.assertTrue(math.isclose(flow, expected, rel_tol=8e-14))

    def test_regulator_flow_preserves_head_smaller_than_smallest_binary64(self):
        config = example_feed_system().configuration
        head = Fraction(1, 1 << 1200)
        with localcontext() as context:
            context.prec = 100
            # The relative correction to the square-root head limit is O(head).
            h = Decimal(1) / Decimal(head.denominator)
            expected = (Decimal(config.regulator_discharge_coefficient)
                        * Decimal(config.regulator_maximum_area_m2)
                        * (2 * h / (Decimal(config.gas.specific_gas_constant_j_kg_k)
                                    * Decimal(300))).sqrt())
        try:
            flow, choked = RegulatedFeedSystem(config)._regulator_flow(
                upstream_pressure_pa=1.0, downstream_pressure_pa=1 - head,
                upstream_temperature_k=300.0, opening=1.0,
            )
        except DomainError as error:
            raise AssertionError("A representable square-root flow was refused.") from error
        self.assertFalse(choked)
        self.assertTrue(math.isclose(flow, float(expected), rel_tol=8e-14))

    def test_near_equilibrium_regulator_flow_against_decimal(self):
        for gamma in (math.nextafter(1., math.inf), 1.000001, 1.1, 1.4, 1.667):
            config = replace(example_feed_system().configuration,
                             gas=IdealPressurantGas(300., gamma, "synthetic"))
            system = RegulatedFeedSystem(config)
            upstream = 2_000_000.
            for ratio in (math.nextafter(1., 0.), .999999999, .999, .8, .6, .3):
                downstream = upstream * ratio
                expected, choked = decimal_orifice(config, upstream, downstream, 300., .5)
                actual, actual_choked = system._regulator_flow(
                    upstream_pressure_pa=upstream, downstream_pressure_pa=downstream,
                    upstream_temperature_k=300., opening=.5)
                with self.subTest(gamma=gamma, ratio=ratio):
                    self.assertGreater(actual, 0)
                    self.assertTrue(math.isclose(actual, expected, rel_tol=8e-14, abs_tol=0),
                                    (actual, expected))
                    self.assertEqual(actual_choked, choked)

    def test_gas_thermodynamic_properties_cannot_be_infinite(self):
        with self.assertRaises(DomainError):
            IdealPressurantGas(1e308, 1.01)

    def test_finite_modal_frequency_survives_extreme_line_scaling(self):
        for scale in (1e-200, 1e200):
            line = DynamicFeedLine(replace(example_line().configuration,
                inertance_pa_s2_kg=scale, compliance_kg_pa=scale,
                resistance_pa_s2_kg2=0.))
            result = line.evaluate(DynamicFeedLineState(0., 2.),
                                   upstream_pressure_pa=2., engine_mass_flow_kg_s=0.)
            with self.subTest(scale=scale):
                self.assertTrue(math.isclose(result.natural_angular_frequency_rad_s,
                                             1 / scale, rel_tol=5e-15))

    def test_line_cannot_return_nonfinite_rhs(self):
        line = example_line(resistance=1e308)
        with self.assertRaises(DomainError):
            line.evaluate(DynamicFeedLineState(10., 2e6),
                          upstream_pressure_pa=2e6, engine_mass_flow_kg_s=0.)

    def test_finite_rhs_does_not_allow_an_infinite_line_jacobian(self):
        line = DynamicFeedLine(replace(example_line().configuration,
            compliance_kg_pa=math.nextafter(0., 1.), inertance_pa_s2_kg=1e308,
            resistance_pa_s2_kg2=0.))
        with self.assertRaises(DomainError):
            line.evaluate(DynamicFeedLineState(0., 2.),
                          upstream_pressure_pa=2., engine_mass_flow_kg_s=0.)

    def test_unrepresentable_feed_drop_fails_with_a_domain_error(self):
        system = example_feed_system()
        state = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
            bottle_gas_temperature_k=300., regulator_opening=.5)
        with self.assertRaises(DomainError):
            system.evaluate(state, 1e200)

    def test_near_full_tank_preserves_positive_ullage_and_rejects_adjacent_overfill(self):
        config = replace(example_feed_system().configuration,
            tank_internal_volume_m3=.0023993932475265857,
            propellant_density_kg_m3=404.72024715263893)
        system = RegulatedFeedSystem(config)
        state = RegulatedFeedState(.9710830281553327, .001, 300., .1, 300., .5)
        with localcontext() as context:
            context.prec = 1200
            exact_ullage = Decimal(config.tank_internal_volume_m3) - (
                Decimal(state.propellant_mass_kg) / Decimal(config.propellant_density_kg_m3)
            )
            expected_pressure = (
                Decimal(state.tank_pressurant_mass_kg)
                * Decimal(config.gas.specific_gas_constant_j_kg_k)
                * Decimal(state.tank_gas_temperature_k) / exact_ullage
            )
            expected_mass_partial = (
                expected_pressure / (Decimal(config.propellant_density_kg_m3) * exact_ullage)
            )
            expected_initial_mass = (
                Decimal(2_000_000) * exact_ullage
                / (Decimal(config.gas.specific_gas_constant_j_kg_k) * Decimal(300))
            )
        self.assertGreater(exact_ullage, 0)
        try:
            result = system.evaluate(state, 0.)
        except DomainError as error:
            raise AssertionError("Representable positive ullage was refused.") from error
        self.assertEqual(result.tank_ullage_volume_m3, float(exact_ullage))
        self.assertEqual(result.tank_pressure_pa, float(expected_pressure))
        partials = dict(result.pressure_partials.injector_pressure_wrt_state)
        self.assertEqual(partials["propellant_mass_kg"], float(expected_mass_partial))
        initial = system.initial_state(
            propellant_mass_kg=state.propellant_mass_kg,
            tank_pressure_pa=2e6, tank_gas_temperature_k=300.0,
            bottle_pressure_pa=20e6, bottle_gas_temperature_k=300.0,
        )
        self.assertEqual(initial.tank_pressurant_mass_kg, float(expected_initial_mass))
        overfilled = replace(state, propellant_mass_kg=math.nextafter(
            state.propellant_mass_kg, math.inf,
        ))
        with self.assertRaises(DomainError):
            system.evaluate(overfilled, 0.)


class FeedConservationOracles(unittest.TestCase):
    def test_energy_balance_near_unity_gamma_against_decimal(self):
        for gamma in (math.nextafter(1., math.inf), 1.000001, 1.4):
            config = replace(example_feed_system().configuration,
                             gas=IdealPressurantGas(300., gamma, "synthetic"))
            system = RegulatedFeedSystem(config)
            state = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
                tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
                bottle_gas_temperature_k=300., regulator_opening=.5)
            result = system.evaluate(state, 0.)
            with localcontext() as ctx:
                ctx.prec = 75
                g, gas_r = Decimal(gamma), Decimal.from_float(300.)
                cv = gas_r / (g - 1)
                cp = g * cv
                # Independent open-system first law, before rearranging Tdot.
                mt, mb, t, q = map(Decimal, (state.tank_pressurant_mass_kg,
                    state.bottle_pressurant_mass_kg, 300., result.regulator_mass_flow_kg_s))
                expected_tank = (q * cp * t - cv * t * q) / (mt * cv)
                expected_bottle = (-q * cp * t + cv * t * q) / (mb * cv)
            with self.subTest(gamma=gamma):
                self.assertTrue(math.isclose(result.derivative.tank_gas_temperature_k_s,
                                             float(expected_tank), rel_tol=1e-14, abs_tol=0))
                self.assertTrue(math.isclose(result.derivative.bottle_gas_temperature_k_s,
                                             float(expected_bottle), rel_tol=1e-14, abs_tol=0))

    def test_line_resistance_dissipates_quadratic_storage_in_both_flow_directions(self):
        line = example_line()
        config = line.configuration
        upstream = 2e6
        for flow in (-.2, 0., .2):
            for pressure in (1.9e6, 2e6, 2.1e6):
                result = line.evaluate(DynamicFeedLineState(flow, pressure),
                                      upstream_pressure_pa=upstream, engine_mass_flow_kg_s=0.)
                storage_rate = (config.inertance_pa_s2_kg * flow
                                * result.derivative.mass_flow_kg_s2
                                + config.compliance_kg_pa * (pressure - upstream)
                                * result.derivative.manifold_pressure_pa_s)
                expected = -config.resistance_pa_s2_kg2 * abs(flow)**3
                with self.subTest(flow=flow, pressure=pressure):
                    self.assertLessEqual(storage_rate, 0.)
                    self.assertAlmostEqual(storage_rate, expected, delta=1e-10)

    def test_closed_adiabatic_ullage_follows_pressure_volume_power_law(self):
        system = RegulatedFeedSystem(replace(example_feed_system().configuration,
                                              valve_time_constant_s=100.))
        initial = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=1e5,
            bottle_gas_temperature_k=300., regulator_opening=0.)
        config = system.configuration
        fields = tuple(initial.__dataclass_fields__)
        start = tuple(getattr(initial, field) for field in fields)
        demand, duration = .05, 60.
        def rhs(values):
            return system.evaluate(RegulatedFeedState(*values), demand).derivative.as_tuple()
        v0 = config.tank_internal_volume_m3 - initial.propellant_mass_kg / config.propellant_density_kg_m3
        final_volume = v0 + demand * duration / config.propellant_density_kg_m3
        exact_pressure = 2e6 * (v0 / final_volume) ** config.gas.gamma
        errors = []
        for steps in (20, 40, 80):
            values = start
            for _ in range(steps):
                values = step_rk4(rhs, values, duration / steps)
            final = system.evaluate(RegulatedFeedState(*values), demand)
            errors.append(abs(final.tank_pressure_pa - exact_pressure))
            self.assertEqual(values[1], start[1])
            self.assertEqual(values[3], start[3])
            self.assertEqual(values[4], start[4])
        self.assertLess(errors[-1] / exact_pressure, 2e-11)
        self.assertGreater(errors[0], errors[1] * 10)
        self.assertGreater(errors[1], errors[2] * 10)

    def test_bottle_blowdown_preserves_adiabatic_invariant_and_total_pressurant(self):
        system = example_feed_system()
        initial = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
            bottle_gas_temperature_k=300., regulator_opening=1.)
        start = tuple(getattr(initial, f) for f in initial.__dataclass_fields__)
        def rhs(values):
            return system.evaluate(RegulatedFeedState(*values), 0.).derivative.as_tuple()
        values = start
        for _ in range(400):
            values = step_rk4(rhs, values, .005)
        expected_temperature = start[4] * (values[3] / start[3]) ** (
            system.configuration.gas.gamma - 1)
        self.assertLess(values[3], start[3] * .999)
        self.assertAlmostEqual(values[1] + values[3], start[1] + start[3], delta=1e-14)
        self.assertAlmostEqual(values[4], expected_temperature, delta=2e-10)


if __name__ == "__main__":
    unittest.main()
