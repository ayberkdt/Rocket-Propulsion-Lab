"""Direct inventory helpers: independent ledgers and hostile reserve boundaries."""

import copy
import itertools
import random
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import inf, nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import inventory
from rocket_propulsion.propulsion.burns.models import PropellantRole
from tests.test_burn_constant_numerics import mono


def ledger_fixture(rates=(1.,), loaded=10., reserve=0.):
    definition, point = mono()
    tanks = (replace(definition.tanks[0], loaded_mass_kg=loaded, reserve_mass_kg=reserve),)
    streams = tuple(replace(point.streams[0], stream_id=f"path-{i}", mass_flow_kg_s=rate)
                    for i, rate in enumerate(rates))
    # Only ledger fields are changed; declared engine performance is not an
    # oracle for these direct stream/tank helpers. Keep constructors independent
    # from the arithmetic extremes under test.
    point = copy.copy(point)
    object.__setattr__(point, "streams", streams)
    return tanks, point


def decimal_remaining(tank, rates, time, precision):
    with localcontext() as context:
        context.prec = precision
        return Decimal.from_float(tank.loaded_mass_kg) - sum(
            (Decimal.from_float(rate) for rate in rates), Decimal(0)) * Decimal.from_float(time)


class LegacyInventoryAttacks(unittest.TestCase):
    def test_missing_mapping_is_rejected_by_all_inventory_queries(self):
        _, point = ledger_fixture()
        for query, args in ((inventory.validate_tank_coverage, ()),
                            (inventory.tank_duration_limit, ()),
                            (inventory.tank_masses_at, (0.,))):
            with self.subTest(query=query.__name__), self.assertRaises(DomainError):
                query((), point, *args)

    def test_time_requires_a_finite_nonnegative_real(self):
        tanks, point = ledger_fixture()
        for time in (float('nan'), inf, -inf, -1., True, False, None, '1', []):
            with self.subTest(time=time), self.assertRaises(DomainError):
                inventory.tank_masses_at(tanks, point, time)
        self.assertEqual(inventory.tank_masses_at(tanks, point, -0.), (('tank', 10.),))
        self.assertEqual(inventory.tank_masses_at(tanks, point, 1), (('tank', 9.),))

    def test_reserve_crossing_cannot_hide_in_absolute_or_relative_tolerance(self):
        for loaded, reserve, rate, time in ((1., 1., 1., 1e-14),
                                            (1e16, 1e16, 1., 1.),
                                            (1e-20, 0., 1e-20, nextafter(1., inf))):
            tanks, point = ledger_fixture((rate,), loaded, reserve)
            with self.subTest(loaded=loaded), self.assertRaises(DomainError):
                inventory.tank_masses_at(tanks, point, time)

    def test_shared_tank_small_paths_survive_before_reserve_decision(self):
        tanks, point = ledger_fixture((1e16, 1.), 1e16)
        self.assertEqual(inventory.flow_by_tank(point), {'tank': 1e16})
        with self.assertRaises(DomainError):
            inventory.tank_masses_at(tanks, point, 1.)
        horizon, binding = inventory.tank_duration_limit(tanks, point)
        self.assertEqual(horizon, nextafter(1., 0.))
        self.assertEqual(binding, 'tank')
        self.assertGreater(dict(inventory.tank_masses_at(tanks, point, horizon))['tank'], 0.)

    def test_horizon_is_the_greatest_safe_float_not_nearest_physical_time(self):
        tanks, point = ledger_fixture((10.,), 1.)
        horizon, binding = inventory.tank_duration_limit(tanks, point)
        self.assertEqual((horizon, binding), (nextafter(.1, 0.), 'tank'))
        for precision in (100, 200):
            remaining = decimal_remaining(tanks[0], (10.,), horizon, precision)
            self.assertGreaterEqual(remaining, 0)
            self.assertEqual(inventory.tank_masses_at(tanks, point, horizon),
                             (('tank', float(remaining)),))
            self.assertLess(decimal_remaining(tanks[0], (10.,), nextafter(horizon, inf), precision), 0)
        with self.assertRaises(DomainError):
            inventory.tank_masses_at(tanks, point, nextafter(horizon, inf))

    def test_binding_tank_is_chosen_before_rounding_and_ties_keep_input_order(self):
        tanks, point = ledger_fixture((1.,), 1.)
        extra = replace(tanks[0], tank_id='other')
        # tank's rate is 1 + 2**-54; displayed rates and nearest cutoff times
        # coincide, but its exact reserve boundary is strictly earlier.
        object.__setattr__(point, 'streams', (
            point.streams[0], replace(point.streams[0], stream_id='tiny', mass_flow_kg_s=2.**-54),
            replace(point.streams[0], stream_id='second', tank_id='other')))
        self.assertEqual(inventory.tank_duration_limit((extra, *tanks), point),
                         (nextafter(1., 0.), 'tank'))
        object.__setattr__(point, 'streams', (point.streams[0], point.streams[2]))
        self.assertEqual(inventory.tank_duration_limit((extra, *tanks), point), (1., 'other'))

    def test_stream_order_does_not_change_grouped_export_or_final_ledger(self):
        rates = (1e16, 1., 1.)
        for order in itertools.permutations(rates):
            tanks, point = ledger_fixture(order, 1e16)
            self.assertEqual(inventory.flow_by_tank(point)['tank'], 1e16+2.)
            expected = decimal_remaining(tanks[0], order, .5, 100)
            self.assertEqual(inventory.tank_masses_at(tanks, point, .5), (('tank', float(expected)),))

    def test_tank_duplicates_and_incompatible_identity_cannot_be_overwritten(self):
        tanks, point = ledger_fixture()
        variants = ((tanks[0], replace(tanks[0], loaded_mass_kg=20.)),
                    (replace(tanks[0], propellant_key='different'),),
                    (replace(tanks[0], role=PropellantRole.FUEL),))
        for bad in variants:
            for query, args in ((inventory.validate_tank_coverage, ()),
                                (inventory.tank_duration_limit, ()),
                                (inventory.tank_masses_at, (0.,))):
                with self.subTest(bad=bad, query=query.__name__), self.assertRaises(DomainError):
                    query(bad, point, *args)

    def test_malformed_collection_point_and_tracking_flag_fail_typed(self):
        tanks, point = ledger_fixture()
        for bad in (None, {}, 'tank', (None,)):
            with self.subTest(tanks=bad), self.assertRaises(DomainError):
                inventory.validate_tank_coverage(bad, point)
        for bad in (None, {}, object()):
            with self.subTest(point=bad), self.assertRaises(DomainError):
                inventory.flow_by_tank(bad)
        for flag in ('false', 0, 1, None):
            changed = copy.copy(point)
            stream = copy.copy(point.streams[0])
            object.__setattr__(stream, 'tank_depleting', flag)
            object.__setattr__(changed, 'streams', (stream,))
            with self.subTest(flag=flag), self.assertRaises(DomainError):
                inventory.flow_by_tank(changed)

    def test_unknown_roles_and_nonnumeric_or_boolean_inventory_fail_typed(self):
        tanks, point = ledger_fixture()
        for field, value in (('role', 'unknown'), ('loaded_mass_kg', True),
                             ('reserve_mass_kg', False), ('loaded_mass_kg', '10')):
            bad = copy.copy(tanks[0])
            object.__setattr__(bad, field, value)
            with self.subTest(field=field), self.assertRaises(DomainError):
                inventory.validate_tank_coverage((bad,), point)
        # Existing records may carry recognized enum spellings; compare meaning
        # without mutating caller records or silently accepting unknown roles.
        known = copy.copy(tanks[0])
        object.__setattr__(known, 'role', 'monopropellant')
        inventory.validate_tank_coverage([known], point)
        self.assertIsInstance(known.role, str)
        self.assertNotIsInstance(known.role, PropellantRole)

    def test_zero_no_drain_and_untracked_paths_preserve_inventory(self):
        tanks, point = ledger_fixture((0.,))
        self.assertEqual(inventory.tank_duration_limit(tanks, point), (inf, None))
        self.assertEqual(inventory.tank_masses_at(tanks, point, float_info.max), (('tank', 10.),))
        object.__setattr__(point, 'streams', (replace(point.streams[0], tank_id='',
                              mass_flow_kg_s=1., tank_depleting=False),))
        self.assertEqual(inventory.flow_by_tank(point), {})
        self.assertEqual(inventory.tank_duration_limit((), point), (inf, None))
        self.assertEqual(inventory.tank_masses_at((), point, 10.), ())

    def test_extreme_horizons_retain_binding_and_all_finite_safe_queries(self):
        tiny = nextafter(0., 1.)
        tanks, point = ledger_fixture((tiny,), 1.)
        self.assertEqual(inventory.tank_duration_limit(tanks, point), (inf, 'tank'))
        expected = decimal_remaining(tanks[0], (tiny,), float_info.max, 1500)
        self.assertEqual(expected, decimal_remaining(tanks[0], (tiny,), float_info.max, 1800))
        self.assertEqual(inventory.tank_masses_at(tanks, point, float_info.max), (('tank', float(expected)),))
        tanks, point = ledger_fixture((float_info.max,), tiny)
        self.assertEqual(inventory.tank_duration_limit(tanks, point), (0., 'tank'))
        self.assertEqual(inventory.tank_masses_at(tanks, point, 0.), (('tank', tiny),))
        with self.assertRaises(DomainError):
            inventory.tank_masses_at(tanks, point, tiny)
        tanks, point = ledger_fixture((1.,), float_info.max)
        self.assertEqual(inventory.tank_duration_limit(tanks, point), (float_info.max, 'tank'))
        self.assertEqual(inventory.tank_masses_at(tanks, point, float_info.max), (('tank', 0.),))
        tanks, point = ledger_fixture((1.,), 1., 1.)
        self.assertEqual(inventory.tank_duration_limit(tanks, point), (0., 'tank'))
        self.assertEqual(inventory.tank_masses_at(tanks, point, 0.), (('tank', 1.),))

    def test_real_output_range_failures_are_explicit_but_coverage_needs_no_rate_export(self):
        tanks, point = ledger_fixture((1e308, 1e308), 1.)
        inventory.validate_tank_coverage(tanks, point)
        self.assertGreater(inventory.tank_duration_limit(tanks, point)[0], 0.)
        self.assertEqual(inventory.tank_masses_at(tanks, point, 0.), (('tank', 1.),))
        with self.assertRaises(DomainError) as caught:
            inventory.flow_by_tank(point)
        self.assertIn('tank', str(caught.exception).lower())
        tiny = nextafter(0., 1.)
        tanks, point = ledger_fixture((tiny,), tiny)
        with self.assertRaises(DomainError):
            inventory.tank_masses_at(tanks, point, .75)
        self.assertEqual(inventory.tank_masses_at(tanks, point, 1.), (('tank', 0.),))

    def test_seeded_ledgers_and_adjacent_horizons_against_paired_decimal(self):
        rng = random.Random(65429)
        for _ in range(48):
            scale = 2.**rng.randint(-700, 700)
            rates = tuple(scale*rng.randint(1, 31)/16 for _ in range(rng.randint(1, 5)))
            loaded = scale*rng.randint(1, 31)
            reserve = loaded/4
            tanks, point = ledger_fixture(rates, loaded, reserve)
            horizon, binding = inventory.tank_duration_limit(tanks, point)
            a = decimal_remaining(tanks[0], rates, horizon, 800)
            b = decimal_remaining(tanks[0], rates, horizon, 1200)
            self.assertEqual(float(a), float(b))
            self.assertGreaterEqual(a, Decimal.from_float(reserve))
            self.assertLess(decimal_remaining(tanks[0], rates, nextafter(horizon, inf), 1200),
                            Decimal.from_float(reserve))
            self.assertEqual(binding, 'tank')
            self.assertEqual(inventory.tank_masses_at(tanks, point, horizon), (('tank', float(b)),))
            with self.assertRaises(DomainError):
                inventory.tank_masses_at(tanks, point, nextafter(horizon, inf))

    def test_malformed_stream_records_identifiers_and_rates_fail_typed(self):
        _, point = ledger_fixture()
        for bad in (None, {}, 'streams', (None,)):
            changed = copy.copy(point)
            object.__setattr__(changed, 'streams', bad)
            with self.subTest(streams=bad), self.assertRaises(DomainError):
                inventory.flow_by_tank(changed)
        for field, value in (('stream_id', None), ('tank_id', ''), ('propellant_key', ''),
                             ('role', 'unknown'), ('mass_flow_kg_s', True),
                             ('mass_flow_kg_s', float('nan')), ('mass_flow_kg_s', -1.)):
            stream = copy.copy(point.streams[0])
            object.__setattr__(stream, field, value)
            changed = copy.copy(point)
            object.__setattr__(changed, 'streams', (stream,))
            with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                inventory.flow_by_tank(changed)

    def test_ordinary_multiple_tanks_keep_order_and_leave_idle_tank_untouched(self):
        tanks, point = ledger_fixture((2., 3.), 20., 5.)
        # Exercise list inputs as well as immutable tuples.
        idle = replace(tanks[0], tank_id='idle', loaded_mass_kg=4., reserve_mass_kg=1.)
        self.assertEqual(inventory.tank_duration_limit([idle, *tanks], point), (3., 'tank'))
        self.assertEqual(inventory.tank_masses_at([idle, *tanks], point, 3.),
                         (('idle', 4.), ('tank', 5.)))
        self.assertEqual(tanks[0].loaded_mass_kg, 20.)
        self.assertEqual(tuple(stream.mass_flow_kg_s for stream in point.streams), (2., 3.))
