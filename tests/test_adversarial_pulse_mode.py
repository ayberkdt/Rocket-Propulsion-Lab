"""Try to break pulse-mode planning with random thrusters and nasty requests.

Oracles written in this file:

* Equal-pulse rocket equation: n pulses of impulse I and propellant m from
  mass m0 give Δv = (I/m) ln(m0 / (m0 - n m)) (one exhaust velocity c = I/m
  for the whole train; NASA Glenn, Ideal Rocket Equation,
  https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/).
* Achievable sets by exhaustive enumeration: for every pulse count n, the
  feasible widths are [t_min, t_max(n)] with n m(t_max) <= budget (found here
  by bisection on m), and Δv_n(t) is continuous and increasing in t, so the
  achievable Δv of n pulses is the interval [Δv_n(t_min), Δv_n(t_max(n))].
  Impulse requests use the intervals [n I_min, n I(t_longest)].
* Shape preservation of the tables (monotone, no overshoot) is the PCHIP
  property: F. N. Fritsch and R. E. Carlson, SIAM J. Numer. Anal. 17(2),
  1980 (bibliographic citation).
* Minimum impulse bit as the smallest commandable pulse: J. M. Parker and
  M. J. Wilson, The Minimum Impulse Thruster, NASA NTRS 20070032005,
  https://ntrs.nasa.gov/citations/20070032005
"""

import copy
import json
import pickle
import random
import unittest
from dataclasses import replace
from math import log1p, nextafter

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import PulseModeThruster, PulsePlan

NEAREST, SAFE = "nearest", "no_overshoot"


def random_thruster(generator: random.Random) -> PulseModeThruster:
    nodes = generator.randint(2, 6)
    on_times = [0.005 * (1 + generator.random())]
    for _ in range(nodes - 1):
        on_times.append(on_times[-1] * generator.uniform(1.2, 3.0))
    thrust = generator.uniform(0.5, 25.0)
    isp = generator.uniform(150.0, 320.0)
    flow = thrust / (9.80665 * isp)
    impulses, masses = [], []
    for t in on_times:
        efficiency = generator.uniform(0.2, 1.0)
        impulses.append(max(thrust * t * efficiency, (impulses[-1] if impulses else 0) * 1.01
                            + 1e-6))
        mass = flow * t * generator.uniform(0.9, 1.6)
        # Equal successive propellant values are allowed ("non-decreasing").
        if masses and generator.random() < 0.2:
            mass = masses[-1]
        masses.append(max(mass, masses[-1] if masses else mass))
    longest = None if generator.random() < 0.5 else on_times[-1] * generator.uniform(1.0, 4.0)
    return PulseModeThruster(
        on_times_s=tuple(on_times), impulse_bits_n_s=tuple(impulses),
        propellant_per_pulse_kg=tuple(masses), steady_thrust_n=thrust,
        steady_mass_flow_kg_s=flow, minimum_off_time_s=0.01, source_id="random",
        source_sha256="b" * 64, maximum_on_time_s=longest,
    )


def delta_v(thruster, mass, count, on_time):
    bit = thruster.pulse(on_time)
    # -log1p(-u/m0) = ln(m0/(m0-u)) without cancellation (DLMF 4.6.1).
    return bit.impulse_n_s / bit.propellant_kg * -log1p(-count * bit.propellant_kg / mass)


def widest(thruster, count, budget):
    short, long_ = thruster.minimum_on_time_s, thruster.longest_on_time_s
    if count * thruster.pulse(long_).propellant_kg <= budget:
        return long_
    low, high = short, long_
    for _ in range(200):
        middle = 0.5 * (low + high)
        if middle in (low, high):  # converged to adjacent doubles
            break
        if count * thruster.pulse(middle).propellant_kg <= budget:
            low = middle
        else:
            high = middle
    return low


def delta_v_intervals(thruster, mass, budget):
    most = int(budget // thruster.pulse(thruster.minimum_on_time_s).propellant_kg)
    intervals = [(0, 0.0, 0.0)]
    for count in range(1, most + 1):
        if count * thruster.pulse(thruster.minimum_on_time_s).propellant_kg > budget:
            break
        intervals.append((count, delta_v(thruster, mass, count, thruster.minimum_on_time_s),
                          delta_v(thruster, mass, count, widest(thruster, count, budget))))
    return intervals


# Documented: within 1e-9 relative counts as meeting the request.  The oracle
# demands exactness only well inside that band, so ties at the band edge are
# judged by the nearest-distance check instead.
TOLERANCE = 0.5e-9


def oracle(intervals, request):
    slack = TOLERANCE * request
    exact = [n for n, low, high in intervals if low - slack <= request <= high + slack]
    nearest = min(0.0 if low <= request <= high else min(abs(request - low), abs(request - high))
                  for _, low, high in intervals)
    best_below = max(min(high, request) for _, low, high in intervals if low <= request)
    return (min(exact) if exact else None), nearest, best_below


class DeltaVPlanAttackTests(unittest.TestCase):
    def test_random_thrusters_and_hostile_requests_against_enumeration(self) -> None:
        generator = random.Random(8675309)
        checked = 0
        for case in range(80):
            thruster = random_thruster(generator)
            mass = generator.uniform(5.0, 500.0)
            m_short = thruster.pulse(thruster.minimum_on_time_s).propellant_kg
            budget = min(m_short * generator.uniform(1.5, 400.0), 0.9 * mass)
            intervals = delta_v_intervals(thruster, mass, budget)
            top = max(high for _, _, high in intervals)
            requests = [generator.uniform(0.0, top) for _ in range(8)]
            for count, low, high in generator.sample(intervals[1:], min(4, len(intervals) - 1)):
                requests += [low, high, nextafter(low, 0.0), nextafter(high, 2 * high),
                             0.5 * (low + high)]
            requests += [nextafter(top, 0.0), top]
            for request in requests:
                exact_count, nearest, best_below = oracle(intervals, request)
                scale = max(request, 1e-30)
                for policy in (NEAREST, SAFE):
                    plan = thruster.plan_delta_v(request, mass, propellant_budget_kg=budget,
                                                 quantization=policy)
                    checked += 1
                    label = (case, request, policy, plan.pulse_count, exact_count)
                    self.assertLessEqual(plan.propellant_kg, budget * (1 + 1e-12), label)
                    if exact_count is not None:
                        self.assertFalse(plan.quantized, label)
                        self.assertEqual(plan.pulse_count, exact_count, label)
                        self.assertLessEqual(abs(plan.residual), 1e-9 * scale, label)
                    elif policy == NEAREST:
                        self.assertLessEqual(abs(plan.residual), nearest + 1e-9 * scale, label)
                    else:
                        self.assertGreaterEqual(plan.residual, -1e-12 * scale, label)
                        self.assertGreaterEqual(plan.delivered, best_below - 1e-9 * scale, label)
                    if plan.pulse is not None:
                        recomputed = delta_v(thruster, mass, plan.pulse_count,
                                             plan.pulse.on_time_s)
                        self.assertLessEqual(abs(recomputed - plan.delivered), 1e-12 * scale + 1e-15)
            # Above everything the budget can deliver: refused, never "quantized".
            with self.assertRaises(DomainError):
                thruster.plan_delta_v(top * 1.001 + 1e-9, mass, propellant_budget_kg=budget)
        self.assertGreater(checked, 2000)

    def test_budget_limited_shorter_pulses_are_found(self) -> None:
        # Longest pulse 2.4..3.0 kg: three longest fit in 9.76 kg, four shorter ones too.
        thruster = PulseModeThruster(
            on_times_s=(0.02, 0.05), impulse_bits_n_s=(6000.0, 9000.0),
            propellant_per_pulse_kg=(2.4, 3.0), steady_thrust_n=1e5,
            steady_mass_flow_kg_s=40.0, minimum_off_time_s=0.0, source_id="big",
            source_sha256="b" * 64,
        )
        intervals = delta_v_intervals(thruster, 10.0, 9.76)
        # Three longest pulses reach at most 6908 m/s; four pulses only fit when
        # shortened and then deliver 8047..9300 m/s.  8500 m/s is therefore met
        # exactly by four budget-limited pulses, and 7000 m/s lies in the gap.
        self.assertEqual(oracle(intervals, 8500.0)[0], 4)
        plan = thruster.plan_delta_v(8500.0, 10.0, propellant_budget_kg=9.76)
        self.assertEqual((plan.pulse_count, plan.quantized), (4, False))
        self.assertLessEqual(plan.propellant_kg, 9.76 * (1 + 1e-12))
        for request in (4605.0, 6000.0, 7000.0, 8047.0, 9000.0):
            plan = thruster.plan_delta_v(request, 10.0, propellant_budget_kg=9.76)
            exact_count, nearest, _ = oracle(intervals, request)
            if exact_count is not None:
                self.assertEqual(plan.pulse_count, exact_count, request)
            else:
                self.assertLessEqual(abs(plan.residual), nearest + 1e-9 * request, request)
        self.assertEqual(thruster.plan_delta_v(7000.0, 10.0, propellant_budget_kg=9.76)
                         .pulse_count, 3)

    def test_repeated_requests_preserve_budget_mass_and_model_independence(self) -> None:
        first = PulseModeThruster(
            on_times_s=(0.02, 0.05), impulse_bits_n_s=(6000.0, 9000.0),
            propellant_per_pulse_kg=(2.4, 3.0), steady_thrust_n=1e5,
            steady_mass_flow_kg_s=40.0, minimum_off_time_s=0.0, source_id="big",
            source_sha256="b" * 64)
        second = replace(first, propellant_per_pulse_kg=(2.0, 3.5))
        hashes = (first.model_hash, second.model_hash)
        for _repeat in range(2):
            for thruster in (first, second, first):
                for budget in (9.76, 5.1, 7.3, 9.76):
                    for mass in (10.0, 15.0):
                        intervals = delta_v_intervals(thruster, mass, budget)
                        for fraction in (0.25, 0.6, 0.95):
                            request = fraction * max(high for _, _, high in intervals)
                            exact, nearest, below = oracle(intervals, request)
                            for policy in (NEAREST, SAFE):
                                plan = thruster.plan_delta_v(
                                    request, mass, propellant_budget_kg=budget, quantization=policy)
                                self.assertLessEqual(plan.propellant_kg, budget * (1 + 1e-12))
                                if exact is not None:
                                    self.assertEqual(plan.pulse_count, exact)
                                    self.assertLessEqual(abs(plan.residual), 1e-9 * request)
                                elif policy == NEAREST:
                                    self.assertLessEqual(abs(plan.residual), nearest + 1e-9 * request)
                                else:
                                    self.assertGreaterEqual(plan.residual, -1e-12 * request)
                                    self.assertGreaterEqual(plan.delivered, below - 1e-9 * request)
        self.assertEqual((first.model_hash, second.model_hash), hashes)
        # Fresh construction and normal dataclass copying/serialization must
        # preserve plans; cache execution state is not exchange data.
        reopened = replace(first)
        for model in (reopened, copy.deepcopy(first), pickle.loads(pickle.dumps(first))):
            self.assertEqual(first.plan_delta_v(8500, 10, propellant_budget_kg=9.76),
                             model.plan_delta_v(8500, 10, propellant_budget_kg=9.76))

    def test_budget_must_leave_dry_mass(self) -> None:
        thruster = random_thruster(random.Random(1))
        for budget in (0.0, -1.0, 10.0, 11.0):
            with self.assertRaises(DomainError):
                thruster.plan_delta_v(0.1, 10.0, propellant_budget_kg=budget)


class ImpulsePlanAttackTests(unittest.TestCase):
    def test_random_impulse_requests_against_enumeration(self) -> None:
        generator = random.Random(31337)
        for _ in range(150):
            thruster = random_thruster(generator)
            low_bit = thruster.minimum_impulse_bit_n_s
            high_bit = thruster.pulse(thruster.longest_on_time_s).impulse_n_s
            intervals = [(0, 0.0, 0.0)] + [(n, n * low_bit, n * high_bit) for n in range(1, 400)]
            requests = [generator.uniform(0.0, 50 * high_bit) for _ in range(10)]
            for n in (1, 2, 3, generator.randint(4, 40)):
                requests += [n * low_bit, n * high_bit, nextafter(n * low_bit, 0.0),
                             nextafter(n * high_bit, 1e9)]
            for request in requests:
                exact_count, nearest, best_below = oracle(intervals, request)
                scale = max(request, 1e-30)
                plan = thruster.plan_impulse(request)
                safe = thruster.plan_impulse(request, quantization=SAFE)
                if exact_count is not None:
                    self.assertEqual(plan.pulse_count, exact_count, request)
                    self.assertLessEqual(abs(plan.residual), 1e-9 * scale)
                else:
                    self.assertLessEqual(abs(plan.residual), nearest + 1e-12 * scale)
                    self.assertGreaterEqual(safe.residual, -1e-12 * scale)
                    self.assertGreaterEqual(safe.delivered, best_below - 1e-12 * scale)

    def test_pulse_limit_applies_to_the_returned_plan(self) -> None:
        thruster = PulseModeThruster(
            on_times_s=(0.01, 0.012), impulse_bits_n_s=(1.0, 1.2),
            propellant_per_pulse_kg=(1e-3, 1.1e-3), steady_thrust_n=100.0,
            steady_mass_flow_kg_s=0.05, minimum_off_time_s=0.0, source_id="lim",
            source_sha256="b" * 64,
        )
        # 1.5 N s lies between one longest pulse (1.2) and two minimum bits (2.0).
        for policy in (NEAREST, SAFE):
            plan = thruster.plan_impulse(1.5, quantization=policy, maximum_pulses=1)
            self.assertEqual(plan.pulse_count, 1)
        with self.assertRaises(DomainError):
            thruster.plan_impulse(1.9, maximum_pulses=1)  # nearest is 2 pulses
        with self.assertRaises(DomainError):
            thruster.plan_impulse(10.0, maximum_pulses=3)


class ClaimTests(unittest.TestCase):
    """Properties the module docstring claims, checked by dense sampling."""

    def test_tables_are_shape_preserving_and_bounds_hold(self) -> None:
        generator = random.Random(2718)
        for _ in range(80):
            thruster = random_thruster(generator)
            bound = thruster.exhaust_velocity_bound_m_s
            previous_impulse = previous_mass = 0.0
            short, long_ = thruster.minimum_on_time_s, thruster.longest_on_time_s
            for k in range(801):
                t = short + (long_ - short) * k / 800
                bit = thruster.pulse(t)
                self.assertGreater(bit.impulse_n_s, previous_impulse - 1e-15)
                self.assertGreaterEqual(bit.propellant_kg, previous_mass - 1e-18)
                self.assertLessEqual(bit.exhaust_velocity_m_s, bound * (1 + 1e-12))
                previous_impulse, previous_mass = bit.impulse_n_s, bit.propellant_kg
            self.assertLessEqual(previous_impulse, max(thruster.impulse_bits_n_s[-1],
                                                       thruster.pulse(long_).impulse_n_s) + 1e-12)

    def test_train_delta_v_increases_with_width_and_count(self) -> None:
        generator = random.Random(1618)
        for _ in range(60):
            thruster = random_thruster(generator)
            mass = 1e3
            short, long_ = thruster.minimum_on_time_s, thruster.longest_on_time_s
            for count in (1, 2, 7, 30):
                previous = -1.0
                for k in range(301):
                    t = short + (long_ - short) * k / 300
                    value = delta_v(thruster, mass, count, t)
                    self.assertGreaterEqual(value, previous - 1e-12 * value)
                    previous = value
                self.assertLess(delta_v(thruster, mass, count, short),
                                delta_v(thruster, mass, count + 1, short))


class PlanIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.thruster = random_thruster(random.Random(5))
        self.plan = self.thruster.plan_impulse(10 * self.thruster.minimum_impulse_bit_n_s)

    def test_artifacts_refuse_forged_or_foreign_plans(self) -> None:
        other = random_thruster(random.Random(6))
        foreign = other.plan_impulse(10 * other.minimum_impulse_bit_n_s)
        period = 10.0
        forged = (
            foreign,
            replace(self.plan, pulse_count=2.5),
            replace(self.plan, pulse_count=0),
            replace(self.plan, pulse_count=True),
            replace(self.plan, pulse=replace(self.plan.pulse, impulse_n_s=1e9)),
            replace(self.plan, pulse=replace(self.plan.pulse, on_time_s=0.0)),
            replace(self.plan, pulse=None),
            "plan",
        )
        for plan in forged:
            with self.assertRaises((DomainError, InputError), msg=repr(plan)[:80]):
                self.thruster.pulse_train_artifact(plan, period)
        artifact = self.thruster.pulse_train_artifact(self.plan, period)
        self.assertAlmostEqual(artifact.delivered_total_impulse_n_s,
                               self.plan.delivered_impulse_n_s, places=12)
        self.assertIsInstance(self.plan, PulsePlan)

    def test_every_tampered_field_and_hostile_document_is_refused(self) -> None:
        document = self.thruster.manifest()
        for key in document:
            if key == "model_hash":
                continue
            for value in (None, "1", 1e300, [], {}, True, -1.0):
                if value is None and document[key] is None:
                    continue  # not a change
                tampered = dict(document, **{key: value})
                with self.assertRaises((DomainError, InputError), msg=key):
                    PulseModeThruster.from_dict(tampered)
            missing = {k: v for k, v in document.items() if k != key}
            with self.assertRaises((DomainError, InputError), msg=key):
                PulseModeThruster.from_dict(missing)
        text = self.thruster.to_json()
        hostile = (
            "9" * 5000, "[" * 100_000 + "]" * 100_000, text.replace("0.0", "NaN", 1),
            text.replace("0.0", "Infinity", 1), "﻿" + text, text[: len(text) // 2], "",
        )
        for value in hostile:
            with self.assertRaises((DomainError, InputError)):
                PulseModeThruster.from_json(value)
        # Duplicate keys: JSON keeps the last value, so the hash must catch it.
        duplicated = text.rstrip().rstrip("}") + ', "steady_thrust_n": 1234.5}'
        with self.assertRaises(InputError):
            PulseModeThruster.from_json(duplicated)
        self.assertEqual(PulseModeThruster.from_json(text), self.thruster)
        json.loads(text)  # the canonical document itself is plain JSON


if __name__ == "__main__":
    unittest.main()
