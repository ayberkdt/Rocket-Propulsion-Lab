"""Whole-field radial oracles, conservation audits and hostile marching data."""
import itertools
import unittest
from collections import Counter
from dataclasses import asdict, replace
from decimal import ROUND_DOWN, Inexact, localcontext
from fractions import Fraction
from math import ldexp, nextafter, sqrt

import mpmath as mp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import (
    MocFieldInputs,
    MocFieldResult,
    MocFlowPoint,
    MocWallSegment,
    solve_moc_axis_departure,
    solve_moc_axis_point,
    solve_moc_field,
)
from rocket_propulsion.geometry import moc_field as model
from tests.test_moc_unit_process import exact_mp


def radial_fixture(n=4, *, axisymmetric=True, compressing=False, gamma=1.4, mach=2., precision=110):
    """rho*V*R**d and analytic Mach-line angles; no production field equations.

    Analytic data are supplied only at the inlet and old axis history. Every
    later axis, interior and wall state must come from the computed field.
    """
    with mp.workdps(precision):
        g, m0 = map(exact_mp, (gamma, mach))
        d = 2 if axisymmetric else 1
        h = mp.mpf('.1') / n
        sign = -1 if compressing else 1

        def area(m):
            return (2 / (g + 1) * (1 + (g - 1) * m * m / 2)) ** ((g + 1) / (2 * (g - 1))) / m

        def nu(m):
            return mp.sqrt((g + 1) / (g - 1)) * mp.atan(mp.sqrt((g - 1) / (g + 1) * (m * m - 1))) - mp.atan(mp.sqrt(m * m - 1))

        def polar(m, beta):
            radius = 2 * (area(m) / area(m0)) ** (mp.mpf(1) / d)
            x = radius * mp.cos(beta)
            return MocFlowPoint(float(x if sign == 1 else 4 - x), float(radius * mp.sin(beta)),
                                float(m), float(sign * beta * 180 / mp.pi))

        history = tuple(polar(m0 - sign * k * h, 0) for k in (3, 2))
        incoming = polar(m0 - sign * h, abs(nu(m0) - nu(m0 - sign * h)) / d)
        axis = solve_moc_axis_point(incoming=incoming, axis_history=history, gamma=gamma, axisymmetric=axisymmetric)

        def exact(x, y):
            with mp.workdps(precision):
                x, y = map(exact_mp, (x, y))
                radial_x = x if sign == 1 else 4 - x
                radius = mp.sqrt(radial_x * radial_x + y * y)
                m = mp.findroot(lambda z: mp.log(area(z) / area(m0)) - d * mp.log(radius / 2), (m0 - mp.mpf('.1'), m0 + mp.mpf('.1')))
                return float(m), float(sign * mp.atan2(y, radial_x) * 180 / mp.pi)

        x = axis.point.x_m
        height = (x if sign == 1 else 4 - x) * .05
        inlet = (axis.point,) + tuple(MocFlowPoint(x, height * j / n, *exact(x, height * j / n)) for j in range(1, n + 1))
        wall = MocWallSegment(x, height, x + .8, height + sign * .04)
    return {'axis_result': axis, 'inlet_points': inlet, 'wall_segment': wall, 'marching_steps': n}, exact


def uniform_fixture(*, mach=2., axisymmetric=True):
    with mp.workdps(110):
        k = float(mp.sqrt(exact_mp(mach) ** 2 - 1))
    axis = solve_moc_axis_point(incoming=MocFlowPoint(-k, 1., mach, 0.),
        axis_history=(MocFlowPoint(-3 * k, 0., mach, 0.), MocFlowPoint(-2 * k, 0., mach, 0.)), axisymmetric=axisymmetric)
    x = axis.point.x_m
    return {'axis_result': axis, 'inlet_points': tuple(MocFlowPoint(x, r, mach, 0.) for r in (0., .5, 1.)),
            'wall_segment': MocWallSegment(x, 1., x + 8 * k, 1.), 'marching_steps': 3}


class MocFieldChecks:
    def valid(self, **inputs):
        try:
            return solve_moc_field(**inputs)
        except (DomainError, ArithmeticError) as error:
            self.fail(f'Valid whole-field witness refused: {error}')

    def convergence(self, mode, compress):
        errors = []
        cell_errors = []
        for n in (4, 8, 16):
            inputs, exact = radial_fixture(n, axisymmetric=mode, compressing=compress)
            result = self.valid(**inputs)
            mach_error = angle_error = 0.
            for row in (*result.fronts, *result.half_fronts):
                for p in row:
                    m, theta = exact(p.x_m, p.radius_m)
                    mach_error = max(mach_error, abs(p.mach / m - 1))
                    angle_error = max(angle_error, abs(p.flow_angle_deg - theta))
            errors.append(max(mach_error, angle_error))
            cell_errors.append(result.max_cell_relative_mass_imbalance)
            self.assertGreater(result.max_relative_mass_drift, 0.)
            self.assertLess(result.wall_mass_leakage_relative, 1e-13)
        self.assertTrue(all(3 < a / b < 12 for a, b in itertools.pairwise(errors)), errors)
        self.assertTrue(all(3 < a / b < 12 for a, b in itertools.pairwise(cell_errors)), cell_errors)
        self.assertLess(mach_error, 1e-5)
        self.assertLess(angle_error, 1e-4)
        self.assertLess(result.max_relative_mass_drift, 1e-6)
        fine = self.valid(**radial_fixture(32, axisymmetric=mode, compressing=compress)[0])
        self.assertTrue(3 < cell_errors[-1] / fine.max_cell_relative_mass_imbalance < 12)
        self.assertLess(fine.max_cell_relative_mass_imbalance, 2e-6)
        self.assertLess(fine.max_relative_mass_drift, 1e-6)





class MocFieldAttacks(MocFieldChecks, unittest.TestCase):
    def test_uniform_flow_has_closed_form_fronts_in_both_geometries(self):
        for mode in (False, True):
            args = uniform_fixture(axisymmetric=mode)
            result = self.valid(**args)
            for step, row in enumerate(result.fronts):
                for j, point in enumerate(row):
                    self.assertAlmostEqual(point.x_m - args['axis_result'].point.x_m, step * .5 * sqrt(3), places=13)
                    self.assertEqual((point.radius_m, point.mach, point.flow_angle_deg), (j * .5, 2., 0.))
            self.assertEqual(result.relative_mass_flux, (1.,) * 4)
            self.assertEqual((result.max_relative_mass_drift, result.max_cell_relative_mass_imbalance, result.wall_mass_leakage_relative), (0., 0., 0.))

    def test_uniform_extreme_mach_does_not_reconstruct_from_saturated_speed(self):
        result = self.valid(**uniform_fixture(mach=1e200))
        self.assertTrue(all(p.mach == 1e200 for row in (*result.fronts, *result.half_fronts) for p in row))
        self.assertEqual(result.relative_mass_flux, (1.,) * 4)

    def test_large_translation_keeps_raw_residual_and_coordinate_allowance_visible(self):
        mach, height = 1.1, 1e295
        axis = solve_moc_axis_point(incoming=MocFlowPoint(1e308, height, mach, 0.),
            axis_history=(MocFlowPoint(0., 0., mach, 0.), MocFlowPoint(1., 0., mach, 0.)))
        x = axis.point.x_m
        result = self.valid(axis_result=axis, inlet_points=tuple(MocFlowPoint(x, height * j / 2, mach, 0.) for j in range(3)),
            wall_segment=MocWallSegment(x, height, x + 8 * height, height), marching_steps=3)
        self.assertGreater(result.max_equation_residual, 1e-3)
        self.assertGreater(result.max_coordinate_rounding_allowance, result.max_equation_residual)
        self.assertEqual(result.max_rounding_adjusted_equation_residual, 0.)
        self.assertEqual(result.relative_mass_flux, (1.,) * 4)

    def test_fixture_inputs_are_stable_at_110_and_180_digits(self):
        for mode, compress in itertools.product((False, True), repeat=2):
            self.assertEqual(radial_fixture(axisymmetric=mode, compressing=compress, precision=110)[0],
                             radial_fixture(axisymmetric=mode, compressing=compress, precision=180)[0])

    def test_nonuniform_flux_diagnostics_match_independent_quadrature(self):
        for mode in (False, True):
            result = self.valid(**radial_fixture(axisymmetric=mode)[0])
            with mp.workdps(110):
                g = exact_mp(result.inputs.axis_result.inputs.gamma)
                ref = exact_mp(result.inputs.axis_result.point.mach)
                def flux(row, *, g=g, ref=ref, mode=mode):
                    total = mp.mpf(0)
                    for a, b in itertools.pairwise(row):
                        values = []
                        for p in (a, b):
                            m = exact_mp(p.mach)
                            density = ((2 + (g - 1) * ref * ref) / (2 + (g - 1) * m * m)) ** (1 / (g - 1))
                            speed = m * mp.sqrt((g + 1) / (2 + (g - 1) * m * m))
                            theta = exact_mp(p.flow_angle_deg) * mp.pi / 180
                            values.append((density * speed * mp.cos(theta), density * speed * mp.sin(theta)))
                        ax, ar, bx, br = map(exact_mp, (a.x_m, a.radius_m, b.x_m, b.radius_m))
                        def integrand(t, *, values=values, ax=ax, ar=ar, bx=bx, br=br, mode=mode):
                            u = values[0][0] * (1 - t) + values[1][0] * t
                            v = values[0][1] * (1 - t) + values[1][1] * t
                            weight = ar * (1 - t) + br * t if mode else 1
                            return weight * (u * (br - ar) - v * (bx - ax))
                        total += mp.quad(integrand, [0, 1])
                    return total
                expected = tuple(float(flux(row) / flux(result.fronts[0])) for row in result.fronts)
                cell_error = mp.mpf(0)
                for old, half, new in zip(result.fronts, result.half_fronts, result.fronts[1:]):
                    for cell in model._cells(old, half, new):
                        edges = [flux((a, b)) for a, b in zip(cell, (*cell[1:], cell[0]))]
                        cell_error = max(cell_error, abs(sum(edges)) / sum(map(abs, edges)))
            self.assertTrue(all(abs(a - b) < 3e-15 for a, b in zip(expected, result.relative_mass_flux)))
            self.assertLess(abs(result.max_relative_mass_drift - max(abs(v - 1) for v in expected)), 3e-15)
            self.assertGreater(result.max_relative_mass_drift, 1e-8)
            self.assertLess(abs(result.max_cell_relative_mass_imbalance - float(cell_error)), 3e-15)

    def test_complete_adjacent_parent_mesh_has_no_missing_or_skipped_fronts(self):
        result = self.valid(**radial_fixture()[0])
        self.assertEqual(tuple(map(len, result.fronts)), (5,) * 5)
        self.assertEqual(tuple(map(len, result.half_fronts)), (4,) * 4)
        self.assertEqual(result.fronts[0], result.inputs.inlet_points)
        self.assertTrue(all(row[0].radius_m == row[0].flow_angle_deg == 0. for row in result.fronts))
        self.assertGreater(result.fronts[-1][0].mach, result.fronts[0][0].mach)

    def test_all_cells_tile_the_entire_strip_without_internal_boundary_holes(self):
        for args in (uniform_fixture(), radial_fixture()[0]):
            result = self.valid(**args)
            for old, half, new in zip(result.fronts, result.half_fronts, result.fronts[1:]):
                cells = model._cells(old, half, new)
                n = len(half)
                self.assertEqual(len(cells), 3 * n + 1)
                def area(points):
                    return abs(sum(Fraction(a.x_m) * Fraction(b.radius_m) - Fraction(b.x_m) * Fraction(a.radius_m)
                                   for a, b in zip(points, (*points[1:], points[0])))) / 2
                outline = (old[0], *new, old[-1], *old[-2:0:-1])
                self.assertEqual(sum(map(area, cells)), area(outline))
                edges = Counter((a, b) for cell in cells for a, b in zip(cell, (*cell[1:], cell[0])))
                net = edges - Counter({(b, a): count for (a, b), count in edges.items()})
                boundary = Counter(zip(outline, (*outline[1:], outline[0])))
                self.assertEqual(net, boundary)

    def test_unenforced_axis_departure_defect_stays_separate_from_solved_equations(self):
        args = radial_fixture()[0]
        result = self.valid(**args)
        first = solve_moc_axis_departure(axis_result=args['axis_result'], minus_parent=args['inlet_points'][1])
        self.assertGreater(result.max_axis_departure_plus_defect, 1e-8)
        self.assertGreaterEqual(result.max_axis_departure_plus_defect, first.averaged_plus_compatibility_defect)
        self.assertLessEqual(result.max_rounding_adjusted_equation_residual, 64 * 2**-52)
        self.assertGreaterEqual(result.max_equation_residual, result.max_rounding_adjusted_equation_residual)
        self.assertGreaterEqual(result.max_coordinate_rounding_allowance, 0.)

    def test_gamma_and_cylindrical_mode_come_from_owned_axis(self):
        for mode in (False, True):
            args, exact = radial_fixture(axisymmetric=mode, gamma=1.22, mach=3.)
            result = self.valid(**args)
            errors = [abs(p.mach / exact(p.x_m, p.radius_m)[0] - 1) for row in result.fronts for p in row]
            self.assertLess(max(errors), 1e-4)
            self.assertEqual(result.inputs.axis_result.inputs.gamma, 1.22)

    def test_power_of_two_coordinate_scaling_preserves_states_and_mass_ratios(self):
        args = radial_fixture()[0]
        result = self.valid(**args)
        for exponent in (-900, 900):
            def scale(p, *, exponent=exponent):
                return replace(p, x_m=ldexp(p.x_m, exponent), radius_m=ldexp(p.radius_m, exponent))
            axis = args['axis_result']
            initial = solve_moc_axis_point(incoming=scale(axis.inputs.incoming), axis_history=tuple(scale(p) for p in axis.inputs.axis_history))
            w = args['wall_segment']
            scaled = self.valid(axis_result=initial, inlet_points=tuple(scale(p) for p in args['inlet_points']),
                wall_segment=MocWallSegment(*(ldexp(v, exponent) for v in (w.x0_m, w.radius0_m, w.x1_m, w.radius1_m))), marching_steps=4)
            for original, changed in zip(result.fronts, scaled.fronts):
                self.assertEqual(tuple((p.mach, p.flow_angle_deg) for p in original), tuple((p.mach, p.flow_angle_deg) for p in changed))
                self.assertEqual(tuple((p.x_m, p.radius_m) for p in original), tuple((ldexp(p.x_m, -exponent), ldexp(p.radius_m, -exponent)) for p in changed))
            self.assertEqual(result.relative_mass_flux, scaled.relative_mass_flux)
            self.assertEqual(result.max_relative_mass_drift, scaled.max_relative_mass_drift)

    def test_mass_audit_is_independent_of_callers_decimal_context(self):
        args = radial_fixture()[0]
        expected = self.valid(**args)
        model._solve.cache_clear()
        model._flux_state.cache_clear()
        with localcontext() as ctx:
            ctx.prec = 7
            ctx.rounding = ROUND_DOWN
            ctx.traps[Inexact] = True
            self.assertEqual(expected, self.valid(**args))

    def test_axis_history_is_required_not_a_bare_point_or_wall_result(self):
        args = uniform_fixture()
        for bad in (args['axis_result'].point, {}, None):
            with self.subTest(bad=bad), self.assertRaises(DomainError):
                solve_moc_field(**{**args, 'axis_result': bad})

    def test_inlet_must_match_axis_order_vertical_section_and_wall(self):
        args = uniform_fixture()
        original = args['inlet_points']
        for points in (original[:2], original[::-1], (replace(original[0], mach=3.), *original[1:]),
                       (original[0], replace(original[1], x_m=1.), original[2]),
                       (original[0], original[1], original[1]),
                       (*original[:-1], replace(original[-1], flow_angle_deg=5.))):
            with self.subTest(points=points), self.assertRaises(DomainError):
                solve_moc_field(**{**args, 'inlet_points': points})
        with self.assertRaises(DomainError):
            solve_moc_field(**{**args, 'wall_segment': replace(args['wall_segment'], radius0_m=2.)})

    def test_inlet_resource_read_is_bounded_for_infinite_iterables(self):
        args = uniform_fixture()
        consumed = []
        def infinite():
            for p in itertools.repeat(args['inlet_points'][0]):
                consumed.append(p)
                yield p
        with self.assertRaises(DomainError):
            solve_moc_field(**{**args, 'inlet_points': infinite()})
        self.assertEqual(len(consumed), 66)

    def test_incompatible_inlet_wall_is_refused_before_marching(self):
        args = uniform_fixture()
        points = (*args['inlet_points'][:-1], replace(args['inlet_points'][-1], flow_angle_deg=5.))
        with self.assertRaises(DomainError) as context:
            MocFieldInputs(**{**args, 'inlet_points': points})
        self.assertEqual(context.exception.field, 'inlet_points')

    def test_resource_fields_reject_booleans_floats_and_excessive_work(self):
        args = uniform_fixture()
        for name, values in (('marching_steps', (True, 1., 0, 129)), ('max_iterations', (True, 1., 0, 257))):
            for value in values:
                with self.subTest(name=name, value=value), self.assertRaises(DomainError):
                    MocFieldInputs(**{**args, name: value})
        axis = args['axis_result'].point
        many = tuple(replace(axis, radius_m=j / 64) for j in range(65))
        with self.assertRaises(DomainError):
            MocFieldInputs(**{**args, 'inlet_points': many, 'marching_steps': 128})

    def test_wall_end_and_iteration_budget_do_not_extrapolate_or_silently_pass(self):
        args = radial_fixture()[0]
        wall = args['wall_segment']
        short = replace(wall, x1_m=wall.x0_m + .001, radius1_m=wall.radius0_m + .00005)
        for change, code in (({'max_iterations': 1}, 'moc_departure_not_converged'),
                             ({'wall_segment': short}, 'moc_wall_outside_segment')):
            with self.subTest(change=change):
                with self.assertRaises(ConvergenceError) as context:
                    solve_moc_field(**{**args, **change})
                self.assertEqual(context.exception.code, code)

    def test_cell_orientation_rejects_folds_and_collapses_at_extreme_scales(self):
        base = tuple(MocFlowPoint(x, r, 2., 0.) for x, r in ((0., .5), (.5, 0.), (1., .5), (.5, 1.)))
        for exponent in (-1000, 0, 1000):
            cell = tuple(replace(p, x_m=ldexp(p.x_m, exponent), radius_m=ldexp(p.radius_m, exponent)) for p in base)
            try:
                model._check_cell(cell)
            except DomainError as error:
                self.fail(f'Valid convex field cell refused: {error}')
            for bad in (cell[::-1], (cell[0], cell[0], *cell[2:])):
                with self.subTest(exponent=exponent):
                    with self.assertRaises(ConvergenceError) as context:
                        model._check_cell(bad)
                    self.assertEqual(context.exception.code, 'moc_field_folded_cell')

    def test_front_predicates_reject_axis_wall_crossings_and_radial_disorder(self):
        inputs = MocFieldInputs(**uniform_fixture())
        a, b = MocFlowPoint(1., .25, 2., 0.), MocFlowPoint(1., .75, 2., 0.)
        model._check_front((a, b), inputs, complete=False)
        for points, code in (((b, a), 'moc_field_front_order'),
                             ((replace(a, radius_m=0.), b), 'moc_field_outside_domain'),
                             ((a, replace(b, radius_m=1.1)), 'moc_field_outside_domain')):
            with self.subTest(points=points):
                with self.assertRaises(ConvergenceError) as context:
                    model._check_front(points, inputs, complete=False)
                self.assertEqual(context.exception.code, code)

    def test_noncausal_compressive_inlet_cannot_return_a_blessed_field(self):
        args = uniform_fixture()
        axis = args['axis_result'].point
        profiles = ((2.526999047132291, -1.0610771786951076),
                    (3.7129089430998308, 3.254856083683798),
                    (2.1769920662203983, -6.962304860855625))
        points = (axis, *(replace(axis, radius_m=(j + 1) / 4, mach=m, flow_angle_deg=t)
                          for j, (m, t) in enumerate(profiles)), replace(axis, radius_m=1.))
        with self.assertRaises(ConvergenceError) as context:
            solve_moc_field(**{**args, 'inlet_points': points, 'marching_steps': 6})
        self.assertEqual(context.exception.code, 'moc_unit_noncausal')

    def test_density_ledger_resource_limit_is_explicit_instead_of_erasing_flux(self):
        with self.assertRaises(ConvergenceError) as context:
            model._flux_state(MocFlowPoint(0., 1., 1e200, 0.), nextafter(1., 2.), 2.)
        self.assertEqual(context.exception.code, 'moc_field_density_resource')

    def test_result_mapping_owns_every_nested_front_and_input(self):
        result = self.valid(**radial_fixture()[0])
        raw = asdict(result)
        copied = MocFieldResult(**raw)
        self.assertEqual(copied, result)
        raw['fronts'][1][1]['mach'] = 100.
        raw['inputs']['axis_result']['inputs']['axis_history'][0]['mach'] = 200.
        self.assertNotEqual(copied.fronts[1][1].mach, 100.)
        self.assertNotEqual(copied.inputs.axis_result.inputs.axis_history[0].mach, 200.)

    def test_state_flux_and_cell_diagnostics_must_replay(self):
        result = self.valid(**radial_fixture()[0])
        rows = list(result.fronts)
        rows[-1] = (replace(rows[-1][0], mach=3.), *rows[-1][1:])
        changes = ({'fronts': tuple(rows)}, {'fronts': result.fronts[:-1]}, {'half_fronts': ()},
            {'relative_mass_flux': (1.,) * 5}, {'max_relative_mass_drift': 0.},
            {'max_cell_relative_mass_imbalance': 0.}, {'wall_mass_leakage_relative': 1.},
            {'max_equation_residual': 1.}, {'max_rounding_adjusted_equation_residual': 1.},
            {'max_coordinate_rounding_allowance': 1.}, {'max_axis_departure_plus_defect': 0.})
        for change in changes:
            with self.subTest(change=change), self.assertRaises(DomainError):
                replace(result, **change)

    def test_qualification_and_discretization_flags_cannot_be_forged(self):
        result = self.valid(**uniform_fixture())
        for name in ('physical_model_qualified', 'sonic_initialization_complete', 'full_nozzle_solution'):
            for value in (True, 0, 1):
                with self.subTest(name=name, value=value), self.assertRaises(DomainError):
                    replace(result, **{name: value})
        for change in ({'discretization_error_bound': 0.}, {'numerical_policy': 'certified'}, {'mass_flux_policy': 'conservative_by_construction'}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                replace(result, **change)

    def test_typed_input_record_replays_all_supersonic_states(self):
        args = uniform_fixture()
        copied = MocFieldInputs(**asdict(MocFieldInputs(**args)))
        self.assertEqual(copied.inlet_points, args['inlet_points'])
        with self.assertRaises(DomainError):
            solve_moc_field(**{**args, 'inlet_points': (args['inlet_points'][0], {'x_m': 0., 'radius_m': .5, 'mach': 1., 'flow_angle_deg': 0.}, args['inlet_points'][-1])})
