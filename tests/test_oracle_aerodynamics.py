"""Incompressible aerodynamics against exact results and internal-consistency laws.

* Lifting-line on an elliptic planform has the closed-form answers
  ``CL = a0 alpha / (1 + a0/(pi AR))``, constant downwash, span efficiency 1
  and ``CDi = CL^2/(pi AR)``; any planform must satisfy ``e <= 1``.
* The vortex-panel solution must give the same lift from the circulation
  (Kutta-Joukowski, ``cl = 2 Gamma/(V c)``) and from integrating its own
  pressure coefficients, produce no pressure drag in 2-D potential flow
  (d'Alembert), satisfy the Kutta condition (equal upper and lower
  trailing-edge velocities), converge as the panel count doubles, and tend to
  thin-airfoil theory ``cl = 2 pi alpha`` as thickness goes to zero; a
  symmetric section must give zero lift at zero incidence and antisymmetric
  lift in alpha.
"""

from __future__ import annotations

import itertools
import math
import unittest

from rocket_propulsion.aerodynamics import airfoil_geometry as ag
from rocket_propulsion.aerodynamics import lifting_line as ll
from rocket_propulsion.aerodynamics import panel_methods as pm


class EllipticWingExact(unittest.TestCase):
    def test_closed_form_lift_downwash_and_drag(self) -> None:
        for span, chord, a0 in ((10.0, 1.5, 2 * math.pi), (4.0, 1.0, 5.7), (30.0, 2.0, 6.1)):
            wing = ll.elliptic_wing(span, chord, section_lift_slope=a0)
            for alpha_deg in (-3.0, 2.0, 7.0):
                alpha = math.radians(alpha_deg)
                r = ll.solve_lifting_line(wing, alpha)
                ar = r.aspect_ratio
                cl = a0 * alpha / (1 + a0 / (math.pi * ar))
                with self.subTest(span=span, chord=chord, alpha=alpha_deg):
                    self.assertAlmostEqual(ar, span * span / (math.pi * span * chord / 4), places=10)
                    self.assertLess(abs(r.lift_coefficient - cl), 1e-12)
                    self.assertLess(abs(r.induced_drag_coefficient - cl * cl / (math.pi * ar)), 1e-13)
                    self.assertLess(abs(r.span_efficiency - 1.0), 1e-12)
                    self.assertLess(max(r.induced_angles_rad) - min(r.induced_angles_rad), 1e-15)

    def test_every_planform_has_efficiency_at_most_one(self) -> None:
        for taper in (0.1, 0.3, 0.45, 0.7, 1.0):
            for twist in (-4.0, 0.0, 3.0):
                wing = ll.trapezoidal_wing(12.0, 2.0, taper, tip_twist_rad=math.radians(twist))
                r = ll.solve_lifting_line(wing, math.radians(5.0))
                with self.subTest(taper=taper, twist=twist):
                    self.assertLessEqual(r.span_efficiency, 1.0 + 1e-12)
                    # CDi = CL^2/(pi AR e) is the definition of e
                    self.assertLess(abs(r.induced_drag_coefficient * math.pi * r.aspect_ratio * r.span_efficiency
                                        - r.lift_coefficient ** 2), 1e-12)

    def test_lift_is_linear_and_symmetric_spanwise(self) -> None:
        wing = ll.trapezoidal_wing(9.0, 1.6, 0.4)
        a = ll.solve_lifting_line(wing, math.radians(2.0))
        b = ll.solve_lifting_line(wing, math.radians(6.0))
        self.assertLess(abs(b.lift_coefficient - 3 * a.lift_coefficient), 1e-12)
        g = a.circulation
        self.assertLess(max(abs(x - y) for x, y in zip(g, reversed(g), strict=True)), 1e-14)


def solve(designation: str, alpha_deg: float, panels: int = 120) -> pm.LiftingPanelSolution:
    nodes = ag.naca_section(designation).panel_nodes(panels)
    return pm.solve_lifting_panels(nodes, math.radians(alpha_deg))


class VortexPanelConsistency(unittest.TestCase):
    def test_kutta_joukowski_equals_pressure_integration_and_no_drag(self) -> None:
        for designation in ("0006", "0012", "2412", "4415", "23012"):
            for alpha in (-4.0, 0.0, 3.0, 8.0):
                r = solve(designation, alpha)
                with self.subTest(section=designation, alpha=alpha):
                    self.assertLess(abs(r.lift_coefficient - 2 * r.circulation / r.chord), 1e-12)
                    # pressure integration is O(h); 1 % at 120 panels per surface
                    self.assertLess(abs(r.pressure_lift_coefficient - r.lift_coefficient),
                                    1e-2 * max(1.0, abs(r.lift_coefficient)))
                    self.assertLess(abs(r.pressure_drag_coefficient), 2e-3)

    def test_pressure_lift_gap_shrinks_with_refinement(self) -> None:
        for designation in ("2412", "4415"):
            gaps = [abs(r.pressure_lift_coefficient - r.lift_coefficient)
                    for r in (solve(designation, 8.0, n) for n in (60, 120, 240))]
            with self.subTest(section=designation):
                self.assertLess(gaps[1], 0.75 * gaps[0])
                self.assertLess(gaps[2], 0.75 * gaps[1])

    def test_cambered_sections_match_an_independent_hess_smith_value(self) -> None:
        """Regression for the unsigned-distance defect: NACA 4412 cl0 was 1.44 with 80 panels.

        0.5172 is the zero-incidence lift of an independent numpy Hess-Smith
        code on the same 80-panel-per-surface nodes; thin-airfoil theory gives
        0.4556 and thickness adds about 0.1 t/c of it.
        """
        self.assertLess(abs(solve("4412", 0.0, 80).lift_coefficient - 0.5172), 5e-4)
        slope = (solve("4412", 4.0, 160).lift_coefficient - solve("4412", 0.0, 160).lift_coefficient) / math.radians(4.0)
        self.assertLess(abs(slope - 6.904), 0.02)

    def test_kutta_condition(self) -> None:
        for designation in ("0012", "2412"):
            r = solve(designation, 5.0)
            v = r.surface_velocities
            self.assertLess(abs(abs(v[0]) - abs(v[-1])), 1e-9 * max(abs(v[0]), 1.0))

    def test_symmetric_section_is_antisymmetric_in_alpha(self) -> None:
        self.assertLess(abs(solve("0012", 0.0).lift_coefficient), 1e-12)
        for alpha in (1.0, 5.0, 9.0):
            self.assertLess(abs(solve("0012", alpha).lift_coefficient + solve("0012", -alpha).lift_coefficient), 1e-11)

    def test_panel_refinement_converges(self) -> None:
        values = [solve("2412", 4.0, n).lift_coefficient for n in (40, 80, 160, 320)]
        diffs = [abs(b - a) for a, b in itertools.pairwise(values)]
        self.assertTrue(diffs[1] < diffs[0] and diffs[2] < diffs[1])
        self.assertLess(diffs[-1] / abs(values[-1]), 1e-3)

    def test_thin_limit_approaches_thin_airfoil_theory(self) -> None:
        alpha = math.radians(4.0)
        errors = []
        for designation in ("0024", "0012", "0006", "0003"):
            r = solve(designation, 4.0, 200)
            errors.append(abs(r.lift_coefficient - 2 * math.pi * alpha) / (2 * math.pi * alpha))
        # thickness raises cl by about 0.8 t/c; the gap must shrink toward zero
        self.assertEqual(errors, sorted(errors, reverse=True))
        self.assertLess(errors[-1], 0.035)
