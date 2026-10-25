"""Whole-field planar radial refinement with unchanged accuracy bounds."""
import unittest

from tests import test_moc_field as cases


class MocFieldPlanarAttacks(cases.MocFieldChecks, unittest.TestCase):
    def test_planar_expansion_whole_field_and_mass_refine(self):
        self.convergence(False, False)

    def test_planar_compression_whole_field_and_mass_refine(self):
        self.convergence(False, True)

