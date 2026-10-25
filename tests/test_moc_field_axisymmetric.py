"""Whole-field axisymmetric radial refinement with unchanged accuracy bounds."""
import unittest

from tests import test_moc_field as cases


class MocFieldAxisymmetricAttacks(cases.MocFieldChecks, unittest.TestCase):
    def test_spherical_expansion_whole_field_and_mass_refine(self):
        self.convergence(True, False)

    def test_spherical_compression_whole_field_and_mass_refine(self):
        self.convergence(True, True)

