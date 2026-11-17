"""Rerun the sonic precision suites with the bounded first stage disabled.

The first stage answers most valid inputs, so these subclasses keep the
complete Decimal stage (its admissions, iterations and refusals) under the
same independent oracles. Base classes are referenced through their modules
so that unittest does not collect them twice.
"""
import unittest
from unittest.mock import patch

import tests.test_geometry_precision as geometry
import tests.test_isentropic_area_angle_inverse as area_angle
import tests.test_prandtl_meyer_precision as prandtl_meyer
from rocket_propulsion.compressible import sonic_numerics


class _CompleteStage:
    def setUp(self):
        super().setUp()
        patcher = patch.object(sonic_numerics, '_FAST_PATH', False)
        patcher.start()
        self.addCleanup(patcher.stop)


class CompleteStageGeometryPrecisionTests(_CompleteStage, geometry.GeometryPrecisionTests):
    pass


class CompleteStageAreaAngleInverseTests(_CompleteStage, area_angle.IsentropicAreaAngleInverseTests):
    pass


class CompleteStagePrandtlMeyerTests(_CompleteStage, prandtl_meyer.PrandtlMeyerPrecisionTests):
    pass


if __name__ == '__main__':
    unittest.main()
