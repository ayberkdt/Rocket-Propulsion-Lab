"""Parametric propulsion-component geometry models.

The package provides axisymmetric converging-diverging nozzle contours
(conical, cubic-Hermite preliminary bell [1]_ and Rao-style thrust-optimized
parabola [2]_), a planar minimum-length method-of-characteristics nozzle [3]_, local
characteristic unit processes (interior, axis, wall and axis departure) and a
marched supersonic characteristic field, Hall's transonic throat expansion,
contour comparison metrics and text exports (CSV, SVG, DXF, ASCII STL). Each
public callable documents its own assumptions, validity limits and sources.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] RAO-1958. G. V. R. Rao, "Exhaust Nozzle Contour for Optimum Thrust",
       Jet Propulsion 28(6), 377-382, 1958. doi:10.2514/8.7324
.. [3] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
       compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
"""

from .comparison import ContourComparison, compare_nozzle_contours
from .export import ContourExportBundle, ContourExportInputs, export_nozzle_contour
from .hall_throat import (
    HallSonicInputs,
    HallSonicResult,
    HallThroatInputs,
    HallThroatModel,
    HallThroatResult,
    evaluate_hall_throat,
    locate_hall_sonic_point,
)
from .moc import (
    CharacteristicPoint,
    MocDesign,
    MocResiduals,
    generate_moc_nozzle,
)
from .moc_axis_departure import (
    MocAxisDepartureInputs,
    MocAxisDepartureResult,
    solve_moc_axis_departure,
)
from .moc_boundary import (
    MocAxisInputs,
    MocBoundaryResult,
    MocWallInputs,
    MocWallSegment,
    solve_moc_axis_point,
    solve_moc_wall_point,
)
from .moc_field import MocFieldInputs, MocFieldResult, solve_moc_field
from .moc_unit_process import MocFlowPoint, MocUnitInputs, MocUnitResult, solve_moc_interior_point
from .nozzle_contour import NozzleContour, NozzleStation, generate_nozzle_contour
from .rao import generate_rao_contour

__all__ = [
    "CharacteristicPoint",
    "ContourComparison",
    "ContourExportBundle",
    "ContourExportInputs",
    "HallSonicInputs",
    "HallSonicResult",
    "HallThroatInputs",
    "HallThroatModel",
    "HallThroatResult",
    "MocAxisDepartureInputs",
    "MocAxisDepartureResult",
    "MocAxisInputs",
    "MocBoundaryResult",
    "MocDesign",
    "MocFieldInputs",
    "MocFieldResult",
    "MocFlowPoint",
    "MocResiduals",
    "MocUnitInputs",
    "MocUnitResult",
    "MocWallInputs",
    "MocWallSegment",
    "NozzleContour",
    "NozzleStation",
    "compare_nozzle_contours",
    "evaluate_hall_throat",
    "export_nozzle_contour",
    "generate_moc_nozzle",
    "generate_nozzle_contour",
    "generate_rao_contour",
    "locate_hall_sonic_point",
    "solve_moc_axis_departure",
    "solve_moc_axis_point",
    "solve_moc_field",
    "solve_moc_interior_point",
    "solve_moc_wall_point",
]

