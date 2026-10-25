"""External aerodynamics: lifting flow, viscous and hypersonic aerothermodynamics,
and launch-vehicle aerodynamic coefficients.

Each submodule documents its own equations, validity range and primary
sources: :mod:`.potential_flow`, :mod:`.airfoil_geometry`, :mod:`.thin_airfoil`,
:mod:`.panel_methods`, :mod:`.lifting_line` (incompressible and lifting flow);
:mod:`.boundary_layer`, :mod:`.couette`, :mod:`.aerothermal`, :mod:`.entry`
(viscous flow and heating); :mod:`.slender_body`, :mod:`.launch_vehicle`,
:mod:`.aero_database` (ascent-vehicle aerodynamics).
"""

from . import (
    aero_database,
    aerothermal,
    airfoil_geometry,
    boundary_layer,
    couette,
    entry,
    launch_vehicle,
    lifting_line,
    panel_methods,
    potential_flow,
    slender_body,
    thin_airfoil,
)

__all__ = [
    "aero_database",
    "aerothermal",
    "airfoil_geometry",
    "boundary_layer",
    "couette",
    "entry",
    "launch_vehicle",
    "lifting_line",
    "panel_methods",
    "potential_flow",
    "slender_body",
    "thin_airfoil",
]
