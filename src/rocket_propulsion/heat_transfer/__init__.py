"""Heat transfer: conduction, convection and radiation for propulsion hardware.

Each submodule documents its own equations, validity range and primary
sources: :mod:`.conduction` (steady and transient walls, fins),
:mod:`.convection` (internal, external, natural, boiling, film cooling and
rocket coolant correlations), :mod:`.radiation` (black-body, view factors,
gray enclosures, gas radiation) and :mod:`.regen_chamber` (coupled
gas-wall-coolant regenerative cooling).
"""

from . import conduction, convection, radiation, regen_chamber

__all__ = ["conduction", "convection", "radiation", "regen_chamber"]
