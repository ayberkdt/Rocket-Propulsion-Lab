"""Checked adapter direction outputs and stable serialized unit coordinates."""

from math import hypot, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import vector3
from rocket_propulsion.propulsion.burns.direction import _direction_geometry

_ADAPTER_DIRECTION_POLICY = "scaled_hypot_exact_binary_projection_v1"


def _directed_vector(direction, magnitude, *, label):
    """Project the final quantity, without rounding a unit vector first."""
    vector = vector3(direction, label)
    return _direction_geometry(vector).project(magnitude, label=label)


def _stored_unit_direction(direction, *, label):
    """Validate numerical unit-vector evidence without changing its coordinates.

    Four ulp of one is a binary64 representation check, not a physical pointing
    allowance. Repeated normalization can change the bytes covered by a hash.
    """
    vector = vector3(direction, label)
    if abs(hypot(*vector) - 1.0) > 4 * ulp(1.0):
        raise DomainError("Stored direction must already be a numerical unit vector.", field=label)
    return vector
