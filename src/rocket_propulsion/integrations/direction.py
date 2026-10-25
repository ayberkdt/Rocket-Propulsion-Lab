"""Checked direction outputs for the astrodynamics adapters.

Adapters project a final quantity (a delta-v or a unit direction) along a
validated direction without rounding an intermediate unit vector, and accept
stored unit directions only when their serialized coordinates are already a
numerical unit vector, so that hashed coordinates are never rewritten.

References
----------
.. [1] PYTHON-3.12-HYPOT. Python 3.12, ``math.hypot``,
       https://docs.python.org/3.12/library/math.html#math.hypot
"""

from math import hypot, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import vector3
from rocket_propulsion.propulsion.burns.direction import _direction_geometry

_ADAPTER_DIRECTION_POLICY = "scaled_hypot_exact_binary_projection_v1"


def _directed_vector(direction, magnitude, *, label):
    """Return ``magnitude * d / |d|`` as three finite floats.

    Parameters
    ----------
    direction : object
        Three finite real components, not all zero; any length.
    magnitude : float
        Quantity to project (for example a delta-v magnitude, or 1 for a
        unit direction).
    label : str
        Name used in error messages.

    Returns
    -------
    tuple of float
        The projected vector; each component is rounded once from an exact
        rational expression, without first rounding a unit vector.

    Raises
    ------
    DomainError
        If ``direction`` is malformed or zero, or a component is not a
        finite float.
    """
    vector = vector3(direction, label)
    return _direction_geometry(vector).project(magnitude, label=label)


def _stored_unit_direction(direction, *, label):
    """Accept a stored unit direction only if it is already numerically unit length.

    Parameters
    ----------
    direction : object
        Three finite real components read from a serialized document.
    label : str
        Name used in error messages.

    Returns
    -------
    tuple of float
        The coordinates exactly as stored.

    Raises
    ------
    DomainError
        If the components are malformed or ``abs(hypot(*d) - 1)`` exceeds four
        units in the last place of 1.0 [1]_.

    Notes
    -----
    The tolerance of four ulp of one is a binary64 representation check, not
    a physical pointing allowance. Coordinates are never renormalized,
    because renormalization can change the bytes covered by an integrity
    hash.

    References
    ----------
    .. [1] PYTHON-3.12-HYPOT. Python 3.12, ``math.hypot``,
           https://docs.python.org/3.12/library/math.html#math.hypot
    """
    vector = vector3(direction, label)
    if abs(hypot(*vector) - 1.0) > 4 * ulp(1.0):
        raise DomainError("Stored direction must already be a numerical unit vector.", field=label)
    return vector
