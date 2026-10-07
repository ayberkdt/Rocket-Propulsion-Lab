"""Body-frame thruster geometry and propellant-dependent mass properties.

Translational propagation only needs the net force along a commanded axis.
Attitude-coupled (6-DOF) simulation, momentum unloading and thrust-vector
alignment also need *where* each thruster pushes and *where* the centre of
mass is while propellant drains.  This module supplies both in the vehicle
body frame:

``ThrusterLayout``
    Per-unit-thrust influence matrix (force and torque about the centre of
    mass), resultant force and torque for given thrust levels, cosine (axial)
    efficiency, and the torque derivative with respect to the centre of mass.
``PropellantMassModel``
    Total mass, centre of mass and inertia about the centre of mass from a dry
    rigid body plus point-mass propellant stations whose masses are the
    propagated tank states, with analytic derivatives of the centre of mass
    and inertia with respect to every tank mass.

Attitude dynamics and frame rotation stay in the consuming framework.  All
vectors here are body-frame and SI.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, real_mapping, real_tuple, text, vector3

Vector = tuple[float, float, float]
Matrix3 = tuple[Vector, Vector, Vector]

_ZERO: Vector = (0.0, 0.0, 0.0)


def _vector(value: tuple[float, ...], label: str) -> Vector:
    return vector3(value, label)


def _cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _sub(a: Vector, b: Vector) -> Vector:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _skew(a: Vector) -> Matrix3:
    """Return ``[a]x`` such that ``[a]x b = a x b``."""

    return ((0.0, -a[2], a[1]), (a[2], 0.0, -a[0]), (-a[1], a[0], 0.0))


def _shift_operator(r: Vector) -> Matrix3:
    """Return ``(r.r) I - r r^T`` (parallel-axis operator per unit mass)."""

    rr = _dot(r, r)
    return tuple(
        tuple((rr if i == j else 0.0) - r[i] * r[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def _matrix_add(*terms: tuple[float, Matrix3]) -> Matrix3:
    return tuple(
        tuple(sum(weight * matrix[i][j] for weight, matrix in terms) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class ThrusterMount:
    """One thruster's body-frame application point and force direction.

    ``thrust_axis`` is the direction of the force *on the vehicle* (opposite
    to the exhaust); it is normalized on construction.

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    thruster_id: str
    position_m: Vector
    thrust_axis: Vector

    def __post_init__(self) -> None:
        text(self.thruster_id, "Thruster identifier")
        object.__setattr__(self, "position_m", _vector(self.position_m, "Thruster position"))
        axis = _vector(self.thrust_axis, "Thrust axis")
        norm = sqrt(_dot(axis, axis))
        if norm <= 0.0:
            raise DomainError("Thrust axis must be non-zero.")
        object.__setattr__(self, "thrust_axis", tuple(value / norm for value in axis))


@dataclass(frozen=True, slots=True)
class ForceTorque:
    """Resultant body-frame force and torque about the centre of mass.

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    force_n: Vector
    torque_n_m: Vector
    total_thrust_n: float
    torque_wrt_center_of_mass: Matrix3


@dataclass(frozen=True, slots=True)
class ThrusterLayout:
    """Immutable set of body-mounted thrusters.

    For thruster ``i`` at ``p_i`` with axis ``u_i`` and thrust ``T_i``, about
    the centre of mass ``c``::

        F = sum T_i u_i,     tau = sum (p_i - c) x T_i u_i
        d tau / d c = [F]x

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    mounts: tuple[ThrusterMount, ...]

    def __post_init__(self) -> None:
        if isinstance(self.mounts, (str, bytes)) or not hasattr(self.mounts, "__iter__"):
            raise DomainError("Thruster mounts must be a sequence.")
        object.__setattr__(self, "mounts", tuple(self.mounts))
        if any(not isinstance(mount, ThrusterMount) for mount in self.mounts):
            raise DomainError("Thruster mounts must be ThrusterMount records.")
        if not self.mounts:
            raise DomainError("A thruster layout requires at least one thruster.")
        identifiers = [mount.thruster_id for mount in self.mounts]
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Thruster identifiers must be unique.")

    @property
    def thruster_ids(self) -> tuple[str, ...]:
        return tuple(mount.thruster_id for mount in self.mounts)

    def influence_matrix(self, center_of_mass_m: Vector) -> tuple[tuple[float, ...], ...]:
        """Return the 6xN force/torque produced by one newton on each thruster.

        Rows 0-2 are force, rows 3-5 torque about ``center_of_mass_m``;
        columns follow :attr:`thruster_ids`.  This is the matrix a thrust
        allocator inverts subject to non-negative, bounded thrust.

        References
        ----------
        B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
        """

        center = _vector(center_of_mass_m, "Centre of mass")
        columns = [
            (*mount.thrust_axis, *_cross(_sub(mount.position_m, center), mount.thrust_axis))
            for mount in self.mounts
        ]
        return tuple(tuple(column[row] for column in columns) for row in range(6))

    def resultant(
        self, thrusts_n: Mapping[str, float], center_of_mass_m: Vector
    ) -> ForceTorque:
        """Return net force, torque about the centre of mass and ``d tau/d c``.

        Thrusters not listed are off.  Thrust values must be finite and
        non-negative.

        References
        ----------
        B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
        """

        center = _vector(center_of_mass_m, "Centre of mass")
        thrusts_n = real_mapping(thrusts_n, "Thruster thrusts")
        unknown = set(thrusts_n) - set(self.thruster_ids)
        if unknown:
            raise DomainError("Unknown thrusters: " + ", ".join(sorted(unknown)))
        force = [0.0, 0.0, 0.0]
        torque = [0.0, 0.0, 0.0]
        total = 0.0
        for mount in self.mounts:
            thrust = thrusts_n.get(mount.thruster_id, 0.0)
            if thrust < 0.0:
                raise DomainError("Thruster thrust must be finite and non-negative.")
            if thrust == 0.0:
                continue
            total += thrust
            arm = _cross(_sub(mount.position_m, center), mount.thrust_axis)
            for index in range(3):
                force[index] += thrust * mount.thrust_axis[index]
                torque[index] += thrust * arm[index]
        net: Vector = (force[0], force[1], force[2])
        return ForceTorque(
            force_n=net,
            torque_n_m=(torque[0], torque[1], torque[2]),
            total_thrust_n=total,
            torque_wrt_center_of_mass=_skew(net),
        )

    def axial_efficiency(self, thrusts_n: Mapping[str, float], axis: Vector) -> float:
        """Return ``(F . a) / sum T_i``: the cosine loss of a canted cluster.

        References
        ----------
        NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
        """

        direction = _vector(axis, "Axis")
        norm = sqrt(_dot(direction, direction))
        if norm <= 0.0:
            raise DomainError("Axis must be non-zero.")
        result = self.resultant(thrusts_n, _ZERO)
        if result.total_thrust_n <= 0.0:
            raise DomainError("Axial efficiency requires at least one firing thruster.")
        return _dot(result.force_n, direction) / (norm * result.total_thrust_n)


@dataclass(frozen=True, slots=True)
class PropellantStation:
    """Point-mass location of one propagated tank state in the body frame.

    The propellant is lumped at a fixed centroid; slosh and centroid travel
    inside the tank are not represented.

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    state_name: str
    position_m: Vector

    def __post_init__(self) -> None:
        text(self.state_name, "Propellant station state name")
        object.__setattr__(self, "position_m", _vector(self.position_m, "Station position"))


@dataclass(frozen=True, slots=True)
class MassPropertiesEvaluation:
    """Mass, centre of mass, inertia and their tank-mass derivatives.

    ``inertia_kg_m2`` is about the current centre of mass in body axes.

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    mass_kg: float
    center_of_mass_m: Vector
    inertia_kg_m2: Matrix3
    center_of_mass_wrt_state: tuple[tuple[str, Vector], ...]
    inertia_wrt_state: tuple[tuple[str, Matrix3], ...]


@dataclass(frozen=True, slots=True)
class PropellantMassModel:
    """Dry rigid body plus point-mass propellant stations.

    With dry mass ``m_d`` at ``d``, dry inertia ``J_d`` about ``d``, and tank
    masses ``m_k`` at ``p_k``; ``M = m_d + sum m_k`` and
    ``c = (m_d d + sum m_k p_k) / M``.  Using ``S(r) = (r.r) I - r r^T``::

        I_c = J_d + m_d S(d) + sum m_k S(p_k) - M S(c)
        dc/dm_k   = (p_k - c) / M
        dI_c/dm_k = S(p_k) - S(c) - [2 (c.e) I - (e c^T + c e^T)],  e = p_k - c

    References
    ----------
    B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
    """

    dry_mass_kg: float
    dry_center_of_mass_m: Vector
    dry_inertia_kg_m2: Matrix3
    stations: tuple[PropellantStation, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "dry_mass_kg", real(self.dry_mass_kg, "Dry mass"))
        if self.dry_mass_kg <= 0.0:
            raise DomainError("Dry mass must be finite and greater than zero.")
        object.__setattr__(
            self, "dry_center_of_mass_m", _vector(self.dry_center_of_mass_m, "Dry centre of mass")
        )
        raw = self.dry_inertia_kg_m2
        if isinstance(raw, (str, bytes)) or not hasattr(raw, "__iter__"):
            raise DomainError("Dry inertia must be a 3x3 matrix.")
        inertia = tuple(real_tuple(row, "Dry inertia") for row in raw)
        if len(inertia) != 3 or any(len(row) != 3 for row in inertia):
            raise DomainError("Dry inertia must be a 3x3 matrix.")
        object.__setattr__(self, "dry_inertia_kg_m2", inertia)
        for i in range(3):
            for j in range(i):
                if abs(inertia[i][j] - inertia[j][i]) > 1e-9 * max(1.0, abs(inertia[i][j])):
                    raise DomainError("Dry inertia must be symmetric.")
        principal = (inertia[0][0], inertia[1][1], inertia[2][2])
        if any(value < 0.0 for value in principal) or any(
            principal[i] > principal[(i + 1) % 3] + principal[(i + 2) % 3] + 1e-9
            for i in range(3)
        ):
            raise DomainError("Dry inertia diagonal violates the triangle inequality.")
        if isinstance(self.stations, (str, bytes)) or not hasattr(self.stations, "__iter__"):
            raise DomainError("Propellant stations must be a sequence.")
        object.__setattr__(self, "stations", tuple(self.stations))
        if any(not isinstance(station, PropellantStation) for station in self.stations):
            raise DomainError("Propellant stations must be PropellantStation records.")
        names = [station.state_name for station in self.stations]
        if len(names) != len(set(names)):
            raise DomainError("Each tank state may have only one propellant station.")

    def evaluate(self, additional_states: Mapping[str, float]) -> MassPropertiesEvaluation:
        """Return mass properties for the current tank-state masses.

        References
        ----------
        B. Wie, Space Vehicle Dynamics and Control, 2nd ed., AIAA, 2008.
        """

        states = real_mapping(additional_states, "Tank states")
        masses: list[float] = []
        for station in self.stations:
            if station.state_name not in states:
                raise DomainError(f"Missing tank state {station.state_name!r}.")
            value = states[station.state_name]
            if value < 0.0:
                raise DomainError("Tank-state masses must be finite and non-negative.")
            masses.append(value)
        total = self.dry_mass_kg + sum(masses)
        moment = [self.dry_mass_kg * value for value in self.dry_center_of_mass_m]
        for mass, station in zip(masses, self.stations, strict=True):
            for index in range(3):
                moment[index] += mass * station.position_m[index]
        center: Vector = (moment[0] / total, moment[1] / total, moment[2] / total)
        terms: list[tuple[float, Matrix3]] = [
            (1.0, self.dry_inertia_kg_m2),
            (self.dry_mass_kg, _shift_operator(self.dry_center_of_mass_m)),
            (-total, _shift_operator(center)),
        ]
        terms.extend(
            (mass, _shift_operator(station.position_m))
            for mass, station in zip(masses, self.stations, strict=True)
        )
        inertia = _matrix_add(*terms)
        center_partials: list[tuple[str, Vector]] = []
        inertia_partials: list[tuple[str, Matrix3]] = []
        shift_center = _shift_operator(center)
        for station in self.stations:
            offset = _sub(station.position_m, center)
            center_partials.append(
                (station.state_name, tuple(value / total for value in offset))  # type: ignore[arg-type]
            )
            correction = tuple(
                tuple(
                    (2.0 * _dot(center, offset) if i == j else 0.0)
                    - (offset[i] * center[j] + center[i] * offset[j])
                    for j in range(3)
                )
                for i in range(3)
            )
            inertia_partials.append(
                (
                    station.state_name,
                    _matrix_add(
                        (1.0, _shift_operator(station.position_m)),
                        (-1.0, shift_center),
                        (-1.0, correction),  # type: ignore[arg-type]
                    ),
                )
            )
        return MassPropertiesEvaluation(
            mass_kg=total,
            center_of_mass_m=center,
            inertia_kg_m2=inertia,
            center_of_mass_wrt_state=tuple(center_partials),
            inertia_wrt_state=tuple(inertia_partials),
        )
