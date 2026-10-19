"""Body-frame force and moment composition for mounted, gimballed engines.

This module owns propulsion hardware geometry only.  It does not propagate
attitude, select a guidance direction, transform reference frames, or evaluate
spacecraft inertia.  An astrodynamics/6-DOF host supplies the current
body-frame reference point (normally the instantaneous centre of mass) and
uses the returned body-frame wrench in its own equations of motion.

The two gimbal rotations are deliberately explicit and order-sensitive.  For
each mount, axis 1 is supplied in body coordinates and axis 2 is derived as
``nominal_direction x axis_1``.  The realised direction is
``R(axis_2, angle_2) R(axis_1, angle_1) nominal_direction``.  This is a
body-fixed, extrinsic two-axis convention; adapters must not silently replace
it with another Euler/Cardan convention.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from math import cos, pi, sin, sqrt
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document, real, records, text, vector3

from .serialization import canonical_sha256


THRUST_VECTORING_SCHEMA = "rocket_propulsion_thrust_vectoring_v1"
THRUST_VECTORING_REFERENCE_IDS = (
    "NASA-SP-8059",
    "NASA-TM-20240012525",
)

Vector3 = tuple[float, float, float]


def _vector3(value: tuple[float, float, float], label: str) -> Vector3:
    return vector3(value, label)


def _dot(left: Vector3, right: Vector3) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def _cross(left: Vector3, right: Vector3) -> Vector3:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _add(left: Vector3, right: Vector3) -> Vector3:
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _subtract(left: Vector3, right: Vector3) -> Vector3:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _scale(scale: float, vector: Vector3) -> Vector3:
    return tuple(scale * component for component in vector)  # type: ignore[return-value]


def _norm(vector: Vector3) -> float:
    return sqrt(_dot(vector, vector))


def _normalized(vector: Vector3, label: str) -> Vector3:
    magnitude = _norm(vector)
    if magnitude <= 0.0:
        raise DomainError(f"{label} must be non-zero.")
    return _scale(1.0 / magnitude, vector)


def _rotate(vector: Vector3, axis: Vector3, angle_rad: float) -> Vector3:
    """Rodrigues rotation around an already-normalised, body-fixed axis."""

    cosine = cos(angle_rad)
    sine = sin(angle_rad)
    return _add(
        _add(_scale(cosine, vector), _scale(sine, _cross(axis, vector))),
        _scale((1.0 - cosine) * _dot(axis, vector), axis),
    )


@dataclass(frozen=True, slots=True)
class BodyFrameEngineMount:
    """Immutable engine location and two-axis gimbal geometry.

    ``nominal_direction_body`` is the direction of the force applied to the
    vehicle, not the exhaust direction.  ``application_point_body_m`` and both
    axes use the vehicle body frame.  Axis 1 must be perpendicular to the
    nominal direction; axis 2 is the normalised cross product
    ``nominal_direction x axis_1``.  Limits apply to the realised commanded
    angles and may be zero for a fixed engine.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    engine_id: str
    application_point_body_m: Vector3
    nominal_direction_body: Vector3
    gimbal_axis_1_body: Vector3
    gimbal_axis_1_limit_rad: float = 0.0
    gimbal_axis_2_limit_rad: float = 0.0

    def __post_init__(self) -> None:
        text(self.engine_id, "Engine mount identifier")
        point = _vector3(self.application_point_body_m, "Engine application point")
        direction = _normalized(
            _vector3(self.nominal_direction_body, "Nominal thrust direction"),
            "Nominal thrust direction",
        )
        axis_1 = _normalized(
            _vector3(self.gimbal_axis_1_body, "First gimbal axis"),
            "First gimbal axis",
        )
        if abs(_dot(direction, axis_1)) > 1e-10:
            raise DomainError("First gimbal axis must be perpendicular to nominal thrust.")
        for value, label in (
            (self.gimbal_axis_1_limit_rad, "First gimbal limit"),
            (self.gimbal_axis_2_limit_rad, "Second gimbal limit"),
        ):
            value = real(value, label)
            if value < 0.0 or value >= 0.5 * pi:
                raise DomainError(f"{label} must be finite and in [0, pi/2).")
        object.__setattr__(self, "application_point_body_m", point)
        object.__setattr__(self, "nominal_direction_body", direction)
        object.__setattr__(self, "gimbal_axis_1_body", axis_1)
        object.__setattr__(
            self,
            "gimbal_axis_1_limit_rad",
            real(self.gimbal_axis_1_limit_rad, "First gimbal limit"),
        )
        object.__setattr__(
            self,
            "gimbal_axis_2_limit_rad",
            real(self.gimbal_axis_2_limit_rad, "Second gimbal limit"),
        )

    @property
    def gimbal_axis_2_body(self) -> Vector3:
        """Derived second body-fixed gimbal axis."""

        return _normalized(
            _cross(self.nominal_direction_body, self.gimbal_axis_1_body),
            "Derived second gimbal axis",
        )

    def direction_and_partials(
        self,
        gimbal_axis_1_rad: float,
        gimbal_axis_2_rad: float,
    ) -> tuple[Vector3, Vector3, Vector3]:
        """Return realised direction and analytic derivatives by both angles."""

        angle_1 = real(gimbal_axis_1_rad, "First gimbal angle")
        angle_2 = real(gimbal_axis_2_rad, "Second gimbal angle")
        tolerance = 1e-14
        if abs(angle_1) > self.gimbal_axis_1_limit_rad + tolerance:
            raise DomainError(f"Engine {self.engine_id!r} first gimbal angle exceeds its limit.")
        if abs(angle_2) > self.gimbal_axis_2_limit_rad + tolerance:
            raise DomainError(f"Engine {self.engine_id!r} second gimbal angle exceeds its limit.")
        axis_1 = self.gimbal_axis_1_body
        axis_2 = self.gimbal_axis_2_body
        after_axis_1 = _rotate(self.nominal_direction_body, axis_1, angle_1)
        direction = _rotate(after_axis_1, axis_2, angle_2)
        derivative_1 = _rotate(_cross(axis_1, after_axis_1), axis_2, angle_2)
        derivative_2 = _cross(axis_2, direction)
        return direction, derivative_1, derivative_2


@dataclass(frozen=True, slots=True)
class EngineVectorCommand:
    """Instantaneous thrust magnitude and two body-fixed gimbal commands.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    engine_id: str
    thrust_n: float
    gimbal_axis_1_rad: float = 0.0
    gimbal_axis_2_rad: float = 0.0
    enabled: bool = True

    def __post_init__(self) -> None:
        text(self.engine_id, "Engine command identifier")
        thrust = real(self.thrust_n, "Commanded thrust")
        angle_1 = real(self.gimbal_axis_1_rad, "First commanded gimbal angle")
        angle_2 = real(self.gimbal_axis_2_rad, "Second commanded gimbal angle")
        if thrust < 0.0:
            raise DomainError("Commanded thrust must be finite and non-negative.")
        if not isinstance(self.enabled, bool):
            raise DomainError("Engine command enabled flag must be boolean.")
        object.__setattr__(self, "thrust_n", thrust)
        object.__setattr__(self, "gimbal_axis_1_rad", angle_1)
        object.__setattr__(self, "gimbal_axis_2_rad", angle_2)


@dataclass(frozen=True, slots=True)
class EngineWrenchPartials:
    """Analytic engine force/moment derivatives for estimation and 6-DOF STM use.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    force_wrt_thrust: Vector3
    moment_wrt_thrust: Vector3
    force_wrt_gimbal_axis_1_rad: Vector3
    force_wrt_gimbal_axis_2_rad: Vector3
    moment_wrt_gimbal_axis_1_rad: Vector3
    moment_wrt_gimbal_axis_2_rad: Vector3


@dataclass(frozen=True, slots=True)
class EngineWrench:
    """One mounted engine's body-frame force and moment contribution.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    engine_id: str
    active: bool
    thrust_n: float
    direction_body: Vector3
    force_body_n: Vector3
    moment_body_n_m: Vector3
    moment_arm_body_m: Vector3
    partials: EngineWrenchPartials


@dataclass(frozen=True, slots=True)
class ClusterWrenchPartials:
    """Named analytic derivatives of the resultant body-frame wrench.

    Each tuple is keyed by engine identifier.  Gimbal entries contain axis-1
    and axis-2 derivatives in that order.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    force_wrt_thrust: tuple[tuple[str, Vector3], ...]
    moment_wrt_thrust: tuple[tuple[str, Vector3], ...]
    force_wrt_gimbal_rad: tuple[tuple[str, tuple[Vector3, Vector3]], ...]
    moment_wrt_gimbal_rad: tuple[tuple[str, tuple[Vector3, Vector3]], ...]


@dataclass(frozen=True, slots=True)
class ClusterWrenchEvaluation:
    """Body-frame force/moment result about a caller-supplied reference point.

    ``scalar_thrust_sum_n - resultant_thrust_n`` is reported as vector
    cancellation loss.  It is not propulsive efficiency or plume loss; it is
    solely the geometric loss caused by summing non-collinear force vectors.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    reference_point_body_m: Vector3
    force_body_n: Vector3
    moment_body_n_m: Vector3
    resultant_thrust_n: float
    resultant_direction_body: Vector3 | None
    scalar_thrust_sum_n: float
    vector_cancellation_loss_n: float
    engines: tuple[EngineWrench, ...]
    partials: ClusterWrenchPartials
    model_sha256: str


@dataclass(frozen=True, slots=True)
class MountedEngineCluster:
    """Hash-identified body-frame geometry for a mounted engine cluster.

    The evaluation uses ``moment = (application_point - reference_point) x
    force``.  Passing the instantaneous body-frame centre of mass therefore
    yields the propulsion torque required by a host 6-DOF propagator.  The
    cluster neither computes that centre of mass nor integrates attitude.

    References
    ----------
    NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
    NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
    """

    mounts: tuple[BodyFrameEngineMount, ...]
    model_id: str = "mounted-engine-cluster"
    schema: str = THRUST_VECTORING_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "mounts", records(self.mounts, BodyFrameEngineMount, "Mounted engines")
        )
        text(self.model_id, "Mounted-cluster model identifier")
        text(self.schema, "Thrust-vectoring schema")
        if self.schema != THRUST_VECTORING_SCHEMA:
            raise DomainError("Unsupported thrust-vectoring schema.")
        if not self.mounts:
            raise DomainError("A mounted engine cluster requires at least one mount.")
        identifiers = tuple(mount.engine_id for mount in self.mounts)
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Mounted engine identifiers must be unique.")
        object.__setattr__(
            self, "mounts", tuple(sorted(self.mounts, key=lambda item: item.engine_id))
        )

    def _payload(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "model_id": self.model_id,
            "reference_ids": list(THRUST_VECTORING_REFERENCE_IDS),
            "rotation_convention": "direction=R_axis2(angle2)@R_axis1(angle1)@nominal",
            "mounts": [
                {
                    "engine_id": mount.engine_id,
                    "application_point_body_m": list(mount.application_point_body_m),
                    "nominal_direction_body": list(mount.nominal_direction_body),
                    "gimbal_axis_1_body": list(mount.gimbal_axis_1_body),
                    "gimbal_axis_1_limit_rad": mount.gimbal_axis_1_limit_rad,
                    "gimbal_axis_2_limit_rad": mount.gimbal_axis_2_limit_rad,
                }
                for mount in self.mounts
            ],
        }

    @property
    def model_sha256(self) -> str:
        """Canonical SHA-256 identity of the geometry and convention."""

        return canonical_sha256(self._payload())

    def model_manifest(self) -> dict[str, Any]:
        """Return the integration-facing provenance and semantics manifest."""

        payload = self._payload()
        return {
            "schema": self.schema,
            "model_id": self.model_id,
            "model_sha256": self.model_sha256,
            "reference_ids": payload["reference_ids"],
            "rotation_convention": payload["rotation_convention"],
            "force_frame": "vehicle_body",
            "moment_frame": "vehicle_body",
            "moment_reference": "caller_supplied_body_frame_point",
            "owns_attitude_propagation": False,
        }

    def to_dict(self) -> dict[str, Any]:
        """Serialize geometry with a recomputable semantic hash."""

        payload = self._payload()
        payload["model_sha256"] = self.model_sha256
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize geometry as deterministic, finite JSON."""

        return json.dumps(
            self.to_dict(),
            allow_nan=False,
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> MountedEngineCluster:
        """Re-open serialized geometry and refuse schema or hash tampering."""

        if not isinstance(payload, Mapping):
            raise InputError("Mounted engine cluster root must be an object.")
        if payload.get("schema") != THRUST_VECTORING_SCHEMA:
            raise InputError("Unsupported thrust-vectoring schema.")
        if payload.get("reference_ids") != list(THRUST_VECTORING_REFERENCE_IDS):
            raise InputError("Thrust-vectoring reference identity does not match schema v1.")
        if payload.get("rotation_convention") != (
            "direction=R_axis2(angle2)@R_axis1(angle1)@nominal"
        ):
            raise InputError("Thrust-vectoring rotation convention does not match schema v1.")
        try:
            mounts = tuple(
                BodyFrameEngineMount(
                    engine_id=item["engine_id"],
                    application_point_body_m=item["application_point_body_m"],
                    nominal_direction_body=item["nominal_direction_body"],
                    gimbal_axis_1_body=item["gimbal_axis_1_body"],
                    gimbal_axis_1_limit_rad=item["gimbal_axis_1_limit_rad"],
                    gimbal_axis_2_limit_rad=item["gimbal_axis_2_limit_rad"],
                )
                for item in payload["mounts"]
            )
            cluster = cls(
                mounts=mounts,
                model_id=payload["model_id"],
                schema=payload["schema"],
            )
        except (KeyError, TypeError, ValueError, DomainError) as exc:
            raise InputError("Malformed mounted engine cluster document.") from exc
        if payload.get("model_sha256") != cluster.model_sha256:
            raise InputError("Mounted engine cluster semantic hash mismatch.")
        return cluster

    @classmethod
    def from_json(cls, document: str) -> MountedEngineCluster:
        """Parse and integrity-check a JSON cluster document."""

        return read_json_document(document, "Mounted engine cluster", cls.from_dict)  # type: ignore[return-value]

    def evaluate(
        self,
        commands: tuple[EngineVectorCommand, ...],
        *,
        reference_point_body_m: Vector3 = (0.0, 0.0, 0.0),
    ) -> ClusterWrenchEvaluation:
        """Compose all engine forces, moments, and analytic command partials.

        Every configured engine must have exactly one command, including
        disabled or zero-thrust engines.  This fail-closed rule prevents an
        omitted engine from being mistaken for an intentional shutdown.

        References
        ----------
        NASA-SP-8059: https://ntrs.nasa.gov/citations/19710016722
        NASA/TM-20240012525: https://ntrs.nasa.gov/citations/20240012525
        """

        commands = records(commands, EngineVectorCommand, "Engine vector commands")
        reference = _vector3(reference_point_body_m, "Wrench reference point")
        command_map = {command.engine_id: command for command in commands}
        if len(command_map) != len(commands):
            raise DomainError("Each mounted engine requires at most one vector command.")
        mount_ids = {mount.engine_id for mount in self.mounts}
        if set(command_map) != mount_ids:
            missing = sorted(mount_ids - set(command_map))
            unknown = sorted(set(command_map) - mount_ids)
            details: list[str] = []
            if missing:
                details.append("missing=" + ",".join(missing))
            if unknown:
                details.append("unknown=" + ",".join(unknown))
            raise DomainError("Engine vector commands do not match mounts: " + "; ".join(details))

        total_force: Vector3 = (0.0, 0.0, 0.0)
        total_moment: Vector3 = (0.0, 0.0, 0.0)
        scalar_sum = 0.0
        engine_results: list[EngineWrench] = []
        force_wrt_thrust: list[tuple[str, Vector3]] = []
        moment_wrt_thrust: list[tuple[str, Vector3]] = []
        force_wrt_gimbal: list[tuple[str, tuple[Vector3, Vector3]]] = []
        moment_wrt_gimbal: list[tuple[str, tuple[Vector3, Vector3]]] = []

        for mount in self.mounts:
            command = command_map[mount.engine_id]
            direction, direction_axis_1, direction_axis_2 = mount.direction_and_partials(
                command.gimbal_axis_1_rad,
                command.gimbal_axis_2_rad,
            )
            thrust = command.thrust_n if command.enabled else 0.0
            force = _scale(thrust, direction)
            arm = _subtract(mount.application_point_body_m, reference)
            moment = _cross(arm, force)
            force_axis_1 = _scale(thrust, direction_axis_1)
            force_axis_2 = _scale(thrust, direction_axis_2)
            moment_thrust = _cross(arm, direction)
            moment_axis_1 = _cross(arm, force_axis_1)
            moment_axis_2 = _cross(arm, force_axis_2)
            partials = EngineWrenchPartials(
                force_wrt_thrust=direction if command.enabled else (0.0, 0.0, 0.0),
                moment_wrt_thrust=moment_thrust if command.enabled else (0.0, 0.0, 0.0),
                force_wrt_gimbal_axis_1_rad=force_axis_1,
                force_wrt_gimbal_axis_2_rad=force_axis_2,
                moment_wrt_gimbal_axis_1_rad=moment_axis_1,
                moment_wrt_gimbal_axis_2_rad=moment_axis_2,
            )
            engine_results.append(
                EngineWrench(
                    engine_id=mount.engine_id,
                    active=bool(command.enabled and thrust > 0.0),
                    thrust_n=thrust,
                    direction_body=direction,
                    force_body_n=force,
                    moment_body_n_m=moment,
                    moment_arm_body_m=arm,
                    partials=partials,
                )
            )
            total_force = _add(total_force, force)
            total_moment = _add(total_moment, moment)
            scalar_sum += thrust
            force_wrt_thrust.append((mount.engine_id, partials.force_wrt_thrust))
            moment_wrt_thrust.append((mount.engine_id, partials.moment_wrt_thrust))
            force_wrt_gimbal.append(
                (
                    mount.engine_id,
                    (
                        partials.force_wrt_gimbal_axis_1_rad,
                        partials.force_wrt_gimbal_axis_2_rad,
                    ),
                )
            )
            moment_wrt_gimbal.append(
                (
                    mount.engine_id,
                    (
                        partials.moment_wrt_gimbal_axis_1_rad,
                        partials.moment_wrt_gimbal_axis_2_rad,
                    ),
                )
            )

        resultant = _norm(total_force)
        cancellation = scalar_sum - resultant
        tolerance = max(1e-12, scalar_sum * 1e-12)
        if cancellation < -tolerance:
            raise DomainError("Resultant thrust exceeds scalar thrust sum.")
        cancellation = max(0.0, cancellation)
        resultant_direction = None if resultant <= 1e-15 else _scale(1.0 / resultant, total_force)
        return ClusterWrenchEvaluation(
            reference_point_body_m=reference,
            force_body_n=total_force,
            moment_body_n_m=total_moment,
            resultant_thrust_n=resultant,
            resultant_direction_body=resultant_direction,
            scalar_thrust_sum_n=scalar_sum,
            vector_cancellation_loss_n=cancellation,
            engines=tuple(engine_results),
            partials=ClusterWrenchPartials(
                force_wrt_thrust=tuple(force_wrt_thrust),
                moment_wrt_thrust=tuple(moment_wrt_thrust),
                force_wrt_gimbal_rad=tuple(force_wrt_gimbal),
                moment_wrt_gimbal_rad=tuple(moment_wrt_gimbal),
            ),
            model_sha256=self.model_sha256,
        )


__all__ = [
    "THRUST_VECTORING_REFERENCE_IDS",
    "THRUST_VECTORING_SCHEMA",
    "BodyFrameEngineMount",
    "ClusterWrenchEvaluation",
    "ClusterWrenchPartials",
    "EngineVectorCommand",
    "EngineWrench",
    "EngineWrenchPartials",
    "MountedEngineCluster",
    "Vector3",
]
