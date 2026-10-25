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
it with another Euler/Cardan convention.  Order-sensitive two-axis gimbal
geometry is derived in [2]_; thrust-vector pointing, misalignment and
centre-of-mass moment arms are treated in [1]_.

References
----------
.. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
       Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
.. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
       Vector Control Actuator-Force / Gimbal-Torque Transformation
       Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from math import cos, pi, sin, sqrt
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import real, records, text, vector3

from .exchange_documents import _read_document
from .serialization import canonical_sha256

THRUST_VECTORING_SCHEMA = "rocket_propulsion_thrust_vectoring_v1"
THRUST_VECTORING_REFERENCE_IDS = (
    "NASA-SP-8059",
    "NASA-TM-20240012525",
)

Vector3 = tuple[float, float, float]


def _vector3(value: tuple[float, float, float], label: str) -> Vector3:
    """Validate and return a finite three-component vector."""
    return vector3(value, label)


def _dot(left: Vector3, right: Vector3) -> float:
    """Return the Euclidean dot product."""
    return sum(a * b for a, b in zip(left, right, strict=True))


def _cross(left: Vector3, right: Vector3) -> Vector3:
    """Return the right-handed cross product ``left x right``."""
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _add(left: Vector3, right: Vector3) -> Vector3:
    """Return the component-wise sum."""
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _subtract(left: Vector3, right: Vector3) -> Vector3:
    """Return the component-wise difference ``left - right``."""
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _scale(scale: float, vector: Vector3) -> Vector3:
    """Return the vector multiplied by a scalar."""
    return tuple(scale * component for component in vector)  # type: ignore[return-value]


def _norm(vector: Vector3) -> float:
    """Return the Euclidean norm ``sqrt(v . v)``."""
    return sqrt(_dot(vector, vector))


def _normalized(vector: Vector3, label: str) -> Vector3:
    """Return the unit vector; raise `DomainError` for a zero vector."""
    magnitude = _norm(vector)
    if magnitude <= 0.0:
        raise DomainError(f"{label} must be non-zero.")
    return _scale(1.0 / magnitude, vector)


def _rotate(vector: Vector3, axis: Vector3, angle_rad: float) -> Vector3:
    """Rotate a vector about a unit, body-fixed axis by the Rodrigues formula.

    ``R(a, t) v = cos(t) v + sin(t) (a x v) + (1 - cos(t)) (a . v) a`` for a
    unit axis ``a``; the rotation is right-handed about ``a``.
    """

    cosine = cos(angle_rad)
    sine = sin(angle_rad)
    return _add(
        _add(_scale(cosine, vector), _scale(sine, _cross(axis, vector))),
        _scale((1.0 - cosine) * _dot(axis, vector), axis),
    )


@dataclass(frozen=True, slots=True)
class BodyFrameEngineMount:
    """Immutable engine location and two-axis gimbal geometry.

    Parameters
    ----------
    engine_id : str
        Non-blank engine identifier.
    application_point_body_m : Vector3
        Thrust application point in the vehicle body frame, m.
    nominal_direction_body : Vector3
        Non-zero direction of the force applied to the vehicle (not the
        exhaust direction) at zero gimbal angles; stored normalised.
    gimbal_axis_1_body : Vector3
        Non-zero first gimbal axis in body coordinates, perpendicular to the
        nominal direction within ``1e-10``; stored normalised.
    gimbal_axis_1_limit_rad, gimbal_axis_2_limit_rad : float, default 0.0
        Symmetric angle limits in ``[0, pi/2)``, rad; zero gives a fixed
        engine about that axis.

    Raises
    ------
    DomainError
        If a field has the wrong type, a vector is zero or non-finite, axis 1
        is not perpendicular to the nominal direction, or a limit is outside
        ``[0, pi/2)``.

    Notes
    -----
    Axis 2 is the normalised cross product ``nominal_direction x axis_1``.
    The realised direction applies the axis-1 rotation first and the axis-2
    rotation second, both about body-fixed axes; this order matters for
    finite angles [2]_.  Thrust-vector pointing and misalignment design
    criteria are given in [1]_.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
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
        """Second body-fixed gimbal axis, ``unit(nominal_direction x axis_1)`` [1]_.

        References
        ----------
        .. [1] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
               Vector Control Actuator-Force / Gimbal-Torque Transformation
               Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
        """

        return _normalized(
            _cross(self.nominal_direction_body, self.gimbal_axis_1_body),
            "Derived second gimbal axis",
        )

    def direction_and_partials(
        self,
        gimbal_axis_1_rad: float,
        gimbal_axis_2_rad: float,
    ) -> tuple[Vector3, Vector3, Vector3]:
        """Return the realised thrust direction and its derivatives by both angles.

        Parameters
        ----------
        gimbal_axis_1_rad, gimbal_axis_2_rad : float
            Gimbal angles, rad, each within its limit (tolerance ``1e-14``).

        Returns
        -------
        tuple of Vector3
            ``(d, dd/dangle_1, dd/dangle_2)`` in body coordinates.

        Raises
        ------
        DomainError
            If an angle is non-finite or exceeds its limit.

        Notes
        -----
        ``d = R(a2, angle_2) R(a1, angle_1) d0`` with body-fixed axes [1]_.
        Because ``dR(a, t)v/dt = a x R(a, t)v`` for a unit axis, the
        derivatives are ``dd/dangle_1 = R(a2, angle_2) (a1 x R(a1, angle_1) d0)``
        and ``dd/dangle_2 = a2 x d``.

        References
        ----------
        .. [1] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
               Vector Control Actuator-Force / Gimbal-Torque Transformation
               Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
        """

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

    The angles drive the order-sensitive two-axis gimbal geometry of [2]_;
    thrust-vector pointing during thrusting manoeuvres is treated in [1]_.

    Parameters
    ----------
    engine_id : str
        Identifier of the commanded mount.
    thrust_n : float
        Non-negative commanded thrust magnitude, N.
    gimbal_axis_1_rad, gimbal_axis_2_rad : float, default 0.0
        Commanded gimbal angles, rad; limits are checked at evaluation.
    enabled : bool, default True
        ``False`` commands zero thrust while keeping the engine accounted for.

    Raises
    ------
    DomainError
        If a field has the wrong type, the thrust is negative or non-finite,
        or ``enabled`` is not a boolean.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
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
    """Analytic derivatives of one engine's force and moment.

    Intended for estimation and six-degree-of-freedom state-transition use.
    For a disabled engine the thrust derivatives are zero and the gimbal
    derivatives vanish with the zero thrust.  Force and moment follow the
    thrust-vector and moment-arm relations of [1]_ and the gimbal geometry of
    [2]_.

    Attributes
    ----------
    force_wrt_thrust : Vector3
        ``dF/dT``, the realised unit direction.
    moment_wrt_thrust : Vector3
        ``dM/dT = r x d``, m.
    force_wrt_gimbal_axis_1_rad, force_wrt_gimbal_axis_2_rad : Vector3
        ``dF/dangle``, N/rad.
    moment_wrt_gimbal_axis_1_rad, moment_wrt_gimbal_axis_2_rad : Vector3
        ``dM/dangle = r x dF/dangle``, N m/rad.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
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

    Attributes
    ----------
    engine_id : str
        Mount identifier.
    active : bool
        ``True`` when enabled with positive thrust.
    thrust_n : float
        Applied thrust magnitude (zero when disabled), N.
    direction_body : Vector3
        Realised unit thrust direction after both gimbal rotations [2]_.
    force_body_n : Vector3
        Force ``T d``, N.
    moment_body_n_m : Vector3
        Moment ``r x F`` about the reference point, N m [1]_, [3]_.
    moment_arm_body_m : Vector3
        Arm ``r`` from the reference point to the application point, m.
    partials : EngineWrenchPartials
        Analytic derivatives by thrust and gimbal angles.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
    .. [3] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
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

    Each tuple is keyed by engine identifier.  Gimbal entries contain the
    axis-1 and axis-2 derivatives in that order.  The derivatives are those of
    the summed thrust-vector forces and moments [1]_ under the gimbal geometry
    of [2]_.

    Attributes
    ----------
    force_wrt_thrust, moment_wrt_thrust : tuple of (str, Vector3)
        Derivatives of the resultant force and moment by each engine's thrust.
    force_wrt_gimbal_rad, moment_wrt_gimbal_rad : tuple of (str, (Vector3, Vector3))
        Derivatives of the resultant force and moment by each engine's two
        gimbal angles.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
    """

    force_wrt_thrust: tuple[tuple[str, Vector3], ...]
    moment_wrt_thrust: tuple[tuple[str, Vector3], ...]
    force_wrt_gimbal_rad: tuple[tuple[str, tuple[Vector3, Vector3]], ...]
    moment_wrt_gimbal_rad: tuple[tuple[str, tuple[Vector3, Vector3]], ...]


@dataclass(frozen=True, slots=True)
class ClusterWrenchEvaluation:
    """Body-frame force and moment about a caller-supplied reference point.

    Attributes
    ----------
    reference_point_body_m : Vector3
        Moment reference point, m.
    force_body_n : Vector3
        Resultant force, N.
    moment_body_n_m : Vector3
        Resultant moment about the reference point, N m.
    resultant_thrust_n : float
        Norm of the resultant force, N.
    resultant_direction_body : Vector3 or None
        Unit resultant direction; ``None`` when the resultant is at most
        ``1e-15`` N.
    scalar_thrust_sum_n : float
        Sum of the individual thrust magnitudes, N.
    vector_cancellation_loss_n : float
        ``scalar_thrust_sum_n - resultant_thrust_n``, N.
    engines : tuple of EngineWrench
        Per-engine contributions, sorted by identifier.
    partials : ClusterWrenchPartials
        Resultant-wrench derivatives by every command.
    model_sha256 : str
        Hash of the evaluated cluster geometry.

    Notes
    -----
    Per-engine forces and moments follow the thrust-vector and moment-arm
    relations of [1]_ and the gimbal geometry of [2]_, and are summed.  The
    vector cancellation loss is not a propulsive efficiency or plume loss; it
    is solely the geometric loss caused by summing non-collinear force
    vectors.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
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

    Parameters
    ----------
    mounts : tuple of BodyFrameEngineMount
        At least one mount with a unique identifier; stored sorted by
        identifier.
    model_id : str, default "mounted-engine-cluster"
        Model identifier recorded in the manifest.
    schema : str, default THRUST_VECTORING_SCHEMA
        Document schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If a field has the wrong type, the schema is unsupported, no mount is
        given, or identifiers repeat.

    Notes
    -----
    The evaluation uses ``moment = (application_point - reference_point) x
    force`` [1]_, [3]_.  Passing the instantaneous body-frame centre of mass
    therefore yields the propulsion torque required by a host 6-DOF
    propagator.  The cluster neither computes that centre of mass nor
    integrates attitude.  Gimbal rotations follow the order-sensitive
    convention of [2]_.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
           Vector Control Actuator-Force / Gimbal-Torque Transformation
           Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
    .. [3] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
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
        """Return the hash-covered geometry, convention and reference content."""
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
        """Canonical SHA-256 identity of the geometry, convention and references."""

        return canonical_sha256(self._payload())

    def model_manifest(self) -> dict[str, Any]:
        """Return the integration-facing provenance and semantics manifest.

        Returns
        -------
        dict
            Schema, model identity and hash, reference identifiers, rotation
            convention, force and moment frames, moment reference and the
            statement that attitude propagation is owned by the host.
        """

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
        """Return the geometry document with its recomputable semantic hash."""

        payload = self._payload()
        payload["model_sha256"] = self.model_sha256
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the geometry as deterministic, key-sorted, finite JSON.

        Parameters
        ----------
        indent : int or None, default 2
            JSON indentation; ``None`` gives the compact form.

        Returns
        -------
        str
            JSON text without NaN or infinity tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return json.dumps(
            self.to_dict(),
            allow_nan=False,
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> MountedEngineCluster:
        """Reopen serialized geometry and reject schema, convention or hash changes.

        Parameters
        ----------
        payload : Mapping
            Document produced by `to_dict`.

        Returns
        -------
        MountedEngineCluster
            The reconstructed, validated cluster.

        Raises
        ------
        InputError
            If the root is not an object, the schema, reference identifiers or
            rotation convention differ from the current schema, a field is
            malformed, or ``model_sha256`` does not match the contents.
        """

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
        """Parse and integrity-check a JSON cluster document.

        Parameters
        ----------
        document : str
            JSON text produced by `to_json`.

        Returns
        -------
        MountedEngineCluster
            The reconstructed, validated cluster.

        Raises
        ------
        InputError
            If the text is not a conforming JSON document [1]_ or fails the
            checks of `from_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return _read_document(document, "Mounted engine cluster", cls.from_dict)

    def evaluate(
        self,
        commands: tuple[EngineVectorCommand, ...],
        *,
        reference_point_body_m: Vector3 = (0.0, 0.0, 0.0),
    ) -> ClusterWrenchEvaluation:
        """Compose all engine forces, moments and analytic command partials.

        Parameters
        ----------
        commands : tuple of EngineVectorCommand
            Exactly one command per mount, including disabled or zero-thrust
            engines.
        reference_point_body_m : Vector3, default (0, 0, 0)
            Moment reference point in the body frame, normally the
            instantaneous centre of mass, m.

        Returns
        -------
        ClusterWrenchEvaluation
            Resultant force and moment, per-engine contributions, analytic
            partials and the geometry hash.

        Raises
        ------
        DomainError
            If an argument has the wrong type, commands repeat, are missing or
            name unknown mounts, a gimbal angle exceeds its limit, or the
            resultant exceeds the scalar thrust sum beyond
            ``max(1e-12, 1e-12 * sum)``.

        Notes
        -----
        Every configured engine requires exactly one command, so an omitted
        engine is never mistaken for an intentional shutdown.  Each engine
        contributes ``F_i = T_i d_i`` and ``M_i = r_i x F_i`` with
        ``r_i = p_i - reference_point`` [1]_, [3]_; the realised direction
        ``d_i`` follows the two-axis gimbal convention of [2]_.  Force and
        moment derivatives are ``dF/dT = d``, ``dM/dT = r x d``,
        ``dF/dangle = T dd/dangle`` and ``dM/dangle = r x dF/dangle``.

        References
        ----------
        .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
               Maneuvers", 1971. https://ntrs.nasa.gov/citations/19710016722
        .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust
               Vector Control Actuator-Force / Gimbal-Torque Transformation
               Matrix", 2024. https://ntrs.nasa.gov/citations/20240012525
        .. [3] B. Wie, "Space Vehicle Dynamics and Control", 2nd ed., AIAA, 2008.
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
