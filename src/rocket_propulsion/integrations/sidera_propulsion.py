"""Primary propulsion-provider hand-off for Sidera numerical propagation.

This module is the stable boundary between Rocket Propulsion Lab and Sidera.
It deliberately contains no orbit integrator: Sidera owns epochs, frames,
attitude, event localization, and state integration; this package owns thrust,
propellant drains, operating-domain limits, propulsion events, and their
analytic partial derivatives. The split mirrors the force, mass-derivative,
parameter and event contract of a propagator propulsion model [1]_.

References
----------
.. [1] Orekit 13.1.5, "PropulsionModel" API,
       https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
.. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket
       Engine Systems", https://ntrs.nasa.gov/citations/20040000363
.. [3] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
       502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from math import isclose

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.core.validation import (
    real,
    records,
    strings,
    text,
    vector3,
)
from rocket_propulsion.propulsion.burns import (
    ClusterWrenchEvaluation,
    EngineVectorCommand,
    MountedEngineCluster,
    PropagatorPropulsionBridge,
    PropulsionDynamicsEvaluation,
)
from rocket_propulsion.propulsion.burns.exchange_documents import _read_document
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256

from .contracts import SideraCapabilities

SIDERA_PROPULSION_HANDOFF_SCHEMA = "rocket_propulsion_sidera_provider_handoff_v1"
SIDERA_PROPULSION_HANDOFF_REFERENCE_IDS = (
    "OREKIT-13.1.5-PROPULSION-MODEL",
    "NASA-20040000363",
    "NASA-SP-8059",
    "NASA-TM-20240012525",
    "CCSDS-502.0-B-3",
)
SIDERA_PROPULSION_REQUIRED_CAPABILITIES = (
    "propulsion_provider_available",
    "independent_mass_flow_available",
    "additional_states_available",
    "event_detectors_available",
    "analytic_partials_available",
)


@dataclass(frozen=True, slots=True)
class SideraBodyPropulsionEvaluation:
    """Coupled propulsion state rates and body-frame force and moment.

    Attributes
    ----------
    dynamics : PropulsionDynamicsEvaluation
        Scalar-thrust bridge evaluation probed along body +X.
    wrench : ClusterWrenchEvaluation
        Force and moment of the mounted engine cluster about the reference
        point, in the body frame.
    force_closure_residual_n : float
        Sum of enabled command thrusts minus the bridge's axial thrust, N.

    Notes
    -----
    ``dynamics.acceleration_m_s2`` is only the bridge's scalar-thrust probe
    along body +X. A 6-DOF Sidera consumer must use
    ``wrench.force_body_n / vehicle_mass`` after transforming that force with
    its attitude state. Keeping both records exposes, rather than hides, the
    closure between scalar performance and mounted hardware: thrust pointing
    and moment arms [1]_ and the gimbal transformation [2]_.

    References
    ----------
    .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
           Maneuvers", https://ntrs.nasa.gov/citations/19710016722
    .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust Vector
           Control Actuator-Force / Gimbal-Torque Transformation Matrix",
           https://ntrs.nasa.gov/citations/20240012525
    """

    dynamics: PropulsionDynamicsEvaluation
    wrench: ClusterWrenchEvaluation
    force_closure_residual_n: float


@dataclass(frozen=True, slots=True)
class SideraPropulsionProviderHandoff:
    """Hash-bound, propagation-ready Rocket Propulsion Lab provider contract.

    Attributes
    ----------
    bridge : PropagatorPropulsionBridge
        Canonical propulsion bridge; re-opened from its serialization on
        construction.
    ignition_epoch_s : float
        Sidera epoch at which the bridge's relative burn time is zero, s.
    additional_state_order : tuple of str
        Order of the propellant states in Sidera's additional-state vector;
        empty means the sorted bound state names.
    engine_cluster : MountedEngineCluster or None
        Optional mounted cluster for the 6-DOF body-wrench contract.
    provider_id : str
        Provider identifier.
    schema : str
        Must equal ``SIDERA_PROPULSION_HANDOFF_SCHEMA``.

    Raises
    ------
    DomainError
        If the schema differs, the bridge or cluster has the wrong type,
        the epoch or identifiers are malformed, or the state order repeats
        names or does not cover every bound propellant state exactly once.

    Notes
    -----
    ``ignition_epoch_s`` maps Sidera's absolute propagation time to the
    bridge's relative burn time. The optional engine cluster adds a body
    wrench without making propulsion own attitude propagation. The contract
    follows the force, mass-derivative, parameter-driver and event split of
    [1]_, with transient-modelling context from [2]_ and versioned,
    unit-explicit exchange [3]_. It does not require an installed Sidera
    provider: it can be serialized, inspected, replayed and evaluated
    independently, which supports the traceability expected by [4]_.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel" API,
           https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket
           Engine Systems", https://ntrs.nasa.gov/citations/20040000363
    .. [3] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
           502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
    .. [4] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    bridge: PropagatorPropulsionBridge
    ignition_epoch_s: float
    additional_state_order: tuple[str, ...] = ()
    engine_cluster: MountedEngineCluster | None = None
    provider_id: str = "rocket-propulsion-lab"
    schema: str = SIDERA_PROPULSION_HANDOFF_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != SIDERA_PROPULSION_HANDOFF_SCHEMA:
            raise DomainError("Unsupported Sidera propulsion hand-off schema.")
        if not isinstance(self.bridge, PropagatorPropulsionBridge):
            raise DomainError("Sidera propulsion hand-off requires a propagator bridge.")
        object.__setattr__(
            self,
            "ignition_epoch_s",
            real(self.ignition_epoch_s, "Sidera propulsion ignition epoch"),
        )
        text(self.provider_id, "Sidera propulsion provider identifier")
        text(self.schema, "Sidera propulsion hand-off schema")
        if self.engine_cluster is not None and not isinstance(
            self.engine_cluster, MountedEngineCluster
        ):
            raise DomainError("Sidera engine cluster must be a MountedEngineCluster.")

        bound_states = tuple(
            sorted({binding.state_name for binding in self.bridge.stream_bindings})
        )
        supplied_order = strings(
            self.additional_state_order, "Sidera propulsion additional-state order"
        )
        if not supplied_order:
            supplied_order = bound_states
        if len(supplied_order) != len(set(supplied_order)):
            raise DomainError("Sidera additional-state names must be unique.")
        if set(supplied_order) != set(bound_states):
            raise DomainError(
                "Sidera additional-state order must cover every bound propulsion state "
                "exactly once."
            )
        object.__setattr__(self, "additional_state_order", supplied_order)
        # Re-open nested public artifacts at the trust boundary.
        PropagatorPropulsionBridge.from_dict(self.bridge.to_dict())
        if self.engine_cluster is not None:
            MountedEngineCluster.from_dict(self.engine_cluster.to_dict())

    @property
    def required_operating_conditions(self) -> tuple[str, ...]:
        """tuple of str: Operating-condition axis names in evaluation order.

        Empty when the bridge has no performance surface.
        """

        surface = self.bridge.performance_surface or self.bridge.stream_performance_surface
        return () if surface is None else tuple(axis.name for axis in surface.axes)

    @property
    def handoff_hash(self) -> str:
        """str: Canonical SHA-256 of the manifest without its own hash."""

        return canonical_sha256(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the state, event, frame and capability semantics of the hand-off.

        Parameters
        ----------
        include_hash : bool, optional
            Add ``handoff_hash`` (the canonical SHA-256 of the other fields).

        Returns
        -------
        dict
            Schema, provider id, ignition epoch, nested model hashes, state
            layout, operating conditions, estimable parameters, event ids,
            force/mass/time semantics, required Sidera capabilities and
            reference identifiers.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel" API,
               https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "provider_id": self.provider_id,
            "ignition_epoch_s": self.ignition_epoch_s,
            "bridge_model_hash": self.bridge.model_hash,
            "engine_cluster_hash": (
                None if self.engine_cluster is None else self.engine_cluster.model_sha256
            ),
            "state_layout": {
                "vehicle_mass": "host_state",
                "additional_states": list(self.additional_state_order),
            },
            "required_operating_conditions": list(self.required_operating_conditions),
            "estimable_parameters": [item.name for item in self.bridge.parameters],
            "event_ids": [item.event_id for item in self.bridge.event_surfaces()],
            "operating_domain_event_ids": [
                item.event_id for item in self.bridge.operating_condition_surfaces()
            ],
            "interpolation_grid_event_ids": [
                item.event_id for item in self.bridge.operating_condition_grid_surfaces()
            ],
            "force_frame": (
                "host_resolved_direction" if self.engine_cluster is None else "vehicle_body"
            ),
            "mass_semantics": "simultaneous_vehicle_and_named_stream_derivatives",
            "time_semantics": "absolute_sidera_epoch_minus_ignition_epoch",
            "attitude_owned_by": "sidera",
            "required_sidera_capabilities": [
                "propulsion_provider",
                "independent_mass_flow",
                "additional_states",
                "event_detectors",
                "analytic_partials",
            ]
            + (["body_frame_force_moment"] if self.engine_cluster is not None else []),
            "reference_ids": list(SIDERA_PROPULSION_HANDOFF_REFERENCE_IDS),
        }
        if include_hash:
            payload["handoff_hash"] = canonical_sha256(payload)
        return payload

    def to_dict(self) -> dict[str, object]:
        """Return the manifest together with the serialized bridge and engine cluster."""

        payload = self.manifest()
        payload["bridge"] = self.bridge.to_dict()
        payload["engine_cluster"] = (
            None if self.engine_cluster is None else self.engine_cluster.to_dict()
        )
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the hand-off as deterministic finite JSON.

        Parameters
        ----------
        indent : int or None, optional
            Indentation passed to :func:`json.dumps`.

        Returns
        -------
        str
            JSON text with sorted keys and no non-finite numbers [1]_.

        Raises
        ------
        ValueError
            If a value is not finite (``allow_nan=False``).

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format", https://www.rfc-editor.org/rfc/rfc8259
        """

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: object) -> SideraPropulsionProviderHandoff:
        """Rebuild the provider from a decoded document and verify its hash.

        Parameters
        ----------
        payload : object
            Decoded hand-off document; it must be a mapping.

        Returns
        -------
        SideraPropulsionProviderHandoff
            The verified hand-off.

        Raises
        ------
        InputError
            If the root is not a mapping, a field is missing or malformed, or
            ``handoff_hash`` does not match the rebuilt contents.
        """

        if not isinstance(payload, Mapping):
            raise InputError("Sidera propulsion hand-off root must be an object.")
        try:
            cluster_payload = payload["engine_cluster"]
            state_layout = payload["state_layout"]
            handoff = cls(
                bridge=PropagatorPropulsionBridge.from_dict(payload["bridge"]),
                ignition_epoch_s=payload["ignition_epoch_s"],
                additional_state_order=state_layout["additional_states"],
                engine_cluster=(
                    None
                    if cluster_payload is None
                    else MountedEngineCluster.from_dict(cluster_payload)
                ),
                provider_id=payload["provider_id"],
                schema=payload["schema"],
            )
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed Sidera propulsion hand-off document.") from error
        if payload.get("handoff_hash") != handoff.handoff_hash:
            raise InputError("Sidera propulsion hand-off hash does not match its contents.")
        return handoff

    @classmethod
    def from_json(cls, document: str) -> SideraPropulsionProviderHandoff:
        """Parse JSON text and integrity-check the provider hand-off.

        Parameters
        ----------
        document : str
            JSON document text.

        Returns
        -------
        SideraPropulsionProviderHandoff
            The verified hand-off.

        Raises
        ------
        InputError
            If the text violates the strict JSON reading rules of
            :func:`~rocket_propulsion.core.validation.read_json_document`
            [1]_, [2]_ or fails :meth:`from_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format", https://www.rfc-editor.org/rfc/rfc8259
        .. [2] PYTHON-3.12-JSON. Python 3.12, "json - JSON encoder and decoder",
               https://docs.python.org/3.12/library/json.html
        """

        return _read_document(document, "Sidera propulsion hand-off", cls.from_dict)

    def evaluate(
        self,
        *,
        epoch_s: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float] | None = None,
        operating_conditions: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        ambient_pressure_pa: float | None = None,
        side: str = "right",
        apply_inventory_limits: bool = True,
        boundary_tolerance_s: float | None = None,
    ) -> PropulsionDynamicsEvaluation:
        """Evaluate propulsion at a Sidera epoch through the canonical bridge.

        Parameters
        ----------
        epoch_s : float
            Absolute Sidera epoch, s; the bridge is evaluated at
            ``epoch_s - ignition_epoch_s``.
        vehicle_mass_kg : float
            Current vehicle mass, kg.
        direction : tuple of float
            Host-resolved thrust direction.
        additional_states : mapping, optional
            Propellant state values by name.
        operating_conditions : mapping, optional
            Operating-condition values by axis name.
        parameter_overrides : mapping, optional
            Estimable-parameter overrides.
        ambient_pressure_pa : float, optional
            Ambient pressure for pressure-thrust correction, Pa.
        side : {"right", "left"}, optional
            One-sided limit at a discontinuity.
        apply_inventory_limits : bool, optional
            Apply the bridge's inventory clamps.
        boundary_tolerance_s : float, optional
            Tolerance for snapping to a profile boundary, s.

        Returns
        -------
        PropulsionDynamicsEvaluation
            Thrust, acceleration, mass derivative and named state
            derivatives.

        Raises
        ------
        DomainError
            If the epoch or direction is malformed or the bridge refuses the
            evaluation.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel" API,
               https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
        """

        epoch = real(epoch_s, "Sidera propagation epoch")
        return self.bridge.evaluate(
            relative_time_s=epoch - self.ignition_epoch_s,
            vehicle_mass_kg=vehicle_mass_kg,
            direction=vector3(direction, "Sidera resolved thrust direction"),
            additional_states=additional_states,
            operating_conditions=operating_conditions,
            parameter_overrides=parameter_overrides,
            ambient_pressure_pa=ambient_pressure_pa,
            side=side,
            apply_inventory_limits=apply_inventory_limits,
            boundary_tolerance_s=boundary_tolerance_s,
        )

    def ordered_additional_state_derivatives(
        self, evaluation: PropulsionDynamicsEvaluation
    ) -> tuple[float, ...]:
        """Pack named state derivatives in the declared Sidera state-vector order.

        Parameters
        ----------
        evaluation : PropulsionDynamicsEvaluation
            Result of :meth:`evaluate`.

        Returns
        -------
        tuple of float
            Derivatives ordered as ``additional_state_order``.

        Raises
        ------
        DomainError
            If ``evaluation`` has the wrong type or its state names differ
            from the declared layout.
        """

        if not isinstance(evaluation, PropulsionDynamicsEvaluation):
            raise DomainError("Expected a PropulsionDynamicsEvaluation.")
        derivatives = dict(evaluation.additional_state_derivatives)
        if set(derivatives) != set(self.additional_state_order):
            raise DomainError("Propulsion evaluation does not match the Sidera state layout.")
        return tuple(derivatives[name] for name in self.additional_state_order)

    def evaluate_body_frame(
        self,
        *,
        epoch_s: float,
        vehicle_mass_kg: float,
        commands: tuple[EngineVectorCommand, ...],
        reference_point_body_m: tuple[float, float, float],
        additional_states: Mapping[str, float] | None = None,
        operating_conditions: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        ambient_pressure_pa: float | None = None,
        side: str = "right",
        apply_inventory_limits: bool = True,
        boundary_tolerance_s: float | None = None,
        closure_tolerance_n: float = 1.0e-9,
    ) -> SideraBodyPropulsionEvaluation:
        """Evaluate scalar propulsion and close it to a mounted body wrench.

        Parameters
        ----------
        epoch_s : float
            Absolute Sidera epoch, s.
        vehicle_mass_kg : float
            Current vehicle mass, kg.
        commands : tuple of EngineVectorCommand
            Per-engine thrust and gimbal commands.
        reference_point_body_m : tuple of float
            Moment reference point in the body frame, m.
        additional_states, operating_conditions, parameter_overrides : mapping, optional
            As for :meth:`evaluate`.
        ambient_pressure_pa : float, optional
            Ambient pressure, Pa.
        side : {"right", "left"}, optional
            One-sided limit at a discontinuity.
        apply_inventory_limits : bool, optional
            Apply the bridge's inventory clamps.
        boundary_tolerance_s : float, optional
            Boundary snapping tolerance, s.
        closure_tolerance_n : float, optional
            Absolute force-closure tolerance, N; nonnegative. A relative
            tolerance of 1e-12 also applies.

        Returns
        -------
        SideraBodyPropulsionEvaluation
            Scalar dynamics, body wrench and force-closure residual.

        Raises
        ------
        DomainError
            If the hand-off has no engine cluster, the commands or tolerance
            are invalid, or the enabled command thrusts do not close to the
            bridge's axial thrust.

        Notes
        -----
        The sum of enabled engine-command thrust magnitudes must reproduce
        the bridge's instantaneous scalar axial thrust. This prevents
        guidance or an adapter from silently creating or losing thrust while
        distributing it across a gimballed cluster [1]_, [2]_.

        References
        ----------
        .. [1] NASA SP-8059. NASA, "Spacecraft Attitude Control During Thrusting
               Maneuvers", https://ntrs.nasa.gov/citations/19710016722
        .. [2] NASA/TM-20240012525. NASA, "MSFC-TVC TB-03: Derivation of Thrust Vector
               Control Actuator-Force / Gimbal-Torque Transformation Matrix",
               https://ntrs.nasa.gov/citations/20240012525
        """

        if self.engine_cluster is None:
            raise DomainError("This Sidera hand-off has no mounted engine cluster.")
        commands = records(commands, EngineVectorCommand, "Sidera engine commands")
        tolerance = real(closure_tolerance_n, "Sidera force-closure tolerance")
        if tolerance < 0.0:
            raise DomainError("Sidera force-closure tolerance must be non-negative.")
        dynamics = self.evaluate(
            epoch_s=epoch_s,
            vehicle_mass_kg=vehicle_mass_kg,
            direction=(1.0, 0.0, 0.0),
            additional_states=additional_states,
            operating_conditions=operating_conditions,
            parameter_overrides=parameter_overrides,
            ambient_pressure_pa=ambient_pressure_pa,
            side=side,
            apply_inventory_limits=apply_inventory_limits,
            boundary_tolerance_s=boundary_tolerance_s,
        )
        commanded_scalar_thrust = sum(command.thrust_n for command in commands if command.enabled)
        residual = commanded_scalar_thrust - dynamics.axial_thrust_n
        if not isclose(
            commanded_scalar_thrust,
            dynamics.axial_thrust_n,
            rel_tol=1.0e-12,
            abs_tol=tolerance,
        ):
            raise DomainError(
                "Mounted-engine command thrust does not close to bridge axial thrust."
            )
        wrench = self.engine_cluster.evaluate(
            commands,
            reference_point_body_m=vector3(
                reference_point_body_m, "Sidera body-wrench reference point"
            ),
        )
        return SideraBodyPropulsionEvaluation(dynamics, wrench, residual)


def missing_sidera_propulsion_capabilities(
    capabilities: SideraCapabilities,
    handoff: SideraPropulsionProviderHandoff,
) -> tuple[str, ...]:
    """Return the public Sidera capabilities that ``handoff`` needs but lacks.

    Parameters
    ----------
    capabilities : SideraCapabilities
        Observed capabilities.
    handoff : SideraPropulsionProviderHandoff
        Provider contract; an engine cluster additionally requires
        ``body_frame_force_moment_available``.

    Returns
    -------
    tuple of str
        Names of missing capability flags; empty when none is missing.

    Raises
    ------
    DomainError
        If either argument has the wrong type.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel" API,
           https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    """

    if not isinstance(capabilities, SideraCapabilities):
        raise DomainError("Expected SideraCapabilities.")
    if not isinstance(handoff, SideraPropulsionProviderHandoff):
        raise DomainError("Expected a SideraPropulsionProviderHandoff.")
    required = list(SIDERA_PROPULSION_REQUIRED_CAPABILITIES)
    if handoff.engine_cluster is not None:
        required.append("body_frame_force_moment_available")
    return tuple(name for name in required if not getattr(capabilities, name))


def require_sidera_propulsion_capabilities(
    capabilities: SideraCapabilities,
    handoff: SideraPropulsionProviderHandoff,
) -> None:
    """Fail closed unless Sidera can consume the provider hand-off exactly.

    Parameters
    ----------
    capabilities : SideraCapabilities
        Observed capabilities.
    handoff : SideraPropulsionProviderHandoff
        Provider contract.

    Raises
    ------
    FeatureUnavailableError
        If a required capability is missing; ``details`` lists the missing
        capabilities, the hand-off hash and the Sidera version.
    DomainError
        If either argument has the wrong type.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel" API,
           https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    missing = missing_sidera_propulsion_capabilities(capabilities, handoff)
    if missing:
        raise FeatureUnavailableError(
            "Sidera cannot yet consume the propulsion-provider hand-off exactly.",
            details={
                "missing_capabilities": list(missing),
                "handoff_hash": handoff.handoff_hash,
                "sidera_version": capabilities.version,
            },
        )
