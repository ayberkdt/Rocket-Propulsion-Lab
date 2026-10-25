"""Weighted tabulated-burn hand-off for a Sidera ensemble consumer.

This module binds propulsion artifacts to execution time, direction and frame
without importing or modifying Sidera. It is an adapter contract, not an orbit
propagator. Net axial thrust is the force channel and total tank flow is the
mass channel, so canted clusters and open-cycle engines remain physically
closed. Weighted realizations carry scenario, conditional and joint
probabilities in the sense of a probabilistic ensemble [3]_; the document is
versioned and unit-explicit [1]_ and hash-bound for traceability [2]_.

References
----------
.. [1] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
       502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
.. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
.. [3] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures Guide
       for NASA Managers and Practitioners",
       https://ntrs.nasa.gov/citations/20120001369
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from math import isclose, isfinite
from typing import Any, Literal

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.core.validation import (
    integer,
    read_json_document,
    real,
    records,
    strings,
    text,
)
from rocket_propulsion.propulsion.burns import HierarchicalPropulsionEnsemble
from rocket_propulsion.propulsion.burns.ensemble_documents import _handoff_document
from rocket_propulsion.propulsion.burns.ensemble_serialization import (
    HIERARCHICAL_ENSEMBLE_SCHEMA,
    hierarchical_ensemble_from_dict,
    hierarchical_ensemble_hash,
    hierarchical_ensemble_to_dict,
)
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational, _fraction
from rocket_propulsion.propulsion.burns.serialization import canonical_json_bytes
from rocket_propulsion.propulsion.burns.tabulated import TabulatedBurnArtifact

from .contracts import SideraCapabilities
from .direction import _ADAPTER_DIRECTION_POLICY, _directed_vector, _stored_unit_direction

SIDERA_ENSEMBLE_HANDOFF_SCHEMA = "rocket_propulsion_sidera_ensemble_handoff_v1"
SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES = (
    "tabulated_finite_burn",
    "independent_mass_flow",
    "segment_boundary_splitting",
    "weighted_maneuver_ensemble",
)
_FRAMES = {"inertial", "ric", "vnb"}


def _sha256_hex(value: str, field: str) -> str:
    """Return a lower-case 64-digit hexadecimal SHA-256 digest or raise :class:`DomainError`."""
    normalized = text(value, field).lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise DomainError(f"{field} must be a SHA-256 hex digest.")
    return normalized


def _unit_direction(values: tuple[float, float, float]) -> tuple[float, float, float]:
    """Return the unit vector of ``values`` rounded once per component."""
    return _directed_vector(values, 1., label="Sidera ensemble unit direction")


@dataclass(frozen=True, slots=True)
class SideraEnsembleBurnJob:
    """One weighted Sidera propagation request backed by a full burn artifact.

    Attributes
    ----------
    scenario_id : str
        Scenario identifier.
    conditional_index : int
        Nonnegative index of the realization within its scenario.
    scenario_probability, conditional_probability, joint_probability : float
        Probabilities in [0, 1]; ``joint = scenario * conditional`` within an
        absolute and relative tolerance of 1e-12.
    t_start_s : float
        Absolute solver time at which the propulsion history starts, s.
    direction : tuple of float
        Unit thrust direction in ``frame``; normalized on construction.
    frame : {"inertial", "ric", "vnb"}
        Direction frame.
    artifact : TabulatedBurnArtifact
        Complete tabulated burn history.

    Raises
    ------
    DomainError
        If a field is malformed, a probability is outside [0, 1], the joint
        probability does not factor, or the frame is unsupported.

    Notes
    -----
    The probability decomposition follows the conditional-probability
    treatment of [2]_; time and frame metadata follow [1]_.

    References
    ----------
    .. [1] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
           502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
    .. [2] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures Guide
           for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369
    """

    scenario_id: str
    conditional_index: int
    scenario_probability: float
    conditional_probability: float
    joint_probability: float
    t_start_s: float
    direction: tuple[float, float, float]
    frame: Literal["inertial", "ric", "vnb"]
    artifact: TabulatedBurnArtifact

    def __post_init__(self) -> None:
        object.__setattr__(self, "scenario_id", text(self.scenario_id, "Sidera ensemble scenario ID"))
        object.__setattr__(
            self, "conditional_index", integer(self.conditional_index, "Conditional index")
        )
        for name in ("scenario_probability", "conditional_probability", "joint_probability",
                     "t_start_s"):
            object.__setattr__(self, name, real(getattr(self, name), f"Sidera job {name}"))
        if not isinstance(self.artifact, TabulatedBurnArtifact):
            raise DomainError("Sidera ensemble jobs require a TabulatedBurnArtifact.")
        if self.conditional_index < 0:
            raise DomainError("Sidera conditional index must be nonnegative.")
        probabilities = (
            self.scenario_probability,
            self.conditional_probability,
            self.joint_probability,
        )
        if any(
            not isfinite(value) or not 0.0 <= value <= 1.0
            for value in probabilities
        ):
            raise DomainError("Sidera ensemble probabilities must lie in [0, 1].")
        if not isclose(
            self.joint_probability,
            self.scenario_probability * self.conditional_probability,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise DomainError("Sidera joint probability does not factor correctly.")
        if not isfinite(self.t_start_s):
            raise DomainError("Sidera ensemble start time must be finite.")
        frame = text(self.frame, "Sidera ensemble frame").strip().lower()
        if frame not in _FRAMES:
            raise DomainError("Sidera ensemble frame must be inertial, ric, or vnb.")
        object.__setattr__(self, "frame", frame)
        object.__setattr__(self, "direction", _unit_direction(self.direction))

    @property
    def t_end_s(self) -> float:
        """float: Absolute solver time at the end of the propulsion history, s.

        ``t_start_s + artifact.duration_s``, summed exactly and rounded once.

        References
        ----------
        .. [1] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
               502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        return _finite_rational(_fraction(self.t_start_s) + _fraction(self.artifact.duration_s), "Sidera job end time")


@dataclass(frozen=True, slots=True)
class SideraEnsembleHandoff:
    """Integrity-checked weighted transient jobs for a Sidera ensemble consumer.

    Attributes
    ----------
    source_ensemble_hash : str
        SHA-256 of the source hierarchical ensemble (lower-case hex).
    jobs : tuple of SideraEnsembleBurnJob
        At least one job; ``(scenario_id, conditional_index)`` pairs are
        unique.
    joint_probability_sum : float
        Declared sum of joint probabilities; must reproduce and equal one
        within 1e-12.
    warnings : tuple of str
        Non-blank consumer obligations.
    required_capabilities : tuple of str
        Must equal ``SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES``.
    schema : str
        Must equal ``SIDERA_ENSEMBLE_HANDOFF_SCHEMA``.

    Raises
    ------
    DomainError
        If the schema or capability contract differs, the hash is malformed,
        identities repeat, conditional indices are incomplete, conditional
        probabilities of a scenario do not sum to one, a scenario
        probability varies within its scenario, or the probability sums do
        not close.

    Notes
    -----
    Probability closure checks follow the ensemble treatment of [2]_; the
    hand-off is hash-bound for traceability of model inputs and outputs [1]_.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures Guide
           for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369
    """

    source_ensemble_hash: str
    jobs: tuple[SideraEnsembleBurnJob, ...]
    joint_probability_sum: float
    warnings: tuple[str, ...]
    required_capabilities: tuple[str, ...] = SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES
    schema: str = SIDERA_ENSEMBLE_HANDOFF_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(self, "schema", text(self.schema, "Sidera handoff schema"))
        object.__setattr__(self, "jobs", records(self.jobs, SideraEnsembleBurnJob, "Sidera jobs"))
        object.__setattr__(self, "warnings", strings(self.warnings, "Sidera warnings"))
        object.__setattr__(self, "required_capabilities", strings(self.required_capabilities, "Sidera capabilities"))
        object.__setattr__(self, "joint_probability_sum", real(self.joint_probability_sum, "Sidera probability sum"))
        if self.schema != SIDERA_ENSEMBLE_HANDOFF_SCHEMA:
            raise DomainError("Unsupported Sidera ensemble hand-off schema.")
        object.__setattr__(
            self,
            "source_ensemble_hash",
            _sha256_hex(self.source_ensemble_hash, "Source ensemble hash"),
        )
        if not self.jobs:
            raise DomainError("Sidera ensemble hand-off requires at least one job.")
        identities = tuple(
            (item.scenario_id, item.conditional_index) for item in self.jobs
        )
        if len(identities) != len(set(identities)):
            raise DomainError("Sidera ensemble job identities must be unique.")
        computed_sum = sum(item.joint_probability for item in self.jobs)
        if not isclose(computed_sum, 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise DomainError("Sidera ensemble joint probabilities must close to one.")
        if not isclose(
            self.joint_probability_sum,
            computed_sum,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise DomainError("Declared Sidera ensemble probability sum does not reproduce.")
        for scenario_id in dict.fromkeys(item.scenario_id for item in self.jobs):
            selected = tuple(item for item in self.jobs if item.scenario_id == scenario_id)
            if tuple(sorted(item.conditional_index for item in selected)) != tuple(range(len(selected))):
                raise DomainError("Sidera conditional indices must be complete.", field="conditional_index")
            if not isclose(
                sum(item.conditional_probability for item in selected),
                1.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise DomainError("Sidera conditional probabilities must close to one.")
            if len({item.scenario_probability for item in selected}) != 1:
                raise DomainError("Sidera scenario probability must be conditionally constant.")
        if self.required_capabilities != SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES:
            raise DomainError("Sidera ensemble capability contract does not match schema v1.")
        if any(not warning.strip() for warning in self.warnings):
            raise DomainError("Sidera ensemble warnings must not be empty.")

    @property
    def handoff_hash(self) -> str:
        """str: SHA-256 of the canonical JSON document without ``handoff_hash``.

        The digest covers all jobs, weights, semantics and nested artifacts.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return hashlib.sha256(
            canonical_json_bytes(self.to_dict(include_hash=False))
        ).hexdigest()

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the canonical hand-off document.

        Parameters
        ----------
        include_hash : bool, optional
            Add ``handoff_hash`` to the document.

        Returns
        -------
        dict
            Schema, source hash, force/mass/direction semantics, required
            capabilities, probability sum, jobs (each with its end time,
            artifact hash and full artifact) and warnings.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        payload: dict[str, Any] = {
            "schema": self.schema,
            "source_ensemble_schema": HIERARCHICAL_ENSEMBLE_SCHEMA,
            "source_ensemble_hash": self.source_ensemble_hash,
            "force_channel": "axial_thrust_n",
            "mass_flow_channel": "total_tank_mass_flow_kg_s",
            "direction_semantics": "unit_vector_in_declared_frame",
            "direction_numerics": _ADAPTER_DIRECTION_POLICY,
            "required_capabilities": list(self.required_capabilities),
            "joint_probability_sum": self.joint_probability_sum,
            "jobs": [
                {
                    "scenario_id": item.scenario_id,
                    "conditional_index": item.conditional_index,
                    "scenario_probability": item.scenario_probability,
                    "conditional_probability": item.conditional_probability,
                    "joint_probability": item.joint_probability,
                    "t_start_s": item.t_start_s,
                    "t_end_s": item.t_end_s,
                    "direction": list(item.direction),
                    "frame": item.frame,
                    "artifact_hash": item.artifact.artifact_hash,
                    "artifact": item.artifact.to_dict(),
                }
                for item in self.jobs
            ],
            "warnings": list(self.warnings),
        }
        if include_hash:
            payload["handoff_hash"] = self.handoff_hash
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize this hand-off as deterministic JSON.

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
        .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            indent=indent,
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> SideraEnsembleHandoff:
        """Rebuild a hand-off from a decoded document and verify every hash.

        Parameters
        ----------
        payload : dict
            Decoded hand-off document.

        Returns
        -------
        SideraEnsembleHandoff
            The verified hand-off.

        Raises
        ------
        InputError
            If the document shape, schema or physical semantics differ, a
            value is malformed, or a nested artifact hash, job end time or
            the outer hand-off hash does not reproduce.

        Notes
        -----
        Stored directions must already be numerical unit vectors and are
        kept bit-for-bit, because they are covered by the hash.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        _handoff_document(payload)
        if payload.get("schema") != SIDERA_ENSEMBLE_HANDOFF_SCHEMA:
            raise InputError("Unsupported Sidera ensemble hand-off schema.")
        expected_semantics = {
            "source_ensemble_schema": HIERARCHICAL_ENSEMBLE_SCHEMA,
            "force_channel": "axial_thrust_n",
            "mass_flow_channel": "total_tank_mass_flow_kg_s",
            "direction_semantics": "unit_vector_in_declared_frame",
            "direction_numerics": _ADAPTER_DIRECTION_POLICY,
        }
        if any(payload.get(key) != value for key, value in expected_semantics.items()):
            raise InputError("Sidera ensemble hand-off physical semantics changed.")
        try:
            jobs = []
            for item in payload["jobs"]:
                # Serialized coordinates are the evidence covered by the hash.
                # Validate them as unit vectors, then preserve their exact bits.
                direction = _stored_unit_direction(item["direction"], label="Stored Sidera direction")
                job = SideraEnsembleBurnJob(
                    scenario_id=str(item["scenario_id"]),
                    conditional_index=int(item["conditional_index"]),
                    scenario_probability=float(item["scenario_probability"]),
                    conditional_probability=float(item["conditional_probability"]),
                    joint_probability=float(item["joint_probability"]),
                    t_start_s=float(item["t_start_s"]),
                    direction=direction,
                    frame=str(item["frame"]),
                    artifact=TabulatedBurnArtifact.from_dict(item["artifact"]),
                )
                object.__setattr__(job, "direction", direction)
                jobs.append(job)
            handoff = cls(
                source_ensemble_hash=str(payload["source_ensemble_hash"]),
                jobs=tuple(jobs),
                joint_probability_sum=float(payload["joint_probability_sum"]),
                warnings=tuple(str(item) for item in payload["warnings"]),
                required_capabilities=tuple(
                    str(item) for item in payload["required_capabilities"]
                ),
                schema=str(payload["schema"]),
            )
        except InputError:
            raise
        except DomainError as error:
            raise InputError("Malformed Sidera ensemble hand-off values.") from error
        except (KeyError, TypeError, ValueError, IndexError) as error:
            raise InputError("Malformed Sidera ensemble hand-off.") from error
        for declared, job in zip(payload["jobs"], handoff.jobs, strict=True):
            if declared.get("artifact_hash") != job.artifact.artifact_hash:
                raise InputError("Sidera job artifact hash does not reproduce.")
            if declared["t_end_s"] != job.t_end_s:
                raise InputError("Sidera job end time does not reproduce.")
        if payload.get("handoff_hash") != handoff.handoff_hash:
            raise InputError("Sidera ensemble hand-off hash does not match its contents.")
        return handoff

    @classmethod
    def from_json(cls, text: str) -> SideraEnsembleHandoff:
        """Parse JSON text and integrity-check the hand-off document.

        Parameters
        ----------
        text : str
            JSON document text.

        Returns
        -------
        SideraEnsembleHandoff
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
        .. [3] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return read_json_document(text, "Sidera ensemble hand-off", cls.from_dict)


def export_sidera_ensemble_handoff(
    ensemble: HierarchicalPropulsionEnsemble,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
) -> SideraEnsembleHandoff:
    """Bind an integrity-checked propulsion ensemble to Sidera execution inputs.

    Parameters
    ----------
    ensemble : HierarchicalPropulsionEnsemble
        Source ensemble; it is round-tripped through its serialization and
        verified before use.
    t_start_s : float
        Absolute solver time at which every realization starts, s.
    direction : tuple of float
        Thrust direction in ``frame``; normalized per job.
    frame : {"inertial", "ric", "vnb"}, optional
        Direction frame.

    Returns
    -------
    SideraEnsembleHandoff
        One job per realization, with the source ensemble hash and the
        consumer obligations as warnings.

    Raises
    ------
    DomainError
        If the ensemble fails verification or a job is invalid.

    Notes
    -----
    The function does not construct native Sidera objects. Native
    consumption requires the public capabilities listed in
    ``SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES`` (see
    :func:`require_sidera_ensemble_capabilities`); this artifact is the
    producer side of that boundary. Joint probabilities are applied exactly
    once per realization [2]_.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures Guide
           for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369
    """

    verified = hierarchical_ensemble_from_dict(
        hierarchical_ensemble_to_dict(ensemble)
    )
    jobs = tuple(
        SideraEnsembleBurnJob(
            scenario_id=item.scenario_id,
            conditional_index=item.conditional_index,
            scenario_probability=item.scenario_probability,
            conditional_probability=item.conditional_probability,
            joint_probability=item.joint_probability,
            t_start_s=t_start_s,
            direction=direction,
            frame=frame,
            artifact=item.artifact,
        )
        for item in verified.realizations
    )
    return SideraEnsembleHandoff(
        source_ensemble_hash=hierarchical_ensemble_hash(verified),
        jobs=jobs,
        joint_probability_sum=sum(item.joint_probability for item in jobs),
        warnings=(
            "Sidera must split propagation at every artifact segment boundary.",
            "Use axial thrust for force and total tank flow for spacecraft mass rate.",
            "Apply each joint probability exactly once and retain zero-thrust jobs.",
            "Epoch, frame, direction, attitude, guidance, and trajectory remain Sidera responsibilities.",
        ),
    )


def missing_sidera_ensemble_capabilities(
    capabilities: SideraCapabilities,
) -> tuple[str, ...]:
    """Return the public Sidera capabilities still missing for exact consumption.

    Parameters
    ----------
    capabilities : SideraCapabilities
        Observed capabilities, e.g. from
        :func:`~rocket_propulsion.integrations.sidera.inspect_sidera_capabilities`.

    Returns
    -------
    tuple of str
        Missing entries of ``SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES`` among
        tabulated burns, independent mass flow and weighted ensembles; empty
        when none is missing.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    missing: list[str] = []
    if not capabilities.tabulated_burn_available:
        missing.append("tabulated_finite_burn")
    if not capabilities.independent_mass_flow_available:
        missing.append("independent_mass_flow")
    if not capabilities.weighted_maneuver_ensemble_available:
        missing.append("weighted_maneuver_ensemble")
    return tuple(missing)


def require_sidera_ensemble_capabilities(
    capabilities: SideraCapabilities,
) -> None:
    """Fail closed unless Sidera can consume the complete hand-off exactly.

    Parameters
    ----------
    capabilities : SideraCapabilities
        Observed capabilities.

    Raises
    ------
    FeatureUnavailableError
        If any required capability is missing; ``details`` lists the missing
        and observed capabilities.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    missing = missing_sidera_ensemble_capabilities(capabilities)
    if missing:
        raise FeatureUnavailableError(
            "Sidera cannot yet consume a weighted tabulated propulsion ensemble exactly.",
            details={
                "missing_capabilities": list(missing),
                "observed_capabilities": asdict(capabilities),
            },
        )

