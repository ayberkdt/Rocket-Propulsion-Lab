"""Future-facing weighted tabulated-burn hand-off for Sidera.

This module binds propulsion artifacts to execution time, direction, and frame
without importing or modifying Sidera. It is an adapter contract, not an orbit
propagator. Net axial thrust is the force channel; total tank flow is the mass
channel, so canted clusters and open-cycle engines remain physically closed.

References
----------
CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from math import isclose, isfinite, sqrt
from typing import Any, Literal

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.propulsion.burns import HierarchicalPropulsionEnsemble
from rocket_propulsion.propulsion.burns.ensemble_serialization import (
    HIERARCHICAL_ENSEMBLE_SCHEMA,
    hierarchical_ensemble_from_dict,
    hierarchical_ensemble_hash,
    hierarchical_ensemble_to_dict,
)
from rocket_propulsion.propulsion.burns.serialization import canonical_json_bytes
from rocket_propulsion.propulsion.burns.tabulated import TabulatedBurnArtifact

from .contracts import SideraCapabilities

SIDERA_ENSEMBLE_HANDOFF_SCHEMA = "rocket_propulsion_sidera_ensemble_handoff_v1"
SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES = (
    "tabulated_finite_burn",
    "independent_mass_flow",
    "segment_boundary_splitting",
    "weighted_maneuver_ensemble",
)
_FRAMES = {"inertial", "ric", "vnb"}


def _sha256_hex(value: str, field: str) -> str:
    normalized = str(value).lower()
    if len(normalized) != 64 or any(
        character not in "0123456789abcdef" for character in normalized
    ):
        raise DomainError(f"{field} must be a SHA-256 hex digest.")
    return normalized


def _unit_direction(values: tuple[float, float, float]) -> tuple[float, float, float]:
    vector = tuple(float(value) for value in values)
    if len(vector) != 3 or any(not isfinite(value) for value in vector):
        raise DomainError("Sidera ensemble direction must contain three finite values.")
    norm = sqrt(sum(value * value for value in vector))
    if norm <= 0.0:
        raise DomainError("Sidera ensemble direction must be non-zero.")
    return tuple(value / norm for value in vector)


@dataclass(frozen=True, slots=True)
class SideraEnsembleBurnJob:
    """One weighted Sidera propagation request backed by a full burn artifact.

    References
    ----------
    CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
        if not self.scenario_id.strip():
            raise DomainError("Sidera ensemble scenario ID must not be empty.")
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
        frame = str(self.frame).strip().lower()
        if frame not in _FRAMES:
            raise DomainError("Sidera ensemble frame must be inertial, ric, or vnb.")
        object.__setattr__(self, "frame", frame)
        object.__setattr__(self, "direction", _unit_direction(self.direction))

    @property
    def t_end_s(self) -> float:
        """Absolute solver time at the end of the propulsion history.

        References
        ----------
        CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        return self.t_start_s + self.artifact.duration_s


@dataclass(frozen=True, slots=True)
class SideraEnsembleHandoff:
    """Integrity-checked weighted transient jobs for a future Sidera consumer.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    source_ensemble_hash: str
    jobs: tuple[SideraEnsembleBurnJob, ...]
    joint_probability_sum: float
    warnings: tuple[str, ...]
    required_capabilities: tuple[str, ...] = SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES
    schema: str = SIDERA_ENSEMBLE_HANDOFF_SCHEMA

    def __post_init__(self) -> None:
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
        """SHA-256 identity covering all jobs, weights, and nested artifacts.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return hashlib.sha256(
            canonical_json_bytes(self.to_dict(include_hash=False))
        ).hexdigest()

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the canonical future-Sidera hand-off document.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        payload: dict[str, Any] = {
            "schema": self.schema,
            "source_ensemble_schema": HIERARCHICAL_ENSEMBLE_SCHEMA,
            "source_ensemble_hash": self.source_ensemble_hash,
            "force_channel": "axial_thrust_n",
            "mass_flow_channel": "total_tank_mass_flow_kg_s",
            "direction_semantics": "unit_vector_in_declared_frame",
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

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
        """Reopen a hand-off and verify nested artifact and outer hashes.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        if not isinstance(payload, dict):
            raise InputError("Sidera ensemble hand-off root must be an object.")
        if payload.get("schema") != SIDERA_ENSEMBLE_HANDOFF_SCHEMA:
            raise InputError("Unsupported Sidera ensemble hand-off schema.")
        expected_semantics = {
            "source_ensemble_schema": HIERARCHICAL_ENSEMBLE_SCHEMA,
            "force_channel": "axial_thrust_n",
            "mass_flow_channel": "total_tank_mass_flow_kg_s",
            "direction_semantics": "unit_vector_in_declared_frame",
        }
        if any(payload.get(key) != value for key, value in expected_semantics.items()):
            raise InputError("Sidera ensemble hand-off physical semantics changed.")
        try:
            jobs = tuple(
                SideraEnsembleBurnJob(
                    scenario_id=str(item["scenario_id"]),
                    conditional_index=int(item["conditional_index"]),
                    scenario_probability=float(item["scenario_probability"]),
                    conditional_probability=float(item["conditional_probability"]),
                    joint_probability=float(item["joint_probability"]),
                    t_start_s=float(item["t_start_s"]),
                    direction=tuple(float(value) for value in item["direction"]),
                    frame=str(item["frame"]),
                    artifact=TabulatedBurnArtifact.from_dict(item["artifact"]),
                )
                for item in payload["jobs"]
            )
            handoff = cls(
                source_ensemble_hash=str(payload["source_ensemble_hash"]),
                jobs=jobs,
                joint_probability_sum=float(payload["joint_probability_sum"]),
                warnings=tuple(str(item) for item in payload["warnings"]),
                required_capabilities=tuple(
                    str(item) for item in payload["required_capabilities"]
                ),
                schema=str(payload["schema"]),
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, IndexError) as error:
            raise InputError("Malformed Sidera ensemble hand-off.") from error
        for declared, job in zip(payload["jobs"], handoff.jobs, strict=True):
            if declared.get("artifact_hash") != job.artifact.artifact_hash:
                raise InputError("Sidera job artifact hash does not reproduce.")
            if not isclose(
                float(declared.get("t_end_s")),
                job.t_end_s,
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                raise InputError("Sidera job end time does not reproduce.")
        if payload.get("handoff_hash") != handoff.handoff_hash:
            raise InputError("Sidera ensemble hand-off hash does not match its contents.")
        return handoff

    @classmethod
    def from_json(cls, text: str) -> SideraEnsembleHandoff:
        """Parse and integrity-check a future-Sidera hand-off document.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        try:
            payload = json.loads(text)
        except json.JSONDecodeError as error:
            raise InputError("Sidera ensemble hand-off is not valid JSON.") from error
        return cls.from_dict(payload)


def export_sidera_ensemble_handoff(
    ensemble: HierarchicalPropulsionEnsemble,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
) -> SideraEnsembleHandoff:
    """Bind an integrity-checked propulsion ensemble to Sidera execution inputs.

    The function deliberately does not construct native Sidera objects. Current
    Sidera lacks the required public weighted tabulated-burn contract; this
    artifact is the isolated, reviewed producer side of that future boundary.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """Return public capabilities still missing for exact native consumption.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

