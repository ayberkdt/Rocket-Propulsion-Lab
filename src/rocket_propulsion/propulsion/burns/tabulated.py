"""Target-neutral, discontinuity-safe tabulated burn exchange artifacts.

The artifact in this module is the isolated contract intended for eventual
Sidera integration.  It contains propulsion force and tank-drain histories,
but intentionally contains no epoch, direction, frame, orbit, or propagator
state.  Segment endpoints are stored independently so instantaneous command
jumps are not smeared or represented by duplicate timestamps.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from dataclasses import dataclass
from enum import Enum
from itertools import pairwise
from math import isfinite
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .models import ProfiledBurnResult
from .profiles import resolved_schedule_intervals
from .serialization import canonical_json_bytes

_ARTIFACT_SCALE_NAMES = {
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
}


class InterpolationPolicy(str, Enum):
    """Interpolation applied inside each declared exchange segment.

    References
    ----------
    CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    LINEAR = "linear"
    HOLD = "hold"


def _nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _flow_pairs(value: tuple[tuple[str, float], ...], label: str) -> None:
    identifiers = [identifier for identifier, _ in value]
    if any(not identifier.strip() for identifier in identifiers):
        raise DomainError(f"{label} identifiers must not be empty.")
    if len(identifiers) != len(set(identifiers)):
        raise DomainError(f"{label} identifiers must be unique.")
    for _, flow in value:
        _nonnegative(flow, label)


@dataclass(frozen=True, slots=True)
class TabulatedBurnSegment:
    """One time segment with independent left/right propulsion states.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    start_time_s: float
    end_time_s: float
    start_delivered_thrust_n: float
    end_delivered_thrust_n: float
    start_axial_thrust_n: float
    end_axial_thrust_n: float
    start_total_mass_flow_kg_s: float
    end_total_mass_flow_kg_s: float
    start_stream_mass_flows_kg_s: tuple[tuple[str, float], ...] = ()
    end_stream_mass_flows_kg_s: tuple[tuple[str, float], ...] = ()
    phase: str = "tabulated"

    def __post_init__(self) -> None:
        if not isfinite(self.start_time_s) or self.start_time_s < 0.0:
            raise DomainError("Tabulated segment start time must be finite and non-negative.")
        if not isfinite(self.end_time_s) or self.end_time_s <= self.start_time_s:
            raise DomainError("Tabulated segment end time must be after its start time.")
        for value, label in (
            (self.start_delivered_thrust_n, "Start delivered thrust"),
            (self.end_delivered_thrust_n, "End delivered thrust"),
            (self.start_axial_thrust_n, "Start axial thrust"),
            (self.end_axial_thrust_n, "End axial thrust"),
            (self.start_total_mass_flow_kg_s, "Start total mass flow"),
            (self.end_total_mass_flow_kg_s, "End total mass flow"),
        ):
            _nonnegative(value, label)
        if self.start_axial_thrust_n > self.start_delivered_thrust_n * (1.0 + 1e-12):
            raise DomainError("Start axial thrust cannot exceed delivered thrust.")
        if self.end_axial_thrust_n > self.end_delivered_thrust_n * (1.0 + 1e-12):
            raise DomainError("End axial thrust cannot exceed delivered thrust.")
        _flow_pairs(self.start_stream_mass_flows_kg_s, "Start stream mass flow")
        _flow_pairs(self.end_stream_mass_flows_kg_s, "End stream mass flow")
        self._validate_stream_closure(
            self.start_stream_mass_flows_kg_s,
            self.start_total_mass_flow_kg_s,
            "start",
        )
        self._validate_stream_closure(
            self.end_stream_mass_flows_kg_s,
            self.end_total_mass_flow_kg_s,
            "end",
        )
        if not self.phase.strip():
            raise DomainError("Tabulated segment phase must not be empty.")

    @staticmethod
    def _validate_stream_closure(
        streams: tuple[tuple[str, float], ...], total: float, side: str
    ) -> None:
        if not streams:
            return
        actual = sum(flow for _, flow in streams)
        if abs(actual - total) > max(1e-12, total * 1e-10):
            raise DomainError(f"Tabulated {side} stream flows do not close to total flow.")

    @property
    def duration_s(self) -> float:
        return self.end_time_s - self.start_time_s

    @staticmethod
    def _integral(start: float, end: float, duration: float) -> float:
        return 0.5 * (start + end) * duration

    @staticmethod
    def _local_first_moment(start: float, end: float, duration: float) -> float:
        return duration * duration * (start / 2.0 + (end - start) / 3.0)

    @property
    def delivered_impulse_n_s(self) -> float:
        return self._integral(
            self.start_delivered_thrust_n, self.end_delivered_thrust_n, self.duration_s
        )

    @property
    def axial_impulse_n_s(self) -> float:
        return self._integral(self.start_axial_thrust_n, self.end_axial_thrust_n, self.duration_s)

    @property
    def propellant_mass_kg(self) -> float:
        return self._integral(
            self.start_total_mass_flow_kg_s,
            self.end_total_mass_flow_kg_s,
            self.duration_s,
        )

    @property
    def delivered_thrust_first_moment_n_s2(self) -> float:
        return self.start_time_s * self.delivered_impulse_n_s + self._local_first_moment(
            self.start_delivered_thrust_n,
            self.end_delivered_thrust_n,
            self.duration_s,
        )

    def state_at(self, time_s: float) -> dict[str, Any]:
        """Linearly interpolate the propulsion state inside this segment."""

        if not isfinite(time_s) or not self.start_time_s <= time_s <= self.end_time_s:
            raise DomainError("Tabulated evaluation time lies outside the segment.")
        fraction = (time_s - self.start_time_s) / self.duration_s

        def interpolate(start: float, end: float) -> float:
            return start + fraction * (end - start)

        start_streams = dict(self.start_stream_mass_flows_kg_s)
        end_streams = dict(self.end_stream_mass_flows_kg_s)
        identifiers = sorted(start_streams.keys() | end_streams.keys())
        return {
            "time_s": time_s,
            "delivered_thrust_n": interpolate(
                self.start_delivered_thrust_n, self.end_delivered_thrust_n
            ),
            "axial_thrust_n": interpolate(self.start_axial_thrust_n, self.end_axial_thrust_n),
            "total_mass_flow_kg_s": interpolate(
                self.start_total_mass_flow_kg_s, self.end_total_mass_flow_kg_s
            ),
            "stream_mass_flows_kg_s": tuple(
                (
                    identifier,
                    interpolate(
                        start_streams.get(identifier, 0.0), end_streams.get(identifier, 0.0)
                    ),
                )
                for identifier in identifiers
            ),
            "phase": self.phase,
        }


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


@dataclass(frozen=True, slots=True)
class TabulatedBurnArtifact:
    """Immutable propulsion history ready for a future dynamics consumer.

    ``reference_ids`` records the reviewed engineering sources supporting the
    producer model. It is included in the content hash, so citations cannot be
    changed without changing artifact identity.

    References
    ----------
    CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    segments: tuple[TabulatedBurnSegment, ...]
    source_id: str
    source_sha256: str
    interpolation: InterpolationPolicy = InterpolationPolicy.LINEAR
    reference_ids: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    schema: str = "rocket_propulsion_tabulated_burn_v1"

    def __post_init__(self) -> None:
        if self.schema != "rocket_propulsion_tabulated_burn_v1":
            raise DomainError("Unsupported tabulated burn schema.")
        if not self.segments:
            raise DomainError("A tabulated burn requires at least one segment.")
        if abs(self.segments[0].start_time_s) > 1e-12:
            raise DomainError("A tabulated burn must begin at zero relative time.")
        for previous, current in zip(self.segments, self.segments[1:], strict=False):
            if abs(previous.end_time_s - current.start_time_s) > 1e-12:
                raise DomainError("Tabulated burn segments must be contiguous.")
        if self.interpolation is InterpolationPolicy.HOLD:
            for segment in self.segments:
                if any(
                    abs(start - end) > max(1e-12, abs(start) * 1e-12)
                    for start, end in (
                        (
                            segment.start_delivered_thrust_n,
                            segment.end_delivered_thrust_n,
                        ),
                        (segment.start_axial_thrust_n, segment.end_axial_thrust_n),
                        (
                            segment.start_total_mass_flow_kg_s,
                            segment.end_total_mass_flow_kg_s,
                        ),
                    )
                ) or dict(segment.start_stream_mass_flows_kg_s) != dict(
                    segment.end_stream_mass_flows_kg_s
                ):
                    raise DomainError(
                        "Hold interpolation requires constant values inside each segment."
                    )
        if not self.source_id.strip():
            raise DomainError("Tabulated burn source identifier must not be empty.")
        if len(self.source_sha256) != 64 or any(
            character not in "0123456789abcdef" for character in self.source_sha256.lower()
        ):
            raise DomainError("Tabulated burn source hash must be a SHA-256 hex digest.")
        if any(not reference_id.strip() for reference_id in self.reference_ids):
            raise DomainError("Tabulated burn reference identifiers must not be empty.")
        if len(self.reference_ids) != len(set(self.reference_ids)):
            raise DomainError("Tabulated burn reference identifiers must be unique.")

    @property
    def duration_s(self) -> float:
        return self.segments[-1].end_time_s

    @property
    def delivered_total_impulse_n_s(self) -> float:
        return sum(segment.delivered_impulse_n_s for segment in self.segments)

    @property
    def axial_total_impulse_n_s(self) -> float:
        return sum(segment.axial_impulse_n_s for segment in self.segments)

    @property
    def consumed_propellant_kg(self) -> float:
        return sum(segment.propellant_mass_kg for segment in self.segments)

    @property
    def thrust_centroid_time_s(self) -> float:
        impulse = self.delivered_total_impulse_n_s
        return (
            sum(segment.delivered_thrust_first_moment_n_s2 for segment in self.segments) / impulse
            if impulse > 0.0
            else 0.0
        )

    @property
    def system_equivalent_specific_impulse_s(self) -> float | None:
        mass = self.consumed_propellant_kg
        return (
            self.delivered_total_impulse_n_s / (STANDARD_GRAVITY_M_S2 * mass)
            if mass > 0.0
            else None
        )

    @property
    def artifact_hash(self) -> str:
        return hashlib.sha256(_canonical_json_bytes(self.to_dict(include_hash=False))).hexdigest()

    def state_at(self, time_s: float, *, side: str = "right") -> dict[str, Any]:
        """Evaluate the artifact with explicit left/right discontinuity semantics."""

        if side not in {"left", "right"}:
            raise DomainError("Tabulated boundary side must be 'left' or 'right'.")
        if not isfinite(time_s) or not 0.0 <= time_s <= self.duration_s:
            raise DomainError("Tabulated burn time lies outside the artifact.")
        tolerance = max(1e-12, self.duration_s * 1e-12)
        for index, segment in enumerate(self.segments):
            if segment.start_time_s <= time_s < segment.end_time_s:
                if side == "left" and index > 0 and abs(time_s - segment.start_time_s) <= tolerance:
                    return self.segments[index - 1].state_at(time_s)
                return segment.state_at(time_s)
        return self.segments[-1].state_at(self.duration_s)

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return a deterministic JSON-compatible artifact."""

        payload: dict[str, Any] = {
            "schema": self.schema,
            "source_id": self.source_id,
            "source_sha256": self.source_sha256,
            "interpolation": self.interpolation.value,
            "segments": [
                {
                    "start_time_s": segment.start_time_s,
                    "end_time_s": segment.end_time_s,
                    "start_delivered_thrust_n": segment.start_delivered_thrust_n,
                    "end_delivered_thrust_n": segment.end_delivered_thrust_n,
                    "start_axial_thrust_n": segment.start_axial_thrust_n,
                    "end_axial_thrust_n": segment.end_axial_thrust_n,
                    "start_total_mass_flow_kg_s": segment.start_total_mass_flow_kg_s,
                    "end_total_mass_flow_kg_s": segment.end_total_mass_flow_kg_s,
                    "start_stream_mass_flows_kg_s": [
                        list(item) for item in segment.start_stream_mass_flows_kg_s
                    ],
                    "end_stream_mass_flows_kg_s": [
                        list(item) for item in segment.end_stream_mass_flows_kg_s
                    ],
                    "phase": segment.phase,
                }
                for segment in self.segments
            ],
            "summary": {
                "duration_s": self.duration_s,
                "delivered_total_impulse_n_s": self.delivered_total_impulse_n_s,
                "axial_total_impulse_n_s": self.axial_total_impulse_n_s,
                "consumed_propellant_kg": self.consumed_propellant_kg,
                "thrust_centroid_time_s": self.thrust_centroid_time_s,
                "system_equivalent_specific_impulse_s": (self.system_equivalent_specific_impulse_s),
            },
            "warnings": list(self.warnings),
        }
        if self.reference_ids:
            payload["reference_ids"] = list(self.reference_ids)
        if include_hash:
            payload["artifact_hash"] = self.artifact_hash
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(
            self.to_dict(), ensure_ascii=False, allow_nan=False, sort_keys=True, indent=indent
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> TabulatedBurnArtifact:
        """Reopen an artifact and verify its content hash and derived summary."""

        if not isinstance(payload, dict):
            raise InputError("Tabulated burn artifact root must be an object.")
        try:
            artifact = cls(
                schema=str(payload["schema"]),
                source_id=str(payload["source_id"]),
                source_sha256=str(payload["source_sha256"]),
                interpolation=InterpolationPolicy(payload["interpolation"]),
                reference_ids=tuple(str(item) for item in payload.get("reference_ids", ())),
                segments=tuple(
                    TabulatedBurnSegment(
                        start_time_s=float(item["start_time_s"]),
                        end_time_s=float(item["end_time_s"]),
                        start_delivered_thrust_n=float(item["start_delivered_thrust_n"]),
                        end_delivered_thrust_n=float(item["end_delivered_thrust_n"]),
                        start_axial_thrust_n=float(item["start_axial_thrust_n"]),
                        end_axial_thrust_n=float(item["end_axial_thrust_n"]),
                        start_total_mass_flow_kg_s=float(item["start_total_mass_flow_kg_s"]),
                        end_total_mass_flow_kg_s=float(item["end_total_mass_flow_kg_s"]),
                        start_stream_mass_flows_kg_s=tuple(
                            (str(pair[0]), float(pair[1]))
                            for pair in item["start_stream_mass_flows_kg_s"]
                        ),
                        end_stream_mass_flows_kg_s=tuple(
                            (str(pair[0]), float(pair[1]))
                            for pair in item["end_stream_mass_flows_kg_s"]
                        ),
                        phase=str(item["phase"]),
                    )
                    for item in payload["segments"]
                ),
                warnings=tuple(str(item) for item in payload.get("warnings", ())),
            )
        except (KeyError, TypeError, ValueError, IndexError) as error:
            raise InputError("Malformed tabulated burn artifact.") from error
        if payload.get("artifact_hash") != artifact.artifact_hash:
            raise InputError("Tabulated burn artifact hash does not match its contents.")
        declared_summary = payload.get("summary")
        if declared_summary != artifact.to_dict(include_hash=False)["summary"]:
            raise InputError("Tabulated burn artifact summary does not reproduce.")
        return artifact

    @classmethod
    def from_json(cls, text: str) -> TabulatedBurnArtifact:
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as error:
            raise InputError("Tabulated burn artifact is not valid JSON.") from error
        return cls.from_dict(payload)


def scale_tabulated_artifact(
    artifact: TabulatedBurnArtifact,
    scales: dict[str, float],
    *,
    source_id: str,
    provenance: Any,
    reference_ids: tuple[str, ...] = (),
    warnings: tuple[str, ...] = (),
) -> TabulatedBurnArtifact:
    """Apply physically coupled positive scales to a propulsion artifact.

    Delivered thrust scales with ``thrust_scale``; tank and named stream flows
    scale with ``thrust_scale / specific_impulse_scale``; time scales every
    segment boundary; and axial efficiency affects axial force without changing
    delivered thrust or tank drain. The resulting source hash binds the base
    artifact, normalized scales, and caller-supplied provenance.

    References
    ----------
    NASA rocket thrust equation:
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    unknown = set(scales) - _ARTIFACT_SCALE_NAMES
    if unknown:
        raise DomainError(f"Unsupported artifact scales: {', '.join(sorted(unknown))}.")
    normalized = {name: 1.0 for name in _ARTIFACT_SCALE_NAMES}
    for name, value in scales.items():
        try:
            numeric = float(value)
        except (TypeError, ValueError) as error:
            raise DomainError(f"Artifact scale {name!r} must be numeric.") from error
        if not isfinite(numeric) or numeric <= 0.0:
            raise DomainError(f"Artifact scale {name!r} must be positive and finite.")
        normalized[name] = numeric
    if not source_id.strip():
        raise DomainError("Scaled artifact source identifier must not be empty.")
    thrust_scale = normalized["thrust_scale"]
    isp_scale = normalized["specific_impulse_scale"]
    time_scale = normalized["time_scale"]
    axial_scale = normalized["axial_efficiency_scale"]
    maximum_axial_ratio = max(
        (
            axial / delivered
            for segment in artifact.segments
            for axial, delivered in (
                (segment.start_axial_thrust_n, segment.start_delivered_thrust_n),
                (segment.end_axial_thrust_n, segment.end_delivered_thrust_n),
            )
            if delivered > 0.0
        ),
        default=0.0,
    )
    if maximum_axial_ratio * axial_scale > 1.0 + 1e-12:
        raise DomainError("Artifact scale can make axial thrust exceed delivered thrust.")
    flow_scale = thrust_scale / isp_scale

    def scaled_streams(
        streams: tuple[tuple[str, float], ...],
    ) -> tuple[tuple[str, float], ...]:
        return tuple((tank_id, flow * flow_scale) for tank_id, flow in streams)

    segments = tuple(
        TabulatedBurnSegment(
            start_time_s=segment.start_time_s * time_scale,
            end_time_s=segment.end_time_s * time_scale,
            start_delivered_thrust_n=segment.start_delivered_thrust_n * thrust_scale,
            end_delivered_thrust_n=segment.end_delivered_thrust_n * thrust_scale,
            start_axial_thrust_n=(segment.start_axial_thrust_n * thrust_scale * axial_scale),
            end_axial_thrust_n=(segment.end_axial_thrust_n * thrust_scale * axial_scale),
            start_total_mass_flow_kg_s=(segment.start_total_mass_flow_kg_s * flow_scale),
            end_total_mass_flow_kg_s=(segment.end_total_mass_flow_kg_s * flow_scale),
            start_stream_mass_flows_kg_s=scaled_streams(segment.start_stream_mass_flows_kg_s),
            end_stream_mass_flows_kg_s=scaled_streams(segment.end_stream_mass_flows_kg_s),
            phase=segment.phase,
        )
        for segment in artifact.segments
    )
    source_sha256 = hashlib.sha256(
        canonical_json_bytes(
            {
                "base_artifact_hash": artifact.artifact_hash,
                "scales": normalized,
                "provenance": provenance,
            }
        )
    ).hexdigest()
    return TabulatedBurnArtifact(
        segments=segments,
        source_id=source_id,
        source_sha256=source_sha256,
        interpolation=artifact.interpolation,
        reference_ids=tuple(dict.fromkeys(artifact.reference_ids + reference_ids)),
        warnings=artifact.warnings + warnings,
    )


def tabulated_artifact_from_profiled_burn(result: ProfiledBurnResult) -> TabulatedBurnArtifact:
    """Convert an analytic L1 result without flattening discontinuities.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    NASA rocket thrust equation:
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    """

    tank_rates: dict[str, float] = {}
    for stream in result.operating_point.streams:
        if stream.tank_depleting:
            tank_rates[stream.tank_id] = tank_rates.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s

    def flows(throttle: float) -> tuple[tuple[str, float], ...]:
        return tuple((tank_id, rate * throttle) for tank_id, rate in sorted(tank_rates.items()))

    point = result.operating_point
    cant = result.definition.cant_efficiency
    segments = tuple(
        TabulatedBurnSegment(
            start_time_s=interval.start_time_s,
            end_time_s=interval.end_time_s,
            start_delivered_thrust_n=point.delivered_thrust_n * interval.start_throttle,
            end_delivered_thrust_n=point.delivered_thrust_n * interval.end_throttle,
            start_axial_thrust_n=point.delivered_thrust_n * cant * interval.start_throttle,
            end_axial_thrust_n=point.delivered_thrust_n * cant * interval.end_throttle,
            start_total_mass_flow_kg_s=point.total_tank_flow_kg_s * interval.start_throttle,
            end_total_mass_flow_kg_s=point.total_tank_flow_kg_s * interval.end_throttle,
            start_stream_mass_flows_kg_s=flows(interval.start_throttle),
            end_stream_mass_flows_kg_s=flows(interval.end_throttle),
            phase=interval.phase,
        )
        for interval in resolved_schedule_intervals(
            result.schedule, stop_time_s=result.summary.duration_s
        )
    )
    artifact = TabulatedBurnArtifact(
        segments=segments,
        source_id=f"profiled-burn:{result.definition.name}",
        source_sha256=result.result_hash,
        reference_ids=("NASA-SP-125", "NASA-GRC-THRUST-EQUATION"),
        warnings=result.warnings,
    )
    for actual, expected, label in (
        (
            artifact.delivered_total_impulse_n_s,
            result.summary.delivered_total_impulse_n_s,
            "delivered impulse",
        ),
        (
            artifact.consumed_propellant_kg,
            result.summary.consumed_propellant_kg,
            "propellant mass",
        ),
    ):
        if abs(actual - expected) > max(1e-10, abs(expected) * 1e-11):
            raise DomainError(f"Tabulated L1 conversion failed {label} closure.")
    return artifact


def tabulated_artifact_from_csv(
    text: str,
    *,
    source_id: str,
    interpolation: InterpolationPolicy = InterpolationPolicy.LINEAR,
) -> TabulatedBurnArtifact:
    """Import an SI trace without smoothing or silently changing its samples.

    Required columns are ``time_s``, ``delivered_thrust_n`` and
    ``total_mass_flow_kg_s``. ``axial_thrust_n`` defaults to delivered thrust.
    Any column named ``stream:<tank-id>`` is preserved as an explicit drain.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA-CR-140800: https://ntrs.nasa.gov/citations/19750003988
    """

    if not isinstance(text, str) or not text.strip():
        raise InputError("Imported burn trace must not be empty.")
    source_hash = hashlib.sha256(text.encode("utf-8")).hexdigest()
    reader = csv.DictReader(io.StringIO(text))
    required = {"time_s", "delivered_thrust_n", "total_mass_flow_kg_s"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise InputError(
            "Imported trace requires time_s, delivered_thrust_n, and total_mass_flow_kg_s columns."
        )
    stream_columns = sorted(name for name in reader.fieldnames if name.startswith("stream:"))
    points: list[dict[str, Any]] = []
    try:
        for row in reader:
            delivered = float(row["delivered_thrust_n"])
            points.append(
                {
                    "time": float(row["time_s"]),
                    "delivered": delivered,
                    "axial": float(row.get("axial_thrust_n") or delivered),
                    "flow": float(row["total_mass_flow_kg_s"]),
                    "streams": tuple(
                        (column.removeprefix("stream:"), float(row[column] or 0.0))
                        for column in stream_columns
                    ),
                    "phase": str(row.get("phase") or "imported-trace"),
                }
            )
    except (TypeError, ValueError) as error:
        raise InputError("Imported trace contains a non-numeric physical value.") from error
    if len(points) < 2:
        raise InputError("Imported trace requires at least two samples.")
    if abs(points[0]["time"]) > 1e-12:
        raise InputError("Imported trace must begin at time_s = 0.")
    if any(right["time"] <= left["time"] for left, right in pairwise(points)):
        raise InputError("Imported trace times must be strictly increasing.")

    segments: list[TabulatedBurnSegment] = []
    for start, end in pairwise(points):
        held = interpolation is InterpolationPolicy.HOLD
        segments.append(
            TabulatedBurnSegment(
                start_time_s=start["time"],
                end_time_s=end["time"],
                start_delivered_thrust_n=start["delivered"],
                end_delivered_thrust_n=(start["delivered"] if held else end["delivered"]),
                start_axial_thrust_n=start["axial"],
                end_axial_thrust_n=(start["axial"] if held else end["axial"]),
                start_total_mass_flow_kg_s=start["flow"],
                end_total_mass_flow_kg_s=(start["flow"] if held else end["flow"]),
                start_stream_mass_flows_kg_s=start["streams"],
                end_stream_mass_flows_kg_s=(start["streams"] if held else end["streams"]),
                phase=start["phase"],
            )
        )
    return TabulatedBurnArtifact(
        segments=tuple(segments),
        source_id=source_id,
        source_sha256=source_hash,
        interpolation=interpolation,
        warnings=("Imported trace samples were preserved without smoothing.",),
    )
