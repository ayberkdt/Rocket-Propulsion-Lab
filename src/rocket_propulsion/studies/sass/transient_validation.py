"""Replay holdout measurements from a byte-identified acquisition export."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from json import dumps
from math import ulp
from typing import Literal, cast

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, read_json_document, records, text

from .command_life import ThrusterCommandHistory, ValveCommand
from .evidence import CitedScalar
from .project import scalar_from_dict
from .rational_numerics import outward_endpoint
from .shutdown_observation import ShutdownWindow, review_shutdown_impulse
from .transient_observation import TransientSample, TransientTrace, review_transient
from .valve_response import (
    ObservedValveTrace,
    ValveLatencyLimits,
    ValveResponseReview,
    review_valve_responses,
)


@dataclass(frozen=True, slots=True)
class PressureReferent:
    """Physical pressure location and reference shared by measurement and model.

    Differential pressure is location minus reference_location. A gauge
    reference identifies the ambient or other gauge reference for the run.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, identity of referents.
    .. [2] BIPM-SI, pressure units.
    """

    location: str
    reference: Literal["absolute", "gauge", "differential"]
    reference_location: str | None = None

    def __post_init__(self) -> None:
        """Require explicit identities for the positive and any reference location.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1 and 4.2.5, input and referent traceability.
        """
        text(self.location, "pressure measurement location")
        if self.reference not in ("absolute", "gauge", "differential"):
            raise DomainError("Pressure referent requires an explicit reference.")
        if self.reference == "absolute":
            if self.reference_location is not None:
                raise DomainError("Absolute pressure does not take a second location.")
        else:
            reference = text(self.reference_location, "pressure reference location")
            if reference == self.location:
                raise DomainError("Pressure locations must distinguish the two sides.")


def _pressure_referent_from_payload(raw: object) -> PressureReferent | None:
    """Read the exact referent schema without adding a default pressure datum.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, source representation.
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or raw.keys() != {"location", "reference", "reference_location"}:
        raise DomainError("Pressure referent has an invalid schema.")
    return PressureReferent(**raw)


@dataclass(frozen=True, slots=True)
class ArchivedTransientChannel:
    """One acquisition channel from a named experimental run.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, identity of experimental referents.
    """

    condition_record: str
    trace: TransientTrace
    pressure_referent: PressureReferent | None = None

    def __post_init__(self) -> None:
        """Require an unambiguous run identity and typed trace.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, input configuration records.
        """
        text(self.condition_record, "archived transient condition")
        if not isinstance(self.trace, TransientTrace):
            raise DomainError("An archived channel requires a transient trace.")
        if self.pressure_referent is not None and (
            self.trace.unit != "Pa" or not isinstance(self.pressure_referent, PressureReferent)
        ):
            raise DomainError("An archived pressure referent requires a Pa channel.")


@dataclass(frozen=True, slots=True)
class ArchivedValveAcquisition:
    """Issued and observed valve histories from one condition in the export.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, identity of measured referents.
    """

    condition_record: str
    commanded: ThrusterCommandHistory
    observed: ObservedValveTrace

    def __post_init__(self) -> None:
        """Require a named run and typed command and observation records.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, data configuration and pedigree.
        """
        text(self.condition_record, "archived valve condition")
        if not isinstance(self.commanded, ThrusterCommandHistory) or not isinstance(
            self.observed, ObservedValveTrace
        ):
            raise DomainError("An archived valve acquisition requires command and observation histories.")


@dataclass(frozen=True, slots=True)
class ArchivedTransientRun:
    """Canonical acquisition export for one article and its holdout runs.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, input pedigree and referents.
    """

    test_article_id: str
    channels: tuple[ArchivedTransientChannel, ...]
    valve_acquisitions: tuple[ArchivedValveAcquisition, ...] = ()

    def __post_init__(self) -> None:
        """Reject duplicate run/channel pairs and freeze deterministic order.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, traceable data identifiers.
        """
        text(self.test_article_id, "archived transient test article ID")
        channels: tuple[ArchivedTransientChannel, ...] = records(
            self.channels, ArchivedTransientChannel, "archived transient channels"
        )
        keys = {(row.condition_record, row.trace.channel_id) for row in channels}
        if len(keys) != len(channels):
            raise DomainError("An acquisition export requires unique run and channel pairs.")
        valves: tuple[ArchivedValveAcquisition, ...] = records(
            self.valve_acquisitions, ArchivedValveAcquisition, "archived valve acquisitions"
        )
        if not channels and not valves:
            raise DomainError("An acquisition export requires measurement or valve records.")
        if len({row.condition_record for row in valves}) != len(valves):
            raise DomainError("An acquisition export requires unique valve condition records.")
        object.__setattr__(self, "channels", tuple(sorted(channels, key=_channel_key)))
        object.__setattr__(self, "valve_acquisitions", tuple(sorted(valves, key=_valve_condition)))


def _channel_key(channel: ArchivedTransientChannel) -> tuple[str, str]:
    """Identify a channel within a physical test run.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, identity of empirical referents.
    """
    return channel.condition_record, channel.trace.channel_id


def _valve_condition(row: ArchivedValveAcquisition) -> str:
    """Order the valve acquisitions by experiment identity.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, data configuration.
    """
    return row.condition_record


def _valve_archive_row(row: ArchivedValveAcquisition) -> dict:
    """Encode both ledgers, including their independent cited time bases.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.7, acquisition pedigree and timing.
    """
    data = asdict(row)
    data["commanded"]["tick_period"]["value"] = row.commanded.tick_period.value.hex()
    data["observed"]["history"]["tick_period"]["value"] = row.observed.history.tick_period.value.hex()
    return data


def transient_archive_bytes(run: ArchivedTransientRun) -> bytes:
    """Encode the interpreted acquisition export that can be replayed.

    Preserve the original DAQ files and the reviewed DAQ-to-export mapping.
    The export uses sorted JSON keys, ASCII escapes and exact binary64 hex
    strings for every numeric measurement and bound. The top-level schema
    also identifies the test article. Hash these exact bytes as raw data
    in a validation evidence record using acquisition replay.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, input data pedigree and referents.
    """
    if not isinstance(run, ArchivedTransientRun):
        raise DomainError("Expected a typed acquisition export.")
    channels = []
    for channel in run.channels:
        trace = channel.trace
        period = asdict(trace.tick_period)
        period["value"] = trace.tick_period.value.hex()
        slew = None
        if trace.max_absolute_slew is not None:
            slew = asdict(trace.max_absolute_slew)
            slew["value"] = trace.max_absolute_slew.value.hex()
        fields = {
            "channel_id": trace.channel_id,
            "unit": trace.unit,
            "tick_period": period,
            "max_absolute_slew": slew,
            "source_record": trace.source_record,
            "acquisition_review_record": trace.acquisition_review_record,
            "max_gap_ticks": trace.max_gap_ticks,
            "samples": [[sample.tick, sample.value.hex(), sample.absolute_error.hex()]
                        for sample in trace.samples],
        }
        channel_data: dict[str, object] = {"condition_record": channel.condition_record, "trace": fields}
        if channel.pressure_referent is not None:
            channel_data["pressure_referent"] = asdict(channel.pressure_referent)
        channels.append(channel_data)
    data: dict[str, object] = {"schema": "rocket-propulsion-transient-run-v1",
                               "test_article_id": run.test_article_id, "channels": channels}
    if run.valve_acquisitions:
        data["schema"] = "rocket-propulsion-transient-run-v2"
        data["valve_acquisitions"] = [_valve_archive_row(row) for row in run.valve_acquisitions]
    if any(channel.pressure_referent is not None for channel in run.channels):
        data["schema"] = "rocket-propulsion-transient-run-v3"
    return dumps(data, sort_keys=True, ensure_ascii=True, allow_nan=False,
                 separators=(",", ":")).encode("ascii")


def _scalar_from_hex(raw: object) -> CitedScalar:
    """Read a cited scalar with the archive's hexadecimal value encoding.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, source records and units.
    """
    if not isinstance(raw, dict):
        raise DomainError("Archived acquisition scalar requires a cited envelope.")
    return scalar_from_dict({**raw, "value": float.fromhex(raw["value"])})


def _read_valve_history(raw: dict) -> ThrusterCommandHistory:
    """Decode a valve state ledger through the same typed admission rules.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, data representation and input checks.
    """
    if not isinstance(raw, dict) or not isinstance(raw.get("commands"), list):
        raise DomainError("An archived valve history requires a command array.")
    return ThrusterCommandHistory(
        **{**raw, "tick_period": _scalar_from_hex(raw["tick_period"]),
           "commands": tuple(ValveCommand(**row) for row in raw["commands"])}
    )


def _read_valve_acquisition(raw: dict) -> ArchivedValveAcquisition:
    """Reconstruct paired histories and their coverage and alignment review.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, scope and clock error.
    """
    if not isinstance(raw, dict) or not isinstance(raw.get("observed"), dict):
        raise DomainError("An archived valve acquisition requires an observation record.")
    observed = raw["observed"]
    return ArchivedValveAcquisition(
        raw["condition_record"], _read_valve_history(raw["commanded"]),
        ObservedValveTrace(**{**observed, "history": _read_valve_history(observed["history"])}),
    )


def _transient_run_from_payload(payload: object) -> ArchivedTransientRun:
    """Build the typed acquisition export from decoded JSON objects.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, traceable data and representation.
    """
    if not isinstance(payload, dict) or payload.get("schema") not in (
        "rocket-propulsion-transient-run-v1", "rocket-propulsion-transient-run-v2",
        "rocket-propulsion-transient-run-v3",
    ):
        raise DomainError("Unknown acquisition archive schema.")
    channels = []
    if not isinstance(payload["channels"], list):
        raise DomainError("Acquisition channels must be an array.")
    for row in payload["channels"]:
        if not isinstance(row, dict) or not isinstance(row["trace"], dict):
            raise DomainError("Each archived channel requires a trace object.")
        raw = row["trace"]
        if not isinstance(raw["samples"], list):
            raise DomainError("Archived samples must be an array.")
        samples = []
        for item in raw["samples"]:
            if not isinstance(item, list) or len(item) != 3:
                raise DomainError("Archived samples require tick, value and absolute error.")
            samples.append(TransientSample(item[0], float.fromhex(item[1]), float.fromhex(item[2])))
        slew = None if raw["max_absolute_slew"] is None else _scalar_from_hex(raw["max_absolute_slew"])
        trace = TransientTrace(
            raw["channel_id"], raw["unit"], _scalar_from_hex(raw["tick_period"]), tuple(samples),
            raw["source_record"], raw["acquisition_review_record"], raw["max_gap_ticks"], slew,
        )
        channels.append(ArchivedTransientChannel(
            row["condition_record"], trace, _pressure_referent_from_payload(row.get("pressure_referent"))
        ))
    valves = []
    if "valve_acquisitions" in payload:
        if not isinstance(payload["valve_acquisitions"], list):
            raise DomainError("Valve acquisitions must be an array.")
        valves = [_read_valve_acquisition(row) for row in payload["valve_acquisitions"]]
    return ArchivedTransientRun(payload["test_article_id"], tuple(channels), tuple(valves))


def transient_archive_from_bytes(data: bytes) -> ArchivedTransientRun:
    """Reopen a canonical acquisition export for independent replay.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, input identity and referent reproducibility.
    """
    if not isinstance(data, bytes):
        raise DomainError("An acquisition archive must be read as bytes.")
    try:
        source = data.decode("ascii")
    except UnicodeError as error:
        raise DomainError("An acquisition archive must use canonical ASCII JSON.") from error
    run = cast(ArchivedTransientRun, read_json_document(
        source, "Acquisition archive", _transient_run_from_payload
    ))
    if transient_archive_bytes(run) != data:
        raise DomainError("Acquisition archive is not in canonical form.")
    return run


@dataclass(frozen=True, slots=True)
class TransientValidationTarget:
    """Reviewed channel, window and statistic for a required holdout sample.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, metric selection and measurement uncertainty.
    """

    channel_id: str
    start_tick: int
    end_tick: int
    metric: Literal["peak", "signed_integral", "positive_impulse"]
    window_review_record: str
    pressure_referent: PressureReferent | None = None

    def __post_init__(self) -> None:
        """Validate the recipe without substituting a missing window endpoint.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1 and 4.2.5, reviewed scope of referents.
        """
        text(self.channel_id, "planned transient channel")
        text(self.window_review_record, "planned transient window review")
        start = integer(self.start_tick, "planned transient start")
        end = integer(self.end_tick, "planned transient end")
        if end <= start:
            raise DomainError("Planned transient end must be after its start.")
        if self.metric not in ("peak", "signed_integral", "positive_impulse"):
            raise DomainError("Unsupported transient validation metric.")
        if self.pressure_referent is not None and (
            self.metric == "positive_impulse" or not isinstance(self.pressure_referent, PressureReferent)
        ):
            raise DomainError("A pressure referent requires a pressure metric.")
        object.__setattr__(self, "start_tick", start)
        object.__setattr__(self, "end_tick", end)


@dataclass(frozen=True, slots=True)
class ShutdownValidationTarget:
    """Reviewed recipe for forward impulse after the final observed closure.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, metric scope and timing uncertainty.
    """

    channel_id: str
    close_commands: tuple[ValveCommand, ...]
    limits: ValveLatencyLimits
    window: ShutdownWindow

    def __post_init__(self) -> None:
        """Bind the selected closing events, approved limits and alignment record.

        References
        ----------
        .. [1] NASA-7009B, 4.1.1.5 and 4.2.1, limits and reviewed inputs.
        """
        text(self.channel_id, "planned shutdown force channel")
        commands: tuple[ValveCommand, ...] = records(
            self.close_commands, ValveCommand, "planned closing commands"
        )
        if (not commands or len({cmd.head_id for cmd in commands}) != len(commands)
                or any(cmd.opened for cmd in commands)):
            raise DomainError("A shutdown target requires one closing command per selected valve.")
        if not isinstance(self.limits, ValveLatencyLimits) or not isinstance(self.window, ShutdownWindow):
            raise DomainError("A shutdown target requires reviewed latency limits and window.")
        object.__setattr__(self, "close_commands", commands)


@dataclass(frozen=True, slots=True)
class TransientMeasurement:
    """An outward symmetric representation for a replayed holdout measurement.

    References
    ----------
    .. [1] NASA-7009B, 4.2.4 and 4.2.7, numerical and measurement uncertainty.
    """

    unit: str
    value: float
    absolute_error: float
    enclosure: tuple[float, float]
    final_closure_ticks: tuple[int, int] | None = None
    valve_timing_status: str | None = None


def transient_measurement(
    trace: TransientTrace, target: TransientValidationTarget,
) -> TransientMeasurement | str:
    """Replay the metric or report the reason it is unusable in a holdout.

    Use the rounded midpoint as the submitted measurement. Round the larger
    exact endpoint distance outward as the absolute uncertainty. A singleton
    uses the smallest positive binary64 radius because holdout measurements
    require positive absolute uncertainty.

    References
    ----------
    .. [1] NASA-7009B, 4.2.4 and 4.2.7, conservative treatment of uncertainty.
    """
    if not isinstance(trace, TransientTrace) or not isinstance(target, TransientValidationTarget):
        raise DomainError("A transient measurement requires a trace and a reviewed target.")
    if trace.channel_id != target.channel_id:
        return "transient_channel_mismatch"
    if target.metric == "positive_impulse" and trace.unit != "N":
        return "transient_metric_unit_mismatch"
    review = review_transient(trace, target.start_tick, target.end_tick,
                              window_review_record=target.window_review_record)
    if review.status != "bounded":
        return "transient_" + review.status
    slew = trace.max_absolute_slew
    assert slew is not None  # A bounded reduction requires a slew bound.
    if trace.tick_period.assumption_notes or slew.assumption_notes:
        return "transient_assumed_time_or_slew"
    unit = trace.unit if target.metric == "peak" else f"{trace.unit}*s"
    interval = {
        "peak": review.peak, "signed_integral": review.signed_integral,
        "positive_impulse": review.positive_impulse_ns,
    }[target.metric]
    assert interval is not None  # A bounded, dimensionally valid metric has an enclosure.
    return _centered_measurement(unit, interval)


def shutdown_valve_timing(
    target: ShutdownValidationTarget, valve: ArchivedValveAcquisition,
) -> ValveResponseReview | str:
    """Replay the approved valve timing check without a force measurement.

    A declared shutdown must identify real issued closes. The command and
    observed clocks must share a reviewed scale. Force-to-command alignment
    is relevant to impulse, while valve timing uses valve-to-command alignment.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, evidence and relative timing uncertainty.
    """
    if not isinstance(target, ShutdownValidationTarget) or not isinstance(valve, ArchivedValveAcquisition):
        raise DomainError("A valve timing replay requires a shutdown target and archived valve records.")
    if not set(target.close_commands).issubset(valve.commanded.commands):
        return "shutdown_close_not_archived"
    command_clock = valve.commanded.tick_period
    observation_clock = valve.observed.history.tick_period
    if command_clock.value != observation_clock.value:
        return "shutdown_valve_clock_scale_mismatch"
    if command_clock.assumption_notes or observation_clock.assumption_notes:
        return "shutdown_assumed_valve_clock"
    return review_valve_responses(valve.commanded, valve.observed, target.limits)


def shutdown_measurement(
    trace: TransientTrace, target: ShutdownValidationTarget, valve: ArchivedValveAcquisition,
) -> TransientMeasurement | str:
    """Replay a closure-relative measurement using the selected archived valve acquisition.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, experiment replay and timing.
    """
    if not isinstance(trace, TransientTrace) or not isinstance(target, ShutdownValidationTarget):
        raise DomainError("A shutdown measurement requires a trace and a reviewed shutdown target.")
    if not isinstance(valve, ArchivedValveAcquisition):
        raise DomainError("A shutdown measurement requires an archived valve acquisition.")
    if trace.channel_id != target.channel_id:
        return "transient_channel_mismatch"
    if trace.unit != "N":
        return "transient_metric_unit_mismatch"
    if not set(target.close_commands).issubset(valve.commanded.commands):
        return "shutdown_close_not_archived"
    if not trace.tick_period.value == valve.commanded.tick_period.value == valve.observed.history.tick_period.value:
        return "shutdown_clock_scale_mismatch"
    review = review_shutdown_impulse(
        valve.commanded, valve.observed, target.limits, trace, target.close_commands, target.window
    )
    if review.status != "bounded":
        return "shutdown_" + review.status
    slew = trace.max_absolute_slew
    assert slew is not None
    if any(item.assumption_notes for item in (
        trace.tick_period, slew, valve.commanded.tick_period, valve.observed.history.tick_period,
    )):
        return "transient_assumed_time_or_slew"
    assert review.thrust_review is not None and review.thrust_review.positive_impulse_ns is not None
    return _centered_measurement(
        "N*s", review.thrust_review.positive_impulse_ns,
        final_closure_ticks=review.final_closure_ticks, valve_timing_status=review.valve_review.status,
    )


def _centered_measurement(
    unit: str, interval: tuple[float, float], *,
    final_closure_ticks: tuple[int, int] | None = None, valve_timing_status: str | None = None,
) -> TransientMeasurement:
    """Represent both endpoints by a binary64 center and outward absolute radius.

    References
    ----------
    .. [1] NASA-7009B, 4.2.4 and 4.2.7, conservative treatment of numerical error.
    """
    lo, hi = map(Fraction, interval)
    center = float((lo + hi) / 2)
    radius = outward_endpoint(max(Fraction(center) - lo, hi - Fraction(center)), upward=True)
    return TransientMeasurement(
        unit, center, max(radius, ulp(0.0)), interval, final_closure_ticks, valve_timing_status
    )

