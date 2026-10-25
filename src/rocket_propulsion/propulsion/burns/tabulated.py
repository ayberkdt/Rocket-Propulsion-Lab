"""Target-neutral, discontinuity-safe tabulated burn exchange artifacts.

The artifact defined here is the exchange contract between propulsion models
and a downstream trajectory consumer.  It carries propulsion force and
tank-drain histories in relative time and contains no epoch, direction,
frame, orbit or propagator state.  Each segment stores its left and right
endpoint states independently, so an instantaneous command jump is neither
smoothed nor represented by duplicate timestamps.

Integrals of the piecewise-linear histories (impulse, consumed mass and the
thrust first moment) are accumulated in exact rational arithmetic and
rounded once on export.  Total impulse and equivalent specific impulse follow
the definitions ``I = integral F dt`` and ``Isp = I / (g0 * m)`` [1]_, with
``g0`` the registered ``STANDARD_GRAVITY_M_S2``.  [2]_ gives the context for
versioned, unit-explicit exchange documents.

References
----------
.. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
.. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
       https://ccsds.org/Pubs/502x0b3e1.pdf
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from functools import lru_cache
from itertools import pairwise
from math import isfinite
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    finite_float_text,
    name_value_pairs,
    real,
    records,
    strings,
    text,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .exchange_documents import _read_document
from .models import ProfiledBurnResult
from .numeric_outputs import _finite_rational, _fraction
from .profiles import resolved_schedule_intervals
from .serialization import burn_result_to_dict, canonical_json_bytes
from .tabulated_documents import _utf8, _validate_document
from .time_coordinates import _TIME_POLICY

_TABULATED_LEDGER_NUMERICS = "exact_binary_segment_aggregate_ledger_v1"


@lru_cache(maxsize=1024)
def _segment_ledger(segment):
    """Return the exact rational integrals of one linear segment.

    For a quantity varying linearly from ``q0`` to ``q1`` over a segment of
    duration ``h`` starting at ``t0``, the integral is ``h * (q0 + q1) / 2``
    and the first moment of thrust about time zero is
    ``t0 * I + h**2 * (F0 + 2 * F1) / 6``.  Total impulse is the time integral
    of thrust [1]_.

    Returns
    -------
    tuple of Fraction
        ``(duration, delivered_impulse, axial_impulse, propellant_mass,
        thrust_first_moment)``, all unrounded.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
    duration = _fraction(segment.end_time_s) - _fraction(segment.start_time_s)
    f0, f1 = _fraction(segment.start_delivered_thrust_n), _fraction(segment.end_delivered_thrust_n)
    impulse = duration * (f0 + f1) / 2
    axial = (
        duration
        * (_fraction(segment.start_axial_thrust_n) + _fraction(segment.end_axial_thrust_n))
        / 2
    )
    mass = (
        duration
        * (
            _fraction(segment.start_total_mass_flow_kg_s)
            + _fraction(segment.end_total_mass_flow_kg_s)
        )
        / 2
    )
    moment = _fraction(segment.start_time_s) * impulse + duration * duration * (f0 + 2 * f1) / 6
    return duration, impulse, axial, mass, moment


@lru_cache(maxsize=256)
def _artifact_ledger(segments):
    """Sum the exact segment ledgers before any public scalar is rounded.

    Returns
    -------
    tuple of Fraction
        ``(delivered_impulse, axial_impulse, propellant_mass,
        thrust_first_moment)`` over all segments.
    """
    rows = tuple(_segment_ledger(segment) for segment in segments)
    return tuple(sum((row[index] for row in rows), Fraction(0)) for index in range(1, 5))


_ARTIFACT_SCALE_NAMES = {
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
}


class InterpolationPolicy(str, Enum):
    """Interpolation applied inside each declared exchange segment.

    ``LINEAR`` varies every state linearly between the segment endpoints;
    ``HOLD`` requires equal endpoint states, so each segment is constant.
    The policy is declared explicitly in the exchange document, in keeping
    with the unit- and version-explicit exchange context of [1]_.

    References
    ----------
    .. [1] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
           https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    LINEAR = "linear"
    HOLD = "hold"


def _nonnegative(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and non-negative."""
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _normalized_pairs(value: object, label: str) -> tuple[tuple[str, float], ...]:
    """Coerce ``value`` to a tuple of ``(name, real)`` pairs."""
    return name_value_pairs(value, label)


def _flow_pairs(value: tuple[tuple[str, float], ...], label: str) -> None:
    """Require non-blank, unique, UTF-8 identifiers and non-negative flows."""
    identifiers = [identifier for identifier, _ in value]
    if any(not identifier.strip() for identifier in identifiers):
        raise DomainError(f"{label} identifiers must not be empty.")
    if len(identifiers) != len(set(identifiers)):
        raise DomainError(f"{label} identifiers must be unique.")
    for _, flow in value:
        _nonnegative(flow, label)
    for identifier, _ in value:
        _utf8(identifier, label)


@dataclass(frozen=True, slots=True)
class TabulatedBurnSegment:
    """One time segment with independent left and right propulsion states.

    Parameters
    ----------
    start_time_s, end_time_s : float
        Segment limits in relative burn time, s; ``0 <= start < end``.
    start_delivered_thrust_n, end_delivered_thrust_n : float
        Non-negative delivered thrust magnitude at each endpoint, N.
    start_axial_thrust_n, end_axial_thrust_n : float
        Non-negative axial thrust at each endpoint, N; at most the delivered
        thrust (relative allowance ``1e-12``).
    start_total_mass_flow_kg_s, end_total_mass_flow_kg_s : float
        Non-negative total tank outflow at each endpoint, kg/s.
    start_stream_mass_flows_kg_s, end_stream_mass_flows_kg_s : tuple of (str, float)
        Optional named tank drains.  When present they must sum to the total
        flow within ``max(1e-12, 1e-10 * total)`` kg/s.
    phase : str, default "tabulated"
        Non-blank, UTF-8 encodable phase label.

    Raises
    ------
    DomainError
        If a value is non-finite or outside its domain, the times are not
        ordered, axial thrust exceeds delivered thrust, a stream identifier is
        blank or repeated, or the stream flows do not close to the total.

    Notes
    -----
    Every state varies linearly inside the segment.  Integrals are exact
    trapezoids of the endpoint states, and total impulse is the time integral
    of thrust [1]_.  [2]_ and [3]_ give the engine transient-modelling
    context for such force and flow histories.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    .. [3] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
        for name in (
            "start_time_s",
            "end_time_s",
            "start_delivered_thrust_n",
            "end_delivered_thrust_n",
            "start_axial_thrust_n",
            "end_axial_thrust_n",
            "start_total_mass_flow_kg_s",
            "end_total_mass_flow_kg_s",
        ):
            object.__setattr__(self, name, real(getattr(self, name), f"Tabulated segment {name}"))
        for name in ("start_stream_mass_flows_kg_s", "end_stream_mass_flows_kg_s"):
            object.__setattr__(
                self, name, _normalized_pairs(getattr(self, name), f"Tabulated segment {name}")
            )
        object.__setattr__(self, 'phase', _utf8(self.phase, "Tabulated segment phase"))
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
        if _fraction(self.start_axial_thrust_n) > _fraction(
            self.start_delivered_thrust_n
        ) * _fraction(1.0 + 1e-12):
            raise DomainError("Start axial thrust cannot exceed delivered thrust.")
        if _fraction(self.end_axial_thrust_n) > _fraction(self.end_delivered_thrust_n) * _fraction(
            1.0 + 1e-12
        ):
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
        """Require stream flows, when given, to sum to ``total`` within tolerance."""
        if not streams:
            return
        actual = sum((_fraction(flow) for _, flow in streams), Fraction(0))
        if abs(actual - _fraction(total)) > max(
            _fraction(1e-12), _fraction(total) * _fraction(1e-10)
        ):
            raise DomainError(f"Tabulated {side} stream flows do not close to total flow.")

    @property
    def duration_s(self) -> float:
        """Segment duration ``end_time_s - start_time_s``, s."""
        return self.end_time_s - self.start_time_s

    @staticmethod
    def _integral(start: float, end: float, duration: float) -> float:
        """Return the trapezoid integral ``h * (q0 + q1) / 2`` of a linear quantity."""
        return _finite_rational(
            (_fraction(start) + _fraction(end)) * _fraction(duration) / 2, "Tabulated integral"
        )

    @staticmethod
    def _local_first_moment(start: float, end: float, duration: float) -> float:
        """Return ``h**2 * (q0 + 2 * q1) / 6``, the first moment about the segment start."""
        return _finite_rational(
            _fraction(duration) ** 2 * (_fraction(start) + 2 * _fraction(end)) / 6,
            "Tabulated local first moment",
        )

    @property
    def delivered_impulse_n_s(self) -> float:
        """Delivered impulse ``integral F dt`` over the segment, N s [1]_.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        return _finite_rational(_segment_ledger(self)[1], "Tabulated delivered impulse")

    @property
    def axial_impulse_n_s(self) -> float:
        """Axial impulse, the time integral of axial thrust over the segment, N s [1]_.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        return _finite_rational(_segment_ledger(self)[2], "Tabulated axial impulse")

    @property
    def propellant_mass_kg(self) -> float:
        """Propellant drained over the segment, the integral of total flow, kg."""
        return _finite_rational(_segment_ledger(self)[3], "Tabulated propellant mass")

    @property
    def delivered_thrust_first_moment_n_s2(self) -> float:
        """First moment ``integral t F dt`` of delivered thrust about burn start, N s^2."""
        return _finite_rational(_segment_ledger(self)[4], "Tabulated global first moment")

    def _fraction_at(self, time_s):
        """Return the exact normalized position of ``time_s`` inside the segment.

        Raises `DomainError` if the time lies outside ``[start, end]``.
        """
        time_s = _fraction(time_s, "Tabulated segment time")
        if not self.start_time_s <= time_s <= self.end_time_s:
            raise DomainError("Tabulated evaluation time lies outside the segment.")
        return (_fraction(time_s) - _fraction(self.start_time_s)) / _segment_ledger(self)[0]

    def _axial_at(self, time_s):
        """Return the exact linearly interpolated axial thrust at ``time_s``."""
        fraction = self._fraction_at(time_s)
        return (1 - fraction) * _fraction(self.start_axial_thrust_n) + fraction * _fraction(
            self.end_axial_thrust_n
        )

    def state_at(self, time_s: float) -> dict[str, Any]:
        """Return the linearly interpolated propulsion state at ``time_s``.

        Parameters
        ----------
        time_s : float
            Relative burn time inside ``[start_time_s, end_time_s]``, s.

        Returns
        -------
        dict
            ``time_s``, ``delivered_thrust_n``, ``axial_thrust_n``,
            ``total_mass_flow_kg_s``, ``stream_mass_flows_kg_s`` (sorted by
            identifier; a stream absent at one endpoint counts as zero there)
            and ``phase``.  Every interpolated value is checked to be finite.

        Raises
        ------
        DomainError
            If ``time_s`` is invalid or lies outside the segment.
        """
        fraction = self._fraction_at(time_s)
        time_s = real(time_s, "Tabulated segment time")

        def interpolate(start: float, end: float) -> float:
            return _finite_rational(
                (1 - fraction) * _fraction(start) + fraction * _fraction(end),
                "Tabulated interpolated state",
            )

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


def epoch_time_tolerance(*times_s: float) -> float:
    """Return the rounding tolerance for a time formed from absolute epochs.

    Parameters
    ----------
    *times_s : float
        Epoch magnitudes, s, from which the relative time is formed.

    Returns
    -------
    float
        ``16 * eps * max(1, |t_1|, ..., |t_n|)``, s, with binary64 machine
        epsilon ``eps = 2**-52``.

    Raises
    ------
    DomainError
        If any epoch is not a finite real number.

    Notes
    -----
    A time built by adding or subtracting a few binary64 epochs of magnitude
    ``|t|`` carries an error of a few units in the last place [1]_, well
    within ``16 * eps * |t|`` (about 3.6e-8 s at 1e7 s and 3.6e-15 s near
    zero).

    References
    ----------
    .. [1] D. Goldberg, "What Every Computer Scientist Should Know About
           Floating-Point Arithmetic", ACM Computing Surveys 23(1), 1991.
    """

    scale = 1.0
    for value in times_s:
        scale = max(scale, abs(real(value, "Epoch")))
    return 16.0 * 2.220446049250313e-16 * scale


def _segment_start(segment: TabulatedBurnSegment) -> float:
    """Bisection key: the segment start time."""
    return segment.start_time_s


def _segment_end(segment: TabulatedBurnSegment) -> float:
    """Bisection key: the segment end time."""
    return segment.end_time_s


def _canonical_json_bytes(value: Any) -> bytes:
    """Encode ``value`` as sorted, compact, NaN-free UTF-8 JSON for hashing."""
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


@dataclass(frozen=True, slots=True)
class TabulatedBurnArtifact:
    """Immutable, hash-bound propulsion history for a dynamics consumer.

    Parameters
    ----------
    segments : tuple of TabulatedBurnSegment
        Non-empty, exactly contiguous segments starting at relative time zero.
    source_id : str
        Non-blank identifier of the producing model or data set.
    source_sha256 : str
        SHA-256 hex digest of the producer's source content.
    interpolation : InterpolationPolicy, default InterpolationPolicy.LINEAR
        Interpolation inside each segment; ``HOLD`` requires constant
        segments.
    reference_ids : tuple of str, default ()
        Unique, non-blank source identifiers of the producer model.
    warnings : tuple of str, default ()
        Model-scope warnings carried with the history.
    schema : str, default "rocket_propulsion_tabulated_burn_v1"
        Document schema; no other value is accepted.

    Raises
    ------
    DomainError
        If the schema or interpolation policy is unsupported, the segments
        are empty, do not start at zero or are not contiguous, a ``HOLD``
        segment varies, the source identifier is blank, the source hash is not
        a SHA-256 digest, or a reference identifier is blank or repeated.

    Notes
    -----
    ``reference_ids`` records caller-supplied source identifiers and is
    included in the artifact hash.  Its presence does not establish source
    review or physical model validity.  Derived totals follow the impulse and
    specific-impulse definitions of [1]_; [2]_ gives the context for
    versioned, unit-explicit exchange documents and [3]_ the
    transient-modelling context.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
           https://ccsds.org/Pubs/502x0b3e1.pdf
    .. [3] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    """

    segments: tuple[TabulatedBurnSegment, ...]
    source_id: str
    source_sha256: str
    interpolation: InterpolationPolicy = InterpolationPolicy.LINEAR
    reference_ids: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    schema: str = "rocket_propulsion_tabulated_burn_v1"

    def __post_init__(self) -> None:
        object.__setattr__(self, 'schema', _utf8(self.schema, 'Tabulated schema'))
        if self.schema != "rocket_propulsion_tabulated_burn_v1":
            raise DomainError("Unsupported tabulated burn schema.")
        object.__setattr__(
            self, "segments", records(self.segments, TabulatedBurnSegment, "Tabulated segments")
        )
        try:
            object.__setattr__(self, "interpolation", InterpolationPolicy(
                _utf8(self.interpolation, 'Tabulated interpolation')))
        except ValueError as error:
            raise DomainError("Unknown tabulated interpolation policy.") from error
        object.__setattr__(self, 'source_id', _utf8(self.source_id, "Tabulated burn source identifier"))
        object.__setattr__(self, 'source_sha256', _utf8(self.source_sha256, "Tabulated burn source hash"))
        for name in ("reference_ids", "warnings"):
            object.__setattr__(
                self, name, strings(getattr(self, name), f"Tabulated burn {name}", allow_blank=True)
            )
            for item in getattr(self, name):
                _utf8(item, name, allow_blank=True)
        if not self.segments:
            raise DomainError("A tabulated burn requires at least one segment.")
        if self.segments[0].start_time_s != 0.0:
            raise DomainError("A tabulated burn must begin at zero relative time.")
        for previous, current in zip(self.segments, self.segments[1:], strict=False):
            if previous.end_time_s != current.start_time_s:
                raise DomainError("Tabulated burn segments must be contiguous.")
        if self.interpolation is InterpolationPolicy.HOLD:
            for segment in self.segments:
                if any(
                    start != end
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
        """Burn duration: the end time of the last segment, s."""
        return self.segments[-1].end_time_s

    @property
    def delivered_total_impulse_n_s(self) -> float:
        """Delivered total impulse ``integral F dt`` over the burn, N s [1]_.

        The exact segment integrals are summed before one final rounding.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        return _finite_rational(
            _artifact_ledger(self.segments)[0], "Tabulated delivered total impulse"
        )

    @property
    def axial_total_impulse_n_s(self) -> float:
        """Axial total impulse, the burn integral of axial thrust, N s [1]_.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        return _finite_rational(_artifact_ledger(self.segments)[1], "Tabulated axial total impulse")

    @property
    def consumed_propellant_kg(self) -> float:
        """Propellant consumed over the burn, the integral of total tank flow, kg."""
        return _finite_rational(_artifact_ledger(self.segments)[2], "Tabulated consumed propellant")

    @property
    def thrust_centroid_time_s(self) -> float:
        """Impulse-weighted mean time ``integral t F dt / integral F dt``, s.

        Zero when the delivered impulse is zero.
        """
        impulse, _, _, moment = _artifact_ledger(self.segments)
        return _finite_rational(moment / impulse, "Tabulated thrust centroid") if impulse else 0.0

    @property
    def system_equivalent_specific_impulse_s(self) -> float | None:
        """System-equivalent specific impulse ``I / (g0 * m)``, s, or ``None``.

        ``I`` is the delivered total impulse, ``m`` the consumed propellant and
        ``g0`` the registered ``STANDARD_GRAVITY_M_S2`` [1]_.  ``None`` is
        returned when no propellant is consumed.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        impulse, _, mass, _ = _artifact_ledger(self.segments)
        return (
            _finite_rational(
                impulse / (_fraction(STANDARD_GRAVITY_M_S2) * mass),
                "Tabulated equivalent specific impulse",
            )
            if mass
            else None
        )

    @property
    def artifact_hash(self) -> str:
        """SHA-256 hex digest of the canonical document without its hash field.

        The digest identifies content; it is not authentication.
        """
        return hashlib.sha256(_canonical_json_bytes(self.to_dict(include_hash=False))).hexdigest()

    def state_at(self, time_s: float, *, side: str = "right") -> dict[str, Any]:
        """Evaluate the propulsion state with explicit discontinuity semantics.

        Parameters
        ----------
        time_s : float
            Relative burn time in ``[0, duration_s]``, s.
        side : {"right", "left"}, default "right"
            Limit taken at an interior segment boundary; ``"left"`` returns the
            state of the segment ending there, ``"right"`` the state of the
            segment starting there.

        Returns
        -------
        dict
            The owning segment's state, as returned by
            `TabulatedBurnSegment.state_at`.

        Raises
        ------
        DomainError
            If ``side`` is invalid or ``time_s`` lies outside the artifact.
        """

        if side not in ("left", "right"):
            raise DomainError("Tabulated boundary side must be 'left' or 'right'.")
        time_s = real(time_s, "Tabulated burn time")
        if not 0.0 <= time_s <= self.duration_s:
            raise DomainError("Tabulated burn time lies outside the artifact.")
        segment = self.segments[self.segment_index_at(time_s, side=side)]
        return segment.state_at(min(max(time_s, segment.start_time_s), segment.end_time_s))

    @property
    def boundary_tolerance_s(self) -> float:
        """Default boundary tolerance, s: zero, as boundaries are exact binary values."""

        return 0.0

    def segment_index_at(self, time_s: float, *, side: str = "right") -> int:
        """Return the index of the segment that owns ``time_s``, in O(log n).

        Parameters
        ----------
        time_s : float
            Relative burn time, s.
        side : {"right", "left"}, default "right"
            Ownership of a time exactly on an interior segment start.

        Returns
        -------
        int
            Index of the owning segment.

        Raises
        ------
        DomainError
            If ``side`` is invalid or ``time_s`` is not a finite real number.

        Notes
        -----
        The owning segment satisfies ``start <= t < end``.  With
        ``side="left"`` a time exactly equal to an interior segment start
        belongs to the preceding segment, so a jump keeps its left limit.  A
        time at or beyond the final end belongs to the last segment.  Times
        before zero map to the first segment; the caller decides whether such
        a time is inside the firing interval.  Segments are exactly
        contiguous, so every representable interior time has one owner.
        [1]_ gives the engine transient-modelling context for start, steady
        and shutdown histories with such phase boundaries.

        References
        ----------
        .. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
               https://ntrs.nasa.gov/citations/19970010379
        """

        if side not in ("left", "right"):
            raise DomainError("Tabulated boundary side must be 'left' or 'right'.")
        time_s = real(time_s, "Tabulated burn time")
        return self._segment_index_at(_fraction(time_s), side=side)

    def _segment_index_at(self, time_s: Fraction, *, side: str) -> int:
        """Return the owning segment index for an exact rational time."""
        segments = self.segments
        probe: Fraction | float = time_s  # exact time bisected against float keys
        index = max(0, bisect_right(segments, probe, key=_segment_start) - 1)
        if (
            side == "left"
            and index > 0
            and time_s < segments[index].end_time_s
            and time_s == segments[index].start_time_s
        ):
            return index - 1
        return index

    @property
    def shortest_segment_s(self) -> float:
        """Duration of the shortest segment, s."""

        return min(segment.duration_s for segment in self.segments)

    def snap_to_boundary(self, time_s: float, tolerance_s: float | None = None) -> float:
        """Snap a time to the nearest boundary within a tolerance.

        Parameters
        ----------
        time_s : float
            Relative burn time, s.
        tolerance_s : float or None, default None
            Non-negative snapping window, s, strictly below half the shortest
            segment duration.  ``None`` means zero: only exact boundaries snap.

        Returns
        -------
        float
            The nearest boundary (zero or a segment end) if it lies within
            ``tolerance_s`` of ``time_s``, otherwise ``time_s`` unchanged.

        Raises
        ------
        DomainError
            If ``time_s`` or ``tolerance_s`` is invalid, ``tolerance_s`` is
            negative, or it is not below half the shortest segment duration.

        Notes
        -----
        There is no implicit tolerance floor.  Distances and the window check
        are evaluated in exact rational arithmetic before any binary rounding
        [1]_.  A consumer working in absolute epochs uses the propagation
        bridge's epoch entry point, which names clock events explicitly in the
        manner of a propagator event detector [2]_, rather than a
        scale-dependent tolerance.

        References
        ----------
        .. [1] D. Goldberg, "What Every Computer Scientist Should Know About
               Floating-Point Arithmetic", ACM Computing Surveys 23(1), 1991.
        .. [2] Orekit 13.1.5, "PropulsionModel / EventDetector" API
               documentation, ``org.orekit.forces.maneuvers.propulsion``
               (URL listed in docs/REFERENCES.md).
        """

        time_s = real(time_s, "Tabulated burn time")
        return _finite_rational(self._snap_to_boundary(_fraction(time_s), tolerance_s),
                                'Snapped profile time')

    def _snap_to_boundary(self, time_s: Fraction, tolerance_s: float | None) -> Fraction:
        """Snap an exact rational time; see `snap_to_boundary` for the contract."""
        tolerance = Fraction(0)
        if tolerance_s is not None:
            requested = real(tolerance_s, "Boundary tolerance")
            if requested < 0.0:
                raise DomainError("Boundary tolerance must be non-negative.")
            shortest = min(_segment_ledger(segment)[0] for segment in self.segments)
            if 2 * _fraction(requested) >= shortest:
                raise DomainError(
                    "Boundary tolerance must be below half the shortest segment duration."
                )
            tolerance = _fraction(requested)
        probe: Fraction | float = time_s  # exact time bisected against float keys
        index = bisect_left(self.segments, probe, key=_segment_end)
        candidates = [0.0]
        if index < len(self.segments):
            candidates.append(self.segments[index].end_time_s)
        if index > 0:
            candidates.append(self.segments[index-1].end_time_s)
        exact_time = _fraction(time_s)
        nearest = min(candidates, key=lambda value: abs(_fraction(value)-exact_time))
        if abs(_fraction(nearest)-exact_time) <= tolerance:
            return _fraction(nearest)
        return time_s

    def is_boundary_time(self, time_s: float) -> bool:
        """Return whether ``time_s`` equals zero or a segment end, in O(log n).

        Parameters
        ----------
        time_s : float
            Relative burn time, s.

        Returns
        -------
        bool
            ``True`` only for an exact boundary; no tolerance is applied.

        Raises
        ------
        DomainError
            If ``time_s`` is not a finite real number.
        """

        time_s = real(time_s, "Tabulated burn time")
        return self._is_boundary_time(_fraction(time_s))

    def _is_boundary_time(self, time_s: Fraction) -> bool:
        """Return whether an exact rational time equals zero or a segment end."""
        if time_s == 0.0:
            return True
        probe: Fraction | float = time_s  # exact time bisected against float keys
        index = bisect_left(self.segments, probe, key=_segment_end)
        return index < len(self.segments) and self.segments[index].end_time_s == time_s

    def _history_dict(self) -> dict[str, Any]:
        """Return the source samples and metadata, without derived summary fields."""

        payload: dict[str, Any] = {
            "schema": self.schema,
            "integral_numerics": _TABULATED_LEDGER_NUMERICS,
            "time_numerics": _TIME_POLICY,
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
            "warnings": list(self.warnings),
        }
        if self.reference_ids:
            payload["reference_ids"] = list(self.reference_ids)
        return payload

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the JSON-compatible artifact document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add ``artifact_hash``, the SHA-256 digest of the
            canonical document without that field.

        Returns
        -------
        dict
            Source samples and metadata, the numerical-policy identifiers, and
            a ``summary`` of the derived totals, each rounded once from its
            exact value.
        """
        payload = self._history_dict()
        payload["summary"] = {
            "duration_s": self.duration_s,
            "delivered_total_impulse_n_s": self.delivered_total_impulse_n_s,
            "axial_total_impulse_n_s": self.axial_total_impulse_n_s,
            "consumed_propellant_kg": self.consumed_propellant_kg,
            "thrust_centroid_time_s": self.thrust_centroid_time_s,
            "system_equivalent_specific_impulse_s": self.system_equivalent_specific_impulse_s,
        }
        if include_hash:
            payload["artifact_hash"] = hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize `to_dict` as sorted JSON without NaN or Infinity tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        """
        return json.dumps(
            self.to_dict(), ensure_ascii=False, allow_nan=False, sort_keys=True, indent=indent
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> TabulatedBurnArtifact:
        """Reconstruct an artifact and verify its hash and derived summary.

        Parameters
        ----------
        payload : dict
            Document produced by `to_dict`.

        Returns
        -------
        TabulatedBurnArtifact
            The reconstructed artifact.

        Raises
        ------
        InputError
            If the document structure is invalid, a field is missing or
            unsupported, a numerical-policy identifier differs, a field value
            violates the artifact or segment contract, the content hash does
            not match, or the summary does not reproduce.
        DomainError
            If a text field of the document is blank or not UTF-8 encodable.
        """

        if not isinstance(payload, dict):
            raise InputError("Tabulated burn artifact root must be an object.")
        _validate_document(payload)
        if payload.get("integral_numerics") != _TABULATED_LEDGER_NUMERICS:
            raise InputError(
                "Tabulated integral numerical policy is missing or unsupported.",
                field="integral_numerics",
            )
        if payload.get("time_numerics") != _TIME_POLICY:
            raise InputError("Tabulated time numerical policy is missing or unsupported.",
                             field="time_numerics")
        try:
            artifact = cls(
                schema=payload["schema"],
                source_id=payload["source_id"],
                source_sha256=payload["source_sha256"],
                interpolation=payload["interpolation"],
                reference_ids=payload.get("reference_ids", ()),
                segments=tuple(
                    TabulatedBurnSegment(
                        start_time_s=item["start_time_s"],
                        end_time_s=item["end_time_s"],
                        start_delivered_thrust_n=item["start_delivered_thrust_n"],
                        end_delivered_thrust_n=item["end_delivered_thrust_n"],
                        start_axial_thrust_n=item["start_axial_thrust_n"],
                        end_axial_thrust_n=item["end_axial_thrust_n"],
                        start_total_mass_flow_kg_s=item["start_total_mass_flow_kg_s"],
                        end_total_mass_flow_kg_s=item["end_total_mass_flow_kg_s"],
                        start_stream_mass_flows_kg_s=item["start_stream_mass_flows_kg_s"],
                        end_stream_mass_flows_kg_s=item["end_stream_mass_flows_kg_s"],
                        phase=item["phase"],
                    )
                    for item in payload["segments"]
                ),
                warnings=payload.get("warnings", ()),
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
        """Parse a strict JSON document [1]_ and reconstruct it with `from_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        """
        return _read_document(text, "Tabulated burn artifact", cls.from_dict)


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

    Parameters
    ----------
    artifact : TabulatedBurnArtifact
        Base history.
    scales : dict of str to float
        Positive finite factors keyed by ``thrust_scale``,
        ``specific_impulse_scale``, ``time_scale`` and
        ``axial_efficiency_scale``; an omitted key means 1.
    source_id : str
        Non-blank identifier of the scaled artifact.
    provenance : object
        JSON-serializable description of why the scales apply; bound into the
        source hash.
    reference_ids : tuple of str, default ()
        Additional source identifiers, appended without duplicates.
    warnings : tuple of str, default ()
        Additional warnings, appended to the base warnings.

    Returns
    -------
    TabulatedBurnArtifact
        Scaled artifact with the base interpolation policy.

    Raises
    ------
    DomainError
        If the inputs have the wrong type, a scale name is unknown, a factor
        is not positive and finite, the source identifier is blank, or the
        axial-efficiency factor would make axial thrust exceed delivered
        thrust.

    Notes
    -----
    Delivered thrust scales with ``thrust_scale``.  Because mass flow is
    ``F / (g0 * Isp)`` [1]_, total and named stream flows scale with
    ``thrust_scale / specific_impulse_scale``.  ``time_scale`` multiplies
    every segment boundary, and ``axial_efficiency_scale`` multiplies axial
    thrust without changing delivered thrust or tank drain.  The source hash
    binds the base history, the normalized scales and the provenance.  [2]_
    gives the model-credibility context for such derived models.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    if not isinstance(artifact, TabulatedBurnArtifact) or not isinstance(scales, dict):
        raise DomainError(
            "Artifact scaling requires an artifact and a scale mapping.", field="Artifact scales"
        )
    if any(not isinstance(name, str) for name in scales):
        raise DomainError("Artifact scale names must be strings.", field="Artifact scales")
    unknown = set(scales) - _ARTIFACT_SCALE_NAMES
    if unknown:
        raise DomainError(f"Unsupported artifact scales: {', '.join(sorted(unknown))}.")
    normalized = {name: 1.0 for name in _ARTIFACT_SCALE_NAMES}
    for name, value in scales.items():
        numeric = real(value, f"Artifact scale {name}")
        if not isfinite(numeric) or numeric <= 0.0:
            raise DomainError(f"Artifact scale {name!r} must be positive and finite.")
        normalized[name] = numeric
    text(source_id, "Scaled artifact source identifier")
    thrust_scale = normalized["thrust_scale"]
    isp_scale = normalized["specific_impulse_scale"]
    time_scale = normalized["time_scale"]
    axial_scale = normalized["axial_efficiency_scale"]
    maximum_axial_ratio = max(
        (
            _fraction(axial) / _fraction(delivered)
            for segment in artifact.segments
            for axial, delivered in (
                (segment.start_axial_thrust_n, segment.start_delivered_thrust_n),
                (segment.end_axial_thrust_n, segment.end_delivered_thrust_n),
            )
            if delivered > 0.0
        ),
        default=0.0,
    )
    if maximum_axial_ratio * _fraction(axial_scale) > _fraction(1.0 + 1e-12):
        raise DomainError("Artifact scale can make axial thrust exceed delivered thrust.")
    flow_scale = _fraction(thrust_scale) / _fraction(isp_scale)

    def scaled(value, factor, label):
        """Multiply exactly and round once, checking finiteness."""
        return _finite_rational(_fraction(value) * _fraction(factor), label)

    def scaled_streams(
        streams: tuple[tuple[str, float], ...],
    ) -> tuple[tuple[str, float], ...]:
        return tuple(
            (tank_id, scaled(flow, flow_scale, f"Scaled stream {tank_id}"))
            for tank_id, flow in streams
        )

    segments = tuple(
        TabulatedBurnSegment(
            start_time_s=scaled(segment.start_time_s, time_scale, "Scaled start time"),
            end_time_s=scaled(segment.end_time_s, time_scale, "Scaled end time"),
            start_delivered_thrust_n=scaled(
                segment.start_delivered_thrust_n, thrust_scale, "Scaled delivered thrust"
            ),
            end_delivered_thrust_n=scaled(
                segment.end_delivered_thrust_n, thrust_scale, "Scaled delivered thrust"
            ),
            start_axial_thrust_n=scaled(
                segment.start_axial_thrust_n,
                _fraction(thrust_scale) * _fraction(axial_scale),
                "Scaled axial thrust",
            ),
            end_axial_thrust_n=scaled(
                segment.end_axial_thrust_n,
                _fraction(thrust_scale) * _fraction(axial_scale),
                "Scaled axial thrust",
            ),
            start_total_mass_flow_kg_s=scaled(
                segment.start_total_mass_flow_kg_s, flow_scale, "Scaled total flow"
            ),
            end_total_mass_flow_kg_s=scaled(
                segment.end_total_mass_flow_kg_s, flow_scale, "Scaled total flow"
            ),
            start_stream_mass_flows_kg_s=scaled_streams(segment.start_stream_mass_flows_kg_s),
            end_stream_mass_flows_kg_s=scaled_streams(segment.end_stream_mass_flows_kg_s),
            phase=segment.phase,
        )
        for segment in artifact.segments
    )
    source_sha256 = hashlib.sha256(
        canonical_json_bytes(
            {
                "base_history_sha256": hashlib.sha256(
                    canonical_json_bytes(artifact._history_dict())
                ).hexdigest(),
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
    """Convert an analytic L1 profiled burn without smoothing discontinuities.

    Parameters
    ----------
    result : ProfiledBurnResult
        Constant-performance burn with a throttle schedule.

    Returns
    -------
    TabulatedBurnArtifact
        One linear segment per resolved schedule interval, with source hash
        ``result.result_hash``.

    Raises
    ------
    DomainError
        If ``result`` is not a `ProfiledBurnResult`, it does not serialize, or
        the artifact's delivered impulse or propellant mass differs from the
        analytic summary by more than ``max(1e-10, 1e-11 * |expected|)``.

    Notes
    -----
    Delivered thrust, axial thrust (delivered thrust times the cant
    efficiency) and every tank-depleting flow scale linearly with the realized
    throttle at each interval endpoint, so specific impulse is constant over
    the burn [1]_, [2]_.  Each schedule interval keeps its own endpoint
    values, so throttle steps remain discontinuities.

    References
    ----------
    .. [1] NASA SP-125. NASA, "Design of Liquid Propellant Rocket Engines".
           https://ntrs.nasa.gov/citations/19710019929
    .. [2] NASA Glenn Research Center, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    """

    if not isinstance(result, ProfiledBurnResult):
        raise DomainError("L1 conversion requires a ProfiledBurnResult.", field="result")
    burn_result_to_dict(result)
    tank_rates: dict[str, Fraction] = {}
    for stream in result.operating_point.streams:
        if stream.tank_depleting:
            tank_rates[stream.tank_id] = tank_rates.get(stream.tank_id, Fraction(0)) + _fraction(
                stream.mass_flow_kg_s
            )

    def flows(throttle: float) -> tuple[tuple[str, float], ...]:
        """Return per-tank drains at a throttle, sorted by tank identifier."""
        return tuple(
            (tank_id, _finite_rational(rate * _fraction(throttle), f"Profile stream {tank_id}"))
            for tank_id, rate in sorted(tank_rates.items())
        )

    point = result.operating_point
    cant = result.definition.cant_efficiency
    segments = tuple(
        TabulatedBurnSegment(
            start_time_s=interval.start_time_s,
            end_time_s=interval.end_time_s,
            start_delivered_thrust_n=_finite_rational(
                _fraction(point.delivered_thrust_n) * _fraction(interval.start_throttle),
                "Profile delivered thrust",
            ),
            end_delivered_thrust_n=_finite_rational(
                _fraction(point.delivered_thrust_n) * _fraction(interval.end_throttle),
                "Profile delivered thrust",
            ),
            start_axial_thrust_n=_finite_rational(
                _fraction(point.delivered_thrust_n)
                * _fraction(cant)
                * _fraction(interval.start_throttle),
                "Profile axial thrust",
            ),
            end_axial_thrust_n=_finite_rational(
                _fraction(point.delivered_thrust_n)
                * _fraction(cant)
                * _fraction(interval.end_throttle),
                "Profile axial thrust",
            ),
            start_total_mass_flow_kg_s=_finite_rational(
                _fraction(point.total_tank_flow_kg_s) * _fraction(interval.start_throttle),
                "Profile total flow",
            ),
            end_total_mass_flow_kg_s=_finite_rational(
                _fraction(point.total_tank_flow_kg_s) * _fraction(interval.end_throttle),
                "Profile total flow",
            ),
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


def _csv_point(row, stream_columns):
    """Parse and validate one CSV row, including the terminal row of a held trace.

    Returns
    -------
    dict
        ``time``, ``delivered``, ``axial`` (delivered thrust when the column
        is absent), ``flow``, ``streams`` and ``phase`` (``"imported-trace"``
        when the column is absent).

    Raises
    ------
    InputError
        If a numeric cell is not finite decimal text or axial thrust exceeds
        delivered thrust.
    DomainError
        If a value is negative, a stream identifier is invalid, or the stream
        flows do not close to the total flow.
    """
    def number(name):
        return finite_float_text(row[name], f'CSV {name}')

    delivered = number('delivered_thrust_n')
    axial = number('axial_thrust_n') if 'axial_thrust_n' in row else delivered
    flow = number('total_mass_flow_kg_s')
    streams = tuple((column.removeprefix('stream:'), number(column)) for column in stream_columns)
    for value, label in ((delivered, 'Delivered thrust'), (axial, 'Axial thrust'),
                         (flow, 'Total flow')):
        _nonnegative(value, label)
    if _fraction(axial) > _fraction(delivered)*_fraction(1.0+1e-12):
        raise InputError('CSV axial thrust cannot exceed delivered thrust.')
    _flow_pairs(streams, 'CSV streams')
    TabulatedBurnSegment._validate_stream_closure(streams, flow, 'CSV')
    phase = _utf8(row['phase'], 'CSV phase') if 'phase' in row else 'imported-trace'
    return {'time': number('time_s'), 'delivered': delivered, 'axial': axial,
            'flow': flow, 'streams': streams, 'phase': phase}


def tabulated_artifact_from_csv(
    text: str,
    *,
    source_id: str,
    interpolation: InterpolationPolicy = InterpolationPolicy.LINEAR,
) -> TabulatedBurnArtifact:
    """Import a strictly declared SI trace with explicit interpolation semantics.

    Parameters
    ----------
    text : str
        Non-empty, UTF-8 encodable CSV text with a header row.
    source_id : str
        Non-blank identifier of the trace source.
    interpolation : InterpolationPolicy, default InterpolationPolicy.LINEAR
        Interpolation between consecutive rows.

    Returns
    -------
    TabulatedBurnArtifact
        One segment per pair of consecutive rows; the source hash is the
        SHA-256 digest of the literal UTF-8 text.

    Raises
    ------
    InputError
        If the policy is unknown, the text is empty or not UTF-8 encodable,
        the header lacks a required column or has duplicate or unsupported
        columns, a record is ragged or malformed, fewer than two rows are
        given, the first time is not zero, times are not strictly increasing,
        or a cell is invalid.
    DomainError
        If a row violates the segment contract.

    Notes
    -----
    Required columns are ``time_s``, ``delivered_thrust_n`` and
    ``total_mass_flow_kg_s``.  ``axial_thrust_n`` defaults to delivered
    thrust and ``phase`` to ``"imported-trace"``; defaults apply only to
    absent columns, and every declared cell must be present and valid.  Each
    column named ``stream:<tank-id>`` is preserved as an explicit tank drain.
    ``HOLD`` applies row ``i`` to ``[t_i, t_(i+1))``; the terminal row is
    validated and hashed but defines no further interval.  Samples are not
    smoothed.  Measured start, steady and shutdown histories of the kind
    compared in [1]_ and the ignition transients of [2]_ are typical sources.

    References
    ----------
    .. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    .. [2] NASA CR-140800. NASA, "Investigation of the Starting Transients of
           High Performance Solid-Propellant Motors".
           https://ntrs.nasa.gov/citations/19750003988
    """

    try:
        interpolation = InterpolationPolicy(_utf8(interpolation, 'CSV interpolation'))
    except (ValueError, TypeError) as error:
        raise InputError("Unknown CSV interpolation policy.", field="interpolation") from error
    if not isinstance(text, str) or not text.strip():
        raise InputError("Imported burn trace must not be empty.")
    try:
        source_hash = hashlib.sha256(text.encode("utf-8")).hexdigest()
    except UnicodeError as error:
        raise InputError('Imported CSV must be UTF-8 encodable.') from error
    reader = csv.DictReader(io.StringIO(text), strict=True)
    required = {"time_s", "delivered_thrust_n", "total_mass_flow_kg_s"}
    points: list[dict[str, Any]] = []
    try:
        names = reader.fieldnames
        if names is None or not required.issubset(names):
            raise InputError('Imported trace requires time_s, delivered_thrust_n, and total_mass_flow_kg_s.')
        if len(names) != len(set(names)):
            raise InputError('CSV column names must be unique.')
        stream_columns = sorted(name for name in names if name.startswith('stream:'))
        if set(names)-required-{'axial_thrust_n', 'phase'}-set(stream_columns):
            raise InputError('CSV contains unsupported columns; declare supported SI fields explicitly.')
        for column in stream_columns:
            _utf8(column.removeprefix('stream:'), 'CSV stream identifier')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise InputError('CSV record width does not match its header.')
            points.append(_csv_point(row, stream_columns))
    except csv.Error as error:
        raise InputError('Malformed CSV quoting or record.', field='CSV') from error
    if len(points) < 2:
        raise InputError("Imported trace requires at least two samples.")
    if points[0]["time"] != 0.0:
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
        warnings=("Imported SI samples were validated without smoothing; interval values follow declared interpolation.",),
    )
