"""Replay selected engineering-unit CSV rows against an interpreted trace."""

from __future__ import annotations

from collections.abc import Iterator
from csv import Error as CsvError
from csv import reader
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from hashlib import sha256
from math import isfinite
from pathlib import Path
from re import fullmatch
from typing import BinaryIO, Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, text

from .transient_observation import TransientTrace

_DECIMAL = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"


@dataclass(frozen=True, slots=True, order=True)
class NativeCsvTrace:
    """Reviewed row range, channel and integer clock offset in a native CSV.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, documented transformation of referent data.
    """

    condition_record: str
    channel_id: str
    native_sha256: str
    value_column: str
    unit: Literal["N", "Pa"]
    first_index: int
    past_end_index: int
    tick_offset: int
    conversion_review_record: str
    absolute_error_column: str | None = None
    tick_period_s: float | None = None

    def __post_init__(self) -> None:
        """Require an explicit native column, range, unit and time offset.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, input representation.
        """
        for field in ("condition_record", "channel_id", "value_column", "conversion_review_record"):
            text(getattr(self, field), field)
        if not isinstance(self.native_sha256, str) or not fullmatch(r"[a-f0-9]{64}", self.native_sha256):
            raise DomainError("Native CSV requires a lowercase SHA-256 digest.")
        if self.unit not in ("Pa", "N") or self.value_column == "tick":
            raise DomainError("Native CSV requires a force or pressure column distinct from tick.")
        if self.absolute_error_column is not None:
            text(self.absolute_error_column, "absolute_error_column")
            if self.absolute_error_column in ("tick", self.value_column):
                raise DomainError("Native CSV absolute error requires its own column.")
        if self.tick_period_s is not None:
            period = real(self.tick_period_s, "reviewed native CSV tick period in seconds")
            if period <= 0:
                raise DomainError("Native CSV tick period must be positive.")
            object.__setattr__(self, "tick_period_s", period)
        first = integer(self.first_index, "first native CSV index")
        end = integer(self.past_end_index, "past-end native CSV index")
        offset = integer(self.tick_offset, "native CSV tick offset")
        if not 0 <= first < end:
            raise DomainError("Native CSV requires a nonempty half-open row range.")
        object.__setattr__(self, "first_index", first)
        object.__setattr__(self, "past_end_index", end)
        object.__setattr__(self, "tick_offset", offset)


def _hashed_lines(stream: BinaryIO, digest) -> Iterator[str]:
    """Feed CSV parsing and the digest from the same stream.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, traceable input bytes.
    """
    for line in stream:
        digest.update(line)
        yield line.decode("ascii")


def _exact_decimal(token: str) -> Fraction:
    """Convert a bounded decimal token to a rational without decimal-context rounding.

    References
    ----------
    .. [1] NASA-7009B, 4.2.7, treatment of numerical uncertainty.
    """
    if len(token) > 1024 or not fullmatch(_DECIMAL, token):
        raise ValueError("Invalid native decimal token.")
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError("Invalid native decimal exponent.") from exc
    if value.is_zero():
        return Fraction(0)
    if not -1075 <= value.adjusted() <= 308:
        raise ValueError("Native decimal outside the supported exact exponent range.")
    return Fraction(value)


def native_csv_trace_status(path: Path, scope: NativeCsvTrace, trace: TransientTrace) -> str:
    """Compare the selected ticks and binary64 sample centers without resampling.

    For a declared absolute error column, the binary64 interval must contain
    the exact decimal interval in every selected row. A declared tick period
    binds the numerical time scale to the approved native clock conversion.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, referents and measurement uncertainty.
    """
    if trace.unit != scope.unit:
        return "native_csv_unit_mismatch"
    if scope.tick_period_s is not None and trace.tick_period.value != scope.tick_period_s:
        return "native_csv_tick_period_mismatch"
    if len(trace.samples) != scope.past_end_index - scope.first_index:
        return "native_csv_sample_count_mismatch"
    digest = sha256()
    mismatch = None
    count = 0
    try:
        if not path.is_file():
            return "native_csv_unreadable"
        with path.open("rb") as stream:
            rows = reader(_hashed_lines(stream, digest), strict=True)
            header = next(rows, ())
            if not header or len(header) != len(set(header)) or "tick" not in header or scope.value_column not in header:
                return "native_csv_invalid_header"
            tick_column, value_column = header.index("tick"), header.index(scope.value_column)
            error_column = None
            if scope.absolute_error_column is not None:
                if scope.absolute_error_column not in header:
                    return "native_csv_invalid_header"
                error_column = header.index(scope.absolute_error_column)
            for index, row in enumerate(rows):
                if len(row) != len(header):
                    return "native_csv_invalid_row"
                if not scope.first_index <= index < scope.past_end_index:
                    continue
                tick_text, value_text = row[tick_column], row[value_column]
                if not fullmatch(r"[+-]?[0-9]+", tick_text) or not fullmatch(_DECIMAL, value_text):
                    return "native_csv_invalid_sample"
                tick, value = int(tick_text), float(value_text)
                if not isfinite(value):
                    return "native_csv_invalid_sample"
                sample = trace.samples[index - scope.first_index]
                if tick + scope.tick_offset != sample.tick or value.hex() != sample.value.hex():
                    mismatch = "native_csv_sample_mismatch"
                elif error_column is not None:
                    try:
                        center = _exact_decimal(value_text)
                    except ValueError:
                        return "native_csv_invalid_sample"
                    try:
                        radius = _exact_decimal(row[error_column])
                    except ValueError:
                        return "native_csv_invalid_absolute_error"
                    if radius < 0:
                        return "native_csv_invalid_absolute_error"
                    required = radius + abs(Fraction(sample.value) - center)
                    if Fraction(sample.absolute_error) < required:
                        mismatch = mismatch or "native_csv_uncertainty_narrowed"
                count += 1
    except (OSError, UnicodeError, CsvError, ValueError):
        return "native_csv_unreadable"
    if digest.hexdigest() != scope.native_sha256:
        return "native_csv_digest_mismatch"
    if count != len(trace.samples):
        return "native_csv_range_not_covered"
    return mismatch or "verified"
