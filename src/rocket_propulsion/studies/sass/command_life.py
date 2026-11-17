"""Per-head commanded valve use over an explicitly bounded observation window."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import groupby
from math import isfinite, nextafter
from operator import attrgetter
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, records, strings, text

from .catalog import ThrusterRecord
from .evidence import CitedScalar, quantity


@dataclass(frozen=True, slots=True)
class ValveCommand:
    """One state transition at an integer clock tick for a specific head.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6, technical data.
    """

    tick: int
    head_id: str
    opened: bool

    def __post_init__(self) -> None:
        """Validate timestamp, head identity and a binary command state.

        References
        ----------
        .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6, input traceability.
        """
        tick = integer(self.tick, "valve command tick")
        if tick < 0 or not isinstance(self.opened, bool):
            raise DomainError("Valve command requires a nonnegative tick and an explicit boolean state.")
        object.__setattr__(self, "tick", tick)
        object.__setattr__(self, "head_id", text(self.head_id, "valve command head ID"))


@dataclass(frozen=True, slots=True)
class ThrusterCommandHistory:
    """Ordered transitions and head states immediately before tick zero.

    Include all installed heads, including those with no commands. Events
    at the end tick are admitted for cycle closure and right censoring.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.5 and 6.6.
    """

    source_record: str
    tick_period: CitedScalar
    duration_ticks: int
    head_ids: tuple[str, ...]
    initially_open: tuple[str, ...]
    commands: tuple[ValveCommand, ...]

    def __post_init__(self) -> None:
        """Freeze the declared identities, timing and command order.

        References
        ----------
        .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6, controlled data.
        """
        tick_period = quantity(self.tick_period, "s", strictly_positive=True)
        duration = integer(self.duration_ticks, "command history duration ticks")
        heads = strings(self.head_ids, "installed command head IDs")
        initial = strings(self.initially_open, "initially open head IDs")
        events = records(self.commands, ValveCommand, "valve commands")
        if duration <= 0 or not heads or len(set(heads)) != len(heads):
            raise DomainError("Command history requires a positive duration and unique installed head IDs.")
        real(Fraction(tick_period) * duration, "command history duration seconds")
        if len(set(initial)) != len(initial) or not set(initial).issubset(heads):
            raise DomainError("Initially open heads must be unique installed IDs.")
        object.__setattr__(self, "source_record", text(self.source_record, "command history source record"))
        object.__setattr__(self, "duration_ticks", duration)
        object.__setattr__(self, "head_ids", heads)
        object.__setattr__(self, "initially_open", initial)
        object.__setattr__(self, "commands", events)


@dataclass(frozen=True, slots=True)
class HeadCommandUsage:
    """Commanded open time and completed command cycles for one head.

    Censoring flags identify open intervals intersecting a window boundary.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6, retained input scope.
    """

    head_id: str
    observed_open_ticks: int
    longest_open_ticks: int
    completed_cycles: int
    left_censored: bool
    right_censored: bool


@dataclass(frozen=True, slots=True)
class CommandLifeComparison:
    """Compare observed command exposure against a catalogue reference.

    References
    ----------
    .. [1] ARIANE-1N. ArianeGroup 1 N thruster, catalogue life table.
    .. [2] ARIANE-10N. ArianeGroup 10 N family, catalogue life table.
    """

    head_id: str
    metric: str
    demand: float | int
    unit: str
    published_limit: CitedScalar | None
    status: Literal["exceeds_reference", "no_exceedance_identified", "unassessed"]


@dataclass(frozen=True, slots=True)
class CommandLifeReview:
    """Ledger and comparisons for the supplied command window.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6 and 6.8.
    """

    source_record: str
    tick_period: CitedScalar
    duration_ticks: int
    peak_simultaneous_open: int
    usage: tuple[HeadCommandUsage, ...]
    comparisons: tuple[CommandLifeComparison, ...]


def review_thruster_command_life(
    product: ThrusterRecord, history: ThrusterCommandHistory,
) -> CommandLifeReview:
    """Count full command cycles and bound open times using exact tick arithmetic.

    Commands on different heads at the same tick take effect together.
    Contradictory transitions for a head or redundant repeated states are
    rejected. A cycle needs both its opening and closing inside the
    declared window. Censored durations contribute their observed portion.
    The result describes commanded states and leaves mechanical response,
    material exposure and catalogue applicability to separate evidence.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6 and 6.8.
    .. [2] ARIANE-1N. ArianeGroup 1 N thruster, catalogue life references.
    """
    if not isinstance(product, ThrusterRecord) or not isinstance(history, ThrusterCommandHistory):
        raise DomainError("Command life review requires a catalogue record and a command history.")
    usage, peak = _command_usage(history)
    return _compare_command_life(product, history, usage, peak)


def _command_usage(history: ThrusterCommandHistory) -> tuple[tuple[HeadCommandUsage, ...], int]:
    """Validate the declared sequence and aggregate its observed intervals.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6.
    """
    total = dict.fromkeys(history.head_ids, 0)
    longest = dict.fromkeys(history.head_ids, 0)
    cycles = dict.fromkeys(history.head_ids, 0)
    left_censored = set(history.initially_open)
    active = dict.fromkeys(history.initially_open, 0)
    opened_in_window: set[str] = set()
    previous_tick = peak = 0
    for tick, simultaneous in groupby(history.commands, key=attrgetter("tick")):
        if tick < previous_tick or tick > history.duration_ticks:
            raise DomainError("Valve commands must be ordered inside the declared history window.")
        if tick > previous_tick:
            peak = max(peak, len(active))
        seen = set()
        for event in simultaneous:
            head = event.head_id
            if head not in total or head in seen:
                raise DomainError("A tick may contain at most one command per installed head.")
            seen.add(head)
            if event.opened == (head in active):
                raise DomainError("Valve command must change the declared state of its head.")
            if event.opened:
                active[head] = tick
                opened_in_window.add(head)
            else:
                duration = tick - active.pop(head)
                total[head] += duration
                longest[head] = max(longest[head], duration)
                if head in opened_in_window:
                    cycles[head] += 1
                    opened_in_window.remove(head)
        previous_tick = tick
    if previous_tick < history.duration_ticks:
        peak = max(peak, len(active))
    for head, start in active.items():
        duration = history.duration_ticks - start
        total[head] += duration
        longest[head] = max(longest[head], duration)
    usage = tuple(HeadCommandUsage(head, total[head], longest[head], cycles[head],
                                   head in left_censored, head in active) for head in history.head_ids)
    return usage, peak


def _compare_command_life(
    product: ThrusterRecord, history: ThrusterCommandHistory,
    usage: tuple[HeadCommandUsage, ...], peak: int,
) -> CommandLifeReview:
    """Compare command exposure with the catalogue values for the chosen record.

    References
    ----------
    .. [1] ARIANE-1N. ArianeGroup 1 N thruster, catalogue life references.
    """
    period = Fraction(history.tick_period.value)
    # Bracket decimal-to-binary conversion of the declared clock period
    # conservatively. Counts remain exact and need no numerical guard.
    period_lower = Fraction(nextafter(history.tick_period.value, 0.0))
    next_period = nextafter(history.tick_period.value, float("inf"))
    period_upper = Fraction(next_period) if isfinite(next_period) else None
    fields = {field.name: field for field in product.fields}
    checks = []
    for metric, unit, attribute in (
        ("single_burn_life", "s", "longest_open_ticks"),
        ("accumulated_burn_life", "s", "observed_open_ticks"),
        ("cycle_life", "1", "completed_cycles"),
    ):
        field = fields.get(metric)
        limit = None if field is None else field.upper
        if field is not None:
            if field.kind != "upper_bound":
                raise DomainError(f"Catalogue {metric} must be a published upper reference.")
            quantity(field.upper, unit, strictly_positive=True)
        for item in usage:
            count = getattr(item, attribute)
            exact = Fraction(count) * period if unit == "s" else Fraction(count)
            demand = real(exact, metric + " command demand") if unit == "s" else count
            lower = period_lower * count if unit == "s" else exact
            upper = (exact if unit != "s" or count == 0 else
                     None if period_upper is None else period_upper * count)
            status: Literal["exceeds_reference", "no_exceedance_identified", "unassessed"] = "unassessed"
            if limit is not None:
                threshold = Fraction(limit.value)
                if lower > threshold:
                    status = "exceeds_reference"
                elif upper is not None and upper <= threshold and not (
                    item.left_censored or item.right_censored
                ):
                    status = "no_exceedance_identified"
            checks.append(CommandLifeComparison(item.head_id, metric, demand, unit, limit, status))
    return CommandLifeReview(history.source_record, history.tick_period, history.duration_ticks,
                             peak, usage, tuple(checks))
