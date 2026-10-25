"""Shared model provenance, validity and warning value objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeGuard

WarningSeverity = Literal["info", "warning", "critical"]
"""Admissible :attr:`WarningMessage.severity` values."""

_WARNING_SEVERITIES: tuple[WarningSeverity, ...] = ("info", "warning", "critical")


def is_warning_severity(value: object) -> TypeGuard[WarningSeverity]:
    """Return whether ``value`` is an admissible warning severity.

    Parameters
    ----------
    value : object
        Candidate severity.

    Returns
    -------
    bool
        True if ``value`` is one of ``"info"``, ``"warning"`` or ``"critical"``.
    """
    return value in _WARNING_SEVERITIES


@dataclass(frozen=True, slots=True)
class ValidityRange:
    """Supported interval of one input quantity of a model.

    Attributes
    ----------
    quantity : str
        Name of the bounded input quantity.
    minimum, maximum : float or None
        Interval endpoints; ``None`` means unbounded on that side.
    unit : str
        Unit of the endpoints (``"1"`` for dimensionless quantities).
    inclusive_minimum, inclusive_maximum : bool
        Whether each finite endpoint belongs to the supported interval.

    Notes
    -----
    A validity range documents the declared model domain. It is not, by
    itself, evidence of experimental validation inside that domain.
    """

    quantity: str
    minimum: float | None = None
    maximum: float | None = None
    unit: str = "1"
    inclusive_minimum: bool = True
    inclusive_maximum: bool = True


@dataclass(frozen=True, slots=True)
class WarningMessage:
    """Stable warning record shared by the Python, JSON and UI layers.

    Attributes
    ----------
    code : str
        Stable machine-readable identifier.
    message : str
        Human-readable explanation.
    severity : {"info", "warning", "critical"}
        Importance of the warning.
    field : str or None
        Input field the warning refers to, if any.
    """

    code: str
    message: str
    severity: WarningSeverity = "warning"
    field: str | None = None


@dataclass(frozen=True, slots=True)
class ModelMetadata:
    """Identify an equation model, its assumptions, units and fidelity limits.

    Attributes
    ----------
    model : str
        Short model name.
    version : str
        Model version string.
    assumptions : tuple of str
        Modelling assumptions and declared limitations.
    units : dict of str to str
        Output-field name patterns mapped to their units.
    validity : tuple of ValidityRange
        Declared input domains, if any.
    references : tuple of str
        Short source labels; the full records are indexed in
        ``docs/REFERENCES.md``.
    """

    model: str
    version: str
    assumptions: tuple[str, ...]
    units: dict[str, str]
    validity: tuple[ValidityRange, ...] = ()
    references: tuple[str, ...] = ()
