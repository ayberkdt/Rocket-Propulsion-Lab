"""Shared model provenance, validity, and warning value objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True, slots=True)
class ValidityRange:
    """Document a model's supported interval for one input quantity."""

    quantity: str
    minimum: float | None = None
    maximum: float | None = None
    unit: str = "1"
    inclusive_minimum: bool = True
    inclusive_maximum: bool = True


@dataclass(frozen=True, slots=True)
class WarningMessage:
    """A stable warning suitable for Python, JSON, and user interfaces."""

    code: str
    message: str
    severity: Literal["info", "warning", "critical"] = "warning"
    field: str | None = None


@dataclass(frozen=True, slots=True)
class ModelMetadata:
    """Identify the equation model, assumptions, units, and fidelity limits."""

    model: str
    version: str
    assumptions: tuple[str, ...]
    units: dict[str, str]
    validity: tuple[ValidityRange, ...] = ()
    references: tuple[str, ...] = ()

