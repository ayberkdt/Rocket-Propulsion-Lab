"""Shared validation, metadata, unit, and numerical utilities."""

from .errors import (
    ConvergenceError,
    DomainError,
    FeatureUnavailableError,
    InputError,
    RocketPropulsionError,
)
from .metadata import ModelMetadata, ValidityRange, WarningMessage
from .numerics import RootResult, bisect_monotonic, solve_bisect_monotonic
from .units import UNIT_REGISTRY, UnitDefinition, from_si, get_unit, to_si

__all__ = [
    "UNIT_REGISTRY",
    "ConvergenceError",
    "DomainError",
    "FeatureUnavailableError",
    "InputError",
    "ModelMetadata",
    "RocketPropulsionError",
    "RootResult",
    "UnitDefinition",
    "ValidityRange",
    "WarningMessage",
    "bisect_monotonic",
    "from_si",
    "get_unit",
    "solve_bisect_monotonic",
    "to_si",
]

