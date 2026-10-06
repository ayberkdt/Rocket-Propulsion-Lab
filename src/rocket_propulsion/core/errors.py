"""Typed exceptions shared by calculation and transport layers."""

from __future__ import annotations

from typing import Any


class RocketPropulsionError(ValueError):
    """Base error with a stable machine code and optional structured context."""

    default_code = "calculation_error"

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        field: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.code = code or self.default_code
        self.field = field
        self.details = details or {}

    def as_dict(self) -> dict[str, Any]:
        """Return the public, JSON-safe error contract."""

        result: dict[str, Any] = {"code": self.code, "message": self.message}
        if self.field is not None:
            result["field"] = self.field
        if self.details:
            result["details"] = self.details
        return result


class DomainError(RocketPropulsionError):
    """Raised when a value is outside the physical domain of an equation."""

    default_code = "physical_domain"


class InputError(DomainError):
    """Raised when a request field is absent, malformed, or unsupported."""

    default_code = "invalid_input"


class ConvergenceError(DomainError):
    """Raised when a valid numerical problem does not converge."""

    default_code = "non_convergence"


class FeatureUnavailableError(RocketPropulsionError):
    """Raised when an optional external capability is not installed."""

    default_code = "feature_unavailable"

