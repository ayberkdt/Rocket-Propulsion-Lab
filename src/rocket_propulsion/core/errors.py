"""Typed exceptions shared by the calculation and transport layers.

Every library error derives from :class:`RocketPropulsionError`, which carries
a stable machine-readable ``code``, an optional request ``field`` and optional
structured ``details``. The HTTP layer maps the subclasses to status codes and
serializes them with :meth:`RocketPropulsionError.as_dict`.
"""

from __future__ import annotations

from typing import Any


class RocketPropulsionError(ValueError):
    """Base error with a stable machine code and optional structured context.

    Parameters
    ----------
    message : str
        Human-readable description of the failure.
    code : str, optional
        Stable machine-readable error code. Defaults to the class attribute
        ``default_code``.
    field : str, optional
        Name of the input field the error refers to, when there is one.
    details : dict, optional
        Additional JSON-serializable context. ``None`` is stored as an empty
        dictionary.

    Attributes
    ----------
    message : str
        The message passed to the constructor.
    code : str
        The resolved machine-readable code.
    field : str or None
        The offending field, if any.
    details : dict
        Structured context; empty when none was supplied.

    Notes
    -----
    The class derives from :class:`ValueError`, so callers that only handle
    built-in value errors still catch library errors.
    """

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
        """Return the public error contract as a JSON-ready dictionary.

        Returns
        -------
        dict
            Always contains ``code`` and ``message``; ``field`` is present when
            set and ``details`` when non-empty.

        Notes
        -----
        The dictionary is JSON-safe provided that ``details`` holds only
        JSON-serializable values.
        """

        result: dict[str, Any] = {"code": self.code, "message": self.message}
        if self.field is not None:
            result["field"] = self.field
        if self.details:
            result["details"] = self.details
        return result


class DomainError(RocketPropulsionError):
    """Raised when a value is outside the physical domain of an equation.

    The default code is ``"physical_domain"``.
    """

    default_code = "physical_domain"


class InputError(DomainError):
    """Raised when a request field is absent, malformed or unsupported.

    The default code is ``"invalid_input"``. Being a :class:`DomainError`,
    it is caught by handlers written for domain violations.
    """

    default_code = "invalid_input"


class ConvergenceError(DomainError):
    """Raised when a well-posed numerical problem does not converge.

    The default code is ``"non_convergence"``. It also covers numerical
    enclosures that cannot be resolved within their precision budget.
    """

    default_code = "non_convergence"


class FeatureUnavailableError(RocketPropulsionError):
    """Raised when an optional external capability is not installed.

    The default code is ``"feature_unavailable"``.
    """

    default_code = "feature_unavailable"
