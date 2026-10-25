"""Explicit SI conversion registry used at application boundaries.

The domain layer is SI-only. These helpers serve request, display, import
and export adapters, so that unit conversion never enters the physics
equations.

References
----------
.. [1] BIPM, "The International System of Units (SI Brochure)", 9th ed., 2019,
       https://www.bipm.org/en/publications/si-brochure
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from .constants import CELSIUS_OFFSET_K
from .errors import InputError


@dataclass(frozen=True, slots=True)
class UnitDefinition:
    """Affine conversion between one display unit and its SI unit.

    Attributes
    ----------
    symbol : str
        Display symbol, e.g. ``"kPa"`` or ``"degC"``.
    quantity : str
        Physical quantity, e.g. ``"pressure"``; conversions are only valid
        between units of the same quantity.
    scale_to_si : float
        Multiplicative factor from the display unit to SI.
    offset_to_si : float
        Additive offset applied after scaling (non-zero only for
        temperature scales with a shifted zero).
    display_precision : int
        Suggested number of significant digits for display.

    Notes
    -----
    ``value_si = value * scale_to_si + offset_to_si``. The Celsius offset is
    the registered constant ``CELSIUS_OFFSET_K`` defined by the SI
    Brochure [1]_.

    References
    ----------
    .. [1] BIPM, "The International System of Units (SI Brochure)", 9th ed.,
           2019, https://www.bipm.org/en/publications/si-brochure
    """

    symbol: str
    quantity: str
    scale_to_si: float
    offset_to_si: float = 0.0
    display_precision: int = 6

    def to_si(self, value: float) -> float:
        """Convert a finite value in this unit to the quantity's SI unit.

        Parameters
        ----------
        value : float
            Finite value expressed in this unit.

        Returns
        -------
        float
            ``value * scale_to_si + offset_to_si``.

        Raises
        ------
        InputError
            If ``value`` is not finite.
        """

        _finite(value)
        return value * self.scale_to_si + self.offset_to_si

    def from_si(self, value: float) -> float:
        """Convert a finite SI value to this display unit.

        Parameters
        ----------
        value : float
            Finite value expressed in the SI unit of the quantity.

        Returns
        -------
        float
            ``(value - offset_to_si) / scale_to_si``.

        Raises
        ------
        InputError
            If ``value`` is not finite.
        """

        _finite(value)
        return (value - self.offset_to_si) / self.scale_to_si


def _finite(value: float) -> None:
    """Raise :class:`InputError` unless ``value`` is a finite number."""
    if not isfinite(value):
        raise InputError("Unit conversion value must be finite.", field="value")


_UNITS = (
    UnitDefinition("Pa", "pressure", 1.0, display_precision=3),
    UnitDefinition("kPa", "pressure", 1_000.0, display_precision=6),
    UnitDefinition("MPa", "pressure", 1_000_000.0, display_precision=6),
    UnitDefinition("bar", "pressure", 100_000.0, display_precision=6),
    UnitDefinition("psi", "pressure", 6_894.757293168, display_precision=6),
    UnitDefinition("K", "temperature", 1.0, display_precision=4),
    UnitDefinition("degC", "temperature", 1.0, CELSIUS_OFFSET_K, display_precision=4),
    UnitDefinition("degF", "temperature", 5.0 / 9.0, 255.3722222222222, 4),
    UnitDefinition("m", "length", 1.0, display_precision=6),
    UnitDefinition("mm", "length", 0.001, display_precision=4),
    UnitDefinition("in", "length", 0.0254, display_precision=6),
    UnitDefinition("m2", "area", 1.0, display_precision=8),
    UnitDefinition("cm2", "area", 0.0001, display_precision=6),
    UnitDefinition("mm2", "area", 0.000001, display_precision=4),
    UnitDefinition("m/s", "velocity", 1.0, display_precision=5),
    UnitDefinition("ft/s", "velocity", 0.3048, display_precision=5),
    UnitDefinition("kg/s", "mass_flow", 1.0, display_precision=6),
    UnitDefinition("lbm/s", "mass_flow", 0.45359237, display_precision=6),
    UnitDefinition("kg", "mass", 1.0, display_precision=6),
    UnitDefinition("lbm", "mass", 0.45359237, display_precision=6),
    UnitDefinition("N", "force", 1.0, display_precision=4),
    UnitDefinition("kN", "force", 1_000.0, display_precision=6),
    UnitDefinition("lbf", "force", 4.4482216152605, display_precision=6),
    UnitDefinition("J/kg/K", "specific_heat", 1.0, display_precision=5),
    UnitDefinition("kJ/kg/K", "specific_heat", 1_000.0, display_precision=6),
    UnitDefinition("s", "time", 1.0, display_precision=5),
    UnitDefinition("deg", "angle", 1.0, display_precision=5),
    UnitDefinition("1", "dimensionless", 1.0, display_precision=8),
)

UNIT_REGISTRY = {definition.symbol: definition for definition in _UNITS}


def get_unit(symbol: str, *, quantity: str | None = None) -> UnitDefinition:
    """Return a registered unit, optionally checking its physical quantity.

    Parameters
    ----------
    symbol : str
        Registered unit symbol (a key of ``UNIT_REGISTRY``).
    quantity : str, optional
        Expected physical quantity; when given, a unit of another quantity
        is refused.

    Returns
    -------
    UnitDefinition
        The registered definition.

    Raises
    ------
    InputError
        With code ``"unsupported_unit"`` if the symbol is not registered, or
        ``"unit_quantity_mismatch"`` if it belongs to another quantity.
    """

    try:
        definition = UNIT_REGISTRY[symbol]
    except KeyError as error:
        raise InputError(
            f"Unsupported unit: {symbol!r}.",
            code="unsupported_unit",
            field="unit",
            details={"supported": sorted(UNIT_REGISTRY)},
        ) from error
    if quantity is not None and definition.quantity != quantity:
        raise InputError(
            f"Unit {symbol!r} is not valid for {quantity}.",
            code="unit_quantity_mismatch",
            field="unit",
            details={"expected_quantity": quantity, "actual_quantity": definition.quantity},
        )
    return definition


def to_si(value: float, symbol: str, *, quantity: str | None = None) -> float:
    """Convert a value from the unit ``symbol`` to the registry's SI unit.

    Parameters
    ----------
    value : float
        Finite value expressed in ``symbol``.
    symbol : str
        Registered unit symbol.
    quantity : str, optional
        Expected physical quantity of ``symbol``.

    Returns
    -------
    float
        The value in SI units.

    Raises
    ------
    InputError
        If the unit is unknown, belongs to another quantity, or ``value`` is
        not finite.
    """

    return get_unit(symbol, quantity=quantity).to_si(value)


def from_si(value: float, symbol: str, *, quantity: str | None = None) -> float:
    """Convert a value from the registry's SI unit to the unit ``symbol``.

    Parameters
    ----------
    value : float
        Finite value in SI units.
    symbol : str
        Registered target unit symbol.
    quantity : str, optional
        Expected physical quantity of ``symbol``.

    Returns
    -------
    float
        The value expressed in ``symbol``.

    Raises
    ------
    InputError
        If the unit is unknown, belongs to another quantity, or ``value`` is
        not finite.
    """

    return get_unit(symbol, quantity=quantity).from_si(value)
