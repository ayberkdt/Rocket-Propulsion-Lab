"""One-dimensional steady conduction through a regenerative chamber wall.

The hot-gas-side wall temperature follows from the planar conduction
relation of NASA SP-125, equation (4-19) [1]_,

``q = (k / t) (T_wg - T_wc)``  so  ``T_wg = T_wc + q t / k``,

with the coolant-side wall temperature ``T_wc`` either prescribed or
obtained from the coolant-side film relation ``q = h_c (T_wc - T_co)``,
equation (4-20) of the same source. The result is compared with the wall
material's sourced temperature limit from
:mod:`rocket_propulsion.thermodynamics.materials`.

The model is a constant-property, planar, single-layer estimate. It does not
include wall curvature, rib/fin conduction, two-dimensional spreading,
coatings, contact resistance, thermal stress, creep or fatigue, and results
carry ``physical_model_qualified=False``.

References
----------
.. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
       Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
       printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [2] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
       718", Publication SMC-045, 2007,
       https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real, records, strings, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational

from .heat_transfer import _owned
from .materials import get_wall_material

WALL_MODEL = "planar_steady_constant_conductivity_wall_v1"
_ASSUMPTIONS = (
    "steady one-dimensional planar conduction through a single wall layer",
    "constant wall conductivity at the catalog reference condition or a user value",
    "prescribed local heat flux; no curvature, rib/fin, coating or contact resistance",
    "limit check against a sourced temperature only; no stress, creep or fatigue",
)
_OPTIONAL = ("coolant_side_wall_temperature_k", "coolant_side_heat_transfer_coefficient_w_m2_k",
             "coolant_bulk_temperature_k", "wall_thermal_conductivity_w_m_k")


@dataclass(frozen=True, slots=True)
class WallConductionInputs:
    """Owned, validated inputs of the wall-conduction step.

    Attributes
    ----------
    wall_material : str
        Key of :func:`~rocket_propulsion.thermodynamics.materials.get_wall_material`;
        normalized to the catalog key on construction.
    wall_thickness_m : float
        Wall thickness between the hot-gas and coolant surfaces, m.
    heat_flux_w_m2 : float
        Local heat flux through the wall, W/m^2 (for example a gas-side
        Bartz station value).
    coolant_side_wall_temperature_k : float or None
        Prescribed coolant-side wall temperature ``T_wc``, K.
    coolant_side_heat_transfer_coefficient_w_m2_k : float or None
        Coolant-side film coefficient ``h_c``; then
        ``T_wc = T_co + q / h_c`` [1]_.
    coolant_bulk_temperature_k : float or None
        Coolant bulk temperature ``T_co``, K; required with ``h_c`` only.
    wall_thermal_conductivity_w_m_k : float or None
        User conductivity replacing the catalog value, W/(m K).

    Raises
    ------
    InputError
        If the material key is unknown, not exactly one of ``T_wc`` and
        ``h_c`` is given (code ``"wall_coolant_side_ambiguous"``), the bulk
        temperature is missing for or given without ``h_c``, or no
        conductivity is available (code ``"wall_conductivity_unavailable"``).
    DomainError
        If a quantity is not a finite positive real.

    References
    ----------
    .. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
           Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
           printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """

    wall_material: str
    wall_thickness_m: float
    heat_flux_w_m2: float
    coolant_side_wall_temperature_k: float | None = None
    coolant_side_heat_transfer_coefficient_w_m2_k: float | None = None
    coolant_bulk_temperature_k: float | None = None
    wall_thermal_conductivity_w_m_k: float | None = None

    def __post_init__(self) -> None:
        record = get_wall_material(self.wall_material)
        object.__setattr__(self, "wall_material", record.key)
        for name in ("wall_thickness_m", "heat_flux_w_m2", *_OPTIONAL):
            value = getattr(self, name)
            if value is None and name in _OPTIONAL:
                continue
            value = real(value, name)
            if value <= 0:
                raise DomainError("Wall-conduction quantities must be positive.", field=name)
            object.__setattr__(self, name, value)
        prescribed = self.coolant_side_wall_temperature_k is not None
        film = self.coolant_side_heat_transfer_coefficient_w_m2_k is not None
        if prescribed == film:
            raise InputError("Give exactly one of the coolant-side wall temperature and the "
                             "coolant-side heat-transfer coefficient.",
                             code="wall_coolant_side_ambiguous",
                             field="coolant_side_wall_temperature_k")
        if film != (self.coolant_bulk_temperature_k is not None):
            raise InputError("The coolant bulk temperature is used with, and only with, the "
                             "coolant-side heat-transfer coefficient.",
                             field="coolant_bulk_temperature_k")
        if self.wall_thermal_conductivity_w_m_k is None and record.thermal_conductivity is None:
            raise InputError(f"{record.key} has no verified catalog conductivity; supply "
                             "wall_thermal_conductivity_w_m_k.",
                             code="wall_conductivity_unavailable",
                             field="wall_thermal_conductivity_w_m_k")


def _optional_float(value: Fraction, label: str) -> float | None:
    """Round an optional diagnostic; ``None`` when it is not representable."""
    try:
        return _finite_rational(value, label)
    except DomainError:
        return None


def _wall_values(inputs: WallConductionInputs) -> dict[str, object]:
    """Return every result field of the wall step for ``inputs`` [1]_.

    All temperatures are formed exactly from the binary64 inputs and rounded
    once; the limit decision compares exact values.

    References
    ----------
    .. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
           Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
           printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    record = get_wall_material(inputs.wall_material)
    catalog = record.thermal_conductivity
    user = inputs.wall_thermal_conductivity_w_m_k
    conductivity = user if user is not None else catalog.value
    flux = Fraction(inputs.heat_flux_w_m2)
    if inputs.coolant_side_wall_temperature_k is not None:
        coolant_wall = Fraction(inputs.coolant_side_wall_temperature_k)
    else:
        coolant_wall = (Fraction(inputs.coolant_bulk_temperature_k)
                        + flux / Fraction(inputs.coolant_side_heat_transfer_coefficient_w_m2_k))
    drop = flux * Fraction(inputs.wall_thickness_m) / Fraction(conductivity)
    hot = coolant_wall + drop
    limit, basis = record.temperature_limit()
    margin = Fraction(limit.value) - hot
    exported_margin = _optional_float(margin, "Wall limit margin")
    hot_out = _finite_rational(hot, "Hot-wall temperature")
    wall_out = _finite_rational(coolant_wall, "Coolant-side wall temperature")
    warnings = [WarningMessage(
        "wall_model_unqualified",
        "One-dimensional steady planar conduction with constant conductivity; curvature, "
        "ribs, coatings, contact resistance and thermal stress are not modelled.")]
    if user is None:
        where = ("at an unstated temperature" if catalog.temperature_k is None
                 else f"at {catalog.temperature_k:.6g} K")
        warnings.append(WarningMessage(
            "wall_conductivity_reference_condition",
            f"Catalog conductivity of {record.key} ({conductivity:.6g} W/(m K)) is stated "
            f"{where}; it is held constant between {float(wall_out):.6g} K and "
            f"{float(hot_out):.6g} K.", "warning", "wall_material"))
        if catalog.upper_value is not None:
            warnings.append(WarningMessage(
                "wall_conductivity_range_lower_bound",
                f"The source gives a range up to {catalog.upper_value:.6g} W/(m K); its lower "
                "end is used, which gives the larger wall temperature drop.", "info",
                "wall_material"))
    else:
        warnings.append(WarningMessage(
            "wall_conductivity_user_value",
            "User conductivity replaces the catalog value; its source is the caller's.", "info",
            "wall_thermal_conductivity_w_m_k"))
    if basis == "solidus_temperature":
        warnings.append(WarningMessage(
            "wall_limit_is_melting_point",
            f"No sourced service temperature for {record.key}; the limit is its melting "
            "point, which is not a strength or life limit.", "warning", "wall_material"))
    if margin < 0:
        warnings.append(WarningMessage(
            "wall_temperature_limit_exceeded",
            f"Hot-gas-side wall temperature exceeds the {basis.replace('_', ' ')} of "
            f"{record.key}.", "critical", "wall_material"))
    return {
        "wall_thermal_conductivity_w_m_k": conductivity,
        "conductivity_source": "user" if user is not None else "catalog",
        "coolant_side_wall_temperature_k": wall_out,
        "wall_temperature_drop_k": _finite_rational(drop, "Wall temperature drop"),
        "hot_wall_temperature_k": hot_out,
        "wall_limit_temperature_k": limit.value,
        "wall_limit_basis": basis,
        "wall_limit_exceeded": margin < 0,
        "wall_limit_margin_k": exported_margin,
        "wall_limit_margin_unrepresentable": exported_margin is None,
        "warnings": tuple(warnings),
        "assumptions": _ASSUMPTIONS,
        "physical_model_qualified": False,
        "model": WALL_MODEL,
    }


def _owned_warnings(value: object, label: str) -> tuple[WarningMessage, ...]:
    """Return warning records (or mappings of their fields) as typed records."""
    owned = []
    for item in records(value, (WarningMessage, Mapping), label):
        owned.append(_owned(item, WarningMessage, label))
    return tuple(owned)


@dataclass(frozen=True, slots=True)
class WallConductionResult:
    """Replayed wall-conduction step with an explicit material-limit flag.

    Attributes
    ----------
    inputs : WallConductionInputs
        Owned inputs from which the record is replayed.
    wall_thermal_conductivity_w_m_k : float
        Conductivity used, W/(m K).
    conductivity_source : {"catalog", "user"}
        Origin of the conductivity.
    coolant_side_wall_temperature_k : float
        ``T_wc``, K.
    wall_temperature_drop_k : float
        ``q t / k``, K.
    hot_wall_temperature_k : float
        ``T_wg = T_wc + q t / k``, K.
    wall_limit_temperature_k : float
        Sourced material limit, K.
    wall_limit_basis : {"maximum_service_temperature", "solidus_temperature"}
        Catalog field that supplied the limit.
    wall_limit_exceeded : bool
        Exact comparison of ``T_wg`` with the limit.
    wall_limit_margin_k : float or None
        Limit minus exact ``T_wg``; ``None`` when not representable.
    wall_limit_margin_unrepresentable : bool
        True when the margin is ``None`` for that reason.
    warnings : tuple of WarningMessage
        Model, property-condition and limit warnings.
    assumptions : tuple of str
        Declared modelling assumptions.
    physical_model_qualified : bool
        Always False.
    model : str
        Model identifier.

    Raises
    ------
    DomainError
        If a field is malformed or differs from the replay of ``inputs``.

    References
    ----------
    .. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
           Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
           printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """

    inputs: WallConductionInputs
    wall_thermal_conductivity_w_m_k: float
    conductivity_source: str
    coolant_side_wall_temperature_k: float
    wall_temperature_drop_k: float
    hot_wall_temperature_k: float
    wall_limit_temperature_k: float
    wall_limit_basis: str
    wall_limit_exceeded: bool
    wall_limit_margin_k: float | None
    wall_limit_margin_unrepresentable: bool
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...]
    physical_model_qualified: bool = False
    model: str = WALL_MODEL

    def __post_init__(self) -> None:
        source = _owned(self.inputs, WallConductionInputs, "Wall inputs")
        object.__setattr__(self, "inputs", source)
        for name in ("wall_limit_exceeded", "wall_limit_margin_unrepresentable",
                     "physical_model_qualified"):
            if type(getattr(self, name)) is not bool:
                raise DomainError("Wall flags must be actual booleans.", field=name)
        for name in ("conductivity_source", "wall_limit_basis", "model"):
            object.__setattr__(self, name, text(getattr(self, name), name))
        object.__setattr__(self, "warnings", _owned_warnings(self.warnings, "Wall warnings"))
        object.__setattr__(self, "assumptions", strings(self.assumptions, "Wall assumptions"))
        if any(getattr(self, name) != value for name, value in _wall_values(source).items()):
            raise DomainError("Wall result does not replay from owned inputs.",
                              field="Wall conduction result")


def calculate_wall_conduction(
    *, wall_material: str, wall_thickness_m: float, heat_flux_w_m2: float,
    coolant_side_wall_temperature_k: float | None = None,
    coolant_side_heat_transfer_coefficient_w_m2_k: float | None = None,
    coolant_bulk_temperature_k: float | None = None,
    wall_thermal_conductivity_w_m_k: float | None = None,
) -> WallConductionResult:
    """Estimate the hot-gas-side wall temperature and check it against the material limit.

    Parameters
    ----------
    wall_material : str
        Wall-material catalog key, for example ``"grcop-84"``.
    wall_thickness_m : float
        Wall thickness, m.
    heat_flux_w_m2 : float
        Local heat flux through the wall, W/m^2.
    coolant_side_wall_temperature_k : float, optional
        Prescribed ``T_wc``, K. Give this or ``h_c``.
    coolant_side_heat_transfer_coefficient_w_m2_k : float, optional
        Coolant film coefficient ``h_c``, W/(m^2 K); requires
        ``coolant_bulk_temperature_k``.
    coolant_bulk_temperature_k : float, optional
        Coolant bulk temperature ``T_co``, K.
    wall_thermal_conductivity_w_m_k : float, optional
        Conductivity replacing the catalog value, W/(m K).

    Returns
    -------
    WallConductionResult
        Wall temperatures, the material limit, a limit flag and warnings.

    Raises
    ------
    InputError
        If the material is unknown, the coolant-side specification is
        ambiguous or incomplete, or no conductivity is available.
    DomainError
        If a quantity is not finite and positive or an output overflows.

    Notes
    -----
    ``T_wg = T_wc + q t / k`` from ``q = (k/t)(T_wg - T_wc)``, and, with a
    film coefficient, ``T_wc = T_co + q / h_c`` from ``q = h_c (T_wc - T_co)``
    [1]_. Conductivity is taken at its catalog reference condition; for
    alloys whose conductivity rises with temperature (as tabulated for
    Inconel 718 [2]_), a room-temperature value overstates the drop. The
    limit is the sourced maximum service temperature when the catalog has
    one and otherwise the melting point; the latter is not a strength limit.

    References
    ----------
    .. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
           Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
           printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
           718", Publication SMC-045, 2007, Table 5,
           https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
    """
    source = WallConductionInputs(
        wall_material, wall_thickness_m, heat_flux_w_m2, coolant_side_wall_temperature_k,
        coolant_side_heat_transfer_coefficient_w_m2_k, coolant_bulk_temperature_k,
        wall_thermal_conductivity_w_m_k)
    return WallConductionResult(inputs=source, **_wall_values(source))
