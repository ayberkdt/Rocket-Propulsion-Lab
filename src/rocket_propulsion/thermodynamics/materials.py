"""Field-cited coolant and chamber-wall material catalogs for cooling estimates.

Every catalog number is a single published value read from a named source,
stored with that source's locator (table, page, row or query) and with the
reference condition at which the source states it. A value converted from
U.S. customary units also cites the conversion factor that was applied, and
keeps the value exactly as printed in ``source_value``. Properties that
could not be verified against a primary source are absent (``None``); the
catalog never fills a gap with an estimate.

The records are reference-condition snapshots, not property models. They do
not interpolate in temperature or pressure, and a value stated at one
condition is not evidence of the property at another.

References
----------
.. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
       "Thermophysical Properties of Fluid Systems", NIST Chemistry WebBook,
       SRD 69, https://webbook.nist.gov/chemistry/fluid/
.. [2] NISTIR-6646-RP1. J. W. Magee et al., "Thermophysical Properties
       Measurements and Models for Rocket Propellant RP-1: Phase I",
       NISTIR 6646, 2007,
       https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf
.. [3] AURUBIS-C10200. Aurubis, "Technical Datasheet C10200 (Cu-OF)",
       revision 18 08 US,
       https://www.aurubis.be/en/dam/jcr:90029b00-b968-486e-bde0-ce56f16022cd/c10200-cu-of-us.pdf
.. [4] NASA-TM-2005-213566. D. L. Ellis, "GRCop-84: A High-Temperature
       Copper Alloy for High-Heat-Flux Applications", NASA/TM-2005-213566,
       https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf
.. [5] NASA-CR-2000-210055. D. L. Ellis and D. J. Keller, "Thermophysical
       Properties of GRCop-84", NASA/CR-2000-210055,
       https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf
.. [6] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
       718", Publication SMC-045, 2007,
       https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
.. [7] AKSTEEL-304-PDB. AK Steel, "Type 304/304L Stainless Steel" product
       data bulletin, 2018, https://Protolabs.com/media/1020110/304-304l-stainless.pdf
.. [8] NIST-SP811-B8. NIST SP 811, Appendix B.8, factors for units listed
       alphabetically,
       https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
.. [9] NIST-SP811-B9. NIST SP 811, Appendix B.9, factors for units listed by
       kind of quantity,
       https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from fractions import Fraction
from types import MappingProxyType

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import real, strings, text

ACCESSED = "2026-10-25"

MATERIAL_SOURCES = MappingProxyType({
    "NIST-WEBBOOK-FLUID": (
        "E. W. Lemmon, M. O. McLinden and D. G. Friend, Thermophysical Properties of Fluid "
        "Systems, NIST Chemistry WebBook SRD 69", "https://webbook.nist.gov/chemistry/fluid/"),
    "NISTIR-6646-RP1": (
        "J. W. Magee et al., Thermophysical Properties Measurements and Models for Rocket "
        "Propellant RP-1: Phase I, NISTIR 6646, 2007",
        "https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf"),
    "AURUBIS-C10200": (
        "Aurubis, Technical Datasheet C10200 (Cu-OF), revision 18 08 US",
        "https://www.aurubis.be/en/dam/jcr:90029b00-b968-486e-bde0-ce56f16022cd/"
        "c10200-cu-of-us.pdf"),
    "NASA-TM-2005-213566": (
        "D. L. Ellis, GRCop-84: A High-Temperature Copper Alloy for High-Heat-Flux "
        "Applications, NASA/TM-2005-213566",
        "https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf"),
    "NASA-CR-2000-210055": (
        "D. L. Ellis and D. J. Keller, Thermophysical Properties of GRCop-84, "
        "NASA/CR-2000-210055",
        "https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf"),
    "SMC-045-INCONEL-718": (
        "Special Metals Corporation, INCONEL alloy 718, Publication SMC-045, 2007",
        "https://www.specialmetals.com/documents/technical-bulletins/inconel/"
        "inconel-alloy-718.pdf"),
    "AKSTEEL-304-PDB": (
        "AK Steel, Type 304/304L Stainless Steel product data bulletin, 2018",
        "https://Protolabs.com/media/1020110/304-304l-stainless.pdf"),
    "NIST-SP811-B8": (
        "NIST SP 811, Appendix B.8: Factors for Units Listed Alphabetically",
        "https://www.nist.gov/pml/special-publication-811/"
        "nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8"),
    "NIST-SP811-B9": (
        "NIST SP 811, Appendix B.9: Factors for Units Listed by Kind of Quantity",
        "https://www.nist.gov/pml/special-publication-811/"
        "nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9"),
})

_COOLANT_UNITS = {
    "density": "kg/m^3", "specific_heat": "J/(kg K)", "dynamic_viscosity": "Pa s",
    "kinematic_viscosity": "m^2/s", "thermal_conductivity": "W/(m K)",
    "critical_temperature": "K", "critical_pressure": "Pa",
}
_WALL_UNITS = {
    "thermal_conductivity": "W/(m K)", "density": "kg/m^3", "solidus_temperature": "K",
    "liquidus_temperature": "K", "maximum_service_temperature": "K",
}


@dataclass(frozen=True, slots=True)
class PropertyCitation:
    """A registered material source with a claim-specific locator.

    Attributes
    ----------
    source_id : str
        Key of :data:`MATERIAL_SOURCES`; the same key indexes
        ``docs/REFERENCES.md``.
    locator : str
        Table, page, row or query that states the value.

    Raises
    ------
    InputError
        If ``source_id`` is not registered (code ``"unknown_material_source"``).
    DomainError
        If the locator is empty or not text.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    """

    source_id: str
    locator: str

    def __post_init__(self) -> None:
        if not isinstance(self.source_id, str) or self.source_id not in MATERIAL_SOURCES:
            raise InputError("Unregistered material source.", code="unknown_material_source",
                             field="source_id")
        object.__setattr__(self, "locator", text(self.locator, "citation locator"))


@dataclass(frozen=True, slots=True)
class CitedProperty:
    """One published material property with its source and reference condition.

    Attributes
    ----------
    name : str
        Property name, for example ``"density"``.
    value : float
        Value in SI ``unit``. When ``upper_value`` is set, ``value`` is the
        lower end of a published range.
    unit : str
        SI unit symbol.
    citation : PropertyCitation
        Source and locator of the printed value.
    source_value : str
        The value exactly as printed in the source, with its unit.
    condition : str
        Reference condition in words (state, temper, sample, basis).
    temperature_k : float or None
        Reference temperature, K; ``None`` when the source does not state it.
    pressure_pa : float or None
        Reference pressure, Pa; ``None`` when not stated or not applicable.
    conversion : PropertyCitation or None
        Source of the unit-conversion factor applied to ``source_value``;
        ``None`` when only SI prefixes or the Celsius offset were applied.
    upper_value : float or None
        Upper end of a published range, in ``unit``.

    Raises
    ------
    DomainError
        If a value is not finite and positive, a reference condition is
        non-positive, or ``upper_value`` is below ``value``.

    Notes
    -----
    A property value is meaningful only at its stated condition; using it
    elsewhere is a modelling assumption that callers must declare.

    References
    ----------
    .. [1] NIST-SP811-B8. NIST SP 811, Appendix B.8, factors for units listed
           alphabetically,
           https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
    """

    name: str
    value: float
    unit: str
    citation: PropertyCitation
    source_value: str
    condition: str
    temperature_k: float | None = None
    pressure_pa: float | None = None
    conversion: PropertyCitation | None = None
    upper_value: float | None = None

    def __post_init__(self) -> None:
        for label in ("name", "unit", "source_value", "condition"):
            object.__setattr__(self, label, text(getattr(self, label), label))
        if not isinstance(self.citation, PropertyCitation):
            raise DomainError("A catalog property requires a typed citation.", field=self.name)
        if self.conversion is not None and not isinstance(self.conversion, PropertyCitation):
            raise DomainError("A conversion source must be a typed citation.", field=self.name)
        for label in ("value", "temperature_k", "pressure_pa", "upper_value"):
            item = getattr(self, label)
            if item is None and label != "value":
                continue
            item = real(item, f"{self.name}.{label}")
            if item <= 0:
                raise DomainError("Catalog values and conditions must be positive.",
                                  field=f"{self.name}.{label}")
            object.__setattr__(self, label, item)
        if self.upper_value is not None and self.upper_value < self.value:
            raise DomainError("Published range is reversed.", field=self.name)


def _check_record(record, units, label):
    """Validate text fields, per-field units and limitations of a catalog record."""
    for name in ("key", "name", "identity"):
        object.__setattr__(record, name, text(getattr(record, name), f"{label} {name}"))
    for name, unit in units.items():
        item = getattr(record, name)
        if item is None:
            continue
        if not isinstance(item, CitedProperty) or item.name != name or item.unit != unit:
            raise DomainError(f"{label} field must be a cited {unit} property.", field=name)
    limitations = strings(record.limitations, f"{label} limitations")
    if not limitations:
        raise DomainError(f"{label} records must state their limitations.", field="limitations")
    object.__setattr__(record, "limitations", limitations)


@dataclass(frozen=True, slots=True)
class CoolantRecord:
    """Reference-condition properties of one regenerative coolant.

    Attributes
    ----------
    key : str
        Catalog key, for example ``"liquid-methane"``.
    name : str
        Display name.
    identity : str
        Fluid identity as given by the source (for example the NIST fluid ID).
    density, specific_heat, dynamic_viscosity, thermal_conductivity : CitedProperty or None
        Liquid properties at the condition stated on each property.
    kinematic_viscosity : CitedProperty or None
        Published kinematic viscosity where the source gives it instead of
        the dynamic viscosity.
    critical_temperature, critical_pressure : CitedProperty or None
        Critical-point data of the source equation of state.
    limitations : tuple of str
        Evidence limitations of the record.

    Raises
    ------
    DomainError
        If a text field is empty, a property has the wrong name or unit, or
        no limitation is stated.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    .. [2] NISTIR-6646-RP1. J. W. Magee et al., "Thermophysical Properties
           Measurements and Models for Rocket Propellant RP-1: Phase I",
           NISTIR 6646, 2007,
           https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf
    """

    key: str
    name: str
    identity: str
    density: CitedProperty | None
    specific_heat: CitedProperty | None
    dynamic_viscosity: CitedProperty | None
    kinematic_viscosity: CitedProperty | None
    thermal_conductivity: CitedProperty | None
    critical_temperature: CitedProperty | None
    critical_pressure: CitedProperty | None
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        _check_record(self, _COOLANT_UNITS, "Coolant")

    def properties(self) -> tuple[CitedProperty, ...]:
        """Return every published property of the record, in field order.

        Returns
        -------
        tuple of CitedProperty
            Absent (unverified) properties are omitted, never zero-filled.

        References
        ----------
        .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G.
               Friend, "Thermophysical Properties of Fluid Systems", NIST
               Chemistry WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
        """
        return tuple(getattr(self, f.name) for f in fields(self)
                     if f.name in _COOLANT_UNITS and getattr(self, f.name) is not None)


@dataclass(frozen=True, slots=True)
class WallMaterialRecord:
    """Reference-condition properties and a sourced temperature limit of a wall alloy.

    Attributes
    ----------
    key : str
        Catalog key, for example ``"inconel-718"``.
    name : str
        Display name.
    identity : str
        Alloy designation (UNS number or source designation).
    thermal_conductivity : CitedProperty or None
        Thermal conductivity at its stated temperature; a published range
        stores its lower end in ``value``.
    density : CitedProperty or None
        Density as published.
    solidus_temperature, liquidus_temperature : CitedProperty or None
        Lower and upper end of the published melting range; a single
        published melting point is stored as the solidus.
    maximum_service_temperature : CitedProperty or None
        Sourced maximum service or design temperature; its ``condition``
        states the source's criterion.
    limitations : tuple of str
        Evidence limitations of the record.

    Raises
    ------
    DomainError
        If a text field is empty, a property has the wrong name or unit, or
        no limitation is stated.

    References
    ----------
    .. [1] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
           718", Publication SMC-045, 2007,
           https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
    .. [2] NASA-CR-2000-210055. D. L. Ellis and D. J. Keller, "Thermophysical
           Properties of GRCop-84", NASA/CR-2000-210055,
           https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf
    """

    key: str
    name: str
    identity: str
    thermal_conductivity: CitedProperty | None
    density: CitedProperty | None
    solidus_temperature: CitedProperty | None
    liquidus_temperature: CitedProperty | None
    maximum_service_temperature: CitedProperty | None
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        _check_record(self, _WALL_UNITS, "Wall material")

    def properties(self) -> tuple[CitedProperty, ...]:
        """Return every published property of the record, in field order.

        Returns
        -------
        tuple of CitedProperty
            Absent (unverified) properties are omitted, never zero-filled.

        References
        ----------
        .. [1] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
               718", Publication SMC-045, 2007,
               https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
        """
        return tuple(getattr(self, f.name) for f in fields(self)
                     if f.name in _WALL_UNITS and getattr(self, f.name) is not None)

    def temperature_limit(self) -> tuple[CitedProperty, str]:
        """Return the wall temperature limit and the basis it was taken from.

        Returns
        -------
        limit : CitedProperty
            The sourced maximum service temperature when the catalog has one,
            otherwise the solidus (or single melting point).
        basis : {"maximum_service_temperature", "solidus_temperature"}
            Which field supplied the limit. A solidus limit only bounds
            melting; it is not a strength or life limit.

        Raises
        ------
        DomainError
            If the record has neither field (code ``"wall_limit_unavailable"``).

        References
        ----------
        .. [1] NASA-CR-2000-210055. D. L. Ellis and D. J. Keller,
               "Thermophysical Properties of GRCop-84", NASA/CR-2000-210055,
               https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf
        """
        if self.maximum_service_temperature is not None:
            return self.maximum_service_temperature, "maximum_service_temperature"
        if self.solidus_temperature is not None:
            return self.solidus_temperature, "solidus_temperature"
        raise DomainError(f"{self.key} has no sourced temperature limit.",
                          code="wall_limit_unavailable", field="wall_material")


def _exact(text_value: str, factor: str = "1", offset: str = "0", divisor: str = "1") -> float:
    """Return ``(text_value + offset) * factor / divisor`` rounded once to binary64."""
    value = (Fraction(text_value) + Fraction(offset)) * Fraction(factor) / Fraction(divisor)
    return float(value)


def _celsius(text_value: str) -> float:
    """Convert a printed Celsius temperature with ``T/K = t/degC + 273.15`` [1]_.

    References
    ----------
    .. [1] NIST-SP811-B8. NIST SP 811, Appendix B.8, degree Celsius entry,
           https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
    """
    return _exact(text_value, offset="273.15")


def _fahrenheit(text_value: str) -> float:
    """Convert a printed Fahrenheit temperature with ``T/K = (t/degF + 459.67)/1.8`` [1]_.

    References
    ----------
    .. [1] NIST-SP811-B8. NIST SP 811, Appendix B.8, degree Fahrenheit entry,
           https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
    """
    return _exact(text_value, offset="459.67", divisor="1.8")


def _cite(source_id: str, locator: str) -> PropertyCitation:
    """Return a citation whose locator carries the access date of this review."""
    return PropertyCitation(source_id, f"{locator}; accessed {ACCESSED}")


_SP811_CELSIUS = PropertyCitation(
    "NIST-SP811-B8", f"degree Celsius: T/K = t/degC + 273.15; accessed {ACCESSED}")
_SP811_FAHRENHEIT = PropertyCitation(
    "NIST-SP811-B8", f"degree Fahrenheit: T/K = (t/degF + 459.67)/1.8; accessed {ACCESSED}")
_SP811_GCC = PropertyCitation("NIST-SP811-B8",
                              f"gram per cubic centimeter: 1.0 E+03 kg/m^3; accessed {ACCESSED}")
_SP811_LB_IN3 = PropertyCitation(
    "NIST-SP811-B9", f"Density: pound per cubic inch, 2.767 990 E+04 kg/m^3; accessed {ACCESSED}")
_SP811_BTU_IN = PropertyCitation(
    "NIST-SP811-B9", "Thermal conductivity: Btu_IT in/(h ft^2 degF), 1.442 279 E-01 W/(m K); "
    f"accessed {ACCESSED}")


def _nist_fluid(fluid_id, temperature, density, cp, viscosity, conductivity, critical_t,
                critical_p, condition):
    """Return the NIST WebBook liquid row and critical point as cited properties [1]_.

    The liquid row is an isobaric query at 1.0000 MPa returned with five
    digits (``Digits=5``) in K, MPa, kg/m3, J/(g K), Pa s and W/(m K); the
    critical point is the final row of the full saturation table (``Type=SatT``)
    of the same fluid. ``J/(g K)`` is multiplied by 1000 (SI prefix only).

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    """
    row = (f"fluid.cgi ID={fluid_id}, Type=IsoBar, P=1.0000 MPa, row T={temperature} K, "
           "liquid phase")
    sat = f"fluid.cgi ID={fluid_id}, Type=SatT full saturation table, final (critical) row"
    t, p = float(temperature), 1.0e6

    def liquid(name, value, unit, printed):
        return CitedProperty(name, value, unit, _cite("NIST-WEBBOOK-FLUID", row), printed,
                             condition, t, p)

    return dict(
        density=liquid("density", _exact(density), "kg/m^3", f"{density} kg/m3"),
        specific_heat=liquid("specific_heat", _exact(cp, "1000"), "J/(kg K)", f"{cp} J/g*K"),
        dynamic_viscosity=liquid("dynamic_viscosity", _exact(viscosity), "Pa s",
                                 f"{viscosity} Pa*s"),
        kinematic_viscosity=None,
        thermal_conductivity=liquid("thermal_conductivity", _exact(conductivity), "W/(m K)",
                                    f"{conductivity} W/m*K"),
        critical_temperature=CitedProperty(
            "critical_temperature", _exact(critical_t), "K", _cite("NIST-WEBBOOK-FLUID", sat),
            f"{critical_t} K", "critical point of the WebBook equation of state"),
        critical_pressure=CitedProperty(
            "critical_pressure", _exact(critical_p, "1000000"), "Pa",
            _cite("NIST-WEBBOOK-FLUID", sat), f"{critical_p} MPa",
            "critical point of the WebBook equation of state"),
    )


COOLANTS = (
    CoolantRecord(
        "liquid-hydrogen", "Liquid parahydrogen", "NIST WebBook fluid ID B5000001 (Parahydrogen)",
        **_nist_fluid("B5000001", "20.000", "72.292", "9.1868", "1.4779e-05", "0.10243",
                      "32.938", "1.2858", "compressed liquid parahydrogen"),
        limitations=(
            "Parahydrogen equation-of-state values; ortho/para composition of a real load "
            "is not represented.",
            "Single compressed-liquid point at 20 K and 1 MPa; properties vary strongly with "
            "temperature near the critical point.",
        ),
    ),
    CoolantRecord(
        "liquid-methane", "Liquid methane", "NIST WebBook fluid ID C74828 (Methane)",
        **_nist_fluid("C74828", "111.00", "424.16", "3.4660", "0.00011994", "0.18563",
                      "190.56", "4.5992", "compressed liquid methane"),
        limitations=(
            "Pure methane; liquefied-natural-gas impurities are not represented.",
            "Single compressed-liquid point at 111 K and 1 MPa.",
        ),
    ),
    CoolantRecord(
        "water", "Water", "NIST WebBook fluid ID C7732185 (Water)",
        **_nist_fluid("C7732185", "300.00", "996.96", "4.1781", "0.00085366", "0.61000",
                      "647.10", "22.064", "compressed liquid water"),
        limitations=("Single compressed-liquid point at 300 K and 1 MPa.",),
    ),
    CoolantRecord(
        "rp-1", "RP-1 kerosene", "RP-1, NIST original sample (NISTIR 6646)",
        density=CitedProperty(
            "density", _exact("799.01"), "kg/m^3",
            _cite("NISTIR-6646-RP1", "Table 8 (printed p. 35, PDF p. 41), original sample, "
                  "row 23.283 degC / 83.72 kPa"),
            "799.01 kg/m3", "liquid RP-1 (original sample) under nitrogen",
            _celsius("23.283"), _exact("83.72", "1000")),
        specific_heat=None,
        dynamic_viscosity=None,
        kinematic_viscosity=CitedProperty(
            "kinematic_viscosity", _exact("1.867", "1e-6"), "m^2/s",
            _cite("NISTIR-6646-RP1", "Table 11 (printed p. 61, PDF p. 67), original sample, "
                  "first 298.15 K row, capillary I"),
            "1.867 mm2/s", "liquid RP-1 (original sample) at atmospheric pressure "
            "(section 7.1); pressure value not stated in the table", _exact("298.15")),
        thermal_conductivity=CitedProperty(
            "thermal_conductivity", _exact("0.11122"), "W/(m K)",
            _cite("NISTIR-6646-RP1", "Table 10 (printed p. 48, PDF p. 54), point 2056, "
                  "radiation-corrected lambda_c at Te = 302.028 K, Pe = 0.2179 MPa"),
            "0.11122 W/(m K)", "liquid RP-1, transient hot wire, radiation corrected",
            _exact("302.028"), _exact("0.2179", "1000000")),
        critical_temperature=None,
        critical_pressure=None,
        limitations=(
            "Measured on specific NIST samples; RP-1 composition varies between lots.",
            "Density, viscosity and conductivity are at three different stated conditions.",
            "No primary heat capacity table was found in the source (unverified, absent).",
            "Dynamic viscosity is not published here; only the kinematic viscosity is stored.",
        ),
    ),
)

WALL_MATERIALS = (
    WallMaterialRecord(
        "ofhc-copper", "Oxygen-free copper (Cu-OF)", "UNS C10200 (Aurubis C102)",
        thermal_conductivity=CitedProperty(
            "thermal_conductivity", _exact("391"), "W/(m K)",
            _cite("AURUBIS-C10200", "p. 1, physical properties table, thermal conductivity "
                  "at 20 degC"),
            "391 W/mK (226 Btu/ft h degF)", "datasheet value; the sheet notes that the "
            "specified conductivity applies to the soft condition", _celsius("20")),
        density=CitedProperty(
            "density", _exact("8.90", "1000"), "kg/m^3",
            _cite("AURUBIS-C10200", "p. 1, physical properties table, density"),
            "8.90 g/cm3 (0.323 lb/in3)", "datasheet value; temperature not stated",
            conversion=_SP811_GCC),
        solidus_temperature=CitedProperty(
            "solidus_temperature", _celsius("1083"), "K",
            _cite("AURUBIS-C10200", "p. 1, physical properties table, melting point"),
            "1083 degC (1981 degF)", "single published melting point stored as solidus",
            conversion=_SP811_CELSIUS),
        liquidus_temperature=None,
        maximum_service_temperature=None,
        limitations=(
            "Manufacturer typical values, not specification minima.",
            "No maximum service temperature is published in the source; the limit check "
            "falls back to the melting point, which is not a strength limit.",
        ),
    ),
    WallMaterialRecord(
        "grcop-84", "GRCop-84 (Cu-8 at.% Cr-4 at.% Nb)", "NASA GRCop-84",
        thermal_conductivity=CitedProperty(
            "thermal_conductivity", _exact("305"), "W/(m K)",
            _cite("NASA-TM-2005-213566", "printed p. 8 (PDF p. 12), Thermal Conductivity: "
                  "between 305 and 320 W/m K"),
            "305 to 320 W/m K (176 to 185 BTU/h ft degF)",
            "range over the operating temperature range of an MCC liner; the temperature "
            "is not stated", upper_value=_exact("320")),
        density=None,
        solidus_temperature=None,
        liquidus_temperature=None,
        maximum_service_temperature=CitedProperty(
            "maximum_service_temperature", _exact("973"), "K",
            _cite("NASA-CR-2000-210055", "printed p. 13 (PDF p. 17): using a yield strength "
                  "of 100 MPa as the criterion, GRCop-84 can be used up to 973 K"),
            "973 K", "maximum operating temperature at a 100 MPa yield-strength criterion"),
        limitations=(
            "Conductivity is a published range without a stated temperature; the catalog "
            "value is its lower end (larger wall temperature drop).",
            "Density and solidus were not verified from a primary source and are absent.",
        ),
    ),
    WallMaterialRecord(
        "inconel-718", "INCONEL alloy 718", "UNS N07718 / W.Nr. 2.4668",
        thermal_conductivity=CitedProperty(
            "thermal_conductivity", _exact("79", "0.1442279"), "W/(m K)",
            _cite("SMC-045-INCONEL-718", "Table 5 Thermal Properties (PDF p. 2), 70 degF row, "
                  "annealed + aged column"),
            "79 BTU in/ft2 h degF", "annealed 1800 degF/1 hr and aged; the source notes "
            "conductivity was calculated from resistivity", _fahrenheit("70"),
            conversion=_SP811_BTU_IN),
        density=CitedProperty(
            "density", _exact("0.297", "27679.90"), "kg/m^3",
            _cite("SMC-045-INCONEL-718", "Table 2 Physical Constants (PDF p. 1), density, "
                  "annealed and aged"),
            "0.297 lb/in3", "annealed and aged; temperature not stated",
            conversion=_SP811_LB_IN3),
        solidus_temperature=CitedProperty(
            "solidus_temperature", _celsius("1260"), "K",
            _cite("SMC-045-INCONEL-718", "Table 2 Physical Constants (PDF p. 1), melting "
                  "range 1260-1336 degC, lower end"),
            "1260 degC (2300 degF)", "lower end of the published melting range",
            conversion=_SP811_CELSIUS),
        liquidus_temperature=CitedProperty(
            "liquidus_temperature", _celsius("1336"), "K",
            _cite("SMC-045-INCONEL-718", "Table 2 Physical Constants (PDF p. 1), melting "
                  "range 1260-1336 degC, upper end"),
            "1336 degC (2437 degF)", "upper end of the published melting range",
            conversion=_SP811_CELSIUS),
        maximum_service_temperature=CitedProperty(
            "maximum_service_temperature", _fahrenheit("1300"), "K",
            _cite("SMC-045-INCONEL-718", "p. 1 introduction: material used at -423 to "
                  "1300 degF"),
            "1300 degF", "upper end of the stated use range; not a coupled stress limit",
            conversion=_SP811_FAHRENHEIT),
        limitations=(
            "Typical values, not suitable for specification purposes (source statement).",
            "Btu is taken as the International Table Btu of the cited conversion factor.",
        ),
    ),
    WallMaterialRecord(
        "stainless-304", "Type 304/304L stainless steel", "AISI 304/304L (AK Steel)",
        thermal_conductivity=CitedProperty(
            "thermal_conductivity", _exact("16.2"), "W/(m K)",
            _cite("AKSTEEL-304-PDB", "PDF p. 4, Physical Properties, thermal conductivity "
                  "212 degF (100 degC) row"),
            "9.4 BTU/hr/ft/degF (16.2 W/(m K))", "product data bulletin typical value",
            _celsius("100")),
        density=CitedProperty(
            "density", _exact("8.03", "1000"), "kg/m^3",
            _cite("AKSTEEL-304-PDB", "PDF p. 4, Physical Properties, density"),
            "0.29 lbs/in3 (8.03 g/cm3)", "temperature not stated", conversion=_SP811_GCC),
        solidus_temperature=CitedProperty(
            "solidus_temperature", _celsius("1399"), "K",
            _cite("AKSTEEL-304-PDB", "PDF p. 4, Physical Properties, melting range "
                  "1399-1454 degC, lower end"),
            "1399 degC (2550 degF)", "lower end of the published melting range",
            conversion=_SP811_CELSIUS),
        liquidus_temperature=CitedProperty(
            "liquidus_temperature", _celsius("1454"), "K",
            _cite("AKSTEEL-304-PDB", "PDF p. 4, Physical Properties, melting range "
                  "1399-1454 degC, upper end"),
            "1454 degC (2650 degF)", "upper end of the published melting range",
            conversion=_SP811_CELSIUS),
        maximum_service_temperature=CitedProperty(
            "maximum_service_temperature", _celsius("899"), "K",
            _cite("AKSTEEL-304-PDB", "PDF p. 4, Oxidation Resistance: maximum temperature "
                  "for continuous exposure without appreciable scaling"),
            "1650 degF (899 degC)", "continuous exposure without appreciable scaling "
            "(oxidation criterion, not a strength limit)", conversion=_SP811_CELSIUS),
        limitations=(
            "Typical product-bulletin values for 304/304L; 347 is not covered by this source.",
            "The service limit is an oxidation-scaling criterion in air.",
        ),
    ),
)


def _normalized(key: object, label: str) -> str:
    """Return a lower-case, hyphenated catalog key."""
    if not isinstance(key, str):
        raise InputError(f"{label} key must be text.", field=label)
    return key.strip().lower().replace("_", "-").replace(" ", "-")


def coolant_catalog() -> tuple[CoolantRecord, ...]:
    """Return the immutable coolant catalog.

    Returns
    -------
    tuple of CoolantRecord
        Liquid parahydrogen, liquid methane, water and RP-1 records.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    .. [2] NISTIR-6646-RP1. J. W. Magee et al., "Thermophysical Properties
           Measurements and Models for Rocket Propellant RP-1: Phase I",
           NISTIR 6646, 2007,
           https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf
    """
    return COOLANTS


def get_coolant(key: str) -> CoolantRecord:
    """Resolve a coolant record by key.

    Parameters
    ----------
    key : str
        Catalog key; case, surrounding spaces and underscores are normalized
        (``"Liquid_Methane"`` resolves ``"liquid-methane"``).

    Returns
    -------
    CoolantRecord

    Raises
    ------
    InputError
        If the key is not text or not in the catalog (code
        ``"unknown_coolant"``).

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    """
    wanted = _normalized(key, "coolant")
    for record in COOLANTS:
        if record.key == wanted:
            return record
    raise InputError(f"Unknown coolant {key!r}; known: {', '.join(r.key for r in COOLANTS)}.",
                     code="unknown_coolant", field="coolant")


def wall_material_catalog() -> tuple[WallMaterialRecord, ...]:
    """Return the immutable chamber-wall material catalog.

    Returns
    -------
    tuple of WallMaterialRecord
        Oxygen-free copper, GRCop-84, Inconel 718 and 304 stainless records.

    References
    ----------
    .. [1] AURUBIS-C10200. Aurubis, "Technical Datasheet C10200 (Cu-OF)",
           revision 18 08 US,
           https://www.aurubis.be/en/dam/jcr:90029b00-b968-486e-bde0-ce56f16022cd/c10200-cu-of-us.pdf
    .. [2] NASA-TM-2005-213566. D. L. Ellis, "GRCop-84: A High-Temperature
           Copper Alloy for High-Heat-Flux Applications", NASA/TM-2005-213566,
           https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf
    .. [3] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
           718", Publication SMC-045, 2007,
           https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
    .. [4] AKSTEEL-304-PDB. AK Steel, "Type 304/304L Stainless Steel" product
           data bulletin, 2018, https://Protolabs.com/media/1020110/304-304l-stainless.pdf
    """
    return WALL_MATERIALS


def get_wall_material(key: str) -> WallMaterialRecord:
    """Resolve a chamber-wall material record by key.

    Parameters
    ----------
    key : str
        Catalog key; case, surrounding spaces and underscores are normalized.

    Returns
    -------
    WallMaterialRecord

    Raises
    ------
    InputError
        If the key is not text or not in the catalog (code
        ``"unknown_wall_material"``).

    References
    ----------
    .. [1] SMC-045-INCONEL-718. Special Metals Corporation, "INCONEL alloy
           718", Publication SMC-045, 2007,
           https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
    """
    wanted = _normalized(key, "wall_material")
    for record in WALL_MATERIALS:
        if record.key == wanted:
            return record
    raise InputError(
        f"Unknown wall material {key!r}; known: {', '.join(r.key for r in WALL_MATERIALS)}.",
        code="unknown_wall_material", field="wall_material")
