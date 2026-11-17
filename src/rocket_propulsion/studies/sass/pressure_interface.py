"""Screen reviewed feed outlet conditions against exact thruster inlet references."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import text

from .catalog import CatalogField, ThrusterRecord
from .evidence import CitedScalar, quantity
from .rational_numerics import outward_endpoint

InletPort = Literal["inlet", "fuel_inlet", "oxidizer_inlet"]

# NBSIR 82-2401, Appendix 4, table A4: exact standard gravity times the
# avoirdupois pound in kilograms, divided by one square inch in square metres.
_PA_PER_PSI = Fraction("9.80665") * Fraction("0.45359237") / Fraction("0.0254") ** 2


@dataclass(frozen=True, slots=True)
class InletPressureWindow:
    """Pressure bounds at one thruster port for the declared operating window.

    For gauge pressure, the optional reference interval is the absolute
    pressure at the reference side of the measurement, over the same window.
    The reference record identifies the reviewed pairing and pressure datum.

    References
    ----------
    .. [1] BIPM-SI, pressure in Pa.
    .. [2] NASA-SE-6, 6.3 and 6.6, interface and data traceability.
    """

    port: InletPort
    reference: Literal["absolute", "gauge", "differential"]
    lower: CitedScalar
    upper: CitedScalar
    condition_record: str
    gauge_reference_lower: CitedScalar | None = None
    gauge_reference_upper: CitedScalar | None = None
    gauge_reference_record: str | None = None

    def __post_init__(self) -> None:
        """Require named ports and complete, coherent pressure bounds.

        References
        ----------
        .. [1] NASA-SE-6, 6.3 and 6.6, interface definition and data scope.
        """
        text(self.condition_record, "inlet pressure condition")
        if self.port not in ("inlet", "fuel_inlet", "oxidizer_inlet"):
            raise DomainError("Pressure must be declared at a thruster inlet port.")
        if self.reference not in ("absolute", "gauge", "differential"):
            raise DomainError("An inlet pressure window requires an explicit pressure reference.")
        low = quantity(self.lower, "Pa")
        high = quantity(self.upper, "Pa")
        if low > high or (self.reference == "absolute" and low < 0):
            raise DomainError("Inlet pressure bounds are reversed or below absolute zero pressure.")
        reference_fields = (
            self.gauge_reference_lower, self.gauge_reference_upper, self.gauge_reference_record
        )
        if any(value is not None for value in reference_fields):
            if (self.reference != "gauge" or self.gauge_reference_lower is None
                    or self.gauge_reference_upper is None or self.gauge_reference_record is None):
                raise DomainError("Gauge conversion requires both absolute reference bounds and their review record.")
            text(self.gauge_reference_record, "gauge reference pressure review")
            ref_low = quantity(self.gauge_reference_lower, "Pa", minimum=0)
            ref_high = quantity(self.gauge_reference_upper, "Pa", minimum=0)
            if ref_low > ref_high or Fraction(low) + Fraction(ref_low) < 0:
                raise DomainError("Gauge reference bounds are reversed or imply negative absolute pressure.")


@dataclass(frozen=True, slots=True)
class CatalogueMeopCheck:
    """Comparison with the published component MEOP for this pressure port.

    References
    ----------
    .. [1] MOOG-COLD, published component pressure ratings by model.
    .. [2] NASA-SE-6, 6.3, component and system interfaces.
    """

    status: Literal["within_meop_reference", "crosses_meop_reference", "exceeds_meop_reference",
                    "meop_unit_disagreement", "unassessed"]
    published_limit: CatalogField
    reason: str


@dataclass(frozen=True, slots=True)
class CataloguePressureCheck:
    """A comparison of one submitted interface window with published evidence.

    References
    ----------
    .. [1] NASA-SE-6, 6.3 and 6.6, scope and pedigree of interface comparisons.
    """

    port: InletPort
    status: Literal["matches_nominal", "within_published_range", "outside_reference",
                    "overlaps_reference", "unassessed"]
    absolute_window_pa: tuple[float, float] | None
    published: CatalogField | None
    reason: str
    meop: CatalogueMeopCheck | None = None


def inlet_ports(record: ThrusterRecord) -> tuple[InletPort, ...]:
    """Identify the separate propellant inlet ports for a catalogue record.

    References
    ----------
    .. [1] ARIANE-10N, separate fuel and oxidizer flow calibration and propellants.
    .. [2] MOOG-COLD, 058E151 inlet.
    """
    return ("fuel_inlet", "oxidizer_inlet") if record.family == "bipropellant" else ("inlet",)


def _absolute_window(window: InletPressureWindow) -> tuple[Fraction, Fraction] | None:
    """Convert a reviewed gauge window by adding the two interval endpoints.

    References
    ----------
    .. [1] BIPM-SI, retain Pa after the reference conversion.
    .. [2] NASA-SE-6, 6.6, preserve the documented measurement reference.
    """
    low, high = Fraction(window.lower.value), Fraction(window.upper.value)
    if window.reference == "absolute":
        return low, high
    if window.reference == "gauge" and window.gauge_reference_lower is not None:
        assert window.gauge_reference_upper is not None
        return low + Fraction(window.gauge_reference_lower.value), high + Fraction(window.gauge_reference_upper.value)
    return None


def _meop_check(
    fields: tuple[CatalogField, ...], absolute: tuple[Fraction, Fraction] | None,
) -> CatalogueMeopCheck | None:
    """Compare the inlet window against every expression of the operating rating.

    References
    ----------
    .. [1] MOOG-COLD, MEOP, proof and burst rows for the same hardware.
    .. [2] NASA-SE-6, 6.3 and 6.6, scope and pedigree of reference data.
    """
    limit = next((field for field in fields if field.name == "maximum_expected_operating_pressure"), None)
    if limit is None:
        return None
    if absolute is None or limit.pressure_reference != "absolute" or limit.kind != "upper_bound":
        return CatalogueMeopCheck("unassessed", limit, "No comparable absolute inlet window and operating rating.")
    low, high = absolute
    metric = Fraction(quantity(limit.upper, "Pa", strictly_positive=True))
    alternate = limit.reported_psia
    imperial = metric if alternate is None else Fraction(alternate.value) * _PA_PER_PSI
    lower_limit, upper_limit = sorted((metric, imperial))
    if low > upper_limit:
        return CatalogueMeopCheck("exceeds_meop_reference", limit, "The entire inlet window exceeds published MEOP.")
    if high <= lower_limit:
        return CatalogueMeopCheck("within_meop_reference", limit, "The submitted inlet window is within published MEOP.")
    if low > lower_limit and high <= upper_limit:
        return CatalogueMeopCheck("meop_unit_disagreement", limit,
                                  "Only one printed MEOP value admits the entire inlet window.")
    return CatalogueMeopCheck("crosses_meop_reference", limit,
                              "The inlet window crosses at least one printed MEOP value.")


def _inlet_check(
    record: ThrusterRecord, port: InletPort, window: InletPressureWindow | None,
) -> CataloguePressureCheck:
    """Compare one compatible port and pressure reference without extrapolation.

    References
    ----------
    .. [1] NASA-SE-6, 6.3, interface requirements and assumptions.
    .. [2] MOOG-COLD, published nominal inlet condition.
    """
    fields = tuple(field for field in record.fields if field.pressure_port == port)
    ratings = ("maximum_expected_operating_pressure", "proof_pressure", "burst_pressure")
    references = tuple(field for field in fields if field.name not in ratings)
    source = references[0] if len(references) == 1 else None
    absolute = None if window is None else _absolute_window(window)
    meop = _meop_check(fields, absolute)
    if window is None:
        return CataloguePressureCheck(port, "unassessed", None, source, "Inlet window not supplied.", meop)
    if absolute is None:
        return CataloguePressureCheck(port, "unassessed", None, source,
                                      "No reviewed conversion to absolute inlet pressure.", meop)
    low, high = absolute
    reported = (outward_endpoint(low, upward=False), outward_endpoint(high, upward=True))
    if (source is None or source.pressure_reference != "absolute" or source.lower is None
            or source.kind not in ("nominal", "range")):
        return CataloguePressureCheck(port, "unassessed", reported, source,
                                      "No comparable published absolute pressure for this inlet port.", meop)
    if source.kind == "nominal":
        pair = record.nominal_pair()
        if (pair is None or source.nominal_point_id is None
                or source.nominal_point_id != record.field("thrust").nominal_point_id
                or not set.intersection(*({c.source_id for c in value.citations}
                                          for value in (*pair, source.lower, source.upper)))):
            return CataloguePressureCheck(port, "unassessed", reported, source,
                                          "Published pressure is not linked to the nominal performance point.", meop)
    ref_low = Fraction(quantity(source.lower, "Pa", minimum=0))
    ref_high = Fraction(quantity(source.upper, "Pa", minimum=0))
    if high < ref_low or low > ref_high:
        return CataloguePressureCheck(port, "outside_reference", reported, source,
                                      "The submitted absolute window is disjoint from the published reference.", meop)
    if low < ref_low or high > ref_high:
        return CataloguePressureCheck(port, "overlaps_reference", reported, source,
                                      "Part of the submitted window extends beyond the published reference.", meop)
    return CataloguePressureCheck(
        port, "matches_nominal" if source.kind == "nominal" else "within_published_range",
        reported, source, "Pressure reference comparison for preliminary sizing.", meop
    )


def _catalogue_pressure_checks(
    record: ThrusterRecord, windows: tuple[InletPressureWindow, ...],
) -> tuple[CataloguePressureCheck, ...]:
    """Check submitted ports and retain missing companion inlets for bipropellant hardware.

    References
    ----------
    .. [1] NASA-SE-6, 6.3, complete interface coverage.
    .. [2] ARIANE-10N, separately adjusted propellant flow paths.
    """
    if not windows:
        return ()
    by_port = {window.port: window for window in windows}
    return tuple(_inlet_check(record, port, by_port.get(port)) for port in inlet_ports(record))

