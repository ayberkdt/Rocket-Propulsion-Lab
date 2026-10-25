"""Centrifugal-pump similarity, suction and head-flow characteristic relations.

The module complements the ideal shaft-power balance of
:mod:`rocket_propulsion.propulsion.cycles`, which has no pump map or
cavitation model, with four preliminary-design relations:

* speed and impeller-trim scaling of flow, head and shaft power;
* specific speed ``Ns = N Q^(1/2) / H^(3/4)`` and suction specific speed
  ``Ss = N Q^(1/2) / NPSH^(3/4)`` in the U.S. customary units of the
  source monographs (rpm, U.S. gpm, ft) [1]_, [3]_;
* net positive suction head ``NPSH = (p_t - p_v) / (rho g)`` and a
  cavitation-margin check against a required NPSH [1]_, [3]_;
* a normalized quadratic, homologous head-flow characteristic
  ``H / H_d = A n^2 + B n q + C q^2`` (``q = Q/Q_d``, ``n = N/N_d``) [2]_.

These are similarity and bookkeeping relations. They do not predict a pump's
characteristic, efficiency or required NPSH; those are inputs. Fluid
compressibility and density change break the incompressible similarity
relations [1]_, and the results carry ``physical_model_qualified=False``.

References
----------
.. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
       Turbopumps", NASA SP-8109, December 1973,
       https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
.. [2] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
       equation (53), and section 5.4.4,
       https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
.. [3] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
       Pumping System Performance: A Sourcebook for Industry", 2nd ed., 2006,
       https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
.. [4] NIST-SP811-B9. NIST SP 811, Appendix B.9, factors for units listed by
       kind of quantity,
       https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

from .cycles import EngineCycleResult

#: U.S. gallon per minute in m^3/s (NIST SP 811 B.9, 6.309 020 E-05).
GPM_M3_S = 6.309020e-05
#: International foot in m (NIST SP 811 B.9, 3.048 E-01).
FOOT_M = 0.3048
#: Revolutions per minute per rad/s; one revolution is 2 pi rad.
RPM_PER_RAD_S = 60.0 / (2.0 * pi)
#: Lower and upper specific speeds of flight-proven centrifugal pumps, SP-8109 p. 7.
SP8109_FLIGHT_SPECIFIC_SPEED_RANGE = (450.0, 2100.0)
#: Impeller diameter ratio below which DOE notes trimming is rarely done.
DOE_MINIMUM_TRIM_RATIO = 0.7
_MAP_TOLERANCE = 1e-12


def _positive(value: object, label: str) -> float:
    """Return ``value`` as a finite positive float or raise :class:`DomainError`."""
    number = real(value, label)
    if number <= 0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _finite(value: float, label: str) -> float:
    """Reject a non-finite derived value."""
    if not isfinite(value):
        raise DomainError(f"{label} is outside the finite range.", field=label)
    return value


def pump_specific_speed(*, shaft_speed_rad_s: float, flow_m3_s: float, head_m: float) -> float:
    """Return the pump specific speed in U.S. customary units.

    Parameters
    ----------
    shaft_speed_rad_s : float
        Shaft speed, rad/s.
    flow_m3_s : float
        Volumetric flow, m^3/s.
    head_m : float
        Head rise per stage, m.

    Returns
    -------
    float
        ``Ns = N Q^(1/2) / H^(3/4)`` with ``N`` in rpm, ``Q`` in U.S. gpm and
        ``H`` in ft, the convention of [1]_ and [2]_.

    Raises
    ------
    DomainError
        If an input is not finite and positive.

    Notes
    -----
    SP-8109 equation (1) gives the definition and reports flight-proven
    centrifugal pumps between 450 and 2100 [1]_. The SI inputs are
    converted with ``1 gpm = 6.309020e-5 m^3/s`` and ``1 ft = 0.3048 m``
    [3]_. The value is not dimensionless and is not comparable with
    specific speeds in other unit systems.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 6 (PDF p. 18),
           equation (1); printed p. 7 (PDF p. 19),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [2] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
           Pumping System Performance", 2nd ed., 2006, printed p. 39,
           https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
    .. [3] NIST-SP811-B9. NIST SP 811, Appendix B.9, factors for units
           listed by kind of quantity,
           https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9
    """
    speed = _positive(shaft_speed_rad_s, "shaft_speed_rad_s") * RPM_PER_RAD_S
    flow = _positive(flow_m3_s, "flow_m3_s") / GPM_M3_S
    head = _positive(head_m, "head_m") / FOOT_M
    return _finite(speed * flow ** 0.5 / head ** 0.75, "Specific speed")


def pump_suction_specific_speed(*, shaft_speed_rad_s: float, flow_m3_s: float,
                                npsh_m: float) -> float:
    """Return the suction specific speed in U.S. customary units.

    Parameters
    ----------
    shaft_speed_rad_s : float
        Shaft speed, rad/s.
    flow_m3_s : float
        Volumetric flow, m^3/s.
    npsh_m : float
        Net positive suction head, m; usually the required NPSH.

    Returns
    -------
    float
        ``Ss = N Q^(1/2) / NPSH^(3/4)`` with ``N`` in rpm, ``Q`` in U.S. gpm
        and ``NPSH`` in ft [1]_.

    Raises
    ------
    DomainError
        If an input is not finite and positive.

    Notes
    -----
    SP-8109 reports suction specific speed capabilities from 15 000 to over
    40 000 for rocket pumps on fluids with vapor pressures similar to cold
    water, and notes that the pumped fluid strongly affects suction
    performance [1]_.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 11 (PDF p. 23),
           equation (8),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [2] NIST-SP811-B9. NIST SP 811, Appendix B.9, factors for units
           listed by kind of quantity,
           https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9
    """
    speed = _positive(shaft_speed_rad_s, "shaft_speed_rad_s") * RPM_PER_RAD_S
    flow = _positive(flow_m3_s, "flow_m3_s") / GPM_M3_S
    npsh = _positive(npsh_m, "npsh_m") / FOOT_M
    return _finite(speed * flow ** 0.5 / npsh ** 0.75, "Suction specific speed")


def net_positive_suction_head(*, inlet_total_pressure_pa: float, vapor_pressure_pa: float,
                              density_kg_m3: float) -> float:
    """Return the available net positive suction head at the pump inlet.

    Parameters
    ----------
    inlet_total_pressure_pa : float
        Pump-inlet total pressure, Pa; positive.
    vapor_pressure_pa : float
        Liquid vapor pressure at the inlet temperature, Pa; non-negative.
    density_kg_m3 : float
        Liquid density, kg/m^3; positive.

    Returns
    -------
    float
        ``(p_t - p_v) / (rho g_n)``, m. A negative value means the inlet
        total pressure is below the vapor pressure.

    Raises
    ------
    DomainError
        If an input is out of domain.

    Notes
    -----
    SP-8109 defines ``NPSH = (P_t - P_v) / rho`` with ``rho`` a weight
    density in its customary units [1]_; with a mass density the head is
    divided by ``g``. The conventional ``g_n`` of the SI Brochure is used
    [2]_.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, Appendix A glossary, printed
           p. 88 (PDF p. 100),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [2] BIPM-SI. BIPM, "SI Brochure, 9th edition", conventional standard
           gravity, https://www.bipm.org/en/publications/si-brochure
    """
    total = _positive(inlet_total_pressure_pa, "inlet_total_pressure_pa")
    vapor = real(vapor_pressure_pa, "vapor_pressure_pa")
    if vapor < 0:
        raise DomainError("Vapor pressure must be non-negative.", field="vapor_pressure_pa")
    density = _positive(density_kg_m3, "density_kg_m3")
    return _finite((total - vapor) / (density * STANDARD_GRAVITY_M_S2), "NPSH")


@dataclass(frozen=True, slots=True)
class PumpAffinityResult:
    """Pump operating point scaled by speed and impeller-trim relations.

    Attributes
    ----------
    flow_m3_s, head_m, shaft_power_w : float
        Scaled flow (m^3/s), head (m) and shaft power (W).
    speed_ratio, diameter_ratio : float
        New-to-reference shaft speed and impeller diameter ratios.
    warnings : tuple of WarningMessage
        Similarity-scope warnings.
    physical_model_qualified : bool
        Always False.

    References
    ----------
    .. [1] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 5.4.4,
           printed p. 110 (PDF p. 150),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
           Pumping System Performance", 2nd ed., 2006, printed pp. 49-50,
           https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
    """

    flow_m3_s: float
    head_m: float
    shaft_power_w: float
    speed_ratio: float
    diameter_ratio: float
    warnings: tuple[WarningMessage, ...]
    physical_model_qualified: bool = False


def scale_pump_affinity(*, flow_m3_s: float, head_m: float, shaft_power_w: float,
                        speed_ratio: float, diameter_ratio: float = 1.0) -> PumpAffinityResult:
    """Scale a pump operating point to a new shaft speed and trimmed impeller diameter.

    Parameters
    ----------
    flow_m3_s : float
        Reference volumetric flow, m^3/s.
    head_m : float
        Reference head rise, m.
    shaft_power_w : float
        Reference shaft power, W.
    speed_ratio : float
        ``N2 / N1``; positive.
    diameter_ratio : float, default 1.0
        Trimmed-to-original impeller diameter ``D2 / D1``; positive.

    Returns
    -------
    PumpAffinityResult
        ``Q2 = Q1 n d``, ``H2 = H1 n^2 d^2`` and ``P2 = P1 n^3 d^3`` with
        ``n = speed_ratio`` and ``d = diameter_ratio``.

    Raises
    ------
    DomainError
        If an input is not finite and positive or an output overflows.

    Notes
    -----
    Speed scaling is the homologous form used by GFSSP's turbopump option,
    where ``Head/Speed^2`` and ``Torque/(Density Speed^2)`` are functions of
    ``flow rate/speed`` [1]_; at a homologous point and fixed density this
    gives ``Q ~ N``, ``H ~ N^2`` and shaft power ``torque x speed ~ N^3``.
    Diameter scaling is the constant-speed impeller-trimming relation
    ``Q ~ D``, ``H ~ D^2``, ``bhp ~ D^3`` [2]_, which the source calls an
    approximation for small changes; trims are rarely below 70 % of the
    original diameter [2]_. The trimming relation is not the geometric
    similarity of differently sized pumps. SP-8109 cautions that density
    and compressibility changes alter the incompressible affinity
    relationships [3]_.

    References
    ----------
    .. [1] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 5.4.4,
           printed p. 110 (PDF p. 150),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
           Pumping System Performance", 2nd ed., 2006, printed pp. 49-50
           (PDF pp. 53-54),
           https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
    .. [3] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 83 (PDF p. 95),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    """
    flow = _positive(flow_m3_s, "flow_m3_s")
    head = _positive(head_m, "head_m")
    power = _positive(shaft_power_w, "shaft_power_w")
    n = _positive(speed_ratio, "speed_ratio")
    d = _positive(diameter_ratio, "diameter_ratio")
    warnings = [WarningMessage(
        "pump_similarity_incompressible",
        "Ideal incompressible similarity at fixed density and efficiency; density and "
        "compressibility changes alter the relations (SP-8109).", "info")]
    if d > 1:
        warnings.append(WarningMessage(
            "impeller_enlargement_outside_trim_relation",
            "The trimming relation describes diameter reduction at constant speed; an "
            "enlarged impeller is outside its basis.", "warning", "diameter_ratio"))
    elif d < DOE_MINIMUM_TRIM_RATIO:
        warnings.append(WarningMessage(
            "impeller_trim_below_doe_practice",
            "Impeller diameters are rarely reduced below 70 % of the original (DOE "
            "sourcebook); the trimming relation is an approximation for small changes.",
            "warning", "diameter_ratio"))
    return PumpAffinityResult(
        flow_m3_s=_finite(flow * n * d, "Scaled flow"),
        head_m=_finite(head * (n * d) ** 2, "Scaled head"),
        shaft_power_w=_finite(power * (n * d) ** 3, "Scaled shaft power"),
        speed_ratio=n, diameter_ratio=d, warnings=tuple(warnings))


@dataclass(frozen=True, slots=True)
class PumpSuctionResult:
    """Available versus required NPSH with a cavitation-margin decision.

    Attributes
    ----------
    npsh_available_m : float
        ``(p_t - p_v) / (rho g_n)``, m.
    npsh_required_m : float
        Required NPSH supplied by the caller (pump test or vendor data), m.
    npsh_margin_m : float
        Available minus required NPSH, m.
    npsh_margin_ratio : float
        Available divided by required NPSH.
    required_margin_ratio : float
        Ratio the available NPSH must reach.
    cavitation_margin_adequate : bool
        ``npsh_available >= required_margin_ratio * npsh_required``.
    suction_specific_speed : float or None
        ``Ss`` at the required NPSH in U.S. customary units, when speed and
        flow are given.
    warnings : tuple of WarningMessage
        Margin and scope warnings.
    physical_model_qualified : bool
        Always False.

    References
    ----------
    .. [1] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
           Pumping System Performance", 2nd ed., 2006, printed p. 38,
           https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
    """

    npsh_available_m: float
    npsh_required_m: float
    npsh_margin_m: float
    npsh_margin_ratio: float
    required_margin_ratio: float
    cavitation_margin_adequate: bool
    suction_specific_speed: float | None
    warnings: tuple[WarningMessage, ...]
    physical_model_qualified: bool = False


def evaluate_pump_suction(
    *, inlet_total_pressure_pa: float, vapor_pressure_pa: float, density_kg_m3: float,
    npsh_required_m: float, required_margin_ratio: float = 1.0,
    shaft_speed_rad_s: float | None = None, flow_m3_s: float | None = None,
) -> PumpSuctionResult:
    """Compare available and required NPSH and flag an insufficient cavitation margin.

    Parameters
    ----------
    inlet_total_pressure_pa : float
        Pump-inlet total pressure, Pa.
    vapor_pressure_pa : float
        Liquid vapor pressure at the inlet temperature, Pa.
    density_kg_m3 : float
        Liquid density, kg/m^3.
    npsh_required_m : float
        Required NPSH of the pump at the operating flow, m; positive. Its
        head-loss criterion is the caller's (SP-8109 plots pump data at 2 %
        and inducer data at 10 % head loss [1]_).
    required_margin_ratio : float, default 1.0
        Minimum ``NPSHA / NPSHR``; at least 1. The default applies only the
        bare criterion ``NPSHA >= NPSHR``. A DOE design rule for industrial
        systems uses 1.25 [2]_; it is not applied by default.
    shaft_speed_rad_s, flow_m3_s : float, optional
        Given together, the suction specific speed at ``npsh_required_m``
        is reported.

    Returns
    -------
    PumpSuctionResult

    Raises
    ------
    DomainError
        If an input is out of domain or the margin ratio is below one.
    InputError
        If only one of speed and flow is given.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 11 (PDF p. 23) and
           Appendix A, printed p. 88 (PDF p. 100),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [2] DOE-PUMP-SOURCEBOOK. U.S. DOE and Hydraulic Institute, "Improving
           Pumping System Performance", 2nd ed., 2006, printed p. 38
           (PDF p. 42),
           https://www.energy.gov/sites/prod/files/2014/05/f16/pump.pdf
    """
    available = net_positive_suction_head(
        inlet_total_pressure_pa=inlet_total_pressure_pa, vapor_pressure_pa=vapor_pressure_pa,
        density_kg_m3=density_kg_m3)
    required = _positive(npsh_required_m, "npsh_required_m")
    ratio_required = real(required_margin_ratio, "required_margin_ratio")
    if ratio_required < 1:
        raise DomainError("Required NPSH margin ratio must be at least one.",
                          field="required_margin_ratio")
    if (shaft_speed_rad_s is None) != (flow_m3_s is None):
        raise InputError("Give shaft speed and flow together for suction specific speed.",
                         field="shaft_speed_rad_s")
    suction = None if flow_m3_s is None else pump_suction_specific_speed(
        shaft_speed_rad_s=shaft_speed_rad_s, flow_m3_s=flow_m3_s, npsh_m=required)
    adequate = available >= ratio_required * required
    warnings = [WarningMessage(
        "npsh_required_is_input",
        "Required NPSH is caller data; thermodynamic suppression head and fluid effects on "
        "suction performance are not modelled.", "info", "npsh_required_m")]
    if available <= 0:
        warnings.append(WarningMessage(
            "inlet_below_vapor_pressure",
            "Inlet total pressure does not exceed the vapor pressure.", "critical",
            "inlet_total_pressure_pa"))
    if not adequate:
        warnings.append(WarningMessage(
            "insufficient_npsh_margin",
            "Available NPSH is below the required margin; cavitation-induced head loss is "
            "expected.", "critical", "npsh_required_m"))
    return PumpSuctionResult(
        npsh_available_m=available, npsh_required_m=required,
        npsh_margin_m=_finite(available - required, "NPSH margin"),
        npsh_margin_ratio=_finite(available / required, "NPSH margin ratio"),
        required_margin_ratio=ratio_required, cavitation_margin_adequate=adequate,
        suction_specific_speed=suction, warnings=tuple(warnings))


@dataclass(frozen=True, slots=True)
class PumpCharacteristic:
    """Normalized homologous quadratic head-flow characteristic of a pump.

    ``H / H_d = A n^2 + B n q + C q^2`` with ``q = Q / Q_d`` and
    ``n = N / N_d``. Passing through the design point requires
    ``A + B + C = 1``.

    Attributes
    ----------
    design_flow_m3_s, design_head_m, design_speed_rad_s : float
        Design (reference) flow, head and shaft speed.
    shutoff_coefficient : float
        ``A``: normalized shut-off head at design speed.
    linear_coefficient : float
        ``B``.
    quadratic_coefficient : float
        ``C``.

    Raises
    ------
    DomainError
        If a design quantity is not finite and positive, a coefficient is
        not finite, ``A`` is not positive, or ``|A + B + C - 1|`` exceeds a
        1e-12 numerical tolerance.

    Notes
    -----
    GFSSP represents a pump by the curve fit ``dp = A0 + B0 m + C0 m^2``
    of pressure rise against mass flow [1]_ and, in its turbopump option,
    by ``Head/Speed^2`` as a function of ``flow rate/speed`` [2]_. Writing
    that homologous function as the same quadratic and normalizing by the
    design point gives the form above. The coefficients come from pump test
    or vendor data; none are supplied by this library.

    References
    ----------
    .. [1] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
           printed p. 51 (PDF p. 91), equation (53),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 5.4.4,
           printed p. 110 (PDF p. 150),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    """

    design_flow_m3_s: float
    design_head_m: float
    design_speed_rad_s: float
    shutoff_coefficient: float
    linear_coefficient: float
    quadratic_coefficient: float

    def __post_init__(self) -> None:
        for name in ("design_flow_m3_s", "design_head_m", "design_speed_rad_s"):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        for name in ("shutoff_coefficient", "linear_coefficient", "quadratic_coefficient"):
            object.__setattr__(self, name, real(getattr(self, name), name))
        if self.shutoff_coefficient <= 0:
            raise DomainError("Shut-off head coefficient must be positive.",
                              field="shutoff_coefficient")
        total = self.shutoff_coefficient + self.linear_coefficient + self.quadratic_coefficient
        if abs(total - 1.0) > _MAP_TOLERANCE:
            raise DomainError("Characteristic must pass through the design point (A+B+C=1).",
                              code="characteristic_not_normalized", field="shutoff_coefficient")

    def head_m(self, flow_m3_s: float, shaft_speed_rad_s: float) -> float:
        """Return the head rise at a flow and shaft speed.

        Parameters
        ----------
        flow_m3_s : float
            Volumetric flow, m^3/s; non-negative.
        shaft_speed_rad_s : float
            Shaft speed, rad/s; positive.

        Returns
        -------
        float
            ``H_d (A n^2 + B n q + C q^2)``, m; may be negative outside the
            characteristic's useful range.

        Raises
        ------
        DomainError
            If an input is out of domain.

        References
        ----------
        .. [1] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
               equation (53), and section 5.4.4,
               https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
        """
        q, n = self._normalized(flow_m3_s, shaft_speed_rad_s)
        return _finite(self.design_head_m * (self.shutoff_coefficient * n * n
                       + self.linear_coefficient * n * q + self.quadratic_coefficient * q * q),
                       "Characteristic head")

    def _normalized(self, flow_m3_s: float, shaft_speed_rad_s: float) -> tuple[float, float]:
        """Return ``(Q/Q_d, N/N_d)`` after domain checks."""
        flow = real(flow_m3_s, "flow_m3_s")
        if flow < 0:
            raise DomainError("Flow must be non-negative.", field="flow_m3_s")
        speed = _positive(shaft_speed_rad_s, "shaft_speed_rad_s")
        return flow / self.design_flow_m3_s, speed / self.design_speed_rad_s


@dataclass(frozen=True, slots=True)
class PumpMapPoint:
    """One evaluated point of a :class:`PumpCharacteristic`.

    Attributes
    ----------
    characteristic : PumpCharacteristic
        Evaluated characteristic.
    flow_m3_s, shaft_speed_rad_s : float
        Operating flow (m^3/s) and speed (rad/s).
    normalized_flow, normalized_speed, normalized_head : float
        ``q``, ``n`` and ``H / H_d``.
    head_m : float
        Head rise, m.
    head_slope_s_m2 : float
        ``dH/dQ`` at constant speed, m per (m^3/s).
    zero_slope_flow_m3_s : float or None
        Flow of zero slope at this speed, ``-B n Q_d / (2 C)``; ``None`` when
        ``C = 0``.
    negative_slope : bool
        Whether head decreases with increasing flow at the point.
    meets_sp8109_zero_slope_rule : bool or None
        For ``C < 0``: whether the zero-slope flow is at least 10 % below the
        operating flow (SP-8109 recommendation, applied to this point as the
        lowest required flow). ``None`` when ``C >= 0`` (no head maximum).
    specific_speed : float or None
        U.S. customary specific speed at the point; ``None`` when flow or
        head is not positive.
    warnings : tuple of WarningMessage
        Stability, range and scope warnings.
    physical_model_qualified : bool
        Always False.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 67 (PDF p. 79),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    """

    characteristic: PumpCharacteristic
    flow_m3_s: float
    shaft_speed_rad_s: float
    normalized_flow: float
    normalized_speed: float
    normalized_head: float
    head_m: float
    head_slope_s_m2: float
    zero_slope_flow_m3_s: float | None
    negative_slope: bool
    meets_sp8109_zero_slope_rule: bool | None
    specific_speed: float | None
    warnings: tuple[WarningMessage, ...]
    physical_model_qualified: bool = False


def evaluate_pump_map(characteristic: PumpCharacteristic, *, flow_m3_s: float,
                      shaft_speed_rad_s: float) -> PumpMapPoint:
    """Evaluate head, slope and stability indicators of a pump characteristic.

    Parameters
    ----------
    characteristic : PumpCharacteristic
        Normalized homologous characteristic.
    flow_m3_s : float
        Operating flow, m^3/s; non-negative.
    shaft_speed_rad_s : float
        Operating shaft speed, rad/s; positive.

    Returns
    -------
    PumpMapPoint

    Raises
    ------
    DomainError
        If ``characteristic`` is not a :class:`PumpCharacteristic` or an
        input is out of domain.

    Notes
    -----
    SP-8109 recommends a characteristic of decreasing head rise with
    increasing flow, with the zero-slope point of the H-Q curve at least
    10 % below the lowest required flow and the stall point at least 15 %
    below it [1]_; the stall point is not represented by the quadratic and
    is not checked. Flight-proven centrifugal pumps span specific speeds of
    450 to 2100 [2]_; points outside give an information warning.

    References
    ----------
    .. [1] NASA-SP8109-PUMP. NASA, "Liquid Rocket Engine Centrifugal Flow
           Turbopumps", NASA SP-8109, 1973, printed p. 67 (PDF p. 79),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [2] NASA-SP8109-PUMP. NASA SP-8109, printed p. 7 (PDF p. 19),
           https://ntrs.nasa.gov/api/citations/19740020848/downloads/19740020848.pdf
    .. [3] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
           equation (53), and section 5.4.4,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    """
    if not isinstance(characteristic, PumpCharacteristic):
        raise DomainError("A PumpCharacteristic is required.", field="characteristic")
    c = characteristic
    q, n = c._normalized(flow_m3_s, shaft_speed_rad_s)
    head = c.head_m(flow_m3_s, shaft_speed_rad_s)
    scale = c.design_head_m / c.design_flow_m3_s
    slope = _finite(scale * (c.linear_coefficient * n + 2 * c.quadratic_coefficient * q),
                    "Head slope")
    zero = None
    if c.quadratic_coefficient != 0:
        zero = _finite(-c.linear_coefficient * n * c.design_flow_m3_s
                       / (2 * c.quadratic_coefficient), "Zero-slope flow")
    rule = None if c.quadratic_coefficient >= 0 else zero <= 0.9 * float(flow_m3_s)
    warnings = [WarningMessage(
        "pump_map_is_caller_data",
        "Characteristic coefficients are caller data; efficiency, stall and cavitation are "
        "not represented by the quadratic.", "info")]
    specific = None
    if head > 0 and flow_m3_s > 0:
        specific = pump_specific_speed(shaft_speed_rad_s=shaft_speed_rad_s,
                                       flow_m3_s=flow_m3_s, head_m=head)
        low, high = SP8109_FLIGHT_SPECIFIC_SPEED_RANGE
        if not low <= specific <= high:
            warnings.append(WarningMessage(
                "specific_speed_outside_sp8109_flight_range",
                "Specific speed is outside the 450-2100 range of flight-proven centrifugal "
                "pumps reported by SP-8109.", "info"))
    if head <= 0:
        warnings.append(WarningMessage(
            "pump_map_nonpositive_head",
            "The characteristic gives no positive head at this point; it is outside the "
            "useful map.", "critical", "flow_m3_s"))
    if slope >= 0:
        warnings.append(WarningMessage(
            "pump_head_slope_not_negative",
            "Head does not decrease with flow here; SP-8109 recommends a decreasing "
            "head-flow characteristic for stable operation.", "warning", "flow_m3_s"))
    if rule is False:
        warnings.append(WarningMessage(
            "pump_zero_slope_margin_below_sp8109",
            "Zero-slope flow is less than 10 % below this flow (SP-8109 recommendation).",
            "warning", "flow_m3_s"))
    return PumpMapPoint(
        characteristic=c, flow_m3_s=float(flow_m3_s), shaft_speed_rad_s=float(shaft_speed_rad_s),
        normalized_flow=q, normalized_speed=n, normalized_head=head / c.design_head_m,
        head_m=head, head_slope_s_m2=slope, zero_slope_flow_m3_s=zero,
        negative_slope=slope < 0, meets_sp8109_zero_slope_rule=rule, specific_speed=specific,
        warnings=tuple(warnings))


@dataclass(frozen=True, slots=True)
class CyclePumpOperatingPoint:
    """Main-pump duty of an engine-cycle result placed on a pump characteristic.

    Attributes
    ----------
    required_head_m : float
        Cycle pump pressure rise divided by ``rho g_n``, m.
    volumetric_flow_m3_s : float
        Cycle pump mass flow divided by density, m^3/s.
    map_point : PumpMapPoint
        Characteristic evaluated at that flow and the given speed.
    head_margin_m : float
        Characteristic head minus required head, m.
    head_sufficient : bool
        Whether the characteristic head meets the cycle requirement.
    warnings : tuple of WarningMessage
        Margin and scope warnings (the map point keeps its own).
    physical_model_qualified : bool
        Always False.

    References
    ----------
    .. [1] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
           equation (53),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    """

    required_head_m: float
    volumetric_flow_m3_s: float
    map_point: PumpMapPoint
    head_margin_m: float
    head_sufficient: bool
    warnings: tuple[WarningMessage, ...]
    physical_model_qualified: bool = False


def evaluate_cycle_pump_point(cycle: EngineCycleResult, characteristic: PumpCharacteristic, *,
                              shaft_speed_rad_s: float) -> CyclePumpOperatingPoint:
    """Check whether a pump characteristic delivers a pumped cycle's head at its flow.

    Parameters
    ----------
    cycle : EngineCycleResult
        Result of :func:`~rocket_propulsion.propulsion.cycles.calculate_engine_cycle`
        for a gas-generator or staged-combustion cycle.
    characteristic : PumpCharacteristic
        Main-pump characteristic.
    shaft_speed_rad_s : float
        Pump shaft speed, rad/s.

    Returns
    -------
    CyclePumpOperatingPoint

    Raises
    ------
    InputError
        If ``cycle`` is not an :class:`EngineCycleResult` or is pressure-fed
        (code ``"cycle_has_no_pump"``).
    DomainError
        If an input is out of domain.

    Notes
    -----
    The cycle's single equivalent density converts its pump pressure rise
    and mass flow to head and volumetric flow, with the conventional
    ``g_n`` [1]_. The cycle balance itself is not changed: the shaft power
    and turbine flow still follow its constant-efficiency balance, not the
    characteristic.

    References
    ----------
    .. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition", conventional standard
           gravity, https://www.bipm.org/en/publications/si-brochure
    .. [2] NASA-GFSSP6-PUMP. NASA GFSSP v6 user manual, section 3.1.7.14,
           equation (53),
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    """
    if not isinstance(cycle, EngineCycleResult):
        raise InputError("An EngineCycleResult is required.", field="cycle")
    if cycle.cycle == "pressure-fed" or cycle.pump_mass_flow_kg_s <= 0:
        raise InputError("The cycle has no pump duty to place on a characteristic.",
                         code="cycle_has_no_pump", field="cycle")
    density = cycle.inputs.propellant_density_kg_m3
    required = _finite(cycle.pump_pressure_rise_pa / (density * STANDARD_GRAVITY_M_S2),
                       "Required head")
    flow = _finite(cycle.pump_mass_flow_kg_s / density, "Volumetric flow")
    point = evaluate_pump_map(characteristic, flow_m3_s=flow, shaft_speed_rad_s=shaft_speed_rad_s)
    margin = _finite(point.head_m - required, "Head margin")
    warnings = [WarningMessage(
        "cycle_pump_point_not_coupled",
        "Characteristic head is compared with the cycle requirement only; the cycle shaft "
        "balance is not re-solved on the characteristic.", "info")]
    if margin < 0:
        warnings.append(WarningMessage(
            "pump_head_below_cycle_requirement",
            "The characteristic head at the cycle flow and given speed is below the required "
            "pump head.", "critical", "shaft_speed_rad_s"))
    return CyclePumpOperatingPoint(required_head_m=required, volumetric_flow_m3_s=flow,
                                   map_point=point, head_margin_m=margin,
                                   head_sufficient=margin >= 0, warnings=tuple(warnings))
