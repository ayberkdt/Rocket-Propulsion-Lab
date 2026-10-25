"""Slender-body normal force, viscous cross-flow and Barrowman fin relations.

The module describes a pointed, axisymmetric launch-vehicle body as a chain
of continuous segments (cone or tangent-ogive nose, cylinders, conical
shoulders and boattails) and evaluates

* the slender-body (Munk) normal-force slope and centre of pressure in the
  form derived by Barrowman [1]_ (equations 3-57 to 3-89);
* the Allen-Perkins viscous cross-flow correction for inclined bodies of
  revolution [2]_ (the cross-force expression on printed page 8 and
  equations 14 and 15);
* Barrowman's subsonic normal-force slope and centre of pressure of a set of
  three or four trapezoidal fins [1]_ (equations 3-12 and 3-30), and the
  slender-body fin-body interference factors of Pitts, Nielsen and
  Kaattari [3]_ (equations 14 and 21).

Validity
--------
Barrowman's body and fin relations assume an angle of attack near zero,
steady irrotational flow, a rigid body and a sharp nose tip [1]_ (section
2.40). Allen and Perkins compared their method with tests of two
ogive-cylinders of fineness ratio 21.1 and 13.1 at Mach 1.98 up to about
20 degrees angle of attack [2]_ (figures 7 and 8); the cross-flow drag
coefficient is an input that the user must supply (the report's appendix
uses 1.2 for a circular-cylinder section at low cross-flow Mach number times
a finite-length ratio of 0.755 or 0.692).

References
----------
.. [1] BARROWMAN-1967. J. S. Barrowman, "The Practical Calculation of the
       Aerodynamic Characteristics of Slender Finned Vehicles", M.S. thesis,
       Catholic University of America, 1967; NTRS 20010047838.
.. [2] NACA-TR-1048. H. J. Allen and E. W. Perkins, "A Study of Effects of
       Viscosity on Flow over Slender Inclined Bodies of Revolution", NACA
       Report 1048, 1951.
.. [3] NACA-TR-1307. W. C. Pitts, J. N. Nielsen and G. E. Kaattari, "Lift
       and Center of Pressure of Wing-Body-Tail Combinations at Subsonic,
       Transonic, and Supersonic Speeds", NACA Report 1307, 1959.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import atan, cos, isfinite, pi, radians, sin, sqrt

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

_SEGMENT_KINDS = ("cone", "tangent_ogive", "cylinder", "frustum")
_NOSE_KINDS = ("cone", "tangent_ogive")
_SIMPSON_INTERVALS = 2000
_CONTINUITY_TOLERANCE = 1.0e-9


def _require_positive(value: float, field: str) -> float:
    """Return ``value`` as float if finite and strictly positive.

    Raises
    ------
    InputError
        If the value is not finite or not strictly positive.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise InputError(f"{field} must be a number.", field=field)
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise InputError(f"{field} must be finite and positive.", field=field)
    return number


def _require_non_negative(value: float, field: str) -> float:
    """Return ``value`` as float if finite and not negative.

    Raises
    ------
    InputError
        If the value is not finite or negative.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise InputError(f"{field} must be a number.", field=field)
    number = float(value)
    if not isfinite(number) or number < 0.0:
        raise InputError(f"{field} must be finite and non-negative.", field=field)
    return number


@dataclass(frozen=True, slots=True)
class BodySegment:
    """One axisymmetric body component of Barrowman's figure 3-7 [1]_.

    Attributes
    ----------
    kind : {"cone", "tangent_ogive", "cylinder", "frustum"}
        Component shape. ``"frustum"`` covers both a conical shoulder
        (radius growing aft) and a conical boattail (radius shrinking aft).
    length_m : float
        Axial length of the component, m.
    fore_radius_m, aft_radius_m : float
        Radii at the forward and aft faces, m. Noses start at zero radius.

    Notes
    -----
    A tangent ogive is the circular arc through the tip that is tangent to the
    cylinder of radius ``aft_radius_m`` at the shoulder; its generating radius
    is ``rho = (L^2 + R^2) / (2 R)`` and its local radius
    ``r(x) = sqrt(rho^2 - (L - x)^2) + R - rho``. This construction reproduces
    the body dimensions of Allen and Perkins [2]_ (figure 6: generating radius
    25.000 in, nose length 4.314 in, base radius 0.375 in).

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           figure 3-7 (body component shapes).
    .. [2] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951,
           figure 6.
    """

    kind: str
    length_m: float
    fore_radius_m: float
    aft_radius_m: float

    def __post_init__(self) -> None:
        if self.kind not in _SEGMENT_KINDS:
            raise InputError(
                f"Segment kind must be one of {', '.join(_SEGMENT_KINDS)}.", field="kind"
            )
        _require_positive(self.length_m, "length_m")
        _require_non_negative(self.fore_radius_m, "fore_radius_m")
        _require_non_negative(self.aft_radius_m, "aft_radius_m")
        if self.kind in _NOSE_KINDS:
            if self.fore_radius_m != 0.0:
                raise InputError("A nose segment must start at zero radius.", field="fore_radius_m")
            if self.aft_radius_m <= 0.0:
                raise InputError("A nose segment needs a positive base radius.", field="aft_radius_m")
        elif self.kind == "cylinder":
            if self.fore_radius_m != self.aft_radius_m or self.aft_radius_m <= 0.0:
                raise InputError(
                    "A cylinder needs equal, positive fore and aft radii.", field="aft_radius_m"
                )
        elif self.fore_radius_m <= 0.0 or self.aft_radius_m <= 0.0:
            raise InputError("A frustum needs positive fore and aft radii.", field="fore_radius_m")

    def radius_at(self, local_x_m: float) -> float:
        """Return the local radius at ``local_x_m`` from the segment front, m.

        Parameters
        ----------
        local_x_m : float
            Axial position measured from the forward face, ``0 <= x <= L``.

        Returns
        -------
        float
            Local body radius, m.

        Raises
        ------
        DomainError
            If ``local_x_m`` lies outside the segment.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
               figure 3-7.
        """

        x = float(local_x_m)
        if not isfinite(x) or x < -1.0e-12 or x > self.length_m * (1.0 + 1.0e-12):
            raise DomainError("Position lies outside the segment.", field="local_x_m")
        x = min(max(x, 0.0), self.length_m)
        if self.kind == "tangent_ogive":
            length = self.length_m
            base = self.aft_radius_m
            rho = (length * length + base * base) / (2.0 * base)
            return max(sqrt(max(rho * rho - (length - x) ** 2, 0.0)) + base - rho, 0.0)
        return self.fore_radius_m + (self.aft_radius_m - self.fore_radius_m) * x / self.length_m

    def radius_slope_at(self, local_x_m: float) -> float:
        """Return ``dr/dx`` at ``local_x_m`` from the segment front.

        Parameters
        ----------
        local_x_m : float
            Axial position measured from the forward face, ``0 <= x <= L``.

        Returns
        -------
        float
            Radius slope (dimensionless).

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
               figure 3-7.
        """

        x = min(max(float(local_x_m), 0.0), self.length_m)
        if self.kind == "tangent_ogive":
            length = self.length_m
            base = self.aft_radius_m
            rho = (length * length + base * base) / (2.0 * base)
            return (length - x) / sqrt(rho * rho - (length - x) ** 2)
        return (self.aft_radius_m - self.fore_radius_m) / self.length_m


def cone_nose(length_m: float, base_radius_m: float) -> BodySegment:
    """Return a conical nose segment of Barrowman's figure 3-7(a) [1]_.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    return BodySegment("cone", length_m, 0.0, base_radius_m)


def tangent_ogive_nose(length_m: float, base_radius_m: float) -> BodySegment:
    """Return a tangent-ogive nose segment of Barrowman's figure 3-7(a) [1]_.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    return BodySegment("tangent_ogive", length_m, 0.0, base_radius_m)


def cylinder_section(length_m: float, radius_m: float) -> BodySegment:
    """Return a cylindrical section of Barrowman's figure 3-7(b) [1]_.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    return BodySegment("cylinder", length_m, radius_m, radius_m)


def conical_frustum(length_m: float, fore_radius_m: float, aft_radius_m: float) -> BodySegment:
    """Return a conical shoulder or boattail of Barrowman's figure 3-7(c) [1]_.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    return BodySegment("frustum", length_m, fore_radius_m, aft_radius_m)


def _simpson(function, lower: float, upper: float, intervals: int = _SIMPSON_INTERVALS) -> float:
    """Composite Simpson integral of ``function`` on ``[lower, upper]``.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    step = (upper - lower) / intervals
    total = function(lower) + function(upper)
    for index in range(1, intervals):
        total += (4.0 if index % 2 else 2.0) * function(lower + index * step)
    return total * step / 3.0


@dataclass(frozen=True, slots=True)
class SegmentGeometry:
    """Integrated geometry of one body segment.

    Attributes
    ----------
    start_x_m, end_x_m : float
        Axial stations of the segment faces measured from the nose tip, m.
    fore_area_m2, aft_area_m2 : float
        Cross-section areas ``A(x0)`` and ``A(x1)``, m^2.
    volume_m3 : float
        Enclosed volume, m^3.
    planform_area_m2 : float
        Projected planform area ``integral 2 r dx``, m^2.
    planform_moment_m3 : float
        First moment of planform area about the nose tip ``integral 2 r x dx``,
        m^3.
    wetted_area_m2 : float
        Lateral surface area ``integral 2 pi r sqrt(1 + r'^2) dx``, m^2.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    start_x_m: float
    end_x_m: float
    fore_area_m2: float
    aft_area_m2: float
    volume_m3: float
    planform_area_m2: float
    planform_moment_m3: float
    wetted_area_m2: float


def _segment_geometry(segment: BodySegment, start_x: float) -> SegmentGeometry:
    """Integrate one segment; linear profiles are integrated in closed form.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    length = segment.length_m
    r0 = segment.fore_radius_m
    r1 = segment.aft_radius_m
    if segment.kind == "tangent_ogive":
        volume = _simpson(lambda x: pi * segment.radius_at(x) ** 2, 0.0, length)
        planform = _simpson(lambda x: 2.0 * segment.radius_at(x), 0.0, length)
        moment_local = _simpson(lambda x: 2.0 * segment.radius_at(x) * x, 0.0, length)
        wetted = _simpson(
            lambda x: 2.0 * pi * segment.radius_at(x) * sqrt(1.0 + segment.radius_slope_at(x) ** 2),
            0.0,
            length,
        )
    else:
        volume = pi * length * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0
        planform = length * (r0 + r1)
        moment_local = length * length * (r0 + 2.0 * r1) / 3.0
        wetted = pi * (r0 + r1) * sqrt(length * length + (r1 - r0) ** 2)
    return SegmentGeometry(
        start_x_m=start_x,
        end_x_m=start_x + length,
        fore_area_m2=pi * r0 * r0,
        aft_area_m2=pi * r1 * r1,
        volume_m3=volume,
        planform_area_m2=planform,
        planform_moment_m3=moment_local + start_x * planform,
        wetted_area_m2=wetted,
    )


@dataclass(frozen=True, slots=True)
class AxisymmetricBody:
    """A pointed body of revolution built from continuous segments [1]_.

    Attributes
    ----------
    segments : tuple of BodySegment
        Segments from the nose tip aft. The first must be a cone or tangent
        ogive; only the first may be a nose; adjacent radii must match
        (Barrowman requires no discontinuities between components, section
        3.21 [1]_).

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           sections 2.40 and 3.21.
    """

    segments: tuple[BodySegment, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.segments, tuple) or not self.segments:
            raise InputError("A body needs a non-empty tuple of segments.", field="segments")
        for index, segment in enumerate(self.segments):
            if not isinstance(segment, BodySegment):
                raise InputError("Every segment must be a BodySegment.", field="segments")
            if (segment.kind in _NOSE_KINDS) != (index == 0):
                raise InputError(
                    "The first segment, and only the first, must be a cone or tangent-ogive nose.",
                    field="segments",
                )
            if index:
                previous = self.segments[index - 1].aft_radius_m
                if abs(previous - segment.fore_radius_m) > _CONTINUITY_TOLERANCE * max(previous, 1.0):
                    raise InputError(
                        f"Radius jumps between segments {index - 1} and {index}.", field="segments"
                    )

    @property
    def length_m(self) -> float:
        """Total body length, m.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
        """

        return sum(segment.length_m for segment in self.segments)

    def _locate(self, x_m: float) -> tuple[BodySegment, float]:
        """Return the segment containing ``x_m`` and the local coordinate.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
        """

        x = float(x_m)
        if not isfinite(x) or x < 0.0 or x > self.length_m * (1.0 + 1.0e-12):
            raise DomainError("Station lies outside the body.", field="x_m")
        start = 0.0
        for segment in self.segments:
            if x <= start + segment.length_m or segment is self.segments[-1]:
                return segment, min(x - start, segment.length_m)
            start += segment.length_m
        raise DomainError("Station lies outside the body.", field="x_m")  # pragma: no cover

    def radius_at(self, x_m: float) -> float:
        """Return the body radius at station ``x_m`` from the nose tip, m.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
        """

        segment, local = self._locate(x_m)
        return segment.radius_at(local)

    def area_slope_at(self, x_m: float) -> float:
        """Return ``dA/dx = 2 pi r dr/dx`` at station ``x_m``, m.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
               equation 3-61.
        """

        segment, local = self._locate(x_m)
        return 2.0 * pi * segment.radius_at(local) * segment.radius_slope_at(local)


@dataclass(frozen=True, slots=True)
class BodyGeometry:
    """Integrated geometry of an :class:`AxisymmetricBody`.

    Attributes
    ----------
    length_m : float
        Total length ``L``, m.
    base_area_m2 : float
        Base cross-section area ``S_b = A(L)``, m^2.
    volume_m3 : float
        Body volume ``Q``, m^3.
    planform_area_m2 : float
        Planform area ``A_p``, m^2.
    planform_centroid_m : float
        Distance from the nose tip to the planform-area centroid, m.
    wetted_area_m2 : float
        Lateral (wetted) surface area excluding the base, m^2.
    fineness_ratio : float
        ``L`` divided by the maximum diameter.
    segments : tuple of SegmentGeometry
        Per-segment integrals.

    References
    ----------
    .. [1] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951,
           symbols A_p, S_b, Q, L of equation 14.
    """

    length_m: float
    base_area_m2: float
    volume_m3: float
    planform_area_m2: float
    planform_centroid_m: float
    wetted_area_m2: float
    fineness_ratio: float
    segments: tuple[SegmentGeometry, ...]


def body_geometry(body: AxisymmetricBody) -> BodyGeometry:
    """Integrate the area, volume, planform and wetted area of a body.

    Parameters
    ----------
    body : AxisymmetricBody
        The body to integrate.

    Returns
    -------
    BodyGeometry
        Integrated geometry.

    Notes
    -----
    Cones, frustums and cylinders are integrated in closed form; tangent
    ogives with composite Simpson quadrature on 2000 intervals. Barrowman's
    closed-form tangent-ogive volume (equation 3-90 [1]_) is not used: with
    the normalisation read as ``d_N A(L)`` its right-hand side equals one
    quarter of the integrated volume for every fineness ratio checked, i.e.
    it matches ``V = pi d_N^3 (RHS)``; see the module test.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equation 3-90.
    .. [2] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951.
    """

    if not isinstance(body, AxisymmetricBody):
        raise InputError("body must be an AxisymmetricBody.", field="body")
    start = 0.0
    pieces = []
    for segment in body.segments:
        pieces.append(_segment_geometry(segment, start))
        start += segment.length_m
    planform = sum(piece.planform_area_m2 for piece in pieces)
    max_radius = max(max(seg.fore_radius_m, seg.aft_radius_m) for seg in body.segments)
    return BodyGeometry(
        length_m=start,
        base_area_m2=pieces[-1].aft_area_m2,
        volume_m3=sum(piece.volume_m3 for piece in pieces),
        planform_area_m2=planform,
        planform_centroid_m=sum(piece.planform_moment_m3 for piece in pieces) / planform,
        wetted_area_m2=sum(piece.wetted_area_m2 for piece in pieces),
        fineness_ratio=start / (2.0 * max_radius),
        segments=tuple(pieces),
    )


@dataclass(frozen=True, slots=True)
class NormalForceComponent:
    """Normal-force slope and centre of pressure of one vehicle component.

    Attributes
    ----------
    name : str
        Component label.
    normal_force_slope_per_rad : float
        ``C_N_alpha`` referred to the vehicle reference area, 1/rad.
    center_of_pressure_m : float or None
        Distance from the nose tip to the component centre of pressure, m;
        ``None`` when the component carries no normal force (constant-area
        sections in slender-body theory).

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-106 and 3-107.
    """

    name: str
    normal_force_slope_per_rad: float
    center_of_pressure_m: float | None


def slender_body_normal_force(
    body: AxisymmetricBody, *, reference_area_m2: float
) -> tuple[NormalForceComponent, ...]:
    """Return Barrowman's subsonic slender-body ``C_N_alpha`` and CP per segment.

    Parameters
    ----------
    body : AxisymmetricBody
        Pointed body of revolution.
    reference_area_m2 : float
        Reference area ``A_r``, m^2.

    Returns
    -------
    tuple of NormalForceComponent
        One entry per segment. Combine them with
        :func:`combine_normal_force_components`.

    Notes
    -----
    Barrowman [1]_ integrates the slender-body running load
    ``n(x) = rho V^2 alpha dA/dx`` (equations 3-57 to 3-61) to obtain
    ``C_N_alpha = 2 [A(x1) - A(x0)] / A_r`` per radian (equation 3-65) and,
    integrating ``x dA/dx`` by parts (equations 3-83 to 3-85), the pitching
    moment about the nose tip ``(2 / (A_r L_r)) [x1 A(x1) - x0 A(x0) - V]``.
    Their ratio gives the component centre of pressure
    ``X = [x1 A(x1) - x0 A(x0) - V] / [A(x1) - A(x0)]``; for a whole body
    with ``A(0) = 0`` this is equation 3-89, ``X_B = l_0 - V_B / A_B``.
    Constant-area sections carry no normal force. Barrowman notes that the
    analysis ignores compressibility and that ``C_N_alpha`` rises near Mach
    one.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           section 3.21-3.22, equations 3-57 to 3-89.
    """

    reference_area = _require_positive(reference_area_m2, "reference_area_m2")
    geometry = body_geometry(body)
    components = []
    for index, (segment, piece) in enumerate(zip(body.segments, geometry.segments, strict=True)):
        delta_area = piece.aft_area_m2 - piece.fore_area_m2
        slope = 2.0 * delta_area / reference_area
        if abs(delta_area) <= 1.0e-15 * max(piece.aft_area_m2, 1.0e-300):
            components.append(NormalForceComponent(f"{index}:{segment.kind}", 0.0, None))
            continue
        moment = (
            piece.end_x_m * piece.aft_area_m2
            - piece.start_x_m * piece.fore_area_m2
            - piece.volume_m3
        )
        components.append(
            NormalForceComponent(f"{index}:{segment.kind}", slope, moment / delta_area)
        )
    return tuple(components)


@dataclass(frozen=True, slots=True)
class CombinedNormalForce:
    """Total normal-force slope and centre of pressure of a vehicle.

    Attributes
    ----------
    normal_force_slope_per_rad : float
        Sum of the component slopes, 1/rad (Barrowman equation 3-106).
    center_of_pressure_m : float
        Slope-weighted centre of pressure from the nose tip, m (equation
        3-107).
    components : tuple of NormalForceComponent
        The inputs.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-106 and 3-107.
    """

    normal_force_slope_per_rad: float
    center_of_pressure_m: float
    components: tuple[NormalForceComponent, ...]


def combine_normal_force_components(
    components: Sequence[NormalForceComponent],
) -> CombinedNormalForce:
    """Sum component normal-force slopes and take the moment-weighted CP.

    Parameters
    ----------
    components : sequence of NormalForceComponent
        Components referred to the same reference area.

    Returns
    -------
    CombinedNormalForce
        ``C_N_alpha = sum(C_N_alpha_i)`` and
        ``X = sum(X_i C_N_alpha_i) / C_N_alpha``.

    Raises
    ------
    DomainError
        If the total slope is zero, so the centre of pressure is undefined.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           section 3.40, equations 3-106 and 3-107.
    """

    items = tuple(components)
    if not items or not all(isinstance(item, NormalForceComponent) for item in items):
        raise InputError("components must be NormalForceComponent objects.", field="components")
    total = 0.0
    moment = 0.0
    for item in items:
        if not isfinite(item.normal_force_slope_per_rad):
            raise InputError("Component slopes must be finite.", field="components")
        if item.normal_force_slope_per_rad == 0.0:
            continue
        if item.center_of_pressure_m is None or not isfinite(item.center_of_pressure_m):
            raise InputError("A loaded component needs a finite CP.", field="components")
        total += item.normal_force_slope_per_rad
        moment += item.normal_force_slope_per_rad * item.center_of_pressure_m
    if total == 0.0:
        raise DomainError("Total normal-force slope is zero; the CP is undefined.")
    return CombinedNormalForce(total, moment / total, items)


@dataclass(frozen=True, slots=True)
class ViscousCrossflowResult:
    """Allen-Perkins force and moment coefficients of an inclined body [1]_.

    Attributes
    ----------
    angle_of_attack_deg : float
        Angle of attack, degrees.
    normal_force_coefficient : float
        Force normal to the body axis divided by ``q A`` (integral of the
        cross force ``f`` of printed page 8 [1]_).
    lift_coefficient : float
        ``C_L`` of equation 14 (or 15).
    foredrag_increment_coefficient : float
        ``C_Dp - C_Dp(alpha=0) cos^3(alpha)`` of equation 14 (or the
        small-angle ``Delta C_Dp`` of equation 15).
    pitching_moment_coefficient : float
        ``C_M`` about ``x_m`` (positive nose-up), equation 14 (or 15).
    center_of_pressure_m : float
        ``x_m - C_M X / C_N``, distance from the nose tip, m.
    potential_normal_force_coefficient : float
        The potential (Munk) part of ``normal_force_coefficient``.
    viscous_normal_force_coefficient : float
        The viscous cross-flow part of ``normal_force_coefficient``.
    small_angle : bool
        True when equation 15 (sines replaced by angles) was used.
    warnings : tuple of WarningMessage
        Validity notes.

    References
    ----------
    .. [1] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951.
    """

    angle_of_attack_deg: float
    normal_force_coefficient: float
    lift_coefficient: float
    foredrag_increment_coefficient: float
    pitching_moment_coefficient: float
    center_of_pressure_m: float
    potential_normal_force_coefficient: float
    viscous_normal_force_coefficient: float
    small_angle: bool
    warnings: tuple[WarningMessage, ...]


def viscous_crossflow_coefficients(
    body: AxisymmetricBody,
    *,
    angle_of_attack_deg: float,
    crossflow_drag_coefficient: float,
    reference_area_m2: float,
    reference_length_m: float,
    moment_reference_m: float = 0.0,
    small_angle: bool = False,
) -> ViscousCrossflowResult:
    """Evaluate the Allen-Perkins potential-plus-viscous cross-flow model.

    Parameters
    ----------
    body : AxisymmetricBody
        Pointed body of revolution.
    angle_of_attack_deg : float
        Angle of attack in degrees, ``0 <= alpha < 90``.
    crossflow_drag_coefficient : float
        ``C_d(alpha=90)``: the planform-averaged drag coefficient of the body
        in cross flow, already including any finite-length correction. The
        report's appendix [1]_ estimates it as a circular-cylinder section
        coefficient at the cross-flow Mach number ``M_c = M sin(alpha)`` and
        Reynolds number (1.2 at low ``M_c``) times the finite-length ratio
        (0.755 and 0.692 for its two bodies).
    reference_area_m2 : float
        Coefficient reference area ``A``, m^2.
    reference_length_m : float
        Moment reference length ``X``, m.
    moment_reference_m : float, default 0.0
        Moment centre ``x_m`` from the nose tip, m.
    small_angle : bool, default False
        Use equation 15 (cosines replaced by one and sines by the angle).

    Returns
    -------
    ViscousCrossflowResult
        Coefficients and centre of pressure.

    Notes
    -----
    The cross force per unit length normal to the axis is [1]_ (printed
    p. 8) ``f = q (dS/dx) sin(2 alpha) cos(alpha/2) + 2 R c_d q sin^2(alpha)``:
    Munk's potential cross force (inclined midway between the normals to the
    axis and to the wind, footnote 7) plus the drag of a cylinder element in
    the cross flow ``V sin(alpha)``. Integrating with a uniform ``c_d`` (so
    that ``x(alpha=90)`` is the planform centroid, as the appendix assumes)
    gives equation 14::

        C_L  = (S_b/A) sin2a cos(a/2) + C_d (A_p/A) sin^2 a cos a
        dC_Dp = (S_b/A) sin2a sin(a/2) + C_d (A_p/A) sin^3 a
        C_M  = [(Q - S_b (L - x_m)) / (A X)] sin2a cos(a/2)
               + C_d (A_p/A) ((x_m - x_p)/X) sin^2 a

    and ``C_N = (S_b/A) sin2a cos(a/2) + C_d (A_p/A) sin^2 a``. The
    small-angle forms are equation 15. Allen and Perkins validated the method
    for slender bodies of fineness ratio 13.1 and 21.1 at Mach 1.98 and
    found it adequate for very high fineness ratio, less so for the lower.

    References
    ----------
    .. [1] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951,
           printed pp. 8-9, equations 14 and 15, appendix (p. 13).
    """

    alpha_deg = float(angle_of_attack_deg)
    if not isfinite(alpha_deg) or not 0.0 <= alpha_deg < 90.0:
        raise DomainError("Angle of attack must lie in [0, 90) degrees.", field="angle_of_attack_deg")
    crossflow = _require_non_negative(crossflow_drag_coefficient, "crossflow_drag_coefficient")
    area = _require_positive(reference_area_m2, "reference_area_m2")
    ref_length = _require_positive(reference_length_m, "reference_length_m")
    x_m = float(moment_reference_m)
    if not isfinite(x_m):
        raise InputError("moment_reference_m must be finite.", field="moment_reference_m")
    geometry = body_geometry(body)
    s_b = geometry.base_area_m2
    a_p = geometry.planform_area_m2
    x_p = geometry.planform_centroid_m
    volume = geometry.volume_m3
    length = geometry.length_m
    alpha = radians(alpha_deg)
    if small_angle:
        potential_cn = 2.0 * (s_b / area) * alpha
        viscous_cn = crossflow * (a_p / area) * alpha * alpha
        lift = potential_cn + viscous_cn
        foredrag = (s_b / area) * alpha * alpha + crossflow * (a_p / area) * alpha**3
        moment = 2.0 * ((volume - s_b * (length - x_m)) / (area * ref_length)) * alpha + crossflow * (
            a_p / area
        ) * ((x_m - x_p) / ref_length) * alpha * alpha
        normal = lift
    else:
        potential_cn = (s_b / area) * sin(2.0 * alpha) * cos(0.5 * alpha)
        viscous_cn = crossflow * (a_p / area) * sin(alpha) ** 2
        normal = potential_cn + viscous_cn
        lift = potential_cn + viscous_cn * cos(alpha)
        foredrag = (s_b / area) * sin(2.0 * alpha) * sin(0.5 * alpha) + crossflow * (
            a_p / area
        ) * sin(alpha) ** 3
        moment = ((volume - s_b * (length - x_m)) / (area * ref_length)) * sin(2.0 * alpha) * cos(
            0.5 * alpha
        ) + crossflow * (a_p / area) * ((x_m - x_p) / ref_length) * sin(alpha) ** 2
    warnings: list[WarningMessage] = []
    if alpha_deg > 20.0:
        warnings.append(
            WarningMessage(
                "crossflow_alpha_beyond_validation",
                "NACA TR 1048 compared the method with tests up to about 20 degrees.",
                field="angle_of_attack_deg",
            )
        )
    if geometry.fineness_ratio < 13.0:
        warnings.append(
            WarningMessage(
                "crossflow_low_fineness",
                "NACA TR 1048 found the method less accurate for its fineness-ratio 13.1 body "
                "than for 21.1; lower fineness ratios were not tested.",
            )
        )
    if normal == 0.0:
        center = float("nan")
    else:
        center = x_m - moment * ref_length / normal
    return ViscousCrossflowResult(
        angle_of_attack_deg=alpha_deg,
        normal_force_coefficient=normal,
        lift_coefficient=lift,
        foredrag_increment_coefficient=foredrag,
        pitching_moment_coefficient=moment,
        center_of_pressure_m=center,
        potential_normal_force_coefficient=potential_cn,
        viscous_normal_force_coefficient=viscous_cn,
        small_angle=bool(small_angle),
        warnings=tuple(warnings),
    )


def crossflow_line_load(
    body: AxisymmetricBody,
    *,
    x_m: float,
    angle_of_attack_deg: float,
    dynamic_pressure_pa: float,
    crossflow_drag_coefficient: float,
) -> float:
    """Return the Allen-Perkins cross force per unit length at a station, N/m.

    Parameters
    ----------
    body : AxisymmetricBody
        Pointed body of revolution.
    x_m : float
        Station from the nose tip, m.
    angle_of_attack_deg : float
        Angle of attack in degrees, ``0 <= alpha < 90``.
    dynamic_pressure_pa : float
        Free-stream dynamic pressure ``q``, Pa.
    crossflow_drag_coefficient : float
        Local cross-flow drag coefficient ``c_d``.

    Returns
    -------
    float
        ``f = q (dS/dx) sin(2 alpha) cos(alpha/2) + 2 R c_d q sin^2(alpha)``,
        normal to the body axis [1]_ (printed p. 8). With ``c_d = 0`` and
        small ``alpha`` it reduces to Barrowman's running load
        ``n(x) = rho V^2 alpha dA/dx`` [2]_ (equation 3-59).

    References
    ----------
    .. [1] NACA-TR-1048. H. J. Allen and E. W. Perkins, NACA Report 1048, 1951,
           printed p. 8.
    .. [2] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equation 3-59.
    """

    alpha_deg = float(angle_of_attack_deg)
    if not isfinite(alpha_deg) or not 0.0 <= alpha_deg < 90.0:
        raise DomainError("Angle of attack must lie in [0, 90) degrees.", field="angle_of_attack_deg")
    q = _require_non_negative(dynamic_pressure_pa, "dynamic_pressure_pa")
    crossflow = _require_non_negative(crossflow_drag_coefficient, "crossflow_drag_coefficient")
    alpha = radians(alpha_deg)
    radius = body.radius_at(x_m)
    return q * body.area_slope_at(x_m) * sin(2.0 * alpha) * cos(0.5 * alpha) + (
        2.0 * radius * crossflow * q * sin(alpha) ** 2
    )


@dataclass(frozen=True, slots=True)
class TrapezoidalFinSet:
    """Three or four identical trapezoidal fins of Barrowman's figure 3-2 [1]_.

    Attributes
    ----------
    count : int
        Number of fins, 3 or 4 (equation 3-11 is derived for these only).
    root_chord_m, tip_chord_m : float
        ``c_r`` and ``c_t``, m.
    semispan_m : float
        Exposed semispan ``s`` measured from the root chord, m.
    tip_le_offset_m : float
        ``x_t``: axial distance from the root leading edge to the tip leading
        edge, m.
    root_le_station_m : float
        ``l_T``: distance from the nose tip to the root leading edge, m.
    body_radius_m : float
        Body radius at the tail ``r_t``, m.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           figure 3-2 and list of symbols.
    """

    count: int
    root_chord_m: float
    tip_chord_m: float
    semispan_m: float
    tip_le_offset_m: float
    root_le_station_m: float
    body_radius_m: float

    def __post_init__(self) -> None:
        if self.count not in (3, 4) or isinstance(self.count, bool):
            raise InputError("Barrowman's fin relations cover 3 or 4 fins.", field="count")
        _require_positive(self.root_chord_m, "root_chord_m")
        _require_non_negative(self.tip_chord_m, "tip_chord_m")
        _require_positive(self.semispan_m, "semispan_m")
        if not isfinite(self.tip_le_offset_m):
            raise InputError("tip_le_offset_m must be finite.", field="tip_le_offset_m")
        _require_non_negative(self.root_le_station_m, "root_le_station_m")
        _require_positive(self.body_radius_m, "body_radius_m")

    @property
    def panel_area_m2(self) -> float:
        """Planform area of one exposed panel ``A_T = (c_r + c_t) s / 2``, m^2.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
               derivation of equation 3-24.
        """

        return 0.5 * (self.root_chord_m + self.tip_chord_m) * self.semispan_m

    @property
    def midchord_sweep_rad(self) -> float:
        """Mid-chord sweep angle ``Gamma_c`` from the figure 3-2 geometry, rad.

        Notes
        -----
        The mid-chord line joins the root and tip mid-chord points, so
        ``tan(Gamma_c) = (x_t + c_t/2 - c_r/2) / s``.

        References
        ----------
        .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
               figure 3-2.
        """

        return atan(
            (self.tip_le_offset_m + 0.5 * self.tip_chord_m - 0.5 * self.root_chord_m)
            / self.semispan_m
        )


def barrowman_fin_normal_force_slope(
    fins: TrapezoidalFinSet,
    *,
    mach: float,
    aspect_ratio: float,
    reference_area_m2: float,
) -> float:
    """Return the subsonic ``(C_N_alpha)_T`` of the fin set (no body), 1/rad.

    Parameters
    ----------
    fins : TrapezoidalFinSet
        The fin set.
    mach : float
        Free-stream Mach number, ``0 <= M < 1``.
    aspect_ratio : float
        Barrowman's "aspect ratio of exposed fin" (list of symbols). The
        thesis does not state the formula used to compute it, so the caller
        must supply the value consistent with Diederich's correlation.
    reference_area_m2 : float
        Reference area ``A_r``, m^2.

    Returns
    -------
    float
        ``N pi AR (A_f/A_r) / (2 + sqrt(4 + (beta AR / cos(Gamma_c))^2))``
        (equation 3-12, from Diederich's correlation with
        ``C_N_alpha0 = 2 pi / beta``, equations 3-1 to 3-6 and 3-11), with
        ``beta = sqrt(1 - M^2)``.

    Raises
    ------
    DomainError
        For ``M >= 1``: the supersonic value needs the strip theory of
        Barrowman's appendix A, which is not implemented.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-1 to 3-12.
    """

    if not isinstance(fins, TrapezoidalFinSet):
        raise InputError("fins must be a TrapezoidalFinSet.", field="fins")
    mach_number = _require_non_negative(mach, "mach")
    if mach_number >= 1.0:
        raise DomainError(
            "Barrowman's equation 3-12 is subsonic; the supersonic fin method is not implemented.",
            field="mach",
        )
    ar = _require_positive(aspect_ratio, "aspect_ratio")
    reference = _require_positive(reference_area_m2, "reference_area_m2")
    beta = sqrt(1.0 - mach_number * mach_number)
    ratio = beta * ar / cos(fins.midchord_sweep_rad)
    return (
        fins.count * pi * ar * (fins.panel_area_m2 / reference) / (2.0 + sqrt(4.0 + ratio * ratio))
    )


def barrowman_fin_center_of_pressure(fins: TrapezoidalFinSet) -> tuple[float, float]:
    """Return the subsonic fin centre of pressure ``(X_T, Y_T)``, m.

    Parameters
    ----------
    fins : TrapezoidalFinSet
        The fin set.

    Returns
    -------
    tuple of float
        ``X_T`` from the nose tip (equation 3-30)::

            X_T = l_T + (x_t/3) (c_r + 2 c_t)/(c_r + c_t)
                  + (1/6) [c_r + c_t - c_r c_t/(c_r + c_t)]

        and the spanwise ``Y_T = r_t + (s/3)(c_r + 2 c_t)/(c_r + c_t)``
        (equation 3-27), the intersection of the quarter-chord line and the
        mean aerodynamic chord.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           section 3.12, equations 3-24 to 3-30.
    """

    if not isinstance(fins, TrapezoidalFinSet):
        raise InputError("fins must be a TrapezoidalFinSet.", field="fins")
    c_r = fins.root_chord_m
    c_t = fins.tip_chord_m
    total = c_r + c_t
    x_cp = (
        fins.root_le_station_m
        + fins.tip_le_offset_m / 3.0 * (c_r + 2.0 * c_t) / total
        + (total - c_r * c_t / total) / 6.0
    )
    y_cp = fins.body_radius_m + fins.semispan_m / 3.0 * (c_r + 2.0 * c_t) / total
    return x_cp, y_cp


@dataclass(frozen=True, slots=True)
class InterferenceFactors:
    """Slender-body wing(tail)-body interference factors [1]_.

    Attributes
    ----------
    radius_to_semispan : float
        ``r/s`` with ``s`` the semispan measured from the body axis.
    tail_in_presence_of_body : float
        ``K_W(B)`` (Barrowman's ``K_T(B)``), equation 14 of [1]_.
    body_in_presence_of_tail : float
        ``K_B(W)`` (Barrowman's ``K_B(T)``), equation 21 of [1]_.

    References
    ----------
    .. [1] NACA-TR-1307. W. C. Pitts, J. N. Nielsen and G. E. Kaattari, NACA
           Report 1307, 1959, equations 14 and 21.
    """

    radius_to_semispan: float
    tail_in_presence_of_body: float
    body_in_presence_of_tail: float


def interference_factors(*, body_radius_m: float, exposed_semispan_m: float) -> InterferenceFactors:
    """Return the slender-body interference factors ``K_W(B)`` and ``K_B(W)``.

    Parameters
    ----------
    body_radius_m : float
        Body radius ``r`` at the tail, m.
    exposed_semispan_m : float
        Exposed semispan measured from the root chord, m; the semispan from
        the axis is ``s = exposed + r`` (Barrowman's ``tau = (s + r_t)/r_t``).

    Returns
    -------
    InterferenceFactors
        With ``x = r/s``::

            K_W(B) = (2/pi) {(1 + x^4)[0.5 atan(0.5 (1/x - x)) + pi/4]
                     - x^2 [(1/x - x) + 2 atan(x)]} / (1 - x)^2
            K_B(W) = (1 - x^2)^2 / (1 - x)^2 - K_W(B)

    Notes
    -----
    NACA TR 1307 [1]_ (equation 14) gives limits ``K_W(B) = 1`` at
    ``r/s = 0`` and 2 as ``r/s -> 1``, and ``K_B(W) = 0`` and 2 (equation 21).
    Barrowman's transcription (equation 3-96 [2]_) omits the ``+ pi/4`` term;
    without it the factor is negative at ``tau = 2``, so the NACA TR 1307 form
    is used. Barrowman's equation 3-98 for ``K_B(T)`` agrees with equation 21.

    References
    ----------
    .. [1] NACA-TR-1307. W. C. Pitts, J. N. Nielsen and G. E. Kaattari, NACA
           Report 1307, 1959, printed pp. 4 and 6, equations 14 and 21.
    .. [2] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-96 and 3-98.
    """

    radius = _require_positive(body_radius_m, "body_radius_m")
    exposed = _require_positive(exposed_semispan_m, "exposed_semispan_m")
    x = radius / (exposed + radius)
    if x > 0.99999:
        raise DomainError(
            "r/s > 0.99999: equation 14 loses precision by cancellation (limit K = 2).",
            field="exposed_semispan_m",
        )
    inverse_minus = 1.0 / x - x
    k_wb = (
        (2.0 / pi)
        * (
            (1.0 + x**4) * (0.5 * atan(0.5 * inverse_minus) + 0.25 * pi)
            - x * x * (inverse_minus + 2.0 * atan(x))
        )
        / (1.0 - x) ** 2
    )
    k_bw = (1.0 - x * x) ** 2 / (1.0 - x) ** 2 - k_wb
    return InterferenceFactors(x, k_wb, k_bw)


def fin_in_presence_of_body(
    fins: TrapezoidalFinSet,
    *,
    mach: float,
    aspect_ratio: float,
    reference_area_m2: float,
) -> NormalForceComponent:
    """Return ``(C_N_alpha)_T(B) = (C_N_alpha)_T K_T(B)`` at the fin CP.

    Parameters
    ----------
    fins : TrapezoidalFinSet
        The fin set.
    mach : float
        Free-stream Mach number, ``0 <= M < 1``.
    aspect_ratio : float
        See :func:`barrowman_fin_normal_force_slope`.
    reference_area_m2 : float
        Reference area, m^2.

    Returns
    -------
    NormalForceComponent
        Slope (equation 3-93 with ``K_T(B)`` from NACA TR 1307 equation 14)
        at ``X_T(B) = X_T`` (equation 3-97).

    Notes
    -----
    The body-in-presence-of-tail term ``(C_N_alpha)_B(T)`` (equation 3-94)
    needs a centre of pressure that Barrowman gives only as a curve (figure
    3-10, subsonic); it is left to the caller, who can build it with
    :func:`interference_factors` and their own ``X_B(T)``.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-93 and 3-97.
    .. [2] NACA-TR-1307. W. C. Pitts, J. N. Nielsen and G. E. Kaattari, NACA
           Report 1307, 1959, equation 14.
    """

    slope = barrowman_fin_normal_force_slope(
        fins, mach=mach, aspect_ratio=aspect_ratio, reference_area_m2=reference_area_m2
    )
    factors = interference_factors(
        body_radius_m=fins.body_radius_m, exposed_semispan_m=fins.semispan_m
    )
    x_cp, _ = barrowman_fin_center_of_pressure(fins)
    return NormalForceComponent(
        "fins_in_presence_of_body", slope * factors.tail_in_presence_of_body, x_cp
    )
