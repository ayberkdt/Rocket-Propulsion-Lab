"""Taylor-Maccoll supersonic flow over a right circular cone at zero incidence.

The conical flow between an attached straight shock and the cone surface is
irrotational and isentropic; with all properties constant along rays from
the vertex the Euler equations reduce to the Taylor-Maccoll ordinary
differential equation in the polar angle [1]_. The equation is integrated
from the shock, where the oblique-shock relations give the initial state,
inward to the ray on which the normal velocity vanishes (the cone surface),
following the inverse procedure of Anderson section 13.6.3 [1]_. A direct
solve for a given cone angle iterates on the shock angle. Results are for a
calorically perfect gas and are checked against the Sims cone tables [2]_.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of Aerodynamics",
       6th ed., McGraw-Hill, 2017, section 13.6, eqs. (13.65)-(13.81),
       pp. 858-867.
.. [2] NASA-SP-3004. J. L. Sims, "Tables for Supersonic Flow Around Right
       Circular Cones at Zero Angle of Attack", NASA SP-3004, 1964.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, cos, degrees, isfinite, pi, radians, sin, sqrt, tan

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from ._advanced_common import bracketed_root, maximize_unimodal
from .isentropic import calculate_isentropic
from .normal_shock import calculate_normal_shock
from .oblique_shock import calculate_oblique_shock

_STEP_RAD = 2.0e-4
_MIN_STEP_RAD = 1.0e-12
_MARGIN_FRACTION = 0.01


@dataclass(frozen=True, slots=True)
class ConeFlowResult:
    """Taylor-Maccoll solution for a sharp cone at zero incidence [1]_.

    Attributes
    ----------
    free_stream_mach : float
        Upstream Mach number.
    gamma : float
        Ratio of specific heats.
    cone_half_angle_deg, shock_angle_deg : float
        Cone semivertex angle and conical shock-wave angle, degrees.
    post_shock_deflection_deg, post_shock_mach : float
        Flow deflection and Mach number immediately behind the shock.
    surface_mach : float
        Mach number on the cone surface.
    surface_pressure_ratio, surface_density_ratio, surface_temperature_ratio : float
        Surface ``p_c/p_inf``, ``rho_c/rho_inf`` and ``T_c/T_inf``.
    surface_pressure_coefficient : float
        ``(p_c - p_inf)/q_inf``.
    stagnation_pressure_ratio : float
        Total-pressure ratio across the shock (constant behind it).
    maximum_cone_half_angle_deg : float or None
        Detachment cone angle when it was computed, else ``None``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 13.6.
    """

    free_stream_mach: float
    gamma: float
    cone_half_angle_deg: float
    shock_angle_deg: float
    post_shock_deflection_deg: float
    post_shock_mach: float
    surface_mach: float
    surface_pressure_ratio: float
    surface_density_ratio: float
    surface_temperature_ratio: float
    surface_pressure_coefficient: float
    stagnation_pressure_ratio: float
    maximum_cone_half_angle_deg: float | None = None


def _validate(mach: float, gamma: float) -> None:
    if not isfinite(mach) or mach <= 1.0:
        raise DomainError("Free-stream Mach number must be greater than 1.", field="mach")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.", field="gamma")


def _wedge_deflection(mach: float, wave_angle: float, gamma: float) -> float:
    """theta-beta-M relation, Anderson eq. (9.23)."""
    mn2 = (mach * sin(wave_angle)) ** 2
    return atan(2.0 / tan(wave_angle) * (mn2 - 1.0) / (mach * mach * (gamma + cos(2.0 * wave_angle)) + 2.0))


def _nondimensional_speed(mach: float, gamma: float) -> float:
    """``V' = V/V_max = [2/((gamma-1) M^2) + 1]^(-1/2)``, Anderson eq. (13.81)."""
    return (2.0 / ((gamma - 1.0) * mach * mach) + 1.0) ** -0.5


def _mach_from_speed(speed: float, gamma: float) -> float:
    """Invert Anderson eq. (13.81)."""
    return sqrt(2.0 / ((gamma - 1.0) * (1.0 / (speed * speed) - 1.0)))


def _derivative(theta: float, vr: float, vt: float, gamma: float) -> tuple[float, float]:
    """Return ``(dVr/dtheta, dVtheta/dtheta)`` from Anderson eqs. (13.79)-(13.80)."""
    a2 = 0.5 * (gamma - 1.0) * (1.0 - vr * vr - vt * vt)
    denominator = a2 - vt * vt
    if denominator == 0.0:
        raise ConvergenceError("Taylor-Maccoll integration reached a singular point.")
    second = (vr * vt * vt - a2 * (2.0 * vr + vt / tan(theta))) / denominator
    return vt, second


def _rk4(theta: float, vr: float, vt: float, step: float, gamma: float) -> tuple[float, float]:
    k1r, k1t = _derivative(theta, vr, vt, gamma)
    k2r, k2t = _derivative(theta + 0.5 * step, vr + 0.5 * step * k1r, vt + 0.5 * step * k1t, gamma)
    k3r, k3t = _derivative(theta + 0.5 * step, vr + 0.5 * step * k2r, vt + 0.5 * step * k2t, gamma)
    k4r, k4t = _derivative(theta + step, vr + step * k3r, vt + step * k3t, gamma)
    return (
        vr + step * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0,
        vt + step * (k1t + 2.0 * k2t + 2.0 * k3t + k4t) / 6.0,
    )


def _integrate_from_shock(mach: float, shock_angle: float, gamma: float) -> tuple[float, float, float, float]:
    """Return ``(cone_angle, V'_r at surface, deflection, M2)`` for a shock angle.

    Notes
    -----
    Steps 1-4 of Anderson section 13.6.3: oblique-shock state behind the
    shock, its radial and normal velocity components, then fourth-order
    Runge-Kutta integration of eq. (13.80) toward the axis until
    ``V'_theta = 0``; the zero is located inside the last step by bisection
    of a single Runge-Kutta sub-step.
    """

    deflection = _wedge_deflection(mach, shock_angle, gamma)
    normal = calculate_normal_shock(mach * sin(shock_angle), gamma)
    m2 = normal.downstream_mach / sin(shock_angle - deflection)
    speed = _nondimensional_speed(m2, gamma)
    vr = speed * cos(shock_angle - deflection)
    vt = -speed * sin(shock_angle - deflection)
    theta = shock_angle
    while True:
        # Close to a weak shock the normal velocity is nearly sonic and the
        # coefficient a^2 - V_theta^2 of eq. (13.80) is small; the step is
        # limited to a fraction of 1 - (V_theta/a)^2 there.
        a2 = 0.5 * (gamma - 1.0) * (1.0 - vr * vr - vt * vt)
        margin = (a2 - vt * vt) / a2
        step = -min(_STEP_RAD, 0.5 * theta, max(_MIN_STEP_RAD, _MARGIN_FRACTION * margin))
        next_vr, next_vt = _rk4(theta, vr, vt, step, gamma)
        if not (isfinite(next_vr) and isfinite(next_vt)):
            raise ConvergenceError("Taylor-Maccoll integration produced a non-finite state.")
        if next_vt >= 0.0:
            def normal_velocity(sub: float, theta=theta, vr=vr, vt=vt) -> float:
                return _rk4(theta, vr, vt, -sub, gamma)[1]

            sub = bracketed_root(normal_velocity, 0.0, -step, x_tolerance=1.0e-15,
                                 label="cone surface ray")
            surface_vr = _rk4(theta, vr, vt, -sub, gamma)[0]
            return theta - sub, surface_vr, deflection, m2
        theta, vr, vt = theta + step, next_vr, next_vt
        if theta < 1.0e-9:
            raise ConvergenceError("Taylor-Maccoll integration reached the axis.")


def _result(
    mach: float, shock_angle: float, gamma: float, maximum_cone_deg: float | None
) -> ConeFlowResult:
    cone_angle, surface_vr, deflection, m2 = _integrate_from_shock(mach, shock_angle, gamma)
    normal = calculate_normal_shock(mach * sin(shock_angle), gamma)
    surface_mach = _mach_from_speed(surface_vr, gamma)
    free = calculate_isentropic(mach, gamma)
    surface = calculate_isentropic(surface_mach, gamma)
    total_ratio = normal.stagnation_pressure_ratio
    pressure_ratio = surface.pressure_total_ratio * total_ratio / free.pressure_total_ratio
    temperature_ratio = surface.temperature_total_ratio / free.temperature_total_ratio
    return ConeFlowResult(
        free_stream_mach=mach,
        gamma=gamma,
        cone_half_angle_deg=degrees(cone_angle),
        shock_angle_deg=degrees(shock_angle),
        post_shock_deflection_deg=degrees(deflection),
        post_shock_mach=m2,
        surface_mach=surface_mach,
        surface_pressure_ratio=pressure_ratio,
        surface_density_ratio=pressure_ratio / temperature_ratio,
        surface_temperature_ratio=temperature_ratio,
        surface_pressure_coefficient=2.0 / (gamma * mach * mach) * (pressure_ratio - 1.0),
        stagnation_pressure_ratio=total_ratio,
        maximum_cone_half_angle_deg=maximum_cone_deg,
    )


def cone_flow_from_shock_angle(
    mach: float, shock_angle_deg: float, gamma: float = 1.4
) -> ConeFlowResult:
    """Inverse Taylor-Maccoll problem: the cone supported by a given shock [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, greater than 1.
    shock_angle_deg : float
        Conical shock angle in degrees, between the Mach angle and 90.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    ConeFlowResult
        Cone angle and surface conditions.

    Raises
    ------
    DomainError
        If an input is out of range.

    Notes
    -----
    Implements the numerical procedure of Anderson section 13.6.3 [1]_:
    oblique-shock relations behind the shock (theta-beta-M eq. (9.23) and the
    normal-shock jump), ``V'`` from eq. (13.81), velocity components
    ``V'_r = V' cos(theta_s - delta)``, ``V'_theta = -V' sin(theta_s - delta)``,
    and fourth-order Runge-Kutta integration of the Taylor-Maccoll equation
    (13.80) with a ``2e-4`` rad step until ``V'_theta = 0``. Surface pressure
    follows from the isentropic relations with the post-shock total pressure.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 862-866.
    """

    _validate(mach, gamma)
    if not isfinite(shock_angle_deg):
        raise DomainError("Shock angle must be finite.", field="shock_angle_deg")
    shock_angle = radians(shock_angle_deg)
    if not asin(1.0 / mach) < shock_angle < 0.5 * pi:
        raise DomainError("Shock angle must lie between the Mach angle and 90 degrees.",
                          field="shock_angle_deg")
    return _result(mach, shock_angle, gamma, None)


def maximum_cone_half_angle(mach: float, gamma: float = 1.4) -> tuple[float, float]:
    """Largest cone angle with an attached conical shock [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, greater than 1.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    tuple of float
        ``(maximum cone half-angle, corresponding shock angle)`` in degrees.

    Notes
    -----
    The cone angle supported by a shock is maximised over the shock angle by
    golden-section search (Anderson Fig. 13.20 shows the maximum, beyond
    which the shock detaches) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 867.
    """

    _validate(mach, gamma)
    mach_angle = asin(1.0 / mach)
    shock, cone = maximize_unimodal(
        lambda beta: _integrate_from_shock(mach, beta, gamma)[0],
        mach_angle + 1.0e-6,
        0.5 * pi - 1.0e-4,
        iterations=45,
    )
    return degrees(cone), degrees(shock)


def solve_cone_flow(
    mach: float, cone_half_angle_deg: float, gamma: float = 1.4
) -> ConeFlowResult:
    """Direct Taylor-Maccoll problem: weak-shock flow over a given cone [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, greater than 1.
    cone_half_angle_deg : float
        Cone semivertex angle in degrees, positive.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    ConeFlowResult
        Weak-branch shock angle and surface conditions.

    Raises
    ------
    DomainError
        If an input is out of range or the cone is too blunt for an attached
        shock (the message gives the maximum cone angle).

    Notes
    -----
    The shock angle is bracketed between the Mach angle and the attached
    two-dimensional wedge wave angle for the same half-angle (the cone shock
    is weaker than the wedge shock, Anderson section 9.3 [1]_); if the wedge
    is detached the bracket is extended to the maximum-cone shock angle. The
    inverse problem of :func:`cone_flow_from_shock_angle` is then solved for
    the shock angle by bisection.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 634-636, 858-867.
    """

    _validate(mach, gamma)
    if not isfinite(cone_half_angle_deg) or cone_half_angle_deg <= 0.0:
        raise DomainError("Cone half-angle must be finite and positive.",
                          field="cone_half_angle_deg")
    target = radians(cone_half_angle_deg)
    mach_angle = asin(1.0 / mach)
    def cone_error(beta: float) -> float:
        return _integrate_from_shock(mach, beta, gamma)[0] - target

    offset = 1.0e-7
    lower = mach_angle + offset
    while cone_error(lower) > 0.0 and offset > 1.0e-14:
        offset *= 1.0e-2
        lower = mach_angle + offset

    upper = None
    maximum_deg = None
    try:
        wedge = calculate_oblique_shock(
            upstream_mach=mach, turning_angle_deg=cone_half_angle_deg, gamma=gamma
        )
        candidate = radians(wedge.wave_angle_deg)
        if cone_error(candidate) > 0.0:
            upper = candidate
    except DomainError:
        upper = None
    if upper is None:
        maximum_deg, peak_shock_deg = maximum_cone_half_angle(mach, gamma)
        if cone_half_angle_deg > maximum_deg:
            raise DomainError(
                f"Detached shock: maximum cone half-angle is {maximum_deg:.6g} degrees.",
                code="detached_shock",
                details={"maximum_cone_half_angle_deg": maximum_deg},
            )
        upper = radians(peak_shock_deg)
    shock = bracketed_root(cone_error, lower, upper, x_tolerance=1.0e-12,
                           label="cone shock angle")
    return _result(mach, shock, gamma, maximum_deg)
