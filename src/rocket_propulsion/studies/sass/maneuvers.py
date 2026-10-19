"""Analytical mission inputs, with impulsive and low-thrust models kept distinct.

References
----------
MIT-522-L5: orbital energy and thrust acceleration.
MIT-522-L6: equations (6)-(8), circular transfers and their applicability.
"""

from __future__ import annotations

from math import isfinite, pi, sqrt

from rocket_propulsion.core.errors import DomainError

from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity


def hohmann_transfer(initial_radius: CitedScalar, final_radius: CitedScalar,
                     gravitational_parameter: CitedScalar) -> AnalysisResult:
    """Return signed tangential impulses between coplanar circular orbits.

    Radii are measured from the central-body centre, not altitudes. No body
    constants, atmosphere, combined plane change or optimality claim is assumed.
    Equal radii require no transfer; reported coast time is then zero.

    References
    ----------
    MIT-522-L6: equation (7), Hohmann impulses.
    MIT-522-L5: equation (11), specific orbital energy.
    MIT-512-L35: page 5, equation (2), half the elliptical orbital period.
    """
    r0 = quantity(initial_radius, "m", strictly_positive=True)
    r1 = quantity(final_radius, "m", strictly_positive=True)
    mu = quantity(gravitational_parameter, "m^3/s^2", strictly_positive=True)
    a = min(r0, r1) + 0.5 * (max(r0, r1) - min(r0, r1))
    v0, v1 = sqrt(mu) / sqrt(r0), sqrt(mu) / sqrt(r1)
    if not all(isfinite(v) and v > 0 for v in (v0, v1)):
        raise DomainError("Orbital speed is outside the positive finite numerical domain.")
    # Rationalize sqrt(1+x)-1; direct speed subtraction loses adjacent radii.
    # Here x=(r1-r0)/(r1+r0), evaluated without overflowing the radius sum.
    x = ((r1 - r0) / a) * 0.5
    dv0 = v0 * x / (sqrt(r1 / a) + 1)
    dv1 = v1 * x / (sqrt(r0 / a) + 1)
    coast = pi * (a / sqrt(mu)) * sqrt(a) if r0 != r1 else 0.0
    inputs = (initial_radius, final_radius, gravitational_parameter)
    methods = (Citation("MIT-522-L6", "Eq. (7): circular two-impulse transfer"),
               Citation("MIT-522-L5", "Keplerian orbital energy"),
               Citation("MIT-512-L35", "Page 5, Eq. (2): half the transfer ellipse period"))
    values = (
        ("departure_delta_v", dv0, "m/s"), ("arrival_delta_v", dv1, "m/s"),
        ("total_delta_v", abs(dv0) + abs(dv1), "m/s"),
        ("transfer_coast_time", coast, "s"),
        ("initial_orbital_period", 2 * pi * (r0 / v0), "s"),
    )
    return AnalysisResult("coplanar_circular_hohmann", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in values
    ), ("Two-body, instantaneous tangential impulses; not an electric-thruster trajectory.",
        "No atmosphere, perturbations, encounter phasing or disposal-lifetime prediction."))


def circular_spiral_delta_v(initial_radius: CitedScalar, final_radius: CitedScalar,
                           gravitational_parameter: CitedScalar) -> AnalysisResult:
    """Estimate near-circular continuous tangential-thrust velocity expenditure.

    This is a separate low-thrust budget: abs(sqrt(mu/r0)-sqrt(mu/r1)).
    It is not a finite-burn correction factor applied to a Hohmann result.

    References
    ----------
    MIT-522-L6: equations (4)-(8), small tangential acceleration and near-circular orbits.
    """
    r0 = quantity(initial_radius, "m", strictly_positive=True)
    r1 = quantity(final_radius, "m", strictly_positive=True)
    mu = quantity(gravitational_parameter, "m^3/s^2", strictly_positive=True)
    methods = (Citation("MIT-522-L6", "Eq. (6): near-circular tangential spiral"),)
    small, large = min(r0, r1), max(r0, r1)
    speed = sqrt(mu) / sqrt(small)
    # Same rationalization as the cited speed difference, avoiding cancellation.
    difference = speed * ((large - small) / large) / (1 + sqrt(small / large))
    value = derived(difference, "m/s",
                    (initial_radius, final_radius, gravitational_parameter), methods)
    return AnalysisResult("near_circular_low_thrust_spiral", (("total_delta_v", value),), (
        "Approximate, coplanar, near-circular and small tangential acceleration only.",
        "No eclipse, steering, power history, escape or transfer-time solution.",
    ))


def perigee_lowering_delta_v(circular_radius: CitedScalar, target_perigee_radius: CitedScalar,
                           gravitational_parameter: CitedScalar) -> AnalysisResult:
    """Size one retrograde impulse; do not confuse it with a two-burn circularization.

    References
    ----------
    MIT-522-L5: equation (11), vis-viva, with a=(ra+rp)/2.
    MIT-522-L6: equation (7), first tangential impulse, descending transfer.
    """
    ra = quantity(circular_radius, "m", strictly_positive=True)
    rp = quantity(target_perigee_radius, "m", strictly_positive=True)
    mu = quantity(gravitational_parameter, "m^3/s^2", strictly_positive=True)
    if rp > ra:
        raise DomainError("Target perigee must not exceed the starting circular radius.")
    a = rp + .5 * (ra - rp)
    # Only the departure impulse is needed. A period overflow must not
    # invalidate an otherwise representable one-impulse calculation.
    impulse = (sqrt(mu) / sqrt(ra)) * (((ra - rp) / a) * .5) / (sqrt(rp / a) + 1)
    output = derived(impulse, "m/s", (circular_radius, target_perigee_radius,
                                     gravitational_parameter), (
        Citation("MIT-522-L5", "Eq. (11): a=(ra+rp)/2, apogee velocity"),
    ))
    return AnalysisResult("single_impulse_perigee_lowering", (("total_delta_v", output),), (
        "Only an osculating two-body target perigee; no re-entry time or compliance claim.",
    ))


def burn_duration_screen(on_time: CitedScalar, orbital_period: CitedScalar,
                         maximum_period_fraction: CitedScalar) -> AnalysisResult:
    """Compare burn/period to an explicitly chosen impulsive-screen threshold.

    No universal threshold is invented: the analyst supplies and cites it.

    References
    ----------
    MIT-522-L6: page 2, impulsive transfer versus continuous thrust assumptions.
    NASA-SE-6: section 6.8, declare decision criteria and assumptions.
    """
    t = quantity(on_time, "s", minimum=0)
    period = quantity(orbital_period, "s", strictly_positive=True)
    limit = quantity(maximum_period_fraction, "1", strictly_positive=True)
    if limit >= 1:
        raise DomainError("Impulsive-screen fraction must be less than one.")
    inputs = (on_time, orbital_period, maximum_period_fraction)
    methods = (Citation("MIT-522-L6", "High-thrust impulsive versus low-thrust assumption"),
               Citation("NASA-SE-6", "6.8: analyst-defined acceptance threshold"))
    return AnalysisResult("analyst_defined_impulsive_screen", (
        ("burn_to_period_ratio", derived(t / period, "1", inputs, methods)),
        ("fraction_margin", derived(limit - t / period, "1", inputs, methods)),
    ), ("A positive margin passes this chosen screen; it does not verify trajectory accuracy.",))
