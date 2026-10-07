"""Free-space burn metrics for maneuver targeting and impulsive equivalents.

Maneuver targeting usually starts from an impulsive solution and then
replaces it with a finite burn.  The two need to be linked consistently:

* the finite burn's free-space velocity change,
  ``Dv = integral F(t) / m(t) dt``; and
* the time at which an equivalent impulse should be placed, here the
  acceleration-weighted centroid ``t* = integral t a dt / integral a dt``.

For a tabulated artifact, thrust and flow are piecewise linear, so ``m(t)`` is
piecewise quadratic and the integrand is smooth inside each segment.  Each
segment is integrated with Gauss-Legendre quadrature at ``n`` and ``2 n``
points, and the difference is reported as error evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import cos, log, pi

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real

from .tabulated import TabulatedBurnArtifact


@lru_cache(maxsize=16)
def gauss_legendre(points: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return Gauss-Legendre nodes and weights on ``[-1, 1]`` (Newton iteration).

    References
    ----------
    NIST Digital Library of Mathematical Functions, section 3.5(v):
    https://dlmf.nist.gov/3.5#v
    """

    points = integer(points, "Quadrature points")
    if not 1 <= points <= 200:
        raise DomainError("Quadrature points must be between 1 and 200.")
    def legendre(x: float) -> tuple[float, float]:
        previous, current = 1.0, x
        for k in range(2, points + 1):
            previous, current = current, ((2 * k - 1) * x * current - (k - 1) * previous) / k
        return current, points * (x * current - previous) / (x * x - 1.0)

    nodes: list[float] = []
    weights: list[float] = []
    for i in range(1, points + 1):
        x = cos(pi * (i - 0.25) / (points + 0.5))
        for _ in range(100):
            value, derivative = legendre(x)
            step = value / derivative
            x -= step
            if abs(step) < 1e-16:
                break
        _, derivative = legendre(x)
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return tuple(nodes), tuple(weights)


@dataclass(frozen=True, slots=True)
class FreeSpaceBurnMetrics:
    """Free-space (no gravity, no drag) burn summary.

    References
    ----------
    NIST Digital Library of Mathematical Functions, section 3.5(v):
    https://dlmf.nist.gov/3.5#v
    """

    initial_mass_kg: float
    final_mass_kg: float
    consumed_propellant_kg: float
    delta_v_m_s: float
    impulsive_equivalent_time_s: float
    effective_exhaust_velocity_m_s: float | None
    quadrature_points: int
    delta_v_error_estimate_m_s: float


def free_space_burn_metrics(
    artifact: TabulatedBurnArtifact,
    initial_mass_kg: float,
    *,
    thrust_scale: float = 1.0,
    mass_flow_scale: float = 1.0,
    quadrature_points: int = 16,
) -> FreeSpaceBurnMetrics:
    """Return Dv, final mass and the impulsive-equivalent time of a burn.

    Uses the artifact's *axial* thrust and total flow, scaled like the
    propagator bridge drivers.  ``impulsive_equivalent_time_s`` is measured
    from profile time zero.  ``effective_exhaust_velocity_m_s`` is
    ``Dv / ln(m0 / mf)``, the single exhaust velocity that reproduces this Dv.

    References
    ----------
    NIST Digital Library of Mathematical Functions, section 3.5(v):
    https://dlmf.nist.gov/3.5#v
    NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
    """

    if not isinstance(artifact, TabulatedBurnArtifact):
        raise DomainError("Burn metrics require a TabulatedBurnArtifact.")
    initial_mass_kg = real(initial_mass_kg, "Initial mass")
    thrust_scale = real(thrust_scale, "Thrust scale")
    mass_flow_scale = real(mass_flow_scale, "Mass-flow scale")
    if initial_mass_kg <= 0.0 or thrust_scale <= 0.0 or mass_flow_scale <= 0.0:
        raise DomainError("Initial mass and scales must be greater than zero.")
    consumed = artifact.consumed_propellant_kg * mass_flow_scale
    if consumed >= initial_mass_kg:
        raise DomainError("The burn consumes at least the whole initial mass.")

    def integrate(points: int) -> tuple[float, float]:
        nodes, weights = gauss_legendre(points)
        mass_at_start = initial_mass_kg
        delta_v = 0.0
        moment = 0.0
        for segment in artifact.segments:
            duration = segment.duration_s
            f0 = segment.start_axial_thrust_n * thrust_scale
            f1 = segment.end_axial_thrust_n * thrust_scale
            q0 = segment.start_total_mass_flow_kg_s * mass_flow_scale
            q1 = segment.end_total_mass_flow_kg_s * mass_flow_scale
            half = 0.5 * duration
            for node, weight in zip(nodes, weights, strict=True):
                tau = half * (node + 1.0)
                fraction = tau / duration
                thrust = f0 + (f1 - f0) * fraction
                mass = mass_at_start - (q0 * tau + 0.5 * (q1 - q0) * tau * fraction)
                acceleration = thrust / mass
                delta_v += half * weight * acceleration
                moment += half * weight * acceleration * (segment.start_time_s + tau)
            mass_at_start -= 0.5 * (q0 + q1) * duration
        return delta_v, moment

    points = integer(quadrature_points, "Quadrature points")
    if not 1 <= points <= 100:
        raise DomainError("Quadrature points must be between 1 and 100 (refined at 2n).")
    delta_v, _ = integrate(points)
    refined, refined_moment = integrate(2 * points)
    final_mass = initial_mass_kg - consumed
    ratio = log(initial_mass_kg / final_mass)
    return FreeSpaceBurnMetrics(
        initial_mass_kg=initial_mass_kg,
        final_mass_kg=final_mass,
        consumed_propellant_kg=consumed,
        delta_v_m_s=refined,
        impulsive_equivalent_time_s=(refined_moment / refined if refined > 0.0 else 0.0),
        effective_exhaust_velocity_m_s=(refined / ratio if ratio > 0.0 else None),
        quadrature_points=2 * points,
        delta_v_error_estimate_m_s=abs(refined - delta_v),
    )
