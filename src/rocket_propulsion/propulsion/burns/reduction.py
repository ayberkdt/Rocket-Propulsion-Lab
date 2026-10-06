"""Auditable equivalent-constant reductions for propulsion burn artifacts.

References
----------
NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
NASA rocket thrust equation:
https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .models import BurnResult, ProfiledBurnResult


@dataclass(frozen=True, slots=True)
class ConstantBurnReduction:
    """Delivered/axial impulse and propellant preserving constant reduction.

    ``equivalent_thrust_n`` retains the historical delivered-thrust meaning.
    Dynamics consumers must use ``equivalent_axial_thrust_n`` together with
    ``equivalent_rocket_specific_impulse_s`` so hardware cant reduces force
    without incorrectly reducing total tank drain.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    NASA rocket thrust equation:
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    """

    source_schema: str
    source_result_hash: str
    duration_s: float
    equivalent_thrust_n: float
    equivalent_system_specific_impulse_s: float
    equivalent_axial_thrust_n: float
    equivalent_rocket_specific_impulse_s: float
    source_total_impulse_n_s: float
    source_axial_total_impulse_n_s: float
    source_consumed_propellant_kg: float
    thrust_relative_range: float
    mass_flow_relative_range: float
    thrust_centroid_offset_s: float
    impulse_residual_n_s: float
    axial_impulse_residual_n_s: float
    propellant_residual_kg: float
    exact: bool
    warnings: tuple[str, ...]


def _relative_range(values: tuple[float, ...]) -> float:
    maximum = max(values, default=0.0)
    minimum = min(values, default=0.0)
    return 0.0 if maximum <= 0.0 else (maximum - minimum) / maximum


def reduce_to_constant(
    result: BurnResult | ProfiledBurnResult,
    *,
    exact_tolerance: float = 1e-12,
) -> ConstantBurnReduction:
    """Preserve duration, delivered/axial impulse, and propellant.

    The reduction does not claim to preserve the force-time centroid, orbital
    trajectory, gravity loss, steering response, or tank-by-tank transient.
    ``exact`` is true only when both thrust and total mass flow are already
    constant within ``exact_tolerance``.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    NASA rocket thrust equation:
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    """

    if not isfinite(exact_tolerance) or exact_tolerance < 0.0:
        raise DomainError("Reduction exactness tolerance must be finite and non-negative.")
    summary = result.summary
    if summary.duration_s <= 0.0 or summary.consumed_propellant_kg <= 0.0:
        raise DomainError("A positive-duration, propellant-consuming burn is required.")
    equivalent_thrust = summary.delivered_total_impulse_n_s / summary.duration_s
    equivalent_isp = summary.delivered_total_impulse_n_s / (
        STANDARD_GRAVITY_M_S2 * summary.consumed_propellant_kg
    )
    equivalent_axial_thrust = (
        summary.axial_total_impulse_n_s / summary.duration_s
    )
    equivalent_rocket_isp = summary.axial_total_impulse_n_s / (
        STANDARD_GRAVITY_M_S2 * summary.consumed_propellant_kg
    )
    thrust_values = tuple(point.delivered_thrust_n for point in result.profile)
    flow_values = tuple(point.total_mass_flow_kg_s for point in result.profile)
    thrust_range = _relative_range(thrust_values)
    flow_range = _relative_range(flow_values)
    exact = thrust_range <= exact_tolerance and flow_range <= exact_tolerance
    impulse_residual = abs(
        equivalent_thrust * summary.duration_s
        - summary.delivered_total_impulse_n_s
    )
    axial_impulse_residual = abs(
        equivalent_axial_thrust * summary.duration_s
        - summary.axial_total_impulse_n_s
    )
    equivalent_flow = equivalent_thrust / (STANDARD_GRAVITY_M_S2 * equivalent_isp)
    propellant_residual = abs(
        equivalent_flow * summary.duration_s - summary.consumed_propellant_kg
    )
    warnings: tuple[str, ...] = ()
    if not exact:
        warnings = (
            (
                "Equivalent constant burn preserves duration, total impulse, and total "
                "propellant only."
            ),
            (
                "Force-time centroid and trajectory response are not preserved; evaluate "
                "the native profile in a transient-capable dynamics tool."
            ),
        )
    return ConstantBurnReduction(
        source_schema=result.schema,
        source_result_hash=result.result_hash,
        duration_s=summary.duration_s,
        equivalent_thrust_n=equivalent_thrust,
        equivalent_system_specific_impulse_s=equivalent_isp,
        equivalent_axial_thrust_n=equivalent_axial_thrust,
        equivalent_rocket_specific_impulse_s=equivalent_rocket_isp,
        source_total_impulse_n_s=summary.delivered_total_impulse_n_s,
        source_axial_total_impulse_n_s=summary.axial_total_impulse_n_s,
        source_consumed_propellant_kg=summary.consumed_propellant_kg,
        thrust_relative_range=thrust_range,
        mass_flow_relative_range=flow_range,
        thrust_centroid_offset_s=(
            summary.thrust_centroid_time_s - summary.duration_s / 2.0
        ),
        impulse_residual_n_s=impulse_residual,
        axial_impulse_residual_n_s=axial_impulse_residual,
        propellant_residual_kg=propellant_residual,
        exact=exact,
        warnings=warnings,
    )
