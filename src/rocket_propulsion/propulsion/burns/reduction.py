"""Auditable equivalent-constant reductions of propulsion burn results.

A variable-thrust burn is replaced by a constant burn with the same duration,
total impulse and consumed propellant. Equivalent thrust is impulse divided by
duration [1]_, and equivalent specific impulse is impulse divided by ``g0``
times the consumed propellant mass [2]_. Residuals of the reconstruction and
measures of profile nonuniformity are reported with the result.

References
----------
.. [1] NASA Glenn, "Rocket Thrust Equation".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
.. [2] NASA Glenn, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .models import BurnResult, ProfiledBurnResult


@dataclass(frozen=True, slots=True)
class ConstantBurnReduction:
    """Constant burn preserving duration, delivered/axial impulse and propellant.

    ``equivalent_thrust_n`` is the delivered (uncanted) thrust. Dynamics
    consumers must use ``equivalent_axial_thrust_n`` together with
    ``equivalent_rocket_specific_impulse_s``, so that hardware cant reduces the
    axial force without reducing the total tank drain.

    Parameters
    ----------
    source_schema, source_result_hash : str
        Schema and hash of the reduced result.
    duration_s : float
        Burn duration ``t_b`` [s].
    equivalent_thrust_n : float
        ``I / t_b`` with delivered total impulse ``I`` [N] [1]_.
    equivalent_system_specific_impulse_s : float
        ``I / (g0 m_p)`` with consumed propellant ``m_p`` [s] [2]_.
    equivalent_axial_thrust_n : float
        ``I_ax / t_b`` with axial total impulse ``I_ax`` [N].
    equivalent_rocket_specific_impulse_s : float
        ``I_ax / (g0 m_p)`` [s].
    source_total_impulse_n_s, source_axial_total_impulse_n_s : float
        ``I`` and ``I_ax`` of the source [N s].
    source_consumed_propellant_kg : float
        ``m_p`` [kg].
    thrust_relative_range, mass_flow_relative_range : float
        ``(max - min) / max`` over the profile samples; zero when the maximum
        is not positive.
    thrust_centroid_offset_s : float
        Thrust centroid time minus ``t_b / 2`` [s]; the constant burn has its
        centroid at mid-burn.
    impulse_residual_n_s, axial_impulse_residual_n_s : float
        Floating-point reconstruction residuals of the impulses [N s].
    propellant_residual_kg : float
        Reconstruction residual of the propellant mass [kg].
    exact : bool
        ``True`` when thrust and total mass flow are constant within the
        reduction tolerance.
    warnings : tuple of str
        Limits of the reduction when it is not exact.

    References
    ----------
    .. [1] NASA Glenn, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    .. [2] NASA Glenn, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
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
    """Return ``(max - min) / max``, or zero when the maximum is not positive."""
    maximum = max(values, default=0.0)
    minimum = min(values, default=0.0)
    return 0.0 if maximum <= 0.0 else (maximum - minimum) / maximum


def reduce_to_constant(
    result: BurnResult | ProfiledBurnResult,
    *,
    exact_tolerance: float = 1e-12,
) -> ConstantBurnReduction:
    """Reduce a burn result to a constant burn with the same impulse and propellant.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Burn with positive duration and positive consumed propellant.
    exact_tolerance : float, optional
        Finite, nonnegative tolerance on the relative thrust and mass-flow
        ranges for declaring the reduction exact.

    Returns
    -------
    ConstantBurnReduction
        Equivalent thrusts and specific impulses [1]_ [2]_, residuals,
        nonuniformity measures and warnings.

    Raises
    ------
    DomainError
        If the tolerance is invalid or the burn has no positive duration or
        propellant consumption.

    Notes
    -----
    The reduction preserves duration, delivered and axial impulse and consumed
    propellant only. It does not preserve the force-time centroid, orbital
    trajectory, gravity loss, steering response or tank-by-tank transient;
    a non-exact reduction should be checked against the native profile in a
    transient-capable dynamics tool. ``exact`` is ``True`` only when both
    thrust and total mass flow are already constant within
    ``exact_tolerance``.

    References
    ----------
    .. [1] NASA Glenn, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    .. [2] NASA Glenn, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
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
