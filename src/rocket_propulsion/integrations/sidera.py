"""Optional, fail-closed adapter from propulsion burn results to Sidera burns.

The adapter inspects only the installed public Sidera interface and never
searches local source paths. Hand-off documents can be produced without
Sidera; native objects are built only when a compatible installation is
found.

References
----------
.. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
.. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
       502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
"""

from __future__ import annotations

import inspect
from dataclasses import asdict
from importlib import import_module
from typing import Any, Literal

from rocket_propulsion import __version__
from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError
from rocket_propulsion.core.validation import real, text
from rocket_propulsion.propulsion.burns import (
    BurnResult,
    ProfiledBurnResult,
    burn_result_to_dict,
    reduce_to_constant,
)
from rocket_propulsion.propulsion.burns.serialization import burn_result_from_dict

from .contracts import SideraCapabilities
from .direction import _ADAPTER_DIRECTION_POLICY, _directed_vector

_FINITE_BURN_FIELDS = {
    "t_start_s",
    "duration_s",
    "thrust_n",
    "isp_s",
    "direction",
    "frame",
    "throttle",
}
_FRAMES = ("inertial", "ric", "vnb")


def inspect_sidera_capabilities() -> SideraCapabilities:
    """Inspect the installed public Sidera interface and report its capabilities.

    Returns
    -------
    SideraCapabilities
        ``status="unavailable"`` when Sidera cannot be imported,
        ``"incompatible"`` when its public propagation facade cannot be read
        or ``FiniteBurn``/``ManeuverPlan`` do not expose the required
        contract, and ``"supported"`` otherwise.

    Notes
    -----
    Only ``sidera`` and ``sidera.core.propagation`` are imported; local
    paths are never searched. ``FiniteBurn`` is compatible when its
    constructor accepts every field the adapter writes. Provider
    sub-capabilities are inferred from constructor parameter names; a
    provider whose signature cannot be inspected is reported as present
    without sub-capabilities.
    """

    try:
        sidera = import_module("sidera")
        propagation = import_module("sidera.core.propagation")
    except (ImportError, ModuleNotFoundError) as error:
        return SideraCapabilities(
            status="unavailable",
            version=None,
            finite_burn_available=False,
            maneuver_plan_available=False,
            constructor_compatible=False,
            supported_frames=(),
            constant_burn_exact=False,
            tabulated_burn_available=False,
            independent_mass_flow_available=False,
            reason=str(error),
        )

    try:
        finite_burn = getattr(propagation, "FiniteBurn", None)
        maneuver_plan = getattr(propagation, "ManeuverPlan", None)
        tabulated = getattr(propagation, "TabulatedFiniteBurn", None)
        weighted_ensemble = getattr(propagation, "WeightedManeuverEnsemble", None)
        propulsion_provider = getattr(propagation, "PropulsionProvider", None)
    except Exception as error:  # noqa: BLE001 - optional integration boundary
        return SideraCapabilities(
            status="incompatible",
            version=str(getattr(sidera, "__version__", "unknown")),
            finite_burn_available=False,
            maneuver_plan_available=False,
            constructor_compatible=False,
            supported_frames=(),
            constant_burn_exact=False,
            tabulated_burn_available=False,
            independent_mass_flow_available=False,
            reason=f"Sidera public propagation facade could not load: {error}",
        )
    compatible = False
    reason: str | None = None
    if finite_burn is not None:
        try:
            parameters = set(inspect.signature(finite_burn).parameters)
            compatible = _FINITE_BURN_FIELDS <= parameters
            if not compatible:
                reason = "FiniteBurn constructor does not expose the required public fields."
        except (TypeError, ValueError) as error:
            reason = f"FiniteBurn constructor could not be inspected: {error}"
    else:
        reason = "Sidera does not expose FiniteBurn from sidera.core.propagation."

    exact = compatible and maneuver_plan is not None
    independent_mass_flow = (
        tabulated is not None and "mass_flow_kg_s" in inspect.signature(tabulated).parameters
    )
    weighted_ensemble_available = weighted_ensemble is not None
    provider_parameters: set[str] = set()
    if propulsion_provider is not None:
        try:
            provider_parameters = set(inspect.signature(propulsion_provider).parameters)
        except (TypeError, ValueError):
            # Presence is still reported, but unsupported/introspection-hostile
            # providers do not receive inferred sub-capabilities.
            provider_parameters = set()
    return SideraCapabilities(
        status="supported" if exact else "incompatible",
        version=str(getattr(sidera, "__version__", "unknown")),
        finite_burn_available=finite_burn is not None,
        maneuver_plan_available=maneuver_plan is not None,
        constructor_compatible=compatible,
        supported_frames=_FRAMES if compatible else (),
        constant_burn_exact=exact,
        tabulated_burn_available=tabulated is not None,
        independent_mass_flow_available=independent_mass_flow,
        reason=reason,
        weighted_maneuver_ensemble_available=weighted_ensemble_available,
        exact_weighted_tabulated_ensemble=(
            tabulated is not None and independent_mass_flow and weighted_ensemble_available
        ),
        propulsion_provider_available=propulsion_provider is not None,
        additional_states_available="additional_states" in provider_parameters,
        event_detectors_available="event_detectors" in provider_parameters,
        analytic_partials_available="partials" in provider_parameters,
        body_frame_force_moment_available="body_wrench" in provider_parameters,
    )


def _normalized_direction(direction: tuple[float, float, float]) -> tuple[float, float, float]:
    """Return the unit vector of ``direction`` rounded once per component."""
    return _directed_vector(direction, 1., label="Sidera unit direction")


def _validated_export_inputs(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: str,
    require_exact: bool,
) -> tuple[tuple[float, float, float], str]:
    """Validate export inputs and return ``(unit direction, normalized frame)``.

    The burn result is re-opened through its public integrity boundary
    (serialization and hash verification) before use. With
    ``require_exact`` the result must be eligible for exact constant export.
    """
    if not isinstance(result, (BurnResult, ProfiledBurnResult)):
        raise DomainError("Sidera export requires a BurnResult or ProfiledBurnResult.")
    # Re-open the artifact through the public integrity boundary before using it.
    burn_result_from_dict(burn_result_to_dict(result))
    if require_exact and not result.summary.exact_constant_sidera_export_available:
        raise DomainError("Burn result is not eligible for exact constant Sidera export.")
    real(t_start_s, "Sidera start time")
    normalized_frame = text(frame, "Sidera frame").strip().lower()
    if normalized_frame not in _FRAMES:
        raise DomainError("Sidera finite-burn frame must be inertial, ric, or vnb.")
    return _normalized_direction(direction), normalized_frame


def export_sidera_artifact(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
    mode: Literal["exact", "equivalent_constant"] = "exact",
) -> dict[str, Any]:
    """Create a Sidera finite-burn hand-off document without importing Sidera.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Verified burn artifact; it is re-opened through its integrity
        boundary first.
    t_start_s : float
        Burn start time in the consumer's time scale, s.
    direction : tuple of float
        Thrust direction in ``frame``; normalized on export.
    frame : {"inertial", "ric", "vnb"}, optional
        Sidera finite-burn frame.
    mode : {"exact", "equivalent_constant"}, optional
        ``"exact"`` requires a burn that is exactly representable as a
        constant burn; ``"equivalent_constant"`` reduces any burn to a
        constant one and records the reduction.

    Returns
    -------
    dict
        ``{"schema": "sidera_maneuver_plan_v1", "maneuvers": [...],
        "manifest": {...}}`` with one ``finite_burn`` maneuver and a
        manifest of source hashes, semantics and omitted domains.

    Raises
    ------
    DomainError
        If ``mode`` or ``frame`` is unsupported, ``result`` is not a burn
        result or fails re-verification, the direction or start time is
        invalid, or exact export is requested for an ineligible burn.

    Notes
    -----
    Thrust is the net axial thrust after hardware cant. For exact export
    ``thrust = axial total impulse / duration`` and ``isp_s`` is the
    rocket-equivalent specific impulse, so ``isp_s`` with standard gravity
    represents the total tank drain [1]_. Profiled results are always
    reduced. The equivalent-constant reduction preserves duration, total
    impulse and total propellant, not the force-time history or the
    trajectory response. ``throttle`` is written as 1 because the operating
    point already contains the delivered thrust. The manifest follows the
    versioned, unit-explicit exchange practice of [2]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
           502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    if not isinstance(mode, str) or mode not in ("exact", "equivalent_constant"):
        raise DomainError("Sidera export mode must be exact or equivalent_constant.")

    normalized_direction, normalized_frame = _validated_export_inputs(
        result,
        t_start_s=t_start_s,
        direction=direction,
        frame=frame,
        require_exact=mode == "exact",
    )
    reduction = None
    if mode == "equivalent_constant" or isinstance(result, ProfiledBurnResult):
        reduction = reduce_to_constant(result)
        thrust_n = reduction.equivalent_axial_thrust_n
        isp_s = reduction.equivalent_rocket_specific_impulse_s
    else:
        thrust_n = result.summary.axial_total_impulse_n_s / result.summary.duration_s
        isp_s = result.summary.rocket_equivalent_specific_impulse_s
    burn = {
        "kind": "finite_burn",
        "t_start_s": real(t_start_s, "Sidera start time"),
        "duration_s": result.summary.duration_s,
        "thrust_n": thrust_n,
        "isp_s": isp_s,
        "direction": list(normalized_direction),
        "frame": normalized_frame,
        # The operating point already contains actual delivered thrust. Applying
        # its original throttle here would incorrectly scale force and flow twice.
        "throttle": 1.0,
    }
    manifest: dict[str, object] = {
        "source_schema": result.schema,
        "source_input_hash": result.input_hash,
        "source_result_hash": result.result_hash,
        "rocket_propulsion_version": __version__,
        "mapping": "exact_constant" if mode == "exact" else "equivalent_constant",
        "direction_numerics": _ADAPTER_DIRECTION_POLICY,
        "mass_flow_semantics": "total_tank_drain_via_system_isp",
        "force_semantics": "net_axial_thrust_after_hardware_cant",
        "source_delivered_impulse_n_s": result.summary.delivered_total_impulse_n_s,
        "source_axial_impulse_n_s": result.summary.axial_total_impulse_n_s,
        "source_total_tank_flow_kg_s": result.operating_point.total_tank_flow_kg_s,
        "source_commanded_throttle": result.operating_point.commanded_throttle,
        "source_realized_throttle": result.operating_point.realized_throttle,
        "target_throttle": 1.0,
        "omitted_domains": [
            "tank_stream_breakdown",
            "tank_reserves",
            "propulsion_warnings",
        ],
    }
    if reduction is not None:
        manifest["reduction"] = asdict(reduction)
    if mode == "equivalent_constant" and reduction is not None and not reduction.exact:
        manifest["approximation_warning"] = (
            "Equivalent-constant export preserves duration, total impulse, and total "
            "propellant; it does not preserve the force-time history or trajectory response."
        )
    return {
        "schema": "sidera_maneuver_plan_v1",
        "maneuvers": [burn],
        "manifest": manifest,
    }


def to_native_finite_burn(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
) -> Any:
    """Construct Sidera's public ``FiniteBurn`` from an exactly constant burn.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Verified burn artifact eligible for exact constant export.
    t_start_s : float
        Burn start time, s.
    direction : tuple of float
        Thrust direction in ``frame``.
    frame : {"inertial", "ric", "vnb"}, optional
        Sidera finite-burn frame.

    Returns
    -------
    object
        ``sidera.core.propagation.FiniteBurn``.

    Raises
    ------
    FeatureUnavailableError
        If no compatible Sidera installation is found; ``details`` carries
        the observed capabilities.
    DomainError
        If the burn or the export inputs are invalid (see
        :func:`export_sidera_artifact`).

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    capabilities = inspect_sidera_capabilities()
    if not capabilities.constant_burn_exact:
        raise FeatureUnavailableError(
            "A compatible Sidera FiniteBurn installation is unavailable.",
            details={"capabilities": asdict(capabilities)},
        )
    artifact = export_sidera_artifact(result, t_start_s=t_start_s, direction=direction, frame=frame)
    finite_burn = import_module("sidera.core.propagation").FiniteBurn
    return finite_burn(
        **{key: value for key, value in artifact["maneuvers"][0].items() if key != "kind"}
    )


def to_native_maneuver_plan(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
) -> Any:
    """Wrap the exact native finite burn in Sidera's public ``ManeuverPlan``.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Verified burn artifact eligible for exact constant export.
    t_start_s : float
        Burn start time, s.
    direction : tuple of float
        Thrust direction in ``frame``.
    frame : {"inertial", "ric", "vnb"}, optional
        Sidera finite-burn frame.

    Returns
    -------
    object
        ``sidera.core.propagation.ManeuverPlan`` holding one ``FiniteBurn``.

    Raises
    ------
    FeatureUnavailableError
        If no compatible Sidera installation is found.
    DomainError
        If the burn or the export inputs are invalid.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    burn = to_native_finite_burn(result, t_start_s=t_start_s, direction=direction, frame=frame)
    maneuver_plan = import_module("sidera.core.propagation").ManeuverPlan
    return maneuver_plan((burn,))
