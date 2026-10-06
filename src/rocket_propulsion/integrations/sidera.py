"""Optional, fail-closed adapter from propulsion results to Sidera burns."""

from __future__ import annotations

import inspect
from dataclasses import asdict
from importlib import import_module
from math import isfinite, sqrt
from typing import Any, Literal

from rocket_propulsion import __version__
from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError
from rocket_propulsion.propulsion.burns import (
    BurnResult,
    ProfiledBurnResult,
    burn_result_to_dict,
    reduce_to_constant,
)
from rocket_propulsion.propulsion.burns.serialization import burn_result_from_dict

from .contracts import SideraCapabilities

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
    """Inspect installed public Sidera contracts without searching local paths."""

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
        tabulated is not None
        and "mass_flow_kg_s" in inspect.signature(tabulated).parameters
    )
    weighted_ensemble_available = weighted_ensemble is not None
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
            tabulated is not None
            and independent_mass_flow
            and weighted_ensemble_available
        ),
    )


def _normalized_direction(direction: tuple[float, float, float]) -> tuple[float, float, float]:
    if len(direction) != 3 or any(not isfinite(value) for value in direction):
        raise DomainError("Sidera direction must contain exactly three finite components.")
    norm = sqrt(sum(value * value for value in direction))
    if norm <= 0.0:
        raise DomainError("Sidera direction must be non-zero.")
    return tuple(value / norm for value in direction)


def _validated_export_inputs(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: str,
    require_exact: bool,
) -> tuple[tuple[float, float, float], str]:
    # Re-open the artifact through the public integrity boundary before using it.
    burn_result_from_dict(burn_result_to_dict(result))
    if require_exact and not result.summary.exact_constant_sidera_export_available:
        raise DomainError("Burn result is not eligible for exact constant Sidera export.")
    if not isfinite(t_start_s):
        raise DomainError("Sidera start time must be finite.")
    normalized_frame = str(frame).strip().lower()
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
    """Create an exact or explicitly approximate hand-off without importing Sidera."""

    if mode not in {"exact", "equivalent_constant"}:
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
        "t_start_s": float(t_start_s),
        "duration_s": result.summary.duration_s,
        "thrust_n": thrust_n,
        "isp_s": isp_s,
        "direction": list(normalized_direction),
        "frame": normalized_frame,
        # The operating point already contains actual delivered thrust. Applying
        # its original throttle here would incorrectly scale force and flow twice.
        "throttle": 1.0,
    }
    manifest = {
        "source_schema": result.schema,
        "source_input_hash": result.input_hash,
        "source_result_hash": result.result_hash,
        "rocket_propulsion_version": __version__,
        "mapping": "exact_constant" if mode == "exact" else "equivalent_constant",
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
    """Construct Sidera's public ``FiniteBurn`` or fail with a typed error."""

    capabilities = inspect_sidera_capabilities()
    if not capabilities.constant_burn_exact:
        raise FeatureUnavailableError(
            "A compatible Sidera FiniteBurn installation is unavailable.",
            details={"capabilities": asdict(capabilities)},
        )
    artifact = export_sidera_artifact(
        result, t_start_s=t_start_s, direction=direction, frame=frame
    )
    finite_burn = import_module("sidera.core.propagation").FiniteBurn
    return finite_burn(
        **{
            key: value
            for key, value in artifact["maneuvers"][0].items()
            if key != "kind"
        }
    )


def to_native_maneuver_plan(
    result: BurnResult | ProfiledBurnResult,
    *,
    t_start_s: float,
    direction: tuple[float, float, float],
    frame: Literal["inertial", "ric", "vnb"] = "inertial",
) -> Any:
    """Wrap the exact native finite burn in Sidera's public ``ManeuverPlan``."""

    burn = to_native_finite_burn(
        result, t_start_s=t_start_s, direction=direction, frame=frame
    )
    maneuver_plan = import_module("sidera.core.propagation").ManeuverPlan
    return maneuver_plan((burn,))

