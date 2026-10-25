"""Optional, fail-closed native Sidera thrust provider for the propulsion bridge.

The wrapper flies a :class:`PropagatorPropulsionBridge` burn through Sidera's
public numerical propagation entry point (``sidera.core.propagation.propagate``)
without modifying Sidera and without reducing the burn to an equivalent
constant Isp. It uses only three public Sidera contracts:

1. **Dynamics protocol.** Each smooth burn segment is propagated with a
   :class:`SideraBridgeSegmentDynamics`, an object satisfying Sidera's
   ``PropagationDynamics`` protocol that delegates gravity, frames, epochs and
   restart points to the caller's Sidera engine and adds the bridge's thrust
   acceleration ``F / m`` [2]_ to ``dy[3:6]`` and its total mass rate to
   ``dy[6]``, following the force/mass-derivative split of a propagator
   propulsion model [1]_.
2. **Frame convention.** The commanded direction is given on Sidera's
   maneuver axes (``"inertial"``, ``"ric"`` or ``"vnb"``) and resolved from
   the current state at every evaluation with Sidera's public
   ``maneuver_frame_basis``, the frame definitions of CCSDS orbit messages
   [3]_.
3. **Events.** Time surfaces of the bridge (burn start, interior profile
   boundaries, burn stop) are honoured as segment boundaries, exactly as
   Sidera segments its own profiled burns, so no accepted integration step
   contains a jump or kink of the thrust. The protected dry-mass surface is
   registered as a terminal Sidera event ``g(t, y) = y[6] - m_dry`` located
   by SciPy's sign-change root finding [4]_; the burn then stops and the run
   coasts on.

Sidera's propagated state is the 6-state or the 7-state with mass, so named
tank states, tank reserves, operating-condition surfaces and ambient pressure
have no native channel and are refused with :class:`FeatureUnavailableError`
instead of being approximated [5]_.

References
----------
.. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
       EventDetector" API: force, mass-derivative, parameter-driver and
       event contract.
.. [2] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical Astrodynamics",
       equation (19), thrust acceleration F/m,
       https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf
.. [3] CCSDS-502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
       502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
.. [4] SCIPY-1.14.1-IVP. SciPy 1.14.1, "solve_ivp": event sign-change
       detection, terminal and direction attributes.
.. [5] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
"""

from __future__ import annotations

import inspect
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from importlib import import_module
from types import SimpleNamespace
from typing import Any

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError
from rocket_propulsion.core.validation import real, text
from rocket_propulsion.propulsion.burns import BurnEventKind, PropagatorPropulsionBridge
from rocket_propulsion.propulsion.burns.propagator import INVENTORY_LIMIT_TOLERANCE
from rocket_propulsion.propulsion.burns.time_coordinates import (
    _clock_table,
    _profile_boundaries,
)

from .contracts import SideraCapabilities
from .direction import _directed_vector
from .sidera import inspect_sidera_capabilities

SIDERA_NATIVE_FRAMES = ("inertial", "ric", "vnb")
SIDERA_NATIVE_REQUIRED_CAPABILITIES = (
    "propagate_available",
    "dynamics_protocol_available",
    "extra_events_available",
    "frame_basis_available",
)
# Members of Sidera's public ``PropagationDynamics`` protocol that the wrapper
# forwards; their presence on the protocol class is part of the probe.
_PROTOCOL_MEMBERS = (
    "flags",
    "grav",
    "build_rhs",
    "body_context",
    "body_fixed_rotation",
    "ephemeris_window_s",
    "solver_epoch_tdb_s",
    "restart_plan",
    "kernel_provenance",
)
_PROPAGATE_PARAMETERS = {"dynamics", "y0", "cfg", "time_cfg"}
# Sidera's own segment-end acceptance for maneuver segments (absolute, s).
_SEGMENT_END_TOLERANCE_S = 1.0e-9


@dataclass(frozen=True, slots=True)
class SideraNativeCapabilities:
    """Native-wrapper capabilities observed on the installed public Sidera interface.

    Every flag records what inspection of the installed package found; a
    False flag means the capability was not observed, not that it is known
    to be absent upstream.

    Attributes
    ----------
    status : str
        ``"supported"``, ``"incompatible"`` or ``"unavailable"``.
    version : str or None
        Reported Sidera version, ``"unknown"`` if it exposes none, ``None``
        when Sidera is not importable.
    propagate_available : bool
        ``sidera.core.propagation.propagate`` accepts ``dynamics``, ``y0``,
        ``cfg`` and ``time_cfg``.
    dynamics_protocol_available : bool
        ``sidera.core.propagation.dynamics_protocol`` exposes
        ``PropagationDynamics`` with every forwarded member and
        ``RestartPlan``.
    extra_events_available : bool
        ``propagate`` accepts ``extra_events``.
    frame_basis_available : bool
        ``sidera.core.propagation.plans.maneuver_frame_basis`` is callable.
    additional_states_available : bool
        Copied from :class:`SideraCapabilities`; named tank states can be
        carried beside the vehicle mass.
    operating_conditions_available, ambient_pressure_available : bool
        A host channel for performance-surface conditions or ambient
        pressure was observed. No such public channel is probed for, so both
        are always False.
    reason : str or None
        Explanation when a capability is missing or inspection failed.
    provider : SideraCapabilities or None
        The general Sidera capability record the probe started from.
    """

    status: str
    version: str | None
    propagate_available: bool
    dynamics_protocol_available: bool
    extra_events_available: bool
    frame_basis_available: bool
    additional_states_available: bool = False
    operating_conditions_available: bool = False
    ambient_pressure_available: bool = False
    reason: str | None = None
    provider: SideraCapabilities | None = None


@dataclass(frozen=True, slots=True)
class SideraNativeEventMapping:
    """How one bridge event surface is represented in a native Sidera run.

    Attributes
    ----------
    event_id : str
        Bridge event identifier, for example ``"burn:stop"``.
    kind : str
        :class:`BurnEventKind` value.
    mechanism : {"segment_boundary", "terminal_event"}
        ``"segment_boundary"``: a time surface on which one Sidera
        propagation segment ends and the next starts. ``"terminal_event"``:
        a Sidera extra event that stops the burn segment at its root.
    date_s : float or None
        Absolute Sidera time of a time surface, s.
    threshold : float or None
        Limit value of an inventory surface, kg.
    """

    event_id: str
    kind: str
    mechanism: str
    date_s: float | None = None
    threshold: float | None = None


@dataclass(frozen=True, slots=True)
class SideraNativeEventRecord:
    """A bridge event that a native Sidera run honoured.

    Attributes
    ----------
    event_id : str
        Bridge event identifier.
    kind : str
        :class:`BurnEventKind` value.
    time_s : float
        Absolute Sidera time at which the event was applied, s.
    vehicle_mass_kg : float
        Vehicle mass ``y[6]`` at that time, kg.
    """

    event_id: str
    kind: str
    time_s: float
    vehicle_mass_kg: float


@dataclass(frozen=True, slots=True)
class SideraNativeBurnRun:
    """Segments, honoured events and mass ledger of one native Sidera burn run.

    Attributes
    ----------
    segments : tuple of (str, object)
        ``(label, sidera PropagationResult)`` in time order; labels are
        ``"coast:pre-burn"``, ``"burn:segment:<k>"`` and ``"coast:post-burn"``.
    events : tuple of SideraNativeEventRecord
        Bridge events honoured, in time order.
    ignition_date_s, stop_date_s : float
        Planned absolute burn start and stop, s (ignition bias included).
    cutoff_date_s : float or None
        Absolute time at which thrust actually ended, s; ``None`` when the
        run stopped before ignition.
    mass_at_ignition_kg, mass_at_cutoff_kg : float or None
        Vehicle mass at ignition and at cutoff, kg.
    completed : bool
        The burn reached its stop surface.
    stopped_early : bool
        A caller-supplied terminal event (or an impact) ended the run.
    stop_reason : str or None
        Sidera's stop reason when ``stopped_early``.
    bridge_model_hash : str
        Hash of the propulsion bridge that was flown.
    """

    segments: tuple[tuple[str, Any], ...]
    events: tuple[SideraNativeEventRecord, ...]
    ignition_date_s: float
    stop_date_s: float
    cutoff_date_s: float | None
    mass_at_ignition_kg: float | None
    mass_at_cutoff_kg: float | None
    completed: bool
    stopped_early: bool
    stop_reason: str | None
    bridge_model_hash: str

    @property
    def consumed_propellant_kg(self) -> float | None:
        """float or None: Propagated mass at ignition minus mass at cutoff, kg."""

        if self.mass_at_ignition_kg is None or self.mass_at_cutoff_kg is None:
            return None
        return self.mass_at_ignition_kg - self.mass_at_cutoff_kg

    def stitched(self) -> tuple[Any, Any]:
        """Return the run's samples as one time array and one state array.

        Returns
        -------
        tuple of numpy.ndarray
            ``(t, y)`` with shapes ``(n,)`` and ``(n, 7)``. A segment's first
            row is dropped when it repeats the previous segment's last time
            (absolute tolerance 1e-12 s), as Sidera does for its own
            maneuver segments.
        """

        numpy = import_module("numpy")
        times: list[float] = []
        rows: list[Any] = []
        for _, result in self.segments:
            seg_t = numpy.asarray(result.t, dtype=numpy.float64)
            seg_y = numpy.asarray(result.y, dtype=numpy.float64)
            start = 1 if times and seg_t.size and abs(float(seg_t[0]) - times[-1]) <= 1e-12 else 0
            times.extend(float(value) for value in seg_t[start:])
            rows.extend(seg_y[start:])
        return numpy.asarray(times, dtype=numpy.float64), numpy.vstack(rows)


def _unavailable(reason: str | None, provider: SideraCapabilities) -> SideraNativeCapabilities:
    """Return the capability record of an absent Sidera installation."""
    return SideraNativeCapabilities(
        status="unavailable",
        version=None,
        propagate_available=False,
        dynamics_protocol_available=False,
        extra_events_available=False,
        frame_basis_available=False,
        reason=reason,
        provider=provider,
    )


def inspect_sidera_native_capabilities() -> SideraNativeCapabilities:
    """Inspect the installed public Sidera interface for the native thrust wrapper.

    Returns
    -------
    SideraNativeCapabilities
        ``status="unavailable"`` when Sidera cannot be imported,
        ``"incompatible"`` when a required public entry point is missing or
        cannot be inspected, and ``"supported"`` otherwise.

    Notes
    -----
    Only installed modules are imported (``sidera``,
    ``sidera.core.propagation``, its ``dynamics_protocol`` and ``plans``
    modules); local source paths are never searched. Presence is checked by
    attribute and signature inspection; nothing is propagated. The probe
    follows the fail-closed model-applicability practice of [1]_.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    provider = inspect_sidera_capabilities()
    if provider.status == "unavailable":
        return _unavailable(provider.reason, provider)
    try:
        propagation = import_module("sidera.core.propagation")
        protocol_module = import_module("sidera.core.propagation.dynamics_protocol")
        plans = import_module("sidera.core.propagation.plans")
        propagate = getattr(propagation, "propagate", None)
        protocol = getattr(protocol_module, "PropagationDynamics", None)
        restart_plan = getattr(protocol_module, "RestartPlan", None)
        frame_basis = getattr(plans, "maneuver_frame_basis", None)
    except Exception as error:  # noqa: BLE001 - optional integration boundary
        return SideraNativeCapabilities(
            status="incompatible",
            version=provider.version,
            propagate_available=False,
            dynamics_protocol_available=False,
            extra_events_available=False,
            frame_basis_available=False,
            reason=f"Sidera public propagation modules could not load: {error}",
            provider=provider,
        )

    reasons: list[str] = []
    parameters: set[str] = set()
    if callable(propagate):
        try:
            parameters = set(inspect.signature(propagate).parameters)
        except (TypeError, ValueError) as error:
            reasons.append(f"propagate could not be inspected: {error}")
    else:
        reasons.append("sidera.core.propagation does not expose propagate.")
    propagate_available = _PROPAGATE_PARAMETERS <= parameters
    if callable(propagate) and parameters and not propagate_available:
        reasons.append("propagate does not accept dynamics, y0, cfg and time_cfg.")
    extra_events_available = "extra_events" in parameters
    if parameters and not extra_events_available:
        reasons.append("propagate does not accept extra_events.")
    protocol_available = (
        protocol is not None
        and restart_plan is not None
        and all(hasattr(protocol, name) for name in _PROTOCOL_MEMBERS)
    )
    if not protocol_available:
        reasons.append("PropagationDynamics/RestartPlan protocol is not exposed as expected.")
    frame_basis_available = callable(frame_basis)
    if not frame_basis_available:
        reasons.append("plans.maneuver_frame_basis is not exposed.")
    supported = (
        propagate_available
        and protocol_available
        and extra_events_available
        and frame_basis_available
    )
    if not provider.additional_states_available:
        reasons.append(
            "Named additional (tank) states are not observed; propagate carries the "
            "6-state or the 7-state with mass only."
        )
    return SideraNativeCapabilities(
        status="supported" if supported else "incompatible",
        version=provider.version,
        propagate_available=propagate_available,
        dynamics_protocol_available=protocol_available,
        extra_events_available=extra_events_available,
        frame_basis_available=frame_basis_available,
        additional_states_available=provider.additional_states_available,
        reason=" ".join(reasons) or None,
        provider=provider,
    )


def missing_sidera_native_capabilities(
    capabilities: SideraNativeCapabilities,
    bridge: PropagatorPropulsionBridge,
) -> tuple[str, ...]:
    """Return the native-wrapper capabilities that ``bridge`` needs but Sidera lacks.

    Parameters
    ----------
    capabilities : SideraNativeCapabilities
        Observed capabilities.
    bridge : PropagatorPropulsionBridge
        Burn to be flown. Stream bindings or reserves additionally require
        ``additional_states_available``; a performance surface requires
        ``operating_conditions_available``; a nozzle exit area requires
        ``ambient_pressure_available``.

    Returns
    -------
    tuple of str
        Names of missing capability flags; empty when none is missing.

    Raises
    ------
    DomainError
        If either argument has the wrong type.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    """

    if not isinstance(capabilities, SideraNativeCapabilities):
        raise DomainError("Expected SideraNativeCapabilities.")
    if not isinstance(bridge, PropagatorPropulsionBridge):
        raise DomainError("Expected a PropagatorPropulsionBridge.")
    required = list(SIDERA_NATIVE_REQUIRED_CAPABILITIES)
    if bridge.stream_bindings or bridge.state_constraints:
        required.append("additional_states_available")
    if bridge.performance_surface is not None or bridge.stream_performance_surface is not None:
        required.append("operating_conditions_available")
    if bridge.nozzle_exit_area_m2 is not None:
        required.append("ambient_pressure_available")
    return tuple(name for name in required if not getattr(capabilities, name))


def require_sidera_native_capabilities(
    capabilities: SideraNativeCapabilities,
    bridge: PropagatorPropulsionBridge,
) -> None:
    """Fail closed unless Sidera can fly ``bridge`` natively without approximation.

    Parameters
    ----------
    capabilities : SideraNativeCapabilities
        Observed capabilities.
    bridge : PropagatorPropulsionBridge
        Burn to be flown.

    Raises
    ------
    FeatureUnavailableError
        If a required capability is missing; ``details`` lists the missing
        capabilities, the bridge hash, the Sidera status, version and reason.
    DomainError
        If either argument has the wrong type.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    missing = missing_sidera_native_capabilities(capabilities, bridge)
    if missing:
        raise FeatureUnavailableError(
            "Sidera cannot fly this propulsion bridge natively.",
            details={
                "missing_capabilities": list(missing),
                "bridge_model_hash": bridge.model_hash,
                "sidera_status": capabilities.status,
                "sidera_version": capabilities.version,
                "reason": capabilities.reason,
            },
        )


def _sidera_modules() -> SimpleNamespace:
    """Import the public Sidera members the wrapper calls, or fail closed."""
    try:
        propagation = import_module("sidera.core.propagation")
        protocol_module = import_module("sidera.core.propagation.dynamics_protocol")
        plans = import_module("sidera.core.propagation.plans")
        return SimpleNamespace(
            numpy=import_module("numpy"),
            propagate=propagation.propagate,
            protocol=protocol_module.PropagationDynamics,
            frame_basis=plans.maneuver_frame_basis,
        )
    except (ImportError, AttributeError) as error:
        raise FeatureUnavailableError(
            "The public Sidera propagation interface could not be imported.",
            details={"reason": str(error)},
        ) from error


class SideraBridgeSegmentDynamics:
    """Sidera ``PropagationDynamics`` for one smooth segment of a bridge burn.

    Parameters
    ----------
    burn : SideraNativeBurn
        Validated native burn.
    index : int
        Zero-based profile segment, between the canonical dates
        ``burn.segment_dates_s[index]`` and ``[index + 1]``.

    Raises
    ------
    DomainError
        If ``burn`` is not a :class:`SideraNativeBurn` or ``index`` is out of
        range.

    Notes
    -----
    Every protocol member other than ``build_rhs`` is forwarded to the
    wrapped Sidera engine, so the propagator sees the same gravity model,
    epoch, ephemeris window, impact sphere and restart points. ``drag_model``
    is forwarded too, so Sidera's integrator guards still see a runtime drag
    model. No private ``_build_solver_rhs`` is exposed: Sidera's
    ``propagate`` prefers that hook when present and would otherwise bypass
    the thrust, so the wrapper deliberately offers only the public,
    state-validating ``build_rhs``.

    The right-hand side requires the 7-state with mass at ``y[6]``. It
    resolves the direction from the state with Sidera's
    ``maneuver_frame_basis`` [3]_, evaluates the bridge at the absolute
    Sidera time, takes the right-hand limit at the segment start and the
    left-hand limit at the segment end, and adds ``F / m`` [2]_ and the total
    mass rate, as in [1]_. Inside a segment the thrust and flow are linear
    in time, so the right-hand side is smooth there.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [2] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical Astrodynamics",
           equation (19), thrust acceleration F/m,
           https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf
    .. [3] CCSDS-502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
           502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    def __init__(self, burn: SideraNativeBurn, index: int) -> None:
        if not isinstance(burn, SideraNativeBurn):
            raise DomainError("Segment dynamics require a SideraNativeBurn.")
        if isinstance(index, bool) or not isinstance(index, int):
            raise DomainError("Segment index must be an integer.")
        if not 0 <= index < len(burn.segment_dates_s) - 1:
            raise DomainError("Segment index is outside the burn profile.")
        self._burn = burn
        self._base = burn.base_dynamics
        self._index = index
        self._start_s = burn.segment_dates_s[index]
        self._end_s = burn.segment_dates_s[index + 1]

    @property
    def segment_window_s(self) -> tuple[float, float]:
        """tuple of float: Absolute Sidera start and end of this segment, s."""

        return self._start_s, self._end_s

    @property
    def flags(self) -> Any:
        """object: The wrapped engine's perturbation flags."""

        return self._base.flags

    @property
    def grav(self) -> Any:
        """object: The wrapped engine's central gravity model."""

        return self._base.grav

    @property
    def drag_model(self) -> Any:
        """object or None: The wrapped engine's runtime drag model, if any."""

        return getattr(self._base, "drag_model", None)

    def body_context(self) -> Any:
        """Return the wrapped engine's central-body context."""

        return self._base.body_context()

    def body_fixed_rotation(self) -> Any:
        """Return the wrapped engine's inertial-to-body-fixed rotation, or ``None``."""

        return self._base.body_fixed_rotation()

    def ephemeris_window_s(self) -> tuple[float, float] | None:
        """Return the wrapped engine's tabulated-ephemeris window, or ``None``."""

        return self._base.ephemeris_window_s()

    def solver_epoch_tdb_s(self) -> float | None:
        """Return the wrapped engine's TDB epoch at solver time 0, or ``None``."""

        return self._base.solver_epoch_tdb_s()

    def restart_plan(self) -> Any:
        """Return the wrapped engine's restart plan; the burn adds no kink inside a segment."""

        return self._base.restart_plan()

    def kernel_provenance(self) -> Mapping[str, Any] | None:
        """Return the wrapped engine's SPICE-kernel provenance, or ``None``."""

        return self._base.kernel_provenance()

    def thrust(self, t: float, y: Sequence[float]) -> tuple[tuple[float, float, float], float]:
        """Return the bridge's inertial thrust acceleration and mass rate at ``(t, y)``.

        Parameters
        ----------
        t : float
            Absolute Sidera time, s.
        y : sequence of float
            7-state ``[r (m), v (m/s), mass (kg)]`` on the engine's inertial
            axes.

        Returns
        -------
        tuple
            ``(acceleration (m/s**2), dm/dt (kg/s))``; the mass rate is
            negative while firing.

        Raises
        ------
        DomainError
            If the state is not a 7-state with positive mass or the bridge
            refuses the evaluation.

        References
        ----------
        .. [1] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical Astrodynamics",
               equation (19), thrust acceleration F/m,
               https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf
        """

        burn = self._burn
        numpy = burn._modules.numpy
        state = numpy.asarray(y, dtype=numpy.float64).reshape(-1)
        if state.size != 7:
            raise DomainError("A native Sidera burn needs the 7-state with mass at y[6].")
        mass = real(float(state[6]), "Sidera vehicle mass")
        if mass <= 0.0:
            raise DomainError("Sidera vehicle mass must be positive during a burn.")
        time_s = real(float(t), "Sidera propagation time")
        basis = burn._modules.frame_basis(state, burn.frame)
        resolved = basis @ numpy.asarray(burn.direction, dtype=numpy.float64)
        side = "left" if time_s == self._end_s else "right"
        value = burn.bridge.evaluate_at_epoch(
            time_s=time_s,
            ignition_time_s=burn.ignition_epoch_s,
            vehicle_mass_kg=mass,
            direction=(float(resolved[0]), float(resolved[1]), float(resolved[2])),
            parameter_overrides=burn.parameter_overrides,
            side=side,
            apply_inventory_limits=burn.bridge.protected_dry_mass_kg is None,
        )
        return value.acceleration_m_s2, value.mass_derivative_kg_s

    def build_rhs(self, *, force_rebuild: bool = False) -> Callable[[float, Any], Any]:
        """Return the wrapped engine's public right-hand side with the bridge thrust added.

        Parameters
        ----------
        force_rebuild : bool, optional
            Forwarded to the wrapped engine's ``build_rhs``.

        Returns
        -------
        callable
            ``f(t, y)`` returning a new array: the engine's derivative with
            the thrust acceleration added to ``[3:6]`` and the mass rate to
            ``[6]``.
        """

        numpy = self._burn._modules.numpy
        base_rhs = self._base.build_rhs(force_rebuild=force_rebuild)
        thrust = self.thrust

        def rhs(t: float, y: Any) -> Any:
            derivative = numpy.array(base_rhs(t, y), dtype=numpy.float64, copy=True)
            acceleration, mass_rate = thrust(t, y)
            derivative[3] += acceleration[0]
            derivative[4] += acceleration[1]
            derivative[5] += acceleration[2]
            derivative[6] += mass_rate
            return derivative

        return rhs


class SideraNativeBurn:
    """Fly a propagator-bridge burn natively through Sidera's public propagator.

    Parameters
    ----------
    base_dynamics : object
        Sidera engine satisfying ``PropagationDynamics`` (for example from
        ``sidera.core.body_engine.orbits.build_body_dynamics``).
    bridge : PropagatorPropulsionBridge
        Mass-only bridge: no stream bindings, reserves, performance surface
        or nozzle exit area (see :func:`missing_sidera_native_capabilities`).
        ``protected_dry_mass_kg`` is supported.
    ignition_epoch_s : float
        Absolute Sidera time at which the bridge's relative burn time is
        zero, s; the ``ignition_time_bias_s`` driver shifts it.
    direction : tuple of float
        Thrust direction on ``frame`` axes; normalized once.
    frame : {"inertial", "ric", "vnb"}, optional
        Sidera maneuver frame of ``direction``.
    parameter_overrides : mapping, optional
        Values of ``thrust_scale``, ``mass_flow_scale`` and
        ``ignition_time_bias_s``.
    capabilities : SideraNativeCapabilities, optional
        Observed capabilities; probed with
        :func:`inspect_sidera_native_capabilities` when omitted.

    Raises
    ------
    FeatureUnavailableError
        If Sidera is absent or lacks a capability this bridge needs.
    DomainError
        If the bridge, base dynamics, epoch, direction, frame or overrides
        are invalid, or the clock cannot distinguish two profile boundaries.

    Notes
    -----
    The capability gate runs before anything Sidera-specific is touched, so
    an absent installation is always reported as
    :class:`FeatureUnavailableError`. The burn is flown as Sidera flies its
    own profiled burns: one propagation segment per profile segment, bounded
    by the canonical dates ``ignition + bias + t_k`` rounded once from the
    exact rational sum [4]_. Start, boundary and stop surfaces are therefore
    segment boundaries, and the protected dry-mass surface is a terminal
    Sidera event on ``y[6]`` [3]_. The bridge is evaluated with
    ``apply_inventory_limits=False`` when that event is registered, so no
    Runge-Kutta stage sees the clamp. Thrust is the bridge's net axial force
    ``F / m`` [2]_ along the state-resolved direction and the mass rate is
    the bridge's total tank drain [1]_; no equivalent-constant Isp is formed.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [2] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical Astrodynamics",
           equation (19), thrust acceleration F/m,
           https://ocw.mit.edu/courses/16-522-space-propulsion-spring-2015/6127948308ccbe73c1c44751e2509ae0_MIT16_522S15_Lecture5.pdf
    .. [3] SCIPY-1.14.1-IVP. SciPy 1.14.1, "solve_ivp": event sign-change
           detection, terminal and direction attributes.
    .. [4] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """

    def __init__(
        self,
        base_dynamics: Any,
        bridge: PropagatorPropulsionBridge,
        *,
        ignition_epoch_s: float,
        direction: tuple[float, float, float],
        frame: str = "inertial",
        parameter_overrides: Mapping[str, float] | None = None,
        capabilities: SideraNativeCapabilities | None = None,
    ) -> None:
        if not isinstance(bridge, PropagatorPropulsionBridge):
            raise DomainError("A native Sidera burn requires a PropagatorPropulsionBridge.")
        observed = inspect_sidera_native_capabilities() if capabilities is None else capabilities
        if not isinstance(observed, SideraNativeCapabilities):
            raise DomainError("Expected SideraNativeCapabilities.")
        require_sidera_native_capabilities(observed, bridge)
        modules = _sidera_modules()
        if not isinstance(base_dynamics, modules.protocol):
            raise DomainError(
                "Base dynamics must satisfy Sidera's public PropagationDynamics protocol."
            )
        # Re-open the bridge through its public integrity boundary before use.
        PropagatorPropulsionBridge.from_dict(bridge.to_dict())
        normalized_frame = text(frame, "Sidera maneuver frame").strip().lower()
        if normalized_frame not in SIDERA_NATIVE_FRAMES:
            raise DomainError("Sidera maneuver frame must be inertial, ric or vnb.")
        drivers = bridge.resolved_parameters(parameter_overrides)
        ignition = real(ignition_epoch_s, "Sidera ignition epoch")
        dates, _ = _clock_table(
            _profile_boundaries(bridge.artifact.segments, None),
            ignition,
            drivers["ignition_time_bias_s"],
        )
        self._modules = modules
        self._base = base_dynamics
        self._bridge = bridge
        self._ignition = ignition
        self._direction = _directed_vector(direction, 1.0, label="Sidera thrust direction")
        self._frame = normalized_frame
        self._overrides = dict(drivers)
        self._dates = tuple(float(date) for date in dates)
        self._capabilities = observed

    @property
    def base_dynamics(self) -> Any:
        """object: The wrapped Sidera engine."""

        return self._base

    @property
    def bridge(self) -> PropagatorPropulsionBridge:
        """PropagatorPropulsionBridge: The burn being flown."""

        return self._bridge

    @property
    def ignition_epoch_s(self) -> float:
        """float: Absolute nominal ignition time before the bias, s."""

        return self._ignition

    @property
    def direction(self) -> tuple[float, float, float]:
        """tuple of float: Unit thrust direction on :attr:`frame` axes."""

        return self._direction

    @property
    def frame(self) -> str:
        """str: Sidera maneuver frame of :attr:`direction`."""

        return self._frame

    @property
    def parameter_overrides(self) -> dict[str, float]:
        """dict: Every resolved estimable-parameter value (copy)."""

        return dict(self._overrides)

    @property
    def capabilities(self) -> SideraNativeCapabilities:
        """SideraNativeCapabilities: Capabilities the burn was admitted with."""

        return self._capabilities

    @property
    def segment_dates_s(self) -> tuple[float, ...]:
        """tuple of float: Canonical absolute profile boundaries, start and stop included, s."""

        return self._dates

    def event_mapping(self) -> tuple[SideraNativeEventMapping, ...]:
        """Return how every bridge event surface is represented natively.

        Returns
        -------
        tuple of SideraNativeEventMapping
            Time surfaces as ``"segment_boundary"`` at their canonical
            dates; the dry-mass surface as ``"terminal_event"``.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        """

        mappings: list[SideraNativeEventMapping] = []
        boundary = 0
        for surface in self._bridge.event_surfaces():
            if surface.kind is BurnEventKind.START:
                date = self._dates[0]
            elif surface.kind is BurnEventKind.PROFILE_BOUNDARY:
                boundary += 1
                date = self._dates[boundary]
            elif surface.kind is BurnEventKind.STOP:
                date = self._dates[-1]
            else:
                mappings.append(
                    SideraNativeEventMapping(
                        surface.event_id,
                        surface.kind.value,
                        "terminal_event",
                        threshold=surface.threshold,
                    )
                )
                continue
            mappings.append(
                SideraNativeEventMapping(
                    surface.event_id, surface.kind.value, "segment_boundary", date_s=date
                )
            )
        return tuple(mappings)

    def segment_dynamics(self, index: int) -> SideraBridgeSegmentDynamics:
        """Return the Sidera dynamics of profile segment ``index``.

        Parameters
        ----------
        index : int
            Zero-based profile segment.

        Returns
        -------
        SideraBridgeSegmentDynamics
            Dynamics valid on ``[segment_dates_s[index], segment_dates_s[index + 1]]``.
        """

        return SideraBridgeSegmentDynamics(self, index)

    def inventory_events(self) -> tuple[Callable[[float, Any], float], ...]:
        """Return the bridge's inventory surfaces as Sidera extra events.

        Returns
        -------
        tuple of callable
            Empty without ``protected_dry_mass_kg``; otherwise one event
            ``g(t, y) = y[6] - m_dry`` with ``terminal=True``,
            ``direction=-1`` and ``uses_full_state=True`` (Sidera passes
            only ``y[:6]`` to extra events without that flag).

        Notes
        -----
        The event refuses a state without ``y[6]``, so an installation that
        ignored ``uses_full_state`` fails loudly instead of missing the
        root. Roots are located from sign changes [1]_.

        References
        ----------
        .. [1] SCIPY-1.14.1-IVP. SciPy 1.14.1, "solve_ivp": event sign-change
               detection, terminal and direction attributes.
        """

        threshold = self._bridge.protected_dry_mass_kg
        if threshold is None:
            return ()
        limit = float(threshold)

        def dry_mass(_t: float, y: Any) -> float:
            if len(y) != 7:
                raise DomainError("The dry-mass event needs the 7-state with mass at y[6].")
            return float(y[6]) - limit

        dry_mass.__name__ = "burn_dry_mass"
        dry_mass.terminal = True  # type: ignore[attr-defined]
        dry_mass.direction = -1.0  # type: ignore[attr-defined]
        dry_mass.uses_full_state = True  # type: ignore[attr-defined]
        return (dry_mass,)

    def propagate(
        self,
        y0: Sequence[float],
        cfg: Any,
        *,
        time_cfg: Any,
        extra_events: Sequence[Callable[[float, Any], float]] = (),
    ) -> SideraNativeBurnRun:
        """Propagate coast, burn segments and coast with Sidera's public ``propagate``.

        Parameters
        ----------
        y0 : sequence of float
            Initial 7-state ``[r (m), v (m/s), mass (kg)]`` at ``time_cfg.t0_s``.
        cfg : object
            Sidera ``PropagatorConfig``; used unchanged for every segment.
        time_cfg : object
            Sidera ``TimeConfig`` dataclass with an explicit ``output_dt_s``;
            each segment uses a copy with its own ``t0_s`` and ``duration_s``.
        extra_events : sequence of callable, optional
            Caller events, passed to every segment; a terminal one ends the
            run.

        Returns
        -------
        SideraNativeBurnRun
            Segments, honoured events and mass ledger.

        Raises
        ------
        DomainError
            If the state is not a finite 7-state with positive mass, the
            time configuration is unusable, the burn window is outside the
            run, the burn would consume the whole vehicle without a dry-mass
            limit, or a segment does not end at its boundary.

        Notes
        -----
        A dry-mass root stops thrust at the located time; the rest of the
        run coasts on the base dynamics. Burn start, interior boundaries and
        stop are recorded at their canonical dates [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        .. [2] SCIPY-1.14.1-IVP. SciPy 1.14.1, "solve_ivp": event sign-change
               detection, terminal and direction attributes.
        """

        numpy = self._modules.numpy
        propagate = self._modules.propagate
        state = numpy.asarray(y0, dtype=numpy.float64).reshape(-1).copy()
        if state.size != 7 or not bool(numpy.all(numpy.isfinite(state))) or state[6] <= 0.0:
            raise DomainError("A native Sidera burn needs a finite 7-state with positive mass.")
        if getattr(time_cfg, "output_dt_s", None) is None:
            raise DomainError("Native Sidera burn segments need an explicit output_dt_s.")
        t0 = real(getattr(time_cfg, "t0_s", 0.0) or 0.0, "Sidera start time")
        tf = t0 + real(getattr(time_cfg, "duration_s", None), "Sidera duration")
        user_events = list(extra_events)
        ignition, stop = self._dates[0], self._dates[-1]
        if ignition < t0 or stop > tf:
            raise DomainError("The burn window must lie inside the propagation window.")

        segments: list[tuple[str, Any]] = []
        events: list[SideraNativeEventRecord] = []

        def run(dynamics: Any, start: Any, lo: float, hi: float, label: str, hooks: list) -> Any:
            window = replace(time_cfg, t0_s=lo, duration_s=hi - lo)
            result = propagate(dynamics, start, cfg, time_cfg=window, extra_events=hooks or None)
            segments.append((label, result))
            return result

        def finish(**fields: Any) -> SideraNativeBurnRun:
            return SideraNativeBurnRun(
                segments=tuple(segments),
                events=tuple(events),
                ignition_date_s=ignition,
                stop_date_s=stop,
                bridge_model_hash=self._bridge.model_hash,
                **fields,
            )

        def reached(result: Any, hi: float) -> Any:
            end = float(result.t[-1])
            if abs(end - hi) > _SEGMENT_END_TOLERANCE_S:
                raise DomainError(f"Sidera segment ended at t={end} s instead of {hi} s.")
            return numpy.asarray(result.y[-1], dtype=numpy.float64).copy()

        def record(surface_id: str, kind: BurnEventKind, time_s: float, mass: float) -> None:
            events.append(SideraNativeEventRecord(surface_id, kind.value, time_s, mass))

        if ignition > t0:
            result = run(self._base, state, t0, ignition, "coast:pre-burn", user_events)
            if result.stopped_early:
                return finish(
                    cutoff_date_s=None,
                    mass_at_ignition_kg=None,
                    mass_at_cutoff_kg=None,
                    completed=False,
                    stopped_early=True,
                    stop_reason=result.stop_reason,
                )
            state = reached(result, ignition)
        mass_at_ignition = float(state[6])
        dry = self._bridge.protected_dry_mass_kg
        consumption = (
            self._bridge.artifact.consumed_propellant_kg * self._overrides["mass_flow_scale"]
        )
        if dry is None and not mass_at_ignition - consumption > 0.0:
            raise DomainError(
                f"The burn consumes {consumption!r} kg of a {mass_at_ignition!r} kg vehicle "
                "and the bridge has no protected_dry_mass_kg."
            )
        record("burn:start", BurnEventKind.START, ignition, mass_at_ignition)
        inventory = list(self.inventory_events())
        cutoff = stop
        completed = True
        if dry is not None and mass_at_ignition <= dry + INVENTORY_LIMIT_TOLERANCE:
            record("burn:dry-mass", BurnEventKind.DRY_MASS, ignition, mass_at_ignition)
            cutoff, completed = ignition, False
        else:
            last = len(self._dates) - 2
            for index in range(last + 1):
                lo, hi = self._dates[index], self._dates[index + 1]
                result = run(
                    self.segment_dynamics(index),
                    state,
                    lo,
                    hi,
                    f"burn:segment:{index + 1}",
                    user_events + inventory,
                )
                if result.stopped_early:
                    # Extra events are appended after Sidera's own, so the
                    # inventory columns are the last ones; the root time and
                    # state come from the event record itself.
                    first = len(result.t_events) - len(inventory)
                    fired = [
                        column
                        for column in range(first, len(result.t_events))
                        if inventory and numpy.asarray(result.t_events[column]).size
                    ]
                    if not fired:
                        state = numpy.asarray(result.y[-1], dtype=numpy.float64).copy()
                        return finish(
                            cutoff_date_s=float(result.t[-1]),
                            mass_at_ignition_kg=mass_at_ignition,
                            mass_at_cutoff_kg=float(state[6]),
                            completed=False,
                            stopped_early=True,
                            stop_reason=result.stop_reason,
                        )
                    column = fired[0]
                    cutoff = float(numpy.asarray(result.t_events[column]).reshape(-1)[0])
                    state = numpy.asarray(
                        result.y_events[column], dtype=numpy.float64
                    ).reshape(-1, 7)[0].copy()
                    completed = False
                    record("burn:dry-mass", BurnEventKind.DRY_MASS, cutoff, float(state[6]))
                    break
                state = reached(result, hi)
                if index < last:
                    record(
                        f"burn:boundary:{index + 1}",
                        BurnEventKind.PROFILE_BOUNDARY,
                        hi,
                        float(state[6]),
                    )
            if completed:
                record("burn:stop", BurnEventKind.STOP, stop, float(state[6]))
        mass_at_cutoff = float(state[6])
        stopped_early, stop_reason = False, None
        if cutoff < tf:
            result = run(self._base, state, cutoff, tf, "coast:post-burn", user_events)
            stopped_early, stop_reason = bool(result.stopped_early), result.stop_reason
        return finish(
            cutoff_date_s=cutoff,
            mass_at_ignition_kg=mass_at_ignition,
            mass_at_cutoff_kg=mass_at_cutoff,
            completed=completed,
            stopped_early=stopped_early,
            stop_reason=stop_reason,
        )


__all__ = [
    "SIDERA_NATIVE_FRAMES",
    "SIDERA_NATIVE_REQUIRED_CAPABILITIES",
    "SideraBridgeSegmentDynamics",
    "SideraNativeBurn",
    "SideraNativeBurnRun",
    "SideraNativeCapabilities",
    "SideraNativeEventMapping",
    "SideraNativeEventRecord",
    "inspect_sidera_native_capabilities",
    "missing_sidera_native_capabilities",
    "require_sidera_native_capabilities",
]
