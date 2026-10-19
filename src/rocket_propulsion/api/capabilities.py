"""A nonblocking, process-local snapshot of optional Sidera inspection.

The synchronous integration inspector remains available to Python callers.
HTTP clients see an explicit unknown state while the optional modules load.
"""

from collections.abc import Callable
from dataclasses import asdict
from threading import Lock, Thread
from typing import Any

from rocket_propulsion.integrations.contracts import SideraCapabilities


def _failed_probe(error: BaseException) -> SideraCapabilities:
    """Represent an inspection failure without claiming any observed support."""
    return SideraCapabilities(
        status="incompatible", version=None, finite_burn_available=False,
        maneuver_plan_available=False, constructor_compatible=False, supported_frames=(),
        constant_burn_exact=False, tabulated_burn_available=False,
        independent_mass_flow_available=False,
        reason=f"Sidera inspection failed: {type(error).__name__}: {error}",
    )


class SideraCapabilityProbe:
    """Inspect once per process; pending flags remain unknown until completion.

    Installation/configuration changes require restarting the server. A hung
    external import remains checking; requests never wait for its completion
    or start replacement workers. The worker is daemonized for server shutdown.
    """

    def __init__(self, inspector: Callable[[], SideraCapabilities]) -> None:
        """Retain the inspector without importing optional native packages."""
        self._inspector = inspector
        self._lock = Lock()
        self._started = False
        self._result: SideraCapabilities | None = None

    def _run(self) -> None:
        """Publish one immutable result, failing closed on optional-package errors."""
        try:
            result = self._inspector()
            if not isinstance(result, SideraCapabilities):
                raise TypeError("Inspector did not return a SideraCapabilities record.")
        except BaseException as error:  # worker exits must publish a terminal result
            result = _failed_probe(error)
        with self._lock:
            self._result = result

    def snapshot(self) -> dict[str, Any]:
        """Return immediately, without substituting false flags for pending evidence."""
        with self._lock:
            if self._result is None and not self._started:
                self._started = True
                try:
                    Thread(target=self._run, name="sidera-capability-probe", daemon=True).start()
                except RuntimeError as error:
                    self._result = _failed_probe(error)
            if self._result is not None:
                return {**asdict(self._result), "probe_complete": True}
            return {"status": "checking", "probe_complete": False,
                    "reason": "Optional Sidera modules are being inspected; capabilities are unknown."}
