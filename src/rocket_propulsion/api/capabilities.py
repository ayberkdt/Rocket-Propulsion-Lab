"""Nonblocking, process-local snapshot of the optional Sidera inspection.

Python callers can call the synchronous inspector
:func:`~rocket_propulsion.integrations.sidera.inspect_sidera_capabilities`
directly. HTTP clients instead receive an explicit unknown state while the
optional modules load in a background thread [1]_.

References
----------
.. [1] Python 3.12, "threading - Thread-based parallelism" (Thread, Lock and
       daemon semantics), https://docs.python.org/3.12/library/threading.html
"""

from collections.abc import Callable
from dataclasses import asdict
from threading import Lock, Thread
from typing import Any

from rocket_propulsion.integrations.contracts import SideraCapabilities


def _failed_probe(error: BaseException) -> SideraCapabilities:
    """Return an ``"incompatible"`` record with every capability flag false.

    The reason names the exception type and message, so a failed inspection
    never claims observed support.
    """
    return SideraCapabilities(
        status="incompatible", version=None, finite_burn_available=False,
        maneuver_plan_available=False, constructor_compatible=False, supported_frames=(),
        constant_burn_exact=False, tabulated_burn_available=False,
        independent_mass_flow_available=False,
        reason=f"Sidera inspection failed: {type(error).__name__}: {error}",
    )


class SideraCapabilityProbe:
    """Run an inspector once per process and serve its result without blocking.

    Parameters
    ----------
    inspector : callable
        Zero-argument function returning a :class:`SideraCapabilities`
        record; it is called at most once, in a daemon worker thread.

    Notes
    -----
    Until the worker finishes, capability flags are reported as unknown,
    never as false. Installation or configuration changes take effect only
    after a server restart. A hung external import stays in the checking
    state: requests never wait for it and never start replacement workers.
    The worker is a daemon thread so that it cannot block interpreter
    shutdown [1]_.

    References
    ----------
    .. [1] Python 3.12, "threading - Thread-based parallelism" (Thread, Lock and
           daemon semantics), https://docs.python.org/3.12/library/threading.html
    """

    def __init__(self, inspector: Callable[[], SideraCapabilities]) -> None:
        """Store the inspector; no optional package is imported here."""
        self._inspector = inspector
        self._lock = Lock()
        self._started = False
        self._result: SideraCapabilities | None = None

    def _run(self) -> None:
        """Run the inspector once and publish its result under the lock.

        Any exception, including ``SystemExit`` raised in the worker, is
        converted into a failed-inspection record, so a terminal result is
        always published.
        """
        try:
            result = self._inspector()
            if not isinstance(result, SideraCapabilities):
                raise TypeError("Inspector did not return a SideraCapabilities record.")
        except BaseException as error:  # noqa: BLE001 - Even worker SystemExit must publish a terminal result.
            result = _failed_probe(error)
        with self._lock:
            self._result = result

    def snapshot(self) -> dict[str, Any]:
        """Return the current inspection state immediately.

        The first call starts the worker thread. If the thread cannot be
        started, the failure itself becomes the result.

        Returns
        -------
        dict
            The capability record plus ``probe_complete=True`` once the
            inspection has finished; otherwise
            ``{"status": "checking", "probe_complete": False, "reason": ...}``
            with no capability flags.

        References
        ----------
        .. [1] Python 3.12, "threading - Thread-based parallelism" (Thread, Lock and
               daemon semantics), https://docs.python.org/3.12/library/threading.html
        """
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
