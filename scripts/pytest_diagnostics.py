"""Best-effort diagnostic snapshots from a GIL-held Python timer callback.

The parent runner owns the independent hard process deadline. This plugin
does not alter assertions, test selection or that deadline. Unlike a C
faulthandler watchdog, the callback acquires the GIL before inspecting frames.
Upstream background: python/cpython#116008 and #131580.
"""

import faulthandler
import math
import os
import sys
from threading import Timer

import pytest


def pytest_addoption(parser):
    """Configure the diagnostic delay, independently of the process timeout."""
    parser.addoption("--diagnostic-timeout", type=float, default=20.,
                     help="Best-effort Python-thread stack snapshot delay in seconds")


def pytest_configure(config):
    """Keep a live diagnostic descriptor outside pytest's per-test capture."""
    delay = config.getoption("diagnostic_timeout")
    if not math.isfinite(delay) or delay <= 0:
        raise pytest.UsageError("Diagnostic timeout must be finite and positive.")
    if hasattr(sys, "_is_gil_enabled") and not sys._is_gil_enabled():
        raise pytest.UsageError("This diagnostic protocol requires a GIL-enabled runtime.")
    config._rpl_diagnostic_stream = os.fdopen(os.dup(2), "w", encoding="utf-8")


def _snapshot(stream, delay):
    """Run from Python with the GIL; a held native GIL may delay this callback."""
    stream.write(f"Diagnostic snapshot after {delay:g}s (Python thread holds GIL)\n")
    stream.flush()
    faulthandler.dump_traceback(file=stream, all_threads=True)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    """Bound the lifetime of one snapshot timer to one complete test protocol."""
    delay = item.config.getoption("diagnostic_timeout")
    timer = Timer(delay, _snapshot, args=(item.config._rpl_diagnostic_stream, delay))
    timer.daemon = True
    timer.start()
    try:
        yield
    finally:
        timer.cancel()
        timer.join(timeout=1.)


def pytest_unconfigure(config):
    """Close the duplicate after all per-test timers have been cancelled."""
    stream = getattr(config, "_rpl_diagnostic_stream", None)
    if stream is not None:
        stream.close()
