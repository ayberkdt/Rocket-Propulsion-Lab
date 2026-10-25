"""Pulse-mode thruster: impulse bits, minimum impulse bit and Δv quantization.

Attitude-control and small trajectory-correction thrusters are fired in short
electrical pulses.  For short pulses the delivered impulse is not
``F_steady * t_on``: valve opening, chamber fill and tail-off make the impulse
bit (and the specific impulse) a nonlinear function of the electrical on-time
[2]_, [3]_, and the smallest pulse the valve can command sets the **minimum
impulse bit** (MIB) [1]_.  A trajectory tool therefore cannot ask for an
arbitrarily small Δv: it receives an integer number of pulses whose total may
differ from the request.

:class:`PulseModeThruster` takes caller-supplied data, an impulse-bit table
``I(t_on)`` and a propellant-per-pulse table ``m(t_on)``, and never invents
generic pulse coefficients.  Between table nodes it uses shape-preserving
PCHIP interpolation (monotone, no overshoot) [4]_.  Beyond the last node the
pulse is fully developed, and extra on-time adds steady-state impulse and flow
linearly.

The tables are intended to represent characterized hardware.  This class does
not authenticate measurements or establish empirical validity; a model is only
as valid as the table it is given.

For ``n`` identical pulses of on-time ``t`` from mass ``m0`` in free space,
every pulse has the same exhaust velocity ``c = I/m``, so the ideal rocket
equation [5]_ holds exactly for the whole train:

    Δv_n(t) = (I(t) / m(t)) * ln(m0 / (m0 - n m(t))).

In exact arithmetic ``Δv_n`` is strictly increasing in ``t`` (``I`` strictly
increasing, ``m`` non-decreasing) and in ``n``, so a requested Δv is met by
the smallest ``n`` whose longest permitted pulse reaches it, and then by
bisection on ``t``.  When even ``n`` minimum impulse bits overshoot, the
request is quantized and the residual is reported, never hidden.

References
----------
.. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
       NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
.. [2] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
       https://ntrs.nasa.gov/citations/19720011120
.. [3] NASA, "Space Simulation Experiments on Reaction Control System
       Thruster Plumes". https://ntrs.nasa.gov/citations/19730002016
.. [4] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
       Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
.. [5] NASA Glenn Research Center, "Ideal Rocket Equation".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import pairwise
from math import ceil, floor, log1p
from threading import Lock

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    integer,
    real,
    real_tuple,
    strings,
    text,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .exchange_documents import _read_document
from .serialization import _require_numerical_policies, canonical_sha256
from .tabulated import InterpolationPolicy, TabulatedBurnArtifact, TabulatedBurnSegment
from .throttleable import MonotoneCubicCurve

PULSE_MODE_SCHEMA = "rocket_propulsion_pulse_mode_thruster_v1"
PULSE_MODE_REFERENCE_IDS = ("NASA-20070032005",)
QUANTIZATION_POLICIES = ("nearest", "no_overshoot")
# A plan whose delivered value lies within this relative distance of the
# request meets it (``quantized`` is ``False``).
_RELATIVE_TOLERANCE = 1e-9
_MAX_PULSES = 10_000_000
_WIDTH_CACHE_LOCK = Lock()


@dataclass(frozen=True, slots=True)
class PulseBit:
    """Impulse, propellant and their on-time derivatives for one pulse.

    Pulse quantities are indexed by electrical on-time; [1]_ relates the
    minimum impulse bit to valve response and electrical pulse width.

    Attributes
    ----------
    on_time_s : float
        Electrical on-time (pulse width), s.
    impulse_n_s : float
        Impulse bit, N s.
    propellant_kg : float
        Propellant consumed by the pulse, kg.
    impulse_wrt_on_time_n : float
        ``dI/dt_on``, N.
    propellant_wrt_on_time_kg_s : float
        ``dm/dt_on``, kg/s.

    References
    ----------
    .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
           NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
    """

    on_time_s: float
    impulse_n_s: float
    propellant_kg: float
    impulse_wrt_on_time_n: float
    propellant_wrt_on_time_kg_s: float

    @property
    def specific_impulse_s(self) -> float:
        """Pulse-averaged specific impulse ``I / (m g0)``, s.

        ``g0`` is the registered ``STANDARD_GRAVITY_M_S2`` [1]_.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """

        return self.impulse_n_s / (self.propellant_kg * STANDARD_GRAVITY_M_S2)

    @property
    def exhaust_velocity_m_s(self) -> float:
        """Pulse-averaged effective exhaust velocity ``I / m``, m/s [1]_.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """

        return self.impulse_n_s / self.propellant_kg


@dataclass(frozen=True, slots=True)
class PulsePlan:
    """An integer pulse train that meets or quantizes a request.

    Attributes
    ----------
    target : str
        ``"impulse"`` or ``"delta_v"``.
    requested : float
        Requested impulse in N s or Δv in m/s, as given by ``target``.
    delivered : float
        Delivered value in the same unit.
    residual : float
        ``requested - delivered``; zero to solver tolerance unless
        ``quantized`` is ``True``.
    quantized : bool
        ``True`` when ``|residual|`` exceeds a relative ``1e-9`` of the
        request, i.e. no integer train meets it.
    pulse_count : int
        Number of identical pulses.
    pulse : PulseBit or None
        The repeated pulse; ``None`` for an empty train.
    delivered_impulse_n_s : float
        Total impulse of the train, N s.
    propellant_kg : float
        Total propellant of the train, kg.
    initial_mass_kg, final_mass_kg : float or None
        Vehicle mass before and after the train for a Δv plan; ``None`` for
        an impulse plan.

    Notes
    -----
    Quantization arises from the minimum impulse bit [1]_.

    References
    ----------
    .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
           NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
    """

    target: str
    requested: float
    delivered: float
    residual: float
    quantized: bool
    pulse_count: int
    pulse: PulseBit | None
    delivered_impulse_n_s: float
    propellant_kg: float
    initial_mass_kg: float | None
    final_mass_kg: float | None

    @property
    def effective_specific_impulse_s(self) -> float | None:
        """Train-averaged specific impulse, s, or ``None`` without pulses.

        All pulses are identical, so the train average equals the single-pulse
        value ``I / (m g0)`` [1]_.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """

        return None if self.pulse is None else self.pulse.specific_impulse_s


@dataclass(frozen=True, slots=True)
class PulseModeThruster:
    """Caller-supplied pulse-mode characteristics of one thruster.

    Parameters
    ----------
    on_times_s : tuple of float
        Positive, strictly increasing electrical pulse widths, s; at least two
        nodes.  The first is the minimum on-time.
    impulse_bits_n_s : tuple of float
        Positive, strictly increasing impulse bit per node, N s; the first is
        the minimum impulse bit.
    propellant_per_pulse_kg : tuple of float
        Positive, non-decreasing propellant per pulse per node, kg.
    steady_thrust_n : float
        Positive steady-state thrust, N.
    steady_mass_flow_kg_s : float
        Positive steady-state mass flow, kg/s.
    minimum_off_time_s : float
        Non-negative minimum time between pulses, s.
    source_id : str
        Non-blank identifier of the characterization data.
    source_sha256 : str
        SHA-256 hex digest of the characterization data.
    maximum_on_time_s : float or None, default None
        Longest permitted pulse, s, not shorter than the last node; ``None``
        uses the last node.
    model_id : str, default "pulse-mode-thruster"
        Non-blank model identifier.
    reference_ids : tuple of str, default PULSE_MODE_REFERENCE_IDS
        Source identifiers recorded with the model.
    schema : str, default PULSE_MODE_SCHEMA
        Document schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If the schema is unsupported, a field has the wrong type, an
        identifier is blank, the hash is not a SHA-256 digest, or a table or
        scalar violates the ordering and sign requirements above.

    Notes
    -----
    The minimum on-time and its impulse bit characterize the valve-limited
    minimum impulse bit [1]_.  Between nodes, impulse and propellant are
    interpolated by monotone PCHIP [2]_, [3]_.  Beyond the last node a pulse
    is fully developed: ``I = I_N + F (t - t_N)`` and
    ``m = m_N + mdot (t - t_N)`` with the steady thrust and flow, up to the
    longest permitted pulse.  The values are continuous at ``t_N``; the slope
    may change.

    References
    ----------
    .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
           NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
    .. [2] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    .. [3] F. N. Fritsch and J. Butland, SIAM J. Sci. Stat. Comput. 5(2), 1984.
    """

    on_times_s: tuple[float, ...]
    impulse_bits_n_s: tuple[float, ...]
    propellant_per_pulse_kg: tuple[float, ...]
    steady_thrust_n: float
    steady_mass_flow_kg_s: float
    minimum_off_time_s: float
    source_id: str
    source_sha256: str
    maximum_on_time_s: float | None = None
    model_id: str = "pulse-mode-thruster"
    reference_ids: tuple[str, ...] = PULSE_MODE_REFERENCE_IDS
    schema: str = PULSE_MODE_SCHEMA
    _impulse: MonotoneCubicCurve = field(init=False, repr=False, compare=False)
    _propellant: MonotoneCubicCurve = field(init=False, repr=False, compare=False)
    _width_cache: dict[float, float] = field(default_factory=dict, init=False,
                                           repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.schema != PULSE_MODE_SCHEMA:
            raise DomainError("Unsupported pulse-mode thruster schema.")
        for name in ("on_times_s", "impulse_bits_n_s", "propellant_per_pulse_kg"):
            object.__setattr__(self, name, real_tuple(getattr(self, name), f"Pulse {name}"))
        for name in ("steady_thrust_n", "steady_mass_flow_kg_s", "minimum_off_time_s"):
            object.__setattr__(self, name, real(getattr(self, name), f"Pulse {name}"))
        if self.maximum_on_time_s is not None:
            object.__setattr__(
                self, "maximum_on_time_s", real(self.maximum_on_time_s, "Maximum on-time")
            )
        object.__setattr__(self, "model_id", text(self.model_id, "Pulse model id").strip())
        object.__setattr__(self, "source_id", text(self.source_id, "Pulse source id").strip())
        text(self.source_sha256, "Pulse source hash")
        object.__setattr__(
            self, "reference_ids", strings(self.reference_ids, "Pulse reference ids")
        )
        if not self.model_id or not self.source_id:
            raise DomainError("Pulse model and source identifiers must not be empty.")
        if len(self.source_sha256) != 64 or any(
            character not in "0123456789abcdef" for character in self.source_sha256.lower()
        ):
            raise DomainError("Pulse source hash must be a SHA-256 hex digest.")
        times, impulses, masses = self.on_times_s, self.impulse_bits_n_s, self.propellant_per_pulse_kg
        if len(times) < 2 or len(times) != len(impulses) or len(times) != len(masses):
            raise DomainError("Pulse tables need at least two nodes and equal lengths.")
        if times[0] <= 0.0 or any(b <= a for a, b in pairwise(times)):
            raise DomainError("On-times must be positive and strictly increasing.")
        if impulses[0] <= 0.0 or any(b <= a for a, b in pairwise(impulses)):
            raise DomainError("Impulse bits must be positive and strictly increasing.")
        if masses[0] <= 0.0 or any(b < a for a, b in pairwise(masses)):
            raise DomainError("Propellant per pulse must be positive and non-decreasing.")
        if self.steady_thrust_n <= 0.0 or self.steady_mass_flow_kg_s <= 0.0:
            raise DomainError("Steady thrust and mass flow must be greater than zero.")
        if self.minimum_off_time_s < 0.0:
            raise DomainError("Minimum off-time must be non-negative.")
        if self.maximum_on_time_s is not None and self.maximum_on_time_s < times[-1]:
            raise DomainError("Maximum on-time must not be shorter than the last table node.")
        object.__setattr__(self, "_impulse", MonotoneCubicCurve(times, impulses))
        object.__setattr__(self, "_propellant", MonotoneCubicCurve(times, masses))

    def _solve_budget_width(self, share: float) -> float:
        """Return the largest on-time whose propellant per pulse fits ``share``.

        Bisection on the non-decreasing propellant curve.  The result depends
        only on this immutable model and ``share`` and is memoized in a
        thread-safe cache of at most 512 entries per model that is excluded
        from manifests, hashes and equality.
        """
        with _WIDTH_CACHE_LOCK:
            cached = self._width_cache.get(share)
            if cached is not None:
                return cached
        shortest, longest = self.minimum_on_time_s, self.longest_on_time_s
        low, high = shortest, longest
        last = self.on_times_s[-1]
        for _ in range(200):
            middle = 0.5 * (low + high)
            if middle in (low, high):
                break
            if middle <= last:
                mass = self._propellant._value_checked(middle)
            else:
                mass = (self.propellant_per_pulse_kg[-1]
                        + self.steady_mass_flow_kg_s * (middle - last))
            if mass <= share:
                low = middle
            else:
                high = middle
        with _WIDTH_CACHE_LOCK:
            if len(self._width_cache) >= 512:
                self._width_cache.pop(next(iter(self._width_cache)))
            self._width_cache[share] = low
        return low

    @property
    def minimum_on_time_s(self) -> float:
        """Minimum electrical on-time, the first table node, s."""
        return self.on_times_s[0]

    @property
    def longest_on_time_s(self) -> float:
        """Longest permitted pulse, ``maximum_on_time_s`` or the last node, s."""

        return self.on_times_s[-1] if self.maximum_on_time_s is None else self.maximum_on_time_s

    @property
    def minimum_impulse_bit_n_s(self) -> float:
        """Minimum impulse bit [1]_, the impulse at the minimum on-time, N s.

        References
        ----------
        .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        """
        return self.impulse_bits_n_s[0]

    @property
    def exhaust_velocity_bound_m_s(self) -> float:
        """Upper bound of the exhaust velocity ``I(t)/m(t)`` over every permitted width.

        PCHIP is monotone between monotone nodes [1]_, so on ``[t_k, t_k+1]``
        ``I <= I_k+1`` and ``m >= m_k``, and the bound there is
        ``I_k+1 / m_k``.  Beyond the last node the ratio moves monotonically
        from ``I_N/m_N`` towards ``F/mdot``.

        References
        ----------
        .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        """

        impulses, masses = self.impulse_bits_n_s, self.propellant_per_pulse_kg
        bound = max(impulses[k + 1] / masses[k] for k in range(len(masses) - 1))
        bound = max(bound, impulses[-1] / masses[-1])
        if self.longest_on_time_s > self.on_times_s[-1]:
            bound = max(bound, self.steady_thrust_n / self.steady_mass_flow_kg_s)
        return bound

    def pulse(self, on_time_s: float) -> PulseBit:
        """Return impulse, propellant and their on-time derivatives for one pulse.

        Parameters
        ----------
        on_time_s : float
            Electrical on-time, s, between the minimum and the longest
            permitted pulse; values within ``1e-12 * max(1, t_longest)`` of
            either limit are clamped to it.

        Returns
        -------
        PulseBit
            PCHIP-interpolated values [2]_ at or below the last node, the
            fully developed linear extension beyond it.

        Raises
        ------
        DomainError
            If the on-time is non-finite or outside the permitted range,
            including below the minimum on-time that sets the minimum impulse
            bit [1]_.

        References
        ----------
        .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        .. [2] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        """

        on_time = real(on_time_s, "On-time")
        tolerance = 1e-12 * max(1.0, self.longest_on_time_s)
        if on_time < self.minimum_on_time_s - tolerance:
            raise DomainError(
                f"On-time {on_time} s is below the minimum on-time {self.minimum_on_time_s} s."
            )
        if on_time > self.longest_on_time_s + tolerance:
            raise DomainError(
                f"On-time {on_time} s exceeds the longest permitted pulse "
                f"{self.longest_on_time_s} s."
            )
        on_time = min(max(on_time, self.minimum_on_time_s), self.longest_on_time_s)
        return PulseBit(on_time, *self._raw(on_time))

    def _raw(self, on_time: float) -> tuple[float, float, float, float]:
        """Return impulse, propellant and slopes at an already validated on-time.

        PCHIP values and derivatives [1]_ up to the last node; beyond it, the
        linear steady-state extension ``I_N + F (t - t_N)``,
        ``m_N + mdot (t - t_N)`` with slopes ``F`` and ``mdot``.

        References
        ----------
        .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        """

        last = self.on_times_s[-1]
        if on_time <= last:
            impulse, impulse_slope = self._impulse._evaluate_checked(on_time)
            mass, mass_slope = self._propellant._evaluate_checked(on_time)
            return impulse, mass, impulse_slope, mass_slope
        extra = on_time - last
        return (
            self.impulse_bits_n_s[-1] + self.steady_thrust_n * extra,
            self.propellant_per_pulse_kg[-1] + self.steady_mass_flow_kg_s * extra,
            self.steady_thrust_n,
            self.steady_mass_flow_kg_s,
        )

    def _raw_values(self, on_time: float) -> tuple[float, float]:
        """Return impulse and propellant only, at a validated on-time.

        Used inside planning searches; the final `PulseBit` of a plan is
        evaluated by `pulse`, which also returns both slopes.
        """
        last = self.on_times_s[-1]
        if on_time <= last:
            return self._impulse._value_checked(on_time), self._propellant._value_checked(on_time)
        extra = on_time - last
        return (self.impulse_bits_n_s[-1] + self.steady_thrust_n * extra,
                self.propellant_per_pulse_kg[-1] + self.steady_mass_flow_kg_s * extra)

    def on_time_for_impulse(self, impulse_n_s: float) -> float:
        """Return the on-time that delivers a given impulse bit.

        Parameters
        ----------
        impulse_n_s : float
            Requested impulse bit, N s, between the minimum impulse bit and the
            impulse of the longest permitted pulse (relative tolerance
            ``1e-12``).

        Returns
        -------
        float
            On-time, s, found by bisection on the strictly increasing impulse
            curve.

        Raises
        ------
        DomainError
            If the impulse is non-finite, below the minimum impulse bit [1]_,
            or above the longest pulse's impulse.

        References
        ----------
        .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        """

        impulse = real(impulse_n_s, "Impulse bit")
        low, high = self.minimum_on_time_s, self.longest_on_time_s
        if impulse < self.minimum_impulse_bit_n_s * (1.0 - 1e-12):
            raise DomainError("Impulse is below the minimum impulse bit.")
        if impulse > self.pulse(high).impulse_n_s * (1.0 + 1e-12):
            raise DomainError("Impulse exceeds the longest permitted pulse.")
        return _bisect(lambda t: self.pulse(t).impulse_n_s - impulse, low, high)

    def train_delta_v_m_s(self, pulse_count: int, on_time_s: float, initial_mass_kg: float) -> float:
        """Return the free-space Δv of ``pulse_count`` identical pulses.

        Parameters
        ----------
        pulse_count : int
            Non-negative number of pulses.
        on_time_s : float
            On-time of every pulse, s.
        initial_mass_kg : float
            Positive vehicle mass before the train, kg.

        Returns
        -------
        float
            ``(I/m) * ln(m0 / (m0 - n m))``, m/s, evaluated as
            ``-(I/m) * log1p(-n m / m0)``; zero for ``pulse_count == 0``.

        Raises
        ------
        DomainError
            If the count is negative, the mass is not positive, the on-time is
            outside the permitted range, or the train consumes at least the
            whole initial mass.

        Notes
        -----
        Every pulse has the same effective exhaust velocity ``I/m``, so the
        ideal rocket equation [1]_ is exact for the whole train.  Gravity,
        drag and the time between pulses are ignored.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
        """

        count = integer(pulse_count, "Pulse count")
        mass = real(initial_mass_kg, "Initial mass")
        if count < 0 or mass <= 0.0:
            raise DomainError("Pulse count must be non-negative and initial mass positive.")
        if count == 0:
            return 0.0
        bit = self.pulse(on_time_s)
        consumed = count * bit.propellant_kg
        if consumed >= mass:
            raise DomainError("The pulse train consumes at least the whole initial mass.")
        return bit.exhaust_velocity_m_s * -log1p(-consumed / mass)

    def plan_impulse(
        self,
        impulse_n_s: float,
        *,
        quantization: str = "nearest",
        maximum_pulses: int | None = None,
    ) -> PulsePlan:
        """Return the fewest equal pulses whose total impulse meets a request.

        Parameters
        ----------
        impulse_n_s : float
            Non-negative requested total impulse, N s.
        quantization : {"nearest", "no_overshoot"}, default "nearest"
            Choice of train when no integer train meets the request: the
            nearest achievable total, or the largest one not above the request.
        maximum_pulses : int or None, default None
            Upper bound on the pulse count; ``None`` uses the module limit.

        Returns
        -------
        PulsePlan
            Pulse count, pulse and delivered impulse with residual.

        Raises
        ------
        DomainError
            If the request is negative, the policy or limit is invalid, or the
            plan needs more pulses than the limit.

        Notes
        -----
        The count is the smallest ``n`` with ``n * I(t_longest)`` reaching the
        request; the on-time then solves ``I(t) = request / n``.  A request
        between ``(n - 1)`` longest pulses and ``n`` minimum impulse bits [1]_
        cannot be met exactly and is quantized.

        References
        ----------
        .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        """

        request = real(impulse_n_s, "Requested impulse")
        if request < 0.0:
            raise DomainError("Requested impulse must be non-negative.")
        limit = self._pulse_limit(maximum_pulses)
        policy = _policy(quantization)
        longest = self.pulse(self.longest_on_time_s)
        minimum = self.pulse(self.minimum_on_time_s)
        if request == 0.0:
            return self._plan("impulse", 0.0, 0, None, None, limit)
        count = max(1, ceil(request * (1.0 - _RELATIVE_TOLERANCE) / longest.impulse_n_s))
        if request / count >= minimum.impulse_n_s * (1.0 - _RELATIVE_TOLERANCE):
            on_time = self.on_time_for_impulse(
                min(max(request / count, minimum.impulse_n_s), longest.impulse_n_s)
            )
            return self._plan("impulse", request, count, self.pulse(on_time), None, limit)
        # The request lies in the gap between (count - 1) longest pulses and
        # count minimum impulse bits; no integer train meets it exactly.
        below = (count - 1) * longest.impulse_n_s
        above = count * minimum.impulse_n_s
        if policy == "no_overshoot" or request - below < above - request:
            return self._plan(
                "impulse", request, count - 1, longest if count > 1 else None, None, limit
            )
        return self._plan("impulse", request, count, minimum, None, limit)

    def plan_delta_v(
        self,
        delta_v_m_s: float,
        initial_mass_kg: float,
        *,
        propellant_budget_kg: float,
        quantization: str = "nearest",
        maximum_pulses: int | None = None,
        maximum_search_pulses: int = 100_000,
    ) -> PulsePlan:
        """Return the fewest equal pulses whose free-space Δv meets a request.

        Parameters
        ----------
        delta_v_m_s : float
            Non-negative requested Δv, m/s.
        initial_mass_kg : float
            Positive vehicle mass before the train, kg.
        propellant_budget_kg : float
            Propellant available to the train, kg, in ``(0, initial_mass_kg)``.
        quantization : {"nearest", "no_overshoot"}, default "nearest"
            Choice of train when no integer train meets the request.
        maximum_pulses : int or None, default None
            Upper bound on the pulse count; ``None`` uses the module limit.
        maximum_search_pulses : int, default 100000
            Upper bound on the number of budget-limited counts scanned.

        Returns
        -------
        PulsePlan
            Pulse count, pulse, delivered Δv with residual, and the initial
            and final masses.

        Raises
        ------
        DomainError
            If an argument is outside its domain, the budget is below one
            minimum pulse, the request exceeds the Δv the budget can deliver,
            the scan exceeds ``maximum_search_pulses``, or the plan needs more
            pulses than the limit.

        Notes
        -----
        Every pulse train fits in ``propellant_budget_kg``, which is below
        ``initial_mass_kg`` because the dry structure is not propellant.  For
        ``n`` pulses the feasible widths are ``[t_min, t_max(n)]`` with
        ``t_max(n) = min(t_longest, m^-1(budget / n))``, and the achievable Δv
        from the ideal rocket equation [1]_ is the interval
        ``[Δv_n(t_min), Δv_n(t_max(n))]`` (``Δv_n`` increases with the width).
        The lower ends increase with ``n``; up to the count the budget allows
        at the longest width the upper ends do too, so that part is bisected,
        and the budget-limited counts beyond it are scanned.  A request inside
        some interval is met exactly with the smallest such ``n``; otherwise it
        is quantized, as the minimum impulse bit [2]_ requires, to the nearest
        achievable train (or the largest one not above it with
        ``no_overshoot``).  A request above everything the budget can deliver
        is refused.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
        .. [2] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        """

        request = real(delta_v_m_s, "Requested delta-v")
        mass = real(initial_mass_kg, "Initial mass")
        budget = real(propellant_budget_kg, "Propellant budget")
        if request < 0.0 or mass <= 0.0:
            raise DomainError("Delta-v must be non-negative and initial mass positive.")
        if not 0.0 < budget < mass:
            raise DomainError("The propellant budget must be positive and below the initial mass.")
        limit = self._pulse_limit(maximum_pulses)
        scan_cap = integer(maximum_search_pulses, "Maximum search pulses")
        if scan_cap < 1:
            raise DomainError("Maximum search pulses must be at least one.")
        policy = _policy(quantization)
        if request == 0.0:
            return self._plan("delta_v", 0.0, 0, None, mass, limit)
        if request > self.exhaust_velocity_bound_m_s * -log1p(-budget / mass):
            raise DomainError(
                "The request exceeds the delta-v the propellant budget can deliver."
            )
        shortest, longest = self.minimum_on_time_s, self.longest_on_time_s
        shortest_mass = self.pulse(shortest).propellant_kg
        longest_mass = self.pulse(longest).propellant_kg
        most = floor(budget / shortest_mass * (1.0 + 1e-12))
        if most < 1:
            raise DomainError("The propellant budget is below one minimum pulse.")
        at_longest = min(most, floor(budget / longest_mass * (1.0 + 1e-12)))

        # Full pulse values are reused only within this planning call.
        raw = lru_cache(maxsize=256)(self._raw_values)

        def delta_v(count: int, on_time: float) -> float:
            impulse, propellant = raw(on_time)
            used = min(count * propellant, budget)
            return impulse / propellant * -log1p(-used / mass)

        def widest(count: int) -> float:
            if count <= at_longest:
                return longest
            return self._solve_budget_width(budget / count)

        def low_end(count: int) -> float:
            return delta_v(count, shortest)

        def high_end(count: int) -> float:
            return delta_v(count, widest(count))

        exact = None
        below: tuple[float, int, float] = (0.0, 0, shortest)  # (Δv, count, width)
        above: tuple[float, int] | None = None                 # (Δv, count) at t_min
        reach = request * (1.0 - _RELATIVE_TOLERANCE)
        if at_longest >= 1 and high_end(at_longest) >= reach:
            # Monotone part: smallest count whose longest pulse reaches the request.
            lower, upper = 0, at_longest
            while upper - lower > 1:
                middle = (lower + upper) // 2
                if high_end(middle) >= reach:
                    upper = middle
                else:
                    lower = middle
            if low_end(upper) <= request * (1.0 + _RELATIVE_TOLERANCE):
                exact = upper
            else:
                if upper > 1:
                    below = (high_end(upper - 1), upper - 1, longest)
                above = (low_end(upper), upper)
        else:
            if at_longest >= 1:
                below = (high_end(at_longest), at_longest, longest)
            count = at_longest + 1
            scanned = 0
            while count <= most:
                scanned += 1
                if scanned > scan_cap:
                    raise DomainError(
                        "The budget-limited pulse search exceeded maximum_search_pulses."
                    )
                if low_end(count) > request * (1.0 + _RELATIVE_TOLERANCE):
                    above = (low_end(count), count)
                    break
                top = high_end(count)
                if top >= reach:
                    exact = count
                    break
                if top > below[0]:
                    below = (top, count, widest(count))
                count += 1
            # At the very top of the range the best train may fall short of the
            # request by rounding only; that is a met request.
            if (exact is None and above is None
                    and request - below[0] > _RELATIVE_TOLERANCE * request):
                raise DomainError(
                    "The request exceeds the delta-v the propellant budget can deliver."
                )
        if exact is not None:
            top = widest(exact)
            on_time = _bisect(lambda t: delta_v(exact, t) - request, shortest, top)
            return self._plan("delta_v", request, exact, self.pulse(on_time), mass, limit)
        if above is None or policy == "no_overshoot" or request - below[0] < above[0] - request:
            return self._plan(
                "delta_v", request, below[1],
                self.pulse(below[2]) if below[1] else None, mass, limit,
            )
        return self._plan("delta_v", request, above[1], self.pulse(shortest), mass, limit)

    def pulse_train_artifact(
        self, plan: PulsePlan, period_s: float, *, phase: str = "pulse"
    ) -> TabulatedBurnArtifact:
        """Return a propagator-ready tabulated burn for a pulse train.

        Parameters
        ----------
        plan : PulsePlan
            Plan with at least one pulse produced by this model.
        period_s : float
            Pulse repetition period, s; longer than the on-time and leaving at
            least the minimum off-time (tolerance ``1e-12`` s).
        phase : str, default "pulse"
            Phase label of the on-segments.

        Returns
        -------
        TabulatedBurnArtifact
            Alternating on and ``"off"`` segments with hold interpolation and
            a source hash bound to this model's hash, the pulse count, the
            on-time and the period.

        Raises
        ------
        DomainError
            If ``plan`` is not a `PulsePlan`, has no pulse, carries a pulse not
            produced by this model, has more than 100000 pulses, or the period
            is too short.

        Notes
        -----
        Each pulse is a rectangle over its electrical on-time with constant
        thrust ``I/t_on`` and flow ``m/t_on``, so impulse, propellant and
        exhaust velocity per pulse are exact [1]_; the internal pulse shape
        (rise and tail-off [2]_) is not represented.  Component transient
        models such as [3]_ resolve that shape.

        References
        ----------
        .. [1] J. M. Parker and M. J. Wilson, "The Minimum Impulse Thruster", 2005,
               NTRS 20070032005. https://ntrs.nasa.gov/citations/20070032005
        .. [2] NASA CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
               https://ntrs.nasa.gov/citations/19720011120
        .. [3] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
               https://ntrs.nasa.gov/citations/19970010379
        """

        if not isinstance(plan, PulsePlan):
            raise DomainError("A PulsePlan is required.")
        count = integer(plan.pulse_count, "Pulse count")
        if count < 1 or not isinstance(plan.pulse, PulseBit):
            raise DomainError("A pulse train needs at least one pulse.")
        # Only a pulse this thruster actually produces may be stamped with its hash.
        if plan.pulse != self.pulse(plan.pulse.on_time_s):
            raise DomainError("The plan's pulse was not produced by this thruster model.")
        if count > 100_000:
            raise DomainError("Pulse trains above 100000 pulses are not tabulated.")
        period = real(period_s, "Pulse period")
        on_time = plan.pulse.on_time_s
        if period < on_time + self.minimum_off_time_s - 1e-12 or period <= on_time:
            raise DomainError("The pulse period leaves less than the minimum off-time.")
        phase = text(phase, "Pulse phase")
        thrust = plan.pulse.impulse_n_s / on_time
        flow = plan.pulse.propellant_kg / on_time
        segments: list[TabulatedBurnSegment] = []
        for index in range(count):
            start = index * period
            segments.append(
                TabulatedBurnSegment(start, start + on_time, thrust, thrust, thrust, thrust,
                                     flow, flow, phase=phase)
            )
            if index + 1 < count:
                segments.append(
                    TabulatedBurnSegment(start + on_time, (index + 1) * period,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phase="off")
                )
        source = canonical_sha256(
            {
                "model_hash": self.model_hash,
                "pulse_count": count,
                "on_time_s": on_time,
                "period_s": period,
            }
        )
        return TabulatedBurnArtifact(
            segments=tuple(segments),
            source_id=f"{self.model_id}/pulse-train",
            source_sha256=source,
            interpolation=InterpolationPolicy.HOLD,
            reference_ids=self.reference_ids,
        )

    def _pulse_limit(self, maximum_pulses: int | None) -> int:
        """Validate a caller pulse limit; ``None`` gives the module limit."""
        if maximum_pulses is None:
            return _MAX_PULSES
        limit = integer(maximum_pulses, "Maximum pulses")
        if limit < 1:
            raise DomainError("Maximum pulses must be at least one.")
        return limit

    def _plan(
        self, target: str, request: float, count: int, bit: PulseBit | None,
        initial_mass: float | None, limit: int,
    ) -> PulsePlan:
        """Assemble a `PulsePlan`, computing delivered value and residual.

        The delivered Δv uses the ideal rocket equation [1]_ evaluated with
        ``log1p``; the plan is marked quantized when the residual exceeds a
        relative ``1e-9`` of the request.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
        """
        if count > limit:
            raise DomainError(f"The plan needs {count} pulses, above the limit {limit}.")
        impulse = 0.0 if bit is None else count * bit.impulse_n_s
        propellant = 0.0 if bit is None else count * bit.propellant_kg
        if target == "impulse":
            delivered = impulse
        elif bit is None:
            delivered = 0.0
        elif initial_mass is None:  # every delta-v caller passes the validated mass
            raise DomainError("A delta-v pulse plan requires the initial mass.")
        else:
            delivered = bit.exhaust_velocity_m_s * -log1p(-propellant / initial_mass)
        residual = request - delivered
        scale = max(abs(request), 1e-300)
        return PulsePlan(
            target=target,
            requested=request,
            delivered=delivered,
            residual=residual,
            quantized=abs(residual) > _RELATIVE_TOLERANCE * scale,
            pulse_count=count,
            pulse=bit if count else None,
            delivered_impulse_n_s=impulse,
            propellant_kg=propellant,
            initial_mass_kg=initial_mass,
            final_mass_kg=None if initial_mass is None else initial_mass - propellant,
        )

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the complete deterministic model document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add the ``model_hash`` field.

        Returns
        -------
        dict
            Schema, tables, steady values, limits, provenance, reference
            identifiers and the curve numerical policy.
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "model_id": self.model_id,
            "on_times_s": list(self.on_times_s),
            "curve_numerics": "rational_pchip_and_hermite_export_v1",
            "impulse_bits_n_s": list(self.impulse_bits_n_s),
            "propellant_per_pulse_kg": list(self.propellant_per_pulse_kg),
            "steady_thrust_n": self.steady_thrust_n,
            "steady_mass_flow_kg_s": self.steady_mass_flow_kg_s,
            "minimum_off_time_s": self.minimum_off_time_s,
            "maximum_on_time_s": self.maximum_on_time_s,
            "source_id": self.source_id,
            "source_sha256": self.source_sha256,
            "reference_ids": list(self.reference_ids),
        }
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    @property
    def model_hash(self) -> str:
        """SHA-256 content identity of the pulse model document."""

        return canonical_sha256(self.manifest(include_hash=False))

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete model document as key-sorted, finite JSON.

        Parameters
        ----------
        indent : int or None, default 2
            JSON indentation; ``None`` gives the compact form.

        Returns
        -------
        str
            JSON text without NaN or infinity tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return json.dumps(
            self.manifest(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PulseModeThruster:
        """Rebuild a model from its document and verify its hash.

        Parameters
        ----------
        payload : Mapping
            Document produced by `manifest`.

        Returns
        -------
        PulseModeThruster
            The reconstructed, validated model.

        Raises
        ------
        InputError
            If the root is not an object, a field is missing, malformed or
            violates the model contract, the hash does not match the contents,
            or the numerical policies differ from the current ones.
        """

        if not isinstance(payload, Mapping):
            raise InputError("Pulse-mode document root must be an object.")
        try:
            # The document has no separate structural validator, so each member
            # passes the same normaliser the constructor applies to it; any
            # failure is a DomainError, reported below as a malformed document.
            maximum_on_time = payload["maximum_on_time_s"]
            model = cls(
                on_times_s=real_tuple(payload["on_times_s"], "Pulse on_times_s"),
                impulse_bits_n_s=real_tuple(payload["impulse_bits_n_s"], "Pulse impulse_bits_n_s"),
                propellant_per_pulse_kg=real_tuple(
                    payload["propellant_per_pulse_kg"], "Pulse propellant_per_pulse_kg"
                ),
                steady_thrust_n=real(payload["steady_thrust_n"], "Pulse steady_thrust_n"),
                steady_mass_flow_kg_s=real(
                    payload["steady_mass_flow_kg_s"], "Pulse steady_mass_flow_kg_s"
                ),
                minimum_off_time_s=real(payload["minimum_off_time_s"], "Pulse minimum_off_time_s"),
                maximum_on_time_s=(
                    None if maximum_on_time is None else real(maximum_on_time, "Maximum on-time")
                ),
                source_id=text(payload["source_id"], "Pulse source id"),
                source_sha256=text(payload["source_sha256"], "Pulse source hash"),
                model_id=text(payload["model_id"], "Pulse model id"),
                reference_ids=strings(payload["reference_ids"], "Pulse reference ids"),
                schema=text(payload["schema"], "Pulse-mode schema"),
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed pulse-mode document.") from error
        if payload.get("model_hash") != model.model_hash:
            raise InputError("Pulse-mode hash does not match its contents.")
        _require_numerical_policies(payload, model.manifest(include_hash=False))
        return model

    @classmethod
    def from_json(cls, text_value: str) -> PulseModeThruster:
        """Read and verify a JSON model document.

        Parameters
        ----------
        text_value : str
            JSON text produced by `to_json`.

        Returns
        -------
        PulseModeThruster
            The reconstructed, validated model.

        Raises
        ------
        InputError
            If the text is not a conforming JSON document [1]_ or fails the
            checks of `from_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return _read_document(text_value, "Pulse-mode document", cls.from_dict)


def _policy(value: str) -> str:
    """Normalize and validate a quantization policy name."""
    policy = text(value, "Quantization policy").strip().lower()
    if policy not in QUANTIZATION_POLICIES:
        raise DomainError(f"Quantization policy must be one of {QUANTIZATION_POLICIES}.")
    return policy


def _bisect(function, low: float, high: float) -> float:
    """Return a root of a non-decreasing function on ``[low, high]`` by bisection.

    Returns ``low`` when ``function(low) >= 0`` and ``high`` when
    ``function(high) <= 0``; otherwise halves the bracket until it stops
    shrinking in floating point or 200 halvings are done.
    """
    f_low = function(low)
    if f_low >= 0.0:
        return low
    if function(high) <= 0.0:
        return high
    for _ in range(200):
        middle = 0.5 * (low + high)
        if middle in (low, high):
            break
        if function(middle) < 0.0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)
