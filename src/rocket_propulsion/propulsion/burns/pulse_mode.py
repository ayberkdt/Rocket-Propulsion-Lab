"""Pulse-mode thruster: impulse bits, minimum impulse bit and Δv quantization.

Attitude-control and small trajectory-correction thrusters are fired in short
electrical pulses.  For short pulses the delivered impulse is not
``F_steady * t_on``: valve opening, chamber fill and tail-off make the impulse
bit (and the specific impulse) a nonlinear function of the electrical on-time,
and there is a smallest pulse the valve can command, the **minimum impulse
bit** (MIB).  A trajectory tool therefore cannot ask for an arbitrary small
Δv: it gets an integer number of pulses whose total differs from the request.

:class:`PulseModeThruster` takes *measured* data, an impulse-bit table
``I(t_on)`` and a propellant-per-pulse table ``m(t_on)``, and never invents
generic pulse coefficients.  Between table nodes it uses shape-preserving
PCHIP interpolation (monotone, no overshoot).  Beyond the last node the pulse
is fully developed, and extra on-time adds steady-state impulse and flow
linearly.

For ``n`` identical pulses of on-time ``t`` from mass ``m0`` in free space,
every pulse has the same exhaust velocity ``c = I/m``, so the rocket equation
holds exactly for the whole train:

    Δv_n(t) = (I(t) / m(t)) * ln(m0 / (m0 - n m(t))).

``Δv_n`` is strictly increasing in ``t`` (``I`` strictly increasing, ``m``
non-decreasing) and in ``n``, so a requested Δv is met by the smallest ``n``
whose longest permitted pulse reaches it, and then by bisection on ``t``.
When even ``n`` minimum impulse bits overshoot, the request is quantized and
the residual is reported, never hidden.

References
----------
J. M. Parker and M. J. Wilson, The Minimum Impulse Thruster, 2005,
NASA NTRS 20070032005: https://ntrs.nasa.gov/citations/20070032005
F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
SIAM J. Numer. Anal. 17(2), 1980.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from itertools import pairwise
from math import ceil, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import integer, real, real_tuple, strings, text
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .serialization import canonical_sha256
from .tabulated import InterpolationPolicy, TabulatedBurnArtifact, TabulatedBurnSegment
from .throttleable import MonotoneCubicCurve

PULSE_MODE_SCHEMA = "rocket_propulsion_pulse_mode_thruster_v1"
PULSE_MODE_REFERENCE_IDS = ("NASA-20070032005",)
QUANTIZATION_POLICIES = ("nearest", "no_overshoot")
_MAX_PULSES = 10_000_000


@dataclass(frozen=True, slots=True)
class PulseBit:
    """One pulse of a given electrical on-time.

    References
    ----------
    NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
    """

    on_time_s: float
    impulse_n_s: float
    propellant_kg: float
    impulse_wrt_on_time_n: float
    propellant_wrt_on_time_kg_s: float

    @property
    def specific_impulse_s(self) -> float:
        """Pulse-averaged specific impulse ``I / (m g0)``."""

        return self.impulse_n_s / (self.propellant_kg * STANDARD_GRAVITY_M_S2)

    @property
    def exhaust_velocity_m_s(self) -> float:
        """Pulse-averaged effective exhaust velocity ``I / m``."""

        return self.impulse_n_s / self.propellant_kg


@dataclass(frozen=True, slots=True)
class PulsePlan:
    """An integer pulse train meeting (or quantizing) a request.

    ``requested`` is an impulse in N·s or a Δv in m/s, as given by ``target``.
    ``residual`` is ``requested - delivered`` in the same unit; it is zero
    (to solver tolerance) unless ``quantized`` is ``True``.

    References
    ----------
    NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
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
        """Train-averaged specific impulse, or ``None`` without pulses."""

        return None if self.pulse is None else self.pulse.specific_impulse_s


@dataclass(frozen=True, slots=True)
class PulseModeThruster:
    """Measured pulse-mode characteristics of one thruster.

    ``on_times_s`` are electrical pulse widths, strictly increasing; the first
    is the minimum on-time, so ``impulse_bits_n_s[0]`` is the minimum impulse
    bit.  Impulse bits must be strictly increasing and propellant per pulse
    non-decreasing and positive.  Beyond the last node a pulse is fully
    developed: ``I = I_N + F (t - t_N)``, ``m = m_N + mdot (t - t_N)`` with
    the steady thrust and flow, up to ``maximum_on_time_s`` (default: the last
    node).  The value is continuous there; the slope may change.

    References
    ----------
    NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
    F. N. Fritsch and R. E. Carlson, SIAM J. Numer. Anal. 17(2), 1980.
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

    @property
    def minimum_on_time_s(self) -> float:
        return self.on_times_s[0]

    @property
    def longest_on_time_s(self) -> float:
        """Longest permitted pulse (``maximum_on_time_s`` or the last node)."""

        return self.on_times_s[-1] if self.maximum_on_time_s is None else self.maximum_on_time_s

    @property
    def minimum_impulse_bit_n_s(self) -> float:
        return self.impulse_bits_n_s[0]

    def pulse(self, on_time_s: float) -> PulseBit:
        """Return impulse, propellant and their on-time derivatives for one pulse.

        References
        ----------
        NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
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
        last = self.on_times_s[-1]
        if on_time <= last:
            impulse, impulse_slope = self._impulse.evaluate(on_time)
            mass, mass_slope = self._propellant.evaluate(on_time)
        else:
            extra = on_time - last
            impulse = self.impulse_bits_n_s[-1] + self.steady_thrust_n * extra
            mass = self.propellant_per_pulse_kg[-1] + self.steady_mass_flow_kg_s * extra
            impulse_slope, mass_slope = self.steady_thrust_n, self.steady_mass_flow_kg_s
        return PulseBit(on_time, impulse, mass, impulse_slope, mass_slope)

    def on_time_for_impulse(self, impulse_n_s: float) -> float:
        """Invert the impulse-bit curve (bisection on a monotone function).

        References
        ----------
        NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
        """

        impulse = real(impulse_n_s, "Impulse bit")
        low, high = self.minimum_on_time_s, self.longest_on_time_s
        if impulse < self.minimum_impulse_bit_n_s * (1.0 - 1e-12):
            raise DomainError("Impulse is below the minimum impulse bit.")
        if impulse > self.pulse(high).impulse_n_s * (1.0 + 1e-12):
            raise DomainError("Impulse exceeds the longest permitted pulse.")
        return _bisect(lambda t: self.pulse(t).impulse_n_s - impulse, low, high)

    def train_delta_v_m_s(self, pulse_count: int, on_time_s: float, initial_mass_kg: float) -> float:
        """Free-space Δv of ``pulse_count`` identical pulses (exact rocket equation).

        References
        ----------
        NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
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
        return bit.exhaust_velocity_m_s * log(mass / (mass - consumed))

    def plan_impulse(
        self,
        impulse_n_s: float,
        *,
        quantization: str = "nearest",
        maximum_pulses: int | None = None,
    ) -> PulsePlan:
        """Return the fewest equal pulses whose total impulse meets a request.

        References
        ----------
        NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
        """

        request = real(impulse_n_s, "Requested impulse")
        if request < 0.0:
            raise DomainError("Requested impulse must be non-negative.")
        limit = self._pulse_limit(maximum_pulses)
        policy = _policy(quantization)
        longest = self.pulse(self.longest_on_time_s)
        minimum = self.pulse(self.minimum_on_time_s)
        if request == 0.0:
            return self._plan("impulse", 0.0, 0, None, None)
        count = max(1, ceil(request / longest.impulse_n_s - 1e-12))
        if count > limit:
            raise DomainError(f"The request needs {count} pulses, above the limit {limit}.")
        if request / count >= minimum.impulse_n_s * (1.0 - 1e-12):
            on_time = self.on_time_for_impulse(min(request / count, longest.impulse_n_s))
            return self._plan("impulse", request, count, self.pulse(on_time), None)
        # The request lies in the gap between (count - 1) longest pulses and
        # count minimum impulse bits; no integer train meets it exactly.
        below = (count - 1) * longest.impulse_n_s
        above = count * minimum.impulse_n_s
        if policy == "no_overshoot" or request - below < above - request:
            return self._plan("impulse", request, count - 1, longest if count > 1 else None, None)
        return self._plan("impulse", request, count, minimum, None)

    def plan_delta_v(
        self,
        delta_v_m_s: float,
        initial_mass_kg: float,
        *,
        quantization: str = "nearest",
        maximum_pulses: int | None = None,
    ) -> PulsePlan:
        """Return the fewest equal pulses whose free-space Δv meets a request.

        References
        ----------
        NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
        NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
        """

        request = real(delta_v_m_s, "Requested delta-v")
        mass = real(initial_mass_kg, "Initial mass")
        if request < 0.0 or mass <= 0.0:
            raise DomainError("Delta-v must be non-negative and initial mass positive.")
        limit = self._pulse_limit(maximum_pulses)
        policy = _policy(quantization)
        if request == 0.0:
            return self._plan("delta_v", 0.0, 0, None, mass)
        shortest, longest = self.minimum_on_time_s, self.longest_on_time_s
        longest_bit = self.pulse(longest)
        affordable = ceil(mass / longest_bit.propellant_kg) - 1

        def reach(count: int, on_time: float) -> float:
            bit = self.pulse(on_time)
            if count * bit.propellant_kg >= mass:
                return float("inf")
            return bit.exhaust_velocity_m_s * log(mass / (mass - count * bit.propellant_kg))

        # Smallest count whose longest pulse reaches the request.  Δv_n grows
        # at least linearly in n (Δv_n >= n I / m0), so the first-order count
        # is an upper bound, and Δv_n is increasing in n: bisect on n.
        upper = max(1, min(ceil(request * mass / longest_bit.impulse_n_s), affordable))
        if reach(upper, longest) < request:
            raise DomainError("The request cannot be met with the available mass.")
        lower = 0
        while upper - lower > 1:
            middle = (lower + upper) // 2
            if reach(middle, longest) >= request:
                upper = middle
            else:
                lower = middle
        count = upper
        if count > limit:
            raise DomainError(f"The request needs {count} pulses, above the limit {limit}.")
        if reach(count, shortest) <= request:
            on_time = _bisect(lambda t: reach(count, t) - request, shortest, longest)
            return self._plan("delta_v", request, count, self.pulse(on_time), mass)
        # The request lies in the gap between (count - 1) longest pulses and
        # count minimum impulse bits; no integer train meets it exactly.
        below = reach(count - 1, longest) if count > 1 else 0.0
        above = reach(count, shortest)
        if policy == "no_overshoot" or request - below < above - request:
            return self._plan(
                "delta_v", request, count - 1, longest_bit if count > 1 else None, mass
            )
        return self._plan("delta_v", request, count, self.pulse(shortest), mass)

    def pulse_train_artifact(
        self, plan: PulsePlan, period_s: float, *, phase: str = "pulse"
    ) -> TabulatedBurnArtifact:
        """Return a propagator-ready tabulated burn for a pulse train.

        Each pulse is a rectangle over its electrical on-time with constant
        thrust ``I/t_on`` and flow ``m/t_on``, so impulse, propellant and
        exhaust velocity per pulse are exact; the pulse's internal shape
        (rise and tail-off) is not represented.  Pulses repeat every
        ``period_s``, which must leave the minimum off-time.

        References
        ----------
        NASA-20070032005: https://ntrs.nasa.gov/citations/20070032005
        NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
        """

        if not isinstance(plan, PulsePlan):
            raise DomainError("A PulsePlan is required.")
        if plan.pulse_count < 1 or plan.pulse is None:
            raise DomainError("A pulse train needs at least one pulse.")
        if plan.pulse_count > 100_000:
            raise DomainError("Pulse trains above 100000 pulses are not tabulated.")
        period = real(period_s, "Pulse period")
        on_time = plan.pulse.on_time_s
        if period < on_time + self.minimum_off_time_s - 1e-12 or period <= on_time:
            raise DomainError("The pulse period leaves less than the minimum off-time.")
        phase = text(phase, "Pulse phase")
        thrust = plan.pulse.impulse_n_s / on_time
        flow = plan.pulse.propellant_kg / on_time
        segments: list[TabulatedBurnSegment] = []
        for index in range(plan.pulse_count):
            start = index * period
            segments.append(
                TabulatedBurnSegment(start, start + on_time, thrust, thrust, thrust, thrust,
                                     flow, flow, phase=phase)
            )
            if index + 1 < plan.pulse_count:
                segments.append(
                    TabulatedBurnSegment(start + on_time, (index + 1) * period,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phase="off")
                )
        source = canonical_sha256(
            {
                "model_hash": self.model_hash,
                "pulse_count": plan.pulse_count,
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
        if maximum_pulses is None:
            return _MAX_PULSES
        limit = integer(maximum_pulses, "Maximum pulses")
        if limit < 1:
            raise DomainError("Maximum pulses must be at least one.")
        return limit

    def _plan(
        self, target: str, request: float, count: int, bit: PulseBit | None,
        initial_mass: float | None,
    ) -> PulsePlan:
        impulse = 0.0 if bit is None else count * bit.impulse_n_s
        propellant = 0.0 if bit is None else count * bit.propellant_kg
        if target == "impulse":
            delivered = impulse
        else:
            delivered = 0.0 if bit is None else bit.exhaust_velocity_m_s * log(
                initial_mass / (initial_mass - propellant)
            )
        residual = request - delivered
        scale = max(abs(request), 1e-300)
        return PulsePlan(
            target=target,
            requested=request,
            delivered=delivered,
            residual=residual,
            quantized=abs(residual) > 1e-9 * scale,
            pulse_count=count,
            pulse=bit if count else None,
            delivered_impulse_n_s=impulse,
            propellant_kg=propellant,
            initial_mass_kg=initial_mass,
            final_mass_kg=None if initial_mass is None else initial_mass - propellant,
        )

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the complete deterministic model document."""

        payload: dict[str, object] = {
            "schema": self.schema,
            "model_id": self.model_id,
            "on_times_s": list(self.on_times_s),
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
        """Content identity of the pulse model."""

        return canonical_sha256(self.manifest(include_hash=False))

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete model document."""

        return json.dumps(
            self.manifest(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PulseModeThruster:
        """Rebuild a model and verify its hash."""

        if not isinstance(payload, Mapping):
            raise InputError("Pulse-mode document root must be an object.")
        try:
            model = cls(
                on_times_s=payload["on_times_s"],
                impulse_bits_n_s=payload["impulse_bits_n_s"],
                propellant_per_pulse_kg=payload["propellant_per_pulse_kg"],
                steady_thrust_n=payload["steady_thrust_n"],
                steady_mass_flow_kg_s=payload["steady_mass_flow_kg_s"],
                minimum_off_time_s=payload["minimum_off_time_s"],
                maximum_on_time_s=payload["maximum_on_time_s"],
                source_id=payload["source_id"],
                source_sha256=payload["source_sha256"],
                model_id=payload["model_id"],
                reference_ids=payload["reference_ids"],
                schema=payload["schema"],
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed pulse-mode document.") from error
        if payload.get("model_hash") != model.model_hash:
            raise InputError("Pulse-mode hash does not match its contents.")
        return model

    @classmethod
    def from_json(cls, text_value: str) -> PulseModeThruster:
        """Read and verify a JSON model document."""

        if not isinstance(text_value, str):
            raise InputError("Pulse-mode JSON must be text.")
        try:
            payload = json.loads(text_value)
        except json.JSONDecodeError as error:
            raise InputError("Pulse-mode JSON is malformed.") from error
        return cls.from_dict(payload)


def _policy(value: str) -> str:
    policy = text(value, "Quantization policy").strip().lower()
    if policy not in QUANTIZATION_POLICIES:
        raise DomainError(f"Quantization policy must be one of {QUANTIZATION_POLICIES}.")
    return policy


def _bisect(function, low: float, high: float) -> float:
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
